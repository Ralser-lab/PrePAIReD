FROM mambaorg/micromamba:latest

WORKDIR /workflow

COPY environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml -n unifiedenv && \
    micromamba clean --all --yes

RUN micromamba run -n unifiedenv R -e "\
  pkgs <- c('iq','mice','missForest'); \
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]; \
  if(length(to_install)) install.packages(to_install, repos='https://cloud.r-project.org'); \
  cat('Installed:', paste(to_install, collapse=', '), '\n')"

ENV PATH="/opt/conda/envs/unifiedenv/bin:$PATH"
ENV TZ=Europe/Berlin

COPY Initialization/ ./Initialization/
COPY modules_python/ ./modules_python/
COPY modules_R/ ./modules_R/
COPY docs/ ./docs/

# Fix permissions so the (non-root) micromamba user can write build artifacts
USER root
RUN chown -R $MAMBA_USER:$MAMBA_USER /workflow
USER $MAMBA_USER

# Install a C++ compiler into the unifiedenv (do NOT wrap micromamba install in micromamba run)
RUN micromamba install -n unifiedenv -y -c conda-forge cxx-compiler eigen && \
    micromamba clean --all --yes

RUN micromamba run -n unifiedenv bash -lc '\
  set -eux; \
  cd /workflow/modules_python; \
  CXX="$(ls -1 $CONDA_PREFIX/bin/*-conda-linux-gnu-g++ 2>/dev/null | head -n 1 || true)"; \
  [ -z "$CXX" ] && CXX="$(command -v g++ || true)"; \
  [ -z "$CXX" ] && echo "No C++ compiler found" && exit 1; \
  "$CXX" -O3 -fPIC -I"$CONDA_PREFIX/include/eigen3" -c maxlfq_docker.cpp -o maxlfq.o; \
  "$CXX" -shared -o maxlfq_docker.so maxlfq.o \
'

CMD ["bash"]