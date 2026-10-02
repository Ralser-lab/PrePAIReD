FROM mambaorg/micromamba:latest

WORKDIR /workflow

COPY environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml -n unifiedenv && \
    micromamba clean --all --yes

# CRAN packages not taken from conda. 'iq' provides the MaxLFQ implementation
# that replaces the former C++ code. install.packages() only warns on failure,
# so the build is stopped explicitly if any package is still missing.
RUN micromamba run -n unifiedenv R -e "\
  pkgs <- c('iq','mice','missForest'); \
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]; \
  if(length(to_install)) install.packages(to_install, repos='https://cloud.r-project.org'); \
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; \
  if(length(missing)) stop('R packages failed to install: ', paste(missing, collapse=', ')); \
  cat('R packages OK:', paste(pkgs, collapse=', '), '\n')"

ENV PATH="/opt/conda/envs/unifiedenv/bin:$PATH"
ENV TZ=Europe/Berlin

COPY Initialization/ ./Initialization/
COPY modules_python/ ./modules_python/
COPY modules_R/ ./modules_R/
COPY docs/ ./docs/

# Give the (non-root) micromamba user ownership of /workflow
USER root
RUN chown -R $MAMBA_USER:$MAMBA_USER /workflow
USER $MAMBA_USER

CMD ["bash"]