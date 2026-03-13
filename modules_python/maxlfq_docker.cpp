#include <iostream>
#include <vector>
#include<typeinfo>
#include <Eigen/Dense>
#include<thread>

#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

static std::vector<Eigen::MatrixXd> vA(1024);
static std::vector<Eigen::VectorXd> vB(1024);
static std::vector<Eigen::VectorXd> vX(1024);


using namespace std;

extern "C"{
    void free_memory(double* a)
    {
     free(a);
    }
}


extern "C"{
    double* maxlfq (
        int thread_id,  // must be between 0 and 1023; two threads must not use the same thread_id simultaneously
        double* reference_intensities,  // protein intensities calculated using some simple non-MaxLFQ algorithm; vector of size n
        double* ratios,  // matrix of the size n x n, stored row-wise
        double* weights,  // matrix of the size n x n, stored row-wise
        int n) {  // number of runs

    int i, j, k;
    auto& A = vA[thread_id];
    auto& B = vB[thread_id];
    auto& X = vX[thread_id];

    double* intensities_p = (double*)malloc(n*sizeof(double));

    if (intensities_p == NULL)
      {
        cout << "Null pointer has been returned";
      }

    if (B.size() != n) {
        A.resize(n, n);
        B.resize(n);
        X.resize(n);
    }

    for (i = 0; i < n; i++) {
        B(i) = 0.0;
        for (j = 0; j < n; j++) A(i, j) = 0.0;
    }

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            auto ratio = ratios[i * n + j];

            auto weight = weights[i * n + j];

            if (weight > 1e-9) {
                A(i, i) += weight;
                A(j, j) += weight;
                A(i, j) = A(j, i) = -weight;
                B(i) += ratio * weight;
                B(j) -= ratio * weight;
            }
        }
    }


    // small weights towards approximating the 'reference' intensities
    for (i = 0; i < n; i++) {
        double reg = 0.0001 * (A(i, i) >= 1.0 ? A(i, i) : 1.0);
        A(i, i) += reg;
        B(i) += reference_intensities[i] * reg;

    }

    X.noalias() = A.llt().solve(B);

    for (int i = 0; i < n; i++) {
       intensities_p[i] = X[i];
    }

    return intensities_p;
    }
}
