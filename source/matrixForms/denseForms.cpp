#include <iostream>

#include <Eigen/Dense>

#include "matrixForms/denseForms.h"
#include "debug.h"

using namespace Eigen;
using namespace std;

namespace DenseForms {
  void makeStokesMatrix (Ref<MatrixXd> stokesMatrix,
                         const int M,
                         const int N,
                         const double h,
                         const double * viscosityData) {
    if (DEBUG) cerr << "Creating stokesMatrix." << endl;

    stokesMatrix = MatrixXd::Zero (3 * M * N - M - N, 3 * M * N - M - N);

    makeLaplacianXBlock (stokesMatrix.block (0,                 0,                 M * (N - 1), M * (N - 1)), M, N, h, viscosityData);
    makeLaplacianYBlock (stokesMatrix.block (M * (N - 1),       M * (N - 1),       (M - 1) * N, (M - 1) * N), M, N, h, viscosityData);
    makeGradXBlock      (stokesMatrix.block (0,                 2 * M * N - M - N, M * (N - 1), M * N),       M, N, h);
    makeGradYBlock      (stokesMatrix.block (M * (N - 1),       2 * M * N - M - N, (M - 1) * N, M * N),       M, N, h);
    makeDivXBlock       (stokesMatrix.block (2 * M * N - M - N, 0,                  M * N,      M * (N - 1)), M, N, h);
    makeDivYBlock       (stokesMatrix.block (2 * M * N - M - N, M * (N - 1),        M * N,      (M - 1) * N), M, N, h);
  }

  void makeLaplacianXBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosityData) {
    if (DEBUG) cerr << "Creating LaplacianXBlock." << endl;

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        double viscosity = (viscosityData [i * (N + 1) + (j + 1)] +
                            viscosityData [(i + 1) * (N + 1) + (j + 1)]) / 2;

        if (i == 0 || i == (M - 1))
          laplacian (i * (N - 1) + j, i       * (N - 1) + j)       =  viscosity * 5 / (h * h);
        else
          laplacian (i * (N - 1) + j, i       * (N - 1) + j)       =  viscosity * 4 / (h * h);

        if (i > 0)
          laplacian (i * (N - 1) + j, (i - 1) * (N - 1) + j)       = -viscosity / (h * h);

        if (i < (M - 1))
          laplacian (i * (N - 1) + j, (i + 1) * (N - 1) + j)       = -viscosity / (h * h);

        if (j > 0)
          laplacian (i * (N - 1) + j, i       * (N - 1) + (j - 1)) = -viscosity / (h * h);

        if (j < N - 2)
          laplacian (i * (N - 1) + j, i       * (N - 1) + (j + 1)) = -viscosity / (h * h);
      }
    }
  }

  void makeLaplacianYBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosityData) {
    if (DEBUG) cerr << "Creating LaplacianYBlock." << endl;

    for (int i = 0; i < (M - 1); ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityData [(i + 1) * (N + 1) + j] +
                            viscosityData [(i + 1) * (N + 1) + (j + 1)]) / 2;

        if ((j == 0) || (j == (N - 1)))
          laplacian (i * N + j, i       * N + j)       =  viscosity * 5 / (h * h);
        else
          laplacian (i * N + j, i       * N + j)       =  viscosity * 4 / (h * h);

        if (j > 0)
          laplacian (i * N + j, i       * N + (j - 1)) = -viscosity / (h * h);

        if (j < (N - 1))
          laplacian (i * N + j, i       * N + (j + 1)) = -viscosity / (h * h);

        if (i > 0)
          laplacian (i * N + j, (i - 1) * N + j)       = -viscosity / (h * h);

        if (i < (M - 2))
          laplacian (i * N + j, (i + 1) * N + j)       = -viscosity / (h * h);
      }
    }
  }

  void makeGradXBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradXBLock." << endl;
    
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        grad (i * (N - 1) + j, i * N + j)       = -1 / h;
        grad (i * (N - 1) + j, i * N + (j + 1)) =  1 / h;
      }
    }
  }

  void makeGradYBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradYBlock." << endl;
    
    for (int i = 0; i < (M - 1) * N; ++i) {
      grad (i, i)     = -1 / h;
      grad (i, i + N) =  1 / h;
    }
  }

  void makeDivXBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivXBLock." << endl;

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        div (i * N + j,       i * (N - 1) + j) =  1 / h;
        div (i * N + (j + 1), i * (N - 1) + j) = -1 / h;
      }
    }
  }

  void makeDivYBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h) {
    if (DEBUG) cerr << "Creating DivYBlock." << endl;
   
    for (int i = 0; i < (M - 1) * N; ++i) {
      div (i,     i) =  1 / h;
      div (i + N, i) = -1 / h;
    }
  }

  void makeForcingMatrix (Ref<MatrixXd> forcingMatrix,
                          const int M,
                          const int N) {
    if (DEBUG) cerr << "Creating ForcingMatrix." << endl;

    forcingMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M * N - M - N);

    for (int i = 0; i < 2 * M * N - M - N; ++i)
      forcingMatrix (i, i) = 1;
  }

  void makeBoundaryMatrix (Ref<MatrixXd> boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double * viscosityData) {
    if (DEBUG) cerr << "Creating BoundaryMatrix." << endl;
    
    boundaryMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M + 2 * N);

    makeBCLaplacianXBlock (boundaryMatrix.block (0,                 0,     M * (N - 1), 2 * M), M, N, h, viscosityData);
    makeBCLaplacianYBlock (boundaryMatrix.block (M * (N - 1),       2 * M, (M - 1) * N, 2 * N), M, N, h, viscosityData);
    makeBCDivXBlock       (boundaryMatrix.block (2 * M * N - M - N, 0,     M * N,       2 * M), M, N, h);
    makeBCDivYBlock       (boundaryMatrix.block (2 * M * N - M - N, 2 * M, M * N,       2 * N), M, N, h);
  }

  void makeBCLaplacianXBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosityData) {
    if (DEBUG) cerr << "Creating BCLaplacianXBlock." << endl;
   
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < 2; ++j) {
        double viscosity = (viscosityData [i *       (N + 1) + j * N] +
                            viscosityData [(i + 1) * (N + 1) + j * N]) / 2;
        laplacianBC (i * (N - 1) + j * (N - 2), i * 2 + j) = viscosity / (h * h);
      }
    }
  }

  void makeBCLaplacianYBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosityData) {
    if (DEBUG) cerr << "Creating BCLaplacianYBlock." << endl;
  
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityData [(i * M) * (N + 1) + j] +
                            viscosityData [(i * M) * (N + 1) + (j - 1)]) / 2;
        laplacianBC (i * (M - 2) * N + j, i * N + j) = viscosity / (h * h);
      }
    }
  }

  void makeBCDivXBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivXBlock." << endl;
  
    for (int i = 0; i < M; ++i) {
      divBC (i *       N,     i * 2)     =  1 / h;
      divBC ((i + 1) * N - 1, i * 2 + 1) = -1 / h;
    }
  }

  void makeBCDivYBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivYBlock." << endl;
  
    for (int i = 0; i < N; ++i) {
      divBC (              i,     i) =  1 / h;
      divBC ((M - 1) * N + i, N + i) = -1 / h;
    }
  }
}

  
