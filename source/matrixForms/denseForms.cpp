#include <iostream>

#include <Eigen/Dense>

#include "geometry/dataWindow.h"
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
    #ifdef DEBUG
      cout << "<Creating " << 3 * M * N - M - N << "x" << 3 * M * N - M - N << " stokesMatrix>" << endl << endl;
    #endif
    stokesMatrix = MatrixXd::Zero (3 * M * N - M - N, 3 * M * N - M - N);

    makeLaplacianXBlock (stokesMatrix.block (0,                 0,                 M * (N - 1), M * (N - 1)), M, N, h, viscosityData);
    makeLaplacianYBlock (stokesMatrix.block (M * (N - 1),       M * (N - 1),       (M - 1) * N, (M - 1) * N), M, N, h, viscosityData);
    makeGradXBlock      (stokesMatrix.block (0,                 2 * M * N - M - N, M * (N - 1), M * N),       M, N, h);
    makeGradYBlock      (stokesMatrix.block (M * (N - 1),       2 * M * N - M - N, (M - 1) * N, M * N),       M, N, h);
    makeDivXBlock       (stokesMatrix.block (2 * M * N - M - N, 0,                  M * N,      M * (N - 1)), M, N, h);
    makeDivYBlock       (stokesMatrix.block (2 * M * N - M - N, M * (N - 1),        M * N,      (M - 1) * N), M, N, h);
  
    #ifdef DEBUG
      cout << endl;
    #endif
  }

  void makeLaplacianXBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosityData) {
    #ifdef DEBUG 
      cout << "<Creating " << M * (N - 1) << "x" << M * (N - 1) << " LaplacianXBlock>" << endl;
    #endif
    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1); 

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        double viscosity = (viscosityWindow (j + 1, i) + viscosityWindow (j + 1, i + 1)) / 2;

        // First and last rows are non-standard because the laplacian would sample points which
        // do not exist in our gridding.
        if (i == 0 || i == (M - 1))
          laplacian (i * (N - 1) + j, i       * (N - 1) + j)       =  viscosity * 5 / (h * h);
        else
          laplacian (i * (N - 1) + j, i       * (N - 1) + j)       =  viscosity * 4 / (h * h);
        // First and last rows are missing a neighbor in one of two directions
        if (i > 0)
          laplacian (i * (N - 1) + j, (i - 1) * (N - 1) + j)       = -viscosity / (h * h);
        if (i < (M - 1))
          laplacian (i * (N - 1) + j, (i + 1) * (N - 1) + j)       = -viscosity / (h * h);
        // First and last elements of each row are missing a neighbor in one of two directions
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
    #ifdef DEBUG
      cout << "<Creating " << (M - 1) * N << "x" << (M - 1) * N << " LaplacianYBlock>" << endl;
    #endif

    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1);
    for (int i = 0; i < (M - 1); ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityWindow (j, i + 1) + viscosityWindow (j + 1, i + 1)) / 2;

        // The first and last elements of each row are non-standard because the four-point
        // laplacian relies upon points not included in our gridding
        if ((j == 0) || (j == (N - 1)))
          laplacian (i * N + j, i       * N + j)       =  viscosity * 5 / (h * h);
        else
          laplacian (i * N + j, i       * N + j)       =  viscosity * 4 / (h * h);

        // First and last elements of each row are missing a neighbor in one of two directions
        if (j > 0)
          laplacian (i * N + j, i       * N + (j - 1)) = -viscosity / (h * h);
        if (j < (N - 1))
          laplacian (i * N + j, i       * N + (j + 1)) = -viscosity / (h * h);

        // Elements of the first and last rows are missing a neighbor in one of two directions
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
    #ifdef DEBUG
      cout << "<Creating " << M * (N - 1) << "x" << M * N << " GradXBLock>" << endl;
    #endif
    
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
    #ifdef DEBUG
      cout << "<Creating " << (M - 1) * N << "x" << M * N << " GradYBlock>" << endl;
    #endif

    for (int i = 0; i < (M - 1) * N; ++i) {
      grad (i, i)     = -1 / h;
      grad (i, i + N) =  1 / h;
    }
  }

  void makeDivXBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h) {
    #ifdef DEBUG
      cout << "<Creating " << M * N << "x" << M * (N - 1) << " DivXBLock>" << endl;
    #endif

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
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << (M - 1) * N << " DivYBlock>" << endl;
    #endif
   
    for (int i = 0; i < (M - 1) * N; ++i) {
      div (i,     i) =  1 / h;
      div (i + N, i) = -1 / h;
    }
  }

  void makeForcingMatrix (Ref<MatrixXd> forcingMatrix,
                          const int M,
                          const int N) {
    #ifdef DEBUG
      cout << "<Creating " << 3 * M * N - M - N << "x" << 2 * M * N - M - N << " ForcingMatrix>" << endl;
    #endif

    forcingMatrix = MatrixXd::Zero (3 * M * N - M - N, 2 * M * N - M - N);

    for (int i = 0; i < 2 * M * N - M - N; ++i)
      forcingMatrix (i, i) = 1;

    #ifdef DEBUG
      cout << endl;
    #endif
  }

  void makeBoundaryMatrix (Ref<MatrixXd> boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double * viscosityData) {
    #ifdef DEBUG 
      cout << "<Creating " << 3 * M * N - M - N << "x" << 2 * M + 2 * N << " BoundaryMatrix>" << endl;
    #endif
    
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
    #ifdef DEBUG 
      cout << "<Creating " << M * (N - 1) << "x" << 2 * M << " BCLaplacianXBlock>" << endl;
    #endif
   
    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < 2; ++j) {
        double viscosity = (viscosityWindow (j * N, i) + viscosityWindow (j * N, i + 1)) / 2;
        laplacianBC (i * (N - 1) + j * (N - 2), i * 2 + j) = viscosity / (h * h);
      }
    }
  }

  void makeBCLaplacianYBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosityData) {
    #ifdef DEBUG 
      cout << "<Creating " << (M - 1) * N << "x" << 2 * N << " BCLaplacianYBlock>" << endl;
    #endif

    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityWindow (j, i * M) + viscosityWindow (j - 1, i * M)) / 2;
        laplacianBC (i * (M - 2) * N + j, i * N + j) = viscosity / (h * h);
      }
    }
  }

  void makeBCDivXBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << 2 * M << " BCDivXBlock>" << endl;
    #endif
  
    for (int i = 0; i < M; ++i) {
      divBC (i *       N,     i * 2)     =  1 / h;
      divBC ((i + 1) * N - 1, i * 2 + 1) = -1 / h;
    }
  }

  void makeBCDivYBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h) {
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << 2 * N << " BCDivYBlock>" << endl;
    #endif  

    for (int i = 0; i < N; ++i) {
      divBC (              i,     i) =  1 / h;
      divBC ((M - 1) * N + i, N + i) = -1 / h;
    }
  }
}

  
