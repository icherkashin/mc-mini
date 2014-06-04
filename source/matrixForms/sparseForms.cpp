#include <iostream>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "debug.h"

using namespace Eigen;

namespace SparseForms {
  void makeStokesMatrix (SparseMatrix<double>& stokesMatrix,
                         const int M,
                         const int N,
                         const double h,
                         const double * viscosityData) {
    #ifdef DEBUG 
      cout << "<Creating " << 3 * M * N - M - N << "x" << 3 * M * N - M - N << " stokesMatrix>" << endl;
    #endif

    vector<Triplet<double> > tripletList;

    makeLaplacianXBlock (tripletList, 0,                 0,                 M, N, h, viscosityData);
    makeLaplacianYBlock (tripletList, M * (N - 1),       M * (N - 1),       M, N, h, viscosityData);
    makeGradXBlock      (tripletList, 0,                 2 * M * N - M - N, M, N, h);
    makeGradYBlock      (tripletList, M * (N - 1),       2 * M * N - M - N, M, N, h);
    makeDivXBlock       (tripletList, 2 * M * N - M - N, 0,                 M, N, h);
    makeDivYBlock       (tripletList, 2 * M * N - M - N, M * (N - 1),       M, N, h);
    
    stokesMatrix.setFromTriplets (tripletList.begin(), tripletList.end());
    #ifdef DEBUG
      cout << endl;
    #endif
  }

  void makeLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
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
          tripletList.push_back (
              Triplet<double> (M0 + i * (N - 1) + j, 
                               N0 + i * (N - 1) + j, 
                               viscosity * 5 / (h * h)));
        else
          tripletList.push_back (
              Triplet<double> (M0 + i * (N - 1) + j, 
                               N0 + i * (N - 1) + j, 
                               viscosity * 4 / (h * h)));

        // First and last rows are missing a neighbor in one of two directions
        if (i > 0)
          tripletList.push_back (
              Triplet<double> (M0 + i       * (N - 1) + j, 
                               N0 + (i - 1) * (N - 1) + j, 
                               -viscosity / (h * h)));
        if (i < (M - 1))
          tripletList.push_back (
              Triplet<double> (M0 + i       * (N - 1) + j, 
                               N0 + (i + 1) * (N - 1) + j, 
                               -viscosity / (h * h)));

        // First and last elements of each row are missing a neighbor in one of two directions 
        if (j > 0)
          tripletList.push_back (
              Triplet<double> (M0 + i * (N - 1) + j, 
                               N0 + i * (N - 1) + j - 1, 
                               -viscosity / (h * h)));
        if (j < (N - 2))
          tripletList.push_back (
              Triplet<double> (M0 + i * (N - 1) + j, 
                               N0 + i * (N - 1) + j + 1, 
                               -viscosity / (h * h)));
      }
    }
  }

  void makeLaplacianYBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
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
          tripletList.push_back (Triplet<double> (M0 + i * N + j,
                                                  N0 + i * N + j,
                                                  viscosity * 5 / (h * h)));
        else
          tripletList.push_back (Triplet<double> (M0 + i * N + j,
                                                  N0 + i * N + j,
                                                  viscosity * 4 / (h * h)));

        // First and last elements of each row are missing a neighbor in one of two directions
        if (j > 0)
          tripletList.push_back (Triplet<double> (M0 + i * N + j,
                                                  N0 + i * N + (j - 1),
                                                  -viscosity / (h * h)));
        if (j < (N - 1))
          tripletList.push_back (Triplet<double> (M0 + i * N + j,
                                                  N0 + i * N + (j + 1),
                                                  -viscosity / (h * h)));

        // Elements of the first and last rows are missing a neighbor in one of two directions
        if (i > 0)
          tripletList.push_back (Triplet<double> (M0 + i       * N + j,
                                                  N0 + (i - 1) * N + j,
                                                  -viscosity / (h * h)));
        if (i < (M - 2))
          tripletList.push_back (Triplet<double> (M0 + i       * N + j,
                                                  N0 + (i + 1) * N + j,
                                                  -viscosity / (h * h)));
      }
    }
  }

  void makeGradXBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h) {
    #ifdef DEBUG
      cout << "<Creating " << M * (N - 1) << "x" << M * N << " GradXBlock>" << endl;
    #endif

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + j, N0 + i * N + j,     -1 / h));
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + j, N0 + i * N + (j + 1),  1 / h));
      }
    }
  }

  void makeGradYBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h) {
    #ifdef DEBUG
      cout << "<Creating " << (M - 1) * N << "x" << M * N << " GradYBlock>" << endl;
    #endif

    for (int i = 0; i < (M - 1) * N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i, N0 + i,     -1 / h));
      tripletList.push_back (Triplet<double> (M0 + i, N0 + N + i,  1 / h));
    }
  }

  void makeDivXBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << (M - 1) * N << " DivXBlock>" << endl;
    #endif

    for (int i = 0; i < M; ++i) {
      for (int x = 0; x < (N - 1); ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * N + x,     N0 + i * (N - 1) + x,  1 / h));
        tripletList.push_back (Triplet<double> (M0 + i * N + 1 + x, N0 + i * (N - 1) + x, -1 / h));
      }
    }
  }

  void makeDivYBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h) {
    #ifdef DEBUG
      cout << "<Creating " << M * N << "x" << (M - 1) * N << " DivYBlock>" << endl;
    #endif

    for (int i = 0; i < (M - 1) * N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,     N0 + i,  1 / h));
      tripletList.push_back (Triplet<double> (M0 + i + N, N0 + i, -1 / h));
    }
  }

  void makeForcingMatrix (SparseMatrix<double>& forcingMatrix,
                          const int M,
                          const int N) {
    #ifdef DEBUG
      cout << "<Creating " << 3 * M * N - M - N << "x" << 2 * M * N - M - N << " ForcingMatrix>" << endl;
    #endif
    vector<Triplet<double> > tripletList;

    for (int i = 0; i < 2 * M * N - M - N; ++i)
      tripletList.push_back (Triplet<double> (i, i, 1));
  
    forcingMatrix.setFromTriplets (tripletList.begin (), tripletList.end ());
    
    #ifdef DEBUG
      cout << endl;
    #endif
  }

  void makeBoundaryMatrix (SparseMatrix<double>& boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double * viscosityData) {
    #ifdef DEBUG 
      cout << "<Creating " << 3 * M * N - M - N << "x" << 2 * M + 2 * N << " BoundaryMatrix>" << endl;
    #endif

    vector<Triplet<double> > tripletList;

    makeBCLaplacianXBlock (tripletList, 0,                 0,     M, N, h, viscosityData);
    makeBCLaplacianYBlock (tripletList, M * (N - 1),       2 * M, M, N, h, viscosityData);
    makeBCDivXBlock       (tripletList, 2 * M * N - M - N, 0,     M, N, h);
    makeBCDivYBlock       (tripletList, 2 * M * N - M - N, 2 * M, M, N, h);

    boundaryMatrix.setFromTriplets (tripletList.begin(), tripletList.end());

    #ifdef DEBUG
      cout << endl;
    #endif
  }

  void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                              const int      M0,
                              const int      N0,
                              const int      M,
                              const int      N,
                              const double   h,
                              const double * viscosityData) {
    #ifdef DEBUG
      cout << "<Creating " << (M - 1) * N << "x" << 2 * M << " BCLaplacianXBlock>" << endl;
    #endif

    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < 2; ++j) {
        double viscosity = (viscosityWindow (j * N, i) + viscosityWindow (j * N, i + 1)) / 2;
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + j * (N - 2),
                                                N0 + i * 2       + j,
                                                viscosity / (h * h)));
      }
    }
  }

  void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
                              const int      M0,
                              const int      N0,
                              const int      M,
                              const int      N,
                              const double   h,
                              const double * viscosityData) {
    #ifdef DEBUG
      cout << "<Creating " << (M - 1) * N << "x" << 2 * N << " BCLaplacianYBlock>" << endl;
    #endif

    DataWindow<const double> viscosityWindow (viscosityData, N + 1, M + 1);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityWindow (j, i * M) + viscosityWindow (j + 1, i * M)) / 2;
        tripletList.push_back (Triplet<double> (M0 + i * (M - 2) * N + j,  
                                                N0 + i           * N + j, 
                                                viscosity / (h * h)));
      }
    }
  }

  void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h) {
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << 2 * M << " BCDivXBlock>" << endl;
    #endif

    for (int i = 0; i < M; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i * N,           N0 + i * 2,      1 / h));
      tripletList.push_back (Triplet<double> (M0 + (i + 1) * N - 1, N0 + i * 2 + 1, -1 / h));
    }
  }

  void makeBCDivYBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h) {
    #ifdef DEBUG 
      cout << "<Creating " << M * N << "x" << 2 * N << " BCDivYBlock>" << endl;
    #endif

    for (int i = 0; i < N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,               N0 + i,      1 / h));
      tripletList.push_back (Triplet<double> (M0 + (M - 1) * N + i, N0 + N + i, -1 / h));
    }
  }
}
