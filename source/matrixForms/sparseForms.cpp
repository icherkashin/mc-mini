#include <iostream>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"

#define DEBUG true

using namespace Eigen;

namespace SparseForms {
  void makeStokesMatrix (SparseMatrix<double>& stokesMatrix,
                         const int M,
                         const int N,
                         const double h,
                         const double * viscosityData) {
    assert ((M > 0) && (N > 0));
    if (DEBUG) cerr << "Creating stokesMatrix." << endl;

    vector<Triplet<double> > tripletList;

    makeLaplacianXBlock (tripletList, 0,                 0,                 M, N, h, viscosityData);
    makeLaplacianYBlock (tripletList, M * (N - 1),       M * (N - 1),       M, N, h, viscosityData);
    makeGradXBlock      (tripletList, 0,                 2 * M * N - M - N, M, N, h);
    makeGradYBlock      (tripletList, M * (N - 1),       2 * M * N - M - N, M, N, h);
    makeDivXBlock       (tripletList, 2 * M * N - M - N, 0,                 M, N, h);
    makeDivYBlock       (tripletList, 2 * M * N - M - N, M * (N - 1),       M, N, h);
    
    stokesMatrix.setFromTriplets (tripletList.begin(), tripletList.end());
  }

  void makeLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosityData) {
    assert ((M > 0) && (N > 0));
    if (DEBUG) cerr << "Creating LaplacianXBlock." << endl;

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < (N - 1); ++j) {
        double viscosity = (viscosityData[i * (N + 1) + (j + 1)] + 
                            viscosityData[(i + 1) * (N + 1) + (j + 1)]) / 2;
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
    assert ((M > 0) && (N > 0));
    if (DEBUG) cerr << "Creating LaplacianYBlock." << endl;

    for (int i = 0; i < (M - 1); ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityData[(i + 1) * (N + 1) + j] + 
                            viscosityData[(i + 1) * (N + 1) + (j + 1)]) / 2;
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
    if (DEBUG) cerr << "Creating GradXBlock." << endl;
    for (int i = 0; i < M; ++i) {
      for (int x = 0; x < N - 1; ++x) {
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * N + x,     -1 / h));
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + x, N0 + i * N + 1 + x,  1 / h));
      }
    }
  }

  void makeGradYBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h) {
    if (DEBUG) cerr << "Creating GradYBlock." << endl;
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
    if (DEBUG) cerr << "Creating DivXBlock." << endl;
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
    if (DEBUG) cerr << "Creating DivYBlock." << endl;
    for (int i = 0; i < (M - 1) * N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,     N0 + i,  1 / h));
      tripletList.push_back (Triplet<double> (M0 + i + N, N0 + i, -1 / h));
    }
  }

  void makeForcingMatrix (SparseMatrix<double>& forcingMatrix,
                          const int M,
                          const int N) {
    if (DEBUG) cerr << "Creating ForcingMatrix." << endl;
    vector<Triplet<double> > tripletList;

    for (int i = 0; i < 2 * M * N - M - N; ++i)
      tripletList.push_back (Triplet<double> (i, i, 1));
  
    forcingMatrix.setFromTriplets (tripletList.begin (), tripletList.end ());
  }

  void makeBoundaryMatrix (SparseMatrix<double>& boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double * viscosityData) {
    if (DEBUG) cerr << "Creating BoundaryMatrix." << endl;
    vector<Triplet<double> > tripletList;

    makeBCLaplacianXBlock (tripletList, 0,                 0,     M, N, h, viscosityData);
    makeBCLaplacianYBlock (tripletList, M * (N - 1),       2 * M, M, N, h, viscosityData);
    makeBCDivXBlock       (tripletList, 2 * M * N - M - N, 0,     M, N, h);
    makeBCDivYBlock       (tripletList, 2 * M * N - M - N, 2 * M, M, N, h);

    boundaryMatrix.setFromTriplets (tripletList.begin(), tripletList.end());
  }

  void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                              const int      M0,
                              const int      N0,
                              const int      M,
                              const int      N,
                              const double   h,
                              const double * viscosityData) {
    assert ((M > 0) && (N > 0));
    if (DEBUG) cerr << "Creating BCLaplacianXBlock." << endl;

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < 2; ++j) {
        double viscosity = (viscosityData[i       * (N + 1) + j * N] + 
                            viscosityData[(i + 1) * (N + 1) + j * N]) / 2;
        tripletList.push_back (Triplet<double> (M0 + i * (N - 1) + j * (N - 2),
                                                N0 + i * 2       + j,
                                                viscosity / (h * h)));
      }
    }
    /*
    for (int i = 0; i < M; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i * (N - 1),           N0 + i * 2,     viscosityData[0] / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + (i + 1) * (N - 1) - 1, N0 + i * 2 + 1, viscosityData[0] / (h * h)));
    }
    */
  }

  void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
                              const int      M0,
                              const int      N0,
                              const int      M,
                              const int      N,
                              const double   h,
                              const double * viscosityData) {
    if (DEBUG) cerr << "Creating BCLaplacianYBlock." << endl;

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < N; ++j) {
        double viscosity = (viscosityData[(i * M) * (N + 1) + j] +
                            viscosityData[(i * M) * (N + 1) + (j + 1)]) / 2;
        tripletList.push_back (Triplet<double> (M0 + i * (M - 2) * N + j,  
                                                N0 + i           * N + j, 
                                                viscosity / (h * h)));
      }
    }
    /*
    for (int i = 0; i < N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,               N0 + i,     viscosityData[0] / (h * h)));
      tripletList.push_back (Triplet<double> (M0 + (M - 2) * N + i, N0 + N + i, viscosityData[0] / (h * h)));
    }
    */
  }

  void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h) {
    if (DEBUG) cerr << "Creating BCDivXBlock." << endl;
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
    if (DEBUG) cerr << "Creating BCDivYBlock." << endl;
    for (int i = 0; i < N; ++i) {
      tripletList.push_back (Triplet<double> (M0 + i,               N0 + i,      1 / h));
      tripletList.push_back (Triplet<double> (M0 + (M - 1) * N + i, N0 + N + i, -1 / h));
    }
  }
}
