#pragma once

#include <vector>

#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

namespace SparseForms {
  void makeStokesMatrix (SparseMatrix<double>& stokesMatrix,
                         const int M,
                         const int N,
                         const double h,
                         const double * viscosity);

  void makeLaplacianXBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosity);

  void makeLaplacianYBlock (vector<Triplet<double> >& tripletList,
                            const int M0,
                            const int N0,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosity);

  void makeGradXBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h);

  void makeGradYBlock (vector<Triplet<double> >& tripletList,
                       const int M0,
                       const int N0,
                       const int M,
                       const int N,
                       const double h);

  void makeDivXBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h);

  void makeDivYBlock (vector<Triplet<double> >& tripletList,
                      const int M0,
                      const int N0,
                      const int M,
                      const int N,
                      const double h);

  void makeForcingMatrix (SparseMatrix<double>& forcingMatrix,
                          const int M,
                          const int N);

  void makeBoundaryMatrix (SparseMatrix<double>& boundaryMatrix,
                           const int N,
                           const int M,
                           const double h,
                           const double * viscosity);

  void makeBCLaplacianXBlock (vector<Triplet<double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosity);

  void makeBCLaplacianYBlock (vector<Triplet<double> >& tripletList,
                              const int M0,
                              const int N0,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosity);

  void makeBCDivXBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h);

  void makeBCDivYBlock (vector<Triplet<double> >& tripletList,
                        const int M0,
                        const int N0,
                        const int M,
                        const int N,
                        const double h);

  void makeForwardEulerMatrix (SparseMatrix<double> matrix,
                               const int M,
                               const int N,
                               const double dt,
                               const double diffusivity,
                               const double h);
}
