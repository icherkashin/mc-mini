#pragma once

#include <Eigen/Dense>

using namespace Eigen;

namespace DenseForms {
  void makeStokesMatrix (Ref<MatrixXd> stokesMatrix,
                         const int M, 
                         const int N,
                         const double h,
                         const double * viscosity);

  void makeLaplacianXBlock (Ref<MatrixXd> laplacian,
                            const int M, 
                            const int N, 
                            const double h, 
                            const double * viscosity);
  void makeLaplacianYBlock (Ref<MatrixXd> laplacian,
                            const int M,
                            const int N,
                            const double h,
                            const double * viscosity);

  void makeGradXBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h);
  void makeGradYBlock (Ref<MatrixXd> grad,
                       const int M,
                       const int N,
                       const double h);

  void makeDivXBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h);
  void makeDivYBlock (Ref<MatrixXd> div,
                      const int M,
                      const int N,
                      const double h);

  void makeForcingMatrix (Ref<MatrixXd> forcingMatrix,
                          const int M,
                          const int N);

  void makeBoundaryMatrix (Ref<MatrixXd> boundaryMatrix,
                           const int M,
                           const int N,
                           const double h,
                           const double * viscosity);

  void makeBCLaplacianXBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,
                              const double * viscosity);

  void makeBCLaplacianYBlock (Ref<MatrixXd> laplacianBC,
                              const int M,
                              const int N,
                              const double h,                                 
                              const double * viscosity);

  void makeBCDivXBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h);

  void makeBCDivYBlock (Ref<MatrixXd> divBC,
                        const int M,
                        const int N,
                        const double h);
}
