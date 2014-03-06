#include <iostream>
#include <iomanip>
#include <cassert>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "matrixForms/denseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "output/output.h"
#include "parser/parser.h"


using namespace Eigen;
using namespace std;

int main(int argc, char ** argv) {

  if (argc == 1) {
    cerr << "usage: " << string{argv[0]} << " <parameter file>." << endl << endl;
    exit (-1);
  }

  setNbThreads(8);

  ParamParser parser(string{argv[1]});
  GeometryStructure geometry (parser);
  ProblemStructure  problem  (parser, geometry);
  OutputStructure   output   (parser, geometry, problem);

  problem.initializeProblem();
  
  problem.updateForcingTerms();
  problem.solveStokes();

  
  while (problem.advanceTimestep()) {
    output.writeHDF5File (problem.getTimestepNumber());
    
    problem.updateForcingTerms();
    problem.solveStokes();
    problem.recalculateTimestep();
    problem.solveAdvectionDiffusion();

  } while (problem.advanceTimestep());

  output.writeHDF5File();
  
  int M = geometry.getM();
  int N = geometry.getN();
  double h = problem.getH();

  MatrixXd analyticU;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > uSolnMatrix (geometry.getUVelocityData(), M, (N - 1));
  analyticU.resizeLike(uSolnMatrix);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N - 1; ++j)
      analyticU (i, j) = cos ((j + 1) * h) * sin ((i + 0.5) * h);

  MatrixXd analyticV;
  Map<Matrix<double, Dynamic, Dynamic, RowMajor> > vSolnMatrix (geometry.getVVelocityData(), (M - 1), N);
  analyticV.resizeLike(vSolnMatrix);
  for (int i = 0; i < M - 1; ++i)
    for (int j = 0; j < N; ++j)
      analyticV (i, j) = - sin ((j + 0.5) * h) * cos ((i + 1) * h);

  cout << sqrt((uSolnMatrix - analyticU).squaredNorm() * h * h) 
       << "\t" 
       << sqrt((vSolnMatrix - analyticV).squaredNorm() * h * h) << endl;
  
  return 0;
}
