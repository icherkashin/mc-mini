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

  setNbThreads(16);

  ParamParser parser(string{argv[1]});
  GeometryStructure geometry (parser);
  ProblemStructure  problem  (parser, geometry);
  OutputStructure   output   (parser, geometry, problem);

  problem.initializeProblem();
  
  problem.updateForcingTerms();
  problem.solveStokes();

  problem.recalculateTimestep();
  
  do {    
    output.writeHDF5File (problem.getTimestepNumber());
    cerr << "Timestep: " << problem.getTimestepNumber() << "; t = " << problem.getTime() << endl;
    problem.updateForcingTerms();
    problem.solveStokes();
    problem.recalculateTimestep();
    problem.solveAdvectionDiffusion();

  } while (problem.advanceTimestep());

  output.writeHDF5File (problem.getTimestepNumber());
  cerr << "Timestep: " << problem.getTimestepNumber() << "; t = " << problem.getTime() << endl;

  output.writeHDF5File();
  
  return 0;
}
