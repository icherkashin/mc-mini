#include <iostream>

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "output/output.h"
#include "parser/parser.h"

int main(int argc, char ** argv) {

  if (argc == 1) {
    throw std::invalid_argument("usage: " + std::string{argv[0]} + " <parameter file>.");
  }

  ParamParser parser(std::string{argv[1]});
  GeometryStructure geometry (parser);
  ProblemStructure  problem  (parser, geometry);
  OutputStructure   output   (parser, geometry, problem);

  problem.initializeProblem();
  
  do 
  {      
    std::cout << "<Timestep: " << problem.getTimestepNumber() << "; t=" << problem.getTime() << ">" << std::endl;
    problem.updateForcingTerms();
    problem.solveStokes();
    problem.recalculateTimestep();
    problem.solveAdvectionDiffusion();
    output.outputData (problem.getTimestepNumber());
  } 
  while (problem.advanceTimestep());
 
  return 0;
}