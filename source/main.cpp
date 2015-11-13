// Necessary for output in the command line.
#include <iostream>
// Necessary for throwing exceptions.
#include <stdexcept>

// Functions and data structures related to the geometry of the problem.
#include "geometry/geometry.h"
// Functions and data structures related to the implementation of numerical methods.
#include "problem/problem.h"
// Functions and data structures related to the output of the solution data.
#include "output/output.h"
// Functions and data structures related to the parser of parameter files.
#include "parser/parser.h"

int main(int argc, char ** argv) {

// The valid command line usage is "./mc-mini <parameter file>". Otherwise, throw an exception.
  if (argc == 1) {
    throw std::invalid_argument("usage: " + std::string{argv[0]} + " <parameter file>.");
  }

// Initialize the parser with the specified parameter file.
  ParamParser parser(std::string{argv[1]});
// Initialize geometry parameters.
  GeometryStructure geometry (parser);
// Initialize parameters related to the implementation of the employed numerical method.
  ProblemStructure  problem  (parser, geometry);
// Initialize parameters related to output structure.
  OutputStructure   output   (parser, geometry, problem);

// Initialize the initial data for the problem to be solved.
  problem.initializeProblem();

// Main loop where computations are made and data is output for each timestep of the problem.
  do
  {
    // 1. Solve Stokes equations.
    problem.solveStokes();
    // 2. Initialize the right hand side (forcing terms).
    problem.updateForcingTerms();
    // 3. Recalculate time step.
    problem.recalculateTimestep();
    // 4. Output the solution data.
    output.outputData (problem.getTimestepNumber());
    // 5. Solve advection-diffusion equation.
    problem.solveAdvectionDiffusion();
    // 6. Output which time step is being computed.
    std::cout << "<Timestep: " << problem.getTimestepNumber() << "; t=" << problem.getTime() << ">" << std::endl;
  }
  while (problem.advanceTimestep()); // Loop termination criterion: problem.getTimestepNumber() = end_timestep.

  problem.solveStokes(); // Solve Stokes equations.
  output.outputData (problem.getTimestepNumber()); // Output the solution data.

  return 0;
}