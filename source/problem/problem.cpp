#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"

using namespace Eigen;
using namespace std;

ProblemStructure::ProblemStructure 
    (ParamParser&       pp,
     GeometryStructure& gs) :
    parser (pp),
    geometry (gs) {
  M = geometry.getM();
  N = geometry.getN();

  if (parser.push ("problemParams")) {
    parser.getParamDouble   ("cfl",              cfl);

    parser.getParamDouble   ("startTime",        time);
    parser.queryParamDouble ("endTime",          endTime, INT_MAX);
    parser.queryParamInt    ("endStep",          endStep, INT_MAX);
    timestepNumber = 0;

    parser.queryParamDouble ("xExtent",          xExtent, 0.0);
    parser.queryParamDouble ("yExtent",          yExtent, 0.0);
    
    // Exactly one of xExtent or yExtent must be set. The other will be set
    // depending upon the first.
    assert (((xExtent == 0) ^ (yExtent == 0)));
    
    if (xExtent == 0) {
      h      = yExtent / double(M);
      xExtent = h * double(N);
    } else {
      h      = xExtent / double(N);
      yExtent = h * double(M);
    }

    parser.getParamDouble   ("diffusivity",      diffusivity);

    parser.queryParamString ("forcingModel",     forcingModel,     "tauBenchmark");
    parser.queryParamString ("temperatureModel", temperatureModel, "constant");
    parser.queryParamString ("viscosityModel",   viscosityModel,   "constant");
    parser.queryParamString ("boundaryModel",    boundaryModel,    "tauBenchmark");

    parser.queryParamString ("outputFile",       outputFile,       "output.h5");

    parser.pop();
  }
}

ProblemStructure::~ProblemStructure() {
}                                   


bool ProblemStructure::advanceTimestep() {
  time += deltaT;
  timestepNumber++;

  return ((time < endTime) && (timestepNumber < endStep));
}

void ProblemStructure::recalculateTimestep() {
  Map<VectorXd> uVelocityVector (geometry.getUVelocityData(), M *       (N - 1));
  Map<VectorXd> vVelocityVector (geometry.getVVelocityData(), (M - 1) * N);
  Map<VectorXd> uVelocityBoundaryVector (geometry.getUVelocityBoundaryData(), 2 * N);
  Map<VectorXd> vVelocityBoundaryVector (geometry.getVVelocityBoundaryData(), 2 * M);

  double maxInteriorUVelocity = uVelocityVector.maxCoeff();
  double maxBoundaryUVelocity = uVelocityBoundaryVector.maxCoeff();
  double maxUVelocity = (maxInteriorUVelocity > maxBoundaryUVelocity) ? maxInteriorUVelocity : maxBoundaryUVelocity;

  double maxInteriorVVelocity = vVelocityVector.maxCoeff();
  double maxBoundaryVVelocity = vVelocityBoundaryVector.maxCoeff();
  double maxVVelocity = (maxInteriorVVelocity > maxBoundaryVVelocity) ? maxInteriorVVelocity : maxBoundaryVVelocity;
  double velocityDeltaT = cfl * h / sqrt (maxUVelocity * maxUVelocity + maxVVelocity * maxVVelocity);
  if (maxUVelocity == 0 || maxVVelocity == 0)
    velocityDeltaT = INT_MAX;
  double diffusionDeltaT = cfl * h / diffusivity;

  deltaT = (velocityDeltaT < diffusionDeltaT) ? velocityDeltaT : diffusionDeltaT;
}

double ProblemStructure::getH() {
  return h;
}

double ProblemStructure::getTime() {
  return time;
}

int ProblemStructure::getTimestepNumber() {
  return timestepNumber;
}

double ProblemStructure::getEndTime() {
  return endTime;
}
