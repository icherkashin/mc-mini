#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;

// upwind method. Stable but inefficient.
void ProblemStructure::upwindMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  double * uVelocityData = geometry.getUVelocityData();
  double * vVelocityData = geometry.getVVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;
  double leftFlux, rightFlux, bottomFlux, topFlux;

  VectorXd temporaryVector = temperatureVector;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftVelocity   = (j == 0)       ? 
                        uVelocityBoundaryData[2 * i]     : 
                        uVelocityData[i * (N - 1) + (j - 1)];
      rightVelocity  = (j == (N - 1)) ? 
                        uVelocityBoundaryData[2 * i + 1] : 
                        uVelocityData[i * (N - 1) + j];
      bottomVelocity = (i == 0)       ?
                        vVelocityBoundaryData[j]         :
                        vVelocityData[(i - 1) * N + j];
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryData[M + j]     :
                        vVelocityData[i * N + j];
      

      leftFlux = rightFlux = 0;
      topFlux = bottomFlux = 0;

      if (leftVelocity < 0) {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + j) * leftVelocity * deltaT / h;
        }
      } else {
        if (j > 0) {
          leftFlux = temporaryVector(i * N + (j - 1)) * leftVelocity * deltaT / h;
        }
      }

      if (rightVelocity > 0) {
        if (j < (N - 1)) {
          rightFlux = temporaryVector(i * N + j) * rightVelocity * deltaT / h;
        }
      } else {
        if (j < (N - 1)) {
          rightFlux = temporaryVector(i * N + (j + 1)) * rightVelocity * deltaT / h;
        }
      }

      if (bottomVelocity < 0) {
        if (i > 0) {
          bottomFlux = temporaryVector(i * N + j) * bottomVelocity * deltaT / h;
        }
      } else {
        if (i > 0) {
          bottomFlux = temporaryVector((i - 1) * N + j) * bottomVelocity * deltaT / h;
        }
      }

      if (topVelocity > 0) {
        if (i < (M - 1)) {
          topFlux = temporaryVector (i * N + j) * topVelocity * deltaT / h;
        }
      } else {
        if (i < (M - 1)) {
          topFlux = temporaryVector ((i + 1) * N + j) * topVelocity * deltaT / h;
        }
      }

      temperatureVector(i * N + j) = temporaryVector(i * N + j) + leftFlux - rightFlux + bottomFlux - topFlux;
    }
  }
}
