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
    for (int j = 0; j <  N; ++j) {
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryData [2 * i] :
                        uVelocityData [i * (N - 1) + (j - 1)];
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryData [2 * i + 1] :
                        uVelocityData [i * (N - 1) + j];
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryData [j] :
                        vVelocityData [(i - 1) * N + j];
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryData [N + j] :
                        vVelocityData [i * N + j];

      leftFlux = rightFlux = 0;
      topFlux = bottomFlux = 0;

      if (j > 0) {
        if (leftVelocity < 0) {
          leftFlux = temporaryVector (i * N + j) * leftVelocity * deltaT / h;
        } else {
          leftFlux = temporaryVector (i * N + (j - 1)) * leftVelocity * deltaT / h;
        }
      }

      if (j < (N - 1)) {
        if (rightVelocity > 0) {
          rightFlux = temporaryVector (i * N + j) * rightVelocity * deltaT / h;
        } else {
          rightFlux = temporaryVector (i * N + (j + 1)) * rightVelocity * deltaT / h;
        }
      }

      if (i > 0) {
        if (bottomVelocity < 0) {
          bottomFlux = temporaryVector (i * N + j) * bottomVelocity * deltaT / h;
        } else {
          bottomFlux = temporaryVector ((i - 1) * N + j) * bottomVelocity * deltaT / h;
        }
      }

      if (i < (M - 1)) {
        if (topVelocity > 0) {
          topFlux = temporaryVector (i * N + j) * topVelocity * deltaT / h;
        } else {
          topFlux = temporaryVector ((i + 1) * N + j) * topVelocity * deltaT / h;
        }
      }

      temperatureVector (i * N + j) = temporaryVector (i * N + j) + 
                                        leftFlux - rightFlux +
                                        bottomFlux - topFlux;
    }
  }
}

void ProblemStructure::laxWendroff() {
}

void ProblemStructure::frommMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);
  
  double * uVelocityData = geometry.getUVelocityData();
  double * vVelocityData = geometry.getVVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  // Initilize half-time data
  static double * halfTimeTemperatureData;
  static double * halfTimeForcingData;
  static double * halfTimeStokesSolnData;

  // Initialize half-time solver

  static SparseMatrix<double> stokesMatrix   (3 * M * N - M - N, 3 * M * N - M - N);
  static SparseMatrix<double> forcingMatrix  (3 * M * N - M - N, 2 * M * N - M - N);
  static SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  static bool initialized;

  if (!(initialized) || !(viscosityModel=="constant")) {
    halfTimeTemperatureData = new double[M * N];
    halfTimeForcingData = new double[2 * M * N - M - N];
    halfTimeStokesSolnData = new double[3 * M * N - M - N];
    
    double * viscosityData = geometry.getViscosityData();

    SparseForms::makeStokesMatrix (stokesMatrix, M, N, h, viscosityData);
    stokesMatrix.makeCompressed();
    SparseForms::makeForcingMatrix (forcingMatrix, M, N);
    forcingMatrix.makeCompressed();
    SparseForms::makeBoundaryMatrix (boundaryMatrix, M, N, h, viscosityData);
    boundaryMatrix.makeCompressed();

    solver.analyzePattern (stokesMatrix);
    solver.factorize      (stokesMatrix);

    initialized = true;
  }
  
  Map<VectorXd> halfTimeStokesSolnVector (halfTimeStokesSolnData, 3 * M * N - M - N);
  
  // Calculate half-time temperature data
  double leftNeighbor, rightNeighbor, bottomNeighbor, topNeighbor;
  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryData [2 * i] :
                        uVelocityData [i * (N - 1) + (j - 1)];
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryData [2 * i + 1] :
                        uVelocityData [i * (N - 1) + j];
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryData [j] :
                        vVelocityData [(i - 1) * N + j];
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryData [N + j] :
                        vVelocityData [i * N + j];
      
      leftNeighbor   = (j == 0) ?
                        temperatureVector (i * N + j) :
                        temperatureVector (i * N + (j - 1));
      rightNeighbor  = (j == (N - 1)) ?
                        temperatureVector (i * N + j) :
                        temperatureVector (i * N + (j + 1));
      bottomNeighbor = (i == 0) ?
                        temperatureBoundaryVector (j) :
                        temperatureVector ((i - 1) * N + j);
      topNeighbor    = (i == (M - 1)) ?
                        temperatureBoundaryVector (N + j) :
                        temperatureVector ((i + 1) * N + j);
      
      halfTimeTemperatureData [i * N + j] = 
        temperatureVector (i * N + j) - deltaT / 2 * 
        ((leftVelocity + rightVelocity) * (rightNeighbor - leftNeighbor) / (2 * h) + 
         (topVelocity + bottomVelocity) * (topNeighbor - bottomNeighbor) / (2 * h) -
         diffusivity * (topNeighbor + bottomNeighbor + leftNeighbor + rightNeighbor - 4 * temperatureVector (i * N + j)));
    }
  }
  
  // Calculate half-time forcing
  double referenceTemperature;
  double densityConstant;
  double thermalExpansion;

  if (parser.push ("problemParams")) {
    if (parser.push ("buoyancyModelParams")) {
      parser.queryParamDouble ("referenceTemperature", referenceTemperature, 273.15);
      parser.queryParamDouble ("densityConstant",      densityConstant,      100.0);
      parser.queryParamDouble ("thermalExpansion",     thermalExpansion,       1.0);

      parser.pop();
    }
    parser.pop();
  }

  double * halfTimeUForcingData = halfTimeForcingData;
  double * halfTimeVForcingData = halfTimeForcingData + M * (N - 1);

  for (int i = 0; i < M; ++i)
    for (int j = 0; j < (N - 1); ++j)
      halfTimeUForcingData [i * (N - 1) + j] = 0;

  for (int i = 0; i < (M - 1); ++i)
    for (int j = 0; j < N; ++j) 
      halfTimeVForcingData [i * N + j] =-1 * densityConstant *
                                         (1 - thermalExpansion *
                                          ((halfTimeTemperatureData [i * N + j] +
                                            halfTimeTemperatureData [(i + 1) * N + j]) / 2 - 
                                           referenceTemperature));

  // Solve stokes at the half-time to find velocities
  halfTimeStokesSolnVector = solver.solve 
           (forcingMatrix  * Map<VectorXd>(halfTimeForcingData, 2 * M * N - M - N) + 
            boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));

  double * halfTimeUVelocity = halfTimeStokesSolnData;
  double * halfTimeVVelocity = halfTimeStokesSolnData + M * (N - 1);

  // Use the half-time temperature and velocity to find the second-order next timestep.
  VectorXd temporaryVector = temperatureVector;
  
  double leftHalfTimeNeighbor, rightHalfTimeNeighbor, 
         bottomHalfTimeNeighbor, topHalfTimeNeighbor;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      leftHalfTimeNeighbor   = (j == 0) ?
                                halfTimeTemperatureData [i * N + j] :
                                halfTimeTemperatureData [i * N + (j - 1)];
      rightHalfTimeNeighbor  = (j == (N - 1)) ?
                                halfTimeTemperatureData [i * N + j] :
                                halfTimeTemperatureData [i * N + (j + 1)];
      bottomHalfTimeNeighbor = (i == 0) ?
                                temperatureBoundaryVector (j) :
                                halfTimeTemperatureData [(i - 1) * N + j];
      topHalfTimeNeighbor    = (i == (M - 1)) ?
                                temperatureBoundaryVector (N + j) :
                                halfTimeTemperatureData [(i + 1) * N + j];
      
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryData [2 * i] :
                        halfTimeUVelocity [i * (N - 1) + (j - 1)];
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryData [2 * i + 1] :
                        halfTimeUVelocity [i * (N - 1) + j];
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryData [j] :
                        halfTimeVVelocity [(i - 1) * N + j];
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryData [N + j] :
                        halfTimeVVelocity [i * N + j];

      temperatureVector (i * N + j) = temporaryVector (i * N + j) + deltaT / h * 
        (leftVelocity * (leftHalfTimeNeighbor + halfTimeTemperatureData [i * N + j]) / 2 -
         rightVelocity * (rightHalfTimeNeighbor + halfTimeTemperatureData [i * N + j]) / 2 +
         bottomVelocity * (bottomHalfTimeNeighbor + halfTimeTemperatureData [i * N + j]) / 2 -
         topVelocity * (topHalfTimeNeighbor + halfTimeTemperatureData [i * N + j]) / 2);
    }
  }

}

void ProblemStructure::frommVanLeer() {
}
