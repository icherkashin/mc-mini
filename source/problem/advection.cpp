#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "debug.h"

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
  double * temperatureData = geometry.getTemperatureData();
  double * temperatureBoundaryData = geometry.getTemperatureBoundaryData();
  
  double * uVelocityData = geometry.getUVelocityData();
  double * vVelocityData = geometry.getVVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  // Initilize half-time data
  static double * halfTimeTemperatureData = new double[M * N];
  static double * halfTimeUOffsetTemperatureData = new double[M * (N - 1)];
  static double * halfTimeVOffsetTemperatureData = new double[(M - 1) * N];
  static double * halfTimeForcingData = new double[2 * M * N - M - N];
  static double * halfTimeStokesSolnData = new double[3 * M * N - M - N];
  // Cell-centered velocities at time t^n
  static double * cellCenteredVelocityData = new double[2 * M * N];
  static double * cellCenteredUVelocity = cellCenteredVelocityData;
  static double * cellCenteredVVelocity = cellCenteredVelocityData + M * N;

  Map<VectorXd> halfTimeStokesSolnVector (halfTimeStokesSolnData, 3 * M * N - M - N);
  
  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;
  double leftNeighborT, rightNeighborT, bottomNeighborT, topNeighborT;
  
  // Calculate cell-centered velocities
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
      
      cellCenteredVVelocity[i * N + j] = 0.5 * (topVelocity + bottomVelocity);
      cellCenteredUVelocity[i * N + j] = 0.5 * (leftVelocity + rightVelocity);
    }
  }
  
  // Calculate temperatures at half-time level
  // Loop over U-velocity positions
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < (N - 1); ++j) {
      halfTimeUOffsetTemperatureData[i * (N - 1) + j] = 0;
      // Flag to show which Riemann Problem solution to use
      int riemannFlag;
      if (cellCenteredUVelocity[i * N + j] > 0 && cellCenteredUVelocity[i * N + (j + 1)] > 0) {
        riemannFlag = 0b01;
      } else if (cellCenteredUVelocity[i * N + j] < 0 && cellCenteredUVelocity[i * N + (j + 1)] < 0) {
        riemannFlag = 0b10;
      } else {
        riemannFlag = 0b11;
      }

      if (riemannFlag & 0b01) {
        if (j == 0) {
          leftNeighborT = temperatureData[i * N + j];
        } else {
          leftNeighborT = temperatureData[i * N + (j - 1)];
        }
        rightNeighborT = temperatureData[i * N + (j + 1)];

        halfTimeUOffsetTemperatureData[i * (N - 1) + j] += 
            temperatureData[i * N + j] + 
              (h / 2 - deltaT / 2 * cellCenteredUVelocity[i * N + j]) * (rightNeighborT - leftNeighborT) / (2 * h);
      }
      if (riemannFlag & 0b10) {
        leftNeighborT = temperatureData[i * N + j];
        if (j == (N - 2)) {
          rightNeighborT = temperatureData[i * N + (j + 1)];
        } else {
          rightNeighborT = temperatureData[i * N + (j + 2)];
        }

        halfTimeUOffsetTemperatureData[i * (N - 1) + j] += 
            temperatureData[i * N + (j + 1)] -
              (h / 2 - deltaT / 2 * cellCenteredUVelocity[i * N + (j + 1)]) * (rightNeighborT - leftNeighborT) / (2 * h); 
      
      }
      if (riemannFlag == 0b11) {
        halfTimeUOffsetTemperatureData[i * (N - 1) + j] /= 2;
      }

      if (j == 0) {
        leftNeighborT = temperatureData[i * N + j];
      } else {
        leftNeighborT = temperatureData[i * N + (j - 1)];
      }
      rightNeighborT = temperatureData[i * N + (j + 1)]; 
      if (i == 0) {
        bottomNeighborT = temperatureData[i * N + j];
      } else {
        bottomNeighborT = temperatureData[(i - 1) * N + j];
      }
      if (i == (M - 1)) {
        topNeighborT = temperatureData[i * N + j];
      } else {
        topNeighborT = temperatureData[(i + 1) * N + j];
      }
    }
  }

  // Loop over V-velocity positions
  for (int i = 0; i < (M - 1); ++i) {
    for (int j = 0; j < N; ++j) {
      halfTimeVOffsetTemperatureData[i * N + j] = 0;
      // Flag to show which Riemann Problem solution to use
      int riemannFlag;
      if (cellCenteredVVelocity[i * N + j] > 0 && cellCenteredVVelocity[i * (N + 1) + j] > 0) {
        riemannFlag = 0b01;
      } else if (cellCenteredVVelocity[i * N + j] < 0 && cellCenteredVVelocity[i * (N + 1) + j] > 0) {
        riemannFlag = 0b10;
      } else {
        riemannFlag = 0b11;
      }

      // Calculate T_bottom
      if (riemannFlag & 0b01) {
        if (i == 0) {
          bottomNeighborT = temperatureData[i * N + j];
        } else {
          bottomNeighborT = temperatureData[(i - 1) * N + j];
        }
        topNeighborT = temperatureData[(i + 1) * N + j];

        halfTimeVOffsetTemperatureData[i * N + j] +=
            temperatureData[i * N + j] +
              (h / 2 - deltaT / 2 * cellCenteredVVelocity[i * N + j]) * (topNeighborT - bottomNeighborT) / (2 * h);
      }
      
      // Calculate T_top
      if (riemannFlag & 0b10) {
        bottomNeighborT = temperatureData[i * N + j];
        if (i == (N - 1)) {
          topNeighborT = temperatureData[i * N + j];
        } else {
          topNeighborT = temperatureData[(i + 1) * N + j];
        }

        halfTimeVOffsetTemperatureData[i * N + j] += 
            temperatureData[(i + 1) * N + j] -
              (h / 2 - deltaT / 2 * cellCenteredVVelocity[(i + 1) * N + j]) * (topNeighborT - bottomNeighborT) / (2 * h);

      }
      // Average T_bottom and T_top
      if (riemannFlag == 0b11) {
        halfTimeVOffsetTemperatureData[i * N + j] /= 2;
      }

      if (j == 0) {
        leftNeighborT = temperatureData[i * N + j];
      } else {
        leftNeighborT = temperatureData[i * N + (j - 1)];
      }
      if (j == (N - 1)) {
        rightNeighborT = temperatureData[i * N + j];
      } else {
        rightNeighborT = temperatureData[i * N + (j + 1)];
      }
    }
  }
  
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = temperatureData[i * N + j]; 
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureData[i * (N - 1) + (j - 1)]; 
      }
      if (j == (N - 1)) {
        rightNeighborT = temperatureData[i * N + j]; 
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureData[i * (N - 1) + j];
      }

      if (i == 0) {
        bottomNeighborT = temperatureBoundaryData[j];
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureData[(i - 1) * N + j];
      }

      if (i == (M - 1)) {
        topNeighborT = temperatureBoundaryData[N + j];
      } else {
        topNeighborT = halfTimeVOffsetTemperatureData[i * N + j];
      }

      halfTimeTemperatureData[i * N + j] = (rightNeighborT + leftNeighborT + topNeighborT + bottomNeighborT) / 4;
    }
  }
  
  if (DEBUG) {
    std::cout << "<Half-Time Temperature Data>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (halfTimeTemperatureData, M, N) << std::endl << std::endl;
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
  if (DEBUG) {
    std::cout << "<Half-Time U Forcing Data>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (halfTimeUForcingData, M, N - 1) << std::endl << std::endl;

    std::cout << "<Half-Time V Forcing Data>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (halfTimeVForcingData, M - 1, N) << std::endl << std::endl;
  }

  // Initialize half-time solver
  static SparseMatrix<double> stokesMatrix   (3 * M * N - M - N, 3 * M * N - M - N);
  static SparseMatrix<double> forcingMatrix  (3 * M * N - M - N, 2 * M * N - M - N);
  static SparseMatrix<double> boundaryMatrix (3 * M * N - M - N, 2 * M + 2 * N);
  static SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  static bool initialized;

  if (!(initialized) || !(viscosityModel=="constant")) {
    
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

    #ifdef DEBUG
      cout << "<Viscosity Data>" << endl;
      cout << Map<Matrix <double, Dynamic, Dynamic, RowMajor> > (viscosityData, M, N) << endl << endl;
    #endif
  }

  if (DEBUG) {
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (geometry.getVelocityBoundaryData(), M, 2) << std::endl <<std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (geometry.getVelocityBoundaryData() + 2 * M, 2, N) << std::endl << std::endl;
  }

  // Solve stokes at the half-time to find velocities
  halfTimeStokesSolnVector = solver.solve 
           (forcingMatrix  * Map<VectorXd>(halfTimeForcingData, 2 * M * N - M - N) + 
            boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));

  double * halfTimeUVelocity = halfTimeStokesSolnData;
  double * halfTimeVVelocity = halfTimeStokesSolnData + M * (N - 1);
  for (int i = 0; i < M; i++)
    for (int j = 0; j < (N - 1); j++)
      if (abs(halfTimeUVelocity[i * (N - 1) + j]) < 10E-10)
        halfTimeUVelocity[i * (N - 1) + j] = 0;

  for (int i = 0; i < (M - 1); i++) 
    for (int j = 0; j < N; ++j)
      if (abs (halfTimeVVelocity[i * N + j]) < 10E-10)
        halfTimeVVelocity[i * N + j] = 0;

  if (DEBUG) {
    std::cout << "<Half-Time U Velocity>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (halfTimeUVelocity, M, N - 1) << std::endl << std::endl;

    std::cout << "<Half-Time V Velocity>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (halfTimeVVelocity, M - 1, N) << std::endl << std::endl;
  }

  // Solve for full-time temperature
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = 0; //halfTimeUOffsetTemperatureData[i * (N - 1) + j];
        leftVelocity = 0;
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureData[i * (N - 1) + (j - 1)]; 
        leftVelocity = halfTimeUVelocity[i * (N - 1) + (j - 1)];
      }
      if (j == (N - 1)) {
        rightNeighborT = 0; //halfTimeUOffsetTemperatureData[i * (N - 1) + (j - 1)];
        rightVelocity = 0;
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureData[i * (N - 1) + j];
        rightVelocity = halfTimeUVelocity[i * (N - 1) + j];
      }

      if (i == 0) {
        bottomNeighborT = 0; // halfTimeVOffsetTemperatureData[i * N + j];
        bottomVelocity = 0;
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureData[(i - 1) * N + j];
        bottomVelocity = halfTimeVVelocity[(i - 1) * N + j];
      }

      if (i == (M - 1)) {
        topNeighborT = 0; // halfTimeVOffsetTemperatureData[i * N + j];
        topVelocity = 0;
      } else {
        topNeighborT = halfTimeVOffsetTemperatureData[i * N + j];
        topVelocity = halfTimeVVelocity[i * N + j];
      }
      temperatureData[i * N + j] += deltaT / h * (leftVelocity * leftNeighborT - rightVelocity * rightNeighborT) + deltaT / h * (bottomVelocity * bottomNeighborT - topVelocity * topNeighborT);

      if (std::isnan((double)temperatureData[i * N + j])) {
        if (DEBUG) {
          std::cout << "<NaN at  = " << i << ", j = " << j << ">" << std::endl;
          std::cout << "<Neighbor values were: " << std::endl;
          std::cout << "\tleftNeighborT   = " << leftNeighborT << std::endl;
          std::cout << "\trightNeighborT  = " << rightNeighborT << std::endl;
          std::cout << "\tbottomNeighborT = " << bottomNeighborT << std::endl;
          std::cout << "\ttopNeighborT    = " << topNeighborT << ">" << std::endl;
          std::cout << "<Velocity values were: " << std::endl;
          std::cout << "\tleftVelocity    = " << leftVelocity << std::endl;
          std::cout << "\trightVelocity   = " << rightVelocity << std::endl;
          std::cout << "\tbottomVelocity  = " << bottomVelocity << std::endl;
          std::cout << "\ttopVelocity     = " << topVelocity << ">" << std::endl;
          std::cout << "<Shutting down now!>" << std::endl;
        }
        throw "found NaN";
      }
    }
  }

  if (DEBUG) {
    std::cout << "<Full-Time Temperature Data>" << std::endl;
    std::cout << Map<Matrix<double, Dynamic, Dynamic, RowMajor> > (temperatureData, M, N) << std::endl << std::endl;
  }
}
