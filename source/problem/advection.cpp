#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "debug.h"

using namespace Eigen;
using namespace std;

// upwind method. Stable but inefficient.
void ProblemStructure::upwindMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;
  double leftFlux, rightFlux, bottomFlux, topFlux;

  VectorXd temporaryVector = temperatureVector;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j <  N; ++j) {
      // Find all four edge velocities for the current cell. 
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryWindow (0, i) :
                        uVelocityWindow (j - 1, i);
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryWindow (1, i) :
                        uVelocityWindow (j, i);
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryWindow (j, 0) :
                        vVelocityWindow (j, i - 1);
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryWindow (j, 1) :
                        vVelocityWindow (j, i);

      leftFlux = rightFlux = 0;
      topFlux = bottomFlux = 0;
      
      // Solve the Riemann problem on the neighboring velocities and calculate
      // the fluxes accross each edge.
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
  // Temperature data (MxN cell-centered grid)
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), N, M);
  // Temperature boundary data (2xN transverse boundary grid)
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), N, 2);
  
  // U Velocity Data (Mx(N-1) lateral offset grid)
  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  // V Velocity Data ((M-1)xN transverse offset grid)
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  
  // U Velocity Boundary Data (Mx2 lateral boundary grid)
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  // V Velocity Boundary Data (2xN transverse boundary grid)
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);

  // Initilize half-time data
  // Half-time temperature data (MxN cell-centered grid)
  static DataWindow<double> halfTimeTemperatureWindow (new double[M * N], N, M);
  // Half-time U-Offset temperature data (Mx(N-1) lateral offset grid)
  static DataWindow<double> halfTimeUOffsetTemperatureWindow (new double[M * (N - 1)], N - 1, M);
  // Half-time V-offset temperature data ((M-1)xN transverse offset grid)
  static DataWindow<double> halfTimeVOffsetTemperatureWindow (new double[(M - 1) * N], N, M - 1);

  // Half-time Forcing Data (for use in the Stokes solve)
  static double * halfTimeForcingData = new double[2 * M * N - M - N];
  // Half-time U forcing data (Mx(N-1) lateral offset grid)
  static double * halfTimeUForcingData = halfTimeForcingData;
  static DataWindow<double> halfTimeUForcingWindow (halfTimeUForcingData, N - 1, M); 
  // Half-time V forcing data ((M-1)xN transverse offset grid)
  static double * halfTimeVForcingData = halfTimeForcingData + M * (N - 1);
  static DataWindow<double> halfTimeVForcingWindow (halfTimeVForcingData, N, M - 1);
  
  // Half-time Stokes solution data (for use in the Stokes solve)
  static double * halfTimeStokesSolnData = new double[3 * M * N - M - N];
  // Half-time U velocity data (Mx(N-1) lateral offset grid)
  static double * halfTimeUVelocityData = halfTimeStokesSolnData;
  static DataWindow<double> halfTimeUVelocityWindow (halfTimeUVelocityData, N - 1, M);
  // Half-time V velocity data ((M-1)xN transverse offset grid)
  static double * halfTimeVVelocityData = halfTimeStokesSolnData + M * (N - 1);
  static DataWindow<double> halfTimeVVelocityWindow (halfTimeVVelocityData, N, M - 1);
  // Half-time pressure data (MxN cell-centered grid) (not technically necessary)
  static double * halfTimePressureData = halfTimeStokesSolnData + 2 * M * N - M - N;
  static DataWindow<double> halfTimePressureWindow (halfTimePressureData, N, M);

  // Cell-centered velocities
  static double * cellCenteredVelocityData = new double[2 * M * N];
  // Cell-centered U velocity (MxN cell-centered grid)
  static double * cellCenteredUVelocityData = cellCenteredVelocityData;
  static DataWindow<double> cellCenteredUVelocityWindow (cellCenteredUVelocityData, N, M);
  // Cell-centered V velocity (MxN cell-centered grid)
  static double * cellCenteredVVelocityData = cellCenteredVelocityData + M * N;
  static DataWindow<double> cellCenteredVVelocityWindow (cellCenteredVVelocityData, N, M);

  Map<VectorXd> halfTimeStokesSolnVector (halfTimeStokesSolnData, 3 * M * N - M - N);
  
  double leftNeighborT, rightNeighborT, bottomNeighborT, topNeighborT;
  
  // Calculate cell-centered velocities (MxN cell-centered grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double leftVelocity   = (j == 0) ?
                              uVelocityBoundaryWindow (0, i) :
                              uVelocityWindow (j - 1, i);
      double rightVelocity  = (j == (N - 1)) ?
                              uVelocityBoundaryWindow (1, i) :
                              uVelocityWindow (j, i);
      double bottomVelocity = (i == 0) ?
                              vVelocityBoundaryWindow (j, 0) :
                              vVelocityWindow (j, i - 1);
      double topVelocity    = (i == (M - 1)) ?
                              vVelocityBoundaryWindow (j, 1) :
                              vVelocityWindow (j, i);
      
      cellCenteredVVelocityWindow (j, i) = (topVelocity + bottomVelocity) / 2;
      cellCenteredUVelocityWindow (j, i) = (leftVelocity + rightVelocity) / 2;
      // Check if the cell-centered velocity is machine-epsilon, and make it true zero if so.
      if (abs (cellCenteredUVelocityWindow (j, i)) < 1E-10) 
        cellCenteredUVelocityWindow (j, i) = 0;
      if (abs (cellCenteredVVelocityWindow (j, i)) < 1E-10)
        cellCenteredVVelocityWindow (j, i) = 0;
    }
  }

  #ifdef DEBUG
    cout << "<Half-time Cell-Centered Velocities Calculated>" << std::endl;
    cout << "<Cell-Centered U Velocity>" << endl;
    cout << cellCenteredUVelocityWindow.displayMatrix() << endl;
    cout << "<Cell-Centered V Velocity>" << endl;
    cout << cellCenteredVVelocityWindow.displayMatrix() << endl;
  #endif
  
  // Calculate temperatures at half-time level
  // Calculate half-time U-offset temperatures (Mx(N-1) lateral offset grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < (N - 1); ++j) {
      halfTimeUOffsetTemperatureWindow (j, i) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemmanFlag since the case where velocities do not share the 
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;
      if (cellCenteredUVelocityWindow (j,     i) >= 0 && 
          cellCenteredUVelocityWindow (j + 1, i) >= 0) {
        riemannFlag = 0b01;
      } else if (cellCenteredUVelocityWindow (j,     i) <= 0 &&
                 cellCenteredUVelocityWindow (j + 1, i) <= 0) {
        riemannFlag = 0b10;
      } else {
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the left
      // neighboring cell.
      if (riemannFlag & 0b01) {
        if (j == 0) {
          leftNeighborT = temperatureWindow (j, i);
        } else {
          leftNeighborT = temperatureWindow (j - 1, i);
        }
        rightNeighborT = temperatureWindow (j + 1, i);

        halfTimeUOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j, i) +
              (h / 2 - deltaT / 2 * cellCenteredUVelocityWindow (j, i)) * 
                (rightNeighborT - leftNeighborT) / (2 * h);
      }
      
      // One or both velocities are negative. Take the flux from the right
      // neighboring cell.
      if (riemannFlag & 0b10) {
        leftNeighborT = temperatureWindow (j, i);
        if (j == (N - 2)) {
          rightNeighborT = temperatureWindow (j + 1, i);
        } else {
          rightNeighborT = temperatureWindow (j + 2, i);
        }

        halfTimeUOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j + 1, i) -
              (h / 2 + deltaT / 2 * cellCenteredUVelocityWindow (j + 1, i)) * 
                (rightNeighborT - leftNeighborT) / (2 * h); 
      
      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells
      if (riemannFlag == 0b11) halfTimeUOffsetTemperatureWindow (j, i) /= 2;
    }
  }

  #ifdef DEBUG
    cout << "<Half-Time Offset Temperatures>" << endl;
    cout << "<Half-Time U-Offset Temperature>" << endl;
    cout << halfTimeUOffsetTemperatureWindow.displayMatrix() << endl;
  #endif

  // Calculate half-time V-offset temperatures ((M-1)xN transverse offset grid)
  for (int i = 0; i < (M - 1); ++i) {
    for (int j = 0; j < N; ++j) {
      halfTimeVOffsetTemperatureWindow (j, i) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemannFlag since the case where velocities do not share the
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;

      if (cellCenteredVVelocityWindow (j, i    ) >= 0 && 
          cellCenteredVVelocityWindow (j, i + 1) >= 0) {
        // One or both velocities are positive
        riemannFlag = 0b01;
      } else if (cellCenteredVVelocityWindow (j, i    ) <= 0 && 
                 cellCenteredVVelocityWindow (j, i + 1) <= 0) {
        // One or both velocities are negative
        riemannFlag = 0b10;
      } else {
        // One velocity is positive, one is negative.
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the bottom
      // neighboring cell
      if (riemannFlag & 0b01) {
        bottomNeighborT = temperatureWindow (j, i + 1);
        if (i == 0) {
          topNeighborT = temperatureWindow (j, i);
        } else {
          topNeighborT = temperatureWindow (j, i - 1);
        }

        halfTimeVOffsetTemperatureWindow (j, i) +=
            temperatureWindow (j, i) -
              (h / 2 + deltaT / 2 * cellCenteredVVelocityWindow (j, i)) * 
               (topNeighborT - bottomNeighborT) / (2 * h);
      }
      
      // One or both velocities are negative. Take the flux from the top
      // neighboring cell
      if (riemannFlag & 0b10) {
        if (i == (M - 2)) {
          bottomNeighborT = temperatureWindow (j, i + 1);
        } else {
          bottomNeighborT = temperatureWindow (j, i + 2);
        }
        topNeighborT = temperatureWindow (j, i);

        halfTimeVOffsetTemperatureWindow (j, i) += 
            temperatureWindow (j, i + 1) +
              (h / 2 - deltaT / 2 * cellCenteredVVelocityWindow (j, i + 1)) * 
               (topNeighborT - bottomNeighborT) / (2 * h);
      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells.
      if (riemannFlag == 0b11) {
        halfTimeVOffsetTemperatureWindow (j, i) /= 2;
      }
    }
  }

  #ifdef DEBUG
    cout << "<Half-Time V-Offset Temperature>" << endl;
    cout << halfTimeVOffsetTemperatureWindow.displayMatrix() << endl << endl;
  #endif

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = temperatureWindow (j, i); 
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureWindow (j - 1, i); 
      }
      if (j == (N - 1)) {
        rightNeighborT = temperatureWindow (j, i); 
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureWindow (j, i);
      }

      if (i == 0) {
        bottomNeighborT = temperatureBoundaryWindow (j, 0);
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureWindow (j, i - 1);
      }

      if (i == (M - 1)) {
        topNeighborT = temperatureBoundaryWindow (j, 1);
      } else {
        topNeighborT = halfTimeVOffsetTemperatureWindow (j, i);
      }

      halfTimeTemperatureWindow (j, i) = (rightNeighborT + leftNeighborT + topNeighborT + bottomNeighborT) / 4;
    }
  }
  
  #ifdef DEBUG
    cout << "<Half-Time Temperature Data>" << endl;
    cout << halfTimeTemperatureWindow.displayMatrix() << endl << endl;
  #endif
  
  // Calculate half-time forcing
  if (forcingModel == "tauBenchmark") {
    // Benchmark taken from Tau (1991; JCP Vol. 99)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (j, i) = 3 * cos ((j + 1) * h) * sin ((i + 0.5) * h);

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (j, i) = -sin ((j + 0.5) * h) * cos ((i + 1) * h);

  } else if (forcingModel == "solCXBenchmark" ||
             forcingModel == "solKZBenchmark") {
    // solCX Benchmark taken from Kronbichler et al. (2011)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (j, i) = 0;

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (j, i) = - sin((i + 0.5) * M_PI * h) * cos ((j + 1) * M_PI * h);

  } else if (forcingModel == "buoyancy") {
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

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < (N - 1); ++j)
        halfTimeUForcingWindow (j, i) = 0;

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j) {
        halfTimeVForcingWindow (j, i) =  -1 * densityConstant *
                                          (1 - thermalExpansion * 
                                           ((halfTimeTemperatureWindow (j, i) + 
                                             halfTimeTemperatureWindow (j, i + 1)) / 2 -
                                            referenceTemperature));
      }
  } else {
    cout << "<Unexpected forcing model: \"" << forcingModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Half-Time Forcing Data>" << endl;
    cout << "<Half-Time U Forcing Data>" << endl;
    cout << halfTimeUForcingWindow.displayMatrix() << endl;
    cout << "<Half-Time V Forcing Data>" << endl;
    cout << halfTimeVForcingWindow.displayMatrix() << endl << endl;
  #endif

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
      cout << DataWindow<double> (geometry.getViscosityData(), N + 1, M + 1).displayMatrix() << endl << endl;
    #endif
  }

  #ifdef DEBUG
    cout << "<Velocity Boundary Data>" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    cout << uVelocityBoundaryWindow.displayMatrix() << endl;
    cout << "<V Velocity Boundary Data>" << endl;
    cout << vVelocityBoundaryWindow.displayMatrix() << endl << endl;
  #endif

  // Solve stokes at the half-time to find velocities
  halfTimeStokesSolnVector = solver.solve 
           (forcingMatrix  * Map<VectorXd>(halfTimeForcingData, 2 * M * N - M - N) + 
            boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));

  for (int i = 0; i < M; i++)
    for (int j = 0; j < (N - 1); j++)
      if (abs(halfTimeUVelocityWindow (j, i)) < 10E-10)
        halfTimeUVelocityWindow (j, i) = 0;

  for (int i = 0; i < (M - 1); i++) 
    for (int j = 0; j < N; ++j)
      if (abs (halfTimeVVelocityWindow (j, i)) < 10E-10)
        halfTimeVVelocityWindow (j, i) = 0;

  #ifdef DEBUG
    cout << "<Half-Time Velocities>" << endl;
    cout << "<Half-Time U Velocity>" << endl;
    cout << halfTimeUVelocityWindow.displayMatrix() << endl;
    cout << "<Half-Time V Velocity>" << std::endl;
    cout << halfTimeVVelocityWindow.displayMatrix() << endl << endl;
  #endif

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;

  // Solve for full-time temperature
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = 0; 
        leftVelocity = 0;
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureWindow (j - 1, i);
        leftVelocity = halfTimeUVelocityWindow (j - 1, i);
      }
      if (j == (N - 1)) {
        rightNeighborT = 0;
        rightVelocity = 0;
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureWindow (j, i);
        rightVelocity = halfTimeUVelocityWindow (j, i);
      }

      if (i == 0) {
        bottomNeighborT = 0; 
        bottomVelocity = 0;
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureWindow (j, i - 1);
        bottomVelocity = halfTimeVVelocityWindow (j, i - 1);
      }

      if (i == (M - 1)) {
        topNeighborT = 0;
        topVelocity = 0;
      } else {
        topNeighborT = halfTimeVOffsetTemperatureWindow (j, i);
        topVelocity = halfTimeVVelocityWindow (j, i);
      }
      temperatureWindow (j, i) += deltaT / h * (leftVelocity * leftNeighborT - rightVelocity * rightNeighborT) + deltaT / h * (bottomVelocity * bottomNeighborT - topVelocity * topNeighborT);

      if (std::isnan((double)temperatureWindow (j, i))) {
        #ifdef DEBUG
          cout << "<NaN at " << i << "," << j << ">" << endl;
          cout << "<Neighbor values were: " << endl;
          cout << "\tleftNeighborT   = " << leftNeighborT << endl;
          cout << "\trightNeighborT  = " << rightNeighborT << endl;
          cout << "\tbottomNeighborT = " << bottomNeighborT << endl;
          cout << "\ttopNeighborT    = " << topNeighborT << ">" << endl;
          cout << "<Velocity values were: " << endl;
          cout << "\tleftVelocity    = " << leftVelocity << endl;
          cout << "\trightVelocity   = " << rightVelocity << endl;
          cout << "\tbottomVelocity  = " << bottomVelocity << endl;
          cout << "\ttopVelocity     = " << topVelocity << ">" << endl;
          cout << "<Shutting down now!>" << endl;
        #endif
        throw "found NaN";
      }
    }
  }

  #ifdef DEBUG
    cout << "<Full-Time Temperature Data>" << endl;
    cout << temperatureWindow.displayMatrix() << endl << endl;
  #endif
}
