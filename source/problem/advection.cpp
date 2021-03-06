#include <iostream>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "debug.h"

#include "boost/math/constants/constants.hpp"
const double pi = boost::math::constants::pi<double>();

using namespace Eigen;
using namespace std;

// upwind method. Stable but inefficient.
void ProblemStructure::upwindMethod() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), M, N - 1);
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), M - 1, N);
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), M, 2);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), 2, N);

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;
  double leftFlux, rightFlux, bottomFlux, topFlux;

  VectorXd temporaryVector = temperatureVector;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j <  N; ++j) {
      // Find all four edge velocities for the current cell.
      leftVelocity   = (j == 0) ?
                        uVelocityBoundaryWindow (i, 0) :
                        uVelocityWindow (i, j - 1);
      rightVelocity  = (j == (N - 1)) ?
                        uVelocityBoundaryWindow (i, 1) :
                        uVelocityWindow (i, j);
      bottomVelocity = (i == 0) ?
                        vVelocityBoundaryWindow (0, j) :
                        vVelocityWindow (i - 1, j);
      topVelocity    = (i == (M - 1)) ?
                        vVelocityBoundaryWindow (1, j) :
                        vVelocityWindow (i, j);

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
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), M, N);
  // Temperature boundary data (2xN transverse boundary grid)
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), 2, N);
//
  // U Velocity Data (Mx(N-1) lateral offset grid)
  DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), M, N - 1);
  // V Velocity Data ((M-1)xN transverse offset grid)
  DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), M - 1, N);

  // U Velocity Boundary Data (Mx2 lateral boundary grid)
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), M, 2);
  // V Velocity Boundary Data (2xN transverse boundary grid)
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), 2, N);

  // Initilize half-time data
  // Half-time temperature data (MxN cell-centered grid)
  static DataWindow<double> halfTimeTemperatureWindow (new double[M * N], M, N);
  // Half-time U-Offset temperature data (Mx(N-1) lateral offset grid)
  static DataWindow<double> halfTimeUOffsetTemperatureWindow (new double[M * (N - 1)], M, N - 1);
  // Half-time V-offset temperature data ((M-1)xN transverse offset grid)
  static DataWindow<double> halfTimeVOffsetTemperatureWindow (new double[(M - 1) * N], M - 1, N);

  // Half-time Forcing Data (for use in the Stokes solve)
  static double * halfTimeForcingData = new double[2 * M * N - M - N];
  // Half-time U forcing data (Mx(N-1) lateral offset grid)
  static double * halfTimeUForcingData = halfTimeForcingData;
  static DataWindow<double> halfTimeUForcingWindow (halfTimeUForcingData, M, N - 1);
  // Half-time V forcing data ((M-1)xN transverse offset grid)
  static double * halfTimeVForcingData = halfTimeForcingData + M * (N - 1);
  static DataWindow<double> halfTimeVForcingWindow (halfTimeVForcingData, M - 1, N);

  // Half-time Stokes solution data (for use in the Stokes solve)
  static double * halfTimeStokesSolnData = new double[3 * M * N - M - N];
  // Half-time U velocity data (Mx(N-1) lateral offset grid)
  static double * halfTimeUVelocityData = halfTimeStokesSolnData;
  static DataWindow<double> halfTimeUVelocityWindow (halfTimeUVelocityData, M, N - 1);
  // Half-time V velocity data ((M-1)xN transverse offset grid)
  static double * halfTimeVVelocityData = halfTimeStokesSolnData + M * (N - 1);
  static DataWindow<double> halfTimeVVelocityWindow (halfTimeVVelocityData, M - 1, N);
  // Half-time pressure data (MxN cell-centered grid) (not technically necessary)
  static double * halfTimePressureData = halfTimeStokesSolnData + 2 * M * N - M - N;
  static DataWindow<double> halfTimePressureWindow (halfTimePressureData, M, N);

  // Cell-centered velocities
  static double * cellCenteredVelocityData = new double[2 * M * N];
  // Cell-centered U velocity (MxN cell-centered grid)
  static double * cellCenteredUVelocityData = cellCenteredVelocityData;
  static DataWindow<double> cellCenteredUVelocityWindow (cellCenteredUVelocityData, M, N);
  // Cell-centered V velocity (MxN cell-centered grid)
  static double * cellCenteredVVelocityData = cellCenteredVelocityData + M * N;
  static DataWindow<double> cellCenteredVVelocityWindow (cellCenteredVVelocityData, M, N);

  Map<VectorXd> halfTimeStokesSolnVector (halfTimeStokesSolnData, 3 * M * N - M - N);
  VectorXd temporaryTemperature = Map<VectorXd> (geometry.getTemperatureData(), N * M);

  double leftNeighborT, rightNeighborT, bottomNeighborT, topNeighborT;

  // Calculate cell-centered velocities (MxN cell-centered grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double leftVelocity   = (j == 0) ?
                              uVelocityBoundaryWindow (i, 0) :
                              uVelocityWindow (i, j - 1);
      double rightVelocity  = (j == (N - 1)) ?
                              uVelocityBoundaryWindow (i, 1) :
                              uVelocityWindow (i, j);
      double bottomVelocity = (i == 0) ?
                              vVelocityBoundaryWindow (0, j) :
                              vVelocityWindow (i - 1, j);
      double topVelocity    = (i == (M - 1)) ?
                              vVelocityBoundaryWindow (1, j) :
                              vVelocityWindow (i, j);

      cellCenteredVVelocityWindow (i, j) = (topVelocity + bottomVelocity) / 2;
      cellCenteredUVelocityWindow (i, j) = (leftVelocity + rightVelocity) / 2;
      // Check if the cell-centered velocity is machine-epsilon, and make it true zero if so.
      if (abs (cellCenteredUVelocityWindow (i, j)) < 1E-10)
        cellCenteredUVelocityWindow (i, j) = 0;
      if (abs (cellCenteredVVelocityWindow (i, j)) < 1E-10)
        cellCenteredVVelocityWindow (i, j) = 0;
    }
  }

  #ifdef DEBUG
    cout << "<Half-time Cell-Centered Velocities Calculated>" << std::endl;
    cout << "<Cell-Centered U Velocity>" << endl;
    cellCenteredUVelocityWindow.displayMatrix();
    cout << "<Cell-Centered V Velocity>" << endl;
    cellCenteredVVelocityWindow.displayMatrix();
  #endif

  // Calculate temperatures at half-time level
  // Calculate half-time U-offset temperatures (Mx(N-1) lateral offset grid)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < (N - 1); ++j) {
      halfTimeUOffsetTemperatureWindow (i, j) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemmanFlag since the case where velocities do not share the
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;
      if (cellCenteredUVelocityWindow (i, j) >= 0 &&
          cellCenteredUVelocityWindow (i, j + 1) >= 0) {
        riemannFlag = 0b01;
      } else if (cellCenteredUVelocityWindow (i, j) <= 0 &&
                 cellCenteredUVelocityWindow (i, j + 1) <= 0) {
        riemannFlag = 0b10;
      } else {
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the left
      // neighboring cell.
      if (riemannFlag & 0b01) {
        if (j == 0) {
          leftNeighborT = temperatureWindow (i, j);
        } else {
          leftNeighborT = temperatureWindow (i, j - 1);
        }
        rightNeighborT = temperatureWindow (i, j + 1);

        halfTimeUOffsetTemperatureWindow (i, j) +=
            temperatureWindow (i, j) +
              (h / 2 - deltaT / 2 * cellCenteredUVelocityWindow (i, j)) *
                (rightNeighborT - leftNeighborT) / (2 * h);
      }

      // One or both velocities are negative. Take the flux from the right
      // neighboring cell.
      if (riemannFlag & 0b10) {
        leftNeighborT = temperatureWindow (i, j);
        if (j == (N - 2)) {
          rightNeighborT = temperatureWindow (i, j + 1);
        } else {
          rightNeighborT = temperatureWindow (i, j + 2);
        }

        halfTimeUOffsetTemperatureWindow (i, j) +=
            temperatureWindow (i, j + 1) -
              (h / 2 + deltaT / 2 * cellCenteredUVelocityWindow (i, j + 1)) *
                (rightNeighborT - leftNeighborT) / (2 * h);

      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells
      if (riemannFlag == 0b11) halfTimeUOffsetTemperatureWindow (i, j) /= 2;

    }
  }

  #ifdef DEBUG
    cout << "<Half-Time Offset Temperatures>" << endl;
    cout << "<Half-Time U-Offset Temperature>" << endl;
    halfTimeUOffsetTemperatureWindow.displayMatrix();
  #endif

  // Calculate half-time V-offset temperatures ((M-1)xN transverse offset grid)
  for (int i = 0; i < (M - 1); ++i) {
    for (int j = 0; j < N; ++j) {
      halfTimeVOffsetTemperatureWindow (i, j) = 0;
      // Flag to show which Riemann Problem solution to use. We use a bitflag
      // for the riemannFlag since the case where velocities do not share the
      // same sign is merely the average of the other two cases.
      int riemannFlag = 0b00;

      if (cellCenteredVVelocityWindow (i, j) >= 0 &&
          cellCenteredVVelocityWindow (i + 1, j) >= 0) {
        // One or both velocities are positive
        riemannFlag = 0b01;
      } else if (cellCenteredVVelocityWindow (i, j) <= 0 &&
                 cellCenteredVVelocityWindow (i + 1, j) <= 0) {
        // One or both velocities are negative
        riemannFlag = 0b10;
      } else {
        // One velocity is positive, one is negative.
        riemannFlag = 0b11;
      }

      // One or both velocities are positive. Take the flux from the bottom
      // neighboring cell
      if (riemannFlag & 0b01) {
        bottomNeighborT = temperatureWindow (i + 1, j);
        if (i == 0) {
          topNeighborT = temperatureWindow (i, j);
        } else {
          topNeighborT = temperatureWindow (i - 1, j);
        }

        halfTimeVOffsetTemperatureWindow (i, j) +=
            temperatureWindow (i, j) -
              (h / 2 + deltaT / 2 * cellCenteredVVelocityWindow (i, j)) *
               (topNeighborT - bottomNeighborT) / (2 * h);
      }

      // One or both velocities are negative. Take the flux from the top
      // neighboring cell
      if (riemannFlag & 0b10) {
        if (i == (M - 2)) {
          bottomNeighborT = temperatureWindow (i + 1, j);
        } else {
          bottomNeighborT = temperatureWindow (i + 2, j);
        }
        topNeighborT = temperatureWindow (i, j);

        halfTimeVOffsetTemperatureWindow (i, j) +=
            temperatureWindow (i + 1, j) +
              (h / 2 - deltaT / 2 * cellCenteredVVelocityWindow (i + 1, j)) *
               (topNeighborT - bottomNeighborT) / (2 * h);
      }

      // Velocities are in opposing directions. Take the average of the fluxes
      // from the left and right neighboring cells.
      if (riemannFlag == 0b11) {
        halfTimeVOffsetTemperatureWindow (i, j) /= 2;
      }
    }
  }

  #ifdef DEBUG
    cout << "<Half-Time V-Offset Temperature>" << endl;
    halfTimeVOffsetTemperatureWindow.displayMatrix();
  #endif

  double leftHalfTimeNeighborT, rightHalfTimeNeighborT,
         bottomHalfTimeNeighborT, topHalfTimeNeighborT;

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double diffusionWeight = 4;
      if (j == 0) {
        leftHalfTimeNeighborT = temperatureWindow (i, j);
        leftNeighborT         = temperatureWindow (i, j);
        diffusionWeight = 3;
      } else {
        leftHalfTimeNeighborT = halfTimeUOffsetTemperatureWindow (i, j - 1);
        leftNeighborT         = temperatureWindow (i, j - 1);
      }
      if (j == (N - 1)) {
        rightHalfTimeNeighborT = temperatureWindow (i, j);
        rightNeighborT         = temperatureWindow (i, j);
        diffusionWeight = 3;
      } else {
        rightHalfTimeNeighborT = halfTimeUOffsetTemperatureWindow (i, j);
        rightNeighborT         = temperatureWindow (i, j + 1);
      }

      if (i == 0) {
        bottomHalfTimeNeighborT = temperatureBoundaryWindow (0, j);
        bottomNeighborT         = temperatureBoundaryWindow (0, j);
      } else {
        bottomHalfTimeNeighborT = halfTimeVOffsetTemperatureWindow (i - 1, j);
        bottomNeighborT         = temperatureWindow (i - 1, j);
      }

      if (i == (M - 1)) {
        topHalfTimeNeighborT = temperatureBoundaryWindow (1, j);
        topNeighborT         = temperatureBoundaryWindow (1, j);
      } else {
        topHalfTimeNeighborT = halfTimeVOffsetTemperatureWindow (i, j);
        topNeighborT         = temperatureWindow (i + 1, j);
      }

      halfTimeTemperatureWindow (i, j) =
        (rightHalfTimeNeighborT + leftHalfTimeNeighborT +
         topHalfTimeNeighborT + bottomHalfTimeNeighborT) / 4 -
         (diffusionWeight * temperatureWindow (i, j) +
           leftNeighborT + rightNeighborT + bottomNeighborT + topNeighborT)
         * deltaT * diffusivity / (h * h);
    }
  }

  #ifdef DEBUG
    cout << "<Half-Time Temperature Data>" << endl;
    halfTimeTemperatureWindow.displayMatrix();
  #endif

  // Calculate half-time forcing
  if (forcingModel == "tauBenchmark") {
    // Benchmark taken from Tau (1991; JCP Vol. 99)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (i, j) = 3 * cos ((j + 1) * h) * sin ((i + 0.5) * h);

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (i, j) = -sin ((j + 0.5) * h) * cos ((i + 1) * h);

  } else if (forcingModel == "solCXBenchmark" ||
             forcingModel == "solKZBenchmark") {
    // solCX Benchmark taken from Kronbichler et al. (2011)
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N - 1; ++j)
        halfTimeUForcingWindow (i, j) = 0;

    for (int i = 0; i < M - 1; ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (i, j) = - sin((i + 0.5) * pi * h) * cos ((j + 1) * pi * h);

  } else if (forcingModel == "vorticalFlow") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < (N - 1); ++j)
        halfTimeUForcingWindow (i, j) = 3 * cos ((j + 1) * 2 * h) * sin ((i + 0.5) * 2 * h);

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j)
        halfTimeVForcingWindow (i, j) = -sin ((j + 0.5) * h) * cos ((i + 1) * h);

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
        halfTimeUForcingWindow (i, j) = 0;

    for (int i = 0; i < (M - 1); ++i)
      for (int j = 0; j < N; ++j) {
        halfTimeVForcingWindow (i, j) =  -1 * densityConstant *
                                          (1 - thermalExpansion *
                                           ((halfTimeTemperatureWindow (i, j) +
                                             halfTimeTemperatureWindow (i + 1, j)) / 2 -
                                            referenceTemperature));
      }
  } else {
    throw std::runtime_error("<Unexpected forcing model: \"" + forcingModel + "\" : Shutting down now>");
  }

  #ifdef DEBUG
    cout << "<Half-Time Forcing Data>" << endl;
    cout << "<Half-Time U Forcing Data>" << endl;
    halfTimeUForcingWindow.displayMatrix();
    cout << "<Half-Time V Forcing Data>" << endl;
    halfTimeVForcingWindow.displayMatrix();
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
      DataWindow<double> (geometry.getViscosityData(), M + 1, N + 1).displayMatrix();
    #endif
  }

  #ifdef DEBUG
    cout << "<Velocity Boundary Data>" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    uVelocityBoundaryWindow.displayMatrix();
    cout << "<V Velocity Boundary Data>" << endl;
    vVelocityBoundaryWindow.displayMatrix();
  #endif

  // Solve stokes at the half-time to find velocities
  halfTimeStokesSolnVector = solver.solve
           (forcingMatrix  * Map<VectorXd>(halfTimeForcingData, 2 * M * N - M - N) +
            boundaryMatrix * Map<VectorXd>(geometry.getVelocityBoundaryData(), 2 * M + 2 * N));

  for (int i = 0; i < M; i++)
    for (int j = 0; j < (N - 1); j++)
      if (abs(halfTimeUVelocityWindow (i, j)) < 10E-10)
        halfTimeUVelocityWindow (i, j) = 0;

  for (int i = 0; i < (M - 1); i++)
    for (int j = 0; j < N; ++j)
      if (abs (halfTimeVVelocityWindow (i, j)) < 10E-10)
        halfTimeVVelocityWindow (i, j) = 0;

  #ifdef DEBUG
    cout << "<Half-Time Velocities>" << endl;
    cout << "<Half-Time U Velocity>" << endl;
    halfTimeUVelocityWindow.displayMatrix();
    cout << "<Half-Time V Velocity>" << std::endl;
    halfTimeVVelocityWindow.displayMatrix();
  #endif

  double leftVelocity, rightVelocity, bottomVelocity, topVelocity;

  // Solve for full-time temperature
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        leftNeighborT = 0;
        leftVelocity = 0;
      } else {
        leftNeighborT = halfTimeUOffsetTemperatureWindow (i, j - 1);
        leftVelocity = halfTimeUVelocityWindow (i, j - 1);
      }
      if (j == (N - 1)) {
        rightNeighborT = 0;
        rightVelocity = 0;
      } else {
        rightNeighborT = halfTimeUOffsetTemperatureWindow (i, j);
        rightVelocity = halfTimeUVelocityWindow (i, j);
      }

      if (i == 0) {
        bottomNeighborT = 0;
        bottomVelocity = 0;
      } else {
        bottomNeighborT = halfTimeVOffsetTemperatureWindow (i - 1, j);
        bottomVelocity = halfTimeVVelocityWindow (i - 1, j);
      }

      if (i == (M - 1)) {
        topNeighborT = 0;
        topVelocity = 0;
      } else {
        topNeighborT = halfTimeVOffsetTemperatureWindow (i, j);
        topVelocity = halfTimeVVelocityWindow (i, j);
      }

      // Initialize the flux limiter to point to the desired limiter function.
      static double (ProblemStructure::*limiter) (double,double,double) = NULL;
      if (limiter == NULL) {
        if (fluxLimiter == "minmod") {
          limiter = &ProblemStructure::minmod;
        } else if (fluxLimiter == "superbee") {
          limiter = &ProblemStructure::superbee;
        } else if (fluxLimiter == "vanLeer") {
          limiter = &ProblemStructure::vanLeer;
        } else if (fluxLimiter == "none") {
          limiter = &ProblemStructure::minmod;
        }
      }

      double leftFirstOrderT = 0, rightFirstOrderT = 0,
             bottomFirstOrderT = 0, topFirstOrderT = 0,
             secondLeftFirstOrderT = 0, secondBottomFirstOrderT = 0;
      double leftFirstOrderVelocity, rightFirstOrderVelocity,
             bottomFirstOrderVelocity, topFirstOrderVelocity;
      double lateralFlux = 0, transverseFlux = 0;

      if (fluxLimiter != "none") {
        leftFirstOrderVelocity   = (j == 0) ?
                                   uVelocityBoundaryWindow (i, 0) :
                                   uVelocityWindow (i, j - 1);
        rightFirstOrderVelocity  = (j == (N - 1)) ?
                                   uVelocityBoundaryWindow (i, 1) :
                                   uVelocityWindow (i, j);
        bottomFirstOrderVelocity = (i == 0) ?
                                   vVelocityBoundaryWindow (0, j) :
                                   vVelocityWindow (i - 1, j);
        topFirstOrderVelocity    = (i == (M - 1)) ?
                                   vVelocityBoundaryWindow (1, j) :
                                   vVelocityWindow (i, j);

        double leftFlux = 0, rightFlux = 0;
        double topFlux = 0, bottomFlux = 0;

        if (j > 0) {
          if (leftFirstOrderVelocity < 0) {
            leftFlux = temporaryTemperature (i * N + j) * leftFirstOrderVelocity * deltaT / h;
            leftFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            leftFlux = temporaryTemperature (i * N + (j - 1)) * leftFirstOrderVelocity * deltaT / h;
            leftFirstOrderT = temporaryTemperature (i * N + (j - 1));
          }
        }

        if (j > 1) {
          secondLeftFirstOrderT = temporaryTemperature (i * N + (j - 2));
        } else {
          // Handle the weird symmetric boundary condition using abs magic
          secondLeftFirstOrderT = temporaryTemperature (i * N + abs (1 - j));
        }

        if (j < (N - 1)) {
          if (rightFirstOrderVelocity > 0) {
            rightFlux = temporaryTemperature (i * N + j) * rightFirstOrderVelocity * deltaT / h;
            rightFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            rightFlux = temporaryTemperature (i * N + (j + 1)) * rightFirstOrderVelocity * deltaT / h;
            rightFirstOrderT = temporaryTemperature (i * N + (j + 1));
          }
        }

        if (i > 0) {
          if (bottomFirstOrderVelocity < 0) {
            bottomFlux = temporaryTemperature (i * N + j) * bottomFirstOrderVelocity * deltaT / h;
            bottomFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            bottomFlux = temporaryTemperature ((i - 1) * N + j) * bottomFirstOrderVelocity * deltaT / h;
            bottomFirstOrderT = temporaryTemperature ((i - 1) * N + j);
          }
        }

        if (i > 1) {
          secondBottomFirstOrderT = temporaryTemperature ((i - 2) * N + j);
        } else {
          secondBottomFirstOrderT = temperatureBoundaryWindow (0, j);
        }

        if (i < (M - 1)) {
          if (topFirstOrderVelocity > 0) {
            topFlux = temporaryTemperature (i * N + j) * topFirstOrderVelocity * deltaT / h;
            topFirstOrderT = temporaryTemperature (i * N + j);
          } else {
            topFlux = temporaryTemperature ((i + 1) * N + j) * topFirstOrderVelocity * deltaT / h;
            topFirstOrderT = temporaryTemperature ((i + 1) * N + j);
          }
        }

        double leftPhi = (this->*limiter) (secondLeftFirstOrderT,
                                           leftFirstOrderT,
                                           temporaryTemperature (i * N + j));
        double rightPhi = (this->*limiter) (leftFirstOrderT,
                                            temporaryTemperature (i * N + j),
                                            rightFirstOrderT);
        double bottomPhi = (this->*limiter) (secondBottomFirstOrderT,
                                             bottomFirstOrderT,
                                             temporaryTemperature (i * N + j));
        double topPhi = (this->*limiter) (bottomFirstOrderT,
                                          temporaryTemperature (i * N + j),
                                          topFirstOrderT);

        lateralFlux =
            ((1 - leftPhi) * leftFlux + leftPhi * (leftVelocity * leftNeighborT * deltaT / h)) -
            ((1 - rightPhi) * rightFlux + rightPhi * rightVelocity * rightNeighborT * deltaT / h);
        transverseFlux = ((1 - bottomPhi) * bottomFlux + bottomPhi * bottomVelocity * bottomNeighborT * deltaT / h) -
                         ((1 - topPhi) * topFlux + topPhi * topVelocity * topNeighborT * deltaT / h);

        temperatureWindow (i, j) = temporaryTemperature (i * N + j) + transverseFlux + lateralFlux;
      } else
        temperatureWindow (i, j) = temporaryTemperature (i * N + j) + deltaT / h * (leftVelocity * leftNeighborT - rightVelocity * rightNeighborT) + deltaT / h * (bottomVelocity * bottomNeighborT - topVelocity * topNeighborT);

      if (std::isnan((double)temperatureWindow (i, j))) {
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
    temperatureWindow.displayMatrix();
  #endif
}
