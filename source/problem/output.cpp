#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "H5Cpp.h"

#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;
using namespace H5;

void ProblemStructure::outputPressure() {
  double * pressureData = geometry.getPressureData();

  cout << "Pressure:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(pressureData, M, N).colwise().reverse() << endl
       << endl;
}

void ProblemStructure::outputVelocity() {
  double * uVelocityData = geometry.getUVelocityData();
  double * vVelocityData = geometry.getVVelocityData();

  cout << "U Velocity:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(uVelocityData, M, (N - 1)).colwise().reverse() << endl << endl;

  cout << "V Velocity:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(vVelocityData, (M - 1), N).colwise().reverse() << endl << endl;
}

void ProblemStructure::outputBoundaryVelocity() {
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();

  cout << "U Boundary Velocity:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(uVelocityBoundaryData, M, 2).colwise().reverse() << endl
       << endl;

  cout << "V Boundary Velocity:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(vVelocityBoundaryData, 2, N).colwise().reverse() << endl
       << endl;
}

void ProblemStructure::outputForcing() {
  double * uForcingData = geometry.getUForcingData();
  double * vForcingData = geometry.getVForcingData();

  cout << "U Forcing:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(uForcingData, M, (N - 1)).colwise().reverse() << endl << endl;

  cout << "V Forcing:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(vForcingData, (M - 1), N).colwise().reverse() << endl << endl;
}

void ProblemStructure::outputViscosity() {
  double * viscosityData = geometry.getViscosityData();

  cout << "Viscosity:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(viscosityData, M + 1, N + 1).colwise().reverse() << endl << endl;
}

void ProblemStructure::outputTemperature() {
  double * temperatureData = geometry.getTemperatureData();

  
  cout << "Temperature:" << endl
       << Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(temperatureData, M, N).colwise().reverse() << endl << endl;
}
