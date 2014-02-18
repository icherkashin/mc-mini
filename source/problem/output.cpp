#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "hdf5.h"
#include "H5Cpp.h"

#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;
using namespace H5;

void ProblemStructure::outputH5() {
  H5File file (outputFile.c_str(), H5F_ACC_TRUNC);
  
  hsize_t temperatureDim[] = {M, N};
  DataSpace temperatureSpace (2, temperatureDim);

  hsize_t pressureDim[] =    {M, N};
  DataSpace pressureSpace    (2, pressureDim);

  double * viscosityData = geometry.getViscosityData();
  double * interpolatedViscosity = new double[M * N];

  double * uVelocityData = geometry.getUVelocityData();
  double * uBoundaryVelocityData = geometry.getUVelocityBoundaryData();
  double * vVelocityData = geometry.getVVelocityData();
  double * vBoundaryVelocityData = geometry.getVVelocityBoundaryData();

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      interpolatedViscosity[i * M + j] = (viscosityData[i       * (M + 1) + j] + 
                                          viscosityData[i       * (M + 1) + (j + 1)] + 
                                          viscosityData[(i + 1) * (M + 1) + j] + 
                                          viscosityData[(i + 1) * (M + 1) + (j + 1)]) / 4;
    }
  }

  hsize_t viscosityDim[] = {M, N};
  DataSpace viscositySpace (2, viscosityDim);
  

  DataSet temperatureDataset = file.createDataSet ("Temperature", 
                                                   PredType::NATIVE_DOUBLE, 
                                                   temperatureSpace);

  DataSet pressureDataset =    file.createDataSet ("Pressure",
                                                   PredType::NATIVE_DOUBLE,
                                                   pressureSpace);

  DataSet viscosityDataset =   file.createDataSet ("Viscosity",
                                                   PredType::NATIVE_DOUBLE,
                                                   viscositySpace);

  temperatureDataset.write (geometry.getTemperatureData(), PredType::NATIVE_DOUBLE);
  pressureDataset.write    (geometry.getPressureData(),    PredType::NATIVE_DOUBLE);
  viscosityDataset.write   (interpolatedViscosity,         PredType::NATIVE_DOUBLE);
  
  delete[] interpolatedViscosity;
}

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