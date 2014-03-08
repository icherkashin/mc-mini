#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "H5Cpp.h"

#include "geometry/geometry.h"
#include "problem/problem.h"

using namespace Eigen;
using namespace std;
using namespace H5;

void ProblemStructure::outputH5() {
  hid_t     file_id;
  file_id = H5Fcreate ("output.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  double *x = new double[N + 1];
  double *y = new double[M + 1];
  for (int i = 0; i < N + 1; ++i)
    x[i] = h * i;
  for (int j = 0; j < M + 1; ++j)
    y[j] = h * j;

  double *temperatureData = geometry.getTemperatureData();

  double *pressureData    = geometry.getPressureData();
  
  double *viscosityData   = geometry.getViscosityData();

  double * interpolatedViscosityData = new double[M * N];
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      interpolatedViscosityData[i * M + j] = (viscosityData[i       * (M + 1) + j] + 
                                              viscosityData[i       * (M + 1) + (j + 1)] + 
                                              viscosityData[(i + 1) * (M + 1) + j] + 
                                              viscosityData[(i + 1) * (M + 1) + (j + 1)]) / 4;
    }
  }

  double * uVelocityData = geometry.getUVelocityData();
  double * uBoundaryVelocityData = geometry.getUVelocityBoundaryData();
  double * interpolatedUVelocityData = new double[M * N];
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (j == 0) {
        interpolatedUVelocityData [i * N + j] = (uBoundaryVelocityData [i * 2]     +
                                                 uVelocityData         [i * (N - 1) + j]) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityData [i * N + j] = (uVelocityData         [i * (N - 1) + (j - 1)] + 
                                                 uBoundaryVelocityData [i * 2 + 1]      ) / 2;
      } else {
        interpolatedUVelocityData [i * N + j] = (uVelocityData [i * (N - 1) + (j - 1)] +
                                                 uVelocityData [i * (N - 1) + j]      ) / 2;
      }
    }
  }  

  double * vVelocityData = geometry.getVVelocityData();
  double * vBoundaryVelocityData = geometry.getVVelocityBoundaryData();
  double * interpolatedVVelocityData = new double[M * N];

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == 0) {
        interpolatedVVelocityData[i * N + j] = (vBoundaryVelocityData [        j] + 
                                                vVelocityData         [i * N + j]) / 2; 
      } else if (i == (M - 1)) {
        interpolatedVVelocityData[i * N + j] = (vVelocityData         [(i - 1) * N + j] + 
                                                vBoundaryVelocityData [          N + j]) / 2;
      } else {
        interpolatedVVelocityData[i * N + j] = (vVelocityData [(i - 1) * N + j] +
                                                vVelocityData [i       * N + j]) / 2.0;
      }

      cerr << interpolatedVVelocityData[i * N + j] << endl;
    }
  }

  hid_t   dataset_id, dataspace_id;
  hsize_t dims[3];
  herr_t  status;

  hsize_t xdim[3];
  hsize_t ydim[3];
  hsize_t zdim[3];
  
  xdim[0] = N + 1;
  ydim[0] = M + 1;
  
  dataspace_id = H5Screate_simple(1, xdim, NULL);

    dataset_id = H5Dcreate1 (file_id, 
                             "/X", 
                             H5T_NATIVE_DOUBLE, 
                             dataspace_id, 
                             H5P_DEFAULT);

    status = H5Dwrite (dataset_id, 
                       H5T_NATIVE_DOUBLE, 
                       H5S_ALL, H5S_ALL, 
                       H5P_DEFAULT, 
                       x);

    status = H5Dclose (dataset_id);

  status = H5Sclose (dataspace_id);

  dataspace_id = H5Screate_simple (1, ydim, NULL);
  
    dataset_id = H5Dcreate1 (file_id,
                             "/Y",
                             H5T_NATIVE_DOUBLE,
                             dataspace_id,
                             H5P_DEFAULT);

    status = H5Dwrite (dataset_id,
                       H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL,
                       H5P_DEFAULT,
                       y);

    status = H5Dclose (dataset_id);

  status = H5Sclose (dataspace_id);

  dims[0] = M;
  dims[1] = N;
  dataspace_id = H5Screate_simple (2, dims, NULL);

  dataset_id = H5Dcreate1 (file_id, 
                           "/Temperature", 
                           H5T_NATIVE_DOUBLE, 
                           dataspace_id, 
                           H5P_DEFAULT);
  status = H5Dwrite (dataset_id, 
                     H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, 
                     H5P_DEFAULT, 
                     temperatureData);
  status = H5Dclose (dataset_id);

  dataset_id = H5Dcreate1 (file_id,
                           "/Pressure",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT);
  status = H5Dwrite (dataset_id, 
                     H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, 
                     H5P_DEFAULT, 
                     pressureData);
  status = H5Dclose (dataset_id);

  dataset_id = H5Dcreate1 (file_id,
                           "/Viscosity",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT);
  status = H5Dwrite (dataset_id, 
                     H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, 
                     H5P_DEFAULT, 
                     pressureData);
  status = H5Dclose (dataset_id);

  dataset_id = H5Dcreate1 (file_id,
                           "/UVelocity",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT);
  status = H5Dwrite (dataset_id, 
                     H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, 
                     H5P_DEFAULT, 
                     interpolatedUVelocityData);
  status = H5Dclose (dataset_id);

  dataset_id = H5Dcreate1 (file_id,
                           "/VVelocity",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT);
  status  = H5Dwrite (dataset_id, 
                      H5T_NATIVE_DOUBLE, 
                      H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT,
                      interpolatedVVelocityData);
  status = H5Dclose (dataset_id);

  status = H5Sclose (dataspace_id);

  status = H5Fclose (file_id);

  delete[] interpolatedViscosityData;
  delete[] interpolatedUVelocityData;
  delete[] interpolatedVVelocityData;
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
