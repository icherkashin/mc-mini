#include <iostream>
#include <fstream>

#include "boost/lexical_cast.hpp"

#include "H5Cpp.h"

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "output/output.h"

#include <Eigen/Dense>

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif
using namespace Eigen;
using namespace std;

OutputStructure::OutputStructure (ParamParser&       pp,
                                  GeometryStructure& gs,
                                  ProblemStructure&  ps) :
    parser   (pp),
    geometry (gs),
    problem  (ps) {
  M  = geometry.getM(); 
  N  = geometry.getN();
  dx = problem.getH();

  if (parser.push ("outputParams")) {
    parser.queryParamString ("outputFormat", outputFormat, string("hdf5"));
    parser.queryParamString ("outputPath", outputPath, string("."));
    parser.queryParamString ("outputFilename", outputFilename, string("output"));

    parser.pop();
  }

  char s[128];
  sprintf (s, "test -e %s", outputPath.c_str());
  if (system (s) == 1) {
    sprintf (s, "mkdir %s", outputPath.c_str());
    system (s);
  }

  problemXdmfFile.open ((outputPath + "/" + outputFilename + "-series.xdmf").c_str(), ofstream::out);
  problemXdmfFile << "<?xml version=\"1.0\"?>" << endl
                  << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << endl
                  << "<Xdmf Version=\"2.0\">" << endl
                  << "  <Domain>" << endl
                  << "    <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">" << endl;
}

OutputStructure::~OutputStructure () {
  problemXdmfFile << "    </Grid>" << endl
                  << "  </Domain>" << endl
                  << "</Xdmf>" << endl;
  problemXdmfFile.close();
}

void OutputStructure::writeDefaultFile() {
  if (outputFormat == "hdf5") {
    this->writeHDF5File();
  } else {
    cerr << "Unknown output format " << outputFormat << " specified in parameters." << 
            "Shutting down now." << endl << endl;
    exit (-1);
  }
}

void OutputStructure::writeDefaultFile (const int timestep) {
  if (outputFormat == "hdf5") {
    this->writeHDF5File (timestep);
  } else {
    cerr << "Unknown output format " << outputFormat << " specified in parameters." << 
            "Shutting down now." << endl << endl;
    exit (-1);
  }
}

void OutputStructure::writeHDF5File() {
  
  H5File outputFile (H5std_string (outputPath + "/" +
                                   outputFilename + ".h5"), H5F_ACC_TRUNC);
  
  cout << "#==>> Outputting current data to file \"" << 
          outputPath << "/" << outputFilename << ".h5\"" << endl;
  
  ofstream xdmfFile ((outputPath + "/" + outputFilename + ".xdmf").c_str(), ofstream::out);
  
  xdmfFile << "<?xml version=\"1.0\"?>" << endl
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << endl
           << "<Xdmf Version=\"2.0\">" << endl
           << "  <Domain>" << endl
           << "    <Grid Name=\"mesh\" GridType=\"Uniform\">" << endl
           << "      <Topology TopologyType=\"2DCoRectMesh\" NumberOfElements=\"" << M + 1 << " " << N + 1 << "\"/>" << endl
           << "      <Geometry GeometryType=\"Origin_DxDy\">" << endl
           << "        <DataItem Dimensions=\"2\">" << endl
           << "          0 0" << endl
           << "        </DataItem>" << endl
           << "        <DataItem Dimensions=\"2\">" << endl
           << "          " << dx << " " << dx << endl
           << "        </DataItem>" << endl
           << "      </Geometry>" << endl;

  DataSet dataset;
  hsize_t dimsf[2];
  dimsf[0] = M; dimsf[1] = N;
  
  // Write temperature
  double * temperatureData = geometry.getTemperatureData();
  dataset = outputFile.createDataSet ("Temperature",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (temperatureData, 
                 PredType::NATIVE_DOUBLE);

  xdmfFile << "      <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Temperature" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;

  // Write pressure
  double * pressureData = geometry.getPressureData();
  dataset = outputFile.createDataSet ("Pressure",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (pressureData, 
                 PredType::NATIVE_DOUBLE);

  xdmfFile << "      <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Pressure" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;

  // Write U Velocity
  double * interpolatedUVelocityData = new double [M * N];
  double * uVelocityData = geometry.getUVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityData [i * N + j] = (uVelocityBoundaryData [i * 2]          +
                                                 uVelocityData         [i * (N - 1) + j]) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityData [i * N + j] = (uVelocityData         [i * (N - 1) + (j - 1)] +
                                                 uVelocityBoundaryData [i * 2       + 1      ]) / 2;
      } else {
        interpolatedUVelocityData [i * N + j] = (uVelocityData [i * (N - 1) + (j - 1)] +
                                                 uVelocityData [i * (N - 1) + j]      ) / 2;
      }
  dataset = outputFile.createDataSet ("UVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedUVelocityData,
                 PredType::NATIVE_DOUBLE);
  
  xdmfFile << "      <Attribute Name=\"UVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/UVelocity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
  
  // Write V Velocity
  double * interpolatedVVelocityData = new double [M * N];
  double * vVelocityData             = geometry.getVVelocityData();
  double * vVelocityBoundaryData     = geometry.getVVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityData [i * N + j] = (vVelocityBoundaryData [j]               +
                                                 vVelocityData         [i       * N + j]) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityBoundaryData [N + j]          ) / 2;
      } else {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityData         [i       * N + j]) / 2;
      }
  dataset = outputFile.createDataSet ("VVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedVVelocityData,
                 PredType::NATIVE_DOUBLE);
  
  xdmfFile << "      <Attribute Name=\"VVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/VVelocity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
                                                        

  double * velocityDivergence = new double[M * N];
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double uDivergence, vDivergence;

      if (i == 0) {
        vDivergence = (vVelocityBoundaryData[j] - vVelocityData[i * N + j]) / dx;
      } else if (i == (M - 1)) {
        vDivergence = (vVelocityData[(i - 1) * N + j] - vVelocityBoundaryData[M + j]) / dx;
      } else {
        vDivergence = (vVelocityData[(i - 1) * N + j] - vVelocityData[i * N + j]) / dx;
      }

      if (j == 0) {
        uDivergence = (uVelocityBoundaryData[2 * i] - uVelocityData[i * (N - 1) + j]) / dx;
      } else if (j == (N - 1)) {
        uDivergence = (uVelocityData[i * (N - 1) + (j - 1)] + uVelocityBoundaryData[2 * i + 1]) / dx;
      } else {
        uDivergence = (uVelocityData[i * (N - 1) + (j - 1)] - uVelocityData[i * (N - 1) + j]) / dx;
      }

      velocityDivergence[i * N + j] = uDivergence + vDivergence;
    }
  }

  dataset = outputFile.createDataSet ("Divergence",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (velocityDivergence,
                 PredType::NATIVE_DOUBLE);

  xdmfFile << "      <Attribute Name=\"Divergence\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Divergence" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
 
  dimsf[0] = (M + 1); dimsf[1] = (N + 1);
  
  // Write Viscosity
  double * viscosityData = geometry.getViscosityData();
  dataset = outputFile.createDataSet ("Viscosity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (viscosityData,
                 PredType::NATIVE_DOUBLE);

  xdmfFile << "      <Attribute Name=\"Viscosity\" AttributeType=\"Scalar\" Center=\"Node\">" << endl
           << "        <DataItem Dimensions=\"" << M + 1 << " " << N + 1 << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Viscosity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
  delete[] velocityDivergence;
  delete[] interpolatedUVelocityData;
  delete[] interpolatedVVelocityData;

  xdmfFile << "    </Grid>" << endl
           << "  </Domain>" << endl
           << "</Xdmf>" << endl;

  xdmfFile.close();
}

void OutputStructure::writeHDF5File (const int timestep) {
  H5File outputFile (H5std_string (outputPath + "/" +
                                   outputFilename + "-" +
                                   boost::lexical_cast<std::string> (timestep)  + 
                                   ".h5"),
                     H5F_ACC_TRUNC);

  cout << "#==>> Outputting current data to \"" <<
          outputPath << "/" << outputFilename << "-" << timestep << ".h5\"" << endl;

  problemXdmfFile << "      <Grid Name=\"mesh\" GridType=\"Uniform\">" << endl
                  << "        <Time Value=\"" << problem.getTime() << "\"/>" << endl
                  << "        <Topology TopologyType=\"2DCoRectMesh\" NumberOfElements=\"" << M + 1 << " " << N + 1<< "\"/>" << endl
                  << "        <Geometry GeometryType=\"Origin_DxDy\">" << endl
                  << "          <DataItem Dimensions=\"2\">" << endl
                  << "            0 0" << endl
                  << "          </DataItem>" << endl
                  << "          <DataItem Dimensions=\"2\">" << endl
                  << "            " << dx << " " << dx << endl
                  << "          </DataItem>" << endl
                  << "        </Geometry>" << endl;
  
  DataSet dataset;
  hsize_t dimsf[2];
  dimsf[0] = M; dimsf[1] = N;
  
  // Write temperature
  double * temperatureData = geometry.getTemperatureData();
  dataset = outputFile.createDataSet ("Temperature",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (temperatureData, 
                 PredType::NATIVE_DOUBLE);

  problemXdmfFile << "        <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Temperature" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  // Write pressure
  double * pressureData = geometry.getPressureData();
  dataset = outputFile.createDataSet ("Pressure",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (pressureData, 
                 PredType::NATIVE_DOUBLE);

  problemXdmfFile << "        <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Pressure" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  // Write U Velocity
  double * interpolatedUVelocityData = new double [M * N];
  double * uVelocityData = geometry.getUVelocityData();
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityData [i * N + j] = (uVelocityBoundaryData [i * 2]          +
                                                 uVelocityData         [i * (N - 1) + j]) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityData [i * N + j] = (uVelocityData         [i * (N - 1) + (j - 1)] +
                                                 uVelocityBoundaryData [i * 2       + 1      ]) / 2;
      } else {
        interpolatedUVelocityData [i * N + j] = (uVelocityData [i * (N - 1) + (j - 1)] +
                                                 uVelocityData [i * (N - 1) + j]      ) / 2;
      }
  dataset = outputFile.createDataSet ("UVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedUVelocityData,
                 PredType::NATIVE_DOUBLE);

  problemXdmfFile << "        <Attribute Name=\"UVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/UVelocity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;
  
  // Write V Velocity
  double * interpolatedVVelocityData = new double [M * N];
  double * vVelocityData             = geometry.getVVelocityData();
  double * vVelocityBoundaryData     = geometry.getVVelocityBoundaryData();
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityData [i * N + j] = (vVelocityBoundaryData [j]               +
                                                 vVelocityData         [i       * N + j]) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityBoundaryData [N + j]          ) / 2;
      } else {
        interpolatedVVelocityData [i * N + j] = (vVelocityData         [(i - 1) * N + j] +
                                                 vVelocityData         [i       * N + j]) / 2;
      }
  dataset = outputFile.createDataSet ("VVelocity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (interpolatedVVelocityData,
                 PredType::NATIVE_DOUBLE);
                                                        
  problemXdmfFile << "        <Attribute Name=\"VVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/VVelocity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  double * velocityDivergence = new double[M * N];
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double uDivergence, vDivergence;

      if (i == 0) {
        vDivergence = (vVelocityBoundaryData[j] - vVelocityData[i * N + j]) / dx;
      } else if (i == (M - 1)) {
        vDivergence = (vVelocityData[(i - 1) * N + j] - vVelocityBoundaryData[M + j]) / dx;
      } else {
        vDivergence = (vVelocityData[(i - 1) * N + j] - vVelocityData[i * N + j]) / dx;
      }

      if (j == 0) {
        uDivergence = (uVelocityBoundaryData[2 * i] - uVelocityData[i * (N - 1) + j]) / dx;
      } else if (j == (N - 1)) {
        uDivergence = (uVelocityData[i * (N - 1) + (j - 1)] - uVelocityBoundaryData[2 * i + 1]) / dx;
      } else {
        uDivergence = (uVelocityData[i * (N - 1) + (j - 1)] - uVelocityData[i * (N - 1) + j])  / dx;
      }

      velocityDivergence[i * N + j] = uDivergence + vDivergence;
    }
  }

  dataset = outputFile.createDataSet ("Divergence",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (velocityDivergence,
                 PredType::NATIVE_DOUBLE);

  problemXdmfFile << "        <Attribute Name=\"Divergence\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Divergence" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  dimsf[0] = (M + 1); dimsf[1] = (N + 1);
  
  // Write Viscosity
  double * viscosityData = geometry.getViscosityData();
  dataset = outputFile.createDataSet ("Viscosity",
                                      PredType::NATIVE_DOUBLE,
                                      DataSpace (2, dimsf));
  dataset.write (viscosityData,
                 PredType::NATIVE_DOUBLE);

  problemXdmfFile << "        <Attribute Name=\"Viscosity\" AttributeType=\"Scalar\" Center=\"Node\">" << endl
                  << "          <DataItem Dimensions=\"" << M + 1<< " " << N + 1 << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Viscosity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;
  
  delete[] velocityDivergence;
  delete[] interpolatedUVelocityData;
  delete[] interpolatedVVelocityData;

  problemXdmfFile << "      </Grid>" << endl;
}
