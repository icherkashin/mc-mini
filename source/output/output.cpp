#include <iostream>
#include <fstream>

#include "boost/lexical_cast.hpp"

#include "hdf5.h"

#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "output/output.h"

#include <Eigen/Dense>

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
    if (system (s) == -1) {
      cout << "<Error: couldn't create directory " << outputPath << ">" << endl;
      exit (-1);
    }
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
  
  hid_t outputFile = H5Fcreate((outputPath + "/" + outputFilename + ".h5").c_str(), 
                               H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  cout << "<Outputting current data to file \"" << outputPath << "/" << outputFilename << ".h5\">" << endl;
  
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

  hid_t dataset, datatype, dataspace;
  herr_t status;
  hsize_t dimsf[2];
  dimsf[0] = M; dimsf[1] = N;
  
  dataspace = H5Screate_simple (2, dimsf, NULL);
  datatype = H5Tcopy (H5T_NATIVE_DOUBLE);
  status = H5Tset_order (datatype, H5T_ORDER_LE);

  // Write temperature
  dataset = H5Dcreate (outputFile, "Temperature", datatype, dataspace, 
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getTemperatureData());

  xdmfFile << "      <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Temperature" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;

  // Write pressure
  dataset = H5Dcreate (outputFile, "Pressure", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getPressureData());

  xdmfFile << "      <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Pressure" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;

  // Write U Velocity
  static double * interpolatedUVelocityData = new double [M * N];
  static DataWindow<double> interpolatedUVelocityWindow (interpolatedUVelocityData, N, M);
  static DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  static DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityWindow (j, i) = (uVelocityBoundaryWindow (0, i) +
                                              uVelocityWindow         (j, i)) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityWindow (j, i) = (uVelocityWindow         (j - 1, i) +
                                              uVelocityBoundaryWindow (1, i)) / 2;
      } else {
        interpolatedUVelocityWindow (j, i) = (uVelocityWindow (j - 1, i) +
                                              uVelocityWindow (j,     i)) / 2;
      }
  dataset = H5Dcreate (outputFile, "UVelocity", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, interpolatedUVelocityData);
 
  
  xdmfFile << "      <Attribute Name=\"UVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/UVelocity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
  
  // Write V Velocity
  static double * interpolatedVVelocityData = new double [M * N];
  static DataWindow<double> interpolatedVVelocityWindow (interpolatedVVelocityData, N, M);
  static DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  static DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);  
  
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityWindow (j, i) = (vVelocityBoundaryWindow (j, 0) +
                                              vVelocityWindow         (j, i)) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityWindow (j, i) = (vVelocityWindow         (j, (i - 1)) +
                                              vVelocityBoundaryWindow (j, 1)) / 2;
      } else {
        interpolatedVVelocityWindow (j, i) = (vVelocityWindow         (j, (i - 1)) +
                                              vVelocityWindow         (j, i)) / 2;
      }
  dataset = H5Dcreate (outputFile, "VVelocity", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, interpolatedVVelocityData);
  
  xdmfFile << "      <Attribute Name=\"VVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/VVelocity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
                                                        

  static double * velocityDivergenceData = new double[M * N];
  static DataWindow<double> velocityDivergenceWindow (velocityDivergenceData, N, M);

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double uDivergence, vDivergence;

      if (i == 0) {
        vDivergence = (vVelocityBoundaryWindow (j, 0) - vVelocityWindow (j, i)) / dx;
      } else if (i == (M - 1)) {
        vDivergence = (vVelocityWindow (j, i - 1) - vVelocityBoundaryWindow (j, 1)) / dx;
      } else {
        vDivergence = (vVelocityWindow (j, i - 1) - vVelocityWindow (j, i)) / dx;
      }

      if (j == 0) {
        uDivergence = (uVelocityBoundaryWindow (0, i) - uVelocityWindow (j, i)) / dx;
      } else if (j == (N - 1)) {
        uDivergence = (uVelocityWindow (j - 1, i) - uVelocityBoundaryWindow (1, i)) / dx;
      } else {
        uDivergence = (uVelocityWindow (j - 1, i) - uVelocityWindow (j, i)) / dx;
      }

      velocityDivergenceWindow (j, i) = uDivergence + vDivergence;
    }
  }
  
  dataset = H5Dcreate (outputFile, "Divergence", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, velocityDivergenceData);

  xdmfFile << "      <Attribute Name=\"Divergence\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
           << "        <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Divergence" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl;
 
  dimsf[0] = (M + 1); dimsf[1] = (N + 1);
  
  dataspace = H5Screate_simple (2, dimsf, NULL);
  datatype = H5Tcopy (H5T_NATIVE_DOUBLE);
  status = H5Tset_order (datatype, H5T_ORDER_LE);
  
  // Write Viscosity
  dataset = H5Dcreate (outputFile, "Viscosity", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getViscosityData());
   
  xdmfFile << "      <Attribute Name=\"Viscosity\" AttributeType=\"Scalar\" Center=\"Node\">" << endl
           << "        <DataItem Dimensions=\"" << M + 1 << " " << N + 1 << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
           << "          " << outputFilename + ".h5:/Viscosity" << endl
           << "        </DataItem>" << endl
           << "      </Attribute>" << endl
           << "    </Grid>" << endl
           << "  </Domain>" << endl
           << "</Xdmf>" << endl;

  H5Sclose (dataspace);
  H5Tclose (datatype);
  H5Dclose (dataset);
  H5Fclose (outputFile);

  xdmfFile.close();
}

void OutputStructure::writeHDF5File (const int timestep) {
  hid_t outputFile = H5Fcreate ((outputPath + "/" + outputFilename + "-" + boost::lexical_cast<std::string> (timestep) + ".h5").c_str(),
                     H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  cout << "<Outputting current data to \"" << outputPath << "/" << outputFilename << "-" << timestep << ".h5\">" << endl;

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
  
  hid_t dataset, datatype, dataspace;
  herr_t status;
  hsize_t dimsf[2];
  dimsf[0] = M; dimsf[1] = N;
  
  dataspace = H5Screate_simple (2, dimsf, NULL);
  datatype = H5Tcopy (H5T_NATIVE_DOUBLE);
  status = H5Tset_order (datatype, H5T_ORDER_LE);

  // Write temperature
  dataset = H5Dcreate (outputFile, "Temperature", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getTemperatureData());

  problemXdmfFile << "        <Attribute Name=\"Temperature\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Temperature" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  // Write pressure
  dataset = H5Dcreate (outputFile, "Pressure", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getPressureData());  

  problemXdmfFile << "        <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Pressure" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  // Write U Velocity
  static double * interpolatedUVelocityData = new double [M * N];
  static DataWindow<double> interpolatedUVelocityWindow (interpolatedUVelocityData, N, M);
  static DataWindow<double> uVelocityWindow (geometry.getUVelocityData(), N - 1, M);
  static DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);

  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (j == 0) {
        interpolatedUVelocityWindow (j, i) = (uVelocityBoundaryWindow (0, i) +
                                              uVelocityWindow         (j, i)) / 2;
      } else if (j == (N - 1)) {
        interpolatedUVelocityWindow (j, i) = (uVelocityWindow         (j - 1, i) +
                                              uVelocityBoundaryWindow (1, i)) / 2;
      } else {
        interpolatedUVelocityWindow (j, i) = (uVelocityWindow (j - 1, i) +
                                              uVelocityWindow (j,     i)) / 2;
      }

  dataset = H5Dcreate (outputFile, "UVelocity", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, interpolatedUVelocityData);

  problemXdmfFile << "        <Attribute Name=\"UVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/UVelocity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;
  
  // Write V Velocity
  static double * interpolatedVVelocityData = new double [M * N];
  static DataWindow<double> interpolatedVVelocityWindow (interpolatedVVelocityData, N, M);
  static DataWindow<double> vVelocityWindow (geometry.getVVelocityData(), N, M - 1);
  static DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);
  
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (i == 0) {
        interpolatedVVelocityWindow (j, i) = (vVelocityBoundaryWindow (j, 0) +
                                              vVelocityWindow         (j, i)) / 2;
      } else if (i == (M - 1)) {
        interpolatedVVelocityWindow (j, i) = (vVelocityWindow         (j, i - 1) +
                                              vVelocityBoundaryWindow (j, 1)) / 2;
      } else {
        interpolatedVVelocityWindow (j, i) = (vVelocityWindow (j, i - 1) +
                                              vVelocityWindow (j, i)) / 2;
      }

  dataset = H5Dcreate (outputFile, "VVelocity", datatype, dataspace, 
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, interpolatedVVelocityData);
  
  problemXdmfFile << "        <Attribute Name=\"VVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/VVelocity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  static double * velocityDivergenceData = new double[M * N];
  static DataWindow<double> velocityDivergenceWindow (velocityDivergenceData, N, M);

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double uDivergence, vDivergence;

      if (i == 0) {
        vDivergence = (vVelocityBoundaryWindow (j, 0) - vVelocityWindow (j, i)) / dx;
      } else if (i == (M - 1)) {
        vDivergence = (vVelocityWindow (j, i - 1) - vVelocityBoundaryWindow (j, 1)) / dx;
      } else {
        vDivergence = (vVelocityWindow (j, i - 1) - vVelocityWindow (j, i)) / dx;
      }

      if (j == 0) {
        uDivergence = (uVelocityBoundaryWindow (0, i) - uVelocityWindow (j, i)) / dx;
      } else if (j == (N - 1)) {
        uDivergence = (uVelocityWindow (j - 1, i) - uVelocityBoundaryWindow (1, i)) / dx;
      } else {
        uDivergence = (uVelocityWindow (j - 1, i) - uVelocityWindow (j, i)) / dx;
      }

      velocityDivergenceWindow (j, i) = uDivergence + vDivergence;
    }
  }

  dataset = H5Dcreate (outputFile, "Divergence", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, velocityDivergenceData);

  problemXdmfFile << "        <Attribute Name=\"Divergence\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl
                  << "          <DataItem Dimensions=\"" << M << " " << N << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Divergence" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl;

  dimsf[0] = (M + 1); dimsf[1] = (N + 1);
 
  dataspace = H5Screate_simple (2, dimsf, NULL);
  datatype = H5Tcopy (H5T_NATIVE_DOUBLE);
  status = H5Tset_order (datatype, H5T_ORDER_LE);
 
  // Write Viscosity
  dataset = H5Dcreate (outputFile, "Viscosity", datatype, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, geometry.getViscosityData());

  problemXdmfFile << "        <Attribute Name=\"Viscosity\" AttributeType=\"Scalar\" Center=\"Node\">" << endl
                  << "          <DataItem Dimensions=\"" << M + 1<< " " << N + 1 << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << endl
                  << "            " << outputFilename + "-" + boost::lexical_cast<std::string> (timestep) << ".h5:/Viscosity" << endl
                  << "          </DataItem>" << endl
                  << "        </Attribute>" << endl
                  << "      </Grid>" << endl;

  H5Sclose (dataspace);
  H5Tclose (datatype);
  H5Dclose (dataset);
  H5Fclose (outputFile);  
}
