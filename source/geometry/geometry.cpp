#include "geometry/geometry.h"
#include "parser/parser.h"

GeometryStructure::GeometryStructure (ParamParser& parser) {
  if (parser.push ("geometryParams")) {
    parser.getParamInt ("M", M);
    parser.getParamInt ("N", N);

    parser.pop();
  }
  //                                U             V             P/T        Viscosity
  stokesData              = new double[M * (N - 1) + (M - 1) * N + M * N];
  velocityBoundaryData    = new double[M *  2      +  2      * N];
  forcingData             = new double[M * (N - 1) + (M - 1) * N];
  viscosityData           = new double[                                    (M + 1) * (N + 1)];                      
  temperatureData         = new double[                            M * N];
  temperatureBoundaryData = new double[M * 2       + 2       * N];
}

GeometryStructure::~GeometryStructure() {
  delete[] stokesData;
  delete[] velocityBoundaryData;
  delete[] forcingData;
  delete[] viscosityData;
  delete[] temperatureData;
  delete[] temperatureBoundaryData;
}

int GeometryStructure::getM() {
  return M;
}

int GeometryStructure::getN() {
  return N;
}

// Stokes Data
double * GeometryStructure::getStokesData() {
  return stokesData;
}

// Velocity Data
double * GeometryStructure::getVelocityData() {
  return stokesData;
}

double * GeometryStructure::getUVelocityData() {
  return stokesData;
}

double * GeometryStructure::getVVelocityData() {
  return stokesData + M * (N - 1);
}

// Pressure Data
double * GeometryStructure::getPressureData() {
  return stokesData + M * (N - 1) + (M - 1) * N;
}

// Velocity Boundary Data
double * GeometryStructure::getVelocityBoundaryData() {
  return velocityBoundaryData;
}

double * GeometryStructure::getUVelocityBoundaryData() {
  return velocityBoundaryData;
}

double * GeometryStructure::getVVelocityBoundaryData() {
  return velocityBoundaryData + 2 * M;
}

// Forcing Data
double * GeometryStructure::getForcingData() {
  return forcingData;
}

double * GeometryStructure::getUForcingData() {
  return forcingData;
}

double * GeometryStructure::getVForcingData() {
  return forcingData + M * (N - 1);
}

// Viscosity Data
double * GeometryStructure::getViscosityData() {
  return viscosityData;
}

// Temperature Data
double * GeometryStructure::getTemperatureData() {
  return temperatureData;
}

// Temperature Boundary Data
double * GeometryStructure::getTemperatureBoundaryData() {
  return temperatureBoundaryData;
}

double * GeometryStructure::getUTemperatureBoundaryData() {
  return temperatureBoundaryData;
}

double * GeometryStructure::getVTemperatureBoundaryData() {
  return temperatureBoundaryData + M * 2;
}
