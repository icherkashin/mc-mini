#pragma once

#include "parser/parser.h"

/** \brief A simple wrapper class for geometry-specific data
 *
 *  The GeometryStructure class encapsulates all of the geometry-specific data
 *  for a run, including. Specifically, it manages the size and position of the
 *  raw data arrays underlying the temperature, velocity, pressure, forcing,
 *  and viscosity fields
 */
class GeometryStructure {
  public:
    GeometryStructure(ParamParser& pp);
    ~GeometryStructure();

    // The number of rows in the geometry
    int getM();
    // The number of columns in the geometry
    int getN();

    // The full Stokes data array
    double * getStokesData();
    // The Stokes velocity data array
    double * getVelocityData();
    // The Stokes u-direction velocity data array
    double * getUVelocityData();
    // The Stokes v-direction velocity data array
    double * getVVelocityData();
    // The Stokes pressure data array
    double * getPressureData();

    // The velocity boundary data array
    double * getVelocityBoundaryData();
    // The u-direction velocity boundary data array
    double * getUVelocityBoundaryData();
    // The v-direction velocity boundary data array
    double * getVVelocityBoundaryData();

    // The forcing data array
    double * getForcingData();
    // The u-direction forcing data array
    double * getUForcingData();
    // The v-direction forcing data array
    double * getVForcingData();

    // The viscosity data array
    double * getViscosityData();

    // The temperature data array
    double * getTemperatureData();

    // The full temperature boundary data array
    double * getTemperatureBoundaryData();
    // The u-direction temperature boundary data array
    double * getUTemperatureBoundaryData();
    // The v-direction temperature boundary data array
    double * getVTemperatureBoundaryData();

  private:
    /** @defgroup GeoSizes Domain geometry sizes
     *  @name Domain Geometry Sizes
     *  @{
     */
    /// Number of rows in the domain
    int M;
    /// Number of columns in the domain
    int N;

    /** @} */

    /** @defgroup DataPointers
     *  @name Scalar Data Pointers
     *  @{
     */
    /// Stokes solution data containing velocity and pressure values.
    double * stokesData;
    /// Velocity boundary data containing prescribed velocity boundary values.
    double * velocityBoundaryData;

    /// Forcing data containing stokes equation forcing terms
    double * forcingData;

    /// Viscosity data
    double * viscosityData;

    /// Domain-interior temperature data
    double * temperatureData;
    /// Domain-boundary temperature data
    double * temperatureBoundaryData;
    /** @} */
};
