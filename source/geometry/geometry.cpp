#include "geometry/geometry.h"
#include "parser/parser.h"

/** @brief Constructs the GeometryStructure from parameters.
 *
 *  Constructs a GeometryStructure object from parameters parsed from the
 *  configuration parameter file.
 *
 *  Parameter specification:
 *  Section/Subsection | Name | Type | Description
 *  ------------------ | --------- | ---- | -----------
 *  geometryParams | M | int | The number of rows in the underlying representation of the problem domain (*required*)
 *  geometryParams | N | int | The number of columns in the underlying representation of the problem domain (*required*)
 */
GeometryStructure::GeometryStructure (ParamParser& parser) {
  // Grab the parameters from the required 'geometryParams' section
  if (parser.push ("geometryParams")) {
    // Read the number of rows (M) and columns (N)
    parser.getParamInt ("M", M);
    parser.getParamInt ("N", N);

    parser.pop();
  }


  /** The constructor allocates locations in memory for each of the value
   *  fields in the interior and boundaries of the problem domain. The layout
   *  of each of these problems on the problem domain grid are as follows:
   *  @verbatim
                             /----------------------------------------------------------------------\
                             | LEGEND                                                               |
          2   2   2   2      {======================================================================}
        8-5-8-5-8-5-8-5-8    | 1 : Interior cell-centered values (Temperature, Pressure)            |
      3 7 1 6 1 6 1 6 1 7 3  | 2 : Exterior vertical cell-centered values (Boundary temperature)    |
        8-4-8-4-8-4-8-4-8    | 3 : External horizontal cell-centered values (Boundary temperature)  |
      3 7 1 6 1 6 1 6 1 7 3  | 4 : Interior vertical cell-edge values (Velocity, Forcing)           |
        8-4-8-4-8-4-8-4-8    | 5 : Exterior vertical cell-edge values (Boundary velocity)           |
      3 7 1 6 1 6 1 6 1 7 3  | 6 : Interior horizontal cell-edge values (Velocity, Forcing)         |
        8-5-8-5-8-5-8-5-8    | 7 : Exterior horizontal cell-edge values (Boundary velocity)         |
          2   2   2   2      | 8 : Cell-corner values (Viscosity)                                   |
                             \----------------------------------------------------------------------/
      @endverbatim
   */
  stokesData = new double[M * (N - 1) + (M - 1) * N + M * N];

  velocityBoundaryData = new double[M *  2 + 2 * N];

  forcingData = new double[M * (N - 1) + (M - 1) * N];

  viscosityData = new double[(M + 1) * (N + 1)];

  temperatureData = new double[M * N];

  temperatureBoundaryData = new double[M * 2 + 2 * N];
}

/** @brief Deconstructs the GeometryStructure by deallocating all managed
 *         memory
 */
GeometryStructure::~GeometryStructure() {
  delete[] stokesData;
  delete[] velocityBoundaryData;
  delete[] forcingData;
  delete[] viscosityData;
  delete[] temperatureData;
  delete[] temperatureBoundaryData;
}

/// @brief Returns the number of rows in the problem domain
int GeometryStructure::getM() {
  return M;
}

/// @brief Returns the number of columns in the problem domain
int GeometryStructure::getN() {
  return N;
}

/** @brief Returns a pointer to the Stokes data
 *
 *  The **stokesData** variable contains the pointer to the U/V directional
 *  velocity and pressure data. The positions of these points on the domain
 *  grid are laid out as follows:
 *
 *  @verbatim
    +---+---+---+---+    /---------------------------------\
    | P U P U P U P |    | LEGEND                          |
    +-V-+-V-+-V-+-V-+    {=================================}
    | P U P U P U P |    | U : U-directional velocity data |
    +-V-+-V-+-V-+-V-+    | V : V-directional velocity data |
    | P U P U P U P |    | P : Pressure data               |
    +---+---+---+---+    \---------------------------------/
    @endverbatim
 *
 *  The actual layout of the data in memory is:
 *
 *  @verbatim
      M * (N - 1)   (M - 1) * N   M * N
    +-------------+-------------+----------+
    | U VELOCITY  | V VELOCITY  | PRESSURE |
    +-------------+-------------+----------+
    @endverbatim
 */
double * GeometryStructure::getStokesData() {
  return stokesData;
}

/// Returns a pointer to the velocity region of the Stokes data.
double * GeometryStructure::getVelocityData() {
  return stokesData;
}

/// Returns a pointer to the u-directional velocity region of the Stokes data
double * GeometryStructure::getUVelocityData() {
  return stokesData;
}

/// Returns a pointer to the v-directional velocity region of the Stokes data
double * GeometryStructure::getVVelocityData() {
  return stokesData + M * (N - 1);
}

/// Returns a pointer to the pressure region of the Stokes data
double * GeometryStructure::getPressureData() {
  return stokesData + M * (N - 1) + (M - 1) * N;
}

/** @brief Returns a pointer to the domain-boundary velocity data.
 *
 *  The **velocityBoundaryData** variable contains the pointer to the U/V
 *  directional boundary velocity data. The positions of these points on the
 *  domain grid are laid out as follows:
 *
 *  @verbatim
    +-V-+-V-+-V-+-V-+   /------------------------------------------\
    U   |   |   |   U   | LEGEND                                   |
    +---+---+---+---+   {==========================================}
    U   |   |   |   U   | U : U-directional boundary velocity data |
    +---+---+---+---+   | V : V-directional boundary velocity data |
    U   |   |   |   U   \------------------------------------------/
    +-V-+-V-+-V-+-V-+
    @endverbatim
 *
 *  The actual layout of the data in memory is:
 *  @verbatim
       2 * M                 2 * N
    +---------------------+---------------------+
    | U BOUNDARY VELOCITY | V BOUNDARY VELOCITY |
    +---------------------+---------------------+
    @endverbatim
 */
double * GeometryStructure::getVelocityBoundaryData() {
  return velocityBoundaryData;
}

/** @brief Returns a pointer to the u-directional domain-boundary velocity data.
 *
 *  The pointer is calculated using an offset of \f$0\f$ from the boundary
 *  velocity data pointer.
 */
double * GeometryStructure::getUVelocityBoundaryData() {
  return velocityBoundaryData;
}

/** @brief Returns a pointer to the v-directional domain-boundary velocity data.
 *
 *  The pointer is calculated using an offset of \f$2 * M\f$ from the boundary
 *  velocity data pointer.
 */
double * GeometryStructure::getVVelocityBoundaryData() {
  return velocityBoundaryData + 2 * M;
}

/** @brief Returns a pointer to the forcing data.
 *
 *  The **forcingData** variable contains the pointer to the U/V-directional
 *  forcing data. The positions of these points on the domain grid are laid
 *  out as follows:
 *
 *  @verbatim
    +---+---+---+---+   /--------------------------------\
    |   u   u   u   |   | LEGEND                         |
    +-v-+-v-+-v-+-v-+   {================================}
    |   u   u   u   |   | u : U-directional forcing data |
    +-v-+-v-+-v-+-v-+   | v : V-directional forcing data |
    |   u   u   u   |   \--------------------------------/
    +---+---+---+---+
    @endverbatim
 *
 *  The actual layout of the data in memory is:
 *  @verbatim
      M * (N - 1)   (M - 1) * N
    +-------------+-------------+
    |  U FORCING  |  V FORCING  |
    +-------------+-------------+
    @endverbatim
 */
double * GeometryStructure::getForcingData() {
  return forcingData;
}

/** @brief Returns a pointer to the u-directional forcing data
 *
 *  The pointer is calculated using an offset of \f$0\f$ from the forcing data
 *  pointer.
 */
double * GeometryStructure::getUForcingData() {
  return forcingData;
}

/** @brief Returns a pointer to the v-directional forcing data
 *
 *  The pointer is calculated using an offset of \f$M * (N - 1)\f$ from the forcing
 *  data pointer.
 */
double * GeometryStructure::getVForcingData() {
  return forcingData + M * (N - 1);
}

/** @brief Returns a pointer to the viscosity data.
 *
 *  The **viscosityData** variable contains the pointer to the viscosity
 *  data. The positions of these points on the domain grid are laid out as
 *  follows:
 *
 *  @verbatim
    M---M---M---M---M   /--------------------\
    |   |   |   |   |   | LEGEND             |
    M---M---M---M---M   {====================}
    |   |   |   |   |   | M : Viscosity data |
    M---M---M---M---M   \--------------------/
    |   |   |   |   |
    M---M---M---M---M
    @endverbatim
 *
 *  The actual layout of the data in memory is:
 *  @verbatim
      (M + 1) * (N + 1)
    +-------------------+
    |     VISCOSITY     |
    +-------------------+
    @endverbatim
 */
double * GeometryStructure::getViscosityData() {
  return viscosityData;
}

/** @brief Returns a pointer to the domain-interior temperature data
 *
 *  The **temperatureData** variable contains the pointer to the domain
 *  interior temperature data. The positions of these points on the domain
 *  grid are laid out as follows:
 *
 *  @verbatim
    +---+---+---+---+   /----------------------\
    | T | T | T | T |   | LEGEND               |
    +---+---+---+---+   {======================}
    | T | T | T | T |   | T : Temperature data |
    +---+---+---+---+   \----------------------/
    | T | T | T | T |
    +---+---+---+---+
    @endverbatim
 *
 *  The actual layout of the data in memory is:
 *  @verbatim
      M * N
    +-------------+
    | TEMPERATURE |
    +-------------+
    @endverbatim
 */
double * GeometryStructure::getTemperatureData() {
  return temperatureData;
}

/** @brief Returns a pointer to the domain-boundary temperature data.
 *
 *  The **temperatureBoundaryData** member variable contains the pointer to the
 *  domain boundary temperature data. The positions of these points on
 *  the domain grid are laid out as follows:
 *
 *  @verbatim
        t   t   t   t
      +---+---+---+---+     /-------------------------------\
    t |   |   |   |   | t   | LEGEND                        |
      +---+---+---+---+     {===============================}
    t |   |   |   |   | t   | t : Temperature boundary data |
      +---+---+---+---+     \-------------------------------/
    t |   |   |   |   | t
      +---+---+---+---+
        t   t   t   t
    @endverbatim
 */
double * GeometryStructure::getTemperatureBoundaryData() {
  return temperatureBoundaryData;
}

/** @brief Returns a pointer to the u-directional temperature boundary data.
 *
 *  The pointer is calculated using an offset of \f$0\f$ from the temperature
 *  boundary data pointer.
 */
double * GeometryStructure::getUTemperatureBoundaryData() {
  return temperatureBoundaryData;
}

/** @brief Returns a pointer to the v-directional temperature boundary data.
 *
 *  The pointer is calculated using an offset of \f$2 * M\f$ from the
 *  temperature boundary data pointer.
 */
double * GeometryStructure::getVTemperatureBoundaryData() {
  return temperatureBoundaryData + M * 2;
}
