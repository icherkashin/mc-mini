#include <iostream>
#include <cmath>
#include <stdexcept>

#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "debug.h"

#include "boost/math/constants/constants.hpp"
const double pi = boost::math::constants::pi<double>();

/*
 *
 * Data initialization routines
 *
 */

void ProblemStructure::initializeProblem() {
  initializeTimestep();
  initializeTemperature();
  initializeTemperatureBoundary();
  initializeVelocityBoundary();
  initializeViscosity();
}

void ProblemStructure::initializeTimestep() {
  deltaT = cfl * h / diffusivity;
  int nTimestep = (endTime - time) / deltaT;
  if (abs (nTimestep * deltaT + time - endTime) > 1E-06)
    deltaT = (endTime - time) / ++nTimestep;
  #ifdef DEBUG
    std::cout << "<Timestep initialized to " << deltaT << ">" << std::endl;
  #endif
}

void ProblemStructure::initializeTemperature() {
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), M, N);

  double referenceTemperature;
  double temperatureScale;

  if (parser.push ("problemParams")) {
    if (parser.push ("initialTemperatureParams")) {
      parser.queryParamDouble ("referenceTemperature", referenceTemperature, 273.15);
      parser.queryParamDouble ("temperatureScale",     temperatureScale,     100.0);

      parser.pop();
    }
    parser.pop();
  }

  if (temperatureModel == "constant") {

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        temperatureWindow (i, j) = referenceTemperature;

  } else if (temperatureModel == "sineWave") {
    int xModes;
    int yModes;

    if (parser.push ("problemParams")) {
      if (parser.tryPush ("initialTemperatureParams")) {
        parser.queryParamInt ("xModes", xModes, 2);
        parser.queryParamInt ("yModes", yModes, 2);

        parser.pop();
      }

      parser.pop();
    }

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        temperatureWindow (i, j) = referenceTemperature +
                                   sin ((i + 0.5) * h * xModes * pi / xExtent) *
                                   sin ((j + 0.5) * h * yModes * pi / yExtent) *
                                   temperatureScale;

  } else if (temperatureModel == "squareWave") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j) {
        if ((M / 4 < j && j < 3 * M / 4) && (N / 4 < i && i < 3 * N / 4))
          temperatureWindow (i, j) = referenceTemperature + temperatureScale;
        else
          temperatureWindow (i, j) = referenceTemperature;
      }
  } else if (temperatureModel == "circle") {
     double center_x;
     double center_y;
     double radius;
     if (parser.push ("problemParams")) {
       if (parser.tryPush ("initialTemperatureParams")) {
         parser.getParamDouble ("radius", radius);
         parser.getParamDouble ("xCenter", center_x);
         parser.getParamDouble ("yCenter", center_y);
         parser.pop();
       }
       parser.pop();
     }

     for (int i = 0; i < M; ++i)
       for (int j= 0; j < N; ++j) {
         if ( std::sqrt(std::pow((i*h+h/2)-(center_y),2.0) + std::pow((j*h+h/2)-(center_x),2.0))  < radius )
           temperatureWindow (i, j) = referenceTemperature + temperatureScale;
         else
           temperatureWindow (i, j) = referenceTemperature;
       }
  } else {
    throw std::invalid_argument("<Unexpected temperature model: \"" + temperatureModel + "\" : Shutting down now>");
  }

  #ifdef DEBUG
    cout << "<Initialized temperature model as: \"" << temperatureModel << "\">" << endl;
    cout << "<Temperature Data>" << endl;
    temperatureWindow.displayMatrix();
  #endif
}

void ProblemStructure::initializeTemperatureBoundary() {
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), 2, N);

  double upperTemperature;
  double lowerTemperature;

  if (parser.push ("problemParams")) {
    if (parser.push ("temperatureBoundaryParams")) {
      parser.getParamDouble ("upperBoundaryTemperature", upperTemperature);
      parser.getParamDouble ("lowerBoundaryTemperature", lowerTemperature);

      parser.pop();
    }

    parser.pop();
  }

  for (int j = 0; j < N; ++j) {
    temperatureBoundaryWindow (0, j) = lowerTemperature;
    temperatureBoundaryWindow (1, j) = upperTemperature;
  }
}

void ProblemStructure::initializeVelocityBoundary() {
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), M, 2);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), 2, N);


  if (boundaryModel == "tauBenchmark") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryWindow (i, j) = cos (j * N * h) * sin ((i + 0.5) * h);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryWindow (i, j) = -sin ((j + 0.5) * h) * cos (i * M * h);
  } else if (boundaryModel == "solCXBenchmark" ||
             boundaryModel == "solKZBenchmark" ||
             boundaryModel == "noFlux") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryWindow (i, j) = 0;
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryWindow (i, j) = 0;
  } else {
    throw std::invalid_argument("<Unexpected boundary model: \"" + boundaryModel + "\" : Shutting down now>");
  }

  #ifdef DEBUG
    cout << "<Initialized boundary model as: \"" << boundaryModel << "\">" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    uVelocityBoundaryWindow.displayMatrix();
    cout << "<V Velocity Boundary Data>" << endl;
    vVelocityBoundaryWindow.displayMatrix();
  #endif
}

void ProblemStructure::initializeViscosity() {
  DataWindow<double> viscosityWindow (geometry.getViscosityData(), M + 1, N + 1);

  double viscosity;

  if (viscosityModel == "constant") {
    if (parser.push ("problemParams")) {
      if (parser.tryPush ("initialViscosity")) {
        parser.queryParamDouble ("viscosityScale", viscosity, 1.0);

        parser.pop();
      } else {
        viscosity = 1.0;
      }
      parser.pop();
    }

    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j)
        viscosityWindow (i, j) = viscosity;
  } else if (viscosityModel == "tauBenchmark") {
    viscosity = 1.0;
  } else if (viscosityModel == "solCXBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j)
        viscosityWindow (i, j) = (j <= N / 2) ? 1.0 : 1.0E06;
  } else if (viscosityModel == "solKZBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j)
        viscosityWindow (i, j) = 1.0 + j * h * 1.0E06;
  } else {
    throw std::invalid_argument("Unexpected viscosity model: \"" + viscosityModel + "\" : Shutting down now!");
  }

  #ifdef DEBUG
    cout << "<Viscosity model initialized as: \"" << viscosityModel << "\">" << endl;
    cout << "<Viscosity Data>" << endl;
    viscosityWindow.displayMatrix();
  #endif
}
