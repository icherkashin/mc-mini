#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>

#include "matrixForms/sparseForms.h"
#include "geometry/dataWindow.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"
#include "debug.h"

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
  DataWindow<double> temperatureWindow (geometry.getTemperatureData(), N, M);

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
    
    for (int i = 0; i < N; ++i) 
      for (int j = 0; j < M; ++j) 
        temperatureWindow (j, i) = referenceTemperature;
   
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

    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        temperatureWindow (j, i) = referenceTemperature +
                                   sin ((i + 0.5) * h * xModes * M_PI / xExtent) * 
                                   sin ((j + 0.5) * h * yModes * M_PI / yExtent) * 
                                   temperatureScale;

  } else if (temperatureModel == "squareWave") {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j) { 
        if ((M / 4 < j && j < 3 * M / 4) && (N / 4 < i && i < 3 * N / 4))
          temperatureWindow (j, i) = referenceTemperature + temperatureScale;
        else
          temperatureWindow (j, i) = referenceTemperature;
      }
  } else {
    cerr << "<Unexpected temperature model: \"" << boundaryModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Initialized temperature model as: \"" << temperatureModel << "\">" << endl;
    cout << "<Temperature Data>" << endl;
    cout << temperatureWindow.displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeTemperatureBoundary() {
  DataWindow<double> temperatureBoundaryWindow (geometry.getTemperatureBoundaryData(), N, 2);

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
    temperatureBoundaryWindow (j, 0) = lowerTemperature;
    temperatureBoundaryWindow (j, 1) = upperTemperature;
  }
}

void ProblemStructure::initializeVelocityBoundary() {
  DataWindow<double> uVelocityBoundaryWindow (geometry.getUVelocityBoundaryData(), 2, M);
  DataWindow<double> vVelocityBoundaryWindow (geometry.getVVelocityBoundaryData(), N, 2);


  if (boundaryModel == "tauBenchmark") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryWindow (j, i) = cos (j * N * h) * sin ((i + 0.5) * h);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryWindow (j, i) = -sin ((j + 0.5) * h) * cos (i * M * h);
  } else if (boundaryModel == "solCXBenchmark" ||
             boundaryModel == "solKZBenchmark" ||
             boundaryModel == "noFlux") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryWindow (j, i) = 0;
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryWindow (j, i) = 0;
  } else {
    cerr << "<Unexpected boundary model: \"" << boundaryModel << "\" : Shutting down now>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Initialized boundary model as: \"" << boundaryModel << "\">" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    cout << uVelocityBoundaryWindow.displayMatrix() << endl;
    cout << "<V Velocity Boundary Data>" << endl;
    cout << vVelocityBoundaryWindow.displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeViscosity() {
  DataWindow<double> viscosityWindow (geometry.getViscosityData(), N + 1, M + 1);

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
        viscosityWindow (j, i) = viscosity;
  } else if (viscosityModel == "tauBenchmark") {
    viscosity = 1.0;
  } else if (viscosityModel == "solCXBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j) 
        viscosityWindow (j, i) = (j <= N / 2) ? 1.0 : 1.0E06;
  } else if (viscosityModel == "solKZBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j)
        viscosityWindow (j, i) = 1.0 + j * h * 1.0E06;
  } else {
    cerr << "Unexpected viscosity model: \"" << viscosityModel << "\" : Shutting down now!" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Viscosity model initialized as: \"" << viscosityModel << "\">" << endl;
    cout << "<Viscosity Data>" << endl;
    cout << viscosityWindow.displayMatrix() << endl << endl;
  #endif
}
