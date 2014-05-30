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
  double * temperatureData = geometry.getTemperatureData();

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
        temperatureData[j * N + i] = referenceTemperature;
   
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
        temperatureData[j * N + i] = referenceTemperature +
                                     sin ((i + 0.5) * h * xModes * M_PI / xExtent) * 
                                     sin ((j + 0.5) * h * yModes * M_PI / yExtent) * 
                                     temperatureScale;

  } else if (temperatureModel == "squareWave") {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j) {
        if ((M / 4 < j && j < 3 * M / 4) && (N / 4 < i && i < 3 * N / 4))
          temperatureData[j * N + i] = referenceTemperature + temperatureScale;
        else
          temperatureData[j * N + i] = referenceTemperature - temperatureScale;
      }
  } else {
    std::cerr << "<Unexpected temperature model: \"" << boundaryModel << "\" : Shutting down now!>" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Initialized temperature model as: \"" << temperatureModel << "\">" << endl;
    cout << "<Temperature Data>" << endl;
    cout << DataWindow<double> (temperatureData, N, M).displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeTemperatureBoundary() {
  double * temperatureBoundaryData = geometry.getTemperatureBoundaryData();

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

  for (int i = 0; i < N; ++i) {
    temperatureBoundaryData[i] = lowerTemperature;
    temperatureBoundaryData[N + i] = upperTemperature;
  }
}

void ProblemStructure::initializeVelocityBoundary() {
  double * uVelocityBoundaryData = geometry.getUVelocityBoundaryData();
  double * vVelocityBoundaryData = geometry.getVVelocityBoundaryData();


  if (boundaryModel == "tauBenchmark") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryData [i * 2 + j] = cos (j * N * h) * sin ((i + 0.5) * h);
  
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryData [i * N + j] = -sin ((j + 0.5) * h) * cos (i * M * h);
  } else if (boundaryModel == "solCXBenchmark" ||
             boundaryModel == "solKZBenchmark" ||
             boundaryModel == "noFlux") {
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < 2; ++j)
        uVelocityBoundaryData[i * 2 + j] = 0;
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < N; ++j)
        vVelocityBoundaryData[i * N + j] = 0;
  } else {
    cerr << "Unexpected boundary model: \"" << boundaryModel << "\" : Shutting down now!" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Initialized boundary model as: \"" << boundaryModel << "\">" << endl;
    cout << "<U Velocity Boundary Data>" << endl;
    cout << DataWindow<double> (uVelocityBoundaryData, 2, M).displayMatrix() << endl;
    cout << "<V Velocity Boundary Data>" << endl;
    cout << DataWindow<double> (vVelocityBoundaryData, N, 2).displayMatrix() << endl << endl;
  #endif
}

void ProblemStructure::initializeViscosity() {
  double * viscosityData = geometry.getViscosityData();

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
        viscosityData[i * (N + 1) + j] = viscosity;
  } else if (viscosityModel == "tauBenchmark") {
    viscosity = 1.0;
  } else if (viscosityModel == "solCXBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j) 
        viscosityData[i * (N + 1) + j] = (j <= N / 2) ? 1.0 : 1.0E06;
  } else if (viscosityModel == "solKZBenchmark") {
    for (int i = 0; i < (M + 1); ++i)
      for (int j = 0; j < (N + 1); ++j)
        viscosityData[i * (N + 1) + j] = 1.0 + j * h * 1.0E06;
  } else {
    cerr << "Unexpected viscosity model: \"" << viscosityModel << "\" : Shutting down now!" << endl;
    exit(-1);
  }

  #ifdef DEBUG
    cout << "<Viscosity model initialized as: \"" << viscosityModel << "\">" << endl;
    cout << "<Viscosity Data>" << endl;
    cout << DataWindow<double> (viscosityData, N + 1, M + 1).displayMatrix() << endl << endl;
  #endif
}
