#pragma once

#include "parser/parser.h"

class ProblemStructure {
  public:
    ProblemStructure (ParamParser& pp, GeometryStructure& gs);
    ~ProblemStructure();

    // Initialize problem (set dt, temperature, boundary.)
    void initializeProblem();
    void initializeTimestep();
    void initializeViscosity();
    void initializeTemperature();
    void initializeTemperatureBoundary();
    void initializeVelocityBoundary();

    // Main loop equations.
    void updateForcingTerms();
    void updateViscosity();
    void solveStokes();
    void recalculateTimestep();
    void solveAdvectionDiffusion();
    bool advanceTimestep();

    // Advection methods
    void upwindMethod();
    void laxWendroff();
    void frommMethod();
    void frommVanLeer();

    // Diffusion methods
    void forwardEuler();
    void backwardEuler();
    void crankNicolson();

    void outputPressure();
    void outputVelocity();
    void outputBoundaryVelocity();
    void outputForcing();
    void outputViscosity();
    void outputTemperature();
    void outputBoundaryTemperature();

    double getH();
    double getTime();
    double getEndTime();
    int getTimestepNumber();

  private:
    ParamParser& parser;
    GeometryStructure& geometry;

    string forcingModel;
    string temperatureModel;
    string viscosityModel;
    string boundaryModel;
    string outputFile;

    int M;
    int N;

    double cfl;

    double time;
    double endTime;
    double deltaT;
    int timestepNumber;
    int endStep;

    double xExtent;
    double yExtent;
    double h;

    double viscosity;
    double diffusivity;
};
