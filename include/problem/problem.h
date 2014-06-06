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

    // Flux Limiters
    double minmod (double ub, double u, double uf) {
        double r = (u - ub);
        if ((uf - u) != 0)
          (r /= (uf - u));
        else if (r != 0)
          r = INT_MAX;
        // Using std::max and std::min from <algorithm> allows us to use more
        // than two arguments, but requires us to do all kinds of funky stuff
        // to make types unambiguous.
        return std::max (0.0, std::min (1.0, r));
    }

    double superbee (double ub, double u, double uf) {
        double r = (u - ub);
        if ((uf - u) != 0) 
          r /= (uf - u);
        else if (r != 0)
          r = INT_MAX;
        return std::max (std::max(0.0, std::min (2 * r, 1.0)), std::min (r, 2.0));
    }

    double vanLeer (double ub, double u, double uf) {
        double r = (u - ub);
        if ((uf - u) != 0)
          r /= (uf - u);
        else if (r != 0)
          r = INT_MAX;
        return (r + std::abs (r)) / (1 + std::abs (r));
    }

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
    string advectionMethod;
    string fluxLimiter;
    string diffusionMethod;
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
