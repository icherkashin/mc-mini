#pragma once

#include "parser/parser.h"

/*! ProblemStructure holds all of the specific details and routines for the
 *  current problem which would be out-of-place in the GeometryStructure.
 *  These details include the current time, the x and y-extents of the problem
 *  domain, and certain physical constants such as diffusivity and viscosity.
 */

class ProblemStructure {
  public:
    /** ProblemStructure constructor.
     *  Construct the ProblemStructure object with a reference to the 
     *  ParamParser and GeometryStructure objects.
     */
    ProblemStructure (ParamParser& pp, GeometryStructure& gs);
    ~ProblemStructure();

    /** Initialize the problem. 
     *  This includes initializing the timestep, setting the initial viscosity,
     *  impose the initial temperature condition and temperature and velocity
     *  boundary conditions.
     */
    void initializeProblem();
    /** Initialize the timestep based on the 
     *
     */
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
          r /= (uf - u);
        else 
          if (r <= 0)
            r = 0;
          else
            r = UINT_MAX;
        
        return max (0.0, min (1.0, r));
    }

    double superbee (double ub, double u, double uf) {
        double r = (u - ub);

        if ((uf - u) != 0) 
          r /= (uf - u);
        else
          if (r <= 0)
            r = 0;
          else
            r = UINT_MAX;

        return max (max (0.0, min (2 * r, 1.0)), min (r, 2.0)) / 2; 
    }

    double vanLeer (double ub, double u, double uf) {
        double r = (u - ub);
        if ((uf - u) != 0)
          r /= (uf - u);
        else
          if (r <= 0)
            r = 0;
          else
            r = UINT_MAX;

        return (r + abs (r)) / (1 + abs (r)) / 2;
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
