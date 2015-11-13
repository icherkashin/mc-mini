#pragma once

#include <fstream>

#include "geometry/geometry.h"
#include "problem/problem.h"
#include "parser/parser.h"

class OutputStructure {
  public:
    OutputStructure (ParamParser&       pp,
                     GeometryStructure& gs,
                     ProblemStructure&  ps);

    ~OutputStructure();

    void outputData (const int timestep = 0);

    void writeHDF5File (const int timestep = 0);
  
  private:
    ParamParser&       parser;
    GeometryStructure& geometry;
    ProblemStructure&  problem;

    int M;
    int N;

    double dx;

    string outputFormat;
    string outputPath;
    string outputFilename;

    std::ofstream problemXdmfFile;
};
