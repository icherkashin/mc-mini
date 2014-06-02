#include <iostream>
#include <cassert>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "matrixForms/sparseForms.h"
#include "geometry/geometry.h"
#include "problem/problem.h"
#include "debug.h"

using namespace Eigen;
using namespace std;

// Forward Euler diffusion method. Unstable but fairly efficient.
void ProblemStructure::forwardEuler() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  double mu = deltaT * diffusivity / (h * h);

  SparseMatrix<double> rhs;
  SparseMatrix<double> rhsBoundary;
  rhs.resize (M * N, M * N);
  rhsBoundary.resize (M * N, 2 * N);

  vector<Triplet<double> > tripletList;
  tripletList.reserve (5 * M * N);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; ++j) {
      if ((j == 0) || (j == (N - 1)))
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 - 3 * mu)); 
      else
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 - 4 * mu));
      if (j > 0) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j - 1), mu));
      if (j < (N - 1)) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j + 1), mu));
      if (i > 0)
        tripletList.push_back (Triplet<double> (i * N + j, (i - 1) * N + j, mu));
      if (i < (M - 1))
        tripletList.push_back (Triplet<double> (i * N + j, (i + 1) * N + j, mu));
    }

  rhs.setFromTriplets (tripletList.begin(), tripletList.end());
  rhs.makeCompressed();

  tripletList.clear();
  tripletList.reserve (2 * N);

  for (int j = 0; j < N; ++j) {
    tripletList.push_back (Triplet<double> (j,               j,     mu));
    tripletList.push_back (Triplet<double> ((M - 1) * N + j, N + j, mu));
  }


  rhsBoundary.setFromTriplets (tripletList.begin(), tripletList.end());
  rhsBoundary.makeCompressed();

  VectorXd temporaryVector = temperatureVector;
  temperatureVector = rhs * temporaryVector + rhsBoundary * temperatureBoundaryVector;
}

// Backward Euler Diffusion method. Stable but inefficient.
void ProblemStructure::backwardEuler() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * N);

  double mu = deltaT * diffusivity / (h * h);

  SparseMatrix<double> lhs;
  SparseMatrix<double> rhsBoundary;
  lhs.resize (M * N, M * N);
  rhsBoundary.resize (M * N, 2 * N);

  vector<Triplet<double> > tripletList;
  tripletList.reserve (5 * M * N);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; ++j) {
      if ((j == 0) || (j == (N - 1)))
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 + 3 * mu)); 
      else
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 + 4 * mu));
      if (j > 0) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j - 1), -mu));
      if (j < (N - 1)) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j + 1), -mu));
      if (i > 0)
        tripletList.push_back (Triplet<double> (i * N + j, (i - 1) * N + j, -mu));
      if (i < (M - 1))
        tripletList.push_back (Triplet<double> (i * N + j, (i + 1) * N + j, -mu));
    }

  lhs.setFromTriplets (tripletList.begin(), tripletList.end());
  lhs.makeCompressed();

  tripletList.clear();
  tripletList.reserve (2 * N);

  for (int j = 0; j < N; ++j) {
    tripletList.push_back (Triplet<double> (j,               j,     mu));
    tripletList.push_back (Triplet<double> ((M - 1) * N + j, N + j, mu));
  }

  rhsBoundary.setFromTriplets (tripletList.begin(), tripletList.end());
  rhsBoundary.makeCompressed();

  VectorXd temporaryVector = temperatureVector;

  #ifdef DEBUG
    cout << "<Backward Euler " << lhs.rows() << "x" << lhs.cols() << " LHS Matrix generated>" << endl;
    cout << "<Backward Euler " << rhsBoundary.rows() << "x" << rhsBoundary.cols() << " RHS Boundary Matrix generated>" << endl;
    cout << "<Temperature Boundary Vector has "<< temperatureBoundaryVector.rows() << " elements>" << endl << endl;
  #endif
  
  SimplicialLLT<SparseMatrix<double> > solver;
  solver.compute (lhs);
  temperatureVector =  solver.solve(temporaryVector + rhsBoundary * temperatureBoundaryVector);
}

void ProblemStructure::crankNicolson() {
  Map<VectorXd> temperatureVector (geometry.getTemperatureData(), M * N);
  Map<VectorXd> temperatureBoundaryVector (geometry.getTemperatureBoundaryData(), 2 * M);

  double mu = deltaT * diffusivity / (2 * h * h);

  SparseMatrix<double> rhs (M * N, M * N);
  SparseMatrix<double> rhsBoundary (M * N, 2 * N);

  vector<Triplet<double> > tripletList;
  tripletList.reserve (5 * M * N);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; ++j) {
      if ((j == 0) || (j == (N - 1)))
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 - 3 * mu)); 
      else
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 - 4 * mu));
      if (j > 0) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j - 1), mu));
      if (j < (N - 1)) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j + 1), mu));
      if (i > 0)
        tripletList.push_back (Triplet<double> (i * N + j, (i - 1) * N + j, mu));
      if (i < (M - 1))
        tripletList.push_back (Triplet<double> (i * N + j, (i + 1) * N + j, mu));
    }

  rhs.setFromTriplets (tripletList.begin(), tripletList.end());
  rhs.makeCompressed();

  tripletList.clear();
  tripletList.reserve (2 * N);

  for (int j = 0; j < N; ++j) {
    tripletList.push_back (Triplet<double> (j,               j,     mu));
    tripletList.push_back (Triplet<double> ((M - 1) * N + j, N + j, mu));
  }

  rhsBoundary.setFromTriplets (tripletList.begin(), tripletList.end());
  rhsBoundary.makeCompressed();

  VectorXd temporaryVector = rhs * temperatureVector + rhsBoundary * temperatureBoundaryVector;
  
  SparseMatrix<double> lhs;
  lhs.resize (M * N, M * N);
  rhsBoundary.resize (M * N, 2 * M);

  tripletList.clear();
  tripletList.reserve (5 * M * N);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; ++j) {
      if ((j == 0) || (j == (N - 1)))
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 + 3 * mu)); 
      else
        tripletList.push_back (Triplet<double> (i * N + j, i * N + j, 1 + 4 * mu));
      if (j > 0) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j - 1), -mu));
      if (j < (N - 1)) 
        tripletList.push_back (Triplet<double> (i * N + j, i * N + (j + 1), -mu));
      if (i > 0)
        tripletList.push_back (Triplet<double> (i * N + j, (i - 1) * N + j, -mu));
      if (i < (M - 1))
        tripletList.push_back (Triplet<double> (i * N + j, (i + 1) * N + j, -mu));
    }

  lhs.setFromTriplets (tripletList.begin(), tripletList.end());
  lhs.makeCompressed();
  
  SimplicialLLT<SparseMatrix<double> > solver;
  solver.compute (lhs);
  temperatureVector =  solver.solve(temporaryVector + rhsBoundary * temperatureBoundaryVector);
}
