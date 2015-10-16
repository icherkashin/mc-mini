#include <gtest/gtest.h>

#include "matrixForms/denseForms.h"

TEST(DenseForms, makeLaplacianXBlock) {
  // instantiate test viscosity data
  double *viscosity_data = new double[16];
  for (int i = 0; i < 16; ++i) {
    viscosity_data[i] = 1.0;
  }

  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(6, 6);
  DenseForms::makeLaplacianXBlock(actual_matrix, 3, 3, 1, viscosity_data);

  // clean up after the viscosity data
  delete[] viscosity_data;

  Eigen::MatrixXd expected_matrix(6, 6);
  expected_matrix << 5,-1,-1, 0, 0, 0,
                    -1, 5, 0,-1, 0, 0,
                    -1, 0, 4,-1,-1, 0,
                     0,-1,-1, 4, 0,-1,
                     0, 0,-1, 0, 5,-1,
                     0, 0, 0,-1,-1, 5;

  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeLaplacianYBlock) {
  // instntiate test viscosity data
  double *viscosity_data = new double[16];
  for (int i = 0; i < 16; ++i) { 
    viscosity_data[i] = 1.0; 
  }

  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(6, 6);
  DenseForms::makeLaplacianYBlock(actual_matrix, 3, 3, 1, viscosity_data);

  // clean up after the viscosity data
  delete[] viscosity_data;

  Eigen::MatrixXd expected_matrix(6, 6);
  expected_matrix << 5,-1, 0,-1, 0, 0,
                    -1, 4,-1, 0,-1, 0,
                     0,-1, 5, 0, 0,-1,
                    -1, 0, 0, 5,-1, 0,
                     0,-1, 0,-1, 4,-1,
                     0, 0,-1, 0,-1, 5;

  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeGradXBlock) {
  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(6, 9);
  DenseForms::makeGradXBlock(actual_matrix, 3, 3, 1.0);

  // Create the expected matrix
  Eigen::MatrixXd expected_matrix(6, 9);
  expected_matrix << -1, 1, 0, 0, 0, 0, 0, 0, 0,
                      0,-1, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0,-1, 1, 0, 0, 0, 0,
                      0, 0, 0, 0,-1, 1, 0, 0, 0,
                      0, 0, 0, 0, 0, 0,-1, 1, 0,
                      0, 0, 0, 0, 0, 0, 0,-1, 1;

  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeGradYBlock) {
  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(6, 9);
  DenseForms::makeGradYBlock(actual_matrix, 3, 3, 1.0);

  // Create the expected matrix
  Eigen::MatrixXd expected_matrix(6, 9);
  expected_matrix << -1, 0, 0, 1, 0, 0, 0, 0, 0,
                      0,-1, 0, 0, 1, 0, 0, 0, 0,
                      0, 0,-1, 0, 0, 1, 0, 0, 0,
                      0, 0, 0,-1, 0, 0, 1, 0, 0,
                      0, 0, 0, 0,-1, 0, 0, 1, 0,
                      0, 0, 0, 0, 0,-1, 0, 0, 1;
  
  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeDivXBlock) {
  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(9, 6);
  DenseForms::makeDivXBlock(actual_matrix, 3, 3, 1.0);

  Eigen::MatrixXd expected_matrix(9, 6);
  expected_matrix << 1, 0, 0, 0, 0, 0,
                    -1, 1, 0, 0, 0, 0,
                     0,-1, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0,
                     0, 0,-1, 1, 0, 0,
                     0, 0, 0,-1, 0, 0,
                     0, 0, 0, 0, 1, 0,
                     0, 0, 0, 0,-1, 1,
                     0, 0, 0, 0, 0,-1;

  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeDivYBlock) {
  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(9, 6);
  DenseForms::makeDivYBlock(actual_matrix, 3, 3, 1.0);

  Eigen::MatrixXd expected_matrix(9, 6);
  expected_matrix << 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0,
                    -1, 0, 0, 1, 0, 0,
                     0,-1, 0, 0, 1, 0,
                     0, 0,-1, 0, 0, 1,
                     0, 0, 0,-1, 0, 0,
                     0, 0, 0, 0,-1, 0,
                     0, 0, 0, 0, 0,-1;
  
  ASSERT_EQ(expected_matrix, actual_matrix);
}

TEST(DenseForms, makeForcingMatrix) {
  Eigen::MatrixXd actual_matrix = Eigen::MatrixXd::Zero(8, 4);
  DenseForms::makeForcingMatrix(actual_matrix, 2, 2);

  // Generate the expected matrix
  Eigen::MatrixXd expected_matrix(8, 4);
  expected_matrix << 1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0; 

  ASSERT_EQ(expected_matrix, actual_matrix);
}
