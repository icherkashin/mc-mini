# Specify .cpp files
set(TARGET_NAME ${PROJECT_NAME})
set(SRC 
  geometry/geometry.cpp

  matrixForms/denseForms.cpp
  matrixForms/sparseForms.cpp

  output/output.cpp

  parser/paramTree.cpp
  parser/parser.cpp

  problem/advection.cpp 
  problem/diffusion.cpp
  problem/initialization.cpp
  problem/problem.cpp
  problem/solveRoutines.cpp)

# Build an executable from main.cpp and all the specified source files
add_executable(${PROJECT_NAME} main.cpp ${SRC})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})

# Build a library from all specified source files
# This is required for using Google Test
add_library(${PROJECT_NAME}-lib ${SRC})
target_link_libraries(${PROJECT_NAME}-lib ${LIBRARIES})
set_target_properties(${PROJECT_NAME}-lib PROPERTIES OUPTUT_NAME ${PROJECT_NAME})
