# Specify .cpp files
set(TARGET_NAME ${PROJECT_NAME})
set(SRC advection.cpp denseForms.cpp geometry.cpp paramTree.cpp problem.cpp sparseForms.cpp diffusion.cpp initialization.cpp output.cpp parser.cpp solveRoutines.cpp)

# Build an executable from main.cpp and all the specified source files
add_executable(${PROJECT_NAME} main.cpp ${SRC})
target_link_libraries(${TARGET_NAME} ${LIBRARIES})

# Build a library from all specified source files
# This is required for using Google Test
add_library(${PROJECT_NAME}-lib ${SRC})
target_link_libraries(${PROJECT_NAME}-lib ${LIBRARIES})
set_target_properties(${PROJECT_NAME}-lib PROPERTIES OUPTUT_NAME ${PROJECT_NAME})
