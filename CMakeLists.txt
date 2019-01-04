cmake_minimum_required(VERSION 3.8)
project(heat1D)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
find_package(MPI)
include_directories(
src
${MPI_CXX_INCLUDE_PATH}
)
set(SOURCE_FILES
src/se2.cc
)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})

