cmake_minimum_required(VERSION 3.10)
project(TDSE_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Include JSON
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find ArrayFire (will require system installation or CI setup)
find_package(ArrayFire REQUIRED)

# Source files
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(tdse ${SOURCES})

# Use CUDA or CPU backend (user can change this if needed)
target_link_libraries(tdse PRIVATE ArrayFire::afcuda)
