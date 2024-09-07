cmake_minimum_required(VERSION 3.14)

# Project name and version
project(SimpleMPI VERSION 1.0
	DESCRIPTION "Parallelization of an MD code"
		  LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Add the executable
add_executable(mpimd src/main.cpp)

# Find MPI
find_package(MPI REQUIRED)
target_link_libraries(mpimd PUBLIC MPI::MPI_CXX)

