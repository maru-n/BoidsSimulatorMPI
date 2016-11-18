cmake_minimum_required(VERSION 2.6)

project(MassiveSwarm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package( OpenMP )
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(MassiveSwarmSimCUI
        src/main.cpp
        src/boid.cpp
        src/vector3D.cpp
        src/boid_simulation.cpp
        src/boid_simulation.h)

add_executable(MassiveSwarmSimCUI_legacy
        src_legacy/main.cpp
        src_legacy/boid.cpp
        src_legacy/vector3D.cpp
        )

