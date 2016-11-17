cmake_minimum_required(VERSION 2.6)
project(MassiveSwarm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(
        /usr/local/include/eigen3/
)

add_executable(MassiveSwarmSimCUI
        src/main.cpp
        src/boid.cpp
        src/vector3D.cpp
        )

add_executable(MassiveSwarmSimCUI_legacy
        src_legacy/main.cpp
        src_legacy/boid.cpp
        src_legacy/vector3D.cpp
        )

