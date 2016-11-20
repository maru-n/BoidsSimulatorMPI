#cmake_minimum_required(VERSION 2.6)
cmake_minimum_required(VERSION 3.6)

project(MassiveSwarm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(K_LOCAL "$ENV{HOME}/local")
set(K_C_FLAGS "-Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta -L${K_LOCAL}/lib -I${K_LOCAL}/include")
set(K_CXX_FLAGS "-Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta -L${K_LOCAL}/lib -I${K_LOCAL}/include")
if (${CMAKE_C_COMPILER} MATCHES fccpx)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${K_C_FLAGS}")
endif()
if (${CMAKE_CXX_COMPILER} MATCHES FCCpx)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${K_CXX_FLAGS}")
endif()

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
        src/boid_simulation.h
        src/dtype.h)

add_executable(MassiveSwarmSimCUI_legacy
        src/legacy/main.cpp
        src/legacy/boid.cpp
        src/legacy/vector3D.cpp
        )

