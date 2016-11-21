#cmake_minimum_required(VERSION 2.6)
cmake_minimum_required(VERSION 3.6)

project(MassiveSwarm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if ($ENV{Env_base} MATCHES "K-1.2.0-20-1")
    set(K_LOCAL "$ENV{HOME}/local")
    #set(K_C_CXX_FLAGS "-Kfast -Kparallel -Kopenmp -Koptmsg=2 -V -Nsrc -Nsta -L${K_LOCAL}/lib -I${K_LOCAL}/include")
    set(K_C_CXX_FLAGS   "-w -Kfast -Kparallel -Kopenmp -Koptmsg=1 -L${K_LOCAL}/lib -I${K_LOCAL}/include")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${K_C_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${K_C_CXX_FLAGS}")
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
        src/dtype.h)

if ($ENV{Env_base} MATCHES "K-1.2.0-20-1")
    #find_package( MPI )
    #if (MPI_FOUND)
    add_executable(MassiveSwarmSimCUI_MPI
            src/main.cpp
            src/boid.cpp
            src/vector3D.cpp
            src/boid_simulation.cpp
            src/dtype.h)
    target_compile_options(MassiveSwarmSimCUI_MPI PRIVATE "-DMPI_ENABLE")
    #endif()
endif()

add_executable(MassiveSwarmSimCUI_legacy
        src/legacy/main.cpp
        src/legacy/boid.cpp
        src/legacy/vector3D.cpp
        )

