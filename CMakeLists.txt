cmake_minimum_required(VERSION 3.6)

project(MassiveSwarm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lboost_program_options")

# compiler options on K
if ($ENV{Env_base} MATCHES "K-1.2.0-20-1")
    set(K_C_CXX_FLAGS "-Xg -D_K_ -L$ENV{HOME}/local/lib -I$ENV{HOME}/local/include")
    #set(K_C_CXX_FLAGS "-Xg -D_K_ -L/usr/lib -I/usr/include")
    #set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -Kfast -Kparallel -Koptmsg=2 -V -Nsrc -Nsta")
    #set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -Kfast -Kparallel -Koptmsg=1 -w")
    set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -Kfast -Kparallel -w")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${K_C_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${K_C_CXX_FLAGS}")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MPI)
if (MPI_CXX_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    add_executable(boidsim_mpi
            src/main.cpp
            src/boid.cpp
            src/vector3D.cpp
            src/boid_simulation.cpp
            src/boid_simulation_mpi.cpp)
    target_compile_options(boidsim_mpi PRIVATE "-D_MPI")
endif()

add_executable(boidsim
        src/main.cpp
        src/boid.cpp
        src/vector3D.cpp
        src/boid_simulation.cpp
        src/dtype.h)

add_executable(boidsim_legacy
        src/legacy/main.cpp
        src/legacy/boid.cpp
        src/legacy/vector3D.cpp
        )
