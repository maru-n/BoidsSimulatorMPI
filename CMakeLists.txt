cmake_minimum_required(VERSION 2.6)
project(MassiveSwarm)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}")

# for Debugger
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# compiler options on K
if (DEFINED ENV{Env_base})
    #set(K_C_CXX_FLAGS "-Xg -Kfast -Kparallel -Koptmsg=2 -V -Nsrc -Nsta")
    #set(K_C_CXX_FLAGS "-Xg -Kfast -Kparallel -Koptmsg=1 -w")
    set(K_C_CXX_FLAGS "-Xg -Kfast -Kparallel -w")
    #set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -L/opt/rist/boost-1.53.0/lib -I/opt/rist/boost-1.53.0/include")  # boost1.53.0
    set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -L/opt/klocal/boost-1.65.1-fujitsu/lib -I/opt/klocal/boost-1.65.1-fujitsu/include")  # boost1.65.1
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${K_C_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${K_C_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L/opt/klocal/boost-1.65.1-fujitsu/lib -I/opt/klocal/boost-1.65.1-fujitsu/include")
endif()
# for FX100
#set(K_C_CXX_FLAGS "-Xg -Kfast -Kparallel -w")
#set(K_C_CXX_FLAGS "${K_C_CXX_FLAGS} -L/opt/klocal/boost-1.65.1-fujitsu/lib -I/opt/klocal/boost-1.65.1-fujitsu/include")  # boost1.65.1
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${K_C_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${K_C_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L/opt/klocal/boost-1.65.1-fujitsu/lib -I/opt/klocal/boost-1.65.1-fujitsu/include")

find_package(Boost REQUIRED program_options)
if(Boost_FOUND)
else()
    message(STATUS "!!! BOOST NOT FOUND !!!")
endif()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

execute_process(COMMAND brew --prefix libomp
                OUTPUT_VARIABLE OpenMP_HOME
                OUTPUT_STRIP_TRAILING_WHITESPACE)  
message(STATUS "OpenMP Root : ${OpenMP_HOME}")
# OpenMPのインストールチェック
set(OpenMP_C_LIB_NAMES "omp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "${OpenMP_HOME}/lib/")
# キャッシュ変数なのでFORCEオプションをつける
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -Wno-unused-command-line-argument -I${OpenMP_HOME}/include -lomp -L${OpenMP_omp_LIBRARY}" CACHE STRING "" FORCE) 
# 同じくキャッシュ変数なのでFORCEオプションをつける
set(OpenMP_C_FLAGS "-fopenmp -Wno-unused-command-line-argument -I${OpenMP_HOME}/include -lomp -L${OpenMP_omp_LIBRARY}" CACHE STRING "" FORCE)


find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MPI)
if (MPI_CXX_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    add_executable(boidsim_mpi
            src/main.cpp
            src/args.cpp
            src/boid.cpp
            src/vector3D.cpp
            src/boid_simulation.cpp
            src/boid_simulation_mpi.cpp)
    target_link_libraries(boidsim_mpi ${Boost_LIBRARIES})
    target_compile_definitions(boidsim_mpi PRIVATE _MPI=1)
    install(TARGETS boidsim_mpi DESTINATION bin)

    add_executable(boidsim_mpi_nogrid
            src/main.cpp
            src/args.cpp
            src/boid.cpp
            src/vector3D.cpp
            src/boid_simulation.cpp
            src/boid_simulation_mpi_nogrid.cpp)
    target_link_libraries(boidsim_mpi_nogrid ${Boost_LIBRARIES})            
    target_compile_definitions(boidsim_mpi_nogrid PRIVATE _MPI=1)
    install(TARGETS boidsim_mpi_nogrid DESTINATION bin)
endif()

add_executable(boidsim
        src/main.cpp
        src/args.cpp
        src/boid.cpp
        src/vector3D.cpp
        src/boid_simulation.cpp
        src/dtype.h)
target_link_libraries(boidsim ${Boost_LIBRARIES})
install(TARGETS boidsim DESTINATION bin)

add_executable(boidsim_grid
        src/main.cpp
        src/args.cpp
        src/boid.cpp
        src/vector3D.cpp
        src/boid_simulation_grid.cpp
        src/dtype.h)
target_link_libraries(boidsim_grid ${Boost_LIBRARIES})
install(TARGETS boidsim_grid DESTINATION bin)


add_executable(boidsim_legacy
        src/legacy/main.cpp
        src/legacy/boid.cpp
        #src/legacy/vector3D.cpp
        src/vector3D.cpp
        )

