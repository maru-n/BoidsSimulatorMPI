#!/usr/bin/env bash

export N=1000
export T=50
export build_dir=build
export setting_file="../../settings/mototake.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_name_mpi=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim boidsim_mpi || exit -1

echo -e "calculating normal version... =>" ${data_name}
time ./boidsim -s ${setting_file} -o ${data_name} -T $T > /dev/null || exit -1

echo -e "calculating MPI version... =>" ${data_name_mpi}
time mpiexec -n 8 ./boidsim_mpi -s ${setting_file} -o ${data_name_mpi} -T $T > /dev/null || exit -1

exit 0
