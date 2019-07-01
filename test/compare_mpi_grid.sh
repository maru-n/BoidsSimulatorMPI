#!/usr/bin/env bash

export N=1000
export T=50
export build_dir=build
export setting_file="../../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_name_grid=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim_mpi_nogrid boidsim_mpi || exit -1

echo -e "calculating normal MPI version to " ${data_name}
mpiexec --hostfile ../../hostfile -n 8 ./boidsim_mpi_nogrid -s ${setting_file} -o ${data_name} -T $T -N $N > /dev/null || exit -1

echo -e "calculating grid optimized MPI version to" ${data_name_grid}
mpiexec --hostfile ../../hostfile -n 8 ./boidsim_mpi -s ${setting_file} -o ${data_name_grid} -T $T -N $N > /dev/null || exit -1

cmp ${data_name} ${data_name_grid} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
