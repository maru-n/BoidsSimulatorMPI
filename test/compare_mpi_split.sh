#!/usr/bin/env bash

export N=1000
export T=2
export build_dir=build
#export setting_file="../../settings/test.ini"
export setting_file="../../settings/example.ini"
#export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_dir_mpi=`mktemp -d`
export data_name_mpi=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim boidsim_mpi || exit -1

echo -e "calculating normal version to " ${data_name}
./boidsim -s ${setting_file} -o ${data_name} -T $T > /dev/null || exit -1

echo -e "calculating MPI version to" ${data_dir_mpi}
mpiexec --hostfile ../../hostfile -n 8 ./boidsim_mpi --parallel-output -s ${setting_file} -o ${data_dir_mpi}/out -T $T > /dev/null || exit -1
python ../../tools/merge_cluster_datafiles.py ${data_dir_mpi} ${data_name_mpi}

cmp ${data_name} ${data_name_mpi} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
