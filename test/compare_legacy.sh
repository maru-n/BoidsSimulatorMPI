#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export setting_file="../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_name_mpi=`mktemp`
export data_name_legacy=`mktemp`
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} .. || exit -1
make boidsim boidsim_mpi boidsim_legacy || exit -1

echo -e "calculating normal version."
./boidsim -s ${setting_file} -o ${data_name} -T $T > /dev/null || exit -1

echo -e "calculating MPI version."
mpiexec -n 8 ./boidsim_mpi -s ${setting_file} -o ${data_name_mpi} -T $T > /dev/null || exit -1

echo -e "calculating legacy version."
./boidsim_legacy ${data_name_legacy} $N $T > /dev/null || exit -1

cmp ${data_name} ${data_name_legacy} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "normal:\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "normal:\033[0;31mTest Failed\033[0;39m"
fi

cmp ${data_name_mpi} ${data_name_legacy} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "MPI:   \033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "MPI:   \033[0;31mTest Failed\033[0;39m"
fi

exit 0
