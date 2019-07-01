#!/usr/bin/env bash

#export N=16384
export N=100
export T=30
export build_dir=build
export setting_file="../../settings/test.ini"

export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name_org=`mktemp`
export data_name_grid=`mktemp`

echo -e "\033[0;32mbuildding...\033[0;39m"
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim boidsim_grid || exit -1

echo -e "calculating GRID NEW version to" ${data_name_grid}
time ./boidsim_grid -s ${setting_file} -o ${data_name_grid} -T $T -N $N > /dev/null || exit -1

echo -e "calculating NORMAL version to" ${data_name_org}
time ./boidsim -s ${setting_file} -o ${data_name_org} -T $T -N $N > /dev/null || exit -1

cmp ${data_name_org} ${data_name_grid} 0 0 > /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
