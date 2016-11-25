#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export data_name=data_test.ptcl
export data_name_legacy=data_test_legay.ptcl
export setting_file="../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++"

if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} .. || exit -1
make boidsim boidsim_legacy || exit -1

echo -e "calculating current version."
./boidsim ${data_name} $N $T ${setting_file} > /dev/null || exit -1

echo -e "calculating legacy version."
./boidsim_legacy ${data_name_legacy} $N $T > /dev/null || exit -1

cmp ${data_name} ${data_name_legacy} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
