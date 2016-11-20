#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export data_name=test_data.ptcl
export data_name_legacy=test_data_legay.ptcl
export setting_file="../settings/test.ini"
export CMAKE_OPTIONS="-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++"

if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi

cd ${build_dir}
cmake ${CMAKE_OPTIONS} ..
if [ $? != 0 ]; then
    echo build process failed 1>&2
    exit -1
fi

make
if [ $? != 0 ]; then
    echo build process failed 1>&2
    exit -1
fi

echo -e "calculating current version."
#./MassiveSwarmSimCUI ${data_name} $N $T ${setting_file} 1 0.1 1 0.0001 0.1 1 0.1 0.1 1 0.1 0.003 0.001 > /dev/null
./MassiveSwarmSimCUI ${data_name} $N $T test 12345 1 0.1 1 0.0001 0.1 1 0.1 0.1 1 0.1 0.003 0.001 > /dev/null
if [ $? != 0 ]; then
    echo failed 1>&2
    exit -1
fi

echo -e "calculating legacy version."
./MassiveSwarmSimCUI_legacy ${data_name_legacy} $N $T > /dev/null
if [ $? != 0 ]; then
    echo failed 1>&2
    exit -1
fi

cmp ${data_name} ${data_name_legacy} 0 0> /dev/null

if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

#cd ..
#rm -rf ${build_dir}
exit 0
