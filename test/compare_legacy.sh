#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export data_name=data.ptcl
export data_name_legacy=data_legay.ptcl
export CMAKE_OPTIONS="-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++"

if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi

cd ${build_dir}
cmake ${CMAKE_OPTIONS} ..
if [ $? != 0 ]; then
    echo -e "\033[0;31mbuild process failed\033[0;39m"
    exit -1
fi

make
if [ $? != 0 ]; then
    echo -e "\033[0;31mbuild process failed\033[0;39m"
    exit -1
fi

echo -e "\033[0;32mcalculating current version\033[0;39m"
./MassiveSwarmSimCUI ${data_name} $N $T >> /dev/null

echo -e "\033[0;32mcalculating legacy version\033[0;39m"
./MassiveSwarmSimCUI_legacy ${data_name_legacy} $N $T > /dev/null

cmp ${data_name} ${data_name_legacy} > /dev/null

if [ $? = 0 ]; then
    echo -e "\033[0;32mSucceeded\033[0;39m"
else
    echo -e "\033[0;31mFailed\033[0;39m"
fi

#cd ..
#rm -rf ${build_dir}
exit 0
