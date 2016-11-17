#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export data_name=data.ptcl
export data_name_legacy=data_legay.ptcl

if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi

cd ${build_dir}
cmake ..
make

echo -e "calculating current version."
./MassiveSwarmSimCUI ${data_name} $N $T > /dev/null

echo -e "calculating legacy version."
./MassiveSwarmSimCUI_legacy ${data_name_legacy} $N $T > /dev/null

cmp ${data_name} ${data_name_legacy} > /dev/null

if [ $? = 0 ]; then
    echo -e "\033[0;32mSucceeded\033[0;39m"
else
    echo -e "\033[0;31mFailed\033[0;39m"
fi

#rm -rf build
