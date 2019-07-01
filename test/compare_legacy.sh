#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export setting_file="../../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_name_org=`mktemp`
export data_name_legacy=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim boidsim_legacy || exit -1

echo -e "calculating normal version to" ${data_name_org}
./boidsim -s ${setting_file} -o ${data_name_org} -T $T > /dev/null || exit -1
dd if=${data_name_org} bs=1 skip=50 of=${data_name} 2> /dev/null

echo -e "calculating legacy version to" ${data_name_legacy}
./boidsim_legacy ${data_name_legacy} $N $T > /dev/null || exit -1

cmp ${data_name} ${data_name_legacy} 0 0 > /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
