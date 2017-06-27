#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export setting_file="../../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export data_name=`mktemp`
export data_name_split=`mktemp`
export data_name_split_first=`mktemp`
export data_name_split_last=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim || exit -1

echo -e "calculating normal version..."
./boidsim -s ${setting_file} -T $T -o ${data_name} > /dev/null || exit -1

echo -e "calculating split version..."
t_first=$(($T/2))
t_last=$(( $T-$t_first ))
./boidsim -s ${setting_file} -T $t_first -o ${data_name_split_first} > /dev/null || exit -1
./boidsim -s ${setting_file} -T $t_last  -i ${data_name_split_first} -o ${data_name_split_last} > /dev/null || exit -1
echo $data_name_split_first
echo $data_name_split_last
echo $data_name_split

cmp ${data_name} ${data_name_split} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
