#!/usr/bin/env bash

export N=1000
export T=30
export build_dir=build
export setting_file="../../settings/test.ini"
export CMAKE_OPTIONS="-D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++"

export TOOLS_PATH="../../tools/"

export data_name=`mktemp`

export tmp=`mktemp`
export data_name_split_first=`mktemp`
export data_name_split_first_v=`mktemp`
export data_name_secong_init=`mktemp`
export data_name_split_second=`mktemp`
export data_name_split_merge=`mktemp`

echo -e "build..."
if [ ! -e ${build_dir} ]; then mkdir ${build_dir} ; fi
cd ${build_dir}
cmake ${CMAKE_OPTIONS} ../../ || exit -1
make boidsim || exit -1

echo -e "calculating normal version..."
./boidsim -s ${setting_file} -T $T -o ${data_name} > /dev/null || exit -1


t_first=$(($T/2))
t_last=$(( $T-$t_first+1 ))
echo -e "calculating split version...(1/2)"
./boidsim -s ${setting_file} -T $t_first -o ${data_name_split_first} --velocity-output ${data_name_split_first_v} > /dev/null || exit -1

echo -e "generate next initialization file..."
python ${TOOLS_PATH}generate_initialization_file.py ${data_name_split_first} ${data_name_split_first_v} $(( t_first-1 )) ${data_name_secong_init} || exit -1

echo -e "calculating split version...(2/2)"
./boidsim -s ${setting_file} -T $t_last -o ${data_name_split_second} --initialization ${data_name_secong_init} > /dev/null || exit -1


echo -e "merge split data..."
python ${TOOLS_PATH}split_datafile.py ${data_name_split_second} ${tmp} 1 $((t_last-1))
python ${TOOLS_PATH}merge_datafiles.py ${data_name_split_first} ${tmp} ${data_name_split_merge}


cmp ${data_name} ${data_name_split_merge} 0 0> /dev/null
if [ $? = 0 ]; then
    echo -e "\033[0;32mTest Succeeded\033[0;39m"
else
    echo -e "\033[0;31mTest Failed\033[0;39m"
fi

exit 0
