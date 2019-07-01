#!/usr/bin/env bash

if [[ $1 == '-h' ]]; then
    echo Usage: $0 [target_dir]
    exit
fi

target_dir=$1

if [[ ${target_dir} == '' ]]; then
    cd `dirname $0`
    cd ../
    target_dir='./data'
fi

rsync --exclude clean --exclude old --exclude masumori --exclude trash --exclude backup* --copy-links -ahv -P -e ssh klogin:~/data/* ${target_dir}
