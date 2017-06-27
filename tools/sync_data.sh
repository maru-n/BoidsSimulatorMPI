#!/usr/bin/env bash

if [[ $1 == '-h' ]]; then
    echo Usage: $0 [target_dir]
    exit
fi

target_dir=$1

if [[ ${target_dir} == '' ]]; then
    target_dir='./'
fi

rsync --exclude old --exclude masumori --exclude trash --exclude backup* --copy-links -ahv -P -e ssh klogin:~/data/* ${target_dir}
