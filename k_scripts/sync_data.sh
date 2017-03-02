#!/usr/bin/env bash

target_dir=$1

rsync --exclude old --exclude masumori --exclude trash --exclude backup* --copy-links -ahv -P -e ssh klogin:~/data/* ${target_dir}

