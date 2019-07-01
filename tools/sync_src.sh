#!/usr/bin/env bash

cd `dirname $0`
cd ../

#rsync --exclude testdata --exclude notebook --exclude build --exclude data --exclude bin --exclude test --copy-links -ahv -e ssh . klogin:~/MassiveSwarm
rsync --copy-links -ahv -e ssh CMakeLists.txt klogin:~/MassiveSwarm
rsync --copy-links -ahv -e ssh ./src klogin:~/MassiveSwarm
rsync --copy-links -ahv -e ssh ./job_scripts klogin:~/MassiveSwarm
rsync --copy-links -ahv -e ssh ./settings klogin:~/MassiveSwarm
