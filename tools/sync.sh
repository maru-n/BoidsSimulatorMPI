#!/usr/bin/env bash

cd `dirname $0`
cd ../

rsync --exclude build --exclude data --exclude bin --exclude test --copy-links -ahv -e ssh . klogin:~/MassiveSwarm