#!/bin/bash -x
#
#PJM --rsc-list "node=8"
#PJM --rsc-list "elapse=1:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./MassiveSwarmSimCUI ./"
#PJM --stgout "./data*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name testjob
#

. /work/system/Env_base

MassiveSwarmSimCUI data_<>.ptcl <N> <T>
