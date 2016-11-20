#!/bin/bash -x
#
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=4:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./MassiveSwarmSimCUI ./"
#PJM --stgin "./settings/example.ini ./"
#PJM --stgout "./data*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name testjob
#

. /work/system/Env_base

export AGENT_NUM=10000
export TIME_STEP=1800

MassiveSwarmSimCUI data.ptcl ${AGENT_NUM} ${TIME_STEP} *.ini