#!/bin/bash -x
#
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=1:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./boidsim ./"
#PJM --stgin "./settings/*.ini ./"
#PJM --stgout "./data*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name example
#

. /work/system/Env_base

export AGENT_NUM=1000
export TIME_STEP=1800
export DATA_FILE_NAME=${PJM_JOBID}.ptcl

./boidsim ${DATA_FILE_NAME} ${AGENT_NUM} ${TIME_STEP} settings.ini
