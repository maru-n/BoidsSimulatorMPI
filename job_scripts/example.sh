#!/bin/bash -x
#
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=1:00:00"
#PJM --stg-transfiles all
#PJM --stgin "./boidsim ./"
#PJM --stgin "./settings/*.ini ./"
#PJM --stgout "./*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name example
#



. /work/system/Env_base_1.2.0-20-1
export LD_LIBRARY_PATH=/opt/rist/boost-1.53.0/lib:$LD_LIBRARY_PATH

export TIME_STEP=1800
export AGENT_NUM=1000
export FIELD_SIZE=1.0
export DATA_FILE_NAME=data_${PJM_JOBID}.ptcl

./boidsim -s example.ini -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
