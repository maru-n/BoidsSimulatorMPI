#!/bin/bash -x
#
#PJM --rsc-list "node={{node}}"
#PJM --rsc-list "elapse={{elapse}}"
#PJM --stg-transfiles all
#PJM --stgin "./boidsim ./"
#PJM --stgin "../settings/*.ini ./"
#PJM --stgout "./*.ptcl /data/hp160264/k03378{{data_safix}}/"
#PJM -s
#PJM --name {{name}}
#



. /work/system/Env_base_1.2.0-20-1
export LD_LIBRARY_PATH=/opt/rist/boost-1.53.0/lib:$LD_LIBRARY_PATH

export TIME_STEP={{time_step}}
export AGENT_NUM={{agent_num}}
export FIELD_SIZE={{field_size}}
export DATA_FILE_NAME={{name}}_${PJM_JOBID}.ptcl

./boidsim -s {{setting}} -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
