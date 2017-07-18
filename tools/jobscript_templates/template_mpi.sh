#!/bin/bash -x
#
#PJM --rsc-list "node={{node}}"
#PJM --rsc-list "elapse={{elapse}}"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./boidsim_mpi %r:./"
#PJM --stgin "rank=* ../settings/*.ini %r:./"
#PJM --stgout "rank=* ./*.ptcl* /data/hp160264/k03378{{data_safix}}/"
#PJM --mpi "rank-map-bynode=XYZ"
#PJM --mpi "use-rankdir"
#PJM -s
#PJM --name {{name}}
#

. /work/system/Env_base_1.2.0-20-1
export LD_LIBRARY_PATH=/opt/rist/boost-1.53.0/lib:$LD_LIBRARY_PATH

export TIME_STEP={{time_step}}
export AGENT_NUM={{agent_num}}
export FIELD_SIZE={{field_size}}
export DATA_FILE_NAME={{name}}_${PJM_JOBID}.ptcl

mpiexec ./boidsim_mpi -s {{setting}} -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
#mpiexec -mca mpi_print_stats 1 ./boidsim_mpi -s example.ini -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
