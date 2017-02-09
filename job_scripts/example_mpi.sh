#!/bin/bash -x
#
#PJM --rsc-list "node=2x2x2"
#PJM --rsc-list "elapse=1:00:00"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./boidsim_mpi %r:./"
#PJM --stgin "rank=* ./settings/*.ini %r:./"
#PJM --stgout "rank=0 ./*.ptcl /data/hp160264/k03378/"
#PJM --mpi "rank-map-bynode=XYZ"
#PJM --mpi "use-rankdir"
#PJM -s
#PJM --name example_mpi
#

. /work/system/Env_base_1.2.0-20-1
export LD_LIBRARY_PATH=/opt/rist/boost-1.53.0/lib:$LD_LIBRARY_PATH

export TIME_STEP=1800
export AGENT_NUM=8000
export FIELD_SIZE=2.0
export DATA_FILE_NAME=data_${PJM_JOBID}.ptcl

mpiexec ./boidsim_mpi -s example.ini -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
#mpiexec -mca mpi_print_stats 1 ./boidsim_mpi -s example.ini -T ${TIME_STEP} -N ${AGENT_NUM} -F ${FIELD_SIZE} -o ${DATA_FILE_NAME}
