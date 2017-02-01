#!/bin/bash -x
#
#PJM --rsc-list "node=2x2x2"
#PJM --rsc-list "elapse=1:00:00"
#PJM --mpi "rank-map-bynode=XYZ"
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./boidsim_mpi %r:./"
#PJM --stgin "rank=* ./settings/*.ini %r:./"
#PJM --stgout "rank=0 %r:./*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name example_mpi
#

. /work/system/Env_base
export LD_LIBRARY_PATH=/opt/rist/boost-1.53.0/lib:$LD_LIBRARY_PATH

export AGENT_NUM=10000
export TIME_STEP=1800
export DATA_FILE_NAME=data_${PJM_JOBID}.ptcl

mpiexec ./boidsim_mpi ${DATA_FILE_NAME} ${AGENT_NUM} ${TIME_STEP} example.ini
#mpiexec –mca mpi_print_stats 1 ./boidsim_mpi ${DATA_FILE_NAME} ${AGENT_NUM} ${TIME_STEP} example.ini
