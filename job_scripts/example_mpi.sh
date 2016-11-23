#!/bin/bash -x
#
#PJM --rsc-list "node=2x2x2"
#PJM --rsc-list "elapse=1:00:00"
#PJM --mpi "rank-map-bynode=XYZ"
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./MassiveSwarmSimCUI_MPI %r:./"
#PJM --stgin "rank=* ./settings/example.ini %r:./"
#PJM --stgout "rank=0 %r:./data*.ptcl /data/hp160264/k03378/"
#PJM -s
#PJM --name testjob_mpi
#

. /work/system/Env_base

export AGENT_NUM=10000
export TIME_STEP=1800
export DATA_FILE_NAME=data_${PJM_JOBID}.ptcl

mpiexec ./MassiveSwarmSimCUI_MPI ${DATA_FILE_NAME} ${AGENT_NUM} ${TIME_STEP} example.ini
#mpiexec –mca mpi_print_stats 1 ./MassiveSwarmSimCUI_MPI ${DATA_FILE_NAME} ${AGENT_NUM} ${TIME_STEP} *.ini
