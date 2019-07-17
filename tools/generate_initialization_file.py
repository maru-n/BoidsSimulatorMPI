#!/usr/bin/env python

import sys
from swarm_util import SwarmDataManager


pos_data_file_name = sys.argv[1]
vel_data_file_name = sys.argv[2]
time_step = int(sys.argv[3])
out_data_file_name = sys.argv[4]

sdm = SwarmDataManager(pos_data_file_name)
N = sdm.N
sdm.close()

with open(pos_data_file_name, 'rb') as pf:
    with open(vel_data_file_name, 'rb') as vf:
        pf.seek(50 + N * 4 * 3 * time_step)
        vf.seek(50 + N * 4 * 3 * time_step)
        with open(out_data_file_name, 'wb') as out_f:
            pos_data = pf.read(4*3*N)
            vel_data = vf.read(4*3*N)
            assert len(pos_data) == 4*3*N
            assert len(vel_data) == 4*3*N
            for i in range(N):
                out_f.write(pos_data[4*3*i:4*3*(i+1)])
                out_f.write(vel_data[4*3*i:4*3*(i+1)])
