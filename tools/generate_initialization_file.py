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
            for i in range(N):
                d = pf.read(4*3)
                out_f.write(d)
                d = vf.read(4*3)
                out_f.write(d)