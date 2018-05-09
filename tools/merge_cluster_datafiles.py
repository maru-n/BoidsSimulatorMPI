#!/usr/bin/env python

import numpy as np
import sys, os, struct, glob
from swarm_util import SwarmDataManager

USAGE = """\
USAGE: merge_cluster.py header_file data_file_1 data_file_2 ... output_file\n \
       merge_cluster.py data_file_dir output_file\n \
"""

INITIAL_TIME = 0

if len(sys.argv) < 3:
    print(USAGE)
    sys.exit()

if len(sys.argv) == 3:
    data_dir = sys.argv[1]
    header_fname = glob.glob(os.path.join(data_dir, '*__header'))[0]
    data_fnames = glob.glob(os.path.join(data_dir, '*__c_*'))
    output_fname = sys.argv[2]
else:
    header_fname = sys.argv[1]
    data_fnames = sys.argv[2:-1]
    output_fname = sys.argv[-1]

of = open(output_fname, 'wb')

# get the information from header
sdm = SwarmDataManager(header_fname)
N = sdm.N
sdm.close()

with open(header_fname, 'rb') as f:
    h = f.read()
    of.write(h)

data_files = []
for fname in data_fnames:
    #f = open(fname, 'rb')
    #f = open(fname, 'rb', buffering=4096*1024)
    #f = open(fname, 'rb', buffering=4096*512)
    f = open(fname, 'rb', buffering=4096*256)
    data_files.append(f)

t = 0

X = np.zeros((N, 3), dtype=np.float32)
ids = np.zeros(N)
while True:
    nums = []
    idx = 0
    for f in data_files:
        d = f.read(4)
        n = struct.unpack("<I", d)[0]
        nums.append(n)

        ids[idx:idx+n] = np.ndarray(n, buffer=f.read(4*n), dtype=np.uint32)
        X[idx:idx+n]   = np.ndarray((n, 3), buffer=f.read(4*n*3), dtype=np.float32)
        idx += n


    assert len(np.unique(ids)) == X.shape[0]

    X = X[ids.argsort()]
    X.tofile(of)

    print(t, "(max_num:{} min_num:{})".format(np.max(nums), np.min(nums)))
    t += 1

    # if t == 50:
    #     break
