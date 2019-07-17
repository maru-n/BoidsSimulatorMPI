#!/usr/bin/env python

import numpy as np
import sys, os, struct, glob
from swarm_util import SwarmDataManager

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('input_dir')
parser.add_argument('output_file')
parser.add_argument('--time-start', type=int)
parser.add_argument('--time-end', type=int)
parser.add_argument('--velocity', action="store_true")

args = parser.parse_args()

try:
    header_fname = glob.glob(os.path.join(args.input_dir, '*__header'))[0]
except Exception as e:
    print("can not find header file.")
    sys.exit()

data_fnames = glob.glob(os.path.join(args.input_dir, '*__c_*'))

of = open(args.output_file, 'wb')

# get the information from header
sdm = SwarmDataManager(header_fname)
N = sdm.N
sdm.close()

print(f"N = {N}")

with open(header_fname, 'rb') as f:
    h = f.read()
    of.write(h)

data_files = []
for fname in data_fnames:
    if args.velocity:
        if not 'vel' in fname:
            continue
    else:
        if 'vel' in fname:
            continue
    #f = open(fname, 'rb')
    #f = open(fname, 'rb', buffering=4096*1024)
    #f = open(fname, 'rb', buffering=4096*512)
    f = open(fname, 'rb', buffering=4096*256)
    data_files.append(f)

t = 0

X = np.zeros((N, 3), dtype=np.float32)
ids = np.zeros(N, dtype=np.uint32)
finish = False
while True:
    skip = False
    if args.time_start is not None and t < args.time_start:
        skip = True
    elif args.time_end is not None and t >= args.time_end:
        skip = True
        finish = True

    nums = []
    idx = 0
    for f in data_files:
        #print(f)
        d = f.read(4)
        if(len(d)!=4):
            finish = True
            break
        n = struct.unpack("<I", d)[0]
        nums.append(n)

        if skip:
            f.seek(4*n*4,1)
        else:
            try:
                ids[idx:idx+n] = np.ndarray(n, buffer=f.read(4*n), dtype=np.uint32)
                X[idx:idx+n]   = np.ndarray((n, 3), buffer=f.read(4*n*3), dtype=np.float32)
                idx += n
            except Exception as e:
                finish = True

    if finish:
        print("# merge finished.")
        break

    if not skip:
        assert len(np.unique(ids)) == N

        X = X[ids.argsort()]
        X.tofile(of)

    print(f'{t} (max_num:{np.max(nums)} min_num:{np.min(nums)}) {"(skip)" if skip else "(save)"}')
    t += 1

    # if t == 50:
    #     break
