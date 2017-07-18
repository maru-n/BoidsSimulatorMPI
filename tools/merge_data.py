#!/usr/bin/env python

import numpy as np
import sys
import struct

header_fname = sys.argv[1]
data_fnames = sys.argv[2:-1]
output_fname = sys.argv[-1]

of = open(output_fname, 'wb')

with open(header_fname, 'rb') as f:
    h = f.read()
    of.write(h)

data_files = []
for fname in data_fnames:
    f = open(fname, 'rb')
    data_files.append(f)

t = 0
while True:
    X = np.zeros((0, 3), dtype=np.float32)
    ids = np.zeros(0)
    for f in data_files:
        d = f.read(4)
        n = struct.unpack("<I", d)[0]

        d = f.read(4*n)
        i = np.ndarray(n, buffer=d, dtype=np.uint32)
        d = f.read(4*n*3)
        x = np.ndarray((n, 3), buffer=d, dtype=np.float32)
    
        X = np.r_[X, x]
        ids = np.r_[ids, i]
    assert len(np.unique(ids)) == X.shape[0]

    X = X[ids.argsort()]
    X.tofile(of)

    print(t)
    t += 1
