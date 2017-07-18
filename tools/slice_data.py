#!/usr/bin/env python

import sys
from swarm_util import *
import struct

def slice_data(in_fname, out_fname, start, steps):
    sdm = SwarmDataManager(in_fname)
    N = sdm.N
    sdm.close()
    with open(in_fname, 'rb') as inf, open(out_fname, 'wb') as outf:
        hdr_size = struct.calcsize(HEADER_FORMAT)
        headder = struct.unpack(HEADER_FORMAT, inf.read(hdr_size))
        headder = list(headder)
        headder[5] = steps  # steps
        headder[6] = start  # t_0
        hdr_data = struct.pack(HEADER_FORMAT, *headder)
        outf.write(hdr_data)

        one_step_datasize = N * 12
        inf.seek(hdr_size + one_step_datasize*(start - sdm.t_0))
        for i in range(start, start+steps):
            data = inf.read(one_step_datasize)
            if len(data) != one_step_datasize:
                raise Exception('data length is not suitable. chack the data.')
            outf.write(data)


if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    start = int(sys.argv[3])
    steps = int(sys.argv[4])
    slice_data(infile, outfile, start, steps)
