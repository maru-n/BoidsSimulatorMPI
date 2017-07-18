#!/usr/bin/env python

import sys
from swarm_util import *
import struct

def merge_data(in_fname_list, out_fname, verbose=True):
    for ifn in in_fname_list:
        sdm = SwarmDataManager(ifn)
        print(sdm.N)
    with open(out_fname, 'wb') as of:
        pass

if __name__ == '__main__':
    infile = sys.argv[1:-1]
    outfile = sys.argv[-1]
    merge_data(infile_list, outfile)
