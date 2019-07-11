#!/usr/bin/env python

import sys
from swarm_util import *
import struct, shutil
import fix_datafile

USAGE = """\
USAGE: merge_datafiles.py input_file_1 input_file_2 ... output_file
"""

def merge_data(in_fname_list, out_fname, verbose=True):
    shutil.copy2(in_fname_list[0], out_fname)
    with open(out_fname, 'ab') as outf:
        outf.seek(0, 2)
        for in_fname in in_fname_list[1:]:
            with open(in_fname, 'rb') as inf:
                inf.seek(50)
                while 1:
                    d = inf.read(1)
                    if not d:
                        break
                    outf.write(d)
    fix_datafile.fix_data(out_fname)

if __name__ == '__main__':
    infile_list = sys.argv[1:-1]
    outfile = sys.argv[-1]
    merge_data(infile_list, outfile)
