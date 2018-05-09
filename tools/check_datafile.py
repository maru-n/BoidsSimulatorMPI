#!/usr/bin/env python

import sys
import warnings
from swarm_util import *

USAGE = """\
USAGE: check_datafile.py file1 file2 ...\
"""

if len(sys.argv) < 2:
    print(USAGE)
    sys.exit()

fnames = sys.argv[1:]

for fn in fnames:
    print(fn)
    sdm = SwarmDataManager(fn)
    print('  N:{} t0:{} step:{}'.format(sdm.N, sdm.t_0, sdm.steps), end='')
    try:
        sdm.check_metadata()
    except Exception as e:
        print('\033[91m Invalid metadata: ', end='')
        print(e, end='')
        print('\033[0m', end='')

    print()