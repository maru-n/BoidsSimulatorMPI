#!/usr/bin/env python

import sys
import warnings
from swarm_util import *
import argparse

parser = argparse.ArgumentParser(description='check metadata of swarm simulation data file.')
parser.add_argument('input_files', type=str, nargs='+',
                    help='swarm simulation data files')
parser.add_argument('--verbose', '-v', action='store_true')

args = parser.parse_args()

fnames = args.input_files

for fn in fnames:
    print(fn)
    sdm = SwarmDataManager(fn)
    if args.verbose:
        print('  N:{} t0:{} step:{}'.format(sdm.N, sdm.t_0, sdm.steps), end='')
    else:
        print('  N:{} t0:{} step:{}'.format(sdm.N, sdm.t_0, sdm.steps), end='')
    try:
        sdm.check_metadata()
    except Exception as e:
        print('\033[91m Invalid metadata: ', end='')
        print(e, end='')
        print('\033[0m', end='')

    print()
