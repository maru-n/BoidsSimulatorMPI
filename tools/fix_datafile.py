#!/usr/bin/env python

import sys
from swarm_util import *

def get_overwrite_data(filename):
    sdm = SwarmDataManager(filename)
    steps = sdm.steps
    correct_steps = sdm.valid_steps
    filesize = sdm.file_size
    correct_filesize = sdm.correct_file_size
    sdm.close()
    return steps, correct_steps, filesize, correct_filesize


def ask_user(filename):
    try:
        sdm.check_metadata()
    except Exception as e:
        pass
    else:
        print('Header is correct.')
        return False
    steps, correct_steps, filesize, correct_filesize = get_overwrite_data(filename)

    if steps != correct_steps:
        print('header: step {} -> step {}'.format(steps, correct_steps))
    if filesize != correct_filesize:
        print('filesize: {} -> {}'.format(filesize, correct_filesize))
    choice = input("Overwrite? : ").lower()
    if choice in ['y', 'ye', 'yes']:
        return True
    else:
        return False


def fix_data(filename):
    steps, correct_steps, filesize, correct_filesize = get_overwrite_data(filename)
    with open(filename, 'r+b') as f:
        size = struct.calcsize(HEADER_FORMAT)
        headder = struct.unpack(HEADER_FORMAT, f.read(size))
        headder = list(headder)
        headder[5] = correct_steps
        f.seek(0)
        f.write(struct.pack(HEADER_FORMAT, *headder))
        f.truncate(correct_filesize)

USAGE = """\
clean extra gabage data on the end of file.
USAGE: fix_datafile.py filename\
"""

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(USAGE)
        exit()
    fn = sys.argv[1]
    if ask_user(fn):
        fix_data(fn)
