#!/usr/bin/env python

import argparse
import os.path
from struct import *



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check result of calculation on K from output files. (report files and data file)')
    """
    parser.add_argument('error', metavar='ERROR', type=str,
                        help='error file of calculation on K. ({job_name}.e{job_id})')
    parser.add_argument('info', metavar='INFO', type=str,
                        help='info file of calculation on K. ({job_name}.i{job_id})')
    parser.add_argument('stdout', metavar='STDOUT', type=str,
                        help='stdout file of calculation on K. ({job_name}.s{job_id})')
    parser.add_argument('summary', metavar='SUMMARY', type=str,
                        help='summary file of calculation on K. ({job_name}.s{job_id})')
    parser.add_argument('data', metavar='DATA', type=str,
                        help='data file of calculation on K. (xxx.ptcl)')
    """

    args = parser.parse_args()

    '''
    typedef struct {
        char format_name[4] = {'P', 'T', 'C', 'L'};
        unsigned char major_version = 0;
        unsigned char minor_version = 2;
        unsigned int header_length;
        unsigned int N;
        unsigned int step;
        unsigned int t_0;
        unsigned int fps;
        float x_min;
        float x_max;
        float y_min;
        float y_max;
        float z_min;
        float z_max;
    } data_file_header_v02;
    '''

    header_v002_fmt = ""

    print(os.path.getsize(args.data))
    data_file = open(args.data, 'rb')
    d = data_file.read(10)
    print(d)
