#!/usr/bin/env python

import pystache
import os
from datetime import timedelta
import argparse
import warnings


def flatten_list(nested_list):
    if type(nested_list) is not list:
        return [nested_list]
    res = []
    for l in nested_list:
        res += flatten_list(l)
    return res


def generate_job_name(setting_file, n, t):
    return "{}_N{}_T{}".format(os.path.splitext(os.path.basename(setting_file))[0], str(n).zfill(7), t)


def predict_calculation_time(n, k, step):
    # parameters calculated by notebook.
    if k == 1:
        a, b = [  3.33487120e-08, 1.19010289e-06]
        return (a*n**2 + b*n) * step
    else:
        a, b, c = [  9.17194851e-08,   3.95157428e-06,   1.23766077e+00]
        return (a*(n**2)/(k**6) + b*n + c*(n**(1/3))/k) * step


dirname = os.path.dirname(__file__)
template_single = open(os.path.join(dirname, "jobscript_templates/template.sh")).read()
template_mpi = open(os.path.join(dirname, "jobscript_templates/template_mpi.sh")).read()

def generate_script(setting, population, field_size, time_step, data_safix=""):
    # predict calculation time and if it end on a day, set number of node
    max_node = 6
    max_time = 24*60*60
    for k in range(1, max_node+1):
        t = predict_calculation_time(population, k, time_step)
        t = int(t * 1.5)  # margin
        if k == max_node and t > max_time:
            print("## Warning: this setting not expected end on a day with {}**3 nodes. N:{} time_step:{}".format(max_node, population, time_step))
            t = max_time-1
        if t <= max_time:
            template_data = {
                'name'       : generate_job_name(setting, population, time_step),
                'setting'    : os.path.basename(setting),
                'node'       : '',
                'elapse'     : str(timedelta(seconds=t)),
                'time_step'  : time_step,
                'agent_num'  : population,
                'field_size' : str(field_size),
                'data_safix' : data_safix,
            }
            if k == 1:
                template = template_single
                template_data['node'] = '1'
            else:
                template = template_mpi
                template_data['node'] = '{}x{}x{}'.format(k, k, k)
            return pystache.render(template, template_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate job script for K computer.')
    parser.add_argument('-s', dest='settings', nargs='+', type=str, required=True,
                        help='list of setting files.')
    parser.add_argument('-n', dest='populations', nargs='+', type=str, required=True,
                        help='list of population by python expression, ex) 1000 ,2**10, [2,3] or [i for i in range(10)]')
    parser.add_argument('-d', dest='density', type=float, required=True,
                        help='density of individuals. field size will be calculated by this value.')
    #parser.add_argument('-f', dest='field_size', nargs='+', type=str,
    #                    help='list of population by python expression, ex) 1000 or 2**10.')
    parser.add_argument('-t', dest='time_step', type=int, required=True,
                        help='time step.')
    parser.add_argument('-o', dest='output_dir',
                        help='directory to save output scipts.')


    args = parser.parse_args()

    settings = [s for s in args.settings]
    populations = flatten_list([eval(s) for s in args.populations])
    field_size = [(n/args.density)**(1/3) for n in populations]
    time_step = args.time_step
    output_dir = args.output_dir

    for s in settings:
        for n, fs in zip(populations, field_size):
            script_fname = generate_job_name(s, n, time_step) + ".sh"
            if output_dir:
                script_fname = os.path.join(output_dir, script_fname)
            print(script_fname)
            f = open(script_fname, 'w')
            script = generate_script(s, n, fs, time_step)
            f.write(script)
