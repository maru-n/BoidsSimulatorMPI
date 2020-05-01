#!/usr/bin/env python


from swarm_util import SwarmDataManager
import h5py

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('input_data')
parser.add_argument('output_data')
parser.add_argument('--setting_file', type=str)

args = parser.parse_args()

sdm = SwarmDataManager(args.input_data)

with h5py.File(args.output_data, 'w') as f:
    f.create_dataset('position', [1,2,3])
    f['data'].attr.create('N'{sdm.N}'))

print(f'
print(f'Space_X_Min={sdm.x_min}')
print(f'Space_X_Max={sdm.x_max}')
print(f'Space_Y_Min={sdm.y_min}')
print(f'Space_Y_Max={sdm.y_max}')
print(f'Space_Z_Min={sdm.z_min}')
print(f'Space_Z_Max={sdm.z_max}')

print(f'Separation_ViewingDistance={0.01}')
print(f'Separation_ViewingAngle={0.5}')
print(f'Separation_ForceCoefficient={0.002}')

print(f'Alignment_ViewingDistance={0.05}')
print(f'Alignment_ViewingAngle={0.3333333333333333}')
print(f'Alignment_ForceCoefficient={0.06}')

print(f'Cohesion_ViewingDistance={0.05}')
print(f'Cohesion_ViewingAngle={0.5}')
print(f'Cohesion_ForceCoefficient={0.008}')

print(f'Velocity_Min={0.001}')
print(f'Velocity_Max={0.005}')

print('')
