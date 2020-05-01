import numpy as np

SETTING_STR = '''[Global]
# DENSITY is 16384 (2**14)
POPULATION = 2048
FIELD_SIZE = 0.5
INIT = random_uniform
RANDOM_SEED = 12345

[Separation]
SIGHT_DISTANCE = 0.01
SIGHT_ANGLE= 0.5
FORCE_COEFFICIENT = 0.002

[Alignment]
SIGHT_DISTANCE = 0.05
SIGHT_ANGLE= 0.3333333333333333
FORCE_COEFFICIENT = 0.06

[Cohesion]
SIGHT_DISTANCE = 0.05
SIGHT_ANGLE= 0.5
FORCE_COEFFICIENT = 0.008

[Velocity]
MIN = {:.4f}
MAX = {:.4f}
'''

FNAME = 'vmin{:.4f}_vmax{:.4f}.ini'

for vmin in np.arange(0, 0.006, 0.001):
    for vmax in np.arange(vmin, 0.006, 0.001):
        if vmin == vmax:
            continue
        if vmax == 0.005:
            continue
        fname = FNAME.format(vmin, vmax)
        str = SETTING_STR.format(vmin, vmax)

        with open(fname, 'w') as f:
            f.write(str)
