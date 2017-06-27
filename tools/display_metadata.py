#!/usr/bin/env python

import sys
from swarm_data_manager import SwarmDataManager

if __name__ == '__main__':
    sdm = SwarmDataManager(sys.argv[1])
    print("N    : {}".format(sdm.N))
    print("t_0  : {}".format(sdm.t_0))
    print("step : {}".format(sdm.step))
