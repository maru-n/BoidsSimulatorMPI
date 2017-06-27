#!/usr/bin/env python

import sys
import os
import numpy as np
import struct

class SwarmDataManager(object):
    def __init__(self, filename):
        self.__header_format = "<4sBBIIIIIffffff"
        self.__header_size = struct.calcsize(self.__header_format)
        self.file_name = filename
        self.file = open(filename, 'rb')
        metadata = struct.unpack(self.__header_format, self.file.read(self.__header_size))
        self.file_type = metadata[0]
        self.N = metadata[4]
        self.steps = metadata[5]
        self.t_0 = metadata[6]
        self.fps = metadata[7]
        self.x_min = metadata[8]
        self.x_max = metadata[9]
        self.y_min = metadata[10]
        self.y_max = metadata[11]
        self.z_min = metadata[12]
        self.z_max = metadata[13]
        self.x_size = self.x_max - self.x_min
        self.y_size = self.y_max - self.y_min
        self.z_size = self.z_max - self.z_min
        self.__step_data_format = "<{}f".format(self.N*3)
        self.__step_data_size = struct.calcsize(self.__step_data_format)
        self.file_size = os.path.getsize(self.file_name)
        self.valid_steps = (self.file_size - self.__header_size) // self.__step_data_size

    def is_valid_data(self):
        return (self.valid_steps == self.steps)

    def set_timestep(self, step):
        self.file.seek(self.__header_size + step * self.__step_data_size)

    def read_all(self):
        self.set_timestep(0)
        bd = self.file.read(self.__step_data_size * self.valid_steps)
        data = np.ndarray((self.valid_steps, self.N, 3), buffer=bd, dtype=np.float32)
        return data

    def read_positions(self, start, stop):
        self.set_timestep(start)
        s = stop - start
        bd = self.file.read(self.__step_data_size * s)
        data = np.ndarray((s, self.N, 3), buffer=bd, dtype=np.float32)
        return data

    def read_pos(self, step=None):
        if step is not None:
            self.set_timestep(step)
        bd = self.file.read(self.__step_data_size)
        data = np.ndarray((self.N, 3), buffer=bd, dtype=np.float32)
        return data

    def read_vel(self, step=None):
        if step is not None:
            self.set_timestep(step)
        self.file.seek(-self.__step_data_size, 1)
        bd = self.file.read(self.__step_data_size)
        pre_pos = np.ndarray((self.N, 3), buffer=bd, dtype=np.float32)
        bd = self.file.read(self.__step_data_size)
        pos = np.ndarray((self.N, 3), buffer=bd, dtype=np.float32)
        return pos - pre_pos


if __name__ == '__main__':
    sd = SwarmData(sys.argv[1])
    for i in range(10):
        sd.set_timestep(i)
