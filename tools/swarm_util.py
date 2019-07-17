#!/usr/bin/env python

import struct


HEADER_FORMAT = "<4sBBIIIIIffffff"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)

def load_metadata(filename):
    metadata = {}
    with open(filename, 'rb') as f:
        bindata = struct.unpack(HEADER_FORMAT, f.read(HEADER_SIZE))
        metadata['N'] = bindata[4]
        metadata['steps'] = bindata[5]
        metadata['x_min'] = bindata[8]
        metadata['x_max'] = bindata[9]
        metadata['y_min'] = bindata[10]
        metadata['y_max'] = bindata[11]
        metadata['z_min'] = bindata[12]
        metadata['z_max'] = bindata[13]
    return metadata


def load_data(file_name, time_step, vel_calc_step=1, end_time=None, time_interval=1):
    return load_data_legacy(file_name, time_step, vel_calc_step=1, end_time=None, time_interval=1)

# def load_data(file_name, time_step, vel_calc_step=1, end_time=None, time_interval=1):
#     sdm = SwarmDataManager(filename)
#     if end_time is None:
#         end_time = time + 1
#
#     N = load_metadata(file_name)['N']
#     STEP_DATA_SIZE = struct.calcsize('f') * 3 * N
#
#     X = np.empty((N, 6, (end_time-time)//time_interval))
#     with open(filename, 'rb') as f:
#
#         seek_pos = HEADER_SIZE + time_step * STEP_DATA_SIZE
#         f.seek(seek_pos)
#
#         buffer = f.read(self.__step_data_size * s)
#         data = np.ndarray((s, self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
#         x, v = calc_vel(sdm, t, step=vel_calc_step)
#         X[:,0:3,i] = x
#         X[:,3:6,i] = v
#
#     def set_timestep(self, step):
#         self.file.
    #
    #
    # if end_time is None and time_interval is None:
    #     return calc_vel(sdm, time, step=vel_calc_step)
    # else:
    #     X = np.empty((sdm.N, 6, (end_time-time)//time_interval))
    #     for i, t in enumerate(range(time, end_time, time_interval)):
    #         x =
    #         X[:,0:3,i] = x
    #         X[:,3:6,i] = v
    #     return X



################
# Legacy code
################

import sys
import os
import numpy as np
import struct
import matplotlib.pyplot as plt
import pandas as pd
import glob
import re
import warnings
from mpl_toolkits.mplot3d.axes3d import Axes3D

#DATA_DIM = 3
#DATA_DIM = 12

HEADER_FORMAT = "<4sBBIIIIIffffff"

class SwarmDataManager(object):
    def __init__(self, filename, force=False):
        self.DATA_DIM = 3
        if force:
            self.DATA_DIM = 12
        self.__header_size = struct.calcsize(HEADER_FORMAT)
        self.file_name = filename
        self.file = open(filename, 'rb')
        metadata = struct.unpack(HEADER_FORMAT, self.file.read(self.__header_size))
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
        #self.__step_data_format = "<{}f".format(self.N*3)
        self.__step_data_format = "<{}f".format(self.N*self.DATA_DIM)
        self.__step_data_size = struct.calcsize(self.__step_data_format)
        self.file_size = os.path.getsize(self.file_name)
        self.valid_steps = (self.file_size - self.__header_size) // self.__step_data_size
        self.correct_file_size = self.__header_size + self.__step_data_size * self.valid_steps

    def check_metadata(self):
        if not self.valid_steps == self.steps:
           raise Exception('Header is {} steps, but data size can contain only {} steps'.format(self.steps, self.valid_steps))
        if not (self.file_size - self.__header_size) % self.__step_data_size == 0:
            raise Exception('Data size is {}, but {} step data have to be {}'.format(self.file_size, self.valid_steps, self.correct_file_size))
        return True

    def set_timestep(self, step):
        self.file.seek(self.__header_size + step * self.__step_data_size)

    def read_all(self):
        self.set_timestep(0)
        bd = self.file.read(self.__step_data_size * self.valid_steps)
        data = np.ndarray((self.valid_steps, self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
        return data

    def read_positions(self, start, stop):
        self.set_timestep(start)
        s = stop - start
        bd = self.file.read(self.__step_data_size * s)
        data = np.ndarray((s, self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
        return data

    def read_pos(self, step=None):
        if step is not None:
            self.set_timestep(step)
        bd = self.file.read(self.__step_data_size)
        data = np.ndarray((self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
        return data

    def read_vel(self, step=None):
        warnings.warn('This method is duplicated. Use calc_vel(sdm, time, step)')
        if step is not None:
            self.set_timestep(step)
        self.file.seek(-self.__step_data_size, 1)
        bd = self.file.read(self.__step_data_size)
        pre_pos = np.ndarray((self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
        bd = self.file.read(self.__step_data_size)
        pos = np.ndarray((self.N, self.DATA_DIM), buffer=bd, dtype=np.float32)
        return pos - pre_pos

    def close(self):
        self.file.close()

def load_metadata_legacy(filename):
    sdm = SwarmDataManager(filename)
    metadata = dict(
        file_name = sdm.file_name,
        N = sdm.N,
        step_max = sdm.steps,
        x_min = sdm.x_min,
        x_max = sdm.x_max,
        y_min = sdm.y_min,
        y_max = sdm.y_max,
        z_min = sdm.z_min,
        z_max = sdm.z_max
    )
    return metadata


def load_data_legacy(filename, time, vel_calc_step=1, end_time=None, time_interval=None):
    sdm = SwarmDataManager(filename)
    if end_time is None and time_interval is None:
        return calc_vel(sdm, time, step=vel_calc_step)
    else:
        X = np.empty((sdm.N, 6, (end_time-time)//time_interval))
        for i, t in enumerate(range(time, end_time, time_interval)):
            x, v = calc_vel(sdm, t, step=vel_calc_step)
            X[:,0:3,i] = x
            X[:,3:6,i] = v
        return X

def load_pos(filename, time):
    sdm = SwarmDataManager(filename)
    sdm.set_timestep(time)
    X = sdm.read_pos()
    return X

def calc_vel(sdm, time, step=1):
    sdm.set_timestep(time)
    x0 = sdm.read_pos()

    offset = np.zeros(x0.shape)
    x1 = x0.copy()
    # chech every step data.
    # if particle cross over boundary, add field size value to offset
    for i in range(step):
        x2 = sdm.read_pos()
        v = x2 - x1

        offset -= np.sign(v) * (np.abs(v) > sdm.x_size/2) * sdm.x_size
        """
        rows, cols = np.where(v > sdm.x_size/2)
        for i, j in zip(rows, cols):
            offset[i,j] -= sdm.x_size
        rows, cols = np.where(v < -sdm.x_size/2)
        for i, j in zip(rows, cols):
            offset[i,j] += sdm.x_size
        """
        x1 = x2
    v =  x1 - x0 + offset
    v = v / step

    return x0, v

def calc_vel_legacy(sdm, time, step=1):
    sdm.set_timestep(time)
    x0 = sdm.read_pos()

    offset = np.zeros(x0.shape)
    x1 = x0.copy()
    # chech every step data.
    # if particle cross over boundary, add field size value to offset
    for i in range(step):
        x2 = sdm.read_pos()
        v = x2 - x1
        rows, cols = np.where(v > sdm.x_size/2)
        for i, j in zip(rows, cols):
            offset[i,j] -= sdm.x_size
        rows, cols = np.where(v < -sdm.x_size/2)
        for i, j in zip(rows, cols):
            offset[i,j] += sdm.x_size
        x1 = x2
    v =  x1 - x0 + offset
    v = v / step

    return x0, v

def plt_xv(x, v, dim=(0,1)):
    plt.quiver(x[:,dim[0]], x[:,dim[1]], v[:,dim[0]], v[:,dim[1]], scale=0.25)


def plt_xvc(x, v, c=None, dim=(0,1), type='all', only_class=-1, cmap=plt.cm.prism):
    if type is 'all':
        plt.quiver(x[:,dim[0]], x[:,dim[1]], v[:,dim[0]], v[:,dim[1]], c, scale=0.25, cmap=cmap)
    else:
        x0 = x[c==-1]
        v0 = v[c==-1]
        x1 = x[c!=-1]
        v1 = v[c!=-1]
        c1 = c[c!=-1]
        if type is 'class':
            plt.quiver(x0[:,dim[0]], x0[:,dim[1]], v0[:,dim[0]], v0[:,dim[1]], scale=0.25, color='grey', alpha=0.2)
            plt.quiver(x1[:,dim[0]], x1[:,dim[1]], v1[:,dim[0]], v1[:,dim[1]], c1, scale=0.25, cmap=cmap)
        elif type is 'unclass':
            plt.quiver(x0[:,dim[0]], x0[:,dim[1]], v0[:,dim[0]], v0[:,dim[1]], scale=0.25, color='red')
            plt.quiver(x1[:,dim[0]], x1[:,dim[1]], v1[:,dim[0]], v1[:,dim[1]], scale=0.25, color='gray', alpha=0.2)
        elif type is 'onlyclass':
            x0 = x[c==only_class]
            v0 = v[c==only_class]
            x1 = x[c!=only_class]
            v1 = v[c!=only_class]
            plt.quiver(x0[:,dim[0]], x0[:,dim[1]], v0[:,dim[0]], v0[:,dim[1]], scale=0.25, color='red')
            plt.quiver(x1[:,dim[0]], x1[:,dim[1]], v1[:,dim[0]], v1[:,dim[1]], scale=0.25, color='gray', alpha=0.2)

def plt_3d_xc(x, c=None, cmap=plt.cm.prism, plot_class='all', point_size=1, axis=None):
    if axis is None:
        ax = Axes3D(plt.gcf())
    else:
        ax = axis
    if plot_class=='all':
        ax.scatter3D(x[:,0], x[:, 1], x[:, 2], ',', c=c, cmap=cmap, s=point_size)
    elif plot_class=='class':
        x1 = x[c==-1]
        c1 = c[c==-1]
        ax.scatter3D(x1[:,0], x1[:, 1], x1[:, 2], ',', c='gray', alpha=0.01)
        x0 = x[c!=-1]
        c0 = c[c!=-1]
        ax.scatter3D(x0[:,0], x0[:, 1], x0[:, 2], ',', c=c0, cmap=cmap, s=point_size)
    elif type(plot_class) is int:
        x1 = x[c!=plot_class]
        c1 = c[c!=plot_class]
        ax.scatter3D(x1[:,0], x1[:, 1], x1[:, 2], ',', c='gray', alpha=0.01)
        x0 = x[c==plot_class]
        c0 = c[c==plot_class]
        ax.scatter3D(x0[:,0], x0[:, 1], x0[:, 2], ',', c=c0, cmap=cmap, s=point_size)


def load_param_datas_dataframe(data_dir, PARAM_SET_LIST = ['mototake', 'one_body'], N_LIST = [2**n for n in range(10,20)]):
    index = pd.MultiIndex.from_product([PARAM_SET_LIST, N_LIST], names=['param_set', 'N'])
    sdms = [None] * len(index)
    sdms = pd.Series(sdms, index=index)

    file_list = glob.glob(data_dir + "*ptcl")

    for path in file_list:
        fname = os.path.basename(path)
        r = re.search('(.+)_N([0-9]{7})_T[0-9]+_[0-9]{7}\.ptcl',fname)
        if r is None:
            continue
        ps = r.groups()[0]
        n = int(r.groups()[1])
        try:
            sdms[ps][n] = SwarmDataManager(path)
        except KeyError:
            pass
    return sdms


def gen_metric_func(alpha):
    def metric_func(x0, x1):
        return (alpha * np.linalg.norm(x0[:3]-x1[:3])) + ((1-alpha) * np.linalg.norm(x0[3:]-x1[3:]))
    return metric_func


def divide_class(x, c):
    cls = np.unique(c)
    res = []
    for cc in cls:
        res.append(x[c==cc])
    return res, cls

if __name__ == '__main__':
    sd = SwarmData(sys.argv[1])
    for i in range(10):
        sd.set_timestep(i)
