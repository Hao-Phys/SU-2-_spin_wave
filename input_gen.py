#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : input_gen.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.05.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>


import os 
import numpy as np

hx = 0.0
hy = 0.0
hz = 0.0

cwd = os.getcwd()
dir = cwd + '/gs_info/hx_' + str(hx) \
     + 'hy_' + str(hy) + 'hz_' + str(hz) 

if os.path.isdir(dir):
    pass
else:
    os.mkdir(dir)

os.chdir(dir)
fname = 'input.txt'

if os.path.exists(fname):
    pass
else:
    paras = np.array([hx, hy, hz])
    np.savetxt(fname, paras)
