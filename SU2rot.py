#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : SU2rot.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.15.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import os

num_sub = cf.num_sub
hx = cf.hx
hy = cf.hy
hz = cf.hz
path = cf.path

field = np.sqrt(hx**2 + hy**2 + hz**2)

def R_rot(theta, phi):

    mat = np.array([[-np.sin(phi), -np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta)], \
                    [np.cos(phi), -np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta)], \
                    [0.0, np.sin(theta), np.cos(theta)]])

    return mat


if field == 0:

    theta1 = 0.5*np.pi
    phi1 = 0.0
    phi2 = np.pi
    theta = np.array([theta1, theta1, theta1, theta1, theta1, theta1])
    phi = np.array([phi1, phi1, phi1, phi1, phi2, phi2])

else:
    fname = path + 'opt_angles.txt'
    reason = os.path.isfile(fname)

    if reason == 0:
        print("please generate the optimal order file first")
        exit()

    theta = np.zeros((num_sub))
    phi = np.zeros((num_sub))
    angles = np.loadtxt(fname)

    for flag in range(num_sub):
        theta[flag] = angles[0]
        phi[flag] = angles[1]

R_l2g = np.zeros((num_sub, 3, 3))

for ii in range(num_sub):
    R_l2g[ii, :, :] = R_rot(theta[ii], phi[ii])


