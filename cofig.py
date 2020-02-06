#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : cofig.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.15.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import os
import sys

inFile = sys.argv[1]
paras = np.loadtxt(inFile)

spin = 2.5
num_sub = 6
num_bond = 22

J1 = 1.295
J2 = 0.371
J3 = 0.890
Delta1 = 1.16
Delta2 = 0.0
Delta3 = 0.44

# external field Cartesian frame
hx = paras[0]
hy = paras[1]
hz = paras[2]


print("hx= ", hx)
print("hy= ", hy)
print("hz= ", hz)

cwd = os.getcwd()
path = cwd + '/gs_info/hx_' + str(hx) \
     + 'hy_' + str(hy) + 'hz_' + str(hz) +'/' 
# mu_B = 5.788381E-2

A_mat = np.zeros((2*num_sub, 2*num_sub))

for flag in range(num_sub):
    A_mat[flag, flag] = 1

for flag in range(num_sub, 2*num_sub):
    A_mat[flag, flag] = -1

# the bond information

a1 = np.array([1.0, 0.0, 0.0])
a2 = np.array([0.0, 1.0, 0.0])
a3 = np.array([0.0, 0.0, 1.0])

sub_idx = np.zeros((num_bond, 2))
delta_ij = np.zeros((3, num_bond))

sub_idx[0, 0] = 0
sub_idx[0, 1] = 4
delta_ij[:, 0] = a3/4

sub_idx[1, 0] = 1
sub_idx[1, 1] = 4
delta_ij[:, 1] = -a3/4

sub_idx[2, 0] = 2
sub_idx[2, 1] = 5
delta_ij[:, 2] = -a3/4

sub_idx[3, 0] = 3
sub_idx[3, 1] = 5
delta_ij[:, 3] = a3/4

sub_idx[4, 0] = 2
sub_idx[4, 1] = 0
delta_ij[:, 4] = -a1/3 + a2/3

sub_idx[5, 0] = 2
sub_idx[5, 1] = 0
delta_ij[:, 5] = -a1/3 - 2*a2/3

sub_idx[6, 0] = 2
sub_idx[6, 1] = 0
delta_ij[:, 6] = 2*a1/3 + a2/3

sub_idx[7, 0] = 3
sub_idx[7, 1] = 1
delta_ij[:, 7] = -a1/3+a2/3

sub_idx[8, 0] = 3
sub_idx[8, 1] = 1
delta_ij[:, 8] = 2*a1/3 + a2/3

sub_idx[9, 0] = 3
sub_idx[9, 1] = 1
delta_ij[:, 9] = -a1/3 - 2*a2/3

sub_idx[10, 0] = 3
sub_idx[10, 1] = 4
delta_ij[:, 10] = -a1/3 + a2/3 - a3/4

sub_idx[11, 0] = 3
sub_idx[11, 1] = 4
delta_ij[:, 11] = -a1/3 - 2*a2/3 - a3/4

sub_idx[12, 0] = 3
sub_idx[12, 1] = 4
delta_ij[:, 12] = 2*a1/3 + a2/3 - a3/4

sub_idx[13, 0] = 2
sub_idx[13, 1] = 4
delta_ij[:, 13] = -a1/3 + a2/3 + a3/4

sub_idx[14, 0] = 2
sub_idx[14, 1] = 4
delta_ij[:, 14] = 2*a1/3 + a2/3 + a3/4

sub_idx[15, 0] = 2
sub_idx[15, 1] = 4
delta_ij[:, 15] = -a1/3 - 2*a2/3 + a3/4

sub_idx[16, 0] = 1
sub_idx[16, 1] = 5
delta_ij[:, 16] = a1/3 - a2/3 + a3/4

sub_idx[17, 0] = 1
sub_idx[17, 1] = 5
delta_ij[:, 17] = a1/3 + 2*a2/3 + a3/4

sub_idx[18, 0] = 1
sub_idx[18, 1] = 5
delta_ij[:, 18] = -2*a1/3 - a2/3 + a3/4

sub_idx[19, 0] = 0
sub_idx[19, 1] = 5
delta_ij[:, 19] = a1/3 - a2/3 - a3/4

sub_idx[20, 0] = 0
sub_idx[20, 1] = 5
delta_ij[:, 20] = -2*a1/3 - a2/3 - a3/4

sub_idx[21, 0] = 0
sub_idx[21, 1] = 5
delta_ij[:, 21] = a1/3 + 2*a2/3 - a3/4

sub_idx =sub_idx.astype(int)



