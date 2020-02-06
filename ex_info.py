#!/usr/bin/env python3 
# -*- coding: utf-8 -*- 
# File              : ex_info.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.15.2020
# Last Modified Date: 02.05.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import SU2rot as rot

num_bond = cf.num_bond
num_sub = cf.num_sub
J1 = cf.J1
J2 = cf.J2
J3 = cf.J3
Delta1 = cf.Delta1
Delta2 = cf.Delta2
Delta3 = cf.Delta3
field = cf.field
sub_idx = cf.sub_idx
R_l2g = rot.R_l2g

theta = rot.theta

Jex = np.zeros((num_bond))
Delta = np.zeros((num_bond))

Jex[0:4] = J1
Jex[4:10] = J2
Jex[10:22] = J3

Delta[0:4] = Delta1
Delta[4:10] = Delta2 
Delta[10:22] = Delta3

tHij = np.zeros((num_bond, 3, 3))

for bond in range(num_bond):
    tmp = Jex[bond] * np.diag(np.array([1.0, 1.0, Delta[bond]]))
    sub1 = sub_idx[bond, 0]
    sub2 = sub_idx[bond, 1]

    R1 = R_l2g[sub1, :, :]
    R2 = R_l2g[sub2, :, :]
    tHij[bond, :, :] = R1.T @ tmp @ R2


