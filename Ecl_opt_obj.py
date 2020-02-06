#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : Ecl_opt_obj.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.05.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf

num_bond = cf.num_bond
num_sub = cf.num_sub
spin = cf.spin
J1 = cf.J1
J2 = cf.J2
J3 = cf.J3
Delta1 = cf.Delta1
Delta2 = cf.Delta2
Delta3 = cf.Delta3
hx = cf.hx
hy = cf.hy
hz = cf.hz
sub_idx = cf.sub_idx

Jex = np.zeros((num_bond))
Delta = np.zeros((num_bond))

Jex[0:4] = J1
Jex[4:10] = J2
Jex[10:22] = J3

Delta[0:4] = Delta1
Delta[4:10] = Delta2 
Delta[10:22] = Delta3

def Energ_cl(v):

    funval = 0.0

    for bond in range(num_sub):
        id1 = sub_idx[bond, 0]
        id2 = sub_idx[bond, 1]

        theta1 = v[2*id1]
        phi1   = v[2*id1+1]

        theta2 = v[2*id2]
        phi2   = v[2*id2+1]

        S1 = spin*np.array([np.sin(theta1)*np.cos(phi1), \
                            np.sin(theta1)*np.sin(phi1), \
                            np.cos(theta1)])

        S2 = spin*np.array([np.sin(theta2)*np.cos(phi2), \
                            np.sin(theta2)*np.sin(phi2), \
                            np.cos(theta2)])

        # exchange tensor in global reference frame
        tmp = Jex[bond] * np.diag(np.array([1.0, 1.0, Delta[bond]]))

        funval += S1.T @ tmp @ S2

    for sub_lat in range(num_sub):

        theta = v[2*sub_lat]
        phi   = v[2*sub_lat]
        Si = spin*np.array([np.sin(theta)*np.cos(phi), \
                            np.sin(theta)*np.sin(phi), \
                            np.cos(theta)])

        h_ext = np.array([hx, hy, hz])
        E_Zeeman = - h_ext @ Si
        funval += E_Zeeman

    return funval
