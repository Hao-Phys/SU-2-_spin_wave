#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : GLSW.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.15.2020
# Last Modified Date: 02.05.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import ex_info as exi
import SU2rot as rot

spin = cf.spin
num_bond = cf.num_bond
num_sub = cf.num_sub
delta_ij = cf.delta_ij
sub_idx = cf.sub_idx
spin = cf.spin
A_mat = cf.A_mat
hx = cf.hx
hy = cf.hy
hz = cf.hz

tHij = exi.tHij
theta = np.repeat(rot.theta, 2)

def sw_hamiltonian(q):

    ham11 = np.zeros((num_sub, num_sub), dtype=complex)
    ham22 = np.zeros((num_sub, num_sub), dtype=complex)
    ham12 = np.zeros((num_sub, num_sub), dtype=complex)
    ham21 = np.zeros((num_sub, num_sub), dtype=complex)

    for bond in range(num_bond):

        J_local = tHij[bond, :, :]

        J_a = - spin * J_local[2, 2]
        J_b = 0.5 * spin * (J_local[0, 0] - 1j*J_local[1, 0] \
                           -1j*J_local[0, 1] + J_local[1, 1])
        cJ_b = J_b.conj()
        J_c = 0.5 * spin * (J_local[0, 0] - 1j*J_local[1, 0] \
                           -1j*J_local[0, 1] - J_local[1, 1])
        cJ_c = J_c.conj()

        sub1 = sub_idx[bond, 0]
        sub2 = sub_idx[bond, 1]
        bond_vec = delta_ij[:, bond]
        phase = np.exp(1j*2.0*np.pi*q @ bond_vec)
        cphase = phase.conj()

        ham11[sub1, sub1] += 0.5*J_a
        ham11[sub2, sub2] += 0.5*J_a

        ham22[sub1, sub1] += 0.5*J_a
        ham22[sub2, sub2] += 0.5*J_a

        ham11[sub1, sub2] += 0.5*J_b*phase
        ham11[sub2, sub1] += 0.5*cJ_b*cphase
        ham22[sub1, sub2] += 0.5*cJ_b*phase
        ham22[sub2, sub1] += 0.5*J_b*cphase

        ham21[sub1, sub2] += 0.5*J_c*phase
        ham21[sub2, sub1] += 0.5*J_c*cphase
        ham12[sub1, sub2] += 0.5*cJ_c*phase
        ham12[sub2, sub1] += 0.5*cJ_c*cphase

    for sublat in range(num_sub):

        Ri = R_l2g[sublat, :, :]
        fact = (hx*Ri[0, 2] + hy*Ri[1, 2] + hz*Ri[2, 2])
        ham11[sublat, sublat] += 0.5*fact
        ham22[sublat, sublat] += 0.5*fact

    ham = np.zeros((2*num_sub, 2*num_sub), dtype=complex)

    ham[:num_sub, :num_sub] = ham11
    ham[:num_sub, num_sub:] = ham12
    ham[num_sub:, :num_sub] = ham21
    ham[num_sub:, num_sub:] = ham22


    return ham

def eigensystem(q):
     
    hlsw = sw_hamiltonian(q)
    ham = A_mat @ hlsw
    eigval, eigvec = np.linalg.eig(ham)
    eigval = np.real(eigval)
#    imagpart = np.imag(eigval)
    idx = eigval.argsort()[::-1]
    eigval = eigval[idx]
    eigvec = eigvec[:, idx]
    tmp = eigvec.conj().T @ A_mat @ eigvec

    for kk in range(2*num_sub):
        eigvec[:, kk] = eigvec[:, kk]/np.sqrt(np.abs(tmp[kk, kk]))

    return 2.0*eigval, eigvec

def berry_curvature(q, Nsites):

    dq = 1.0/Nsites    
    hsw0 = sw_hamiltonian(q)
    epsilon, ubov = eigensystem(q)

    q1p = np.array([q[0]+dq, q[1], q[2]])
    q2p = np.array([q[0], q[1]+dq, q[2]])
    q3p = np.array([q[0], q[1], q[2]+dq])

    hsw1 = sw_hamiltonian(q1p)
    hsw2 = sw_hamiltonian(q2p)
    hsw3 = sw_hamiltonian(q3p)

    pH1 = (hsw1 - hsw0)/dq
    pH2 = (hsw2 - hsw0)/dq
    pH3 = (hsw3 - hsw0)/dq



    







