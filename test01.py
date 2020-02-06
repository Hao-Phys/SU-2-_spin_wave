#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : test01.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 01.16.2020
# Last Modified Date: 01.18.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import GLSW
from matplotlib import pyplot as plt

kk = np.arange(100)

lenk = len(kk)
omega = np.zeros((6, lenk))
omega_imag = np.zeros((6, lenk))

for flag in range(lenk):
    k3 = 4.0/lenk*flag
    k = np.array([0.0, 0.0, k3])
    ek, tmp, imagpart = GLSW.eigensystem(k)
    omega[:, flag] = ek[:6]
    omega_imag[:, flag] = imagpart[:6]

print(imagpart)
fig, ax = plt.subplots(2, 1)
for band in range(6):
    ax[0].plot(kk, omega[band, :])
    ax[1].plot(kk, omega_imag[band, :])

plt.ylim(0, 20)
plt.show()
