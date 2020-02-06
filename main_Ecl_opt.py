#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : main_Ecl_opt.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.05.2020
# Last Modified Date: 02.06.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>


import numpy as np
import time
import Ecl_opt_obj as Eobj
import cofig as cf
import scipy.optimize as sciopt

num_sub = cf.num_sub
dirpath = cf.path
spin = cf.spin

x_global_opt = np.zeros(2*num_sub)
e_global_opt = 1000.0

x0 = np.zeros(2*num_sub)
lb0 = np.zeros(2*num_sub)
ub0 = np.zeros(2*num_sub)

lb0 = -0.001

for flag in range(num_sub):
    ub0[flag*2] = np.pi + 0.001
    ub0[flag*2+1] = 2.0*np.pi + 0.001

bb = sciopt.Bounds(lb0, ub0)

N_rand = 100

for flag1 in range(N_rand):
    start_time = time.time()
    np.random.seed(seed=None)

    for flag2 in range(num_sub):
        x0[2*flag2] = np.pi * np.random.rand(1) 
        x0[2*flag2+1] = 2.0 * np.pi * np.random.rand(1)

    res = sciopt.minimize(Eobj.Energ_cl, x0, method='L-BFGS-B', 
                        bounds = bb,
                        options={'ftol': 2.220446049250313e-09, 'disp': 2, 'maxfun': 20000})
    
    x_local_opt = res.x
    e_local_opt = res.fun
    return_flag = res.status
# =============================================================================
#     opt = nlopt.opt(nlopt.LN_NELDERMEAD, 4*num_sub)
#     opt.set_lower_bounds(lb0)
#     opt.set_upper_bounds(ub0)
#     opt.set_min_objective(Eobj.Energy_cl)
#     opt.set_xtol_rel(1e-4)
#     opt.set_maxeval = 1000
#     x_local_opt = opt.optimize(x0)
#     e_local_opt = opt.last_optimum_value()
#     return_flag = opt.last_optimize_result()
#     # nlopt return values
#     # 1: generic success; 2: stopval reached;
#     # 3: ftol_rel (or abs) reached; 4. x_tol reached
#     # 5: maxevl reached; 6: maxtime reached
#     # negative: error, check nlopt document
# =============================================================================

    if e_local_opt < e_global_opt:
        print('find new local minima, updates!')
        x_global_opt = x_local_opt
        e_global_opt = e_local_opt

    print("finishing the %s th random experiment!" % flag1)
    print("time elpase in this experiment is %s seconds" % (time.time() - start_time))
    print("the return flag is %s" % return_flag)


mf_vals = np.zeros((num_sub, 3))

for flag3 in range(num_sub):
    theta  = x_global_opt[flag3*2]
    phi    = x_global_opt[flag3*2+1]

    Si = spin*np.array([np.sin(theta)*np.cos(phi), \
                        np.sin(theta)*np.sin(phi), \
                        np.cos(theta)])

    mf_vals[flag3, 0] = Si[0]
    mf_vals[flag3, 1] = Si[1]
    mf_vals[flag3, 2] = Si[2]

fname = dirpath + 'opt_angles.txt'
np.savetxt(fname, x_global_opt)
fname1 = dirpath + 'mfvals.txt'
np.savetxt(fname1, mf_vals)
