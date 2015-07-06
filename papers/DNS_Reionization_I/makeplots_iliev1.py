# matplotlib-based plotting script
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np
import matplotlib

# update font size
matplotlib.rcParams.update({'font.size': 16})

# iliev test #1 results
#          Grid    Tol   Error      Nt      Wall      RH       HYPRE    Chem
# -------------------------------------------------------------------------
sp_res = np.array([ 
        [16,   1.0e-2,  0.06240,   383,     19.2,     14.7,    12.7,     1.5],
        [16,   1.0e-3,  0.01901,  1065,     47.4,     38.7,    30.3,     7.0],
        [16,   1.0e-4,  0.00225,  7035,    373.0,    330.4,   220.2,    98.9],
        [16,   1.0e-5,  0.00304, 67533,   2300.6,   1895.2,  1208.7,   600.6],
        [32,   1.0e-2,  0.08686,   340,    109.4,    101.0,    87.3,    11.1],
        [32,   1.0e-3,  0.03060,   848,    282.0,    267.3,   209.8,    50.8],
        [32,   1.0e-4,  0.00610,  5992,   1722.5,   1642.4,  1208.1,   388.0],
        [32,   1.0e-5,  0.00129, 56931,  18756.9,  17930.2, 11728.6,  5551.6],
        [64,   1.0e-2,  0.08240,   431,   1753.9,   1710.7,  1520.2,   154.9],
        [64,   1.0e-3,  0.03652,   804,   3806.1,   3732.3,  3011.9,   648.2],
        [64,   1.0e-4,  0.00780,  5335,  21383.5,  20915.0, 15997.3,  4439.9],
        [64,   1.0e-5,  0.00168, 49572, 105700.5, 102716.0, 73277.2, 26584.8] ] )
# -------------------------------------------------------------------------
ex_res = np.array([ 
        [16,   1.0e-2,  0.04274,   377,     19.1,     12.0,    11.6,     0.0],
        [16,   1.0e-3,  0.01086,  1023,     42.9,     29.7,    28.4,     0.0],
        [16,   1.0e-4,  0.00542,  6951,    224.4,    149.3,   141.5,     0.0],
        [16,   1.0e-5,  0.00795, 67703,   2018.1,   1212.5,  1125.5,     0.0],
        [32,   1.0e-2,  0.06367,   334,    126.5,    100.9,    98.0,     0.0],
        [32,   1.0e-3,  0.02183,   783,    240.0,    194.2,   188.0,     0.0],
        [32,   1.0e-4,  0.00328,  5927,   1535.3,   1231.3,  1183.5,     0.0],
        [32,   1.0e-5,  0.00528, 57017,  11553.8,   8697.9,  8172.3,     0.0],
        [64,   1.0e-2,  0.05286,   429,   1810.8,   1559.6,  1523.4,     0.0],
        [64,   1.0e-3,  0.02164,   786,   3293.4,   2856.1,  2784.9,     0.0],
        [64,   1.0e-4,  0.00366,  5334,  18284.7,  15680.4, 15207.8,     0.0],
        [64,   1.0e-5,  0.00442, 49646,  92276.7,  75851.6, 72896.0,     0.0] ] )
# -------------------------------------------------------------------------
grids = ('$16^3$ mesh','$32^3$ mesh','$64^3$ mesh')
ONE = ( 1.0, 1.0, 1.0, 1.0 )

# extract desired arrays from data
nx16_sp = sp_res[0:4,:]
nx32_sp = sp_res[4:8,:]
nx64_sp = sp_res[8:12,:]

nx16_ex = ex_res[0:4,:]
nx32_ex = ex_res[4:8,:]
nx64_ex = ex_res[8:12,:]


# plot error vs average time step size
fig = figure()
loglog(np.divide(ONE,nx16_sp[:,3]), nx16_sp[:,2], 'r-')
loglog(np.divide(ONE,nx32_sp[:,3]), nx32_sp[:,2], 'b--')
loglog(np.divide(ONE,nx64_sp[:,3]), nx64_sp[:,2], 'k-.')
xlabel('$\Delta t_{avg}$', fontsize=20)
ylabel('Error', fontsize=20)
title('Error vs Average Time Step, Split Solver')
legend(grids, loc='upper left', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-error_split.pdf')

fig = figure()
loglog(np.divide(ONE,nx16_ex[:,3]), nx16_ex[:,2], 'r-')
loglog(np.divide(ONE,nx32_ex[:,3]), nx32_ex[:,2], 'b--')
loglog(np.divide(ONE,nx64_ex[:,3]), nx64_ex[:,2], 'k-.')
xlabel('$\Delta t_{avg}$', fontsize=20)
ylabel('Error', fontsize=20)
title('Error vs Average Time Step')
legend(grids, loc='upper left', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-error_enzo.pdf')


# plot runtime vs average time step size
fig = figure()
loglog(np.divide(ONE,nx16_sp[:,3]), nx16_sp[:,4], 'r-')
loglog(np.divide(ONE,nx32_sp[:,3]), nx32_sp[:,4], 'b--')
loglog(np.divide(ONE,nx64_sp[:,3]), nx64_sp[:,4], 'k-.')
xlabel('$\Delta t_{avg}$', fontsize=20)
ylabel('Runtime', fontsize=20)
title('Runtime vs Average Time Step, Split Solver')
legend(grids, loc='upper right', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-runtime_split.pdf')

fig = figure()
loglog(np.divide(ONE,nx16_ex[:,3]), nx16_ex[:,4], 'r-')
loglog(np.divide(ONE,nx32_ex[:,3]), nx32_ex[:,4], 'b--')
loglog(np.divide(ONE,nx64_ex[:,3]), nx64_ex[:,4], 'k-.')
xlabel('$\Delta t_{avg}$', fontsize=20)
ylabel('Runtime', fontsize=20)
title('Runtime vs Average Time Step')
legend(grids, loc='upper right', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-runtime_enzo.pdf')


# plot error vs runtime
fig = figure()
loglog(nx16_sp[:,4], nx16_sp[:,2], 'r-')
loglog(nx32_sp[:,4], nx32_sp[:,2], 'b--')
loglog(nx64_sp[:,4], nx64_sp[:,2], 'k-.')
xlabel('Runtime', fontsize=20)
ylabel('Error', fontsize=20)
title('Error vs Runtime, Split Solver')
legend(grids, loc='upper right', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-efficiency_split.pdf')

fig = figure()
loglog(nx16_ex[:,4], nx16_ex[:,2], 'r-')
loglog(nx32_ex[:,4], nx32_ex[:,2], 'b--')
loglog(nx64_ex[:,4], nx64_ex[:,2], 'k-.')
xlabel('Runtime', fontsize=20)
ylabel('Error', fontsize=20)
title('Error vs Runtime')
legend(grids, loc='upper right', shadow=True)
#axis((1, 100, 1, 8))
grid()
fig.subplots_adjust(bottom=0.15)
savefig('i1-efficiency_enzo.pdf')


