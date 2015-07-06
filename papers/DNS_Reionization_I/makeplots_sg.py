# matplotlib-based plotting script
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np
import matplotlib

# update font size
matplotlib.rcParams.update({'font.size': 16})

# SG test results
#          Grid    Tol   Error      Nt      Wall      RH       HYPRE    Chem
# -------------------------------------------------------------------------
sp_res = np.array([ 
        [16,    1.0e-2,   0.03207,    338,      83.5,      14.9,     13.5,     0.8],
        [16,    1.0e-3,   0.02473,    417,      82.6,      14.1,     11.9,     1.6],
        [16,    1.0e-4,   0.01275,   1116,     231.9,      55.0,     41.6,    11.2],
        [16,    1.0e-5,   0.00632,   8051,    1849.3,     282.5,    206.7,    62.3],
        [32,    1.0e-2,   0.04881,    338,     170.7,      94.9,     86.6,     5.2],
        [32,    1.0e-3,   0.03960,    406,     207.5,     115.4,     99.1,    12.8],
        [32,    1.0e-4,   0.01499,   1199,     588.8,     324.6,    254.0,    59.9],
        [32,    1.0e-5,   0.00603,   8933,    4733.8,    3014.7,   2145.1,   764.7],
        [64,    1.0e-2,   0.06337,    338,    1341.7,    1252.4,   1162.2,    61.6],
        [64,    1.0e-3,   0.06216,    348,    1451.0,    1359.5,   1183.8,   146.0],
        [64,    1.0e-4,   0.02110,   1022,    2941.9,    2678.7,   2142.1,   459.1],
        [64,    1.0e-5,   0.00635,   8594,   22639.0,   20291.0,  15237.6,  4422.9] ] )
# -------------------------------------------------------------------------
ex_res = np.array([ 
        [16,    1.0e-2,   0.02254,    338,     255.1,      10.3,      9.9,     0.0],
        [16,    1.0e-3,   0.01665,    410,     351.9,      19.9,     19.1,     0.0],
        [16,    1.0e-4,   0.00897,   1122,     641.4,      37.8,     35.8,     0.0],
        [16,    1.0e-5,   0.00584,   8035,    3726.3,     159.8,    148.5,     0.0],
        [32,    1.0e-2,   0.03320,    338,    1117.7,      87.5,     84.6,     0.0],
        [32,    1.0e-3,   0.02717,    397,    1202.6,     103.6,    100.2,     0.0],
        [32,    1.0e-4,   0.01131,   1155,    1928.0,     244.2,    234.5,     0.0],
        [32,    1.0e-5,   0.00571,   8907,   10611.0,    1698.8,   1618.2,     0.0],
        [64,    1.0e-2,   0.04249,    338,   11915.4,    1281.9,   1250.9,     0.0],
        [64,    1.0e-3,   0.04164,    349,   11649.0,    1224.0,   1193.5,     0.0],
        [64,    1.0e-4,   0.01495,   1023,   11739.8,    2232.4,   2165.5,     0.0],
        [64,    1.0e-5,   0.00591,   8595,   54755.5,   15538.2,  14945.3,     0.0] ] )
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
savefig('sg-error_split.pdf')

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
savefig('sg-error_enzo.pdf')


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
savefig('sg-runtime_split.pdf')

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
savefig('sg-runtime_enzo.pdf')


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
savefig('sg-efficiency_split.pdf')

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
savefig('sg-efficiency_enzo.pdf')


