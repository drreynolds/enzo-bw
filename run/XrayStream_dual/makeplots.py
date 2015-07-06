# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
execfile("utilities.py")


# set the total number of snapshots
te = 100

# set the graphics output type
pictype = '.png'

# set 1D data arrays for each desired snapshot
tstep=0                                              # set output number
t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep)   # load data
Xr0 = sum(sum(Xr, axis=1), axis=0)                   # collapse into 1D array
tstep=10;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr1  = sum(sum(Xr, axis=1), axis=0)
tstep=20;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr2  = sum(sum(Xr, axis=1), axis=0)
tstep=30;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr3  = sum(sum(Xr, axis=1), axis=0)
tstep=40;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr4  = sum(sum(Xr, axis=1), axis=0)
tstep=50;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr5  = sum(sum(Xr, axis=1), axis=0)
tstep=60;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr6  = sum(sum(Xr, axis=1), axis=0)
tstep=70;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr7  = sum(sum(Xr, axis=1), axis=0)
tstep=80;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr8  = sum(sum(Xr, axis=1), axis=0)
tstep=90;  t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr9  = sum(sum(Xr, axis=1), axis=0)
tstep=100; t, vol, UV, Xr, xHI, xHII, Temp = load_vals(tstep);  Xr10 = sum(sum(Xr, axis=1), axis=0)

# set mesh
x = linspace(0.0, 1.0, Xr0.size)

# generate plot
figure()
plot(x,Xr0,x,Xr1,x,Xr2,x,Xr3,x,Xr4,x,Xr5,x,Xr6,x,Xr7,x,Xr8,x,Xr9,x,Xr10)
xlabel('$x$')
ylabel('X-ray Radiation Energy Density')
title('Propagation of X-ray Front')
grid()
savefig('propagation' + pictype)
