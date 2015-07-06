#!/usr/bin/env python
# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
execfile("utilities.py")

from math import log10
import numpy as na
from matplotlib import rc
rc('font', size=18)
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager, FontProperties
font= FontProperties(size='small')
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()

# set the incident flux
F = 1e-2
#F = 1e-3
#F = 1e-4
#F = 1e-5
#F = 1e-6
#F = 1e-7

# set the total number of computed snapshots
te = 1000

# set parameters for reference solution
redshift = 15.0
Ex = 1e3   # photon energy (eV)
erg_eV = 8.616975441619992e-05
G = 6.673e-8
H0 = 71.0
#yH = 0.76
yH = 1.0
IH = 13.6
h = 6.6260693e-27
mH = 1.67262171e-24
kb = 1.3806504e-16
eV = 1.60217653e-12
yr = 3.1557e7
rhoH = 3 * (H0/3.086e19)**2 / (8 * na.pi * G) * yH * (1.0+redshift)**3
rho = 3 * (H0/3.086e19)**2 / (8 * na.pi * G) * (1.0+redshift)**3
print 'rho = ', rho
print 'Hydrogen fraction = ', rhoH/rho

# read cooling table (primordial abundances -- 2nd column)
cool_table = na.loadtxt("zcool_sd93.dat")[:,0:2]

def CH(T):
    lnT = na.log(T*erg_eV)
    coeff = [-32.71396786, 13.536556, -5.73932875, 1.56315498, -0.2877056,
              3.48255977e-2, -2.63197617e-3, 1.11954395e-4, -2.03914985e-6]
    a = 0.0
    for i in range(len(coeff)):
        a += coeff[i] * lnT**i
    return na.exp(a)

def alphaB(T):
    return 2.59e-13 * (T/1e4)**(-0.7)

def sigmaH(E):
    nuscaled = E / 13.6
    eps = na.sqrt(nuscaled - 1.0)
    return (6.30e-18) * nuscaled**(-4.0) * \
	exp(4.0-4.0*na.arctan(eps)/eps)/(1.0-exp(-2.0*na.pi/eps))
#    return 5.475e-14 * (E / 0.4298 - 1)**2 * (E / 0.4298)**(-4.0185) * \
#        (1 + na.sqrt(E / 14.13))**(-2.963)

def secondary_ion(x):
    return 0.3908 * (1.0 - x**0.4092)**1.7592

def secondary_heat(x):
    return 0.9971 * (1.0 - (1.0 - x**0.2663)**1.3163)

plt.subplots_adjust(left=0.15, right=0.96, bottom=0.07, top=0.9, hspace=1e-3)
labels = {"x": r"$x_e$", "T": "Temperature [K]"}

# compute reference solution
tfinal = 1e9 * yr
dt0 = 1e6 * yr
t = 0.0
x = 1e-5
T = 100.0
Eth = kb * T
results = {"time": [], "x": [], "T": []}
tout = t+dt0
while t < tfinal:
    ne = rhoH*x/mH
    T = Eth / kb
    logT = log10(T)
    kph = secondary_ion(x) * (F / (Ex*eV)) * sigmaH(Ex) * (Ex/IH)
    heat = secondary_heat(x) * (F / (Ex*eV)) * sigmaH(Ex) * (Ex*eV)
    if T < 1e4:
        cool = 0.0
    else:
        cool_idx = na.searchsorted(cool_table[:,0], logT)
        interp_factor = (logT - cool_table[cool_idx,0]) / \
            (cool_table[cool_idx+1,0] - cool_table[cool_idx,0])
        cool = cool_table[cool_idx,1] + \
            interp_factor * (cool_table[cool_idx+1,1] - cool_table[cool_idx,1])
        cool = 10.0**cool
    dx_dt = (1.0 - x) * (kph - ne*CH(T)) - x*ne*alphaB(T)
    dt = min(dt0, 0.01 * abs(x/(dx_dt)),
             0.01 * abs(Eth / (heat - cool)))
    dt = min(dt, tout - t)
    t += dt
    x += dx_dt * dt
    Eth += (heat - cool) * dt
    if ((t - tout)/tout > -1e-8):
        tout = t + dt0
        results["time"].append(t)
        results["x"].append(x)
        results["T"].append(Eth/kb)

for k in results.keys():
    results[k] = na.array(results[k])


# load computed results
Nresults = {"time": [], "x": [], "xHI": [], "T": [], "E": [], "Xray": [], "HIkph": [], "photogamma": []}
for it in range(te+1):
    t, vol, UV, Xr, xHI, xHII, xHeI, xHeII, xHeIII, Temp, E, HIkph, photogamma = load_vals(it)
    Nresults["time"].append(t)
    Nresults["x"].append(na.average(xHII))
    Nresults["xHI"].append(na.average(xHI))
    Nresults["T"].append(na.average(Temp))
    Nresults["E"].append(na.average(E))
    Nresults["Xray"].append(na.average(Xr))
    Nresults["HIkph"].append(na.average(HIkph))
    Nresults["photogamma"].append(na.average(photogamma))

for k in Nresults.keys():
    Nresults[k] = na.array(Nresults[k])


# generate plots
ylims = ( (1e-5, 1), (1e2, 1e5) );
keys = results.keys()
for i in range(2):
    plt.subplot(2,1,i+1)
    p = plt.loglog(results["time"]/yr, results[keys[i]], "b-", label="analytical")
    p = plt.loglog(Nresults["time"]/yr, Nresults[keys[i]], "r-", label="computed")
    #plt.xlim(1e5,1e9)
    plt.xlim(dt0/yr,1e9)
#    plt.ylim(ylims[i][0], ylims[i][1])
    plt.ylabel(labels[keys[i]])
    plt.legend(prop=font, loc=2)
    if i == 0:
        logF = log10(F)
        lab = "$F_x = 10^{%.1f}$" % logF
        plt.title(lab)
        p[0].axes.xaxis.set_major_formatter(nullfmt)
    else:
        plt.xlabel("Time [yr]")

fig = plt.gcf()
fig.set_size_inches(8,10)
plt.savefig("x-evo_comparisons.png")



plt.figure()
plt.subplot(5,1,1)
plt.loglog(Nresults["time"]/yr, Nresults["Xray"])
plt.xlim(dt0/yr,1e9)
plt.ylabel("Xray")
plt.title("Other fields")

plt.subplot(5,1,2)
plt.loglog(Nresults["time"]/yr, Nresults["HIkph"])
plt.xlim(dt0/yr,1e9)
plt.ylabel("HI_kph")

plt.subplot(5,1,3)
plt.loglog(Nresults["time"]/yr, Nresults["photogamma"])
plt.xlim(dt0/yr,1e9)
plt.ylabel("photogamma")

plt.subplot(5,1,4)
plt.loglog(Nresults["time"]/yr, Nresults["E"])
plt.xlim(dt0/yr,1e9)
plt.ylabel("E")

plt.subplot(5,1,5)
plt.plot(Nresults["time"]/yr, Nresults["x"] + Nresults["xHI"])
plt.xlim(dt0/yr,1e9)
plt.ylim(0,2)
plt.ylabel("x + (1-x)")

fig = plt.gcf()
fig.set_size_inches(8,15)
plt.savefig("x-fields.png")
