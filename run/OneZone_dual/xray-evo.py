from math import log10
import numpy as na
from matplotlib import rc
rc('font', size=18)
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager, FontProperties
font= FontProperties(size='small')
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()

redshift = 15.0
allF = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]  # radiation fluxes (erg/cm^2/s)
#allF = [1e-5, 1e-6, 1e-7]  # radiation fluxes (erg/cm^2/s)
Ex = 1e3   # photon energy (eV)

erg_eV = 8.61423e-5
G = 6.673e-8
H0 = 71.0
yH = 0.76
IH = 13.6
mH = 1.67262171e-24
h = 6.626e-27
kb = 1.38e-16
eV = 1.602e-12
yr = 3.1557e7

rhoH = 3 * (H0/3.086e19)**2 / (8 * na.pi * G) * yH * (1.0+redshift)**3

# Read cooling table (primordial abundances -- 2nd column)
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
    return 5.475e-14 * (E / 0.4298 - 1)**2 * (E / 0.4298)**(-4.0185) * \
        (1 + na.sqrt(E / 14.13))**(-2.963)

def secondary_ion(x):
    return 0.3908 * (1.0 - x**0.4092)**1.7592

def secondary_heat(x):
    return 0.9971 * (1.0 - (1.0 - x**0.2663)**1.3163)

plt.subplots_adjust(left=0.15, right=0.96, bottom=0.07, top=0.9, hspace=1e-3)
labels = {"x": r"$x_e$", "T": "Temperature [K]", "dx_dt": r"$\dot{x}_e$ [1/s]", 
          "dE_dt": r"$\dot{E}$ [erg/s]"}
colors = ["k", "b", "r", "g", "m", "c"]
styles = ["-", "--", "-.", ":", "-", "--"]

for ii, F in enumerate(allF):

    tfinal = 1e9 * yr
    dt0 = 1e6 * yr
    t = 0.0
    x = 1e-5
    T = 100.0
    Eth = kb * T
    results = {"time": [], "x": [], "T": [], "dx_dt": [], "dE_dt": []}
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
        t += dt
        x += dx_dt * dt
        Eth += (heat - cool) * dt
        results["time"].append(t)
        results["x"].append(x)
        results["T"].append(Eth/kb)
        results["dx_dt"].append(dx_dt)
        results["dE_dt"].append(heat-cool)

    for k in results.keys():
        results[k] = na.array(results[k])

    keys = results.keys()
    for i in range(4):
        plt.subplot(4,1,i+1)
        if i == 0:
            logF = log10(F)
            label = "$F_x = 10^{%.1f}$" % logF
        else:
            label = None
        p = plt.loglog(results["time"]/yr, results[keys[i]], c=colors[ii],
                       ls=styles[ii], label=label)
#        plt.xlim(1e5,1e9)
        plt.xlim(1e6,1e9)
        plt.ylabel(labels[keys[i]])
        if i == 3:
            plt.xlabel("Time [yr]")
        else:
            p[0].axes.xaxis.set_major_formatter(nullfmt)
        if i == 0:
            plt.legend(loc=8, prop=font, ncol=3, bbox_to_anchor=[0.5,1.0])

fig = plt.gcf()
fig.set_size_inches(8,12)
plt.savefig("x-evo.png")

