#!/usr/bin/env python
# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
execfile("utilities.py")


# set the total number of snapshots
#te = 151
te = 60

# set the graphics output type
pictype = '.png'

# set some constants
Ngammadot = 1.0e50     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
kpc = 3.0857e21        # 1 kpc in cm
megayear = 3.15576e13  # 1 Myr in sec

# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (2, te+1), dtype=float);

 
# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    t, vol, rho, Eg, xHI, xHII, TotE, Vx, Vy, Vz, Temp, Pres, Mach  = load_all_vals(tstep)
    
    # compute volume element
    nx, ny, nz = Eg.shape
    dV = vol/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/megayear
    rdata[1][tstep] = radius
    
    # generate 2D plots at certain times
    if (tstep == 6) or (tstep == 20) or (tstep == 50):
        
        # set time label
        if (tstep == 6):
            Myr = '3'
        elif (tstep == 20):
            Myr = '10'
        else:
            Myr = '25'
        
        # set mesh
        x = linspace(0.0,1.0,nx)
        y = linspace(0.0,1.0,ny)
        X, Y = meshgrid(x,y)
        
        # xHI slice through z=0
        figure()
        sl = log10(xHI[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log HI fraction, t =' + Myr + ' Myr')
        savefig('HIcontour_' + Myr + 'Myr' + pictype)
        
        # xHII slice through z=0
        figure()
        sl = log10(xHII[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log HII fraction, t =' + Myr + ' Myr')
        savefig('HIIcontour_' + Myr + 'Myr' + pictype)
        
        # rho slice through z=0
        figure()
        sl = log10(rho[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Density, t =' + Myr + ' Myr')
        savefig('DensContour_' + Myr + 'Myr' + pictype)
        
        # Mach number slice through z=0
        figure()
        sl = log10(Mach[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Mach number, t =' + Myr + ' Myr')
        savefig('MachContour_' + Myr + 'Myr' + pictype)
        
        # Eg slice through z=0
        figure()
        sl = log10(Eg[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log radiation density, t =' + Myr + ' Myr')
        savefig('Econtour_' + Myr + 'Myr' + pictype)
        
        # Temp slice through z=0
        figure()
        sl = log10(Temp[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Temperature, t =' + Myr + ' Myr')
        savefig('TempContour_' + Myr + 'Myr' + pictype)
        
        # Pres slice through z=0
        figure()
        sl = log10(Pres[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Pressure, t =' + Myr + ' Myr')
        savefig('PressContour_' + Myr + 'Myr' + pictype)
        
        # spherically-averaged profiles for xHI, xHII, Temp, number density, pressure, mach number
        Nradii = nx*3/2
        Hradii = linspace(0.0,sqrt(3.0),Nradii)
        rad_idx = zeros( (nx,ny,nz), dtype=float)
        for k in range(nz):
            zloc = (k+0.5)/nz
            for j in range(ny):
                yloc = (j+0.5)/ny
                for i in range(nx):
                    xloc = (i+0.5)/nx
                    rad_idx[i][j][k] = max(0,floor(sqrt(xloc*xloc + yloc*yloc + zloc*zloc)/sqrt(3.0)*Nradii))
        Hcount = 1.0e-16*ones(Nradii)
        HIprof = zeros(Nradii, dtype=float)
        HIIprof = zeros(Nradii, dtype=float)
        Tprof = zeros(Nradii, dtype=float)
        nProf = zeros(Nradii, dtype=float)
        pProf = zeros(Nradii, dtype=float)
        mProf = zeros(Nradii, dtype=float)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    idx = rad_idx[i][j][k]
                    HIprof[idx] += xHI[i][j][k]
                    HIIprof[idx] += xHII[i][j][k]
                    Tprof[idx] += Temp[i][j][k]
                    nProf[idx] += rho[i][j][k]/mp
                    pProf[idx] += Pres[i][j][k]
                    mProf[idx] += Mach[i][j][k]
                    Hcount[idx] += 1
        HIprof  = log10(HIprof/Hcount)
        HIIprof = log10(HIIprof/Hcount)
        Tprof   = log10(Tprof/Hcount)
        nProf   = log10(nProf/Hcount)
        pProf   = log10(pProf/Hcount)
        mProf   = log10(mProf/Hcount)
        
        # chemistry profiles
        figure()
        plot(Hradii,HIprof,'b-',Hradii,HIIprof,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, t =' + Myr + ' Myr')
        legend( ('xHI','xHII') )
        axis([ 0.0, 1.2, -5.0, 1.0 ])
        savefig('profiles_' + Myr + 'Myr' + pictype)
        
        # Temperature profile
        figure()
        plot(Hradii,Tprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(T) [K]')
        title('Temperature Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, 3.0, 5.0 ])
        savefig('TempProfile_' + Myr + 'Myr' + pictype)
        
        # Number density profile
        figure()
        plot(Hradii,nProf)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(n) [cm$^{-3}$]')
        title('Number density Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, -4.0, 1.0 ])
        savefig('nProfile_' + Myr + 'Myr' + pictype)
        
        # Pressure profile
        figure()
        plot(Hradii,pProf)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(p) [g cm$^{-1}$ s${-2}$]')
        title('Pressure Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, -17.0, -10.0 ])
        savefig('PressProfile_' + Myr + 'Myr' + pictype)
        
        # Mach number profile
        figure()
        plot(Hradii,mProf)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(Mach)')
        title('Mach Number Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, -5.0, 1.0 ])
        savefig('MachProfile_' + Myr + 'Myr' + pictype)


# I-front radius
figure()
plot(rdata[0],rdata[1]/kpc,'b-')
xlabel('$t$ [Myr]')
ylabel('$r_I$ [kpc]')
title('Propagation of HII Region')
axis([ 0.0, 30.0, -0.1, 0.8 ])
grid()
savefig('rad_vs_time' + pictype)
