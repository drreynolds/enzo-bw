# matplotlib-based plotting script for Iliev et al. test #1
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *

# define some helpful functions
def get_params(file):
    """Returns t, vol, gamma, dUnit, tUnit, lUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("InitialTime" in text):
            tval = float(text[len(text)-1])
        elif ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
        elif ("Gamma" in text):
            gamma = float(text[len(text)-1])
        elif ("DomainLeftEdge" in text):
            xL = float(text[len(text)-3])
            yL = float(text[len(text)-2])
            zL = float(text[len(text)-1])
        elif ("DomainRightEdge" in text):
            xR = float(text[len(text)-3])
            yR = float(text[len(text)-2])
            zR = float(text[len(text)-1])
    vol = (xR-xL)*(yR-yL)*(zR-zL)*lUnit*lUnit*lUnit
    tval = tval*tUnit
    return [tval, vol, gamma, dUnit, tUnit, lUnit]

########
def get_params_cosmology(file):
    """Returns z0, z, xR, t0, H0, dUnit, tUnit, lUnit from a parameter files"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("CosmologyInitialRedshift" in text):
            z0 = float(text[len(text)-1])
        elif ("CosmologyCurrentRedshift" in text):
            z = float(text[len(text)-1])
        elif ("InitialTime" in text):
            t0 = float(text[len(text)-1])
        elif ("CosmologyHubbleConstantNow" in text):
            H0 = float(text[len(text)-1])
    f.close()
    f = open(file + '.rtmodule' )
    for line in f:
        text = shlex.split(line)
        if ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
    f.close()
    xR = lUnit
    t0 *= tUnit
    H0 *= 100*1e5/3.0857e24  # H0 units given 100km/s/Mpc, convert to 1/s 
    return [z0, z, xR, t0, H0, dUnit, tUnit, lUnit]
######

def load_vals(tdump):
    """Returns t, vol, UV, Xray, xHI, xHII, xHeI, xHeII, xHeIII, Temp, E, HIkph, photogamma from a given data dump"""
    import h5py
    import numpy as np
    # get general dump information
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    tval, vol, gamma, dUnit, tUnit, lUnit = get_params(pfile)
    # get parallelism information
    dims, nprocs, gextent = get_parallel_decomp(tdump)
    # allocate output arrays
    UV_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    Xr_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHI_tot  = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHII_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHeI_tot  = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHeII_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHeIII_tot  = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    Temp_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    E_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    HIkph_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    photogamma_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    # set some constants
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    # iterate over grids, inserting their values in the appropriate locations
    for i in range(nprocs):
        pstring = repr(i+1).zfill(8)
        hfile = pfile + '.grid' + pstring
        # get data arrays
        f = h5py.File(hfile,'r')
        UV = np.array(f.get('UV_Radiation')) * (dUnit*lUnit*lUnit/tUnit/tUnit)
        Xr = np.array(f.get('Xray_Radiation')) * (dUnit*lUnit*lUnit/tUnit/tUnit)
        energy = np.array(f.get('Total_Energy')) * (lUnit*lUnit/tUnit/tUnit)
        HI = np.array(f.get('HI_Density')) * dUnit
        HII = np.array(f.get('HII_Density')) * dUnit
        HeI = np.array(f.get('HeI_Density')) * dUnit
        HeII = np.array(f.get('HeII_Density')) * dUnit
        HeIII = np.array(f.get('HeIII_Density')) * dUnit
        rho = np.array(f.get('Density')) * dUnit
        HIkph = f.get('HI_kph')
        photogamma = f.get('PhotoGamma')
        # compute derived fields
        H = HI + HII + (1.e-30 * dUnit)
        He = HeI + HeII + HeIII + (1.e-30 * dUnit)
        HIfrac = np.divide(HI,H)
        HIIfrac = np.divide(HII,H)
        HeIfrac = np.divide(HeI,He)
        HeIIfrac = np.divide(HeII,He)
        HeIIIfrac = np.divide(HeIII,He)
        ne = HII + HeII/4.0 + HeIII/2.0
        num_dens = 0.25*He + H + ne
        mu = np.divide(rho, num_dens)
        Temp = np.multiply(energy, mu) * (gamma-1.0)*mp/kb
        # insert into output arrays
        UV_tot[gextent[i][0][0]:gextent[i][0][1]+1,
               gextent[i][1][0]:gextent[i][1][1]+1,
               gextent[i][2][0]:gextent[i][2][1]+1] = UV
        Xr_tot[gextent[i][0][0]:gextent[i][0][1]+1,
               gextent[i][1][0]:gextent[i][1][1]+1,
               gextent[i][2][0]:gextent[i][2][1]+1] = Xr
        xHI_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                gextent[i][1][0]:gextent[i][1][1]+1,
                gextent[i][2][0]:gextent[i][2][1]+1] = HIfrac
        xHII_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                 gextent[i][1][0]:gextent[i][1][1]+1,
                 gextent[i][2][0]:gextent[i][2][1]+1] = HIIfrac
        xHeI_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                 gextent[i][1][0]:gextent[i][1][1]+1,
                 gextent[i][2][0]:gextent[i][2][1]+1] = HeIfrac
        xHeII_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                  gextent[i][1][0]:gextent[i][1][1]+1,
                  gextent[i][2][0]:gextent[i][2][1]+1] = HeIIfrac
        xHeIII_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                   gextent[i][1][0]:gextent[i][1][1]+1,
                   gextent[i][2][0]:gextent[i][2][1]+1] = HeIIIfrac
        Temp_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                 gextent[i][1][0]:gextent[i][1][1]+1,
                 gextent[i][2][0]:gextent[i][2][1]+1] = Temp
        E_tot[gextent[i][0][0]:gextent[i][0][1]+1,
              gextent[i][1][0]:gextent[i][1][1]+1,
              gextent[i][2][0]:gextent[i][2][1]+1] = energy
        HIkph_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                  gextent[i][1][0]:gextent[i][1][1]+1,
                  gextent[i][2][0]:gextent[i][2][1]+1] = HIkph
        photogamma_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                       gextent[i][1][0]:gextent[i][1][1]+1,
                       gextent[i][2][0]:gextent[i][2][1]+1] = photogamma
    # finished, return with results 
    return [tval, vol, UV_tot, Xr_tot, xHI_tot, xHII_tot, xHeI_tot, xHeII_tot, xHeIII_tot, Temp_tot, E_tot, HIkph_tot, photogamma_tot]

########
def load_vals_cosmology(tdump):
    """Returns t, z, xR, nH, UV, Xray, xHI, xHII from a given data dump"""
    import h5py
    import numpy as np
    # get general dump information
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pstring = repr(1).zfill(8)
    hfile = pfile + '.grid' + pstring
    z0, z, xR, tval, H0, dUnit, tUnit, lUnit = get_params_cosmology(pfile)
    # get parallelism information
    dims, nprocs, gextent = get_parallel_decomp(tdump)
    # allocate output arrays
    UV_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    Xr_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHI_tot  = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    xHII_tot = np.zeros( (dims[0],dims[1],dims[2]), dtype=float)
    # set some constants
    mp = 1.67262171e-24
    # iterate over grids, inserting their values in the appropriate locations
    for i in range(nprocs):
        pstring = repr(i+1).zfill(8)
        hfile = pfile + '.grid' + pstring
        f = h5py.File(hfile,'r')
        UV = f.get('UV_Radiation')
        Xr = f.get('Xray_Radiation')
        HI = f.get('HI_Density')
        HII = f.get('HII_Density')
        rho = f.get('Density')
        # add floor values for happy numerics
        HI  = np.add(HI,  1.0e-10)
        HII = np.add(HII, 1.0e-10)
        UV  = np.add(UV,  1.0e-30)
        Xr  = np.add(Xr,  1.0e-30)
        HIfrac = np.divide(HI,rho)
        HIIfrac = np.divide(HII,rho)
        UV = np.add(UV, 1.0e-30)
        UV = np.multiply(UV,dUnit*lUnit*lUnit/tUnit/tUnit)
        Xr = np.add(Xr, 1.0e-30)
        Xr = np.multiply(Xr,dUnit*lUnit*lUnit/tUnit/tUnit)
        mh = 1.67262171e-24        
        nH = rho[0,0,0]*dUnit/mh
        # insert into output arrays
        UV_tot[gextent[i][0][0]:gextent[i][0][1]+1,
               gextent[i][1][0]:gextent[i][1][1]+1,
               gextent[i][2][0]:gextent[i][2][1]+1] = UV
        Xr_tot[gextent[i][0][0]:gextent[i][0][1]+1,
               gextent[i][1][0]:gextent[i][1][1]+1,
               gextent[i][2][0]:gextent[i][2][1]+1] = Xr
        xHI_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                gextent[i][1][0]:gextent[i][1][1]+1,
                gextent[i][2][0]:gextent[i][2][1]+1] = HIfrac
        xHII_tot[gextent[i][0][0]:gextent[i][0][1]+1,
                 gextent[i][1][0]:gextent[i][1][1]+1,
                 gextent[i][2][0]:gextent[i][2][1]+1] = HIIfrac
    # finished, return with results 
    return [tval, z, xR, nH, UV_tot, Xr_tot, xHI_tot, xHII_tot]

########
def get_parallel_decomp(tdump):
    """Returns dims, nprocs, gridextents from a given data dump"""
    import shlex
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.hierarchy'
    nprocs = 0
    # get the total number of grids
    f = open(hfile)
    for line in f:
        text = shlex.split(line)
        if ("Grid" in text):
            nprocs += 1
    f.close()
    # allocate the output arrays
    gridextents = np.zeros( (nprocs,3,2), dtype=int)
    # loop over the grids, filling in data
    igrid = 0
    GridEdges = np.zeros( (nprocs,3,2), dtype=float)
    GridDims  = np.zeros( (nprocs,3), dtype=int)
    GridStart = np.zeros( 3, dtype=int)
    GridEnd   = np.zeros( 3, dtype=int)
    f = open(hfile)
    for line in f:
        text = shlex.split(line)
        if ("GridStartIndex" in text):
            GridStart[0] = int(text[len(text)-3])
            GridStart[1] = int(text[len(text)-2])
            GridStart[2] = int(text[len(text)-1])
        elif ("GridEndIndex" in text):
            GridEnd[0] = int(text[len(text)-3])
            GridEnd[1] = int(text[len(text)-2])
            GridEnd[2] = int(text[len(text)-1])
        elif ("GridLeftEdge" in text):
            GridEdges[igrid,0,0] = float(text[len(text)-3])
            GridEdges[igrid,1,0] = float(text[len(text)-2])
            GridEdges[igrid,2,0] = float(text[len(text)-1])
        elif ("GridRightEdge" in text):
            GridEdges[igrid,0,1] = float(text[len(text)-3])
            GridEdges[igrid,1,1] = float(text[len(text)-2])
            GridEdges[igrid,2,1] = float(text[len(text)-1])
        elif ("NumberOfParticles" in text):
            GridDims[igrid,0] = GridEnd[0]-GridStart[0]+1
            GridDims[igrid,1] = GridEnd[1]-GridStart[1]+1
            GridDims[igrid,2] = GridEnd[2]-GridStart[2]+1
            igrid += 1
    f.close()
    # get spatial mesh size (uniform, so all grids share these)
    dx = float((GridEdges[0,0,1] - GridEdges[0,0,0])/GridDims[0,0])
    dy = float((GridEdges[0,1,1] - GridEdges[0,1,0])/GridDims[0,1])
    dz = float((GridEdges[0,2,1] - GridEdges[0,2,0])/GridDims[0,2])
    # loop over the grids, determining where they sit in the overall layout
    GridLocs = np.zeros( (nprocs,3), dtype=int)
    for i in range(nprocs):
        for j in range(nprocs):
            if ((j != i) & (GridEdges[i,0,0]-GridEdges[j,0,0]>0.5*dx) 
                         & (abs(GridEdges[i,1,0]-GridEdges[j,1,0])<0.5*dy)
                         & (abs(GridEdges[i,2,0]-GridEdges[j,2,0])<0.5*dz)):
                GridLocs[i,0] += 1
                gridextents[i,0,0] += GridDims[j,0]
            if ((j != i) & (abs(GridEdges[i,0,0]-GridEdges[j,0,0])<0.5*dx) 
                         & (GridEdges[i,1,0]-GridEdges[j,1,0]>0.5*dy)
                         & (abs(GridEdges[i,2,0]-GridEdges[j,2,0])<0.5*dz)):
                GridLocs[i,1] += 1
                gridextents[i,1,0] += GridDims[j,1]
            if ((j != i) & (abs(GridEdges[i,0,0]-GridEdges[j,0,0])<0.5*dx) 
                         & (abs(GridEdges[i,1,0]-GridEdges[j,1,0])<0.5*dy)
                         & (GridEdges[i,2,0]-GridEdges[j,2,0]>0.5*dz)):
                GridLocs[i,2] += 1
                gridextents[i,2,0] += GridDims[j,2]
        gridextents[i,0,1] = gridextents[i,0,0] + GridDims[i,0] - 1
        gridextents[i,1,1] = gridextents[i,1,0] + GridDims[i,1] - 1
        gridextents[i,2,1] = gridextents[i,2,0] + GridDims[i,2] - 1
        
    # compute the parallel layout
    xprocs = 1
    yprocs = 1
    zprocs = 1
    for i in range(nprocs):
        xprocs = max(xprocs,GridLocs[i,0]+1)
        yprocs = max(yprocs,GridLocs[i,1]+1)
        zprocs = max(zprocs,GridLocs[i,2]+1)
    # compute the total dimensions of the grid
    dims = np.zeros( 3, dtype=int)
    for j in range(xprocs):
        for i in range(nprocs):
            if ((GridLocs[i,0]==j) & (GridLocs[i,1]==0) & (GridLocs[i,2]==0)):
                dims[0] += GridDims[i,0]
    for j in range(yprocs):
        for i in range(nprocs):
            if ((GridLocs[i,1]==j) & (GridLocs[i,0]==0) & (GridLocs[i,2]==0)):
                dims[1] += GridDims[i,1]
    for j in range(zprocs):
        for i in range(nprocs):
            if ((GridLocs[i,2]==j) & (GridLocs[i,0]==0) & (GridLocs[i,1]==0)):
                dims[2] += GridDims[i,2]
    # account for the fact that HDF5 reverses grid indices
    dimsT = [dims[2], dims[1], dims[0]];
    gridextentsT = np.zeros( (nprocs,3,2), dtype=int)
    for i in range(nprocs):
        gridextentsT[i,2,:] = gridextents[i,0,:]
        gridextentsT[i,1,:] = gridextents[i,1,:]
        gridextentsT[i,0,:] = gridextents[i,2,:]

    return [dimsT, nprocs, gridextentsT]


########## end of file ##########
