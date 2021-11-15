# This script calculates differential field distorsions as a function of zenith angle and wavelength
# for MegaMapper optical designs. It takes as input ZEMAX spot centroid files provided by Robert Barkhouser
# it calculates differential field distorsion fields on field angle units, it creates diagnostic plots,
# and it creates a differential field distorsion object that is pickled, to be used by the mmsim.distsim simulations.

import numpy as np
import astropy.io.ascii as ascii
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mmsim import mkresplot, mkresplot1d
from scipy.interpolate import LinearNDInterpolator
import pickle

# Define the ZEMAX run to be used
auxrun="1536revA"
dir='./zemax/20210817/run1536revA/'
files=glob.glob(dir+'*.txt')

# Make array of wavelengths and zenithal distances 
auxwave=np.empty_like(files)
auxzd=np.empty_like(files)
for i in range(len(files)):
    auxwave[i]=re.split('[_/.]', files[i])[-2]
    auxzd[i]=re.split('[_/.]', files[i])[-3]
wavearr=np.unique(auxwave)
zdarr=np.unique(auxzd)
nwave=len(wavearr)
nzd=len(zdarr)
print("Wavelengths:", wavearr)
print("Zenith Distances:", zdarr)

# Get number of spots in the ZEMAX spot centroid files
infile=dir+'run'+auxrun+'_'+zdarr[0]+'_'+wavearr[0]+'.txt'
tab0=ascii.read(infile, encoding='utf-16')
nspots=len(tab0)

# Crate 3D arrays for all data
fld_x=np.zeros((nwave, nzd, nspots))   # field angle
fld_y=np.zeros((nwave, nzd, nspots))
wave=np.zeros((nwave, nzd, nspots))    # wavelength
chf_x=np.zeros((nwave, nzd, nspots))
chf_y=np.zeros((nwave, nzd, nspots))
cen_x=np.zeros((nwave, nzd, nspots))   # spot centroids
cen_y=np.zeros((nwave, nzd, nspots))
rms=np.zeros((nwave, nzd, nspots))
n_rays=np.zeros((nwave, nzd, nspots))
zd=np.zeros((nwave, nzd, nspots))      # zenith distance
pscale=np.zeros((nwave, nzd, nspots))  # plate scale
dfdx=np.zeros((nwave, nzd, nspots))    # differential field distorsion
dfdy=np.zeros((nwave, nzd, nspots))

# Read data from ZEMAX spot centroid files
for i in range(nwave):
    for j in range(nzd):
        infile=dir+'run'+auxrun+'_'+zdarr[j]+'_'+wavearr[i]+'.txt'
        tab=ascii.read(infile, encoding='utf-16')
        fld_x[i,j,:]=tab['fld_x'].data
        fld_y[i,j,:]=tab['fld_y'].data
        wave[i,j,:]=tab['lambda'].data
        chf_x[i,j,:]=tab['chf_x'].data
        chf_y[i,j,:]=tab['chf_y'].data
        cen_x[i,j,:]=tab['cen_x'].data
        cen_y[i,j,:]=tab['cen_y'].data
        rms[i,j,:]=tab['rms'].data
        n_rays[i,j,:]=tab['n_rays'].data
        zd[i,j,:]=float(zdarr[j])

# Calculate plate scale at field center and global [arcsec/mm] using w=500 nm and zd=0
fr=np.sqrt(fld_x[2,0,:]**2+fld_y[2,0,:]**2).flatten()
rcen=np.sqrt(cen_x[2,0,:]**2+cen_y[2,0,:]**2).flatten()
psarr=np.zeros(len(fr)) # array of platescales out to radius fr
for i in range(len(psarr)):
    psarr[i]=3600*fr[i]/rcen[i]
psarrsort=psarr[fr.argsort()]
pscenter=psarrsort[1]
psglobal=psarrsort[-1]

# Calculate plate scale at each spot sample
for i in range(nwave):
    for j in range(nzd):
        for k in range(nspots):
            fdist=np.sqrt((fld_x[i,j,k]-fld_x[i,j,:])**2+(fld_y[i,j,k]-fld_y[i,j,:])**2)
            cdist=np.sqrt((cen_x[i,j,k]-cen_x[i,j,:])**2+(cen_y[i,j,k]-cen_y[i,j,:])**2)
            auxratio=fdist*3600/cdist
            auxratiosort=auxratio[fdist.argsort()]
            pscale[i,j,k]=np.mean(auxratiosort[1:9])

print("Plate Scale:")
print("ZD, Wave, Platescale (Global/Center/Average) = ", wavearr[2], zdarr[0], psglobal, pscenter, np.mean(pscale[2, 0, :]))

# Calculate diferential field distorsions with respect to zenith at each wavelength
# which platescale shall I use? the one at zenith or the one at the ZD of interest??
# using the latter

for i in range(nwave):
    for j in range(nzd):
        for k in range(nspots):
            dfdx[i,j,k]=(cen_x[i,j,k]-cen_x[i,0,k])*pscale[i,j,k]/3600.
            dfdy[i,j,k]=(cen_y[i,j,k]-cen_y[i,0,k])*pscale[i,j,k]/3600.

# Evaluate scipy interpolation object and pickle it
intdfdx = LinearNDInterpolator(list(zip(wave.flatten(), zd.flatten(), fld_x.flatten(), fld_y.flatten())), dfdx.flatten())
intdfdy = LinearNDInterpolator(list(zip(wave.flatten(), zd.flatten(), fld_x.flatten(), fld_y.flatten())), dfdy.flatten())

with open("./"+auxrun+"_dfd.pkl","wb") as f:
    pickle.dump([intdfdx, intdfdy], f)

# Make platescale map plots
def mkpscaleplot(zd1, w1):
    selzd1=(zdarr.astype(float)==zd1)
    selw1=(wavearr.astype(float)==w1)
    x1=cen_x[selw1,selzd1,:].flatten()
    y1=cen_y[selw1,selzd1,:].flatten()
    pscale1=pscale[selw1,selzd1,:].flatten()

    fig, ax = plt.subplots(figsize=(10,10))
    ax.tricontour(x1, y1, pscale1, linewidths=0.5, colors='k')
    cntr2 = ax.tricontourf(x1, y1, pscale1, levels=14, cmap="RdBu_r")
    fig.colorbar(cntr2, ax=ax)
    ax.set_xlabel('X - [mm]', fontsize=20)
    ax.set_ylabel('Y - [mm]', fontsize=20)
    ax.set_title('Platescale ["/mm]', fontsize=25)
    ax.set_aspect('equal')
    plt.savefig('./plots/run'+auxrun+'_pscale_zd'+zdarr[selzd1][0]+'_w'+wavearr[selw1][0]+'.png')

mkpscaleplot(0,365)
mkpscaleplot(30,365)
mkpscaleplot(50,365)
mkpscaleplot(0,365)
mkpscaleplot(0,560)
mkpscaleplot(0,950)


# Make 2D differential field distorsion plots

def mkdistplot(zd1, w1, zd2, w2):

    selzd1=(zdarr.astype(float)==zd1)
    selw1=(wavearr.astype(float)==w1)
    x1=cen_x[selw1,selzd1,:].flatten()
    y1=cen_y[selw1,selzd1,:].flatten()

    selzd2=(zdarr.astype(float)==zd2)
    selw2=(wavearr.astype(float)==w2)
    x2=cen_x[selw2,selzd2,:].flatten()
    y2=cen_y[selw2,selzd2,:].flatten()

    title='Run '+auxrun+' ; zd='+zdarr[selzd1][0]+' vs zd='+zdarr[selzd2][0]+' ; '+wavearr[selw1][0]+'nm vs '+wavearr[selw2][0]+'nm'
    output='run'+auxrun+'_zd'+zdarr[selzd1][0]+'_zd'+zdarr[selzd2][0]+'_w'+wavearr[selw1][0]+'_w'+wavearr[selw2][0]

    mkresplot(x1,y1,x2,y2, title=title, output=output, factor=100, scale=5, platescale=pscenter)

# Residuals from zd=0 to zd=30/50 at 0.365 microns
mkdistplot(0, 365, 30, 365)
mkdistplot(0, 365, 50, 365)

# Residuals from zd=0 to zd=30/50 at 0.560 microns
mkdistplot(0, 560, 30, 560)
mkdistplot(0, 560, 50, 560)

# Residuals from zd=0 to zd=30/50 at 0.95 microns
mkdistplot(0, 950, 30, 950)
mkdistplot(0, 950, 50, 950)


# Make 1D differential dield distorsion  plots from zd=0 to zd=30/50 at 0.5 microns and fld_x = 0.00/0.75

def mk1ddistplot(zd1, zd2, w1, fld):
    fldsel=fld_x[2,0,:]==fld
    selw=(wavearr.astype(float)==w1)
    selzd1=(zdarr.astype(float)==zd1)
    selzd2=(zdarr.astype(float)==zd2)
    x1=cen_x[selw,selzd1,fldsel]
    y1=cen_y[selw,selzd1,fldsel]
    x2=cen_x[selw,selzd2,fldsel]
    y2=cen_y[selw,selzd2,fldsel]
    title='Run '+auxrun+' ; zd='+zdarr[selzd1][0]+' vs zd='+zdarr[selzd2][0]+' ; '+wavearr[selw][0]+'nm ; fld_x='+str(fld)
    output='run'+auxrun+'_zd'+zdarr[selzd1][0]+'_zd'+zdarr[selzd2][0]+'_w'+wavearr[selw][0]+'_fld'+str(fld)+'_1d'
    mkresplot1d(x1, y1, x2, y2, platescale=pscenter, title=title, output=output)

mk1ddistplot(0, 50, 560, 0)
mk1ddistplot(0, 50, 560, 0.75)
mk1ddistplot(0, 30, 560, 0)
mk1ddistplot(0, 30, 560, 0.75)

