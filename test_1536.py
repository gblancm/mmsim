# This script calculates field distorsions for MegaMapper optical design 1536 
# assuming different observations (dec, HA_start, texp) 

import numpy as np
import astropy.io.ascii as ascii
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mmsim import mkresplot, mkresplot1d


run="1536revA"
dir='./zemax/20210817/run1536revA'
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


# Get number of spots
infile=dir+'run'+run+'_'+zdarr[0]+'_'+wavearr[0]+'.txt'
tab0=ascii.read(infile, encoding='utf-16')
nspots=len(tab0)


# Crate 3D arrays for all data
fld_x=np.zeros((nwave, nzd, nspots))
fld_y=np.zeros((nwave, nzd, nspots))
wave=np.zeros((nwave, nzd, nspots))
chf_x=np.zeros((nwave, nzd, nspots))
chf_y=np.zeros((nwave, nzd, nspots))
cen_x=np.zeros((nwave, nzd, nspots))
cen_y=np.zeros((nwave, nzd, nspots))
rms=np.zeros((nwave, nzd, nspots))
n_rays=np.zeros((nwave, nzd, nspots))
zd=np.zeros((nwave, nzd, nspots))

for i in range(nwave):
    for j in range(nzd):
        infile=dir+'run'+run+'_'+zdarr[j]+'_'+wavearr[i]+'.txt'
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

print("Platescale (Global/Center):", psglobal, pscenter)






# Residuals from zd=0 to zd=30/50 at 0.36 microns
mkresplot(cen_x[0,0,:].flatten(), cen_y[0,0,:].flatten(), cen_x[0,1,:].flatten(), cen_y[0,1,:].flatten(), title='Run '+run+' ; zd=0 vs zd=30 @ 360 nm', output='run'+run+'_zd0_zd30_0360', factor=100, scale=5, platescale=pscenter)
mkresplot(cen_x[0,0,:].flatten(), cen_y[0,0,:].flatten(), cen_x[0,2,:].flatten(), cen_y[0,2,:].flatten(), title='Run '+run+' ; zd=0 vs zd=50 @ 360 nm', output='run'+run+'_zd0_zd50_0360', factor=100, scale=5, platescale=pscenter)

# Residuals from zd=0 to zd=30/50 at 0.5 microns
mkresplot(cen_x[2,0,:].flatten(), cen_y[2,0,:].flatten(), cen_x[2,1,:].flatten(), cen_y[2,1,:].flatten(), title='Run '+run+' ; zd=0 vs zd=30 @ 500 nm', output='run'+run+'_zd0_zd30_0500', factor=100, scale=5, platescale=pscenter)
mkresplot(cen_x[2,0,:].flatten(), cen_y[2,0,:].flatten(), cen_x[2,2,:].flatten(), cen_y[2,2,:].flatten(), title='Run '+run+' ; zd=0 vs zd=50 @ 500 nm', output='run'+run+'_zd0_zd50_0500', factor=100, scale=5, platescale=pscenter)

# Residuals from zd=0 to zd=30/50 at 0.95 microns
mkresplot(cen_x[7,0,:].flatten(), cen_y[7,0,:].flatten(), cen_x[7,1,:].flatten(), cen_y[7,1,:].flatten(), title='Run '+run+' ; zd=0 vs zd=30 @ 950 nm', output='run'+run+'_zd0_zd30_0950', factor=100, scale=5, platescale=pscenter)
mkresplot(cen_x[7,0,:].flatten(), cen_y[7,0,:].flatten(), cen_x[7,2,:].flatten(), cen_y[7,2,:].flatten(), title='Run '+run+' ; zd=0 vs zd=50 @ 950 nm', output='run'+run+'_zd0_zd50_0950', factor=100, scale=5, platescale=pscenter)


# 1D residual plots from zd=0 to zd=30/50 at 0.5 microns and fld_x = 0.00/0.75

fldsel=fld_x[2,0,:]==0.0
x1=cen_x[2,0,fldsel]
y1=cen_y[2,0,fldsel]
x2=cen_x[2,2,fldsel]
y2=cen_y[2,2,fldsel]
mkresplot1d(x1, y1, x2, y2, platescale=pscenter, title='Run '+run+' ; zd=0 vs zd=50 @ 500 nm - fld_x=0.00', output='run'+run+'_zd0_zd50_0500_000_1d')

fldsel=fld_x[2,0,:]==0.75
x1=cen_x[2,0,fldsel]
y1=cen_y[2,0,fldsel]
x2=cen_x[2,2,fldsel]
y2=cen_y[2,2,fldsel]
mkresplot1d(x1, y1, x2, y2, platescale=pscenter, title='Run '+run+' ; zd=0 vs zd=50 @ 500 nm - fld_x=0.75', output='run'+run+'_zd0_zd50_0500_075_1d')

fldsel=fld_x[2,0,:]==0.0
x1=cen_x[2,0,fldsel]
y1=cen_y[2,0,fldsel]
x2=cen_x[2,1,fldsel]
y2=cen_y[2,1,fldsel]
mkresplot1d(x1, y1, x2, y2, platescale=pscenter, title='Run '+run+' ; zd=0 vs zd=30 @ 500 nm - fld_x=0.00', output='run'+run+'_zd0_zd30_0500_000_1d')

fldsel=fld_x[2,0,:]==0.75
x1=cen_x[2,0,fldsel]
y1=cen_y[2,0,fldsel]
x2=cen_x[2,1,fldsel]
y2=cen_y[2,1,fldsel]
mkresplot1d(x1, y1, x2, y2, platescale=pscenter, title='Run '+run+' ; zd=0 vs zd=30 @ 500 nm - fld_x=0.75', output='run'+run+'_zd0_zd30_0500_075_1d')

# Maximum wavelength effect on distortion

x1=cen_x[0,2,:]
y1=cen_y[0,2,:]
x2=cen_x[-1,2,:]
y2=cen_y[-1,2,:]
mkresplot(x1,y1,x2,y2, title='Run '+run+' ; 360 nm vs 950 nm @ zd=50', output='run'+run+'_w360_w950_zd50', factor=100, scale=5, platescale=pscenter)


# Make plots

scale=5
factor=100
platescale=pscenter

fig, ax = plt.subplots(figsize=(30,30), dpi=200)
colorarr=cm.rainbow(np.linspace(0, 1, nwave))
markerarr=['P', 'o', 's']
x0=cen_x[2,0,:].flatten()
y0=cen_y[2,0,:].flatten()
for i in range(nwave):
    for j in range(nzd):
        ax.plot(x0+(cen_x[i,j,:].flatten()-x0)*factor, y0+(cen_y[i,j,:].flatten()-y0)*factor, linestyle='', marker=markerarr[j], color=colorarr[i], alpha=0.5)  
for k in range(nspots):
    circle=plt.Circle((x0[k], y0[k]), 0.5/platescale*factor, fill=False)
    ax.add_artist(circle)  
ax.plot(np.array([np.min(cen_x), np.min(cen_x)+factor*scale/platescale]), np.array([np.max(cen_y),np.max(cen_y)]), color='black', linewidth=3)
ax.text(np.min(cen_x)+factor*scale/platescale/4, np.max(cen_y)*1.05, "{:.0f}".format(scale)+'\"', fontsize=40)
ax.set_xlabel('X - [mm]', fontsize=40)
ax.set_ylabel('Y - [mm]', fontsize=40)
ax.set_title('Run '+run, fontsize=50)
ax.set_aspect('equal')
ax.grid(color='grey',linestyle='--')
ax.tick_params(axis='x', labelsize=40)
ax.tick_params(axis='y', labelsize=40)
plt.savefig('junk.png')

