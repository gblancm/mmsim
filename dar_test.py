import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.coordinates import SkyOffsetFrame, AltAz, SkyCoord
from mmsim import daltatm, daltatmslalib



# Read no-atm-refr file
tab=ascii.read('./zemax/20210425/parax-no-atm-refr_00_0500.txt')
fx=tab['fld_x'].data
fy=tab['fld_y'].data*(-1)
xcen=tab['cen_x'].data
ycen=tab['cen_y'].data*(-1)

xcennd=np.tan(fx*np.pi/180.)*23323.8  
ycennd=np.tan(fy*np.pi/180.)*23323.8  

# Calculate plate scale at field center and global [arcsec/mm]
fr=np.sqrt(fx**2+fy**2)
rcen=np.sqrt(xcen**2+ycen**2)
psarr=np.zeros(len(fr)) # array of platescales out to radius fr

for i in range(len(psarr)):
    psarr[i]=3600*fr[i]/rcen[i]
psarrsort=psarr[fr.argsort()]
pscenter=psarrsort[1]
psglobal=psarrsort[-1]

# Plot platescale maps

fig, ax = plt.subplots(figsize=(10,10))
ax.tricontour(xcen[np.isfinite(psarr)], ycen[np.isfinite(psarr)], psarr[np.isfinite(psarr)], linewidths=0.5, colors='k')
cntr2 = ax.tricontourf(xcen[np.isfinite(psarr)], ycen[np.isfinite(psarr)], psarr[np.isfinite(psarr)], levels=14, cmap="RdBu_r")
fig.colorbar(cntr2, ax=ax)
ax.set_xlabel('X - [mm]', fontsize=20)
ax.set_ylabel('Y - [mm]', fontsize=20)
ax.set_title('Platescale = field_r / spot_r ["/mm]', fontsize=25)
ax.set_aspect('equal')
plt.savefig('platescale.png')


# Make residual plots between two sets of spot plane coordinates

def mkresplot(x1, y1, x2, y2, factor=200, scale=1, platescale=pscenter, output='junk', title='Spot Plane'):

    fig, ax = plt.subplots(figsize=(30,30), dpi=200)
    ax.plot(x1, y1, 'o', markersize=7)
    for i in range(len(x1)):
        if (x1[i]!=x2[i])+(y1[i]!=y2[i]):
            ax.arrow(x1[i], y1[i], (x2[i]-x1[i])*factor, (y2[i]-y1[i])*factor, width=2, head_width=0, length_includes_head=True, head_starts_at_zero=True, color='red')
    auxdist=np.sqrt((x2-x1)**2+(y2-y1)**2)*platescale
    print("Min, Max, RMS = ", np.min(auxdist), np.max(auxdist), np.std(auxdist))
    ax.plot(np.array([np.min(x1), np.min(x1)+factor*scale/platescale]), np.array([np.max(y1),np.max(y1)]), color='black', linewidth=3)
    ax.text(np.min(x1)+factor*scale/platescale/4, np.max(y1)*1.05, "{:.0f}".format(scale)+'\"', fontsize=40)
    ax.set_xlabel('X - [mm]', fontsize=40)
    ax.set_ylabel('Y - [mm]', fontsize=40)
    ax.set_title(title, fontsize=50)
    ax.set_aspect('equal')
    ax.grid(color='grey',linestyle='--')
    ax.tick_params(axis='x', labelsize=40)
    ax.tick_params(axis='y', labelsize=40)
    plt.savefig(output+'.png')


# Plot platescale residuals
mkresplot(fx*3600/pscenter, 1*fy*3600/pscenter, xcen, ycen, title='fx;fy*(pscenter) vs xcen;ycen', output='junk', factor=200)
mkresplot(xcen, ycen, xcennd, ycennd, title='xcen;ycen vs xcennd;ycennd', output='cen_nd', factor=200)



# Read atmos file zd=0
tab=ascii.read('./zemax/20210425/parax-atmospheric_00_0500.txt')
fxa0=tab['fld_x'].data
fya0=tab['fld_y'].data*(-1)
xcena0=tab['cen_x'].data
ycena0=tab['cen_y'].data*(-1)

# Read atmos file zd=30
tab=ascii.read('./zemax/20210425/parax-atmospheric_30_0500.txt')
fxa30=tab['fld_x'].data
fya30=tab['fld_y'].data*(-1)
xcena30=tab['cen_x'].data
ycena30=tab['cen_y'].data*(-1)

# Read atmos file zd=50
tab=ascii.read('./zemax/20210425/parax-atmospheric_50_0500.txt')
fxa50=tab['fld_x'].data
fya50=tab['fld_y'].data*(-1)
xcena50=tab['cen_x'].data
ycena50=tab['cen_y'].data*(-1)


# Read wedge file zd=0
tab=ascii.read('./zemax/20210425/parax-vacuum-wdge_00_0500.txt')
fxw0=tab['fld_x'].data
fyw0=tab['fld_y'].data*(-1)
xcenw0=tab['cen_x'].data
ycenw0=tab['cen_y'].data*(-1)

# Read wedge file zd=30
tab=ascii.read('./zemax/20210425/parax-vacuum-wdge_30_0500.txt')
fxw30=tab['fld_x'].data
fyw30=tab['fld_y'].data*(-1)
xcenw30=tab['cen_x'].data
ycenw30=tab['cen_y'].data*(-1)

# Read wedge file zd=50
tab=ascii.read('./zemax/20210425/parax-vacuum-wdge_50_0500.txt')
fxw50=tab['fld_x'].data
fyw50=tab['fld_y'].data*(-1)
xcenw50=tab['cen_x'].data
ycenw50=tab['cen_y'].data*(-1)



mkresplot(xcenw0, ycenw0, xcenw30, ycenw30, title='Wedge zd=0 vs Wedge zd=30', output='w0_w30', factor=200)
mkresplot(xcenw0, ycenw0, xcenw50, ycenw50, title='Wedge zd=0 vs Wedge zd=50', output='w0_w50', factor=200)
mkresplot(xcena0, ycena0, xcena30, ycena30, title='Atmos zd=0 vs Atmos zd=30', output='a0_a30', factor=10000)
mkresplot(xcena0, ycena0, xcena50, ycena50, title='Atmos zd=0 vs Atmos zd=50', output='a0_a50', factor=10000)




# Compare wedge and atmosphere models to no atmosphere model

mkresplot(xcen, ycen, xcenw0, ycenw0, title='No Atm zd=0 vs Wedge zd=0', output='na0_w0', factor=200)
mkresplot(xcen, ycen, xcenw30, ycenw30, title='No Atm zd=0 vs Wedge zd=30', output='na0_w30', factor=200)
mkresplot(xcen, ycen, xcenw50, ycenw50, title='No Atm zd=0 vs Wedge zd=50', output='na0_w50', factor=200)

mkresplot(xcen, ycen, xcena0, ycena0, title='No Atm zd=0 vs Atmos zd=0', output='na0_a0', factor=200)
mkresplot(xcen, ycen, xcena30, ycena30, title='No Atm zd=0 vs Atmos zd=30', output='na0_a30', factor=200)
mkresplot(xcen, ycen, xcena50, ycena50, title='No Atm zd=0 vs Atmos zd=50', output='na0_a50', factor=200)


#Compare wedge and atmosphere models to atmospheric refraction model from AstroAtmosphere

# Function that applies relative atmospheric refraction to field angles (fx anf fy must be in degrees and contain (0,0))
def atmdistcoord(fx, fy, zd, wave):
    sel=(fx==0)*(fy==0) # origin at which pointing is done
    # Convert field angle to altaz coordinates and calculate zenith distance of each spot
    center = AltAz(180*u.deg, (90-zd)*u.deg, pressure=0)
    telframe = SkyOffsetFrame(origin=center)
    telcoord=SkyCoord(lat=fy*u.degree, lon=fx*u.degree, frame=telframe) 
    altazcoord=telcoord.transform_to(AltAz(pressure=0))
    zdarr=90-altazcoord.alt.value
    # Calculate and apply atmospheric refraction to altaz coordinates
#   dalt=daltatm(zdarr, wave) # Ciddor refraction model
    dalt=daltatmslalib(zdarr, wave) # SLALIB refraction model
    altazcoorddar=SkyCoord(alt=(altazcoord.alt.value+dalt/3600)*u.degree, az=altazcoord.az.value*u.degree, frame='altaz', pressure=0)
    # Convert abck to field angle in telescope frame
    telcoorddar=altazcoorddar.transform_to(telframe)
    # Return refracted field angles assuming pointing is done at origin (0,0), also return zenith distance [deg] and refraction [arcsec]
    return telcoorddar.lon.value-telcoorddar.lon.value[sel], telcoorddar.lat.value-telcoorddar.lat.value[sel], zdarr, dalt

fxdar0, fydar0, zdarr0, dalt0 = atmdistcoord(fx, fy, 0, 5000)
fxdar1, fydar1, zdarr1, dalt1 = atmdistcoord(fx, fy, 1, 5000)
fxdar5, fydar5, zdarr5, dalt5 = atmdistcoord(fx, fy, 5, 5000)
fxdar30, fydar30, zdarr30, dalt30 = atmdistcoord(fx, fy, 30, 5000)
fxdar50, fydar50, zdarr50, dalt50 = atmdistcoord(fx, fy, 50, 5000)

xcendar0=np.tan(fxdar0*np.pi/180.)*23323.8  
ycendar0=np.tan(fydar0*np.pi/180.)*23323.8  
xcendar1=np.tan(fxdar1*np.pi/180.)*23323.8  
ycendar1=np.tan(fydar1*np.pi/180.)*23323.8  
xcendar5=np.tan(fxdar5*np.pi/180.)*23323.8  
ycendar5=np.tan(fydar5*np.pi/180.)*23323.8  
xcendar30=np.tan(fxdar30*np.pi/180.)*23323.8  
ycendar30=np.tan(fydar30*np.pi/180.)*23323.8  
xcendar50=np.tan(fxdar50*np.pi/180.)*23323.8  
ycendar50=np.tan(fydar50*np.pi/180.)*23323.8  

mkresplot(xcen, ycen, xcendar0, ycendar0, title='No Atm zd=0 vs DAR zd=0', output='na0_dar0', factor=200)
mkresplot(xcen, ycen, xcendar1, ycendar1, title='No Atm zd=0 vs DAR zd=1', output='na0_dar1', factor=200)
mkresplot(xcen, ycen, xcendar5, ycendar5, title='No Atm zd=0 vs DAR zd=5', output='na0_dar5', factor=200)
mkresplot(xcen, ycen, xcendar30, ycendar30, title='No Atm zd=0 vs DAR zd=30', output='na0_dar30', factor=200)
mkresplot(xcen, ycen, xcendar50, ycendar50, title='No Atm zd=0 vs DAR zd=50', output='na0_dar50', factor=200)

# Compare wedge to DAR model
mkresplot(xcenw0, ycenw0, xcendar0, ycendar0, title='Wedge zd=0 vs DAR zd=0', output='w0_dar0', factor=6000)
mkresplot(xcenw30, ycenw30, xcendar30, ycendar30, title='Wedge zd=30 vs DAR zd=30', output='w30_dar30', factor=6000)
mkresplot(xcenw50, ycenw50, xcendar50, ycendar50, title='Wedge zd=50 vs DAR zd=50', output='w50_dar50', factor=6000)






