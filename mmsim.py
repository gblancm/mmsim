import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import astropy.units as u
from astropy.time import Time, TimeDelta
import AstroAtmosphere
from AstroAtmosphere.refractivityModels import HohenkerkAndSinclair

import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, SkyOffsetFrame
#from astropy.coordinates import frame_transform_graph
#from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
import pickle

class distsim:
    """
        Contains and computes a MegaMapper field distorsion simulation, including
        atmospheric field distorsion, differential refraction, and instrumental field distorisons.
        Attributes:
        -----------
        site: Class site
            site parameters
        obs: Class obs
            observation parameters
        fov: Class fov
            field-of-view parameters
        wmin: flt
            bluest wavelength to simulate
        wmax: flt
            reddest wavelength to simulate
        Nw: int
            number of samples in wavelength
        Nt: int
            number of samples in time during exposure
        Nf: int
            number of RA/DEC samples from center to edge of telescope FoV
        simtype: str
            'atmos' = simulation of differential distortion from atmospheric refraction 
            'atmoscorr' = simulation of of differential distortion from astmospheric refraction effects with dispersion correction at field center
            'atmostel' = simulation of of differential distortion from astmospheric refraction and telescope optics
        output: str
            Output filename
    """
    def __init__(self, wmin=3650, wmax=9500, Nw=5, Nt=11, Nf=3, dfdfile='1536revA_dfd.pkl', simtype='atmos', output='distsim'):
        self.site=site()
        self.obs=obs()
        self.fov=fov()
        self.wmin=wmin
        self.wmax=wmax
        self.Nw=Nw
        self.Nt=Nt
        self.Nf=Nf
        self.simtype=simtype
        self.output=output
        self.dfdfile=dfdfile
        self.dfdrate=0


    def run(self):
        """
        Method that runs the simulation
        Input:
        ------
        self: Class distsim
            distrosion simulation object
        Output:
        -------
        warr: flt(Nw)
            samples in wavelength for simulation
        tarr: flt(Nt)
            samples in time for the simulation        
        """
        
        # Compute wavelength and time sampling arrays
        self.warr=np.linspace(self.wmin, self.wmax, self.Nw)
        self.tarr=np.linspace(0, self.obs.texp, self.Nt)

        # Define indices for guiding wavelength, centerfield sample, and center of exposure time
        wind=int(self.Nw/2.)
        print("Wavelength Samples: ", self.warr)
        print("Guiding Wavelength is: ", self.warr[wind])
        xind=self.Nf  # center of the field is at index (xind,xind)
        tind=int(self.Nt/2.)

        # Compute RA,DEC of sampling points in the FOV at start of the exposure
            # Compute regular cartesian grid of sampling points on focal plane
        nx=2*self.Nf+1
        arr=np.linspace(-1*self.fov.rfov, self.fov.rfov, nx)*60.*0.99
        grid=np.meshgrid(arr, arr)
        xarr=grid[0]
        yarr=grid[1]
            # Rotate regular grid on focal plane to position angle
        xrot=xarr*np.cos(np.deg2rad(self.fov.posang))-yarr*np.sin(np.deg2rad(self.fov.posang))
        yrot=xarr*np.sin(np.deg2rad(self.fov.posang))+yarr*np.cos(np.deg2rad(self.fov.posang))
            # Compute RA,DEC of rotated grid, assume RA=180 and compute LSTstart from hastart
        ra0=180
        dec0=self.obs.decstart
        lststart=self.obs.hastart-ra0
        decstart=dec0+yrot/60.
        rastart=ra0-xrot/60./np.cos(np.deg2rad(decstart))
        self.radecstart=SkyCoord(ra=rastart*u.degree, dec=decstart*u.degree, frame='icrs')

        # Compute ALT,AZ of sampling points at every time step without distorsions
            # Calculate an aribitraty observing time so the hour angle at start equals hastart
        lco=EarthLocation(lat=self.site.lat*u.deg, lon=self.site.lon*u.deg, height=self.site.height*u.m)
        taux=Time('1993-02-17 00:00:00', scale='utc', location=lco) # my 10th birthday!!
        dtaux=TimeDelta((lststart/15-taux.sidereal_time('apparent').value)*3600.*23.9344696/24, format='sec')
        tobs0=taux+dtaux
        print("Declination: ", self.obs.decstart)
        print("Start Hour Angle: ", tobs0.sidereal_time('apparent').value-ra0/15)
        print("Exposure Time: ", self.obs.texp)
            # Calculate altaz coordinates without distorsion (pressure=0)
        self.altaz=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                dt=TimeDelta(self.tarr[i], format='sec')
                tobs=tobs0+dt
                self.altaz[i,j]=self.radecstart.transform_to(AltAz(obstime=tobs,location=lco, pressure=0)) # pressure=0 implies no atmospheric distorsion correction is applied in coordinate transformation


        # Atmospheric refraction simulations using refraction models in AstroAtmosphere
        if (self.simtype=='atmos')+(self.simtype=='atmoscorr'):
            # Compute atmospheric refraction distorsion for each sampling point at every time and wavelength   
            # Initializing dispersion model
            at=AstroAtmosphere.Observatory()
            # Calculating indices of refraction for each wavelength
            narr=np.ndarray(self.Nw, dtype=float)
            for i in range(self.Nw):
                narr[i]=at.n_tph(l=self.warr[i]/1e4, T=self.site.temp, p=self.site.pres, RH=self.site.rh, xc=self.site.xc)
                # Density of the atmosphere (following CIPM-81/91 equations)
            rho = at.rho(p=self.site.pres, T=self.site.temp, RH=self.site.rh, xc=self.site.xc)
                # Initializing refraction model and setting the reduced height
            ref = AstroAtmosphere.refraction(lat=self.site.lat, h=self.site.height)
            ref.setReducedHeight(p=self.site.pres, rho=rho)
                # Calculating the atmospheric refraction
            self.atmref=np.ndarray((self.Nt, self.Nw), dtype=object)
            for i in range(self.Nt):
                for j in range(self.Nw):
                    self.atmref[i,j]=ref.cassini(narr[j], zenith=90-self.altaz[i,j].alt.value)*3600

            # Initializing dispersion object
            disp = AstroAtmosphere.dispersion(lat=self.site.lat, h=self.site.height)
            disp.setReducedHeight(p=self.site.pres, rho=rho)
            self.atmdis=np.ndarray((self.Nt, self.Nw), dtype=object)
            self.instdist=np.ndarray((self.Nt, self.Nw), dtype=object)
            self.totdist=np.ndarray((self.Nt, self.Nw), dtype=object)
            for i in range(self.Nt):
                for j in range(self.Nw):
                        # Atmospheric dispersion
                    self.atmdis[i,j]=disp.cassini(narr[0], narr[j], zenith=90-self.altaz[i,j].alt.value)*3600
                        # Instrumental distorsion
                    if self.simtype=='atmoscorr':
                        self.instdist[i,j]=(self.atmdis[i,j])[xind,xind] # SET INSTRUMENTAL DISTORSION TO PERFECT ADC CORRECTION AT CENTERFIELD
                    else:
                        self.instdist[i,j]=self.atmref[i,j]-self.atmref[i,j] # SET INSTRUMENTAL DISTORSION TO ZERO
                        # Total is atmosphere + instrumental
                    self.totdist[i,j]=self.atmref[i,j]+self.instdist[i,j]

            # Apply guiding correction at given field position and wavelength to sum of the atmospheric and instrumental distorsion
            self.gcordist=np.ndarray((self.Nt, self.Nw), dtype=object)
            for i in range(self.Nt):
                for j in range(self.Nw):
                    self.gcordist[i,j]=self.totdist[i,j]-(self.totdist[i,wind])[xind,xind]
                    #self.gcordist[i,j]=self.totdist[i,j]   # REMOVE, TEST WITH NO GUIDING

            # Apply guiding corrected total distorsions to ALT/AZ coordinates at every time and wavelength
            self.altazdist=self.altaz 
            for i in range(self.Nt):
                for j in range(self.Nw):
                    dt=TimeDelta(self.tarr[i], format='sec')
                    tobs=tobs0+dt
                    self.altazdist[i,j]=SkyCoord(alt=(self.altazdist[i,j].alt.value+self.gcordist[i,j]/3600)*u.degree, az=self.altazdist[i,j].az.value*u.degree, frame='altaz', obstime=tobs,location=lco, pressure=0)

            
        
        # Do simulation based on ZEMAX model of atmosphere+telescope
        elif self.simtype=='atmostel':

            # Placeholder for function that will read distorsion info from ZEMAX simulations (fx and fy can be arrays)
            intdfdx, intdfdy = pickle.load(open(self.dfdfile, 'rb'))
            #def teldist(zd, wave, fx, fy):
            #    return intdfdx(wave/1e4, zd, fx, fy), intdfdy(wave/1e4, zd, fx, fy)
                
            # Calculate telescope frame field angle coordinates, add distorsion, and transform back to AltAz
            self.telfield=np.ndarray((self.Nt, self.Nw), dtype=object)
            self.telfielddist=np.ndarray((self.Nt, self.Nw), dtype=object)
            self.altazdist=np.ndarray((self.Nt, self.Nw), dtype=object)
            for i in range(self.Nt):
                print(i, self.Nt)
                for j in range(self.Nw):
                    # transform from AltAz to field angle
                    center=self.altaz[i,j][xind,xind]
                    telframe=center.skyoffset_frame()
                    self.telfield[i,j]=self.altaz[i,j].transform_to(telframe)
                    # calculate distorsion offsets and apply them to field angle coordinates
                    wave=self.warr[j]
                    zd=90-self.altaz[i,j].alt.value[xind,xind]
                    fx=self.telfield[i,j].lon.value
                    fy=self.telfield[i,j].lat.value
                    dx=np.empty_like(fx)
                    dy=np.empty_like(fx)
                    #dx,dy=teldist(zd, wave, fx, fy) # evaluate atm+telescope distortion model
                    for k in range(2*self.Nf+1):
                        for l in range(2*self.Nf+1):
                            dx[k,l]=intdfdx(wave/1e4, zd, fx[k,l], fy[k,l])
                            dy[k,l]=intdfdy(wave/1e4, zd, fx[k,l], fy[k,l])
                    self.telfielddist[i,j]=SkyCoord(lon=(fx+dx)*u.deg, lat=(fy+dy)*u.deg, frame=telframe)
                    # transform back from field angle to AltAz 
                    self.altazdist[i,j]=self.altaz[i,j]
                    self.altazdist[i,j]=self.telfielddist[i,j].transform_to(AltAz(obstime=self.altaz[i,j].obstime,location=self.altaz[i,j].location, pressure=0))
        

        else:
            print("'simtype' must be 'atmos', 'atmoscorr' or 'atmostel'")

        
        # Compute RA,DEC corresponding to distorted ALT/AZ coordinates at every time and wavelength
        self.radecdist=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                auxcoord=self.altazdist[i,j].transform_to('icrs')
                self.radecdist[i,j]=self.altazdist[i,j].transform_to('icrs')

        
        # Compute dRA,dDEC offsets [in arcsec] with respect to start of the simulation as a function of time and wavelength
        self.ddec=np.ndarray((self.Nt, self.Nw), dtype=object)
        self.dra=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                self.ddec[i,j]=(self.radecdist[i,j].dec.value-dec0)*3600
                self.dra[i,j]=(self.radecdist[i,j].ra.value-ra0)*np.cos(np.deg2rad(self.radecdist[i,j].dec.value))*3600


        # Compute maximum differential field distortion change rate (arcsec/hr) over exposure
        dtot=np.sqrt((self.ddec[self.Nt-1,wind]-self.ddec[0,wind])**2+(self.dra[self.Nt-1,wind]-self.dra[0,wind])**2)
        rate=(dtot*3600.)/(self.obs.texp/3600.)
        self.dfdrate=np.max(rate[np.isfinite(rate)])

        # Plot results vs time and wavelength growing offsets by "factor"
        factor=700*(self.fov.rfov/1.5)
        fig, ax = plt.subplots(figsize=(40,40), dpi=200)
        # plot sky samples in RA,DEC
        ax.plot(xrot*60, yrot*60, 'o', color='black')
        ax.set_xlim(-1.1*self.fov.rfov*3600, 1.1*self.fov.rfov*3600)
        ax.set_ylim(-1.1*self.fov.rfov*3600, 1.1*self.fov.rfov*3600)
        # plot FOV
        circle=plt.Circle((0,0), self.fov.rfov*3600., fill=False, color='red')
        ax.add_artist(circle)
        # plot sky samples after field distorsion
        colorarr=cm.rainbow(np.linspace(0, 1, self.Nw))
        alphaarr=np.linspace(0.1, 1.0, self.Nt)
        # samples at mid-exposure
        for j in range(self.Nw):
            ax.plot(xrot*60+(-1*self.dra[tind,j]-xrot*60)*factor, yrot*60+(self.ddec[tind,j]-yrot*60)*factor, 'o', color=colorarr[j], alpha=alphaarr[tind])
        # samples at each time step
        for i in range(self.Nt):
            for j in range(self.Nw):
                ax.plot(xrot*60+(-1*self.dra[i,j]-xrot*60)*factor, yrot*60+(self.ddec[i,j]-yrot*60)*factor, '+', color=colorarr[j], alpha=alphaarr[i])
        
        for i in range(2*self.Nf+1):
            for j in range(2*self.Nf+1):
                circle=plt.Circle((xrot[i,j]*60, yrot[i,j]*60), 0.1*factor, fill=False)
                ax.add_artist(circle)
                circle=plt.Circle(((xrot*60+(-1*self.dra[tind,wind]-xrot*60)*factor)[i,j], (yrot*60+(self.ddec[tind,wind]-yrot*60)*factor)[i,j]), 0.5*factor, fill=False)
                ax.add_artist(circle)
                #print(xrot[i,j]*60, yrot[i,j]*60)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\Delta \alpha$ ["]', fontsize=50)
        ax.set_ylabel(r'$\Delta \delta$ ["]', fontsize=50)
        ax.tick_params(axis='x', labelsize=25)
        ax.tick_params(axis='y', labelsize=25)
        ax.set_title(r"$\delta = $"+"{:.6f}".format(self.obs.decstart)+r" deg ; HA($t_{0}$) = "+"{:.1f}".format(self.obs.hastart/15)+r" hr ; $t_{exp}$ = "+"{:.0f}".format(self.obs.texp)+" s ; Circle = 1\"", fontsize=50)
        plt.savefig('./plots/'+self.output+'.png')
        #plt.show()
    
class site:
    """
        Defines an observing site with its parameters
        Deafult = LCO
        Attributes:
        -----------
        lat: flt
            site latitude [deg]
        lon: flt
            site longitude [deg]
        height: flt
            site elevation above sea level [m]
        temp: flt
            site average night-time temperature [K]
        pres: flt
            site average atmospheric pressure [Pa]
        rh: flt
            site average relative humidity [fraction]
        xc: flt
            site averge CO2 density [ppm]
    """
    def __init__(self, lat=-29.028535, lon=-70.685801, height=2380, temp=285.85, pres=76400, rh=0.27, xc=450):
        self.lat=lat
        self.lon=lon
        self.height=height
        self.temp=temp
        self.pres=pres
        self.rh=rh
        self.xc=xc


class obs:
    """
        Defines an observation with its parameters
        Default: start and finish one hour before and after passing through zenit, and 0.65" seeing
                 sample 5 wavelengths between 0.36 and 1 micron, sample 10 time steps during exposure
        Attributes:
        -----------
        hastart: flt
            hour angle of FOV center at start of observation [deg]
        decstart: flt
            declination of FOV center at start of observation [deg]
        texp: flt
            exposure time [s]
        tarr: flt(Ntime)
            time array with 10 time samples during expsure for simulation
        warr: flt(Nwave)
            wavelength array with wavelength samples for simulation
        seeing: flt
            seeing FWHM [arcsec]
    """
    def __init__(self, hastart=-15, decstart=-29.028535, texp=7200, seeing=0.65):
        self.hastart=hastart
        self.decstart=decstart
        self.texp=texp
        self.seeing=seeing


class fov:
    """
        Defines the field of view (FOV) of the telescope and the sampling of objects within it
        Attributes
        ----------
        rfov: flt
            field of view radius [deg]
        posang: flt
            position angle of uniform RA, DEC grid
    """
    def __init__(self, rfov=1.5, posang=0):
        self.rfov=rfov
        self.posang=posang


# Function that returns offset in altitude [arcsec] caused by DAR at LCO for a given wavelength [Angstroms] and zenit distance [deg]
def daltatm(zd, wave):
    at=AstroAtmosphere.Observatory()
    n=at.n_tph(l=wave/1e4, T=273.15+12, p=101325*0.77, RH=0.27, xc=450) # LCO parameters
    rho=at.rho(p=101325*0.77, T=273.15+12, RH=0.27, xc=450)
    ref=AstroAtmosphere.refraction(lat=-29.028535, h=2380)
    ref.setReducedHeight(p=101325*0.77, rho=rho)
    return ref.cassini(n, zenith=zd)*3600

# Function that returns offset in altitude [arcsec] caused by DAR at LCO for a given wavelength [Angstroms] and zenit distance [deg] using SLALIB model
def daltatmslalib(zd, wave):
    daltarr=np.zeros(len(zd))
    for i in range(len(zd)):
        daltarr[i]=AstroAtmosphere.slalib_refraction(wave/10000., zd[i], conditions=None, T=273.15+12, p=101325*0.77, RH=0.27, xc=450, lat=-29.028535, h=2380)*3600.
    return daltarr

# Function that returns Ciddor index of refraction as a function of wavelength, pressure, and temperature
# Assumes RH=0.27 and xc=450
def natm(wave, t, p):
    at=AstroAtmosphere.Observatory()
    n=at.n_tph(l=wave/1e4, T=t, p=p, RH=0.27, xc=450) # LCO parameters
    return n

# Function that returns SLALIB (HohenkerkAndSinclair) index of refraction as a function of wavelength, pressure, and temperature
# Assumes RH=0.27 and xc=450

def natmslalib(wave, t, p):
    n=HohenkerkAndSinclair(l=wave/1e4, T=t, p=p, RH=0.27, lat=-29.028535, h=2380) # LCO parameters
    return n

def mkresplot(x1, y1, x2, y2, factor=200, scale=1, platescale=1, output='junk', title='Spot Plane'):

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
    plt.savefig('./plots/'+output+'.png')



def mkresplot1d(x1,y1,x2,y2,platescale=1,title='1D Residuals', output='junk'):
    dx=x2-x1
    dy=y2-y1
    fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(25, 25))
    ax1.plot(y1, dx, '-o')
    ax2.plot(y1, dy, '-o')
    ax3.plot(y1, np.sqrt(dx**2+dy**2), '-o')
    ax1.set_xlabel(r'Y - [mm]', fontsize=30)
    ax1.set_ylabel(r'$\Delta$ X - [mm]', fontsize=30)
    ax2.set_xlabel(r'Y - [mm]', fontsize=30)
    ax2.set_ylabel(r'$\Delta$ Y - [mm]', fontsize=30)
    ax3.set_xlabel(r'Y - [mm]', fontsize=30)
    ax3.set_ylabel(r'$\sqrt{\Delta X^2 + \Delta Y^2}$ - [mm]', fontsize=30)
    ax1.tick_params(axis='x', labelsize=25)
    ax1.tick_params(axis='y', labelsize=25)
    ax2.tick_params(axis='x', labelsize=25)
    ax2.tick_params(axis='y', labelsize=25)
    ax3.tick_params(axis='x', labelsize=25)
    ax3.tick_params(axis='y', labelsize=25)
    def mm2arcsec(x):
        return x * platescale
    def arcsec2mm(x):
        return x /platescale
    secax1 = ax1.secondary_yaxis('right', functions=(mm2arcsec, arcsec2mm))
    secax2 = ax2.secondary_yaxis('right', functions=(mm2arcsec, arcsec2mm))
    secax3 = ax3.secondary_yaxis('right', functions=(mm2arcsec, arcsec2mm))
    secax1.set_ylabel(r'$\Delta$ X - ["]', fontsize=30)
    secax2.set_ylabel(r'$\Delta$ Y - ["]', fontsize=30)
    secax3.set_ylabel(r'$\sqrt{\Delta X^2 + \Delta Y^2}$ - ["]', fontsize=30)
    secax1.tick_params(axis='y', labelsize=25)
    secax2.tick_params(axis='y', labelsize=25)
    secax3.tick_params(axis='y', labelsize=25)
    ax1.set_title(title, fontsize=50)
    plt.savefig('./plots/'+output+'.png')
