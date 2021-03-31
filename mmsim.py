import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import AstroAtmosphere


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
            number of samples in field angle
        ADC: boolean
            False by default, if True it applies an instrumental distorsion on top of the atmosphere
        output: str
            Output filename
    """
    def __init__(self, wmin=3600, wmax=10000, Nw=5, Nt=101, Nf=3, ADC=False, output='distsim'):
        self.site=site()
        self.obs=obs()
        self.fov=fov()
        self.wmin=wmin
        self.wmax=wmax
        self.Nw=Nw
        self.Nt=Nt
        self.Nf=Nf
        self.ADC=ADC
        self.output=output


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

        # Compute RA,DEC of sampling points in the FOV at start of the exposure
            # Compute regular cartesian grid of sampling points on focal plane
        nx=2*self.Nf+1
        arr=np.linspace(-1*self.fov.rfov, self.fov.rfov, nx)
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
            # Calculate observing time so hour angle at start is hastart
        lco=EarthLocation(lat=self.site.lat*u.deg, lon=self.site.lon*u.deg, height=self.site.height*u.m)
        taux=Time('1993-02-17 00:00:00', scale='utc', location=lco) # my birthday!!
        dtaux=TimeDelta((lststart/15-taux.sidereal_time('apparent').value)*3600.*23.9344696/24, format='sec')
        tobs0=taux+dtaux
        
    
        print("Declination: ", self.obs.decstart)
        print("Start Hour Angle: ", tobs0.sidereal_time('apparent').value-ra0/15)
        print("Exposure Time: ", self.obs.texp)

        # Define index for guiding wavelength, centerfield sample, and center of exposure time
        wind=int(self.Nw/2.)
        print("Wavelength Samples: ", self.warr)
        print("Guiding Wavelength is: ", self.warr[wind])
        xind=self.Nf
        tind=int(self.Nt/2.)

        
            # Calculate altaz coordinates without distorsion (pressure=0)
        self.altaz=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                dt=TimeDelta(self.tarr[i], format='sec')
                tobs=tobs0+dt
                self.altaz[i,j]=self.radecstart.transform_to(AltAz(obstime=tobs,location=lco, pressure=0))

        #colorarr=cm.rainbow(np.linspace(0, 1, self.Nt))
        #for i in range(self.Nt):
        #    plt.plot(self.altaz[i,0].az, self.altaz[i,0].alt, '.', color=colorarr[i])
        #plt.show()

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
       
        #colorarr=cm.rainbow(np.linspace(0, 1, self.Nw))
        #for i in range(self.Nt):
        #    for j in range(self.Nw):
        #        plt.plot(self.altaz[i,j].alt.value, self.atmref[i,j], '.', color=colorarr[j])
        #plt.show()

        # Compute instrumental field distorsions for each sampling point at every time and wavelength
        # In the future will read something based on ZEMAX output, for now just applying a perfect dispersion
        # correction for the center of the field.

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
                if self.ADC:
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

        # Compute RA,DEC corresponding to distorted ALT/AZ coordinates at every time and wavelength
        self.radecdist=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                self.radecdist[i,j]=self.altazdist[i,j].icrs
        
        # Compute dRA,dDEC offsets [in arcsec] with respect to start of the simulation as a function of time and wavelength
        self.ddec=np.ndarray((self.Nt, self.Nw), dtype=object)
        self.dra=np.ndarray((self.Nt, self.Nw), dtype=object)
        for i in range(self.Nt):
            for j in range(self.Nw):
                self.ddec[i,j]=(self.radecdist[i,j].dec.value-dec0)*3600
                self.dra[i,j]=(self.radecdist[i,j].ra.value-ra0)*np.cos(np.deg2rad(self.radecdist[i,j].dec.value))*3600
        
        # Plot results vs time and wavelength

        #colorarr=cm.rainbow(np.linspace(0, 1, self.Nw))
        #for i in range(self.Nt):
        #    for j in range(self.Nw):
        #        plt.plot(self.altaz[i,j].alt.value, self.atmdis[i,j], color=colorarr[j])
        #        plt.plot(self.altaz[i,j].alt.value, self.atmref[i,j], color=colorarr[j])
        #plt.show()


        #fig, ax = plt.subplots()
        #ax.plot(xrot*60, yrot*60, 'o', color='black')
        #circle=plt.Circle((0,0), 1.5*3600., fill=False, color='red')
        #ax.add_artist(circle)
        #tind=int(self.Nt/2.)
        #colorarr=cm.rainbow(np.linspace(0, 1, self.Nw))
        #alphaarr=np.linspace(0.1, 1.0, self.Nt)
        #for j in range(self.Nw):
        #    ax.plot(-1*self.dra[tind,j], self.ddec[tind,j], 'o', color=colorarr[j], alpha=alphaarr[tind])
        #for i in range(self.Nt):
        #    for j in range(self.Nw):
        #        ax.plot(-1*self.dra[i,j], self.ddec[i,j], '.', color=colorarr[j], alpha=alphaarr[i])
        #circle=plt.Circle((0,0), 1.5*3600., fill=False, color='red')
        #ax.add_artist(circle)
        #for i in range(2*self.Nf+1):
        #    for j in range(2*self.Nf+1):
        #        circle=plt.Circle((xrot[i,j]*60, yrot[i,j]*60), 0.5, fill=False)
        #        ax.add_artist(circle)
        #        #print(xrot[i,j]*60, yrot[i,j]*60)
        #plt.show()


        # Plot results vs time and wavelength growing offsets by "factor"
        factor=700*(self.fov.rfov/90)
        fig, ax = plt.subplots(figsize=(40,40), dpi=200)
        ax.plot(xrot*60, yrot*60, 'o', color='black')
        circle=plt.Circle((0,0), self.fov.rfov*60., fill=False, color='red')
        ax.add_artist(circle)
        colorarr=cm.rainbow(np.linspace(0, 1, self.Nw))
        alphaarr=np.linspace(0.1, 1.0, self.Nt)
        for j in range(self.Nw):
            ax.plot(xrot*60+(-1*self.dra[tind,j]-xrot*60)*factor, yrot*60+(self.ddec[tind,j]-yrot*60)*factor, 'o', color=colorarr[j], alpha=alphaarr[tind])
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
        plt.savefig(self.output+'.png')
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
            field of view radius [arcmin]
        posang: flt
            position angle of uniform RA, DEC grid
    """
    def __init__(self, rfov=90, posang=0):
        self.rfov=rfov
        self.posang=posang


