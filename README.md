# mmsim

MegaMapper field distorion and atmospheric refarction/dispersion simulation code:

See test1.py, test2.py, and most of all test3.py to learn how to use it.

Basic usage:

% import mmsim
% mysim = mmsim.distsim()
% mysim.run()

By default (i.e. if no attributes of the "distsim" object are modified) the simulation assumes the following parameters:

lat=-29.028535, lon=-70.685801, height=2380, temp=285.85, pres=76400, rh=0.27, xc=450

# Site parameters for LCO
mysim.site.lat= # latitude [deg]
mysim.site.lon= # longitude [deg]
mysim.site.height=2635 # height above sea level [m]
mysim.site.temp=285.95 # temperature [K]
mysim.site.pres=74300 # pressue [Pa]
mysim.site.rh=0.15 # relative humidity
mysim.site.xc=450 # CO2 concentration [ppm]

# Set field of view parameters
mysim.fov.rfov=8 # FOV radius [arcmin]
mysim.fov.posang=0 # position angle of regular grid of FOV sampling points

# Set observation parameters
mysim.obs.hastart=-30 # starting HA [deg]
mysim.obs.decstart=25 # field center declination [deg]
mysim.obs.texp=14400 # exposure time [s]

# Set simulation parameters
mysim.wmin=3800 # bluest sampled wavelength [A]
mysim.wmax=6500 # reddest sampled wavelength [A]
mysim.Nw=5 # number of wavelength samplings within range; guiding is done at the central wavelngth so needs to be odd
mysim.Nt=31 # number of time samplings within exposure time ; aperture is drawn at center of exposure, so ideally odd
mysim.Nf=3 # number of sampling grid steps per field radius
mysim.ADC=True # apply a perfect ADC correction for the center of the field
mysim.output='distsim_cuby200_1a' # name of the output filenames

