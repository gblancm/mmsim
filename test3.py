# This scripts attempts to reproduce the results for VIMOS in Paranal from Cuby et al. 2000
# "Handling atmospheric dispersion and differential refraction effects in large-field multi-objects spectroscopic observations"

import mmsim

mysim=mmsim.distsim()

# Set location and conditions of Paranal
mysim.site.lat=-24.6275 # latitude [deg]
mysim.site.lon=-70.4044 # longitude [deg]
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

# Run the simulation to reproduce Figure 1a in Cuby et al. 200
mysim.run()

# Now do Figure 1b
mysim.obs.decstart=-75
mysim.output='distsim_cuby200_1b' 
mysim.run()

# Now do Figure 3a
mysim.ADC=False
mysim.Nf=2
mysim.output='distsim_cuby200_3a' 
mysim.run()
