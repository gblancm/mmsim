# mmsim

Field distorion and ADR simulation code:

Codingf plan:

Define a class "distsim" so I can do:

> mysim=distsim()

and then setits attributes. The class has to have attributes related to site location 
(needed for atmospheric refraction models):

mysim.site.latitude = site latitude
mysim.site.elevation = site elevation
mysim.site.temperature = average site temperature during observation
...

the class also has attributes related to the observation (needed to trace dra(lambda),dha(lambda) vs time):

mysim.obs.dec = declination
mysim.obs.ha = hour angle at mid-exposure
mysim.obs.pa = position angle of observation
mysim.obs.texp = exposure time
mysim.obs.seeing = seeing FWHM

...

the class also has attributes related to the sampling points positions (i.e. objects) in the FOV [in deg] at the
begining of the exposure. I will assume the focal plane is at a fixed PA, so it rotates during the exposure.

mysim.obj.coords0 = x(t=0), y(t=0)

the class also has attributes related to the spectroscopic apertures (in this case fibers, but could be slits)

mysim.aper.fiber.r = fiber radius
...

or 

mysim.aper.slit.dx = slit width
mysim.aper.slit.dy = slit height
mysim.aper.slit.pa = slit position angle
...

Finally, the class needs a method to run the simulation:

> mysim.run()

In which the position in the focal plane of each object is calculated as a function of time and wavelength

mysim.obj.coords[i:j]  ; where i,j run over time,wavelength


CODE:

python 3
astropy.Coordinates
AstroAtmosphere



