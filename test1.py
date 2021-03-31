# This script calculates atmospheric and instrumental field distorsion maps
# for three declinations (+1deg, -29deg, -59deg), assuming a 2hr exposure time 
# that starts at an hour angle of -1 hr. It runs simulations with and without ADC
# correction.

import mmsim

mysim=mmsim.distsim()  # initialize distorsion simulation object

mysim.obs.hastart=-15 # set starting hour angle in degrees
mysim.obs.texp=7200 # set exposure time in seconds

# 30 deg South of zenith
mysim.obs.decstart=-59 # set declination of field center
mysim.ADC=True
mysim.output='distsim_south_adc'
mysim.run()

mysim.ADC=False
mysim.output='distsim_south'
mysim.run()

# 30 deg North of zenith
mysim.obs.decstart=1
mysim.ADC=True
mysim.output='distsim_north_adc'
mysim.run()

mysim.ADC=False
mysim.output='distsim_north'
mysim.run()

# Zenith
mysim.obs.decstart=-29
mysim.ADC=True
mysim.output='distsim_zenith_adc'
mysim.run()

mysim.ADC=False
mysim.output='distsim_zenith'
mysim.run()


