# This script calculates atmospheric and instrumental field distorsion maps
# for three declinations (+1deg, -29deg, -59deg), assuming a 2hr exposure time 
# that starts at an hour angle of -1 hr. It runs simulations with and without ADC
# correction.

import mmsim

mysim=mmsim.distsim()  # initialize distorsion simulation object

mysim.obs.hastart=-15 # set starting hour angle in degrees
mysim.obs.texp=7200 # set exposure time in seconds

mysim.dfdfile='1536revA_dfd.pkl' # set differential field distorsion file 



# 30 deg South of zenith
mysim.obs.decstart=-59 # set declination of field center
mysim.simtype='atmos'
mysim.output='distsim_south'
mysim.run()
mysim.simtype='atmoscorr'
mysim.output='distsim_south_adc'
mysim.run()
mysim.simtype='atmostel'
mysim.output='distsim_south_tel'
mysim.run()


# 30 deg North of zenith
mysim.obs.decstart=1
mysim.simtype='atmos'
mysim.output='distsim_north'
mysim.run()
mysim.simtype='atmoscorr'
mysim.output='distsim_north_adc'
mysim.run()
mysim.simtype='atmostel'
mysim.output='distsim_north_tel'
mysim.run()


# Zenith
mysim.obs.decstart=-29
mysim.simtype='atmos'
mysim.output='distsim_zenith'
mysim.run()
mysim.simtype='atmoscorr'
mysim.output='distsim_zenith_adc'
mysim.run()
mysim.simtype='atmostel'
mysim.output='distsim_zenith_tel'
mysim.run()



