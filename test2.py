# This script calculates atmospheric and instrumental field distorsion maps
# for a field at dec=-29deg, assuming a 1hr and 0.5hr exposures
# that starts at hour angles of -60, -50, -40 ,-30 , and -20 degrees. 
# It runs simulations with ADC correction.

import mmsim

mysim=mmsim.distsim()


# Zenith 1h exp at range of starting HA (with ADC)

mysim.obs.decstart=-29
mysim.ADC=True

mysim.obs.texp=1800

mysim.obs.hastart=-60
mysim.output='distsim_m60_1hr'
mysim.run()

mysim.obs.hastart=-50
mysim.output='distsim_m50_1hr'
mysim.run()

mysim.obs.hastart=-40
mysim.output='distsim_m40_1hr'
mysim.run()

mysim.obs.hastart=-30
mysim.output='distsim_m30_1hr'
mysim.run()

mysim.obs.hastart=-20
mysim.output='distsim_m20_1hr'
mysim.run()

mysim.obs.hastart=-10
mysim.output='distsim_m10_1hr'
mysim.run()

mysim.obs.texp=1800

mysim.obs.hastart=-60
mysim.output='distsim_m60_1800s'
mysim.run()

mysim.obs.hastart=-50
mysim.output='distsim_m50_1800s'
mysim.run()

mysim.obs.hastart=-40
mysim.output='distsim_m40_1800s'
mysim.run()

mysim.obs.hastart=-30
mysim.output='distsim_m30_1800s'
mysim.run()

mysim.obs.hastart=-20
mysim.output='distsim_m20_1800s'
mysim.run()

mysim.obs.hastart=-10
mysim.output='distsim_m10_1800s'
mysim.run()

