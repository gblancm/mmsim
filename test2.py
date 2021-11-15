# This script calculates atmospheric and instrumental field distorsion maps
# for a field at dec=-29deg, assuming a 1hr and 0.5hr exposures
# that starts at hour angles of -60, -50, -40 ,-30 , and -20 degrees. 
# It runs simulations with ADC correction.

import mmsim

mysim=mmsim.distsim()

mysim.dfdfile='1536revA_dfd.pkl' # set differential field distorsion file 
mysim.simtype='atmostel'


# 1h exp at dec=-29 at a range of starting HA
mysim.obs.texp=3600
mysim.obs.decstart=-29

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

# 30 min exp at dec=-29 at a range of starting HA
mysim.obs.texp=1800
mysim.obs.decstart=-29

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

