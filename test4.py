import mmsim

mysim=mmsim.distsim()

mysim.obs.hastart=-15 # set starting hour angle in degrees
mysim.obs.texp=7200 # set exposure time in seconds

# 30 deg South of zenith
mysim.obs.decstart=-59 # set declination of field center
mysim.simtype='atmostel'
mysim.dfdfile='1536revA_dfd.pkl'
mysim.output='distsim_south_tel'
mysim.Nw=5
mysim.Nt=11
mysim.run()