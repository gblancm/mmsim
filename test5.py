# This script calculates differential field distorsion change rates [arcsec/hr]
# over a grid in declination and hour angle

import mmsim
import numpy as np
import matplotlib.pyplot as plt

mysim=mmsim.distsim()
mysim.dfdfile='1536revA_dfd.pkl' # set differential field distorsion file 
mysim.simtype='atmostel'
mysim.Nw=3
mysim.Nt=3

decarr=np.linspace(0, -60, 5)
haarr=np.linspace(-30, 30, 5)
grid=np.meshgrid(decarr, haarr)
dec=grid[0]
ha=grid[1]
rate=np.empty_like(dec)


# 1h exp at dec=-29 at a range of starting HA
mysim.obs.texp=900
for i in range(len(decarr)):
    for j in range(len(haarr)):
        mysim.obs.decstart=dec[i,j]
        mysim.obs.hastart=ha[i,j]
        mysim.run()
        rate[i,j]=mysim.dfdrate



fig1,ax1 = plt.subplots(figsize=(15,15))
ax1.contourf(ha, dec, rate) 


