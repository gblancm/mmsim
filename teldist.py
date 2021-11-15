import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
from polynomial2d import polyfit2d
from numpy.polynomial.polynomial import polyval2d
from scipy.optimize import curve_fit


# Read ZEMAX spot field angle and centroid table
tab=ascii.read('./zemax/junk.txt')
fx=tab['col1'].data
fy=tab['col2'].data
xcen=tab['col5'].data
ycen=tab['col6'].data

# Fit transformation xcen(fx,fy) ; ycen(fx, fy)

deg=1
cx=polyfit2d(fx, fy, xcen, deg=deg)
cy=polyfit2d(fx, fy, ycen, deg=deg)
#xcenfit=polyval2d(fx,fy,cx)
#ycenfit=polyval2d(fx,fy,cy)

# Transformation using global plate scale
globscale=np.max(xcen)/np.max(fx) # global plate scale [mm/deg]
globscaley=np.max(ycen)/np.max(fy) # global plate scale [mm/deg]
cxglob=np.array([[0,0],[globscale,0]])
cyglob=np.array([[0,globscale],[0,0]])
xcenfit=polyval2d(fx,fy,cxglob)
ycenfit=polyval2d(fx,fy,cyglob)

# Check radial residuals
rf=np.sqrt(fx**2+fy**2)
rcen=np.sqrt(xcen**2+ycen**2)
rcenfit=np.sqrt(xcenfit**2+ycenfit**2)

def func(x, p0, p1, p2, p3, p4, p5):
  return p0+p1*x+p2*x**2+p3*x**3+p4*x**4+p5*x**5

def oddpol1(x,p1):
    return p1*x
def oddpol3(x,p1,p3):
    return p1*x+p3*x**3
def oddpol5(x,p1,p3,p5):
    return p1*x+p3*x**3+p5*x**5
def oddpol7(x,p1,p3,p5,p7):
    return p1*x+p3*x**3+p5*x**5+p7*x**7
def pol4(x, p1, p2, p3, p4):
  return p1*x+p2*x**2+p3*x**3+p4*x**4

p0=[-1e-6]
cr1,_=curve_fit(oddpol1, rcen, (rcenfit-rcen), p0=p0)
p0=[-1e-6,0]
cr3,_=curve_fit(oddpol3, rcen, (rcenfit-rcen), p0=p0)
p0=[-1e-6,0,0]
cr5,_=curve_fit(oddpol5, rcen, (rcenfit-rcen), p0=p0)
p0=[-1e-6,0,0,0]
cr7,_=curve_fit(oddpol7, rcen, (rcenfit-rcen), p0=p0)
p0=[-1e-6,0,0,0]
fcr4,_=curve_fit(pol4, rcen, (rcenfit-rcen), p0=p0)


#p0=[0,-1e-6,0,0,0,0]
# Odd polynomial fits
#cr1,_=curve_fit(func, rcen, (rcenfit-rcen), p0=p0, bounds=([-1e-30, -np.inf, -1e-30, -1e-30, -1e-30, -1e-30], [1e-30, np.inf, 1e-30, 1e-30, 1e-30, 1e-30]))
#cr3,_=curve_fit(func, rcen, (rcenfit-rcen), p0=p0, bounds=([-1e-30, -np.inf, -1e-30, -np.inf, -1e-30, -1e-30], [1e-30, np.inf, 1e-30, np.inf, 1e-30, 1e-30]))
#cr5,_=curve_fit(func, rcen, (rcenfit-rcen), p0=p0, bounds=([-1e-30, -np.inf, -1e-30, -np.inf, -1e-30, -np.inf], [1e-30, np.inf, 1e-30, np.inf, 1e-30, np.inf]))


# Full polynomial fits
#fcr4,_=curve_fit(func, rcen, (rcenfit-rcen), p0=p0, bounds=([-1e-30, -np.inf, -np.inf, -np.inf, -np.inf, -1e-30], [1e-30, np.inf, np.inf, np.inf, np.inf, 1e-30]))




# Draw grid and transformed grid
aux1fx=np.linspace(-1.5, 1.5, 100)
aux1fy=np.linspace(-1.5, 1.5, 100)
aux2fx=np.linspace(-1.5, 1.5, 19)
aux2fy=np.linspace(-1.5, 1.5, 19)
auxfxx, auxfyx=np.meshgrid(aux2fx, aux1fy)
auxfxy, auxfyy=np.meshgrid(aux1fx, aux2fy)

auxxcenx=polyval2d(auxfxx, auxfyx, cx)
auxycenx=polyval2d(auxfxx, auxfyx, cy)
auxxceny=polyval2d(auxfxy, auxfyy, cx)
auxyceny=polyval2d(auxfxy, auxfyy, cy)

# Make Plots

# Radial residuals
auxr=np.linspace(np.min(rcen), np.max(rcen), 100)
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10, 10), dpi=200)
ax1.plot(rcen, rcenfit-rcen,'o')
#ax1.plot(auxr, auxrres, color='red', label='Fit order = '+"{:.0f}".format()
ax1.plot(auxr, oddpol1(auxr, *cr1), label='Odd Polynomial, order = 1', color='orange')
ax1.plot(auxr, oddpol3(auxr, *cr3), label='Odd Polynomial, order = 3', color='green')
ax1.plot(auxr, oddpol5(auxr, *cr5), label='Odd Polynomial, order = 5', color='red')
ax1.plot(auxr, oddpol7(auxr, *cr7), label='Odd Polynomial, order = 7', color='pink')
ax1.plot(auxr, pol4(auxr, *fcr4), label='Full Polynomial, order = 4', color='brown')
ax1.axhline(0, linestyle='--', color='black', alpha=0.5) 
ax1.set_xlabel('r - Spot Plane [mm]', fontsize=20)
ax1.set_ylabel(r'$\Delta$r - Spot Plane [mm]', fontsize=20)
ax1.legend()


#ax2.plot(rcen, (rcenfit-rcen)-oddpol1(rcen, *cr1), 'o', label='Odd Polynomial, order = 1', color='orange')
#ax2.plot(rcen, (rcenfit-rcen)-oddpol3(rcen, *cr3), 'o', label='Odd Polynomial, order = 3', color='green')
ax2.plot(rcen, (rcenfit-rcen)-oddpol5(rcen, *cr5), 'o', label='Odd Polynomial, order = 5', color='red')
ax2.plot(rcen, (rcenfit-rcen)-oddpol7(rcen, *cr7), 'o', label='Odd Polynomial, order = 7', color='pink')
ax2.plot(rcen, (rcenfit-rcen)-pol4(rcen, *fcr4), 'o', label='Full Polynomial, order = 4', color='brown')
ax2.axhline(0, color='black', alpha=0.5) 
ax2.axhline(0.01*np.max(xcen)/np.max(fx)/3600, linestyle=':', color='black', alpha=0.5) 
ax2.axhline(-0.01*np.max(xcen)/np.max(fx)/3600, linestyle=':', color='black', alpha=0.5) 
ax2.axhline(0.1*np.max(xcen)/np.max(fx)/3600, linestyle=':', color='black', alpha=0.5) 
ax2.axhline(-0.1*np.max(xcen)/np.max(fx)/3600, linestyle=':', color='black', alpha=0.5) 
ax2.set_xlabel('r - Spot Plane [mm]', fontsize=20)
ax2.set_ylabel(r'$\Delta$r - Fit - Spot Plane [mm]', fontsize=20)

plt.savefig('teldist_rres.png')


# Cartesian residuals
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(25, 15), dpi=200)
ax1.plot(xcen, xcenfit-xcen,'o')
ax1.axhline(0, linestyle='--', color='black', alpha=0.5) 
ax1.set_xlabel('X - Spot Plane [mm]', fontsize=20)
ax1.set_ylabel(r'$\Delta$X - Spot Plane [mm]', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax2.plot(xcen, ycenfit-ycen,'o')
ax2.axhline(0, linestyle='--', color='black', alpha=0.5) 
ax2.set_xlabel('X - Spot Plane [mm]', fontsize=20)
ax2.set_ylabel(r'$\Delta$Y - Spot Plane [mm]', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax3.plot(ycen, xcenfit-xcen,'o')
ax3.axhline(0, linestyle='--', color='black', alpha=0.5) 
ax3.set_xlabel('Y - Spot Plane [mm]', fontsize=20)
ax3.set_ylabel(r'$\Delta$X - Spot Plane [mm]', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax4.plot(ycen, ycenfit-ycen,'o')
ax4.axhline(0, linestyle='--', color='black', alpha=0.5) 
ax4.set_xlabel('Y - Spot Plane [mm]', fontsize=20)
ax4.set_ylabel(r'$\Delta$Y - Spot Plane [mm]', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
plt.savefig('teldist2.png')




#factor=2000 # scale factor for residual arrows
#scale=1 # size of scale marker in arcsec

factor=200 # scale factor for residual arrows
scale=10 # size of scale marker in arcsec

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,10), dpi=200)

ax1.plot(auxfxx, auxfyx,'.', color='orange', alpha=0.2)
ax1.plot(auxfxy, auxfyy,'.', color='orange', alpha=0.2)
ax1.plot(fx, fy, 'o', markersize=10)
ax1.set_xlabel(r'$\theta_{x}$ [deg]', fontsize=20)
ax1.set_ylabel(r'$\theta_{y}$ [deg]', fontsize=20)
ax1.set_title('Field Angle', fontsize=25)
ax1.set_aspect('equal')
ax1.grid(color='grey',linestyle='--')
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

ax2.plot(auxxcenx, auxycenx,'.', color='orange', alpha=0.2)
ax2.plot(auxxceny, auxyceny,'.', color='orange', alpha=0.2)
ax2.plot(xcen, ycen, 'o', markersize=10)
#ax2.plot(xcenfit, ycenfit, 'o', color='red', markersize=5)
for i in range(len(xcen)):
    if (xcenfit[i]!=xcen[i])+(ycenfit[i]!=ycen[i]):
        ax2.arrow(xcen[i], ycen[i], (xcenfit[i]-xcen[i])*factor, (ycenfit[i]-ycen[i])*factor, width=5, head_width=20, length_includes_head=True, head_starts_at_zero=True, color='red')
#    ax2.arrow(xcen[i], ycen[i], (xcenfit[i]-xcen[i])*factor, (ycenfit[i]-ycen[i])*factor, width=5, head_width=20, length_includes_head=False, head_starts_at_zero=True, color='red')
ax2.plot(np.array([-500, -500+np.max(xcen)/np.max(fx)/3600*factor*scale]), np.array([550, 550]), color='black', linewidth=3)
ax2.text(-500+np.max(xcen)/np.max(fx)/3600*factor*scale/4, 570, "{:.0f}".format(scale)+'\"', fontsize=20)
ax2.set_xlabel('X - [mm]', fontsize=20)
ax2.set_ylabel('Y - [mm]', fontsize=20)
ax2.set_title('Spot Plane', fontsize=25)
ax2.set_aspect('equal')
ax2.grid(color='grey',linestyle='--')
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
plt.savefig('teldist.png')
