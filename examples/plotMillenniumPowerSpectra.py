#! /usr/bin/env python

import numpy as np
from millenniumpk.plotting.utils import *
from matplotlib.cm import ScalarMappable
from millenniumpk.powerSpectra import MillenniumPowerSpectrum,Pi


MPK = MillenniumPowerSpectrum()
redshifts = MPK.redshift
mask = np.logical_and(redshifts>=0.7,redshifts<=2.2)

carr = colour_array(len(redshifts[mask]),cmap='cool_r')
norm = mpl.colors.Normalize(vmin=redshifts[mask][0],vmax=redshifts[mask][-1])


fig = figure(figsize=(4.5,9))
ax = fig.add_subplot(211,xscale='log',yscale='log')
cax = fig.add_axes([0.35, 0.865, 0.52, 0.02])
cmap = plt.cm.cool
sm = ScalarMappable(norm=norm,cmap=cmap)
sm.set_array([])  
cbar = fig.colorbar(sm,cax=cax,orientation='horizontal')
cbar.set_label("z")
cbar.set_ticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0])
ic = 0
kmin = 9999.9
kmax = -9999.9
for z in redshifts:
    if z >= 0.7 and z <= 2.2: 
        k = MPK.k
        pk = MPK.getPowerSpectrumAtRedshift(z,property="nonLinearPk")
        ax.plot(k,pk,c=carr[ic],label="z = "+sigfig(z,3),lw=1.5)
        ic += 1
        kmin = np.minimum(kmin,k[0])
        kmax = np.maximum(kmax,k[-1])
ax.set_xlim(kmin,kmax)
rlim = 100.0
klim = 2*Pi/rlim
ax.axvline(klim,ls=':',c='k',lw=1.5)
ax.text(klim,0.1,"$r\,=\,"+str(sigfig(rlim,1))+"\,h^{-1}{\\rm Mpc}$",rotation=90.0,ha='right',va='bottom')
ax.set_xlabel("$k\,\,\left [h\,{\\rm Mpc}^{-1}\\right ]$")
ax.set_ylabel("$P(k)\,\,\left [h^3{\\rm Mpc}^{-3}\\right ]$")

r = 10.0**(np.linspace(np.log10(0.1),np.log10(100.0),90))
ic = 0
volumeAveraged = False
dk = 0.01
ax = fig.add_subplot(212,xscale='log',yscale='log')
cax = fig.add_axes([0.35, 0.43, 0.52, 0.02])
cmap = plt.cm.cool
sm = ScalarMappable(norm=norm,cmap=cmap)
sm.set_array([])  
cbar = fig.colorbar(sm,cax=cax,orientation='horizontal')
cbar.set_label("z")
cbar.set_ticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0])
for z in redshifts:
    if z >= 0.7 and z <= 2.2: 
        xi = MPK.getCorrelationFunctionAtRedshift(z,r,snapToSnapshot=True,volumeAveraged=volumeAveraged,N=14,\
                                                       kind='linear',fill_value='extrapolate')
        ax.plot(r,xi,c=carr[ic],lw=1.5)
        ic += 1
ax.set_ylabel("$\\xi_{\\rm mm}(r)$")
ax.set_xlabel("$r\,\,[h^{-1}{\\rm Mpc}]$")
ax.set_xlim(r[0],r[-1])
ax.set_ylim(top=2)
savefig("millenniumPowerSpectraCorrelationFunction.pdf",bbox_inches='tight')
