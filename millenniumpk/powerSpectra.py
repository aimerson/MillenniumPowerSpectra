#! /usr/bin/env python

import sys,glob,fnmatch
import numpy as np
import pkg_resources
from .hdf5 import HDF5
from scipy.interpolate import interp1d


class powerSpectrum(object):

    def __init__(self):
        self.k = None
        self.klow = None
        self.kupp = None
        self.pk = None
        return

    def setWavenumbers(self,kmin=1.0e-4,kmax=1.0e3,N=6,logSpacing=True):
        if kmin == 0:
            logSpacing = False
        if logSpacing:
            self.k = 10.0**(np.linspace(np.log10(kmin),np.log10(kmax),2**N+1))
        else:
            self.k = np.linspace(kmin,kmax,2**N+1)
        return

    def setPowerSpectrum(self,k,pk,**kwargs):
        self.pk = interp1d(k,pk,**kwargs)
        return

    def getFourierWindow(self,R,h0=1.0,dimension=3,trunc=None):
        R /= h0
        result = None
        if dimension == 3:
            result = 3.0*(np.sin(R*self.k) - R*self.k*np.cos(R*self.k))
            result /= (R*self.k)**3
        if dimension == 1:
            result = np.sin(R*self.k)/(self.k*R)
        np.place(result,np.isnan(result),0.0)
        if trunc is not None:
            ktrunc = 2.0*Pi/trunc
            mask = self.k < ktrunc
            np.place(result,mask,0.0)
        return result

   def computeSigma8(self,h0=1.0,verbose=False):
        dk = self.k[1] - self.k[0]
        values = self.pk(self.k)*(self.getFourierWindow(8.0,h0=h0)*self.k)**2
        sigma8 = romb(values,dx=dk,show=verbose)/dk
        sigma8 /= 2.0*(Pi**2)
        return np.sqrt(sigma8)

    def _integratePowerSpectrumAtSeparation(self,isep,dimension,h0=1.0,trunc=None,verbose=False,progressOBJ=None):
        r = self.r[isep]
        dk = self.k[1] - self.k[0]
        values = self.k**2*self.pk(self.k)*self.getFourierWindow(r,h0=h0,dimension=dimension,trunc=trunc)
        result = np.copy(romb(values,dx=dk,show=verbose)/(2.0*(Pi**2)))
        self.xi[isep] = np.copy(result)
        del values
        if progressOBJ is not None:
            progressOBJ.increment()
            progressOBJ.print_status_line("r = "+str(r)+" Mpc")
        return

    def computeXi(self,r,h0=1.0,trunc=None,volumeAveraged=False,verbose=False):
        self.r = None
        self.r = r
        self.xi = np.zeros_like(self.r)
        if volumeAveraged:
            dimension = 3
        else:
            dimension = 1
        PROG = Progress(len(self.r))
        dummy = [self._integratePowerSpectrumAtSeparation(i,dimension,h0=h0,trunc=trunc,verbose=verbose,progressOBJ=PROG) \
                     for i in range(len(self.r))]
        return


    
class MillenniumPowerSpectrum(HDF5):
    
    def __init__(self):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        datafile = pkg_resources.resource_filename(__name__,"datasets/millenniumWMAP1_matterPowerSpectra.hdf5")
        # Initalise HDF5 class
        super(MillenniumPowerSpectra, self).__init__(datafile,'r')
        self.redshift = np.array(self.fileObj["redshift"])[::-1]
        self.snapshot = np.array(self.fileObj["snapshot"])[::-1]
        self.k = np.array(self.fileObj["Snapshots/000/k"])
        self.zPk = np.zeros((len(self.k),len(self.redshift)),dtype=float)
        HDF.close()
        return

    def readSnapshot(self,snapNumber,property):
        hdfdir = "Snapshots/"+str(snapNumber).zfill(3)+"/"        
        if property not in self.lsDatasets(hdfdir):
            raise KeyError("Dataset '"+property+"' not found in snapshot directory!")
        return np.array(self.fileObj[hdfdir+property])

    def storeSnapshot(self,snapNumber,property):        
        iz = np.argmin(np.fabs(self.snapshot-snapNumber))
        data = self.readSnapshot(snapNumber,property)
        self.zPk[:,iz] = np.copy(np.log10(data))
        return

    def read(self,property):
        self.zPk = np.zeros((len(self.k),len(self.redshift)),dtype=float)
        dummy = [self.readSnapshot(iz,property) for iz in self.snapshot]
        HDF.close()
        return
        
    def getPowerSpectrumAtRedshift(self,z,property="nonLinearPk",snapToSnapshot=False,**kwargs):
        if snapToSnapshot:
            iz = np.argmin(np.fabs(self.redshift-z))
            HDF = HDF5(self.datafile,'r')
            pk = self.readSnapshot(HDF,self.snapshot[iz],property)
            HDF.close()
        else:            
            self.read(property)
            f = interp1d(self.redshift,self.zPk,**kwargs)
            pk = 10.0**f(z)
        return pk
    
    def getCorrelationFunctionAtRedshift(self,z,r,snapToSnapshot=True,volumeAveraged=False,N=14,trunc=None,**kwargs):
        PK = powerSpectrum()
        PK.setWavenumbers(self.k.min(),self.k.max(),N=N)
        pk = self.getPowerSpectrumAtRedshift(z,snapToSnapshot=snapToSnapshot)
        PK.setPowerSpectrum(self.k,pk,**kwargs)
        PK.computeXi(r,volumeAveraged=volumeAveraged,trunc=trunc)
        return PK.xi
    

