""" 
Class to handle a set of particles for vortex particles methods

"""
import numpy as np
##

class Particles():
    # constructor
    def __init__(self,nPart,nDim):
        self.nDim  = nDim
        self.nPart = nPart
        self.reset(nPart)
        return
        
    # initializing all particles values to zero
    def reset(self,nPart= None): 
        if nPart is not None:
            self.nPart = nPart
        
        self.P = np.zeros((self.nPart,self.nDim))
        if (self.nDim == 3):
            self.Intensity = np.zeros((self.nPart,self.nDim))
        else:
            self.Intensity = np.zeros(self.nPart)
        
        self.SmoothParam = np.zeros(self.nPart)
        self.SmoothModel = - 1
        self.KernelOrder = - 1
        self.Volume = np.zeros(self.nPart)
        
    def removeZeroVorticity(o):
        I = np.abs(o.Intensity) > 1e-15
        nNew = sum(I)
        nOld = o.nPart
        o.nPart = nNew
        print(o.P.shape)
        print(o.Intensity.shape)
        o.P = o.P[I,:]
        if (o.nDim == 3):
            o.Intensity = o.Intensity[I,:]
        else:
            o.Intensity = o.Intensity[I]
        
        o.Volume = o.Volume[I]
        o.SmoothParam = o.SmoothParam[I]
        print('Removing %d part with zero vorticity\n' % (nOld - nNew))
        
        
    # --------------------------------------------------------------------------------
    # --- Mutator
    # --------------------------------------------------------------------------------
    def setVolume(o,vol):
        if np.asarray(vol).size == 1:
            o.Volume = vol
        elif (vol.shape != o.Volume.shape):
            raise Exception('Wrong Volume size?')
        o.Volume = vol
        
        
    def setIntensity(o,Intensity):
        if Intensity.shape != o.Intensity.shape:
            raise Exception('Wrong Intensity size?')
        
        o.Intensity = Intensity
        
    def setP(o,P): 
        if (P.shape != o.P.shape):
            raise Exception('Wrong P size?')
        
        o.P[:o.nPart,:o.nDim] = P

        #     # Rotating particles
    # #     Xp=Part(1,:);
    # #     Yp=Part(2,:);
    # #     theta=pi/6;
    # #     Part(1,:)=Xp*cos(theta)-Yp*sin(theta);
    # #     Part(2,:)=Xp*sin(theta)+Yp*cos(theta);
