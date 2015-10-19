'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import matplotlib.pyplot as mpl 



class RayStore(object):
    def __init__(self, x = np.array([0,0,0]), xp = np.array([0,0,1]), n = 1.0, l = 564e-9, W = 1.0, ng = 1.0, color = (0,0,0.85)):
        """ Initialize ray with position, direction, wavelength, and energy 
        
        Inputs 
        x: Position, 3 or 4 element numpy array
        xp: direction, 3 or 4 element numpy array
        n: starting refractive index
        l: wavelength
        W: energy
        
        """
        if x.shape[0] == 3:
            x = np.hstack((x,1.0))
        else:
            x[3] = 1.0
        if xp.shape[0] == 3:
            xp = np.hstack((xp,0.0))
        else:
            xp[3] = 0.0
        self.x = [x]
        self.xp = [xp]
        self.W = [W]
        self.n = [n]
        self.ng = [ng]
        self.distance = [0.0]
        self.time = [0.0]
        self.l = l
        
        self.color = color
        
    def addPos(self, newX, newXp, newN, newNg=1.0, t = 1.0, W=1.0, distance = 0.0):
        """ Add new node in the trace of the ray.
        
        Inputs
        newX : New position, 3 or 4 element numpy array
        newXp: New direction, 3 or 4 element numpy array
        newN: New refractive index
        t: Transmission, fraction of the energy remaining
        distance = length of ray segment
        """
        if newX.shape[0] == 3:
            newX = np.hstack((newX,1.0))
        else:
            newX[3] = 1.0
        if newXp.shape[0] == 3:
            newXp = np.hstack((newXp,0.0))
        else:
            newXp[3] = 0.0
        self.x.append(newX)
        self.xp.append(newXp)
        self.W.append(self.W[-1]*t)
        self.n.append(newN)
        self.ng.append(newNg)
        self.distance.append(distance)
        self.time.append(self.time[-1] + distance*newNg/299792458.0)
        
    def drawRay(self, fig):
        mpl.figure(fig)
        x = np.array(self.x)
        mpl.plot(x[:,2], x[:,0], color=self.color)
        
    def getArrayPoints(self):
        return np.array(self.x)
    

        