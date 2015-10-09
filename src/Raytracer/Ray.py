'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import matplotlib.pyplot as mpl 

class Ray(object):
    def __init__(self, x = np.array([0,0,0]), xp = np.array([0,0,1]), n = 1.0, l = 564e-9, W = 1.0):
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
        self.l = l
        
    def addPos(self, newX, newXp, newN, t = 1.0):
        """ Add new node in the trace of the ray.
        
        Inputs
        newX : New position, 3 or 4 element numpy array
        newXp: New direction, 3 or 4 element numpy array
        newN: New refractive index
        t: Transmission, fraction of the energy remaining
        
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
        
    def drawRay(self, fig):
        mpl.figure(fig)
        x = np.array(self.x)
        mpl.plot(x[:,2], x[:,0])

        