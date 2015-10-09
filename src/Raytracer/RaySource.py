'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import Ray
import matplotlib.pyplot as mpl

class RaySource(object):
    def __init__(self, numRays = 1):
        """ Defines a ray source with a set of rays that can be traced.
        Contains a list of rays.
        
        Input
        numRays: Number of rays in the list
        """
        self.numRays = numRays
        self.rays = []
        
    def generateRays(self):
        for rn in range(self.numRays):
            r = Ray.Ray()
            self.rays.append(r)
            
    def drayRays(self, fig):
        for ray in self.rays:
            ray.drawRay(fig)
        
        
class Collimated1DSource(RaySource):
    def __init__(self, numRays = 1000, xDim = 1e-3, l = 567e-9, W = 1):
        """ Defines an 1D ray source with a set of rays that can be traced. The rays
        are collimated in the z direction and initialized along the x axis with
        an extent of xDim.
        Contains a list of rays.
        
        Input
        numRays: Number of rays in the list
        xDim: The rays are evenly spaced around -xDim/2, +xDim/2
        l: Wavelength
        W: Energy
        
        """
        super(Collimated1DSource, self).__init__(numRays)
        self.generateRays(xDim, l, W)
        
    def generateRays(self, xDim = 5e-3, l = 567e-9, W = 1):
        self.rays = []
        for rn in range(self.numRays):
            x = np.array([xDim * (-0.5 + 1.0*rn/self.numRays), 0, 0, 1])
            xp = np.array([0,0,1,0])
            r = Ray.Ray(x, xp, l, W)
            self.rays.append(r)