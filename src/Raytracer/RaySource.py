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
            
    def drawRays(self, fig):
        for ray in self.rays:
            ray.drawRay(fig)
        
    def getRayPoints(self):
        points = []
        for ray in self.rays:
            points.append((ray.getArrayPoints(), ray.color))            
        return points
        
class Collimated1DSource(RaySource):
    def __init__(self, numRays = 1000, xDim = 1e-3, l = 567e-9, W = 1, color = (0,0,0.85)):
        """ Defines an 1D ray source with a set of rays that can be traced. The rays
        are collimated in the z direction and initialized along the x axis with
        an extent of xDim.
        Contains a list of rays.
        
        Input
        numRays: Number of rays in the list
        xDim: The rays are evenly spaced around -xDim/2, +xDim/2
        l: Wavelength
        W: Energy
        color: Drawing color of rays 
        """
        super(Collimated1DSource, self).__init__(numRays)
        self.color = color
        self.generateRays(xDim, l, W, color)
        
    def generateRays(self, xDim = 5e-3, l = 567e-9, W = 1, color = (0,0,0.85)):
        ''' Generate rays matrix and a list of RayStores. A RayStore keeps track
        of a ray as it moves through the optical system.
        
        rays is 3d matrix:
        index0: ray number
        index1: data type (0=x, 1=xp, 2=optical data)
        index2: 4 element data vector
        
        The optical data is encoded as (wavelength, n, ng, W)
        
        Examples: 
        rays[10, 0, :]... position vector of ray 10 
        rays[ :, 1, :]... direction vectors of all rays
        rays[15, 2, 1]... current refractive index for ray 15
        '''
        xx = np.linspace(-xDim/2, xDim/2, self.numRays)
        x = np.column_stack((xx, np.zeros(self.numRays), np.zeros(self.numRays), np.ones(self.numRays)))
        xp = np.column_stack((np.zeros(self.numRays), np.zeros(self.numRays), np.ones(self.numRays), np.zeros(self.numRays)))
        data = np.column_stack((l*np.ones(self.numRays), 1.0*np.ones(self.numRays), 1.0*np.ones(self.numRays), W*np.ones(self.numRays)))
        self.rays = np.dstack((x,xp,data)).swapaxes(1,2)
        self.rayStoreList = []
        for rn in range(self.numRays):
            x=self.rays[rn,0,:]
            xp = self.rays[rn,1,:]
            l = self.rays[rn,2,0]
            n = self.rays[rn,2,1]
            W = self.rays[rn,2,3]
            r = Ray.RayStore(x=x, xp=xp, n=n, l=l, W=W, color=color)
            self.rayStoreList.append(r)
            
    def updateRays(self, newRays):
        ''' Update ray data with new rays (in global coordinates)
        Sets rays matrix to newRays matrix and updates the RayStoresList
        '''
        d = np.sqrt(np.sum((newRays[:,0,0:3]-self.rays[:,0,0:3])**2))
        self.rays = newRays
        for rn in range(self.numRays):
            x=self.rays[rn,0,:]
            xp = self.rays[rn,1,:]
            l = self.rays[rn,2,0]
            n = self.rays[rn,2,1]
            ng = self.rays[rn,2,2]
            W = self.rays[rn,2,3]
            self.rayStoreList[rn].addPos(x=x, xp=xp, l=l, n=n, ng=ng, W=W, distance=d)
            