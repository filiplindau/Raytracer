'''
Created on 13 Oct 2015

@author: Filip Lindau
'''

import numpy as np

class OpticalAperture(object):
    def __init__(self, size=12.7e-3, absorb=True):
        ''' Create an optical aperture. The default implementation is an 
        infinite aperture where the size is specified.
        
        Input:
        size: radius of the circular aperture
        absorb: if True points outside the aperture are absorbed, otherwise just transmitted
        '''
        self.generateEdge()        
        
    def pointInAperture(self, x):
        ''' Calculate if a point is inside of an aperture. Assumes the point is in 
        local coordinates of the surface.
        
        Input:
        x: 4 element position vector 
        '''
        return True
        
    def generateEdge(self):
        self.edge = np.transpose(np.array([np.array([]), np.array([]), np.array([]), np.array([])]))
        
    def getEdge(self):
        return self.edge
        
class CircularAperture(OpticalAperture):
    def __init__(self, size=12.7e-3, absorb=True):
        self.size = size
        self.r2 = size ** 2
        self.absorb = absorb
        OpticalAperture.__init__(self, size, absorb)

    def generateEdge(self):
        nbrPoints = 16
        theta = np.linspace(0, 2 * np.pi, nbrPoints)
        self.edge = np.transpose(np.vstack((self.size * np.cos(theta), self.size * np.sin(theta), np.zeros(nbrPoints), np.ones(nbrPoints))))
        
    def pointInAperture(self, x):
        ''' Calculate if a point is inside of an aperture. Assumes the point is in 
        local coordinates of the surface.
        
        Input:
        x: 4 element position matrix (first index is ray number) 
        '''
        return x[:, 0] * x[:, 0] + x[:, 1] * x[:, 1] < self.r2
        
class RectangularAperture(OpticalAperture):
    def __init__(self, size, absorb=True):
        ''' Create a rectangular aperture.
        
        Input:
        size: Two element list, first element size in x direction, second size in y direction
        '''
        self.x = size[0] / 2
        self.y = size[1] / 2
        self.absorb = absorb
        print "x, y: ", self.x, self.y
        OpticalAperture.__init__(self, size, absorb)

    def generateEdge(self):
        x0 = np.array([-self.x, -self.y, 0, 1])
        x1 = np.array([-self.x, self.y, 0, 1])
        x2 = np.array([self.x, self.y, 0, 1])
        x3 = np.array([self.x, -self.y, 0, 1])
        x4 = np.array([-self.x, -self.y, 0, 1])        
        self.edge = (np.array([x0, x1, x2, x3, x4]))
#        print "edge: ", self.edge
#        raise ValueError
        
    def pointInAperture(self, x):
        return np.logical_and(np.abs(x[:, 0]) < self.x, np.abs(x[:, 1]) < self.y)

class InifiniteAperture(OpticalAperture):
    def __init__(self):
        OpticalAperture.__init__(self)
    
    def pointInAperture(self, x):
        return np.tile(True, x.shape[0])
