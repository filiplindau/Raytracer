'''
Created on 10 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import OpticalElement as oe
import OpticalMaterial as om

class OpticalSystem(object):
    def __init__(self):
        self.elements = []
        self.raySourceList = []
        self.materialLibrary = om.MaterialLibrary()
        self.opticalAxisM = np.identity(4) 
        self.opticalAxisTheta = 0.0
        self.opticalAxisPhi = 0.0
        
    def addElement(self, element):
        element.rotateElement(self.opticalAxisTheta, self.opticalAxisPhi)
        self.elements.append(element)
        
    def rotateOpticalAxis(self, theta, phi):
        self.opticalAxisTheta = theta
        self.opticalAxisPhi = phi
        thM = np.array([[1.0, 0.0,            0.0,           0.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 0.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 1.0], 
                         [0.0,          1.0, 0.0,         1.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 1.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        self.opticalAxisM = np.dot(thM, phM)
        
        
    def addRaySource(self, raySource):
        self.raySourceList.append(raySource)
        
    def traceSystem(self):
        if self.raySourceList != []:
            for raySource in self.raySourceList:
                for (ind, element) in enumerate(self.elements):    
                    raysT = raySource.rays.copy()
                    raysT[:,1,:] = np.transpose(np.dot(self.opticalAxisM, np.transpose(raySource.rays[:,1,:])))        
                    raysList = element.propagateRays(raysT)
                    for rays in raysList:
                        raysT = rays.copy()
                        raysT[:,1,:] = np.transpose(np.dot(np.transpose(self.opticalAxisM), np.transpose(raySource.rays[:,1,:])))
                        raySource.updateRays(raysT) 