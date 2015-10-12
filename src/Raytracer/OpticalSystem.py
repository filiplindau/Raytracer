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
        self.raySource = None
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
        
        
    def setRaySource(self, raySource):
        self.raySource = raySource
        
    def traceSystem(self):
        if self.raySource != None:
            for (ind, element) in enumerate(self.elements):                
                element.propagateRays(self.raySource.rays)