'''
Created on 10 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import OpticalElement as oe

class OpticalSystem(object):
    def __init__(self):
        self.elements = []
        self.raySource = None
        
    def addElement(self, element):
        self.elements.append(element)
        
    def setRaySource(self, raySource):
        self.raySource = raySource
        
    def traceSystem(self):
        if self.raySource != None:
            for element in self.elements:
                element.propagateRays(self.raySource.rays)