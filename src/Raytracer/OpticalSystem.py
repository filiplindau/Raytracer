'''
Created on 10 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import OpticalElement as oe
import OpticalMaterial as om

class OpticalSystem(object):
    ''' Implements an optical axis on which elements can be placed. The 
    axis can be rotated after each element. The element coordinate system
    is with respect to this axis, were then the z coordinate is distance
    along the axis.
    '''
    def __init__(self):
        self.elements = []
        self.raySourceList = []
        self.materialLibrary = om.MaterialLibrary()
        self.opticalAxisXpM = [np.identity(4)] 
        self.opticalAxisXM = [np.identity(4)]
        self.opticalAxisXMT = [np.identity(4)]
        self.opticalAxisTheta = 0.0
        self.opticalAxisPhi = 0.0
        
    def addElement(self, element):
#        element.rotateElement(self.opticalAxisTheta, self.opticalAxisPhi)
        self.elements.append(element)
        x = element.x
        self.opticalAxisXpM.append(self.opticalAxisXpM[-1].copy())
#        self.opticalAxisXM.append(self.opticalAxisXM[-1].copy())
#        self.opticalAxisXMT.append(self.opticalAxisXMT[-1].copy())
        
        xg  = np.dot(self.opticalAxisXMT[-1], np.dot(np.transpose(self.opticalAxisXpM[-1]), np.array([0,0,x[2],1])))
#        xg  = np.dot(np.identity(4), np.dot(self.opticalAxisXMT[-1], np.array([0,0,x[2],1])))
#        self.opticalAxisXM.append(np.identity(4))
#        self.opticalAxisXMT.append(np.identity(4))
        newXM = np.array([[1.0, 0.0, 0.0, -xg[0]],
                             [0.0, 1.0, 0.0, -xg[1]],
                             [0.0, 0.0, 1.0, -xg[2]],
                             [0.0, 0.0, 0.0, 1.0]])

        newXMT = np.array([[1.0, 0.0, 0.0, xg[0]],
                             [0.0, 1.0, 0.0, xg[1]],
                             [0.0, 0.0, 1.0, xg[2]],
                             [0.0, 0.0, 0.0, 1.0]])
        
        self.opticalAxisXM.append(newXM)
        self.opticalAxisXMT.append(newXMT)
        
    def rotateOpticalAxisAfterElement(self, theta, phi, elementNumber):
        self.opticalAxisTheta = theta
        self.opticalAxisPhi = phi
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.opticalAxisXpM[elementNumber + 1] = np.dot(self.opticalAxisXpM[elementNumber + 1], np.dot(thM, phM))
            

    def addRaySource(self, raySource):
        self.raySourceList.append(raySource)
        
    def traceSystem(self):
        if self.raySourceList != []:
            for raySource in self.raySourceList:
                for (ind, element) in enumerate(self.elements):    
                    raysT = raySource.rays.copy()
                    raysT[:, 0, :] = np.transpose(np.dot(self.opticalAxisXpM[ind], np.dot(self.opticalAxisXM[ind], np.transpose(raysT[:, 0, :]))))
                    raysT[:, 1, :] = np.transpose(np.dot(self.opticalAxisXpM[ind], np.transpose(raysT[:, 1, :])))        
                    raysList = element.propagateRays(raysT)
                    for rays in raysList:
                        raysT = rays.copy()
                        raysT[:, 0, :] = np.transpose(np.dot(self.opticalAxisXMT[ind], np.dot(np.transpose(self.opticalAxisXpM[ind]), np.transpose(raysT[:, 0, :]))))
                        raysT[:, 1, :] = np.transpose(np.dot(np.transpose(self.opticalAxisXpM[ind]), np.transpose(raysT[:, 1, :])))
                        raySource.updateRays(raysT) 

    def getRaysFootprint(self, raySourceNumber, elementNumber, surfaceNumber):
        sn = 1
        for el in range(elementNumber):
            sn += self.elements[el].surfaces.__len__()
        if self.raySourceList != []:
            rays = self.raySourceList[raySourceNumber].getRaysPosArray()  
            raysT = np.transpose(np.dot(self.opticalAxisXpM[elementNumber], np.dot(self.opticalAxisXM[elementNumber], np.transpose(rays[:, sn, :]))))
            return self.elements[elementNumber].getRaysFootprint(raysT, surfaceNumber)
            
        
    def getElementEdges(self, elementNumber):
        edges = self.elements[elementNumber].getEdges()
        edgesNew = []
        for edge in edges:
            edgesNew.append(np.transpose(np.dot(self.opticalAxisXMT[elementNumber],
                                       np.dot(np.transpose(self.opticalAxisXpM[elementNumber]), np.transpose(edge)))))
        return edgesNew

