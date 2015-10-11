'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np

from Raytracer.OpticalSurface import Surface
                
class OpticalElement(object):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0]), n = 1.0, thickness = 1.0e-3):
        self.n = n
        self.thickness = thickness
        
        if x.shape[0] == 3:
            self.x = np.hstack((x,1))
        else:
            self.x = x
        if xn.shape[0] == 3:
            self.xn = np.hstack((xn,0))
        else:
            self.xn = xn
        if xt.shape[0] == 3:
            self.xt = np.hstack((xt,0))
        else:
            self.xt = xt
        self.initSurfaces()
        
    def setPosition(self, newPos):
        if newPos.shape[0] == 3:
            self.x = np.hstack((newPos,1))
        else:
            self.x = newPos
        self.initSurfaces()
        
    def setRotation(self, theta, phi):
        thM = np.array([[1.0, 0.0,            0.0,           0.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 0.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 1.0], 
                         [0.0,          1.0, 0.0,         1.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 1.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        xpM = np.dot(thM, phM)
        self.xn = np.dot(xpM, self.xn)
        self.xt = np.dot(xpM, self.xt)
        
        for s in self.surfaces:
            s.setRotationExternal(theta, phi)
#        self.initSurfaces()
        
    def initSurfaces(self):
        self.surfaces = [Surface(self.x, -self.xn, self.xt, self.n)]
        x2 = self.x
        x2[2] += self.thickness
        self.surfaces.append(Surface(x2, self.xn, self.xt))

    
    def propagateRays(self, rays):
        print type(self)
        for ray in rays:
            print "=+=+=+=+ Ray  =+=+=+=+=+=+="
            for surf in self.surfaces:
                
                n0 = ray.n[-1]
                (newX, newXp, newN) = surf.findIntersection(ray.x[-1], ray.xp[-1], n0)
                ray.addPos(newX, newXp, newN, 1.0)
                
class PrismElement(OpticalElement):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0]), n = 1.0, apexAngle = 60*np.pi/180, sideLength = 25e-3):
        self.apexAngle = apexAngle
        self.sideLength = sideLength
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n)
            
    def initSurfaces(self):
        print "prism init surfaces"
        s1 = Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n)
        s1.setRotationInternal(0, self.apexAngle/2)
        s1Pos = self.x+np.array([0,0,-np.sin(self.apexAngle/2)*self.sideLength/2,0])
        s1.setPosition(s1Pos)
        s2 = Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0)
        s2.setRotationInternal(0, -self.apexAngle/2)
        s2Pos = self.x+np.array([0,0,np.sin(self.apexAngle/2)*self.sideLength/2,0])
        s2.setPosition(s2Pos)
        self.surfaces = [s1, s2]
        print s1Pos, s2Pos
        
