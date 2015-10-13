'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np

import Raytracer.OpticalSurface as os
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])
                
class OpticalElement(object):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0]), n = 1.0, thickness = 1.0e-3, material = air):
        self.n = n
        self.thickness = thickness
        self.material = material
        
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

    def rotateElement(self, theta, phi):
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
            s.rotateExternal(theta, phi)
            
    def flipElement(self):
        self.surfaces.reverse()
        for surf in self.surfaces:
            surf.rotateExternal(np.pi, 0)  
        
    def setRotationMatrix(self, xpM):
        self.xn = np.dot(xpM, self.xn)
        self.xt = np.dot(xpM, self.xt)
        for s in self.surfaces:
            s.setRotationExternalMatrix(xpM)
        
    def initSurfaces(self):
        self.surfaces = [os.Surface(self.x, -self.xn, self.xt, self.n, material=self.material)]
        x2 = self.x
        x2[2] += self.thickness
        self.surfaces.append(os.Surface(x2, self.xn, self.xt))

    
    def propagateRays(self, rays):
        print type(self)
        for ray in rays:
            print ""
            print "=+=+=+=+ Ray  =+=+=+=+=+=+="
            for surf in self.surfaces:                
#                n0 = ray.n[-1]
#                (newX, newXp, newN) = surf.findIntersection(ray.x[-1], ray.xp[-1], n0)
                surf.findIntersectionRay(ray)
#                ray.addPos(newX, newXp, newN, 1.0)
                
class PrismElement(OpticalElement):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0]), n = 1.0, apexAngle = 60*np.pi/180, sideLength = 25e-3, material = air):
        self.apexAngle = apexAngle
        self.sideLength = sideLength
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, material = material)
            
    def initSurfaces(self):
        ap = oa.RectangularAperture([self.sideLength, self.sideLength])
        s1 = os.Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1.setRotationInternal(0, self.apexAngle/2)
        s1Pos = self.x+np.array([0,0,-np.sin(self.apexAngle/2)*self.sideLength/2,0])
        s1.setPosition(s1Pos)
        
        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2.setRotationInternal(0, -self.apexAngle/2)
        s2Pos = self.x+np.array([0,0,np.sin(self.apexAngle/2)*self.sideLength/2,0])
        s2.setPosition(s2Pos)
        self.surfaces = [s1, s2]
        
class PCXElement(OpticalElement):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0]), n = 1.0, r = 1.0, thickness = 5e-3, material = air, size = 12.7e-3):
        self.r1 = r
        self.size = size
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, thickness=thickness, material = material)        
        
    def initSurfaces(self):
        ap = oa.CircularAperture(self.size)
        s1 = os.SphericalSurface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, r=self.r1, material=self.material, aperture=ap)
        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2.setPosition(self.x+np.array([0,0,self.thickness,0]))
        self.surfaces = [s1, s2]
        
class ScreenElement(OpticalElement):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,-1,0]), xt = np.array([0,1,0,0]), material = air):
        super(ScreenElement, self).__init__(x=x, xn=xn, xt=xt, n=1.0, material = material)
        
    def initSurfaces(self):
        ap = oa.InifiniteAperture()
        self.surfaces = [os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)]
