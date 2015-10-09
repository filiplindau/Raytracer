'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np

class Surface(object):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,-1,0]), xt = np.array([0,1,0,0]), n = 1.0):
        """ Basic surface. Implements functions for finding intersection with ray and coordinate transforms.
        
        Inputs:
        x: Origin position
        xn: Surface orientation (normal)
        xt: Surface orientation (tangent)
        n: Refractive index after passing surface (None if last surface in element)
        """
        self.x = x
        self.xn = xn / np.sqrt(np.dot(xn,xn))
        self.xt = xt / np.sqrt(np.dot(xt,xt))
        self.n = n
        
        self.generateTransformMatrix()
          
    def generateTransformMatrix(self):
        self.xM = np.array([[1.0, 0.0, 0.0, -self.x[0]], 
                             [0.0, 1.0, 0.0, -self.x[1]],
                             [0.0, 0.0, 1.0, -self.x[2]],
                             [0.0, 0.0, 0.0, 1.0]])

        self.xMT = np.array([[1.0, 0.0, 0.0, self.x[0]], 
                             [0.0, 1.0, 0.0, self.x[1]],
                             [0.0, 0.0, 1.0, self.x[2]],
                             [0.0, 0.0, 0.0, 1.0]])
                
        # Third coordinate axis by cross product:
        xt2 = np.hstack((np.cross(self.xt[0:3], self.xn[0:3]),0))
        self.xpM = np.transpose(np.vstack((xt2, self.xt, self.xn, np.array([0,0,0,1]))))
     
    def findIntersection(self, x, xp, n0 = 1.0):
        """ Reimplement for new surface types. Needs to find ray intersection 
        and surface normal.
        """
        # Transform to local coordinate system:
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))
        xpLocal = np.dot(self.xpM, xp)
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[2]/xpLocal[2]        
        xnLocal = np.array([0.0, 0.0, 1.0, 0.0])
        
        print "=============================="
        print "xnLocal find: ", xnLocal
        xNewLocal = xLocal + t*xpLocal
        xpNewLocal = self.calculateLocalRefraction(xNewLocal, xpLocal, xnLocal, n0)
        xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM),xNewLocal))
        xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)

        print "xp: ", xp
        print "xpLocal: ", xpLocal
        print "xpNewLocal: ", xpNewLocal
        print "xpNew: ", xpNew
        print "------------------------------"        
        print "t: ", t
        print "x: ", x
        print "xLocal: ", xLocal
        print "xNewLocal: ", xNewLocal 
        print "xNew: ", xNew  
        return (xNew, xpNew, self.n)
 
    def calculateLocalRefraction(self, x, xp, xn, n0):
        """ Calculate refraction at surface for local coordinates x and
        local direction xp. Returns new direction and refractive index. 
        
        Inputs:
        x: Local coordinate for ray intersection
        xp: Local coordinates for ray direction
        xn: Local surface normal
        n: refractive index in previous material
        """
        # Surface normal at intersection:        
        xnLocal = xn[0:3]
#         print "xp:", xp
#         print "xn:", xn
#         print "xnLocal:", xnLocal
#         print "x: ", x
        
        n_r = n0/self.n
        nxxp = -np.cross(xnLocal, xp[0:3])
        print "n x xp: ", nxxp
        print "xnLocal: ", xnLocal
        xp_o = n_r*np.cross(xnLocal, nxxp)-xnLocal*np.sqrt(1-n_r**2*np.dot(nxxp, nxxp))                    
        print "n x (-n x xp): ", np.cross(xnLocal, nxxp)
        print "n*sqrt(...): ", xnLocal*np.sqrt(1-n_r**2*np.dot(nxxp, nxxp))
        print "xpNewLocal: ", xp_o
        return np.hstack((xp_o,0))
 
    def calculateGlobalRefraction(self, x, xp, xn, n0):
        """ Calculate refraction at surface for global coordinates x and
        global direction xp. Returns new direction and refractive index. 
        
        Inputs:
        x: Global coordinate for ray intersection
        xp: Global coordinates for ray direction
        n_r is the ratio of
        old and new refractive index
        """
        # Surface normal at intersection:        
        xnLocal = self.xn[0:3]
#         print "xp:", xp
#         print "xn:", xn
#         print "xnLocal:", xnLocal
#         print "x: ", x
        
        n_r = n0/self.n
        nxxp = np.cross(xnLocal, xp[0:3])
        xp_o = n_r*np.cross(xnLocal, -nxxp)-xnLocal*np.sqrt(1-n_r**2*np.dot(nxxp, nxxp))                    
        
        return (np.hstack((xp_o,0)), self.n)

    def setPosition(self, newPos):
        self.x = np.hstack((newPos, 1))
        self.generateTransformMatrix()
 
    def setRotation(self, theta, phi): 
        thM = np.array([[1.0, 0.0,            0.0,           1.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 1.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 1.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 1.0], 
                         [0.0,          1.0, 0.0,         1.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 1.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        self.xpM = np.dot(thM, phM)
                
class OpticalElement(object):
    def __init__(self, n = 1.0, thickness = 1.0e-3, x = np.array([0,0,0,1]), xn = np.array([0,0,1,0]), xt = np.array([0,1,0,0])):
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
        thM = np.array([[1.0, 0.0,            0.0,           1.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 1.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 1.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 1.0], 
                         [0.0,          1.0, 0.0,         1.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 1.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        xpM = np.dot(thM, phM)
        self.xn = np.dot(xpM, self.xn)
        self.xt = np.dot(xpM, self.xt)
        
        self.initSurfaces()
        
    def initSurfaces(self):
        self.surfaces = [Surface(self.x, -self.xn, self.xt, self.n)]
        x2 = self.x
        x2[2] += self.thickness
        self.surfaces.append(Surface(x2, self.xn, self.xt))

    
    def propagateRays(self, rays):
        for ray in rays:
            print "=+=+=+=+ Ray  =+=+=+=+=+=+="
            for surf in self.surfaces:
                
                n0 = ray.n[-1]
                (newX, newXp, newN) = surf.findIntersection(ray.x[-1], ray.xp[-1], n0)
                ray.addPos(newX, newXp, newN, 1.0)
            
