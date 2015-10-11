'''
Created on 11 Oct 2015

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
        self.xpMint = np.identity(4)
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
        self.xpMext = np.transpose(np.vstack((xt2, self.xt, self.xn, np.array([0,0,0,1]))))        
        self.xpM = np.dot(self.xpMint, self.xpMext)
     
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
       
        n_r = n0/self.n
        costh1 = np.dot(xn, xp)
        st1 = xp-costh1*xn
        cos2th2 = 1-n_r**2*(1-costh1**2)
        k2 = (1-cos2th2)/np.dot(st1,st1)
        if k2 > 0:
            xp2 = np.sqrt(cos2th2)*np.sign(costh1)*xn+np.sqrt(k2)*st1
        else:
            xp2 = st1 + 2*costh1*xn
        print "theta_i: ", np.arccos(costh1)*180/np.pi
        print "theta_t: ", np.arccos(np.sqrt(cos2th2))*180/np.pi
        print "n_r: ", n_r
        print "st1: ", st1
        print "cos(theta1)", costh1
        print "cos(theta2)", np.sqrt(cos2th2)
        return xp2

    def calculateLocalRefractionOld(self, x, xp, xn, n0):
        """ Calculate refraction at surface for local coordinates x and
        local direction xp. Returns new direction and refractive index. 
        
        Inputs:
        x: Local coordinate for ray intersection
        xp: Local coordinates for ray direction
        xn: Local surface normal
        n: refractive index in previous material
        """
       
        n_r = n0/self.n
        ndxp = -np.dot(xn, xp)
        sinth2 = n_r**2 * (1-ndxp**2)
        print "theta_i: ", np.arccos(ndxp)*180/np.pi
        print "theta_t: ", np.arcsin(np.sqrt(sinth2))*180/np.pi
        print "ndxp: ", ndxp
        print "n_r: ", n_r
        print "sin(theta)**2: ", sinth2
        print "n_r*ndxp...", (n_r*ndxp-np.sqrt(1-sinth2))*xn
        xp_o = n_r*xp+(n_r*ndxp-np.sqrt(1-sinth2))*xn                    
        return xp_o
 

    def setPosition(self, newPos):
        if newPos.shape[0] == 3:
            self.x = np.hstack((newPos,1))
        else:
            self.x = newPos
        self.generateTransformMatrix()
 
    def setRotationExternal(self, theta, phi): 
        thM = np.array([[1.0, 0.0,            0.0,           0.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 0.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0], 
                         [0.0,          1.0, 0.0,         0.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        self.xpMext = np.dot(thM, phM)
        self.xpM = np.dot(self.xpMint, self.xpMext)
        
    def setRotationInternal(self, theta, phi):
        thM = np.array([[1.0, 0.0,            0.0,           0.0], 
                         [0.0, +np.cos(theta), np.sin(theta), 0.0], 
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0,            0.0,           1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0], 
                         [0.0,          1.0, 0.0,         0.0], 
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0,          0.0, 0.0,         1.0]])
         
        self.xpMint = np.dot(thM, phM)
        self.xpM = np.dot(self.xpMint, self.xpMext)
        
class SphericalSurface(Surface):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,-1,0]), xt = np.array([0,1,0,0]), n = 1.0, r = 1.0):
        Surface.__init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,-1,0]), xt = np.array([0,1,0,0]), n = 1.0)
        