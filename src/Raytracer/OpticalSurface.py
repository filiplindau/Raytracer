'''
Created on 11 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import OpticalMaterial as om
import OpticalAperture as oa
import Ray as r

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])

class Surface(object):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, -1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, material=air, aperture=oa.OpticalAperture(12.7e-3)):
        """ Basic surface. Implements functions for finding intersection with ray and coordinate transforms.
        
        Inputs:
        x: Origin position
        xn: Surface orientation (normal)
        xt: Surface orientation (tangent)
        n: Refractive index after passing surface (None if last surface in element)
        material: Material after passing surface
        aperture: size of the surface
        """
        self.x = x
        self.xn = xn / np.sqrt(np.dot(xn, xn))
        self.xt = xt / np.sqrt(np.dot(xt, xt))
        self.n = n
        self.material = material
        self.xpMint = np.identity(4)
        self.xpM = np.identity(4)
        self.generateTransformMatrix()
        self.aperture = aperture

        xt2 = np.hstack((np.cross(self.xt[0:3], self.xn[0:3]), 0))
        self.xpM = np.transpose(np.vstack((xt2, self.xt, self.xn, np.array([0, 0, 0, 1]))))        
        
        self.generateSurfaceEdge()
          
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
#        self.xpM = np.dot(self.xpMint, self.xpMext)
     
    def findIntersection(self, x, xp, n0=1.0):
        """ This is for intersection with a plane surface. Reimplement for new surface types. 
        Needs to find ray intersection and surface normal.
        """
        # Transform to local coordinate system:
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))        
        xpLocal = np.dot(self.xpM, xp)
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[2] / xpLocal[2]        
        xnLocal = np.array([0.0, 0.0, 1.0, 0.0])
        
        print "=============================="
        print "xnLocal find: ", xnLocal
        xNewLocal = xLocal + t * xpLocal
        xpNewLocal = self.calculateLocalRefraction(xNewLocal, xpLocal, xnLocal, self.n, n0)
        xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))
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
    
    def findIntersectionRay(self, ray, xMel, xMTel, xpMel):
        """ This is for intersection with a plane surface. Re-implement for new surface types. 
        Needs to find ray intersection and surface normal.
        """
        x = np.dot(xpMel, np.dot(xMel, ray.x[-1]))
        xp = np.dot(xpMel, ray.xp[-1])
        n0 = ray.n[-1]
        n = self.material.getRefractiveIndex(ray.l)
        ng = self.material.getGroupRefractiveIndex(ray.l)
        # Transform to local coordinate system:       
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))        
        xpLocal = np.dot(self.xpM, xp)
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[2] / xpLocal[2]        
        xnLocal = np.array([0.0, 0.0, 1.0, 0.0])
        
        xNewLocal = xLocal + t * xpLocal
        xpNewLocal = self.calculateLocalRefraction(xNewLocal, xpLocal, xnLocal, n, n0)
        xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))            
        xNewEl = np.dot(xMTel, np.dot(np.transpose(xpMel), xNew))
        xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
        xpNewEl = np.dot(np.transpose(xpMel), xpNew)
        
        if self.aperture.pointInAperture(xNewLocal) == True:
            ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)
        
    def findIntersectionRays(self, rays):
        """ This is for intersection with a plane surface. Re-implement for new surface types. 
        Needs to find ray intersection and surface normal.
        
        Using the new rays matrix
        """
        x = rays[:, 0, :]
        xp = rays[:, 1, :]
        data = rays[:, 2, :]
        n0 = rays[0, 2, 1] # Assume all rays have the same refractive index
        l = rays[0, 2, 0] # Assume all rays have the same wavelength
        n = self.material.getRefractiveIndex(l)
        ng = self.material.getGroupRefractiveIndex(l)
        data[:, 1] = n
        data[:, 2] = ng
        # Transform to local coordinate system:       
        xLocal = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(x))))        
        xpLocal = np.transpose(np.dot(self.xpM, np.transpose(xp)))
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[:, 2] / xpLocal[:, 2]        
        xnLocal = np.reshape(np.tile(np.array([0.0, 0.0, 1.0, 0.0]), t.shape[0]), (t.shape[0], 4))
        
        xNewLocal = xLocal + np.transpose(np.multiply(np.transpose(xpLocal), t))
        intersectInd = self.aperture.pointInAperture(xNewLocal)
        xpNewLocal = self.calculateLocalRefractions(xNewLocal[intersectInd, :], xpLocal[intersectInd, :], xnLocal[intersectInd, :], n, n0)
        xNew = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(xNewLocal[intersectInd, :]))))            
        xpNew = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(xpNewLocal)))
        
        print "intersectInd: ", intersectInd
        print "xLocal: ", xLocal[0, :]        
        print "xNewLocal: ", xNewLocal[0, :]
        print "xpLocal: ", xpLocal[0, :]
        print "xpNewLocal: ", xpNewLocal
        
        rays[intersectInd, 0, :] = xNew
        rays[intersectInd, 1, :] = xpNew
        rays[intersectInd, 2, 1] = n
        rays[intersectInd, 2, 2] = ng
        return rays
#        return np.dstack((xNew,xpNew,data)).swapaxes(1,2)
        
#         if self.aperture.pointInAperture(xNewLocal) == True:
#             ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)

    def calculateLocalRefractions(self, x, xp, xn, n, n0):
        """ Calculate refraction at surface for local coordinates x and
        local direction xp. Returns new direction and refractive index. 
        
        Inputs:
        x: Local coordinate for ray intersection
        xp: Local coordinates for ray direction
        xn: Local surface normal
        n: refractive index in this material
        n0: refractive index in previous material
        """
       
        n_r = n0 / n
        costh1 = np.sum(np.multiply(xn, xp), 1)
        st1 = xp - np.transpose(np.multiply(np.transpose(xn), costh1))
        cos2th2 = 1 - n_r ** 2 * (1 - costh1 ** 2)
        
        k2 = (1 - cos2th2) / (np.sum(np.multiply(st1, st1), 1) + 1e-10)
#        print "k2: ", k2
#        if k2 >= 0:
        xp2 = np.transpose(np.multiply(np.transpose(xn), np.sqrt(cos2th2) * np.sign(costh1))
                               + np.multiply(np.transpose(st1), np.sqrt(k2)))
#        else:
#            xp2 = st1 - costh1*xn
#         print "theta_i: ", np.arccos(costh1)*180/np.pi
#         print "theta_t: ", np.arccos(np.sqrt(cos2th2))*180/np.pi
#         print "n_r: ", n_r
#         print "st1: ", st1
#         print "cos(theta1)", costh1
#         print "cos(theta2)", np.sqrt(cos2th2)
        return xp2
                
    def calculateLocalRefraction(self, x, xp, xn, n, n0):
        """ Calculate refraction at surface for local coordinates x and
        local direction xp. Returns new direction and refractive index. 
        
        Inputs:
        x: Local coordinate for ray intersection
        xp: Local coordinates for ray direction
        xn: Local surface normal
        n: refractive index in this material
        n0: refractive index in previous material
        """
       
        n_r = n0 / n
        costh1 = np.dot(xn, xp)
        st1 = xp - costh1 * xn
        cos2th2 = 1 - n_r ** 2 * (1 - costh1 ** 2)
        
#         print "xNewLocal: ", x
#         print "costh1: ", costh1
#         print "cos2th2: ", cos2th2
        k2 = (1 - cos2th2) / (np.dot(st1, st1) + 1e-10)
#        print "k2: ", k2
        if k2 >= 0:
            xp2 = np.sqrt(cos2th2) * np.sign(costh1) * xn + np.sqrt(k2) * st1
        else:
            xp2 = st1 - costh1 * xn
#         print "theta_i: ", np.arccos(costh1)*180/np.pi
#         print "theta_t: ", np.arccos(np.sqrt(cos2th2))*180/np.pi
#         print "n_r: ", n_r
#         print "st1: ", st1
#         print "cos(theta1)", costh1
#         print "cos(theta2)", np.sqrt(cos2th2)
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
       
        n_r = n0 / self.n
        ndxp = -np.dot(xn, xp)
        sinth2 = n_r ** 2 * (1 - ndxp ** 2)
        print "-------------------"
        print "theta_i: ", np.arccos(ndxp) * 180 / np.pi
        print "theta_t: ", np.arcsin(np.sqrt(sinth2)) * 180 / np.pi
        print "ndxp: ", ndxp
        print "n_r: ", n_r
        print "sin(theta)**2: ", sinth2
        print "n_r*ndxp...", (n_r * ndxp - np.sqrt(1 - sinth2)) * xn
        xp_o = n_r * xp + (n_r * ndxp - np.sqrt(1 - sinth2)) * xn                    
        return xp_o
 

    def setPosition(self, newPos):
        if newPos.shape[0] == 3:
            self.x = np.hstack((newPos, 1))
        else:
            self.x = newPos
        self.generateTransformMatrix()
 
    def setRotationExternal(self, theta, phi): 
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.xpMext = np.dot(thM, phM)
#        self.xpM = np.dot(self.xpMint, self.xpMext)
        
    def rotateExternal(self, theta, phi):
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.xpMext = np.dot(np.dot(thM, phM), self.xpMext)
#        self.xpM = np.dot(self.xpMint, self.xpMext)
        
    def setRotationExternalMatrix(self, xpMext):          
        self.xpMext = xpMext
#        self.xpM = np.dot(self.xpMint, self.xpMext)
                
    def setRotationInternal(self, theta, phi):
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.xpMint = np.dot(self.xpMint, np.dot(thM, phM))
#        self.xpM = np.dot(self.xpMint, self.xpMext)
        self.xpM = np.dot(self.xpM, np.dot(thM, phM))
        
    def generateSurfaceEdge(self):
        self.surfaceEdge = self.aperture.getEdge()
    
    def getEdges(self):
        se = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(self.surfaceEdge))))
        return se
    
    def getRayFootprint(self, ray):
        xLocal = np.dot(self.xpM, np.dot(self.xM, ray.x))
        return xLocal
        
class SphericalSurface(Surface):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, -1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, r=1.0, material=air, aperture=oa.OpticalAperture(12.7e-3)):
        self.r = r
        Surface.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material, aperture=aperture)        
        
    def findIntersection(self, x, xp, n0=1.0):
        # Transform to local coordinate system:
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))
        xpLocal = np.dot(self.xpM, xp)
        # Spherical surface
        # Intersection where xLocal + t*xpLocal crosses |x-xc|^2 = r^2
        # 
        xDiff = xLocal[0:3] - np.array([0, 0, -self.r])
        
        b = 2 * np.dot(xDiff, xpLocal[0:3])
        c = np.dot(xDiff, xDiff) - self.r ** 2
        sq2 = b * b / 4 - c
        if sq2 > 0:
            # Intersection exists
            sq = np.sqrt(sq2)
            t1 = -b / 2 + sq
            t2 = -b / 2 - sq
            if xLocal[2] + t1 * xpLocal[2] + self.r > 0:
                t = t1
            else:
                t = -b / 2 - sq                
            xNewLocal = xLocal + t * xpLocal
            xnLocal = xNewLocal - np.array([0, 0, -self.r, 1])
            xnLocal = xnLocal / np.sqrt(np.dot(xnLocal, xnLocal))
            
            xpNewLocal = self.calculateLocalRefraction(xNewLocal, xpLocal, xnLocal, self.n, n0)
            xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))
            xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
            nNew = self.n
        else:
            # No intersection
            xNew = None
            xpNew = None
            nNew = n0
                
        return (xNew, xpNew, nNew)

    def findIntersectionRay(self, ray, xMel, xMTel, xpMel):
        x = np.dot(xpMel, np.dot(xMel, ray.x[-1]))
        xp = np.dot(xpMel, ray.xp[-1])
        n0 = ray.n[-1]
        n = self.material.getRefractiveIndex(ray.l)
        ng = self.material.getGroupRefractiveIndex(ray.l)
        # Transform to local coordinate system:       
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))        
        xpLocal = np.dot(self.xpM, xp)
        
        # Spherical surface
        # Intersection where xLocal + t*xpLocal crosses |x-xc|^2 = r^2
        # 
        xDiff = xLocal[0:3] - np.array([0, 0, -self.r])
        
        b = 2 * np.dot(xDiff, xpLocal[0:3])
        c = np.dot(xDiff, xDiff) - self.r ** 2
        sq2 = b * b / 4 - c
        if sq2 > 0:
            # Intersection exists
            sq = np.sqrt(sq2)
            t1 = -b / 2 + sq
            t2 = -b / 2 - sq
            if xLocal[2] + t1 * xpLocal[2] + self.r > 0:
                t = t1
            else:
                t = t2                
            xNewLocal = xLocal + t * xpLocal
            xnLocal = xNewLocal - np.array([0, 0, -self.r, 1])
            xnLocal = xnLocal / np.sqrt(np.dot(xnLocal, xnLocal))
            
            xpNewLocal = self.calculateLocalRefraction(xNewLocal, xpLocal, xnLocal, n, n0)
            xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))            
            xNewEl = np.dot(xMTel, np.dot(np.transpose(xpMel), xNew))
            xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
            xpNewEl = np.dot(np.transpose(xpMel), xpNew)
            if self.aperture.pointInAperture(xNewLocal) == True:
                ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)
        else:
            # No intersection
            xNew = None
            xpNew = None
            nNew = n0

    def findIntersectionRays(self, rays):
        """ This is for intersection with a plane surface. Re-implement for new surface types. 
        Needs to find ray intersection and surface normal.
        
        Using the new rays matrix
        """
        x = rays[:, 0, :]
        xp = rays[:, 1, :]
        data = rays[:, 2, :]
        n0 = rays[0, 2, 1] # Assume all rays have the same refractive index
        l = rays[0, 2, 0] # Assume all rays have the same wavelength
        n = self.material.getRefractiveIndex(l)
        ng = self.material.getGroupRefractiveIndex(l)
        data[:, 1] = n
        data[:, 2] = ng
        nRays = x.shape[0]
        # Transform to local coordinate system:       
        xLocal = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(x))))        
        xpLocal = np.transpose(np.dot(self.xpM, np.transpose(xp)))
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[:, 2] / xpLocal[:, 2]        
        xnLocal = np.reshape(np.tile(np.array([0.0, 0.0, 1.0, 0.0]), t.shape[0]), (t.shape[0], 4))

        # Spherical surface
        # Intersection where xLocal + t*xpLocal crosses |x-xc|^2 = r^2
        # 
        xc = np.reshape(np.tile(np.array([0, 0, -self.r, 1]), nRays), (nRays, 4))
        xDiff = xLocal[:, 0:3] - xc[:, 0:3]
        
        b = 2 * np.sum(np.multiply(xDiff, xpLocal[:, 0:3]), 1)
        c = np.sum(np.multiply(xDiff, xDiff), 1) - self.r ** 2
        sq2 = np.multiply(b, b) / 4 - c
#        print "sq2: ", sq2
        sq2PosInd = sq2 > 0 # Intersection exists for rays with indexes in sq2PosInd
        sq = np.sqrt(sq2[sq2PosInd])
        t = -b[sq2PosInd] / 2 + sq        
        tNegInd = xLocal[sq2PosInd, 2] + np.multiply(t, xpLocal[sq2PosInd, 2]) + self.r < 0
        t[tNegInd] -= 2 * sq[tNegInd]
        
        xNewLocal = xLocal[sq2PosInd, :] + np.transpose(np.multiply(np.transpose(xpLocal[sq2PosInd, :]), t))
        intersectInd = self.aperture.pointInAperture(xNewLocal)
        xnLocal = xNewLocal[intersectInd, :] - xc[sq2PosInd, :][intersectInd, :]
        xnLocal = np.transpose(np.divide(np.transpose(xnLocal) , np.sqrt(np.sum(np.multiply(xnLocal, xnLocal), 1))))
        xpNewLocal = self.calculateLocalRefractions(xNewLocal[intersectInd, :], xpLocal[sq2PosInd, :][intersectInd, :], xnLocal, n, n0)
        xNew = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(xNewLocal[intersectInd, :]))))            
        xpNew = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(xpNewLocal)))
        
#        print "sq2PosInd: ", sq2PosInd.sum()
#        print "tNegInd: ", tNegInd
#        print "t: ", t
#        print "xpLocal: ", xpLocal
#        print "xpNewLocal: ", xpNewLocal
#        print "xnLocal: ", xnLocal
#        print "xp: ", xp
        rays[sq2PosInd[intersectInd], 0, :] = xNew
        rays[sq2PosInd[intersectInd], 1, :] = xpNew
        rays[sq2PosInd[intersectInd], 2, 1] = n
        rays[sq2PosInd[intersectInd], 2, 2] = ng
        return rays
#        return np.dstack((xNew,xpNew,data)).swapaxes(1,2)
        
#         if self.aperture.pointInAperture(xNewLocal) == True:
#             ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)
                
    def generateSurfaceEdge(self):
        nbrPoints = 16
        theta0 = np.arcsin(self.aperture.size / self.r)
        theta = np.linspace(-theta0, theta0, nbrPoints)
        xe = np.vstack((-self.r * np.sin(theta), np.zeros(nbrPoints), -self.r * (1 - np.cos(theta)), np.ones(nbrPoints)))
        ye = np.vstack((np.zeros(nbrPoints), -self.r * np.sin(theta), -self.r * (1 - np.cos(theta)), np.ones(nbrPoints)))
        self.surfaceEdge = np.transpose(np.hstack((xe, ye)))
