'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import numpy as np

import Raytracer.OpticalSurface as os
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import Raytracer.Ray as rr

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])
                
class OpticalElement(object):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, thickness=1.0e-3, material=air, size=12.7e-3):
        self.n = n
        self.thickness = thickness
        self.material = material
        self.size = size
        
        if x.shape[0] == 3:
            self.x = np.hstack((x, 1))
        else:
            self.x = x
        if xn.shape[0] == 3:
            self.xn = np.hstack((xn, 0))
        else:
            self.xn = xn
        if xt.shape[0] == 3:
            self.xt = np.hstack((xt, 0))
        else:
            self.xt = xt
        self.generateTransformMatrix()
        self.initSurfaces()        
        
    def setPosition(self, newPos):
        if newPos.shape[0] == 3:
            self.x = np.hstack((newPos, 1))
        else:
            self.x = newPos
        #self.initSurfaces()
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
        xt2 = np.hstack((np.cross(self.xt[0:3], self.xn[0:3]), 0))
        self.xpM = np.transpose(np.vstack((xt2, self.xt, self.xn, np.array([0, 0, 0, 1]))))     
        
    def setRotation(self, theta, phi):
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.xpM = np.dot(thM, phM)
        self.xn = np.dot(self.xpM, self.xn)
        self.xt = np.dot(self.xpM, self.xt)
        
        
#         for s in self.surfaces:
#             s.setRotationExternal(theta, phi)

    def rotateElement(self, theta, phi):
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.xpM = np.dot(self.xpM, np.dot(thM, phM))
        self.xn = np.dot(self.xpM, self.xn)
        self.xt = np.dot(self.xpM, self.xt)
        
#         for s in self.surfaces:
#             s.rotateExternal(theta, phi)
            
    def flipElement(self):
        nList = [surf.n for surf in self.surfaces]
        nList.reverse()
        matList = [surf.material for surf in self.surfaces]
        matList.reverse()
        for (ind, surf) in enumerate(self.surfaces):
            surf.n = nList[ind]
            surf.material = matList[ind]
        self.surfaces.reverse()
        self.rotateElement(np.pi, 0)
        
    def setRotationMatrix(self, xpM):
        self.xn = np.dot(xpM, self.xn)
        self.xt = np.dot(xpM, self.xt)
        for s in self.surfaces:
            s.setRotationExternalMatrix(xpM)
        
    def getEdges(self):
        edges = []
        for surf in self.surfaces:
            xe = surf.getEdges()
            edges.append(np.dot(self.xMT, np.dot(np.transpose(self.xpM), xe)))
        return edges
            
    def initSurfaces(self):
        ap = oa.CircularAperture(self.size)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=self.xn, xt=self.xt, n=self.n, material=air, aperture=ap)
        self.surfaces = [s1, s2]

    
    def propagateRays(self, rays):
        ''' Propagate rays using the new rays matrix
        The rays are transformed to element local coordinates before 
        being sent to the surfaces.
        '''
        print type(self)
        # raysEl contains the rays transformed to element local coordinates
        raysEl = rays.copy()
        # raysGl contains the returned coordinates from a surface transformed back to global coordinates
        raysGl = rays.copy()
        raysGlList = []
        raysEl[:, 0, :] = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(rays[:, 0, :]))))
        raysEl[:, 1, :] = np.transpose(np.dot(self.xpM, np.transpose(rays[:, 1, :])))
        print "raysEl in: ", raysEl[0, 1, :]
        print "raysGl in: ", raysGl[0, 1, :]
        for surf in self.surfaces:                            
            raysEl = surf.findIntersectionRays(raysEl)
            print "raysEl: ", raysEl[0, 1, :]
            raysGlNew = raysGl.copy()
            raysGlNew[:, 0, :] = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(raysEl[:, 0, :]))))
            raysGlNew[:, 1, :] = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(raysEl[:, 1, :]))) 
            raysGlList.append(raysGlNew)
            print "raysGl new: ", raysGlNew[0, 1, :]
        return raysGlList

    def propagateRaysOld(self, rays):
        ''' Propagate rays using the old rays structure (class with lists)
        '''
        print type(self)
        for ray in rays:
            print ""
            print "=+=+=+=+ Ray  =+=+=+=+=+=+="
            for surf in self.surfaces:                
                print "--------------------------------"
                surf.findIntersectionRay(ray, self.xM, self.xMT, self.xpM)
                
class PrismElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, apexAngle=60 * np.pi / 180, sideLength=25e-3, material=air):
        self.apexAngle = apexAngle
        self.sideLength = sideLength
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material)
            
    def initSurfaces(self):
        ap = oa.RectangularAperture([self.sideLength, self.sideLength])
#        s1 = os.Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn= -self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1.setRotationInternal(0, self.apexAngle / 2)
        s1Pos = np.array([0, 0, -np.sin(self.apexAngle / 2) * self.sideLength / 2, 1])
        s1.setPosition(s1Pos)
        
#        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, 0, 1]), xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2.setRotationInternal(0, -self.apexAngle / 2)
        s2Pos = np.array([0, 0, np.sin(self.apexAngle / 2) * self.sideLength / 2, 1])
        s2.setPosition(s2Pos)
        self.surfaces = [s1, s2]
        
class PCXElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, r=1.0, thickness=5e-3, material=air, size=12.7e-3):
        self.r1 = r
        self.size = size
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, thickness=thickness, material=material)        
        
    def initSurfaces(self):
        ap = oa.CircularAperture(self.size)
#        s1 = os.SphericalSurface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, r=self.r1, material=self.material, aperture=ap)
        s1 = os.SphericalSurface(x=np.array([0, 0, 0, 1]), xn= -self.xn, xt=self.xt, n=self.n, r=self.r1, material=self.material, aperture=ap)
#        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn= -self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
#        s2.setPosition(self.x+np.array([0,0,self.thickness,0]))
        self.surfaces = [s1, s2]
        
class ScreenElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), material=air):
        super(ScreenElement, self).__init__(x=x, xn=xn, xt=xt, n=1.0, material=material, thickness=0.0)
        
    def initSurfaces(self):
        ap = oa.InifiniteAperture()
        self.surfaces = [os.Surface(x=np.array([0, 0, 0, 1]), xn= -self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)]
