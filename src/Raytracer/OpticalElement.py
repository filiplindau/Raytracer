"""
Created on 9 Oct 2015

@author: Filip Lindau
"""

import numpy as np

import Raytracer.OpticalSurface as os
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import Raytracer.Ray as rr

import logging

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])
ml = om.MaterialLibrary()


class OpticalElement(object):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0,
                 thickness=1.0e-3, material=air, size=12.7e-3):
        """
        :param x: Position vector in global coordinates
        :param xn: Element normal vector
        :param xt: Element tangent vector
        :param n: Refractive index if no material is specified
        :param thickness: Element thickness in meters
        :param material: Material as a OpticalMaterial instance
        :param size: Element aperture size in meters
        :returns:
        """
        self.n = n
        self.thickness = thickness
        self.material = material
        self.size = size

        self.logger = logging.getLogger("Element.")
        self.logger.setLevel(logging.INFO)

        self.xM = None              # Transform matrix for position
        self.xMT = None             # Transposed position transform matrix
        self.xpM = None             # Transform matrix for angle
        self.surfaces = list()
        
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
        self.generate_transform_matrix()
        self.init_surfaces()
        
    def set_position(self, new_pos):
        if new_pos.shape[0] == 3:
            self.x = np.hstack((new_pos, 1))
        else:
            self.x = new_pos
        self.generate_transform_matrix()

    def generate_transform_matrix(self):
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
        
    def set_rotation(self, theta, phi):
        """
        Set element rotation to theta, phi angles

        :param theta: Rotation along theta axis
        :param phi: Rotation along phi axis
        """
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

    def rotate_element(self, theta, phi):
        """
        Rotate the element relative to the current rotation

        :param theta: Rotation along theta axis
        :param phi: Rotatino along phi axis
        """
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
            
    def flip_element(self):
        nList = [surf.n for surf in self.surfaces]
        nList.reverse()
        matList = [surf.material for surf in self.surfaces]
        matList.reverse()
        for (ind, surf) in enumerate(self.surfaces):
            surf.n = nList[ind]
            surf.material = matList[ind]
        self.surfaces.reverse()
        self.rotate_element(np.pi, 0)
        
    def reverse_element(self):
        nList = [surf.n for surf in self.surfaces]
        nList.reverse()
        matList = [surf.material for surf in self.surfaces]
        matList.reverse()
        for (ind, surf) in enumerate(self.surfaces):
            surf.n = nList[ind]
            surf.material = matList[ind]
        self.surfaces.reverse()

    def set_rotation_matrix(self, xpM):
        self.xn = np.dot(xpM, self.xn)
        self.xt = np.dot(xpM, self.xt)
        for s in self.surfaces:
            s.set_rotation_external_matrix(xpM)
        
    def get_edges(self):
        edges = []
        for surf in self.surfaces:
            xe = surf.get_edges()
            edges.append(np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(xe)))))
        return edges
            
    def init_surfaces(self):
        ap = oa.CircularAperture(self.size)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                        material=self.material, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=self.xn, xt=self.xt, n=self.n,
                        material=air, aperture=ap)
        self.surfaces = [s1, s2]

    def propagate_rays(self, rays):
        """ Propagate rays using the new rays matrix
        The rays are transformed to element local coordinates before
        being sent to the surfaces.
        """
        # raysEl contains the rays transformed to element local coordinates
        raysEl = rays.copy()
        # raysGl contains the returned coordinates from a surface transformed back to global coordinates
        raysGl = rays.copy()
        raysGlList = []
        raysEl[:, 0, :] = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(rays[:, 0, :]))))
        raysEl[:, 1, :] = np.transpose(np.dot(self.xpM, np.transpose(rays[:, 1, :])))
#        print "raysEl in: ", raysEl[0, 1, :]
#        print "raysGl in: ", raysGl[0, 1, :]
        for surf in self.surfaces:
            raysEl = surf.find_intersection_rays(raysEl)
#            print "raysEl: ", raysEl[0, 1, :]
            raysGlNew = raysEl.copy()

            raysGlNew[:, 0, :] = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM),
                                                                      np.transpose(raysEl[:, 0, :]))))
            raysGlNew[:, 1, :] = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(raysEl[:, 1, :])))
            raysGlList.append(raysGlNew)
#            print "raysGl new: ", raysGlNew[0, 1, :]
        return raysGlList
    
    def get_rays_footprint(self, rays, surface_number):
        raysT = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(rays))))
        return self.surfaces[surface_number].get_rays_footprint(raysT)

    def propagate_rays_old(self, rays):
        """ Propagate rays using the old rays structure (class with lists)
        """
        for ray in rays:
            self.logger.debug("\n=+=+=+=+ Ray  =+=+=+=+=+=+=")
            for surf in self.surfaces:                
                self.logger.debug("--------------------------------")
                surf.find_intersection_ray(ray, self.xM, self.xMT, self.xpM)


class PrismElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0,
                 apex_angle=60 * np.pi / 180, side_length=25e-3, material=air):
        self.apex_angle = apex_angle
        self.side_length = side_length
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material)
            
    def init_surfaces(self):
        ap = oa.RectangularAperture([self.side_length, self.side_length])
#        s1 = os.Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                        material=self.material, aperture=ap)
        s1.set_rotation_internal(0, self.apex_angle / 2)
        s1_pos = np.array([0, 0, -np.sin(self.apex_angle / 2) * self.side_length / 2, 1])
        s1.set_position(s1_pos)
        
#        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, 0, 1]), xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2.set_rotation_internal(0, -self.apex_angle / 2)
        s2_pos = np.array([0, 0, np.sin(self.apex_angle / 2) * self.side_length / 2, 1])
        s2.set_position(s2_pos)
        self.surfaces = [s1, s2]


class PCXElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, r=1.0,
                 thickness=5e-3, material=air, size=12.7e-3):
        self.r1 = r
        self.size = size
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, thickness=thickness, material=material)        
        
    def init_surfaces(self):
        ap = oa.CircularAperture(self.size)
#        s1 = os.SphericalSurface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, r=self.r1,
        #        material=self.material, aperture=ap)
        s1 = os.SphericalSurface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n, r=self.r1,
                                 material=self.material, aperture=ap)
#        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=-self.xn, xt=self.xt, n=1.0,
                        material=air, aperture=ap)
#        s2.setPosition(self.x+np.array([0,0,self.thickness,0]))
        self.surfaces = [s1, s2]


class ScreenElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), material=air):
        super(ScreenElement, self).__init__(x=x, xn=xn, xt=xt, n=1.0, material=material, thickness=0.0)
        
    def init_surfaces(self):
        ap = oa.InifiniteAperture()
        self.surfaces = [os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                                    material=self.material, aperture=ap)]


class GratingElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0,
                 thickness=2e-3, grating_period=1.0/1800e3, m=1, side_length=25e-3, material=ml.get_material("fs")):
        self.grating_period = grating_period
        self.side_length = side_length
        self.thickness = thickness
        self.m = m
        OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material)

    def init_surfaces(self):
        self.logger.info("Init grating surfaces")
        ap = oa.RectangularAperture([self.side_length, self.side_length])
        #        s1 = os.Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                        material=self.material, aperture=ap)
        s1.set_rotation_internal(0, 0)

        #        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=self.xn, xt=self.xt, n=1.0, material=air,
                        aperture=ap, grating_period=self.grating_period, m=self.m)
        s2.set_rotation_internal(0, 0)
        self.surfaces = [s1, s2]


class MirrorElement(OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]),
                 n=1.0, material=air, thickness=5e-5, size=25.4e-3):
        self.size = size
        self.thickness = thickness
        material.reflector = True
        super(MirrorElement, self).__init__(x=x, xn=xn, xt=xt, n=1.0, material=material, thickness=0.0)

    def init_surfaces(self):
        ap = oa.CircularAperture()
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                        material=self.material, aperture=ap)
        s1.set_rotation_internal(0, 0)

        #        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=self.xn, xt=self.xt, n=1.0,
                        material=air, aperture=ap, )
        s2.set_rotation_internal(0, 0)
        self.surfaces = [s1]



