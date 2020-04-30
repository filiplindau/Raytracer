'''
Created on 11 Oct 2015

@author: Filip Lindau
'''

import numpy as np
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import Raytracer.Ray as r

import logging

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])


class Surface(object):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, -1, 0]), xt=np.array([0, 1, 0, 0]),
                 n=1.0, material=air, aperture=oa.OpticalAperture(12.7e-3), grating_period=None, m=0):
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
        self.xM = np.identity(4)
        self.xMT = np.identity(4)
        self.xpMext = np.identity(4)
        self.generate_transform_matrix()
        self.aperture = aperture
        self.surface_edge = None
        self.grating_period = grating_period
        self.m = m                              # Grating diffraction order

        self.logger = logging.getLogger("Surface.")
        self.logger.setLevel(logging.DEBUG)

        xt2 = np.hstack((np.cross(self.xt[0:3], self.xn[0:3]), 0))
        self.xpM = np.transpose(np.vstack((xt2, self.xt, self.xn, np.array([0, 0, 0, 1]))))        
        
        self.generate_surface_edge()
          
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
#        self.xpM = np.dot(self.xpMint, self.xpMext)
     
    def find_intersection(self, x, xp, n0=1.0):
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
        
        self.logger.debug("==============================")
        self.logger.debug("xnLocal find: {0}".format(xnLocal))
        xNewLocal = xLocal + t * xpLocal
        xpNewLocal = self.calculate_local_refraction(xNewLocal, xpLocal, xnLocal, self.n, n0)
        xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))
        xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)

        self.logger.debug("xp: {0}".format(xp))
        self.logger.debug("xpLocal: {0}".format(xpLocal))
        self.logger.debug("xpNewLocal: {0}".format(xpNewLocal))
        self.logger.debug("xpNew: {0}".format(xpNew))
        self.logger.debug("------------------------------")
        self.logger.debug("t: {0}".format(t))
        self.logger.debug("x: {0}".format(x))
        self.logger.debug("xLocal: {0}".format(xLocal))
        self.logger.debug("xNewLocal: {0}".format(xNewLocal))
        self.logger.debug("xNew: {0}".format(xNew))
        return xNew, xpNew, self.n
    
    def find_intersection_ray(self, ray, xM_el, xMT_el, xpM_el):
        """ This is for intersection with a plane surface. Re-implement for new surface types. 
        Needs to find ray intersection and surface normal.
        """
        x = np.dot(xpM_el, np.dot(xM_el, ray.x[-1]))
        xp = np.dot(xpM_el, ray.xp[-1])
        n0 = ray.n[-1]
        n = self.material.get_refractive_index(ray.l)
        ng = self.material.get_group_refractive_index(ray.l)
        # Transform to local coordinate system:       
        xLocal = np.dot(self.xpM, np.dot(self.xM, x))        
        xpLocal = np.dot(self.xpM, xp)
        
        # Plane surface
        # Intersection where xLocal + t*xpLocal crosses z = 0
        # Reimplement here
        t = -xLocal[2] / xpLocal[2]        
        xnLocal = np.array([0.0, 0.0, 1.0, 0.0])
        
        xNewLocal = xLocal + t * xpLocal
        xpNewLocal = self.calculate_local_refraction(xNewLocal, xpLocal, xnLocal, n, n0)
        xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))            
        xNewEl = np.dot(xMT_el, np.dot(np.transpose(xpM_el), xNew))
        xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
        xpNewEl = np.dot(np.transpose(xpM_el), xpNew)
        
        if self.aperture.point_in_aperture(xNewLocal) is True:
            ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)
        
    def find_intersection_rays(self, rays):
        """ This is for intersection with a plane surface. Re-implement for new surface types. 
        Needs to find ray intersection and surface normal.
        
        Using the new rays matrix
        """
        x = rays[:, 0, :]
        xp = rays[:, 1, :]
        data = rays[:, 2, :]
        n0 = rays[0, 2, 1]      # Assume all rays have the same refractive index
        l = rays[0, 2, 0]       # Assume all rays have the same wavelength
        n = self.material.get_refractive_index(l)
        ng = self.material.get_group_refractive_index(l)
        reflector = self.material.get_reflector()

        self.logger.debug("\nIncoming rays in element local coords:\n\nx:\n {0}\n\nxp:\n {1}\n\n"
                          "".format(x, xp))

        # self.logger.info("n0: {0}, n: {1}".format(n0, n))
        data[:, 1] = n
        data[:, 2] = ng
        # Transform to local coordinate system:       
        x_local = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(x))))
        xp_local = np.transpose(np.dot(self.xpM, np.transpose(xp)))
        
        t, xn_local = self.surface_intersect(x_local, xp_local)

        x_new_local = x_local + np.transpose(np.multiply(np.transpose(xp_local), t))
        # Reset positions where the rays where pointing away from the surface:
        dir_wrong_ind = t < 0
        self.logger.debug("dir_wrong_ind: {0}".format(dir_wrong_ind))
        x_new_local[dir_wrong_ind, :] = x_local[dir_wrong_ind, :]
        intersect_ind = self.aperture.point_in_aperture(x_new_local)
        good_ind = np.logical_and(intersect_ind, np.logical_not(dir_wrong_ind))
        self.logger.debug("dir_wrong_ind: {0}, intersect_ind: {1}, x shape: {2}".format(dir_wrong_ind, intersect_ind,
                                                                                        x_local.shape))

        xp_new_local = xp_local.copy()
        xp_new_local[good_ind, :] = self.calculate_local_diffractions(
            x_new_local[good_ind, :], xp_local[good_ind, :],
            xn_local[good_ind, :], n0, n, reflector,
            np.array([1.0, 0, 0, 0]), l, self.grating_period, self.m)
        x_new = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(x_new_local))))
        xp_new = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(xp_new_local)))

        self.logger.debug("\n\nx_local:\n {0}\n\nxp_local:\n {1}\n\nt:\n {2}\n\n"
                          "xn_local:\n {3}\n\nx_new_local:\n {4}\n\nxp_new_local:\n {5}\n\n"
                          "x_new:\n {6}\n\nxp_new:\n {7}"
                          "".format(x_local, xp_local, t, xn_local, x_new_local, xp_new_local, x_new, xp_new))

        rays[:, 0, :] = x_new
        rays[:, 1, :] = xp_new
        rays[intersect_ind, 2, 1] = n
        rays[intersect_ind, 2, 2] = ng
        return rays

    def surface_intersect(self, x_local, xp_local):
        """
        Calculate intersection points and surface normals for rays intersecting plane surface.

        Reimplement for other surface geometries.

        :param x_local: Rays positions in surface local coordinates
        :type x_local: np.array
        :param xp_local: Rays directions in surface local coordinates
        :type xp_local: np.array
        :return: t, xn_local: distances to intersections, normal directions at intersection
        :rtype:
        """
        t = -x_local[:, 2] / xp_local[:, 2]
        xn_local = np.reshape(np.tile(np.array([0.0, 0.0, 1.0, 0.0]), t.shape[0]), (t.shape[0], 4))
        return t, xn_local

    def calculate_local_refractions(self, x, xp, xn, n, n0):
        """ Calculate refraction at surface for local coordinates x and
        local direction xp. Returns new direction.
        
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

    def calculate_local_diffractions(self, x, xp, xn, n1, n2, reflector=False,
                                     q=np.array([1.0, 0.0, 0.0, 0.0]), l0=263e-9, d=None, m=0):
        """ Calculate grating diffraction at surface for local coordinates x and
        local direction xp. Returns new direction and refractive index.

        From: Klein, James "Demystifying the sequential ray-tracing alogrithm", SPIE 2002

        Inputs:
        x: Local coordinate for ray intersection
        xp: Local coordinates for ray direction
        xn: Local surface normal
        n1: Refractive index in the medium where the ray is coming from
        n2: Refractive index in the medium where the ray is exiting to
        reflector: True if the surface is a reflector
        q: Grating unit vector in local coordinates
        l0: Ray wavelength in m
        d: surface grating period in m
        m: diffraction order to calculate
        """
        if reflector:
            if d is None:
                xp_eff = n1 / n2 * xp
            else:
                xp_eff = (n1 * xp + m * l0 / d * q) / n2
            xp2 = xp_eff - 2 * np.sum(np.multiply(xp_eff, xn), 1)[:, np.newaxis] * xn
            self.logger.info("\n-- Reflector --\nxp_eff:\n {0}\n\nxn:\n {1}\n\n"
                             "xp2:\n {2}\n".format(xp_eff, xn, xp2))
        else:
            if d is None:
                xp_eff = n1 / n2 * xp
            else:
                xp_eff = (n1 * xp + m * l0 / d * q) / n2
            xp_dot_n = np.sum(np.multiply(xp_eff, xn), 1)
            quad = xp_dot_n**2 - np.sum(np.multiply(xp_eff, xp_eff), 1) + 1
            if n1 > n2:
                gamma = - xp_dot_n + np.sqrt(quad)
                gamma2 = - xp_dot_n - np.sqrt(quad)
            else:
                gamma = - xp_dot_n - np.sqrt(quad)
                gamma2 = - xp_dot_n + np.sqrt(quad)
            # else:
            #     gamma = (1 - mu**2) / (mu * xp_dot_n - np.sqrt(quad))
            #     gamma2 = (1 - mu**2) / (mu * xp_dot_n + np.sqrt(quad))
            self.logger.debug("\nquad:    {0}\ngamma: {1}".format(quad, gamma))
            self.logger.debug("\nxp_eff: {0} \nxp_dot_n:    {1}\nxn: {2}".format(xp_eff, xp_dot_n, xn))
            self.logger.debug("\nxdotn2: {0} \nx_eff2:    {1}\nn2: {2}".format(xp_dot_n**2, np.sum(np.multiply(xp_eff, xp_eff), 1), n2))
            xp2 = xp_eff / n2 + gamma[:, np.newaxis] * xn
            xp2 = xp_eff + gamma[:, np.newaxis] * xn
            xp2[np.where(quad < 0)] = xp[np.where(quad < 0)]

        xp_norm = xp2 / np.sqrt((xp2**2).sum(1))[:, np.newaxis]
        # xp_norm = xp2
        self.logger.info("\nxp:\n {0}\n\nxp_norm:\n {1}\n\nxn: {2}\n\nn1: {3}\n\nn2: {4}\n".format(xp, xp_norm, xn, n1, n2))
        return xp_norm
                
    def calculate_local_refraction(self, x, xp, xn, n, n0):
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

    def calculate_local_refraction_old(self, x, xp, xn, n0):
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
        self.logger.debug("-------------------")
        self.logger.debug("theta_i: {0}".format(np.arccos(ndxp) * 180 / np.pi))
        self.logger.debug("theta_t: {0}".format(np.arcsin(np.sqrt(sinth2)) * 180 / np.pi))
        self.logger.debug("ndxp: {0}".format(ndxp))
        self.logger.debug("n_r: {0}".format(n_r))
        self.logger.debug("sin(theta)**2: {0}".format(sinth2))
        self.logger.debug("n_r*ndxp... {0}".format((n_r * ndxp - np.sqrt(1 - sinth2)) * xn))
        xp_o = n_r * xp + (n_r * ndxp - np.sqrt(1 - sinth2)) * xn                    
        return xp_o

    def set_position(self, newPos):
        if newPos.shape[0] == 3:
            self.x = np.hstack((newPos, 1))
        else:
            self.x = newPos
        self.generate_transform_matrix()
 
    def set_rotation_external(self, theta, phi):
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
        
    def rotate_external(self, theta, phi):
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
        
    def set_rotation_external_matrix(self, xpMext):
        self.xpMext = xpMext
#        self.xpM = np.dot(self.xpMint, self.xpMext)
                
    def set_rotation_internal(self, theta, phi):
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
        
    def generate_surface_edge(self):
        self.surface_edge = self.aperture.get_edge()
    
    def get_edges(self):
        se = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(self.surface_edge))))
        return se
    
    def get_rays_footprint(self, rays):
        xLocal = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(rays))))
        return xLocal


class SphericalSurface(Surface):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, -1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0, r=1.0, material=air, aperture=oa.OpticalAperture(12.7e-3)):
        self.r = r
        Surface.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material, aperture=aperture)

    def find_intersection(self, x, xp, n0=1.0):
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
            
            xpNewLocal = self.calculate_local_refraction(xNewLocal, xpLocal, xnLocal, self.n, n0)
            xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))
            xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
            nNew = self.n
        else:
            # No intersection
            xNew = None
            xpNew = None
            nNew = n0
                
        return xNew, xpNew, nNew

    def find_intersection_ray(self, ray, xM_el, xMT_el, xpM_el):
        x = np.dot(xpM_el, np.dot(xM_el, ray.x[-1]))
        xp = np.dot(xpM_el, ray.xp[-1])
        n0 = ray.n[-1]
        n = self.material.get_refractive_index(ray.l)
        ng = self.material.get_group_refractive_index(ray.l)
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
            
            xpNewLocal = self.calculate_local_refraction(xNewLocal, xpLocal, xnLocal, n, n0)
            xNew = np.dot(self.xMT, np.dot(np.transpose(self.xpM), xNewLocal))            
            xNewEl = np.dot(xMT_el, np.dot(np.transpose(xpM_el), xNew))
            xpNew = np.dot(np.transpose(self.xpM), xpNewLocal)
            xpNewEl = np.dot(np.transpose(xpM_el), xpNew)
            if self.aperture.point_in_aperture(xNewLocal) == True:
                ray.addPos(xNewEl, xpNewEl, n, ng, 1.0, t)
        else:
            # No intersection
            xNew = None
            xpNew = None
            nNew = n0

#     def find_intersection_rays(self, rays):
#         """ This is for intersection with a plane surface. Re-implement for new surface types.
#         Needs to find ray intersection and surface normal.
#
#         Using the new rays matrix
#         """
#         x = rays[:, 0, :]
#         xp = rays[:, 1, :]
#         data = rays[:, 2, :]
#         n0 = rays[0, 2, 1] # Assume all rays have the same refractive index
#         l = rays[0, 2, 0] # Assume all rays have the same wavelength
#         n = self.material.get_refractive_index(l)
#         ng = self.material.get_group_refractive_index(l)
#         data[:, 1] = n
#         data[:, 2] = ng
#         nRays = x.shape[0]
#
#         # Transform to local coordinate system:
#         xLocal = np.transpose(np.dot(self.xpM, np.dot(self.xM, np.transpose(x))))
#         xpLocal = np.transpose(np.dot(self.xpM, np.transpose(xp)))
#
#         # Spherical surface
#         # Intersection where xLocal + t*xpLocal crosses |x-xc|^2 = r^2
#         #
#         xc = np.reshape(np.tile(np.array([0, 0, -self.r, 1]), nRays), (nRays, 4))
#         xDiff = xLocal[:, 0:3] - xc[:, 0:3]
#
#         b = 2 * np.sum(np.multiply(xDiff, xpLocal[:, 0:3]), 1)
#         c = np.sum(np.multiply(xDiff, xDiff), 1) - self.r ** 2
#         sq2 = np.multiply(b, b) / 4 - c
#         sq2PosInd = sq2 > 0         # Intersection exists for rays with indexes in sq2PosInd
#         sq = np.sqrt(sq2[sq2PosInd])
#         t = -b[sq2PosInd] / 2 + sq
#         tNegInd = xLocal[sq2PosInd, 2] + np.multiply(t, xpLocal[sq2PosInd, 2]) + self.r < 0
#         t[tNegInd] -= 2 * sq[tNegInd]
#
#         xNewLocal = xLocal[sq2PosInd, :] + np.transpose(np.multiply(np.transpose(xpLocal[sq2PosInd, :]), t))
#         intersectInd = self.aperture.point_in_aperture(xNewLocal)
#         xnLocal = xNewLocal[intersectInd, :] - xc[sq2PosInd, :][intersectInd, :]
#         xnLocal = np.transpose(np.divide(np.transpose(xnLocal) , np.sqrt(np.sum(np.multiply(xnLocal, xnLocal), 1))))
#         xpNewLocal = self.calculate_local_diffractions(xNewLocal[intersectInd, :],
#                                                        xpLocal[sq2PosInd, :][intersectInd, :], xnLocal, n, n0)
#         xNew = np.transpose(np.dot(self.xMT, np.dot(np.transpose(self.xpM), np.transpose(xNewLocal[intersectInd, :]))))
#         xpNew = np.transpose(np.dot(np.transpose(self.xpM), np.transpose(xpNewLocal)))
#
# #         self.logger.info("\n\nintersect ind: {0}\nsq2PosInd: {1}".format(intersectInd, sq2PosInd))
#         use_ind = np.logical_and(sq2PosInd, intersectInd)
#         rays[use_ind, 0, :] = xNew
#         rays[use_ind, 1, :] = xpNew
#         rays[use_ind, 2, 1] = n
#         rays[use_ind, 2, 2] = ng
#         return rays

    def surface_intersect(self, x_local, xp_local):
        """
        Calculate intersection points and surface normals for rays intersecting spherical surface.

        :param x_local: Rays positions in surface local coordinates
        :type x_local: np.array
        :param xp_local: Rays directions in surface local coordinates
        :type xp_local: np.array
        :return: t, xn_local: distances to intersections, normal directions at intersection
        :rtype:
        """
        t = -x_local[:, 2] / xp_local[:, 2]
        xn_local = np.reshape(np.tile(np.array([0.0, 0.0, 1.0, 0.0]), t.shape[0]), (t.shape[0], 4))

        # Spherical surface
        # Intersection where xLocal + t*xpLocal crosses |x-xc|^2 = r^2
        #
        n_rays = x_local.shape[0]
        xc = np.reshape(np.tile(np.array([0, 0, -self.r, 1]), n_rays), (n_rays, 4))
        x_diff = x_local[:, 0:3] - xc[:, 0:3]

        b = 2 * np.sum(np.multiply(x_diff, xp_local[:, 0:3]), 1)
        c = np.sum(np.multiply(x_diff, x_diff), 1) - self.r ** 2
        sq2 = np.multiply(b, b) / 4 - c
        sq2_pos_ind = sq2 > 0  # Intersection exists for rays with indexes in sq2_pos_ind
        sq = np.sqrt(sq2[sq2_pos_ind])
        t = -b[sq2_pos_ind] / 2 + sq
        t_neg_ind = x_local[sq2_pos_ind, 2] + np.multiply(t, xp_local[sq2_pos_ind, 2]) + self.r < 0
        t[t_neg_ind] -= 2 * sq[t_neg_ind]

        x_new_local = x_local[sq2_pos_ind, :] + np.transpose(np.multiply(np.transpose(xp_local[sq2_pos_ind, :]), t))
        xn_local = x_new_local - xc[sq2_pos_ind, :]
        xn_local = np.transpose(np.divide(np.transpose(xn_local), np.sqrt(np.sum(np.multiply(xn_local, xn_local), 1))))

        return t, xn_local

    def generate_surface_edge(self):
        nbrPoints = 16
        theta0 = np.arcsin(self.aperture.size / self.r)
        theta = np.linspace(-theta0, theta0, nbrPoints)
        xe = np.vstack((-self.r * np.sin(theta), np.zeros(nbrPoints), -self.r * (1 - np.cos(theta)), np.ones(nbrPoints)))
        ye = np.vstack((np.zeros(nbrPoints), -self.r * np.sin(theta), -self.r * (1 - np.cos(theta)), np.ones(nbrPoints)))
        nan = np.array([np.nan, np.nan, np.nan, 1.0])
        self.surface_edge = np.transpose(np.hstack((xe, nan[:, np.newaxis], ye)))

