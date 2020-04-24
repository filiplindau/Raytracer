"""
Created on 10 Oct 2015

@author: Filip Lindau
"""

import numpy as np
import Raytracer.OpticalElement as oe
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import logging


class OpticalSystem(object):
    """ Implements an optical axis on which elements can be placed. The
    axis can be rotated after each element. The element coordinate system
    is with respect to this axis, were then the z coordinate is distance
    along the axis.
    """
    def __init__(self):
        self.elements = list()
        self.ray_source_list = list()
        self.material_library = om.MaterialLibrary()
        self.optical_axis_xpM = [np.identity(4)]
        self.optical_axis_xM = [np.identity(4)]
        self.optical_axis_xMT = [np.identity(4)]
        self.optical_axis_theta = 0.0
        self.optical_axis_phi = 0.0
        self.logger = logging.getLogger("System.")
        self.logger.setLevel(logging.INFO)
        self.surround_aperture = oa.InsideSphereAperture(1.0)
        
    def add_element(self, element):
        # element.rotateElement(self.opticalAxisTheta, self.opticalAxisPhi)
        self.elements.append(element)
        x = element.x
        self.optical_axis_xpM.append(self.optical_axis_xpM[-1].copy())
#        self.opticalAxisXM.append(self.opticalAxisXM[-1].copy())
#        self.opticalAxisXMT.append(self.opticalAxisXMT[-1].copy())
        
        xg = np.dot(self.optical_axis_xMT[-1], np.dot(np.transpose(self.optical_axis_xpM[-1]), np.array([0, 0, x[2], 1])))
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
        
        self.optical_axis_xM.append(newXM)
        self.optical_axis_xMT.append(newXMT)
        self.surround_aperture.r = 2 * xg
        
    def set_optical_axis_pos(self, element_number):
        x = self.elements[element_number].x
#        self.opticalAxisXpM.append(self.opticalAxisXpM[-1].copy())
        
        xg = np.dot(self.optical_axis_xMT[element_number], np.dot(np.transpose(self.optical_axis_xpM[element_number]), np.array([0, 0, x[2], 1])))
        newXM = np.array([[1.0, 0.0, 0.0, -xg[0]],
                             [0.0, 1.0, 0.0, -xg[1]],
                             [0.0, 0.0, 1.0, -xg[2]],
                             [0.0, 0.0, 0.0, 1.0]])

        newXMT = np.array([[1.0, 0.0, 0.0, xg[0]],
                             [0.0, 1.0, 0.0, xg[1]],
                             [0.0, 0.0, 1.0, xg[2]],
                             [0.0, 0.0, 0.0, 1.0]])
        
        self.optical_axis_xM[element_number + 1] = newXM
        self.optical_axis_xMT[element_number + 1] = newXMT

    def rotate_optical_axis_after_element(self, theta, phi, element_number):
        self.optical_axis_theta = theta
        self.optical_axis_phi = phi
        thM = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, +np.cos(theta), np.sin(theta), 0.0],
                         [0.0, -np.sin(theta), np.cos(theta), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                         [0.0, 1.0, 0.0, 0.0],
                         [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                         [0.0, 0.0, 0.0, 1.0]])
         
        self.optical_axis_xpM[element_number + 1] = np.dot(self.optical_axis_xpM[element_number + 1], np.dot(thM, phM))

    def add_ray_source(self, ray_source):
        self.ray_source_list.append(ray_source)
        
    def trace_system(self):
        if self.ray_source_list:
            for raySource in self.ray_source_list:
                for (ind, element) in enumerate(self.elements):    
                    raysT = raySource.rays.copy()
                    raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.dot(self.optical_axis_xM[ind], np.transpose(raysT[:, 0, :]))))
                    raysT[:, 1, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.transpose(raysT[:, 1, :])))
                    raysList = element.propagate_rays(raysT)
                    for rays in raysList:
                        raysT = rays.copy()
                        raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xMT[ind], np.dot(np.transpose(self.optical_axis_xpM[ind]), np.transpose(raysT[:, 0, :]))))
                        raysT[:, 1, :] = np.transpose(np.dot(np.transpose(self.optical_axis_xpM[ind]), np.transpose(raysT[:, 1, :])))
                        raySource.update_rays(raysT)

    def get_rays_footprint(self, ray_source_number, element_number, surface_number):
        sn = 1
        for el in range(element_number):
            sn += self.elements[el].surfaces.__len__()
        if self.ray_source_list:
            rays = self.ray_source_list[ray_source_number].get_rays_pos_array()
            raysT = np.transpose(np.dot(self.optical_axis_xpM[element_number], np.dot(self.optical_axis_xM[element_number], np.transpose(rays[:, sn, :]))))
            return self.elements[element_number].get_rays_footprint(raysT, surface_number)

    def get_element_edges(self, element_number):
        edges = self.elements[element_number].get_edges()
        edges_new = []
        for edge in edges:
            edges_new.append(np.transpose(np.dot(self.optical_axis_xMT[element_number],
                                                np.dot(np.transpose(self.optical_axis_xpM[element_number]), np.transpose(edge)))))
        return edges_new

