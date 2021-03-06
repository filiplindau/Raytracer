"""
Created on 28 Apr 2020

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
        self.optical_axis_xM_add = [np.identity(4)]
        self.optical_axis_xMT_add = [np.identity(4)]
        self.optical_axis_theta = [0.0]                         # Theta angle compared to previous element
        self.optical_axis_phi = [0.0]                           # Phi angle compared to previous element
        self.optical_axis_pos = [np.array([0.0, 0.0, 0.0])]
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

        xg = np.dot(self.optical_axis_xMT[-1],
                    np.dot(np.transpose(self.optical_axis_xpM[-1]), np.array([0, 0, x[2], 1])))
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
        self.optical_axis_xMT_add.append(newXMT)
        self.optical_axis_xM_add.append(newXM)
        self.optical_axis_theta.append(self.optical_axis_theta[-1])
        self.optical_axis_phi.append(self.optical_axis_phi[-1])
        self.optical_axis_pos.append(element.x)
        self.surround_aperture.r = 2 * xg

    def set_optical_axis_pos(self, element_number):
        x = self.elements[element_number].x
        #        self.opticalAxisXpM.append(self.opticalAxisXpM[-1].copy())

        xg = np.dot(self.optical_axis_xMT[element_number],
                    np.dot(np.transpose(self.optical_axis_xpM[element_number]), np.array([0, 0, x[2], 1])))
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
        self.optical_axis_theta[element_number+1] = theta
        self.optical_axis_phi[element_number+1] = phi
        self.update_transform_matrix_list2()

    def update_transform_matrix_list(self):
        M = np.identity(4)
        for ind in range(len(self.elements)):
            theta = self.optical_axis_theta[ind]
            thM = np.array([[1.0, 0.0, 0.0, 0.0],
                            [0.0, +np.cos(theta), np.sin(theta), 0.0],
                            [0.0, -np.sin(theta), np.cos(theta), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            phi = self.optical_axis_phi[ind]
            phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            M = np.dot(M, np.dot(thM, phM))
            self.optical_axis_xpM[ind] = M

    def update_transform_matrix_list2(self):
        xpM = np.identity(4)
        xM = np.identity(4)
        xMT = np.identity(4)
        for ind in range(len(self.elements)):
            theta = self.optical_axis_theta[ind + 1]
            thM = np.array([[1.0, 0.0, 0.0, 0.0],
                            [0.0, +np.cos(theta), np.sin(theta), 0.0],
                            [0.0, -np.sin(theta), np.cos(theta), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            phi = self.optical_axis_phi[ind + 1]
            phM = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            newXM = np.array([[1.0, 0.0, 0.0, -self.elements[ind].x[0]],
                              [0.0, 1.0, 0.0, -self.elements[ind].x[1]],
                              [0.0, 0.0, 1.0, -self.elements[ind].x[2]],
                              [0.0, 0.0, 0.0, 1.0]])

            newXMT = np.array([[1.0, 0.0, 0.0, self.elements[ind].x[0]],
                               [0.0, 1.0, 0.0, self.elements[ind].x[1]],
                               [0.0, 0.0, 1.0, self.elements[ind].x[2]],
                               [0.0, 0.0, 0.0, 1.0]])

            xpM_tmp = np.dot(thM, phM)
            xM_tmp = np.dot(newXM, xpM_tmp)
            xMT_tmp = np.dot(np.transpose(xpM_tmp), newXMT)

            xpM = np.dot(xpM, xpM_tmp)
            xM = np.dot(xM, xM_tmp)
            xMT = np.dot(xMT, xMT_tmp)
            self.optical_axis_xpM[ind + 1] = xpM
            self.optical_axis_xM_add[ind + 1] = xM
            self.optical_axis_xMT_add[ind + 1] = xMT

    def add_ray_source(self, ray_source):
        self.ray_source_list.append(ray_source)

    def trace_system(self):
        if self.ray_source_list:
            for raySource in self.ray_source_list:
                M = [np.identity(4)]
                for (ind, element) in enumerate(self.elements):
                    M.append(np.dot(np.dot(self.optical_axis_xpM[ind], self.optical_axis_xM[ind]), M[-1]))
                    raysT = raySource.rays.copy()
                    raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.dot(self.optical_axis_xM[ind],
                                                                                            np.transpose(
                                                                                                raysT[:, 0, :]))))
                    raysT[:, 1, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.transpose(raysT[:, 1, :])))
                    raysList = element.propagate_rays(raysT)
                    for rays in raysList:
                        raysT = rays.copy()
                        raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xMT[ind],
                                                             np.dot(np.transpose(self.optical_axis_xpM[ind]),
                                                                    np.transpose(raysT[:, 0, :]))))
                        raysT[:, 1, :] = np.transpose(
                            np.dot(np.transpose(self.optical_axis_xpM[ind]), np.transpose(raysT[:, 1, :])))
                        raySource.update_rays(raysT)

    def trace_system2(self):
        if self.ray_source_list:
            for raySource in self.ray_source_list:
                M = [np.identity(4)]
                for (ind, element) in enumerate(self.elements):
                    M.append(np.dot(np.dot(self.optical_axis_xpM[ind], self.optical_axis_xM[ind]), M[-1]))
                    raysT = raySource.rays.copy()
                    # raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.dot(self.optical_axis_xM[ind],
                    #                                                                         np.transpose(
                    #                                                                             raysT[:, 0, :]))))
                    raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xM_add[ind], np.transpose(raysT[:, 0, :])))
                    raysT[:, 1, :] = np.transpose(np.dot(self.optical_axis_xpM[ind], np.transpose(raysT[:, 1, :])))
                    raysList = element.propagate_rays(raysT)
                    for rays in raysList:
                        raysT = rays.copy()
                        raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xMT_add[ind],
                                                             np.transpose(raysT[:, 0, :])))
                        # raysT[:, 0, :] = np.transpose(np.dot(self.optical_axis_xMT[ind],
                        #                                      np.dot(np.transpose(self.optical_axis_xpM[ind]),
                        #                                             np.transpose(raysT[:, 0, :]))))
                        raysT[:, 1, :] = np.transpose(
                            np.dot(np.transpose(self.optical_axis_xpM[ind]), np.transpose(raysT[:, 1, :])))
                        raySource.update_rays(raysT)

    def aim_ray(self, element_number, pos):
        """
        Aim ray to pass through position on element.

        :param element_number:
        :type element_number:
        :param pos:
        :type pos:
        :return:
        :rtype:
        """

        pass

    def get_rays_footprint(self, ray_source_number, element_number, surface_number):
        sn = 1
        for el in range(element_number):
            sn += self.elements[el].surfaces.__len__()
        if self.ray_source_list:
            rays = self.ray_source_list[ray_source_number].get_rays_pos_array()
            raysT = np.transpose(np.dot(self.optical_axis_xpM[element_number],
                                        np.dot(self.optical_axis_xM[element_number], np.transpose(rays[:, sn, :]))))
            return self.elements[element_number].get_rays_footprint(raysT, surface_number)

    def get_element_edges(self, element_number):
        edges = self.elements[element_number].get_edges()
        edges_new = []
        for edge in edges:
            # edges_new.append(np.transpose(np.dot(self.optical_axis_xMT[element_number],
            #                                      np.dot(np.transpose(self.optical_axis_xpM[element_number]),
            #                                             np.transpose(edge)))))
            # edges_new.append(np.transpose(np.dot(self.optical_axis_xMT_add[element_number],
            #                                      np.dot(np.transpose(self.optical_axis_xpM[element_number]),
            #                                             np.transpose(edge)))))
            edges_new.append(np.transpose(np.dot(self.optical_axis_xMT_add[element_number], np.transpose(edge))))
        return edges_new

