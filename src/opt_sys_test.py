"""
Created on 28 Apr 2020

@author: Filip Lindau
"""

import numpy as np
import Raytracer.OpticalElement as oe
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import logging
import matplotlib.pyplot as mpl


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
        self.optical_axis_theta = [0.0]                         # Theta angle compared to previous element
        self.optical_axis_phi = [0.0]                           # Phi angle compared to previous element
        self.optical_axis_pos = [np.array([0.0, 0.0, 0.0])]
        self.logger = logging.getLogger("System.")
        self.logger.setLevel(logging.INFO)
        self.surround_aperture = oa.InsideSphereAperture(1.0)

    def add_element_dist(self, dist, theta=0.0, phi=0.0):
        self.elements.append(dist)
        self.optical_axis_xpM.append(np.identity(4))

        self.optical_axis_theta.append(theta)
        self.optical_axis_phi.append(phi)
        self.optical_axis_pos.append(dist)

    def trace_system(self):
        r_l = [np.array([0.0, 0.0, 0.0, 1.0])]
        rp_l = [np.array([0.0, 0.0, 1.0, 0.0])]
        r_g = [np.array([0.0, 0.0, 0.0, 1.0])]
        M_rot = [np.identity(4)]
        for ind in range(len(self.elements)):
            self.logger.info("Element {0}/{1}: dist {2}".format(ind, len(self.elements), self.elements[ind]))
            theta = self.optical_axis_theta[ind]
            M_th = np.array([[1.0, 0.0, 0.0, 0.0],
                            [0.0, +np.cos(theta), np.sin(theta), 0.0],
                            [0.0, -np.sin(theta), np.cos(theta), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            phi = self.optical_axis_phi[ind]
            M_ph = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [-np.sin(phi), 0.0, np.cos(phi), 0.0],
                            [0.0, 0.0, 0.0, 1.0]])

            M_rot_l = np.dot(M_th, M_ph)
            M_rot.append(np.dot(M_rot[-1], M_rot_l))
            xg = self.elements[ind]
            M_x = np.array([[1.0, 0.0, 0.0, 0],
                            [0.0, 1.0, 0.0, 0],
                            [0.0, 0.0, 1.0, xg],
                            [0.0, 0.0, 0.0, 1.0]])
            M_xi = np.array([[1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, -xg],
                            [0.0, 0.0, 0.0, 1.0]])

            rp_l.append(np.dot(M_rot[-1], rp_l[-1]))
            # r_l.append(np.dot(M_x, np.dot(M_rot_l.transpose(), r_l[-1])))
            r_l.append(np.dot(M_rot_l, np.array([0, 0, -xg, 1])))
            r_g.append(r_g[-1] + np.dot(M_rot[-1], np.array([0, 0, xg, 1])))
        return np.array(r_l), np.array(rp_l), np.array(r_g)


if __name__ == "__main__":
    opt_sys = OpticalSystem()
    opt_sys.add_element_dist(1, 0, 0)
    opt_sys.add_element_dist(2, 10 * np.pi / 180.0, 0)
    opt_sys.add_element_dist(1, 10 * np.pi / 180.0, 0)
    opt_sys.add_element_dist(1, -10 * np.pi / 180.0, 0)
    opt_sys.add_element_dist(1, 0 * np.pi / 180.0, 0)
    r_l, rp, r_g = opt_sys.trace_system()

    mpl.figure(1)
    mpl.clf()
    mpl.plot(r_g[:, 2], r_g[:, 1], "-d")
    mpl.plot(r_l[:, 2], r_l[:, 1], "-o")
    mpl.axis('equal')
    mpl.grid()
    mpl.draw()
    mpl.show()
