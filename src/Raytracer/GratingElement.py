"""
Created on 9 Oct 2015

@author: Filip Lindau
"""

import numpy as np

import Raytracer.OpticalSurface as os
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import Raytracer.OpticalElement as oe
import Raytracer.Ray as rr

import logging

air = om.OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])
ml = om.MaterialLibrary()


class GratingElement(oe.OpticalElement):
    def __init__(self, x=np.array([0, 0, 0, 1]), xn=np.array([0, 0, 1, 0]), xt=np.array([0, 1, 0, 0]), n=1.0,
                 thickness=2e-3, grating_period=1.0/1800e3, m=1, side_length=25e-3, material=ml.get_material("fs")):
        self.grating_period = grating_period
        self.side_length = side_length
        self.thickness = thickness
        self.m = m
        oe.OpticalElement.__init__(self, x=x, xn=xn, xt=xt, n=n, material=material)

    def init_surfaces(self):
        self.logger.info("Init grating surfaces")
        ap = oa.RectangularAperture([self.side_length, self.side_length])
        #        s1 = os.Surface(x=self.x, xn=-self.xn, xt=self.xt, n=self.n, material=self.material, aperture=ap)
        s1 = os.Surface(x=np.array([0, 0, 0, 1]), xn=-self.xn, xt=self.xt, n=self.n,
                        material=self.material, aperture=ap, grating_period=self.grating_period, m=self.m)
        s1.set_rotation_internal(0, 0)

        #        s2 = os.Surface(x=self.x, xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2 = os.Surface(x=np.array([0, 0, self.thickness, 1]), xn=self.xn, xt=self.xt, n=1.0, material=air, aperture=ap)
        s2.set_rotation_internal(0, 0)
        self.surfaces = [s2, s1]
