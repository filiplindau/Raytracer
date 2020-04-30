"""
Created on 14 Apr 2020

@author: Filip Lindau
"""

import Raytracer as rt
import Raytracer.OpticalElement as oe
import Raytracer.OpticalSystem2 as optsys
import Raytracer.OpticalSurface as optsurf
import Raytracer.OpticalMaterial as om
import Raytracer.OpticalAperture as oa
import Raytracer.RaySource as rs
import Raytracer.Ray as ry
from Raytracer.Ray import RayStore
import numpy as np
import matplotlib.pyplot as mpl
import importlib as imp
import logging

imp.reload(oa)
imp.reload(om)
imp.reload(oe)
imp.reload(rs)
imp.reload(ry)
imp.reload(optsys)
imp.reload(optsurf)

logger = logging.getLogger()
while len(logger.handlers):
    logger.removeHandler(logger.handlers[0])

f = logging.Formatter("%(asctime)s - %(module)s.   %(funcName)s - %(levelname)s - %(message)s")
fh = logging.StreamHandler()
fh.setFormatter(f)
logger.addHandler(fh)
logger.setLevel(logging.INFO)


opt_sys = optsys.OpticalSystem()

fs = opt_sys.material_library.get_material('fs')

grating = oe.GratingElement(x=np.array([0, 0, 50e-3]), grating_period=260e-9, m=1)
# grating = oe.GratingElement(x=np.array([0, 0, 50e-3]), grating_period=30e-8, m=1)
grating.set_rotation(0, -28.5*np.pi/180)

lens0 = oe.PCXElement(x=np.array([0, 0, 100e-3]), r=250e-3, thickness=4.4e-3, material=fs, size=50.8e-3)
lens1 = oe.PCXElement(x=np.array([0, 0, 100e-3]), r=250e-3, thickness=4.4e-3, material=fs, size=50.8e-3)

slab0 = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 100e-3]),
                          material=opt_sys.material_library.get_material('fs'), size=50e-3)
slab0.set_rotation(0, 0 * np.pi / 180)

slab1 = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 100e-3]),
                          material=opt_sys.material_library.get_material('fs'), size=50e-3)
slab1.set_rotation(0, 10 * np.pi / 180)

disp_plane = oe.ScreenElement(x=np.array([0, 0, 500e-3]))
disp_plane.rotate_element(0, 0 * np.pi / 180)

screen = oe.ScreenElement(x=np.array([0, 0, 500e-3]))
screen.rotate_element(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(num_rays=10, x_dim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
r2 = rs.Collimated1DSource(num_rays=3, x_dim=10e-3, l=265e-9, color=(0, 0.5, 0))
theta = 0 * np.pi / 180.0
m = np.array([[1.0, 0.0, 0.0, 0.0],
              [0.0, +np.cos(theta), np.sin(theta), 0.0],
              [0.0, -np.sin(theta), np.cos(theta), 0.0],
              [0.0, 0.0, 0.0, 1.0]])
phi = 0 * np.pi / 180.0
m = np.array([[+np.cos(phi), 0.0, np.sin(phi), 0.0],
              [0.0, 1.0, 0.0, 0.0],
              [-np.sin(phi), 0.0, np.cos(phi), 0.0],
              [0.0, 0.0, 0.0, 1.0]])
r2.rays[:, 1, :] = np.dot(m, r2.rays[:, 1, :].transpose()).transpose()


opt_sys.add_element(grating)
# opt_sys.add_element(slab0)
opt_sys.add_element(lens0)
# opt_sys.add_element(disp_plane)
# opt_sys.add_element(lens1)
opt_sys.add_element(slab0)
opt_sys.add_element(screen)
opt_sys.rotate_optical_axis_after_element(0, -61 * np.pi / 180, 0)

opt_sys.add_ray_source(r2)
# opt_sys.add_ray_source(r2)
logger.info("raytracing...")
opt_sys.trace_system()

# Draw trace:
fig = mpl.figure(1)
mpl.clf()

for ind in range(opt_sys.elements.__len__()):
    ee = opt_sys.get_element_edges(ind)
    for edge in ee:
        mpl.plot(edge[:, 2], edge[:, 0], 'k')

for raySource in opt_sys.ray_source_list:
    p = raySource.get_ray_points()
    for rayP in p:
        mpl.plot(rayP[0][:, 2], rayP[0][:, 0], color=rayP[1])
mpl.axis('equal')
# mpl.gca().set_aspect("equal")
mpl.grid()
mpl.draw()
mpl.show()
