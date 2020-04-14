"""
Created on 9 Oct 2015

@author: Filip Lindau
"""

import Raytracer as rt
import Raytracer.OpticalElement as oe
import Raytracer.OpticalSystem as optsys
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

optSys = optsys.OpticalSystem()

slab = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 100e-3]),
                         material=optSys.material_library.get_material('fs'), size=12.7e-3)
slab.set_rotation(0, 33 * np.pi / 180)

slab2 = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 300e-3]),
                          material=optSys.material_library.get_material('fs'), size=12.7e-3)
slab2.set_rotation(0, 0 * np.pi / 180)

prismMat = om.OpticalMaterial('n', [], [], 1.1)
prismMat = optSys.material_library.get_material('fs')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apex_angle=65 * np.pi / 180, side_length=50e-3,
                        material=prismMat)

prism.set_position(np.array([-5e-3, 0, 200e-3, 1]))
prism.rotate_element(0, 23 * np.pi / 180)


lens = oe.PCXElement(x=np.array([0, 0, 200e-3]), n=1.5, r=0.1, thickness=5e-3,
                     material=optSys.material_library.get_material('fs'), size=12.7e-3)
# lens.rotateElement(0, 47 * np.pi / 180)
lens2 = oe.PCXElement(x=np.array([0, 0, 370e-3]), n=1.5, r=0.1, thickness=5e-3,
                      material=optSys.material_library.get_material('fs'))
lens2.flip_element()
# lens2.rotateElement(0, 47 * np.pi / 180.0)

prism2 = oe.PrismElement(x=np.array([0, 0, 200e-3]), n=1.5, apex_angle=65 * np.pi / 180, side_length=50e-3,
                         material=prismMat)
prism2.rotate_element(0, 23 * np.pi / 180)

screen = oe.ScreenElement(x=np.array([0, 0, 200e-3]))
screen.rotate_element(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(num_rays=20, x_dim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
r2 = rs.Collimated1DSource(num_rays=20, x_dim=10e-3, l=264e-9, color=(0, 0.5, 0))

# optSys.addElement(slab)
# optSys.rotateOpticalAxisAfterElement(0, -10 * np.pi / 180, 0)
# optSys.addElement(slab2)
optSys.add_element(prism)
optSys.rotate_optical_axis_after_element(0, 42.5 * np.pi / 180, 0)
optSys.add_element(lens)
optSys.add_element(lens2)
optSys.add_element(prism2)
optSys.rotate_optical_axis_after_element(0, 42.5 * np.pi / 180, 3)
optSys.add_element(screen)
optSys.add_ray_source(r1)
optSys.add_ray_source(r2)
logger.info("raytracing...")
optSys.trace_system()

# Draw trace:
fig = mpl.figure(1)
mpl.clf()

for ind in range(optSys.elements.__len__()):    
    ee = optSys.get_element_edges(ind)
    for edge in ee:
        mpl.plot(edge[:, 2], edge[:, 0], 'k')

for raySource in optSys.ray_source_list:
    p = raySource.get_ray_points()
    for rayP in p:
        mpl.plot(rayP[0][:, 2], rayP[0][:, 0], color=rayP[1])
mpl.axis('equal')
mpl.grid()
mpl.draw()
mpl.show()
