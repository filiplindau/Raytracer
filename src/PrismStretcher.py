# -*- coding:utf-8 -*-
"""
Created on Oct 21, 2015

@author: Filip Lindau
"""
'''
Created on 19 Oct 2015

@author: Filip Lindau
'''
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
import time

mpl.ion()

reload(oa)
reload(om)
reload(oe)
reload(rs)
reload(ry)
reload(optsys)
reload(optsurf)

optSys = optsys.OpticalSystem()

devAngle = 65 * np.pi / 180
prismAngle = 51 * np.pi / 180

prism1Insert = 0e-3
prism2Insert = 27e-3
prism3Insert = 0e-3
prism4Insert = 0e-3

prismOffset = 7e-3

prismMat = optSys.material_library.get_material('sapphire')
apex = 57 * np.pi / 180
lensMat = optSys.material_library.get_material('fs')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apexAngle=apex, sideLength=50e-3, material=prismMat)

prism.set_position(np.array([-prismOffset + prism1Insert, 0, 200e-3, 1]))
prism.rotate_element(0, prismAngle)

prism2 = oe.PrismElement(x=np.array([-prismOffset + 0e-3 + prism2Insert, 0, 300e-3]), n=1.5, apexAngle=apex, sideLength=50e-3, material=prismMat)
# prism2.reverseElement()
# prism2.rotateElement(0, np.pi - prismAngle)
prism2.rotate_element(0, prismAngle)

prism3 = oe.PrismElement(x=np.array([prismOffset + 3e-3 - prism3Insert, 0, 200e-3]), n=1.5, apexAngle=apex, sideLength=50e-3, material=prismMat)
prism3.reverse_element()
prism3.rotate_element(0, np.pi - prismAngle)

prism4 = oe.PrismElement(x=np.array([-prismOffset - 8e-3 + prism4Insert, 0, 1200e-3]), n=1.5, apexAngle=apex, sideLength=50e-3, material=prismMat)
# prism4.reverseElement()
prism4.rotate_element(0, prismAngle)

lens1 = oe.PCXElement(x=np.array([5e-3, 0, 200e-3]), n=1.5, r=0.1, thickness=2.4e-3, material=lensMat, size=12.7e-3)
lens2 = oe.PCXElement(x=np.array([7e-3, 0, 400e-3]), n=1.5, r=0.1, thickness=2.4e-3, material=lensMat, size=12.7e-3)
lens2.flip_element()
                      
screen = oe.ScreenElement(x=np.array([0, 0, 400e-3]))
screen.rotate_element(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(numRays=5, xDim=5e-3, l=263e-9, color=(0.15, 0.1, 0.75))
r2 = rs.Collimated1DSource(numRays=5, xDim=5e-3, l=265e-9, color=(0, 0.5, 0))

axRot = devAngle
optSys.add_element(prism)
optSys.rotate_optical_axis_after_element(0, axRot, 0)
optSys.add_element(lens1)
optSys.add_element(lens2)
optSys.add_element(prism2)
optSys.rotate_optical_axis_after_element(0, axRot, 3)
# optSys.addElement(prism3)
# optSys.rotateOpticalAxisAfterElement(0, -axRot, 2)
# optSys.addElement(prism4)
# optSys.rotateOpticalAxisAfterElement(0, axRot, 3)
optSys.add_element(screen)
optSys.add_ray_source(r1)
optSys.add_ray_source(r2)
optSys.trace_system()
d = np.linspace(50e-3, 3000e-3, 100)
dt = []
# for dist in d:
#     prism2.setPosition(np.array([-prismOffset + 0e-3 + prism2Insert, 0, dist]))
#     optSys.setOpticalAxisPos(3)
#     r1 = rs.Collimated1DSource(numRays=5, xDim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
#     r2 = rs.Collimated1DSource(numRays=5, xDim=10e-3, l=265e-9, color=(0, 0.5, 0))
#     optSys.raySourceList = [r1, r2]
#     print "raytracing..."
#     optSys.traceSystem()
#     t1 = r1.getRaysTimeOnSurface(-1)
#     t2 = r2.getRaysTimeOnSurface(-1)
#     newDt = np.mean(t2 - t1)
#     dt.append(newDt)
#     print "dt: ", newDt

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

t1 = r1.get_rays_time_on_surface(-1)
t2 = r2.get_rays_time_on_surface(-1)
print t2 - t1
