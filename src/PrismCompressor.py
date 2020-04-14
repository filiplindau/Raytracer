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
import importlib as imp

imp.reload(oa)
imp.reload(om)
imp.reload(oe)
imp.reload(rs)
imp.reload(ry)
imp.reload(optsys)
imp.reload(optsurf)

optSys = optsys.OpticalSystem()

devAngle = 42.73*np.pi/180
prismAngle = 20.4*np.pi/180

prism1Insert = 0e-3
prism2Insert = 0e-3
prism3Insert = 0e-3
prism4Insert = 0e-3

prismOffset = 14e-3

prismMat = optSys.material_library.get_material('fs')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)

prism.set_position(np.array([-prismOffset + prism1Insert, 0, 200e-3, 1]))
prism.rotate_element(0, prismAngle)

prism2 = oe.PrismElement(x=np.array([prismOffset + 9e-3 - prism2Insert, 0, 1200e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism2.reverse_element()
prism2.rotate_element(0, np.pi - prismAngle)

prism3 = oe.PrismElement(x=np.array([prismOffset + 3e-3 - prism3Insert, 0, 200e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism3.reverse_element()
prism3.rotate_element(0, np.pi - prismAngle)

prism4 = oe.PrismElement(x=np.array([-prismOffset -8e-3 + prism4Insert, 0, 1200e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
#prism4.reverseElement()
prism4.rotate_element(0, prismAngle)

screen = oe.ScreenElement(x=np.array([0, 0, 400e-3]))
screen.rotate_element(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(numRays=5, xDim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
r2 = rs.Collimated1DSource(numRays=5, xDim=10e-3, l=264e-9, color=(0, 0.5, 0))

axRot = devAngle
optSys.add_element(prism)
optSys.rotate_optical_axis_after_element(0, axRot, 0)
optSys.add_element(prism2)
optSys.rotate_optical_axis_after_element(0, -axRot, 1)
optSys.add_element(prism3)
optSys.rotate_optical_axis_after_element(0, -axRot, 2)
optSys.add_element(prism4)
optSys.rotate_optical_axis_after_element(0, axRot, 3)
optSys.add_element(screen)
optSys.add_ray_source(r1)
optSys.add_ray_source(r2)
print "raytracing..."
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

t1 = r1.get_rays_time_on_surface(-1)
t2 = r2.get_rays_time_on_surface(-1)
print t2 - t1
