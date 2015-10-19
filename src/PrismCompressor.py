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

reload(oa)
reload(om)
reload(oe)
reload(rs)
reload(ry)
reload(optsys)
reload(optsurf)

optSys = optsys.OpticalSystem()

prismMat = optSys.materialLibrary.getMaterial('fs')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)

prism.setPosition(np.array([-5e-3, 0, 200e-3, 1]))
prism.rotateElement(0, 23 * np.pi / 180)

prism2 = oe.PrismElement(x=np.array([0, 0, 800e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism2.rotateElement(0, 23 * np.pi / 180)

prism3 = oe.PrismElement(x=np.array([0, 0, 200e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism3.flipElement()
prism3.rotateElement(0, 23 * np.pi / 180)

prism4 = oe.PrismElement(x=np.array([0, 0, 800e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism4.rotateElement(0, 203 * np.pi / 180)
prism4.flipElement()

screen = oe.ScreenElement(x=np.array([0, 0, 200e-3]))
screen.rotateElement(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(numRays=20, xDim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
r2 = rs.Collimated1DSource(numRays=20, xDim=10e-3, l=264e-9, color=(0, 0.5, 0))

optSys.addElement(prism)
optSys.rotateOpticalAxisAfterElement(0, 42.5 * np.pi / 180, 0)
optSys.addElement(prism2)
optSys.rotateOpticalAxisAfterElement(0, 42.5 * np.pi / 180, 1)
optSys.addElement(prism3)
optSys.rotateOpticalAxisAfterElement(0, 42.5 * np.pi / 180, 2)
optSys.addElement(prism4)
optSys.rotateOpticalAxisAfterElement(0, 42.5 * np.pi / 180, 3)
optSys.addElement(screen)
optSys.addRaySource(r1)
optSys.addRaySource(r2)
print "raytracing..."
optSys.traceSystem()

# Draw trace:
fig = mpl.figure(1)
mpl.clf()

for ind in range(optSys.elements.__len__()):    
    ee = optSys.getElementEdges(ind)
    for edge in ee:
        mpl.plot(edge[:, 2], edge[:, 0], 'k')

for raySource in optSys.raySourceList:
    p = raySource.getRayPoints()
    for rayP in p:
        mpl.plot(rayP[0][:, 2], rayP[0][:, 0], color=rayP[1])
mpl.axis('equal')
mpl.grid()
mpl.draw()
mpl.show()

