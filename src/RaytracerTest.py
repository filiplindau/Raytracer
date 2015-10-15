'''
Created on 9 Oct 2015

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

slab = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 100e-3]), material=optSys.materialLibrary.getMaterial('fs'))
slab.setRotation(0, 0 * np.pi / 180)

slab2 = oe.OpticalElement(n=1.5, thickness=10e-3, x=np.array([0, 0, 300e-3]), material=optSys.materialLibrary.getMaterial('fs'))
slab2.setRotation(0, 0 * np.pi / 180)

prismMat = om.OpticalMaterial('n', [], [], 1.1)
prismMat = optSys.materialLibrary.getMaterial('fs')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)

prism.setPosition(np.array([-5e-3, 0, 50e-3, 1]))
prism.setRotation(0, 23 * np.pi / 180)


lens = oe.PCXElement(x=np.array([0, 0, 200e-3]), n=1.5, r=0.1, thickness=5e-3, material=optSys.materialLibrary.getMaterial('fs'), size=12.7e-3)
#lens.rotateElement(0, 47 * np.pi / 180)
lens2 = oe.PCXElement(x=np.array([0, 0, 600e-3]), n=1.5, r=0.1, thickness=5e-3, material=optSys.materialLibrary.getMaterial('fs'))
lens2.flipElement()
#lens2.rotateElement(0, 47 * np.pi / 180.0)

prism2 = oe.PrismElement(x=np.array([-605e-3, 0, 715e-3]), n=1.5, apexAngle=65 * np.pi / 180, sideLength=50e-3, material=prismMat)
prism2.rotateElement(0, 2 * 25 * np.pi / 180)

screen = oe.ScreenElement(x=np.array([0, 0, 500e-3]))
screen.rotateElement(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(numRays=2, xDim=10e-3, l=263e-9, color=(0.35, 0, 0.75))
r2 = rs.Collimated1DSource(numRays=20, xDim=10e-3, l=264e-9, color=(0, 0.5, 0))

optSys.addElement(slab)
optSys.rotateOpticalAxisAfterElement(0, -10 * np.pi / 180, 0)
optSys.addElement(slab2)
#optSys.addElement(prism)

#optSys.addElement(lens)
#optSys.addElement(lens2)
#optSys.addElement(prism2)
optSys.addElement(screen)
optSys.addRaySource(r1)
#optSys.addRaySource(r2)
print "raytracing..."
optSys.traceSystem()

# Draw trace:
fig = mpl.figure(1)
mpl.clf()
for element in optSys.elements:
    ee = element.getEdges()
    for edge in ee:
        mpl.plot(edge[2, :], edge[0, :], 'k')

for raySource in optSys.raySourceList:
    p = raySource.getRayPoints()
    for rayP in p:
        mpl.plot(rayP[0][:, 2], rayP[0][:, 0], color=rayP[1])
#mpl.axis('equal')
mpl.grid()
mpl.draw()
mpl.show()
