'''
Created on 20 Oct 2015

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

reload(oa)
reload(om)
reload(oe)
reload(rs)
reload(ry)
reload(optsys)
reload(optsurf)

optSys = optsys.OpticalSystem()

apex = 58 * np.pi / 180
prismMat = optSys.materialLibrary.getMaterial('sapphire')
prism = oe.PrismElement(x=np.array([0, 0, 50e-3]), n=1.5, apexAngle=apex, sideLength=50e-3, material=prismMat)

prism.setPosition(np.array([-10e-3, 0, 200e-3, 1]))

screen = oe.ScreenElement(x=np.array([0, 0, 400e-3]))
screen.rotateElement(0, 0 * np.pi / 180)

r1 = rs.Collimated1DSource(numRays=3, xDim=1e-3, l=263e-9, color=(0.15, 0.1, 0.75))

optSys.addElement(prism)
optSys.addElement(screen)
#optSys.addRaySource(r1)
print "raytracing..."
dev = []

theta_step = 0.1
theta = np.arange(0,45,theta_step)

prism.rotateElement(0, np.min(theta) * np.pi / 180)
for th in theta:
    prism.rotateElement(0, theta_step * np.pi / 180)
    r1 = rs.Collimated1DSource(numRays=20, xDim=10e-3, l=263e-9, color=(0.15, 0.1, 0.75))
    optSys.raySourceList=[r1]
    optSys.traceSystem()
    dev.append(np.arctan(r1.rayStoreList[0].x[-1][0]/0.4)*180/np.pi)

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

fig = mpl.figure(2)
mpl.clf()
mpl.plot(theta, dev)
mpl.grid()
mpl.draw()
mpl.show()
