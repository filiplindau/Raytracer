'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import Raytracer as rt
import Raytracer.OpticalElement as oe
import Raytracer.OpticalSystem as optsys
import Raytracer.OpticalSurface as optsurf
import Raytracer.RaySource as rs
import Raytracer.Ray as ry
import Raytracer.Screen as rsc
from Raytracer.Ray import Ray 
import numpy as np
import matplotlib.pyplot as mpl

reload(oe)
reload(rsc)
reload(rs)
reload(ry)
reload(optsys)
reload(optsurf)

optSys = optsys.OpticalSystem()
slab = oe.OpticalElement(n=1.5, thickness = 5e-3, x=np.array([0,0,10e-3]))
slab.setRotation(0, 10*np.pi/180)
prism = oe.PrismElement(x=np.array([0,0,40e-3]), n=1.5, apexAngle=65*np.pi/180, sideLength=25e-3)
prism.setRotation(0, 23*np.pi/180)
screen = rsc.Screen(x=np.array([0,0,80e-3]))
r = rs.Collimated1DSource(numRays=10, xDim=10e-3)

#optSys.addElement(slab)
optSys.addElement(prism)
optSys.addElement(screen)
optSys.setRaySource(r)
optSys.traceSystem()

fig = mpl.figure(1)
mpl.clf()
optSys.raySource.drayRays(1)
mpl.draw()
mpl.show()