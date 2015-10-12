'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import Raytracer as rt
import Raytracer.OpticalElement as oe
import Raytracer.OpticalSystem as optsys
import Raytracer.OpticalSurface as optsurf
import Raytracer.OpticalMaterial as om
import Raytracer.RaySource as rs
import Raytracer.Ray as ry
from Raytracer.Ray import Ray 
import numpy as np
import matplotlib.pyplot as mpl

reload(om)
reload(oe)
reload(rs)
reload(ry)
reload(optsys)
reload(optsurf)

optSys = optsys.OpticalSystem()

slab = oe.OpticalElement(n=1.5, thickness = 5e-3, x=np.array([0,0,10e-3]), material = optSys.materialLibrary.getMaterial('bk7'))
slab.setRotation(0, 10*np.pi/180)

prism = oe.PrismElement(x=np.array([0,0,30e-3]), n=1.5, apexAngle=65*np.pi/180, sideLength=25e-3, material = om.OpticalMaterial('n',[],[],1.1))
prism.setRotation(0, 23*np.pi/180)

lens = oe.PCXElement(x=np.array([0,0,100e-3]), n=1.5, r=0.025, thickness=10e-3, material = optSys.materialLibrary.getMaterial('bk7'))
lens2 = oe.PCXElement(x=np.array([0,0,200e-3]), n=1.5, r=0.025, thickness=10e-3, material = optSys.materialLibrary.getMaterial('bk7'))

#lens2.setRotation(0*180*np.pi/180, 0.0)
screen = oe.ScreenElement(x=np.array([0,0,500e-3]))

r = rs.Collimated1DSource(numRays=20, xDim=10e-3)

#optSys.addElement(slab)
#optSys.addElement(prism)
#optSys.rotateOpticalAxis(0, 20*np.pi/180)
optSys.addElement(lens)
optSys.addElement(lens2)

optSys.addElement(screen)
optSys.setRaySource(r)
optSys.traceSystem()

fig = mpl.figure(1)
mpl.clf()
optSys.raySource.drayRays(1)
mpl.draw()
mpl.show()