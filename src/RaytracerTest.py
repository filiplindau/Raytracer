'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import Raytracer as rt
import Raytracer.OpticalElement as oe
import Raytracer.RaySource as rs
import Raytracer.Ray as ry
from Raytracer.Ray import Ray 
import numpy as np
import matplotlib.pyplot as mpl

reload(oe)
reload(rs)
reload(ry)

lens = oe.OpticalElement(n=1.0, thickness = 5e-3, x=np.array([0,0,10e-3]))
#lens.setRotation(0, 30*np.pi/180)
r = rs.Collimated1DSource(numRays=3, xDim=10e-3)

lens.propagateRays(r.rays)


fig = mpl.figure(1)
mpl.clf()
r.drayRays(1)
mpl.draw()
mpl.show()