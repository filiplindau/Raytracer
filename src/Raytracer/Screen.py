'''
Created on 9 Oct 2015

@author: Filip Lindau
'''

import Raytracer.OpticalElement as oe
import numpy as np

class Screen(oe.OpticalElement):
    def __init__(self, x = np.array([0,0,0,1]), xn = np.array([0,0,-1,0]), xt = np.array([0,1,0,0])):
        super(Screen, self).__init__(x=x, xn=xn, xt=xt, n=1.0)
        
    def initSurfaces(self):
        self.surfaces = [oe.Surface(x=self.x, xn=self.xn, xt=self.xt, n=self.n)]
        
        