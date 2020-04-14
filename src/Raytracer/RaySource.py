"""
Created on 9 Oct 2015

@author: Filip Lindau
"""

import numpy as np
import Raytracer.Ray as Ray
import matplotlib.pyplot as mpl
import logging


class RaySource(object):
    def __init__(self, num_rays=1):
        """ Defines a ray source with a set of rays that can be traced.
        Contains a list of rays.
        
        Input
        numRays: Number of rays in the list
        """
        self.num_rays = num_rays
        self.rays = None
        self.ray_store_list = list()
        self.logger = logging.getLogger("RaySource.")
        self.logger.setLevel(logging.INFO)
        
    def generate_rays(self):
        """ Reimplement in subclass
        """
        pass
        
    def draw_rays(self, fig):
        for ray in self.rays:
            ray.drawRay(fig)
        
    def get_ray_points(self):
        """ Get list of points including color to draw
        """
        points = list()
        for ray_store in self.ray_store_list:
            points.append((ray_store.getArrayPoints(), ray_store.color))
        return points
    
    def get_rays_pos_array(self):
        """ Get array with ray coordinates. Used for computing.
        """
        return np.array([rl.x for rl in self.ray_store_list])
    
    def get_rays_time_on_surface(self, surface_number):
        """ Get array of group arrival time for the selected surface.
        Use -1 for last surface.
        """
        return np.array([rl.time[surface_number] for rl in self.ray_store_list])


class Collimated1DSource(RaySource):
    def __init__(self, num_rays=1000, x_dim=1e-3, l=567e-9, W=1, color=(0, 0, 0.85)):
        """ Defines an 1D ray source with a set of rays that can be traced. The rays
        are collimated in the z direction and initialized along the x axis with
        an extent of xDim.
        Contains a list of rays.
        
        Input
        numRays: Number of rays in the list
        xDim: The rays are evenly spaced around -xDim/2, +xDim/2
        l: Wavelength
        W: Energy
        color: Drawing color of rays 
        """
        super(Collimated1DSource, self).__init__(num_rays)
        self.color = color
        self.generate_rays(x_dim, l, W, color)
        
    def generate_rays(self, x_dim=5e-3, l=567e-9, W=1, color=(0, 0, 0.85)):
        """ Generate rays matrix and a list of RayStores. A RayStore keeps track
        of a ray as it moves through the optical system.
        
        rays is 3d matrix:
        index0: ray number
        index1: data type (0=x, 1=xp, 2=optical data)
        index2: 4 element data vector
        
        The optical data is encoded as (wavelength, n, ng, W)
        
        Examples: 
        rays[10, 0, :]... position vector of ray 10 
        rays[ :, 1, :]... direction vectors of all rays
        rays[15, 2, 1]... current refractive index for ray 15
        """
        xx = np.linspace(-x_dim / 2, x_dim / 2, self.num_rays)
        x = np.column_stack((xx, np.zeros(self.num_rays), np.zeros(self.num_rays), np.ones(self.num_rays)))
        xp = np.column_stack((np.zeros(self.num_rays), np.zeros(self.num_rays), np.ones(self.num_rays), np.zeros(self.num_rays)))
        data = np.column_stack((l * np.ones(self.num_rays), 1.0 * np.ones(self.num_rays), 1.0 * np.ones(self.num_rays), W * np.ones(self.num_rays)))
        self.rays = np.dstack((x, xp, data)).swapaxes(1, 2)
        self.ray_store_list = list()
        for rn in range(self.num_rays):
            x = self.rays[rn, 0, :]
            xp = self.rays[rn, 1, :]
            l = self.rays[rn, 2, 0]
            n = self.rays[rn, 2, 1]
            W = self.rays[rn, 2, 3]
            r = Ray.RayStore(x=x, xp=xp, n=n, l=l, W=W, color=color)
            self.ray_store_list.append(r)
            
    def update_rays(self, new_rays):
        """ Update ray data with new rays (in global coordinates)
        Sets rays matrix to newRays matrix and updates the RayStoresList
        """
        d = np.sqrt(np.sum((new_rays[:, 0, 0:3] - self.rays[:, 0, 0:3]) ** 2, 1))
        self.rays = new_rays
        for rn in range(self.num_rays):
            x = self.rays[rn, 0, :]
            xp = self.rays[rn, 1, :]
            l = self.rays[rn, 2, 0]
            n = self.rays[rn, 2, 1]
            ng = self.rays[rn, 2, 2]
            W = self.rays[rn, 2, 3]
            self.ray_store_list[rn].addPos(newX=x, newXp=xp, newN=n, newNg=ng, W=W, distance=d[rn])
