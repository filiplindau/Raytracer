'''
Created on 12 Oct 2015

@author: Filip Lindau
'''

import numpy as np

class OpticalMaterial(object):
    def __init__(self, name, sellmeierBList, sellmeierCList, n0 = 1.0):
        ''' A new material defined by the sellmeier coefficients B and C.
        Note that the coefficients should be for a wavelength in um (as is
        customary), but the refractive index is then calculated for a wavelength
        in m.
        '''
        self.name = name
        
        self.generateSellmeierFunction(sellmeierBList, sellmeierCList, n0)
        
    def generateSellmeierFunction(self, BList, CList, n0 = 1.0):
        B = np.array(BList)
        C = np.array(CList)
        def calc(l):
            ll = np.ones(B.shape[0])*l*l*1e12            
            n2 = n0 + np.sum(B*ll/(ll-C))
            return np.sqrt(n2)
        self.sellmeier = calc 
        
    def getRefractiveIndex(self, l):
        ''' Reimplement to implement different optical materials 
        '''        
        n = self.sellmeier(l)
        return n
    
    def getGroupRefractiveIndex(self, l):
        dl = l * 1e-5
        n = self.getRefractiveIndex(l)
        dn = self.getRefractiveIndex(l+dl)-self.getRefractiveIndex(l-dl)
        return n-l*dn/(2*dl)
    
class MaterialLibrary(object):
    def __init__(self):
        self.materials = {}
        self.generateMaterials()
        
    def generateMaterials(self):
        # Data from refractiveindex.info
        air = OpticalMaterial('air', [0.0002433468, 2.927321e-5], [0.00420135, 0.0174331])
        self.materials['air'] = air
        bk7 = OpticalMaterial('bk7', [1.03961212, 0.231792344, 1.01046945], [0.00600069867, 0.0200179144, 103.560653])
        self.materials['bk7'] = bk7
        fs = OpticalMaterial('fs', [0.6961663, 0.4079426, 0.8974794], [0.0684043**2, 0.1162414**2, 9.896161**2])
        self.materials['fs'] = fs
        mgf2 = OpticalMaterial('mgf2', [0.60967, 0.0080, 2.14973], [0.08636**2, 18.0**2, 25.0**2], 1.27620)
        self.materials['mgf2'] = mgf2
        caf2 = OpticalMaterial('caf2', [0.69913, 0.11994, 4.35181], [0.09374**2, 21.18**2, 38.46**2], 1.33973)
        self.materials['caf2'] = caf2
        sf11 = OpticalMaterial('sf11', [1.73759695, 0.313747346, 1.89878101], [0.013188707, 0.0623068142, 155.23629])
        self.materials['sf11'] = sf11
        n = OpticalMaterial('n',[],[],1.5)
        self.materials['n'] = n
        
    def getMaterial(self, name):
        return self.materials[name]
    
    def addMaterial(self, material):
        self.materials[material.name] = material