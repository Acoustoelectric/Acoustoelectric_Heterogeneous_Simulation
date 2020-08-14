# -*- coding: utf-8 -*-
"""
	Read in the diffusion source, and make it a trivial producer. 

"""
import numpy as np 
import haazLib_v2 as hz

infile = "diffusionsource.npz"
data = np.load(infile)
diffusion_source 	= data['diffusion_source']
x 	= data['x']
y 	= data['y']
z 	= data['z']


diffusion_source_producer = hz.visualizeArray(diffusion_source, x,y,z,unit_name="dP/dV",unit="W/m^3",name="Diffusion Source" )

print 'Calculation successful, Diffusion Source scalar field shape: ', diffusion_source.shape









