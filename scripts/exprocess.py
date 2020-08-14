# -*- coding: utf-8 -*-
"""
  Acoustoelectric Effect Simulation: Here we perform the bulk of the computation. 


  Since we are computing a less dense grid into a more dense grid, needed to create an accurate diffusion simulation, 
  we processing a lot of data. This script takes some time if it runs on your computer at all , and your output file, particularly if you have a fine acosutic simulation may be large. 

  	Author: Jean Rintoul  16 July 2020
	Modifications and advice from: Haza Montanaro
	Mathematical advice from: Esra Neufeld
	Thanks to Gails for providing Latte. 

"""
import numpy as np 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os.path

# Change this to the adiabatic compressibility of the medium you are propagating your acoustic wave through. 
adiabatic_compressibility = 1e-9  # inverse of the bulk compression. 

def interpn(*args, **kw):
	"""
	Interpolation on N-D. 
	ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
	where the arrays x, y, z, ... define a rectangular grid
	and a.shape == (len(x), len(y), len(z), ...) are the values
	interpolate at xi, yi, zi, ...

	https://github.com/JohannesBuchner/regulargrid/blob/master/regulargrid/interpn.py
	"""
	
	# method = kw.pop('method', 'cubic')
	# method = kw.pop('method', 'nearest')
	method = kw.pop('method', 'linear')
	# method = kw.pop('method', 'slinear')
	if kw:
		raise ValueError("Unknown arguments: " % kw.keys())
	nd = (len(args)-1)//2
	if len(args) != 2*nd+1:
		raise ValueError("Wrong number of arguments")
	q = args[:nd]
	qi = args[nd+1:]
	a = args[nd]
	for j in range(nd):
		# a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
		# a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False, fill_value=np.nan)(qi[j])
		a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False)(qi[j])
	return a

# p_filename = "p.npz"
# data = np.load(p_filename)
# #print (data.keys())
# p_field_3d 	= data['p_field_3d']
# p_x 		= data['x']
# p_y 		= data['y']
# p_z 		= data['z']

if (os.path.isfile('p_grad.npz')) == False:

	p_filename = "p.npz"
	data = np.load(p_filename)
	# print (data.keys())
	p_field_3d 	= data['p_field_3d']
	p_x 		= data['x']
	p_y 		= data['y']
	p_z 		= data['z']
	# print (type(p_field_3d),p_field_3d.shape, len(p_x),len(p_y),len(p_z))
	# Calculate the gradient of the pressure field. 
	grad_p = np.asarray(np.gradient(p_field_3d, p_x,p_y,p_z  ))
	#grad_p = np.gradient(p_field_3d,p_x,p_y,p_z)
	#print ('gradient calculated 1')
	#np.savez('gradientpreprocessed.npz',grad_p=grad_p,p_x=p_x,p_y=p_y,p_z=p_z)
	#print ('gradient complete')

	#g_filename = "gradientpreprocessed.npz"
	#data = np.load(g_filename)
	#grad_p 	= data['grad_p']
	#p_x 		= data['p_x']
	#p_y 		= data['p_y']
	#p_z 		= data['p_z']
	print ('first grad p',grad_p.shape)
	# flatten grad p ready for dot product. 
	grad_p = grad_p.reshape(3,-1,order="F")
	print ('first grad p',grad_p.shape)
	grad_p = np.moveaxis(grad_p,0,-1)
	outfile = "p_grad.npz"
	np.savez(outfile,grad_p=grad_p,p_x=p_x,p_y=p_y,p_z=p_z)
	print ('gradient complete')
else: # read in the npz. 
	newdata_filename = "p_grad.npz"
	data 	= np.load(newdata_filename)
	grad_p 	= data['grad_p']
	p_x 	= data['p_x']
	p_y 	= data['p_y']
	p_z 	= data['p_z']


if (os.path.isfile('eminpgrid.npz')) == False:
	e_filename = "em.npz"
	data = np.load(e_filename)
	#print (data.keys())
	em_field_3d = data['em_field_3d']
	em_x = data['x']
	em_y = data['y']
	em_z = data['z']
	# interpolate each x,y and z of the E field, into the Pressure grid. 
	# later to do: look into iregular grid interpolation, though for now I have a regular grid. 
	em_inPgrid_x = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,0], p_x, p_y, p_z)
	em_inPgrid_y = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,1], p_x, p_y, p_z)
	em_inPgrid_z = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,2], p_x, p_y, p_z)
	em_inPgrid = np.stack((em_inPgrid_x,em_inPgrid_y,em_inPgrid_z))
	em_inPgrid = np.moveaxis(em_inPgrid,0,-1)
	print (em_inPgrid.shape)
	outfile = "eminpgrid.npz"
	np.savez(outfile,em_inPgrid=em_inPgrid)
else: # read in the npz. 
	newdata_filename = "eminpgrid.npz"
	data = np.load(newdata_filename)
	em_inPgrid = data['em_inPgrid']

# Now we have the interpolated data. 
print ('gradient shape: ', grad_p.shape)  # 609184116, 3 
print ('eminPGrid shape: ', em_inPgrid.shape) # 674 674 1341
print ('single element of grad p: ', grad_p[1000,0]) # Check that the gradient retained complex form. 
em_inPgrid_flattened = em_inPgrid.reshape(len(p_x)*len(p_y)*len(p_z),3,order='F')
print ('flattened em in pgrid shape',em_inPgrid_flattened.shape)

# This matrix multiplication is the key equation. 
diffusion_source = np.matmul(em_inPgrid_flattened[:,None,:], adiabatic_compressibility*grad_p[:,:,None]) [:,0]
print ('completed matmul',diffusion_source.shape)
# Put E source in 3D so that we can visualize it. 
diffusion_source_3d = diffusion_source.reshape( len(p_x),len(p_y),len(p_z), order='F')
print ('reshaped',diffusion_source_3d.shape)
nanless = np.nan_to_num(diffusion_source_3d) # de-nan our diffusion source. 
diffusion_source_3d = nanless

print (diffusion_source_3d.shape, diffusion_source_3d[150,150,90])

diffusion_real = np.real(diffusion_source_3d)
diffusion_imag = np.imag(diffusion_source_3d)
#outfile = "C:\Sim4Life\Jean\ATI_Simulation\Acoustic_Simulation\pgrad.npz"
# Save the result out into a folder.  
outfile = "diffusionsource.npz"
np.savez(outfile,diffusion_source = diffusion_source_3d,x=p_x,y=p_y,z=p_z)

outfile = "diffusion_real.npz"
np.savez(outfile,diffusion_real = diffusion_real,x=p_x,y=p_y,z=p_z)

outfile = "diffusion_imag.npz"
np.savez(outfile,diffusion_imag = diffusion_imag,x=p_x,y=p_y,z=p_z)


# array = np.random.randint(0, 9, size=(100, 100, 100))
# new_array = np.zeros((1000, 100, 100))
# x = np.arange(0, 100, 1)
# x_new = np.arange(0, 100, 0.1)
# 
# for i in x:
#     for j in x:
#         f = interp1d(x, array[:, i, j])
#         new_array[:, i, j] = f(xnew)
# 
# To look at it, I'm going to have to plot a contour from a single slice. 
# Check to see if it is a regular or irregular grid. 
# 
# Find the maximum slice. 
# [resultx,resulty,resultz] = np.where(em_inPgrid_x == np.amax(em_inPgrid_x))
# print (resultx,resulty,resultz)
# [presultx,presulty,presultz] = np.where(p_field_3d == np.amax(p_field_3d))
# print (presultx,presulty,presultz)

# This is the XY slice at the focal point. 
# xy_em = em_inPgrid_x[:,:,presultz]
# xy_p = p_field_3d[:,:,presultz]
# 
# print (xy_em[:,:,0].shape, len(p_x),len(p_y))
# 
# Now plot the slice contour. 
# fig = plt.figure()
# ax = fig.add_subplot(2,1,1)
# CS = ax.contourf(p_x,p_y, xy_em[:,:,0], 10,cmap=plt.cm.bone)
# cbar = fig.colorbar(CS)
# ax.set_title('EM Dipole Field at Focus of Ultrasound')
# ax = fig.add_subplot(2,1,2)
# CS = ax.contourf(p_x,p_y, p_field_3d[:,:,0], 10,cmap=plt.cm.bone)
# cbar = fig.colorbar(CS)
# ax.set_title('Instantaneous Pressure at Focus')
# plt.show()

# print (len(p_z),p_field_3d.shape, p_field_3d[presultx,presulty,:].shape )

# fig = plt.figure()
# ax = fig.add_subplot(3,1,1)
# plt.plot(p_z,np.real(p_field_3d[presultx,presulty,:][0]))
# plt.axvline(x=p_z[presultz],color='k', linestyle='--')
# # plt.plot(p_z,np.real(p_field_3d[presultx,presulty,:][0]),'--bo')
# ax.set_title('Axial P Field showing focus slice line')
# plt.grid()
# ax = fig.add_subplot(3,1,2)
# plt.plot(p_x,np.real(p_field_3d[:,presulty,presultz]))
# ax.set_title('Radial Cross-Section')
# plt.grid()

# ax = fig.add_subplot(3,1,3)
# plt.plot(p_x,np.real(em_inPgrid_x[:,presulty,presultz]))
# ax.set_title('EM Field at focus')
# plt.grid()

# plt.subplots_adjust(wspace=0.6)
# plt.show()


# Now plot a cross section, through the middle in 1D 
# fig = plt.figure()
# ax = fig.add_subplot(2,1,1)
# plt.plot(p_x,np.real(xy_em[:,presulty,0]),'--bo')
# ax.set_title('E Field at Focus')
# plt.grid()
# ax = fig.add_subplot(2,1,2)
# plt.plot(p_x,np.real(p_field_3d[:,presulty,0]),'--bo')
# ax.set_title('Instantaneous Pressure at Focus')
# plt.grid()
# plt.show()