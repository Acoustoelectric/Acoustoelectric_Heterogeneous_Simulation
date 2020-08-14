# -*- coding: utf-8 -*-
import numpy as np 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os.path

# Read in some NPZ files. 
#emfile = "em.npz"
#np.savez(emfile, em_field_3d = em_field_3d, x=em_x_axis2,y=em_y_axis2,z=em_z_axis2)
#pfile = "p.npz"
#np.savez(pfile, p_field_3d = p_field_3d, x=p_x_axis2,y=p_y_axis2,z=p_z_axis2)
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


p_filename = "p.npz"
data = np.load(p_filename)
#print (data.keys())
p_field_3d 	= np.absolute(data['p_field_3d'])
p_x 		= data['x']
p_y 		= data['y']
p_z 		= data['z']


e_filename = "em.npz"
data = np.load(e_filename)
#print (data.keys())
em_field_3d = data['em_field_3d']
em_x = data['x']
em_y = data['y']
em_z = data['z']

# print (p_field_3d.shape, em_field_3d.shape)
# print (em_field_3d[:,:,:,0].shape)

if (os.path.isfile('eminpgrid.npz')) == False:
	# p in em grid worked, but em in p grid is a bit more tricky. 
	# as the dimensions of em are different. 
	em_inPgrid_x = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,0], p_x, p_y, p_z)
	em_inPgrid_y = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,1], p_x, p_y, p_z)
	em_inPgrid_z = interpn(em_x, em_y, em_z, em_field_3d[:,:,:,2], p_x, p_y, p_z)
	em_inPgrid = np.stack((em_inPgrid_x,em_inPgrid_y,em_inPgrid_z))
	print (em_inPgrid.shape)
	em_inPgrid = np.moveaxis(em_inPgrid,0,-1)

	outfile = "eminpgrid.npz"
	np.savez(outfile,em_inPgrid=em_inPgrid)
else: # read in the npz. 
	newdata_filename = "eminpgrid.npz"
	data = np.load(newdata_filename)
	#print (data.keys())
	em_inPgrid_x = data['em_inPgrid_x']

print (em_inPgrid_x.shape)


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
[presultx,presulty,presultz] = np.where(p_field_3d == np.amax(p_field_3d))
print (presultx,presulty,presultz)

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

print (len(p_z),p_field_3d.shape, p_field_3d[presultx,presulty,:].shape )

fig = plt.figure()
ax = fig.add_subplot(3,1,1)
plt.plot(p_z,np.real(p_field_3d[presultx,presulty,:][0]))
plt.axvline(x=p_z[presultz],color='k', linestyle='--')
# plt.plot(p_z,np.real(p_field_3d[presultx,presulty,:][0]),'--bo')
ax.set_title('Axial P Field showing focus slice line')
plt.grid()
ax = fig.add_subplot(3,1,2)
plt.plot(p_x,np.real(p_field_3d[:,presulty,presultz]))
ax.set_title('Radial Cross-Section')
plt.grid()

ax = fig.add_subplot(3,1,3)
plt.plot(p_x,np.real(em_inPgrid_x[:,presulty,presultz]))
ax.set_title('EM Field at focus')
plt.grid()

plt.subplots_adjust(wspace=0.6)
plt.show()


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