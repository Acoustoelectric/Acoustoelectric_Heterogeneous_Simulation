"""
Poisson verification

"""
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib
import matplotlib.cm as cm
from scipy.stats.stats import pearsonr
from scipy import interpolate
from numpy.linalg import norm
from scipy.misc import derivative
from findiff import Gradient, Divergence, Laplacian, Curl
from scipy.interpolate import interp1d

#
# Calculate Gradient two different ways. 
# take the pressure. 
# compare p_grad with the one I compute with central differences here. 

folder_prefix ="gradient"

filename = folder_prefix+"/p_grad.npz" # load the s4l data file. 
data = np.load(filename)
#print (data.keys())
grad_p 	= data['grad_p']
x 	= data['p_x']
y 	= data['p_y']
z 	= data['p_z']
# reshape it properly. 
grad_p = grad_p.reshape( len(x),len(y),len(z),-1, order='F')

filename = folder_prefix+"/p.npz" # load the s4l data file. 
data = np.load(filename)
p_field_3d = data['p_field_3d']
p_x = data['x']
p_y = data['y']
p_z = data['z']
print ('p and grad shape: ',p_field_3d.shape,grad_p.shape,len(p_x),len(p_y),len(p_z))


dx = x[1]-x[0]
dy = y[1]-y[0]
dz = z[1]-z[0]
print ('actual dx/dy/dz: ',dx,dy,dz)

grad_p = np.asarray(np.gradient(p_field_3d, x,y,z ))
grad_p = np.moveaxis(grad_p,0,-1)
# the difference is due to small irregularities in the grid. 

grad_test = np.asarray(np.gradient(p_field_3d, dx,dy,dz ))
# grad_test = gradientO4(p_field_3d, dx=dx)
grad_test = np.moveaxis(grad_test,0,-1)
print ('grad test shape: ', grad_test.shape)
# acc = 1, 398591360.0
# acc = 2, 3862.0
# acc = 3, 630521000.0
# acc = 4, 630521000 .0
# acc = 5, 3862
# acc = 10,3862.0
# Compute the gradient using central differences findiff method. 
grad = Gradient(h=[x, y, z],acc=2)
grad_f = grad(p_field_3d)
grad_f = np.moveaxis(grad_f,0,-1)
print ('gradient comparison shapes: ',grad_f.shape,grad_p.shape)
# laplace_f[laplace_f > 1500] = 0 
# laplace_f[laplace_f < -1500] = 0 


# # Now, remove the area around the electrodes through a mask. 
# # laplace_f[laplace_f > 1500] = 0 
# # laplace_f[laplace_f < -1500] = 0 
# # e1 = [130:152,130:152,img_slice]
# laplace_phi[190:210,110:152,:] = 0
# laplace_phi[55:75,110:152,:] = 0
# source[190:210,110:152,:] = 0
# source[55:75,110:152,:] = 0
grad_original   = grad_p[:,:,:,0]
grad_comparison = grad_f[:,:,:,0] #  grad_f[:,:,:,0]

# cut off all the edges
grad_original = grad_original[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10)]
grad_comparison = grad_comparison[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10)]

# Which slice is the max? 
[resultx,resulty,resultz] = np.where(np.real(grad_comparison) == np.amax(np.real(grad_comparison) ) )
print ('max slice of laplace_f is: ',resultx,resulty,resultz)
# 
[lx,ly,lz] = np.where(np.real(grad_original) == np.amax(np.real(grad_original)))
print ('max slice of source is: ',lx,ly,lz)
img_slice = lz[0]


diff = grad_original - grad_comparison
# # Find the Euclidean distance between matrices. 
distance = np.linalg.norm(diff[:,:,img_slice]) #axis=-1
print ('Euclidean distance between matrices: ',distance)


# Plotting of the result. 
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,4,1)
cs = plt.imshow(np.real(p_field_3d[:,:,img_slice])) 
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('pressure field ')

ax = fig.add_subplot(1,4,2)
cs = plt.imshow(np.real(grad_original[:,:,img_slice])) 
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('original gradient')

ax = fig.add_subplot(1,4,3)
cs = plt.imshow(np.real(grad_comparison[:,:,img_slice]))
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('comparison gradient')

ax = fig.add_subplot(1,4,4)
cs = plt.imshow(np.real(diff[:,:,img_slice]))
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('diff, comparison-original')
fig.suptitle('Euclidean Distance between two gradient calculations: '+str(distance)) # or plt.suptitle('Main title')
fig.subplots_adjust(hspace=0.4,wspace=0.4)
plt.show()



