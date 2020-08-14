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
#folder_prefix ="neumann-0.5conve10-it2000"
#folder_prefix ="neumann-0.5conve09-it2000"
folder_prefix ="neumann-0.5conve08-it1000"
#folder_prefix ="neumann-0.5conve08-it1000"
#folder_prefix ="neumann-0.3gridconve08-it1000"
#folder_prefix ="neumann-0.4conve10-it2000"

filename = folder_prefix+"/ae_data.npz" # load the s4l data file. 
data = np.load(filename)
#print (data.keys())
E 	= data['E']
phi = data['phi']
x 	= data['x']
y 	= data['y']
z 	= data['z']
filename = folder_prefix+"/diffusionsource.npz" # load the s4l data file. 
data = np.load(filename)
diffusion_source = data['diffusion_source']
s_x = data['x']
s_y = data['y']
s_z = data['z']
print ('diffusion shape and phi shape: ',diffusion_source.shape,phi.shape,len(x),len(y),len(z))
print (diffusion_source[20,20,20])

# The laplacian is the divergence of the gradient. 
# KX, KY, KZ = np.meshgrid(x, y, z, indexing='ij')
dx = x[1]-x[0]
dy = y[1]-y[0]
dz = z[1]-z[0]
print ('actual laplace dx/dy/dz: ',dx,dy,dz)
sx=s_x[1]-s_x[0]
sy=s_y[1]-s_y[0]
sz=s_z[1]-s_z[0]
print ('source sx/sy/sz: ',sx,sy,sz)
# acc not specified, 4846
# acc= 2, 4846
# acc= 3, 3907
# acc= 5, 3881
# acc= 8, 3863
# acc= 10, 3858
# acc= 20, 3868
laplace = Laplacian(h=[x,y,z],acc=10)
laplace_f = laplace(phi)
print ('laplace f shape',laplace_f.shape)
# laplace_f[laplace_f > 1500] = 0 
# laplace_f[laplace_f < -1500] = 0 

# 
# It appears like the z-offset is not aligned. 
print ('diffusion source shape: ',diffusion_source.shape)
print  ('laplace shape: ',laplace_f.shape)

# Interpolate the source, into the same grid as phi. 
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

source_inphigrid = interpn(s_x, s_y, s_z, diffusion_source, x, y, z)
print ('source in phi grid shape',source_inphigrid.shape)
source = source_inphigrid
print (source[20,20,20])
# attempt to remove the weird edge effect
laplace_phi = laplace_f[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10)]
source = source[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10)]
print ('maxval comparison laplace/source',np.max(laplace_phi),np.max(source) )

# Now, remove the area around the electrodes through a mask. 
# laplace_f[laplace_f > 1500] = 0 
# laplace_f[laplace_f < -1500] = 0 
# e1 = [130:152,130:152,img_slice]
laplace_phi[190:210,110:152,:] = 0
laplace_phi[55:75,110:152,:] = 0
source[190:210,110:152,:] = 0
source[55:75,110:152,:] = 0



# Which slice is the max? 
[resultx,resulty,resultz] = np.where(np.real(laplace_phi) == np.amax(np.real(laplace_phi) ) )
print ('max slice of laplace_f is: ',resultx,resulty,resultz)
# 
[lx,ly,lz] = np.where(np.real(source) == np.amax(np.real(source) ) )
print ('max slice of source is: ',lx,ly,lz)

# Find the Euclidean distance between matrices. 
distance = np.linalg.norm(laplace_phi-source) #axis=-1
print ('Euclidean distance between matrices: ',distance)
diff = laplace_phi-source

img_slice = lz[0]



## Plotting of the result. 
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,4,1)
cs = plt.imshow(np.real(phi[:,:,img_slice]) )
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('phi(V)')

ax = fig.add_subplot(1,4,2)
cs = plt.imshow(np.real(laplace_phi[:,:,img_slice]) )
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('laplacian of phi')

ax = fig.add_subplot(1,4,3)
cs = plt.imshow(np.real(source[:,:,img_slice]))
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('source')

ax = fig.add_subplot(1,4,4)
cs = plt.imshow(np.real(diff[:,:,img_slice]))
cbar = fig.colorbar(cs)
plt.grid()
ax.set_title('diff')
fig.suptitle('Euclidean Distance between laplacian of phi and source: '+str(distance)) # or plt.suptitle('Main title')
fig.subplots_adjust(hspace=0.4,wspace=0.4)
plt.show()

# Questions: Why is the grid spacing different between source term and phi, 
# 
# It's possible I have a negative the wrong way around somewhere. Check. 

