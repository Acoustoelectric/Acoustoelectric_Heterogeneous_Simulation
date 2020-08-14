"""
Display

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

folder_prefix ="display"

filename = folder_prefix+"/ae_data.npz" # load the s4l data file. 
data = np.load(filename)
#print (data.keys())
E 	= data['E']
phi = data['phi']
x 	= data['x']
y 	= data['y']
z 	= data['z']
print(E.shape,phi.shape,len(x))

# Trim off all the edges 
E = E[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10),:]
phi = phi[10:(len(x)-10),10:(len(y)-10),10:(len(z)-10)]
x = x[10:(len(x)-10)]
y = y[10:(len(y)-10)]
z = z[10:(len(z)-10)]
print ('maxval comparison phi/E',np.max(phi),np.max(E) )
print(E.shape,phi.shape,len(x))
E_mag = np.linalg.norm(E,axis = 3) # get the magnitude.
# Which slice is the max? 
[resultx,resulty,resultz] = np.where(np.real(E_mag) == np.amax(np.real(E_mag) ) )
print ('max slice of laplace_f is: ',resultx,resulty,resultz)
img_slice= resultz[0]
y_slice = resulty[0]
# Plotting of the result. 
fig = plt.figure(figsize=(10,5))
# ax = fig.add_subplot(1,3,1)
# cs = plt.imshow(np.real(phi[:,:,img_slice])) 
# cbar = fig.colorbar(cs)
# plt.grid()
# ax.set_title('phi (V)')


xm, ym, zm = np.meshgrid(x,y,z)


#ax = fig.add_subplot(1,2,2,projection='3d')

#surf = ax.plot_surface(xm[:,:,img_slice],ym[:,:,img_slice],E_mag[:,:,img_slice])
mappable = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
mappable.set_array(E_mag[:,:,img_slice])
#mappable.set_clim(5, 8) # optional

ax2 = fig.add_subplot(1,2,1)
cs = ax2.imshow(E_mag[:,:,img_slice], cmap=mappable.cmap, norm=mappable.norm, extent=(
    np.min(xm[:,:,img_slice]), np.max(xm[:,:,img_slice]), np.min(ym[:,:,img_slice]), np.max(ym[:,:,img_slice])), interpolation='none')
cbar = fig.colorbar(cs)
ax2.set_title('Radial(XY) view E magnitude (V/m)')

print ('xz magnitude: ',E_mag[:,y_slice,:].shape)
print (xm.shape,zm.shape)
print ('x location: ',np.min(xm[:,y_slice,0]),np.max(xm[:,y_slice,0]))
print ('z location: ',np.min(zm[0,y_slice,:]),np.max(zm[0,y_slice,:]) )

print (np.max(E_mag) )
ax3 = fig.add_subplot(1,2,2)
cs = ax3.imshow(E_mag[:,y_slice,:],cmap=mappable.cmap, norm=mappable.norm, extent=(np.min(x),np.max(x),np.min(z),np.max(z)))
# cs = ax3.imshow(E_mag[:,y_slice,:], cmap=mappable.cmap, norm=mappable.norm, extent=(
#      -np.max(xm[:,y_slice,:]), np.max(xm[:,y_slice,:]),  np.min(zm[:,y_slice,:]),np.max(zm[:,y_slice,:])  ), interpolation='none')

cbar = fig.colorbar(cs)
ax3.set_title('Axial(XZ) view E magnitude (V/m)')

# Quiver plot 
#ax3 = fig.add_subplot(1,3,3,projection='3d')
#skip_num = 5
#skip = (slice(None, None, skip_num),  slice(None, None,skip_num),  slice(None, None, skip_num)  )

# The x, y and z coordinates of the arrow locations (default is tail of arrow; see pivot kwarg)
#xx = x
#yy = y
# zz = z[:,:,img_slice]
# E_mag is the scalar value. 
Ex = E[:,:,img_slice,0]
Ey = E[:,:,img_slice,1]
Ez = E[:,:,img_slice,2]

# print ('shapes here:',Ex.shape,Ey.shape,Ez.shape,xx.shape,yy.shape,zz.shape)
# ax.quiver(xx,yy,zz, Ex, Ey, Ez,length=0.1, normalize=True)
# skip3 = (slice(None, None, skip_num), slice(None, None, skip_num), slice(None, None, skip_num))
# # Make the direction data for the arrows. 
# u = np.sin(np.pi * xm) * np.cos(np.pi * ym) * np.cos(np.pi * zm)
# v = -np.cos(np.pi * xm) * np.sin(np.pi * ym) * np.cos(np.pi * zm)
# w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * xm) * np.cos(np.pi * ym) *
#      np.sin(np.pi * zm))

#ax = fig.add_subplot(1,2,2)
#cs = plt.imshow(np.real(E_mag[:,:,img_slice])) 
# ax = fig.add_subplot(1,2,2,projection='3d')
# print (len(xl[::skip_num]),len(yl[::skip_num]),len(zl[::skip_num]), Ex[skip].shape , Ey[skip].shape , Ez[skip].shape)
# I could change the length to the E_mag
# ax.quiver(x[::skip_num], y[::skip_num],z[::skip_num], Ex[skip], Ey[skip], Ez[skip],length=0.1, normalize=True)
#cbar = fig.colorbar(cs)
#plt.grid()
#ax.set_title('E magnitude (V/m)')

fig.suptitle('Acoustoelectric phi and E') # or plt.suptitle('Main title')
fig.subplots_adjust(hspace=0.4,wspace=0.4)
plt.show()

