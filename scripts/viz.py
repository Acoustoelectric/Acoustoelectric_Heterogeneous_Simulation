# -*- coding: utf-8 -*-
# Read in a .mat file and also read in 
# the saved SIM4LIFE file and compare. 
#
# By Jean Rintoul
# 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib
import matplotlib.cm as cm

mdata 			= loadmat('field_test.mat')
pp 				= mdata['pp']
insta_p 		= mdata['insta_p']
pressure_amp 	= mdata['pressure_amp']
matp_x 			= mdata['xx'][0]
matp_z 			= mdata['z'][0]
mat_frequency 	= mdata['f0']

#c0 is speed of sound in water. 
#a is source radius. 
#radius of curvature. 
#x=[0:0.1:10]*1e-3;
#z=[40:0.2:90]*1e-3;
# 'f0','c0','a','roc','x','z','Nmax','termtol','p');
# axes are x and z. 
#print (insta_p.shape)  # 101,251. 

roc=63e-3;  #Radius of curvature
# load the s4l data file. 
filename = "acoustic_test.npz"
data = np.load(filename)
#print (data.keys())
p_s4l = data['p']
x = data['x']
y = data['y']
z = data['z']
#z = z + roc;  #Radius of curvature  # adjust it so 0 is the location of the transducer. 
# Or better yet, move the focus to 0.063
# 
[resultx,resulty,resultz] = np.where(p_s4l == np.amax(p_s4l))
# print (resulty)
#print (p_s4l.shape) # x = 36,y = 34,z = 18 
maxslice = p_s4l[:,resulty[0],:] # 36, 34. 
print (maxslice.shape)
#print (maxslice.shape,len(x),len(z))
#print (x[resultx],z[resultz])
# coordinates of the maximum result. 
print (x[resultx[0]],y[resulty[0]],z[resultz[0]])

z_shift= roc - z[resultz[0]] 
print (z_shift)
z = z + z_shift
#print (insta_p.shape, matp_x.shape,matp_z.shape)
print (insta_p.max(), pressure_amp.max(), pp.max(), maxslice.max(),p_s4l.max())

# dnu0 = -5.34893231317329e-06
# dnu1 = 1.03329524886784e-06
# 

# p maximum is a complex number...
# maxslice max = 29046634.0
# abs(pp)

# Matlab Focusing Gain = G = Pfocus/Psource = 17.154.
# S$L Focusing Gain 
s4lp_source = 1e6 # 1MPA. 
G = maxslice.max()/s4lp_source
pressure_gain = maxslice/s4lp_source

print ('hello', G)

# Pressure (Pascals) in SIm4Life is 3.05e7
# 
# Now plot the slice contour. 
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
CS = ax.contourf(z,x, pressure_gain, 10,cmap=plt.cm.bone,vmin=-20,vmax=20)
cbar = fig.colorbar(CS)
m = plt.cm.ScalarMappable(cmap=plt.cm.bone)
m.set_array(pressure_gain)
m.set_clim(-20, 20)
#plt.colorbar(m, boundaries=np.linspace(-20, 20, 10))
#cbar = clippedcolorbar(CS,extend='neither')
#plt.grid()
plt.ylim(-0.010, 0.010)
plt.xlim(0.040, 0.090)
ax.set_title('Sim4Life Focusing Factor(G):'+str(G))
#ax.set_xlabel('word length anomaly')
#ax.set_ylabel('sentence length anomaly')
# Make a colorbar for the ContourSet returned by the contourf call.

# cbar.ax.set_ylabel('verbosity coefficient')
# end contour drawing. 
# 
ax2 = fig.add_subplot(2,1,2)
# 2.9e7 pa? 
contourp = ax2.contourf(matp_z,matp_x,insta_p, 10, cmap=plt.cm.bone,vmin=-20,vmax=20)
#cbar = clippedcolorbar(contourp,extend='neither')
cbar = fig.colorbar(contourp, ax = ax2)

ax2.set_xlabel('Axial(m)')
ax2.set_title('Matlab Focusing Factor(G):'+str(pressure_amp.max()))
ax.set_ylabel('Radial(m)')
ax2.set_ylabel('Radial(m)')


#c  clim(-20,20)
#plt.grid()
plt.show()
#norm=mpl.colors.Normalize(vmin=-0.5, vmax=1.5)
#plt.xlim(-10, 10)
#plt.ylim(-2, 2)

