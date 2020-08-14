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
from scipy.stats.stats import pearsonr
from scipy import interpolate

mdata 			= loadmat('raw_data_0.5Mhz.mat')
pp 				= mdata['pp']
insta_p 		= mdata['insta_p']
pressure_amp 	= mdata['pressure_amp']
matp_x 			= mdata['xx'][0]
matp_z 			= mdata['z'][0]
mat_frequency 	= mdata['f0']
print (pressure_amp.shape, insta_p.shape, pp.shape)


roc=43e-3;  # Radius of curvature
filename = "acoustic.npz" # load the s4l data file. 
data = np.load(filename)
print (data.keys())
p_s4l_complex = data['p']
p_s4l = np.absolute(p_s4l_complex) # first, calculate the pressure magnitude. 
x = data['x']
y = data['y']
z = data['z']
print ('s4l size:', p_s4l.shape) # 1066 1066 447 
# 
z = z + roc;  #Radius of curvature  # adjust it so 0 is the location of the transducer. 
# Or better yet, move the focus to 0.063
[resultx,resulty,resultz] = np.where(p_s4l == np.amax(p_s4l))
print ('sim4life focus location:',resultx,resulty,resultz) # 532 532 254

[fx,fz] = np.where(pressure_amp == np.amax(pressure_amp)) # 100, 103
print ('focal pt of matlab script',fx,fz)
# Matlab Focusing Gain = G = Pfocus/Psource = 17.154.
mm_multiplier = 1000
print ('offset calculation: ',matp_z[fz],z[resultz], matp_z[fz]-z[resultz])
simulation_offset = mm_multiplier*(matp_z[fz]-z[resultz])[0]

# S4L Focusing Gain 
s4lp_source = 1e6 # 1MPA. 
maxslice = p_s4l.max()
G = maxslice/s4lp_source
pressure_gain = p_s4l/s4lp_source

pressure_phasor = p_s4l_complex/s4lp_source
print ('S4L maxslice and Gain, pressuregain', maxslice, G, pressure_gain.shape)



print (z.shape, pressure_gain[resultx,resulty,:][0].shape)
print ('matlab', matp_z.shape, pressure_amp[fx,:][0].shape )
# Now, let's plot both on one. 
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
plt.plot( mm_multiplier*z, pressure_gain[resultx,resulty,:][0],'r') #70
plt.plot(mm_multiplier*matp_z,pressure_amp[fx,:][0])
#ax.set_title('Axial Component: pearson cross-correlation:'+str(round(pearson_axial,2) ) )
ax.set_title('Pressure Magnitude Z focal offset:'+str(round(simulation_offset,2)) +'mm' )
plt.grid()
ax.legend(['S4l','Matlab'])
#ax.set_xlabel('Distance from base of transducer(mm)')
ax.set_ylabel('Gain')


ax2 = fig.add_subplot(2,1,2)
plt.plot( mm_multiplier*z, pressure_phasor[resultx,resulty,:][0],'r') #70
plt.plot(mm_multiplier*matp_z,pp[fx,:][0])
# # look at the number of pts per wavelength. 
# plt.plot( mm_multiplier*z, pressure_gain[resultx,resulty,:][0],'x') #70
# plt.plot(mm_multiplier*matp_z,insta_p[fx,:][0])
ax2.legend(['S4l','Matlab'])
ax2.set_ylabel('Gain')
ax2.set_xlabel('Distance from base of transducer(mm)')
ax2.set_title('Instantaneous Pressure')
plt.grid()
#plt.title('Pressure magnitude')
plt.show()

#
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
plt.plot(mm_multiplier*x, pressure_gain[:,resulty,resultz],'r') 
plt.plot(mm_multiplier*matp_x,pressure_amp[:,fz])
plt.grid()
ax.legend(['S4l','Matlab'])
#ax.set_xlabel('Radial Distance(mm)')
ax.set_ylabel('Gain')
ax.set_title('Pressure Magnitude')
ax.set_xlabel('Radial Distance(mm)')
#ax2 = fig.add_subplot(2,1,2)
#plt.plot( mm_multiplier*x, pressure_phasor[:,resulty,resultz],'r') #70
#plt.plot(mm_multiplier*matp_x,pp[:,fz])
#ax2.legend(['S4l','Matlab'])

#ax2.set_title('Instantaneous Pressure')
#ax2.set_ylabel('Gain=(focus/source)')
#ax2.set_xlabel('Radial Distance(mm)')
#plt.grid()
plt.subplots_adjust(wspace=0.6)
plt.show()

# 
# Now, plot instantaneous pressure, and determine the number of points. 
# 
fig = plt.figure()
ax2 = fig.add_subplot(1,1,1)
plt.plot( mm_multiplier*z, pressure_phasor[resultx,resulty,:][0],'.r') #70
plt.plot(mm_multiplier*matp_z,pp[fx,:][0])
ax2.legend(['S4l','Matlab'])
ax2.set_ylabel('Gain')
ax2.set_xlabel('Axial Distance(mm)')
ax2.set_title('Instantaneous Pressure: view of irregular grid pattern')
plt.grid()
plt.show()




#print (pressure_gain.shape) # 139, 290
#print (len(x),len(z)) # 139, 290 
#print (insta_p.shape) # 201, 251
# 
# Find the maximum index on the axial and radial locations. 
# 
# z(wher)

# def closest(lst, K): 
#     return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))] 

# s4l_minaxial=closest(z,0.04)
# s4l_maxaxial=closest(z,0.09)
# s4l_minradial=closest(x,-0.01)
# s4l_maxradial=closest(x,0.01)

# #print (s4l_minaxialindex, s4l_maxaxialindex, s4l_minradialindex, s4l_maxradialindex)
# min_axial = z.tolist().index(s4l_minaxial)
# max_axial = z.tolist().index(s4l_maxaxial)
# min_radial = x.tolist().index(s4l_minradial)
# max_radial = x.tolist().index(s4l_maxradial)

# pressures4laxial = pressure_gain[resultx,z.tolist().index(s4l_minaxial):z.tolist().index(s4l_maxaxial) ]
# ztots = z[min_axial:max_axial]
# # interpolate the sim4life so it's the same as matlab. 
# f = interpolate.interp1d(ztots, pressures4laxial,fill_value = 'extrapolate' )
# xnew = matp_z
# ynew = f(xnew)
# pearson_axial = np.corrcoef(ynew, insta_p[100,:])[0, 1]
# print ('axial pearson correlation',  round(pearson_axial,2)  )

# pressures4radial = pressure_gain[x.tolist().index(s4l_minradial):x.tolist().index(s4l_maxradial),resultz]
# xtots = x[min_radial:max_radial]
# xtots= np.reshape(xtots,(len(xtots),1),order="f")
# #print ('lengths',pressures4radial.shape, xtots.shape)
# f = interpolate.interp1d(xtots[:,0], pressures4radial[:,0],fill_value = 'extrapolate' )
# xnew = matp_x
# ynew = f(xnew)
# ynew = np.reshape(ynew,(len(ynew),1),order="f") 
# print ('lengths',ynew.shape, insta_p[:,103].shape)
# pearson_radial = np.corrcoef(ynew[:,0], insta_p[:,103])[0, 1]
# print ('radial pearson correlation',round(pearson_radial,2))

# # Radial component +-0.01
# # Axial component 0.04-0.09
# # compare axial components. 
# newz = np.reshape(z,(1,len(z) ))
# print ('pressure gain shape: ',z.shape, newz.shape, pressure_gain[resultx,:].shape) # 139 ,290 

# fig = plt.figure()
# ax = fig.add_subplot(2,1,1)
# print (newz[0])
# print (pressure_gain[resultx,:][0])
# plt.plot( newz[0], pressure_gain[resultx,:][0] ) #70
# plt.plot(matp_z,insta_p[100,:])
# ax.set_title('Axial Component: pearson cross-correlation:'+str(round(pearson_axial,2) ) )
# plt.grid()
# ax.legend(['S4l','Matlab'])
# ax.set_ylabel('Gain')

# #plt.xlim(0.040, 0.090)

# ax2 = fig.add_subplot(2,1,2)
# plt.plot(x, pressure_gain[:,resultz]) #106
# plt.plot(matp_x,insta_p[:,103])
# ax2.set_title('Radial Component: pearson cross-correlation:' + str(round(pearson_radial,2) ) )
# ax2.legend(['S4l','Matlab'])
# ax2.set_ylabel('Gain')
# plt.xlim(-0.010, 0.010)
# plt.grid()
# plt.show()

# Now plot the slice contour. 
# fig = plt.figure()
# ax = fig.add_subplot(2,1,1)
# CS = ax.contourf(z,x, pressure_gain, 10,cmap=plt.cm.bone,vmin=-20,vmax=20)
# cbar = fig.colorbar(CS)
# m = plt.cm.ScalarMappable(cmap=plt.cm.bone)
# m.set_array(pressure_gain)
# m.set_clim(-20, 20)
# plt.ylim(-0.010, 0.010)
# plt.xlim(0.040, 0.090)
# ax.set_title('Sim4Life Focusing Factor(G):'+str(G))
# ax2 = fig.add_subplot(2,1,2)
# contourp = ax2.contourf(matp_z,matp_x,insta_p, 10, cmap=plt.cm.bone,vmin=-20,vmax=20)
# cbar = fig.colorbar(contourp, ax = ax2)
# ax2.set_xlabel('Axial(m)')
# ax2.set_title('Matlab Focusing Factor(G):'+str(pressure_amp.max()))
# ax.set_ylabel('Radial(m)')
# ax2.set_ylabel('Radial(m)')
# plt.show()


