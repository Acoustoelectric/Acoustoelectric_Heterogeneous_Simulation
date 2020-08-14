import numpy as np 

filename = "potential_field.npz"

data = np.load(filename)
#print (data.keys())
v_field_3d = data['v_field_3d']
x = data['x']
y = data['y']
z = data['z']

# Calculate Gradient. 
g = np.asarray(np.gradient( v_field_3d,x,y,z))
# grad_p = g.reshape(3,-1,order="F")
g = np.moveaxis(g,0,-1)
print('gradient', g.shape)
print ('running code in python 3!')
outfile = "vgrad.npz"
# Save out the gradient into a folder. 
np.savez(outfile,v_grad = g)







