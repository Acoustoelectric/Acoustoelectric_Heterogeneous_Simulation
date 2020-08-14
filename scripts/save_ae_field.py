# -*- coding: utf-8 -*-
'''
 Save the acoustoelectric field after the diffusion step
 ,the pressure field and the lead field. 

'''
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import numpy as np
import subprocess
import sys

outfile = "ae_data.npz"
grad_p = []
sim = s4l.document.AllSimulations['Diffusion ']
if sim.HasResults():
	simulation_extractor = sim.Results()
	sensor_extractor = simulation_extractor["Overall Field"]
	diffusion_output = sensor_extractor.Outputs["T(x,y,z)"]
	
	# Adding a new FieldGradientEvaluator
	inputs = [sensor_extractor.Outputs["T(x,y,z)"]]
	field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
	field_gradient_evaluator.UpdateAttributes()
	gradient_output = field_gradient_evaluator.Outputs["grad(T(x,y,z))"]
	gradient_output.Update()
	grad_p = gradient_output.Data.Field(0)  # add a negative sign. 
	
	diffusion_output.Update()
else: 
	print 'em sim has no results'

v_field = diffusion_output.Data.Field(0)
# print (v_field.shape)
grid_size = diffusion_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) #print field_grid_size # 211,208,158
# print (field_grid_size)

x_axis = diffusion_output.Data.Grid.XAxis
y_axis = diffusion_output.Data.Grid.YAxis
z_axis = diffusion_output.Data.Grid.ZAxis
em_x_axis2 = (x_axis[1:]+x_axis[:-1])/2.0
em_y_axis2 = (y_axis[1:]+y_axis[:-1])/2.0
em_z_axis2 = (z_axis[1:]+z_axis[:-1])/2.0
v_field_3d = v_field.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1,order='F')
print ('v field shape',v_field.shape,type(v_field))
V_source_viz = hz.visualizeArray(v_field_3d, em_x_axis2, em_y_axis2, em_z_axis2,unit_name="Potential field",unit="V",name="Acoustoelectric Voltage Field" )

print ('potential example value: ',v_field_3d[150,150,116])
# put the gradient into a 3D form. 
print ('gradient shape',grad_p.shape,type(grad_p))

# get th emagnitude using linalg.norm
E = np.linalg.norm(grad_p,axis = 1) # get the magnitude. 

g_source_3d = E.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1,order='F')
print ('g_source_3d',g_source_3d.shape,type(g_source_3d))

E_source_viz = hz.visualizeArray(g_source_3d, em_x_axis2, em_y_axis2, em_z_axis2,unit_name="E field",unit="V/m",name="Acoustoelectric E field" )


