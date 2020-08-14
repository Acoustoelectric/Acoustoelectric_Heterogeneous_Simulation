# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import numpy as np
'''
This is a validation test against the matlab focused piston solution of Bessel functions. 

'''

outfile = "acoustic.npz"

acoustic_sim = s4l.document.AllSimulations['Ultrasound Transducer']	
if acoustic_sim.HasResults():	
	simulation_extractor = acoustic_sim.Results()
	acoustic_sensor_extractor = simulation_extractor["Overall Field"]
	p_output = acoustic_sensor_extractor.Outputs["p(x,y,z,f)"]
	p_intensity = acoustic_sensor_extractor.Outputs["Intensity"]
	p_output.Update()
	p_intensity.Update()
else: 
	print 'p sim has no results'
		
# 
# 		
# Get Acoustic Field data
p_field 	= p_output.Data.Field(0) 
p_intensity = p_intensity.Data.Field(0)

# get the real part of the field. 
#p_field = np.real(p_field)
# Note: Field is 1D array
p_grid_size = p_output.Data.Grid.Dimensions
p_field_3d = p_field.reshape( p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')
p_intensity_3d = p_intensity.reshape(p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')
p_x_axis = p_output.Data.Grid.XAxis
p_y_axis = p_output.Data.Grid.YAxis
p_z_axis = p_output.Data.Grid.ZAxis
p_x_axis2 = (p_x_axis[1:]+p_x_axis[:-1])/2.0
p_y_axis2 = (p_y_axis[1:]+p_y_axis[:-1])/2.0
p_z_axis2 = (p_z_axis[1:]+p_z_axis[:-1])/2.0

np.savez(outfile,p=p_field_3d,intensity =p_intensity_3d,x = p_x_axis2, y = p_y_axis2, z = p_z_axis2)

print 'Done!'