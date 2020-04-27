# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import numpy as np
import haazLib_v2 as hz

# Have a look at what dP/dV look like for the pressure field... 
# then have a look at what E_source looks like. 


# Variables to define. 
adiabatic_compressibility = 0.4559 # this is the adiabatic compressibility of water at 20 degrees C. 
real_only = False  # this uses only the real part of the phasor, or if False, the whole complex phasor. 
outfile = "data.npz"

em_sim = s4l.document.AllSimulations['LeadFieldSimulation']
if em_sim.HasResults():
	simulation_extractor = em_sim.Results()
	em_sensor_extractor = simulation_extractor["Overall Field"]
	em_output = em_sensor_extractor.Outputs["EM E(x,y,z,f0)"]
	em_output.Update()
else: 
	print 'em sim has no results'
	
acoustic_sim = s4l.document.AllSimulations['Ultrasound Transducer']	
if acoustic_sim.HasResults():	
	simulation_extractor = acoustic_sim.Results()
	acoustic_sensor_extractor = simulation_extractor["Overall Field"]
	p_output = acoustic_sensor_extractor.Outputs["p(x,y,z,f)"]
	p_output.Update()
else: 
	print 'p sim has no results'
		
# Get Field data
em_field = em_output.Data.Field(0) # 0 corresponds to snapshot number (typically frequency or time)
# Note: Field is 1D array
# Get just the real part... 
if real_only:
	em_field = np.real(em_field)
	#print em_field_real.shape

# Turn to E field into 3D array
grid_size = em_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) #print field_grid_size # 76, 74, 42
em_field_3d = em_field.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1,3,order='F')
x_axis = em_output.Data.Grid.XAxis
y_axis = em_output.Data.Grid.YAxis
z_axis = em_output.Data.Grid.ZAxis
em_x_axis2 = (x_axis[1:]+x_axis[:-1])/2.0
em_y_axis2 = (y_axis[1:]+y_axis[:-1])/2.0
em_z_axis2 = (z_axis[1:]+z_axis[:-1])/2.0
# Note: Careful with Axes since S4L have Fields evaluated at cell centers
# And Axes represent the nodes
# So field_size != x * y * z, it's (x-1) * (y-1) * (z-1)

# Get Acoustic Field data
p_field = p_output.Data.Field(0) 
if real_only:
	p_field = np.real(p_field)
# Note: Field is 1D array
p_grid_size = p_output.Data.Grid.Dimensions
p_field_3d = p_field.reshape( p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')
p_x_axis = p_output.Data.Grid.XAxis
p_y_axis = p_output.Data.Grid.YAxis
p_z_axis = p_output.Data.Grid.ZAxis
p_x_axis2 = (p_x_axis[1:]+p_x_axis[:-1])/2.0
p_y_axis2 = (p_y_axis[1:]+p_y_axis[:-1])/2.0
p_z_axis2 = (p_z_axis[1:]+p_z_axis[:-1])/2.0


print p_field_3d[:,:,47]  # they are all complex numbers. 
print p_field_3d.shape


# interpolate the p_field into the grid of the e_field. p output has to be a 3d array aligned with the grid. 
p_inEMgrid = hz.interpn(p_x_axis2, p_y_axis2, p_z_axis2, p_field_3d, em_x_axis2, em_y_axis2, em_z_axis2)
# Make P_inemgrid  trivial producer, so that we can then pass it to the gradient calculator. 
P_EMgrid_producer = hz.visualizeArray(p_inEMgrid, em_x_axis2, em_y_axis2, em_z_axis2,name="p_inEMgrid",unit_name="dP/dV",unit="W/m^3" )
P_EMgrid_producer.Update()

# Calculate the gradient of the pressure(in the EM grid). 
inputs = [P_EMgrid_producer.Outputs["dP/dV"]]
field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
field_gradient_evaluator.UpdateAttributes()
gradient_output = field_gradient_evaluator.Outputs["grad(dP/dV)"]
gradient_output.Update()
grad_p = gradient_output.Data.Field(0)

# Check the sizes of the matrices ready for E source calculation: 
#print grad_p.shape
#print em_field_real.shape
#print grad_p.shape

# They are both 236208L,3L so use matmul. 
E_source = np.matmul(em_field[:,None,:], adiabatic_compressibility*grad_p[:,:,None]) [:,0]
# then remove the last indice

# Put E source in 3D so that we can visualize it. 
E_source_3d = E_source.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1, order='F')


nanless = np.nan_to_num(E_source_3d)
#print nanless 
#print E_source_3d  # lots of nans... 
E_source_3d = nanless

E_source_producer = hz.visualizeArray(E_source_3d, em_x_axis2, em_y_axis2, em_z_axis2,unit_name="dP/dV",unit="W/m^3",name="E Source" )

np.savez(outfile,E_source_3d=E_source_3d,x = em_x_axis2, y = em_y_axis2, z = em_z_axis2)

print 'Calculation successful, E Source scalar field shape: ', E_source.shape# -*- coding: utf-8 -*-
