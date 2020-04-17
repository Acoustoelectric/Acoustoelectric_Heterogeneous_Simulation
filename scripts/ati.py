# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import numpy as np
import os 
import haazLib_v2 as hz

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
#print em_field.shape # 236208L

# Get just the real part... 
em_field_real = np.real(em_field)
print em_field_real.shape

# Turn to 3D array
grid_size = em_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) 
print field_grid_size # 76, 74, 42

#em_field_3d = em_field.reshape(field_grid_size,order='F')
em_field_3d = em_field_real.reshape( grid_size[0]-1,grid_size[1]-1,-1,order='F')
print em_field_3d.shape  # 76, 74, 126 , this last dimension is a mystery?
# Could it be E, H and ? 
x_axis = em_output.Data.Grid.XAxis
y_axis = em_output.Data.Grid.YAxis
z_axis = em_output.Data.Grid.ZAxis
print len(z_axis)
em_x_axis2 = (x_axis[1:]+x_axis[:-1])/2.0
em_y_axis2 = (y_axis[1:]+y_axis[:-1])/2.0
em_z_axis2 = (z_axis[1:]+z_axis[:-1])/2.0

#To test:
#hz.visualizeArray(em_field_3d[:,:,0:42], x_axis,y_axis,z_axis)

# Note: Careful with Axes since S4L have Fields evaluated at cell centers
# And Axes represent the nodes
# So field_size != x * y * z, it's (x-1) * (y-1) * (z-1)

# Get Field data
p_field = p_output.Data.Field(0) # 0 corresponds to snapshot number (typically frequency or time)
p_field_real = np.real(p_field)
# Note: Field is 1D array

p_x_axis = p_output.Data.Grid.XAxis
p_y_axis = p_output.Data.Grid.YAxis
p_z_axis = p_output.Data.Grid.ZAxis
p_x_axis2 = (p_x_axis[1:]+p_x_axis[:-1])/2.0
p_y_axis2 = (p_y_axis[1:]+p_y_axis[:-1])/2.0
p_z_axis2 = (p_z_axis[1:]+p_z_axis[:-1])/2.0


# Note: Careful with Axes since S4L have Fields evaluated at cell centers
# And Axes represent the nodes
# So field_size != x * y * z, it's (x-1) * (y-1) * (z-1)

# Now, you will need to interpolate the p_field into the grid of the e_field
# Interpolate field into new grid
print len(p_x_axis),len(p_y_axis),len(p_z_axis)
print len(em_x_axis2), len(em_y_axis2), len(em_z_axis2)
print p_field_real.shape
p_grid_size = p_output.Data.Grid.Dimensions
p_field_3d = p_field_real.reshape( p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')
print p_field_3d.shape  # 118 118 99 

#p_inEMgrid = hz.interpn(p_x_axis, p_y_axis, p_z_axis, p_output, em_x_axis2, em_y_axis2, em_z_axis2)

p_inEMgrid = hz.interpn(p_x_axis2, p_y_axis2, p_z_axis2, p_field_3d, em_x_axis2, em_y_axis2, em_z_axis2)


# P output has to be 3D array and must be aligned with grid
print p_inEMgrid.shape
# In theory this now works??? ...who knows which way around everything is... 
adiabatic_compressibility = 0.4559
p_new = adiabatic_compressibility*p_inEMgrid
print p_new.shape

em_leadfield = em_field_3d[:,:,0:42]
print em_leadfield.shape
# You will need to make the appropriate manipulation 
#E_source = 
#one = em_leadfield 
#two = np.gradient( p_new)[0] # why this zero here? 
print np.gradient( p_new)[0].shape
print em_leadfield.shape
E_source= em_leadfield*np.gradient(p_new)[0]
#E_source= np.dot(em_leadfield,np.gradient( p_new)[0])
print E_source.shape
hz.visualizeArray(E_source, x_axis,y_axis,z_axis)
# Adding a new FieldCalculator
#field_calculator = analysis.field.FieldCalculator(inputs=p_output)
#field_calculator.Expression = u"F"
#field_calculator.ValueLocation = field_calculator.ValueLocation.enum.CellCenter
#field_calculator.UpdateAttributes()
#document.AllAlgorithms.Add(field_calculator)

# Adding a new FieldGradientEvaluator
#inputs = [field_calculator.Outputs["Result(x,y,z)"]]
#field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
#field_gradient_evaluator.UpdateAttributes()
#document.AllAlgorithms.Add(field_gradient_evaluator)

#  How to plot E_source in a slice viewer? 
# 
#em_output.Data.Field(0) = E_source

# Add a slice viewer and plot E_source. 
slice_field_viewer = analysis.viewers.SliceFieldViewer(inputs=E_source)
slice_field_viewer.Data.Mode = slice_field_viewer.Data.Mode.enum.QuantityRealPart
slice_field_viewer.Slice.Plane = slice_field_viewer.Slice.Plane.enum.XZ
slice_field_viewer.Slice.Index = 54
	
print (type(slice_field_viewer))
slice_field_viewer.Visualization.Smooth = True
slice_field_viewer.UpdateAttributes()
document.AllAlgorithms.Add(slice_field_viewer)

	
	