# -*- coding: utf-8 -*-
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import s4l_v1.document as document
import s4l_v1.units as units
import s4l_v1.model as model
from s4l_v1 import Unit
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
# print em_field.shape # 236208L

# Get just the real part... 
em_field_real = np.real(em_field)
#print em_field_real.shape

# Turn to 3D array
grid_size = em_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) 
#print field_grid_size # 76, 74, 42

#em_field_3d = em_field.reshape(field_grid_size,order='F')
#em_field_3d = em_field_real.reshape( grid_size[0]-1,grid_size[1]-1,-1,order='F')
em_field_3d = em_field_real.reshape( grid_size[0]-1,grid_size[1]-1,42,3,order='F')
#print em_field_3d.shape  # 76, 74, 126 , this last dimension is a mystery?
# Could it be E, H and ? 
x_axis = em_output.Data.Grid.XAxis
y_axis = em_output.Data.Grid.YAxis
z_axis = em_output.Data.Grid.ZAxis
#print len(z_axis)
em_x_axis2 = (x_axis[1:]+x_axis[:-1])/2.0
em_y_axis2 = (y_axis[1:]+y_axis[:-1])/2.0
em_z_axis2 = (z_axis[1:]+z_axis[:-1])/2.0

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
p_grid_size = p_output.Data.Grid.Dimensions
p_field_3d = p_field_real.reshape( p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')

#print p_field_3d.shape  # 118 118 99 
#p_inEMgrid = hz.interpn(p_x_axis, p_y_axis, p_z_axis, p_output, em_x_axis2, em_y_axis2, em_z_axis2)
p_inEMgrid = hz.interpn(p_x_axis2, p_y_axis2, p_z_axis2, p_field_3d, em_x_axis2, em_y_axis2, em_z_axis2)
# P output has to be 3D array and must be aligned with grid
#print p_inEMgrid.shape 
adiabatic_compressibility = 0.4559 # this is the adiabatic compressibility of water at 20 degrees C. 
#p_scaled = adiabatic_compressibility*p_inEMgrid
#print p_scaled.shape

# Make P_INemgrid  trivial producer, so that we can then pass it to the gradient calculator. 
P_EMgrid_producer = hz.visualizeArray(p_inEMgrid, em_x_axis2, em_y_axis2, em_z_axis2,name="p_inEMgrid",unit_name="dP/dV",unit="W/m^3" )
P_EMgrid_producer.Update()
# the other elements are the magnetic field and the magnetic vector potential field which are hopefully negligible at low frequencies. 
# em_leadfield = em_field_3d[:,:,0:42]  

inputs = [P_EMgrid_producer.Outputs["dP/dV"]]
field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
field_gradient_evaluator.UpdateAttributes()

gradient_output = field_gradient_evaluator.Outputs["grad(dP/dV)"]
gradient_output.Update()
grad_p = gradient_output.Data.Field(0)



#inputs = []
#inputs = [P_EMgrid_producer.Outputs["dP/dV [p_inEMgrid]"]]
# trivial_producer = analysis.extractors.DataCacheImporter(inputs=inputs)
# trivial_producer.Name = "test_grad"
# trivial_producer.FileName = u"c:\\Users\\jeantoul\\Desktop\\sim4life_tests\\ATI\\scripts\\pemgrid_grad.cache"
# trivial_producer.UpdateAttributes()
# document.AllAlgorithms.Add(trivial_producer)

# Adding a new FieldGradientEvaluator
#inputs = [trivial_producer.Outputs["dP/dV"]]
#field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
#field_gradient_evaluator.UpdateAttributes()
#document.AllAlgorithms.Add(field_gradient_evaluator)

#inputs = []
# Adding a new FieldCalculator
#inputs = [P_EMgrid_producer.Outputs["dP/dV [p_inEMgrid]"]]
#field_calculator = analysis.field.FieldCalculator(inputs=inputs)
#field_calculator.UpdateAttributes()
#document.AllAlgorithms.Add(field_calculator)

# Adding a new FieldGradientEvaluator
#inputs = [field_calculator.Outputs["Result(x,y,z)"]]
#field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
#field_gradient_evaluator.UpdateAttributes()
#document.AllAlgorithms.Add(field_gradient_evaluator)
	
# grad(Result(x,y,z))[Gradient]	
# Extract the gradient matrix. 
#gradient_output = field_gradient_evaluator.Outputs["grad(Result(x,y,z))"]
#gradient_output = trivial_producer.Outputs["grad(dP/dV)"]
#gradient_output.Update()
#grad_p = gradient_output.Data.Field(0)

# print p_inEMgrid.shape
#
print grad_p.shape
print em_field_real.shape
# print grad_p.shape
# They are both 236208L,3L
E_test = np.matmul(em_field_real[:,None,:], adiabatic_compressibility*grad_p[:,:,None]) [:,0]
# then remove the last indice

print E_test.shape
# move it to 3D. 
grid_size = em_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) 
print field_grid_size # 76, 74, 42
E_source_3d = E_test.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1, order='F')
E_source_producer = hz.visualizeArray(E_source_3d, em_x_axis2, em_y_axis2, em_z_axis2,unit_name="dP/dV",unit="W/m^3",name="E Source" )

# Problem, this is the gradient calculated on the pressure grid, and we need it on the EM grid.
# 
# Get gradient through the sim4life interface... 
#inputs = [acoustic_sensor_extractor.Outputs["p(x,y,z,f)"]]  # type list [<'Overall Field', Output 0: 'p(x,y,z,f)'>]
# inputs = [P_EMgrid_producer.Outputs[""]]
#inputs = P_EMgrid_producer.Update()
# field_calculator = analysis.field.FieldCalculator(inputs=inputs)
# field_calculator.Expression = u"F"
# field_calculator.ValueLocation = field_calculator.ValueLocation.enum.CellCenter
# field_calculator.UpdateAttributes()
# document.AllAlgorithms.Add(field_calculator)
	
# Adding a new FieldGradientEvaluator
# inputs = [field_calculator.Outputs["Result(x,y,z)"]]
# field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
# field_gradient_evaluator.UpdateAttributes()
# document.AllAlgorithms.Add(field_gradient_evaluator)

# Extract the gradient matrix. 
# gradient_output = field_gradient_evaluator.Outputs["grad(Result(x,y,z))"]
# gradient_output.Update()
# grad_p = gradient_output.Data.Field(0)
# gradient shape is: 1378476, 3L 
# em_field_real.shape  # 236208
# em_field_3d.shape # 76 74 42 3 
# E_test = np.dot(em_field_real, adiabatic_compressibility*grad_p)
# E_test = em_field_real* adiabatic_compressibility*grad_p

# old fashioned numpy calculation of the real part of the acoustoelectric effect. 
# E_source= em_leadfield*np.gradient(p_scaled)[0]
# print E_source.shape

# E_source_producer = hz.visualizeArray(E_source, x_axis,y_axis,z_axis,unit_name="dP/dV",unit="W/m^3",name="E R Source" )

# E_source_p = hz.visualizeArray(E_source, x_axis,y_axis,z_axis,unit_name="dP/dV",unit="W/m^3",name="E Source Test 2" )

# Adding a new SliceFieldViewer
# inputs = [E_source_producer.Outputs[""]]
# slice_field_viewer = analysis.viewers.SliceFieldViewer(inputs=inputs)
# slice_field_viewer.Data.Mode = slice_field_viewer.Data.Mode.enum.QuantityRealPart
# slice_field_viewer.Slice.Index = 20
# slice_field_viewer.GotoMaxSlice()
# slice_field_viewer.UpdateAttributes()
# document.AllAlgorithms.Add(slice_field_viewer)
# slice_field_viewer.Visible = True

# Now save the E_source producer as... 
# inputs = []
# trivial_producer = analysis.extractors.DataCacheImporter(inputs=inputs)
# trivial_producer.Name = "AE Real Source"
# trivial_producer.FileName = u"c:\\Users\\jeantoul\Desktop\\sim4life_tests\\ATI\\scripts\\E_R_source.cache"
# trivial_producer.UpdateAttributes()
# document.AllAlgorithms.Add(trivial_producer)

# Adding a new DataCacheExporter
# inputs = [trivial_producer.Outputs[""]]
# data_cache_exporter = analysis.exporters.DataCacheExporter(inputs=inputs)
# data_cache_exporter.UpdateAttributes()
# document.AllAlgorithms.Add(data_cache_exporter)

# Get gradient through the sim4life interface... 
# inputs = [acoustic_sensor_extractor.Outputs["p(x,y,z,f)"]]
# field_calculator = analysis.field.FieldCalculator(inputs=inputs)
# field_calculator.Expression = u"F"
# field_calculator.ValueLocation = field_calculator.ValueLocation.enum.CellCenter
# field_calculator.UpdateAttributes()
# document.AllAlgorithms.Add(field_calculator)
	
# Adding a new FieldGradientEvaluator
# inputs = [field_calculator.Outputs["Result(x,y,z)"]]
# field_gradient_evaluator = analysis.core.FieldGradientEvaluator(inputs=inputs)
# field_gradient_evaluator.UpdateAttributes()
# document.AllAlgorithms.Add(field_gradient_evaluator)

# Extract the gradient matrix. 
# gradient_output = field_gradient_evaluator.Outputs["grad(Result(x,y,z))"]
# gradient_output.Update()
# grad = gradient_output.Data.Field(0)
