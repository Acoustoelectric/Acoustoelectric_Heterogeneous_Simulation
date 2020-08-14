# -*- coding: utf-8 -*-
"""
Acoustoelectric Effect: 

 This code evaluates the source term ready to solve the Poission equation in 
 the modified thermal solver(by settings Penne's Bioheat Equation to zero in appropriate parts)

 To run this script, you will need python 3 and associated(matpotlib, scipy, numpy) scientific computing libraries installed 
 on your machine, so it can be called through this script. The reason we do this is the numpy version in python 2.7 that
 is currently used in Sim4Life evaluates only the real part of the gradient, and we would prefer it returned as a complex number, 
 otherwise we'll lose all the information as we continue our processing line. 
 

Author: Jean Rintoul  16 July 2020
Modifications and advice from: Haza Montanaro
Mathematical advice from: Esra Neufeld
Thanks to Gails for providing Latte. 


"""
import s4l_v1 as s4l
import s4l_v1.analysis as analysis
import s4l_v1.document as document
import numpy as np
import subprocess
import sys
import haazLib_v2 as hz

# Edit this to your own path: 
basepath = "C:\\Sim4Life\\Jean\\ati-working\\"
codepath = basepath + "scripts\\"
command = 'python '+ codepath + 'exprocess.py'

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
# Turn to E field into 3D array
grid_size = em_output.Data.Grid.Dimensions
field_grid_size = (grid_size[0]-1,grid_size[1]-1, grid_size[2]-1) #print field_grid_size # 76, 74, 42
em_field_3d = em_field.reshape( grid_size[0]-1,grid_size[1]-1,grid_size[2]-1,3,order='F')
x_axis = em_output.Data.Grid.XAxis
y_axis = em_output.Data.Grid.YAxis
z_axis = em_output.Data.Grid.ZAxis
em_x = (x_axis[1:]+x_axis[:-1])/2.0
em_y = (y_axis[1:]+y_axis[:-1])/2.0
em_z = (z_axis[1:]+z_axis[:-1])/2.0
# Note: Careful with Axes since S4L have Fields evaluated at cell centers
# And Axes represent the nodes
# So field_size != x * y * z, it's (x-1) * (y-1) * (z-1)

# Get Acoustic Field data
p_field = p_output.Data.Field(0) 

# Note: Field is 1D array
p_grid_size = p_output.Data.Grid.Dimensions
p_field_3d = p_field.reshape( p_grid_size[0]-1,p_grid_size[1]-1,-1,order='F')
p_x_axis = p_output.Data.Grid.XAxis
p_y_axis = p_output.Data.Grid.YAxis
p_z_axis = p_output.Data.Grid.ZAxis
p_x = (p_x_axis[1:]+p_x_axis[:-1])/2.0
p_y = (p_y_axis[1:]+p_y_axis[:-1])/2.0
p_z = (p_z_axis[1:]+p_z_axis[:-1])/2.0

# save out the raw field data so we can call it through python3. 
emfile = basepath + "em.npz"
np.savez(emfile, em_field_3d = em_field_3d, x=em_x,y=em_y,z=em_z)
pfile = basepath + "p.npz"
np.savez(pfile, p_field_3d = p_field_3d, x=p_x,y=p_y,z=p_z)

# Now, open a subprocess, but read in npz files, 
# otherwise we are piping too much data(it's tens of Gb if you've got a 
# fine resolution grid for the acoustic sim)
child = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
while True:
    out = child.stderr.read(1)
    if out == '' and child.poll() != None:
		# print 'no output or no child'
		break
    if out != '':
		sys.stdout.write(out)
		sys.stdout.flush()
		
print 'External Python 3 Command exprocess Complete'
# Now, read back in the resulting source term. 
# Note here we give it a negative sign, as per the equation. 
diffusion_data = "diffusionsource.npz"
data = np.load(diffusion_data)
diffusion_source = data['diffusion_source']
diffusion_source_producer = hz.visualizeComplexArray(-diffusion_source, p_x, p_y, p_z,unit_name="dP/dV",unit="W/m^3",name="Diffusion Source" )

# Create a source term for just the real part.  
diffusion_data_real = "diffusion_real.npz"
data = np.load(diffusion_data_real)
diffusion_real = data['diffusion_real']
diffusion_real_source_producer = hz.visualizeArray(-diffusion_real, p_x, p_y, p_z,unit_name="dP/dV",unit="W/m^3",name="Real Part Diffusion" )

# Create a source term for just the imaginary part.  
diffusion_data_imag = "diffusion_imag.npz"
data = np.load(diffusion_data_imag)
diffusion_imag = data['diffusion_imag']
diffusion_imag_source_producer = hz.visualizeArray(-diffusion_imag, p_x, p_y, p_z,unit_name="dP/dV",unit="W/m^3",name="Imag Part Diffusion" )

# Automatically add it to the source term. 
inputs = [diffusion_real_source_producer.Outputs["dP/dV"]]
field_snapshot_filter = analysis.field.FieldSnapshotFilter(inputs=inputs)
# field_snapshot_filter.Snapshots.ExistingValues = u"0 s"
# field_snapshot_filter.Snapshots.TargetValue = 0.0, units.Seconds
field_snapshot_filter.UpdateAttributes()
document.AllAlgorithms.Add(field_snapshot_filter)
## Adding a new FieldSnapshotFilter
inputs = [diffusion_imag_source_producer.Outputs["dP/dV"]]
field_snapshot_filter = analysis.field.FieldSnapshotFilter(inputs=inputs)
# field_snapshot_filter.Snapshots.ExistingValues = u"0 s"
# field_snapshot_filter.Snapshots.TargetValue = 0.0, units.Seconds
field_snapshot_filter.UpdateAttributes()
document.AllAlgorithms.Add(field_snapshot_filter)

print 'Calculation successful, Diffusion Source scalar field shape: ', diffusion_source.shape
