# -*- coding: utf-8 -*-
import haazLib_v2 as hz
reload(hz) #TODO: Remove
import s4l_v1 as s4l
import s4l_v1.simulation.acoustic as acoustic
import s4l_v1.analysis as analysis
import extract_phasors as phasors
import time
import struct 
import uuid
import h5py 
import numpy as np

def constantTransform(X, c):
	"""
	def constantTransform(X, c)
	"""
	X = np.array(X)
	out = X*0.+c
	return out

def linearTransform(X, x1,x2,y1,y2):
	"""
	def linearTransform(X, x1,x2,y1,y2)
	"""
	X = np.array(X)
	out = y1 + (y2 - y1) * (X - x1) / (x2 - x1)
	return out

def quadraticTransform(X, x1,x2,x3, y1,y2,y3):
	"""
	def quadraticTransform(X, x1,x2,x3, y1,y2,y3)
	"""
	x = np.array([x1, x2, x3])
	y = np.array([y1, y2, y3])
	z = np.polyfit(x, y, 2)

	p = np.poly1d(z)
	return p(X)

def transformmuToHU(X, mu_ref1=-1000, mu_ref2=534.73, HU_ref1=0, HU_ref2=3000):
	"""
	Convert mu to HU Units using reference points
	transformmuToHU(X, mu_ref1=-1000, mu_ref2=534.73, HU_ref1=0, HU_ref2=3000)
	based on water and cortical bone, and minimum and maximum value of extracted CT data
	X - float / list - HU 
	HUmin - Minimum value of CT data
	HUmax- Maximum value of CT data
	pmin - Density of water
	pmax - Maximum density of cortical bone
	return - float / list - speed of sound in m/s
	"""
	return linearTransform(X, mu_ref1, mu_ref2, HU_ref1, HU_ref2)	

def transformmuToHU_referenceStyrofoamEthanol(X, mu_styrofoam=138.0, mu_ethanol_70=2379.0):
	"""
	transformmuToHU_referenceStyrofoamEthanol(X, mu_styrofoam=138.0, mu_ethanol_70=2379.0)
	"""
	hu_ethanol_70 = -124.5
	styrofoam_density = 32.5
	hu_styrofoam = .95*styrofoam_density-1003
	m = 1.*(hu_styrofoam - hu_ethanol_70) / (mu_styrofoam - mu_ethanol_70)
	b = hu_ethanol_70 - m*mu_ethanol_70
	print "m:",m,"b:",b
	return m*X+b

def transformHUToDensity(X, HU_ref1=-1000., HU_ref2=3000., p_ref1=1000., p_ref2=1908.00):
	"""
	Convert HU Units to Density (kg/m^3)
	transformHUToDensity(X, HU_ref1=-1000., HU_ref2=3000., p_ref1=1000., p_ref2=1908.00)
	based on water and cortical bone, and minimum and maximum value of extracted CT data
	X - float / list - HU 
	HUmin - Minimum value of CT data
	HUmax- Maximum value of CT data
	pmin - Density of water
	pmax - Maximum density of cortical bone
	return - float / list - speed of sound in m/s
	"""		
	return linearTransform(X, HU_ref1, HU_ref2, p_ref1, p_ref2)	
	
def transformHUToSpeed_Linear(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=1498., c_ref2=2813.7):
	"""
	Convert HU Units to Speed of sound (linear function)
	based on water and cortical bone
	transformHUToSpeed_Linear(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=1498., c_ref2=2813.7)
	X - float / list - HU
	cmin - Speed of sound of water
	cmax - Max speed of sound of cortical bone
	pmin - Density of water
	pmax - Maximum density of cortical bone
	return - float / list - speed of sound in m/s
	"""		
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	S = linearTransform(P, p_ref1, p_ref2, c_ref1, c_ref2)

	print P
	print S

	return S
	
def transformHUToSpeed_Pichardo_836(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=None, c_ref2=None):
	"""
	Convert HU Units to Speed of sound (cubic spline function)
	Based on Pichardo paper, values taken for 836 kHz
	based on water and cortical bone
	transformHUToSpeed_Pichardo_836(X)
	X - float / list - HU 
	return - float / list - speed of sound in m/s
	"""		
	# Create Spline
	x_rho_c_836 = [1.0000, 1.4500, 1.6920, 2.5180, 2.6420, 3.1430, 3.4000]
	y_c_836 = [1500.0, 2299.0, 2163.0, 1865.0, 5006.0, 2227.0, 3796.0]
	from interpolate_3 import SmoothSpline
	pp = SmoothSpline(x_rho_c_836, y_c_836, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	return pp(P/1000.0)+0

def transformHUToSpeed_Pichardo_270(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=None, c_ref2=None):
	"""
	Convert HU Units to Speed of sound (cubic spline function)
	Based on Pichardo paper, values taken for 836 kHz
	based on water and cortical bone
	transformHUToSpeed_Pichardo_270(X)
	X - float / list - HU 
	return - float / list - speed of sound in m/s
	"""		
	# Create Spline
	x_rho_c_270 = [1.0000, 1.565, 1.697, 2.028, 2.160, 3.045, 3.4000]
	y_c_270 = [1501.0, 1979.0, 2264.0, 1874.0, 3150.0, 3989.0, 5052.0]
	from interpolate_3 import SmoothSpline
	pp = SmoothSpline(x_rho_c_270, y_c_270, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	return pp(P/1000.0)+0

def transformHUToSpeed_Pichardo_836_Linearized(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=None, c_ref2=None):
	"""
	Convert HU Units to Speed of sound (cubic spline function)
	Based on Pichardo paper, values taken for 836 kHz
	based on water and cortical bone
	transformHUToSpeed_Pichardo_836(X)
	X - float / list - HU 
	return - float / list - speed of sound in m/s
	"""		
	# Create Spline
	x_rho_c_836 = [1.0000, 1.4500, 1.6920, 2.5180, 2.6420, 3.1430, 3.4000]
	y_c_836 = [1500.0, 2299.0, 2163.0, 1865.0, 5006.0, 2227.0, 3796.0]
	from interpolate_3 import SmoothSpline
	pp = SmoothSpline(x_rho_c_836, y_c_836, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	S = linearTransform(P/1000., 1.2, 2.6, pp(1.2), pp(2.6))	

	return S

def transformHUToSpeed_Pichardo_500_Interp(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0, c_ref1=None, c_ref2=None):
	"""
	Convert HU Units to Speed of sound (cubic spline function)
	Based on Pichardo paper, values taken for 836 kHz
	based on water and cortical bone
	transformHUToSpeed_Pichardo_836(X)
	X - float / list - HU 
	return - float / list - speed of sound in m/s
	"""		
	f = 500.
	
	from interpolate_3 import SmoothSpline
	# Create Spline 270
	x_rho_c_270 = [1.0000, 1.565, 1.697, 2.028, 2.160, 3.045, 3.4000]
	y_c_270 = [1501.0, 1979.0, 2264.0, 1874.0, 3150.0, 3989.0, 5052.0]
	pp = SmoothSpline(x_rho_c_270, y_c_270, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	out_270 = pp(P/1000.0)+0

	# Create Spline 836
	x_rho_c_836 = [1.0000, 1.4500, 1.6920, 2.5180, 2.6420, 3.1430, 3.4000]
	y_c_836 = [1500.0, 2299.0, 2163.0, 1865.0, 5006.0, 2227.0, 3796.0]
	pp = SmoothSpline(x_rho_c_836, y_c_836, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	out_836 = pp(P/1000.0)+0

	k = (f - 270.) / (836. - 270.)
	out_500 = k * ( out_836 - out_270 ) + out_270

	# print P
	# print out_270
	# print out_836
	# print out_500 	

	return out_500

	
def transformHUToAttenuation_Pichardo_836(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0):
	"""
	Convert HU Units to Attenuation (cubic spline function)
	Based on Pichardo paper, values taken for 836 kHz
	based on water and cortical bone
	transformHUToAttenuation_Pichardo_836(X)
	X - float / list - HU 
	return - float / list - attenuation in Np/m
	"""		
	# Create Spline
	x_rho_att_836 = [1.000, 1.331, 1.638, 1.963, 2.062, 3.028, 3.400]
	y_att_836 = [456.8, 26.59, 290.2, 226.2, 281.1, 249.1, 306.3]
	from interpolate_3 import SmoothSpline
	pp = SmoothSpline(x_rho_att_836, y_att_836, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)

	#H Todo, turn this into a 'visualize function'
	# print "Att 836"
	# x = np.linspace(1200.,2600.,100)
	# for i in range(x.size):
		# print x[i], pp(x[i]/1000.0)+0

	return pp(P/1000.0)+0
	
def transformHUToAttenuation_Pichardo_270(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0):
	"""
	Convert HU Units to Attenuation (cubic spline function)
	Based on Pichardo paper, values taken for 270 kHz
	based on water and cortical bone
	transformHUToAttenuation_Pichardo_270(X)
	X - float / list - HU 
	return - float / list - attenuation in Np/m
	"""		
	# Create Spline
	x_rho_att_270 = [1.000, 1.608, 1.903, 2.018, 3.200, 3.300, 3.400]
	y_att_270 = [158.5, 60.03, 1.0, 1.0, 250.0, 126.2, 40.12]
	from interpolate_3 import SmoothSpline
	pp = SmoothSpline(x_rho_att_270, y_att_270, p=.99)
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)
	
	print "start att pich 270"
	x = np.linspace(1200.,2600.,100)
	for i in range(x.size):
		print x[i], pp(x[i]/1000.0)+0
	print "end att pich 270"
	
	P = np.clip(P, 1200., 2600.) #

	return pp(P/1000.0)+0
	
def transformHUToAttenuation_Pichardo_270_836_mix(X, HU_ref1 = -1000., HU_ref2 = 3000., p_ref1=1000., p_ref2=1908.0):
	"""
	Convert HU Units to Attenuation (cubic spline function)
	Based on Pichardo paper, values taken as average of 270kHz and 836 kHz functions
	based on water and cortical bone
	transformHUToAttenuation_Pichardo_270_836_mix(X)
	X - float / list - HU 
	return - float / list - attenuation in Np/m
	"""		
	from interpolate_3 import SmoothSpline
	P = transformHUToDensity(X, HU_ref1, HU_ref2, p_ref1, p_ref2)

	# Create Spline for 270 kHz
	x_rho_att_270 = [1.000, 1.608, 1.903, 2.018, 3.200, 3.300, 3.400]
	y_att_270 = [158.5, 60.03, 1.0, 1.0, 250.0, 126.2, 40.12]
	pp_270 = SmoothSpline(x_rho_att_270, y_att_270, p=.99)

	# Create Spline for 836 kHz
	x_rho_att_836 = [1.000, 1.331, 1.638, 1.963, 2.062, 3.028, 3.400]
	y_att_836 = [456.8, 26.59, 290.2, 226.2, 281.1, 249.1, 306.3]
	pp_836 = SmoothSpline(x_rho_att_836, y_att_836, p=.99)

	return ( pp_270(P/1000.0)+0 + pp_836(P/1000.0)+0 ) / 2.0

def transformHUToAttenuation_constant(X, att_ref=1.):
	"""
	transformHUToAttenuation_constant(X, att_ref=1.)
	"""
	return X*0.+att_ref

def setDirectFocusingSources(sim_direct_idx, sim_reverse_idx):	
	"""
	Sets phased sources for direct focusing simulation
	setDirectFocusingSources(sim_direct_idx, sim_reverse_idx)
	sim_direct_idx - int / str - simulation idx for direct (output) simulation 
	sim_reverse_idx - int / str - simulation idx for reverse (input) simulation from which phases are extracted
	"""	
	sim_reverse = hz.getSimulation(sim_reverse_idx)
	sim_direct = hz.getSimulation(sim_direct_idx)

	print "Updating Sources"
	with hz.Timer():
		# Get Phased sources
		phased_sources = phasors.create_phasors(sim_reverse, amplitude=1)
		
		# Remove sources before adding them again
		to_del = []
		for s in sim_direct.AllSettings:
			if isinstance(s,acoustic.SourceSettings):
				to_del.append(s)
		for d in to_del:
			sim_direct.Remove(d)
		
		# Add sources
		for source in phased_sources:
			source_settings = acoustic.SourceSettings()
			source_settings.Amplitude = source['amplitude']
			source_settings.Phase = source['phase']
			sim_direct.Add(source_settings, source['model'])

def visualizeMaxSlice(sim_idx, axis=1):
	"""
	Visualize max slice in postpro slice field viewer for pressure field. 
	Set axis to zero to display all 3 directions at a time.
	visualizeMaxSlice(sim_idx, axis=1)
	sim_idx - int / str - index of simulation whose pressure field will be extracted
	axis - int - 0: 3 Slice Field Viewer in all three directions. 1: XY 2: XZ 3: YZ
	"""
	out = hz.getSimulation(sim_idx).Results()
	p = out["Overall Field"]["p(x,y,z,f)"] 
	
	if (axis == 1 or axis == 0): 
		# Create a slice viewer for the pressure
		slice_field_viewer_pressure = analysis.viewers.SliceFieldViewer()
		slice_field_viewer_pressure.SetInputConnection( p )
		slice_field_viewer_pressure.Plane = analysis.viewers.SliceFieldViewer.ePlane.kXY
		slice_field_viewer_pressure.Update(0)
		slice_field_viewer_pressure.GotoMaxSlice()
		slice_field_viewer_pressure.Update(0)
		analysis.RegisterAlgorithm( slice_field_viewer_pressure )
	
	if (axis == 2 or axis == 0): 
		# Create a slice viewer for the pressure
		slice_field_viewer_pressure = analysis.viewers.SliceFieldViewer()
		slice_field_viewer_pressure.SetInputConnection( p )
		slice_field_viewer_pressure.Plane = analysis.viewers.SliceFieldViewer.ePlane.kXZ
		slice_field_viewer_pressure.Update(0)
		slice_field_viewer_pressure.GotoMaxSlice()
		slice_field_viewer_pressure.Update(0)
		analysis.RegisterAlgorithm( slice_field_viewer_pressure )
	
	if (axis == 3 or axis == 0): 
		# Create a slice viewer for the pressure
		slice_field_viewer_pressure = analysis.viewers.SliceFieldViewer()
		slice_field_viewer_pressure.SetInputConnection( p )
		slice_field_viewer_pressure.Plane = analysis.viewers.SliceFieldViewer.ePlane.kYZ
		slice_field_viewer_pressure.Update(0)
		slice_field_viewer_pressure.GotoMaxSlice()
		slice_field_viewer_pressure.Update(0)
		analysis.RegisterAlgorithm( slice_field_viewer_pressure )
	
def visualizeVoxelProperty(idx, property):
	"""
	Visualize Voxel array with data as stored in input file
	visualizeVoxelProperty(idx, property)
	idx - int - index of simulation
	property - string - name (or beginning of name) of property
		'Mass', 'Speed of Sound', 'Attenuation Coefficient', 'Nonlinearity Parameter (B/A)'
	"""
	h5_f = hz.getInputFileName(idx)

	# Get id_map
	id_map = hz.geth5Dataset(h5_f, 'Meshes', 'id_map')
	# Convert to UUID String
	id_map_keys = []
	for id in id_map:
		id_map_keys.append(str(uuid.UUID(bytes=str(bytearray(id)))))
	
	# List Updated Material Properties
	id_map_values = np.zeros(len(id_map_keys))
	with h5py.File(h5_f, 'r') as f:
		for i in f['AllMaterialMaps']:
			for j in f['AllMaterialMaps/' + i]:
				for k in f['AllMaterialMaps/' + i + '/' + j]:
					if f['AllMaterialMaps/' + i + '/' + j + '/' + k + '/_Object'].attrs['PropertyName'].startswith(property) == True:
						id_map_values[id_map_keys.index(j)] = f['AllMaterialMaps/' + i + '/' + j + '/' + k + '/_Object'].attrs['Value']

	voxels = hz.getVoxels(idx)
	axes = hz.getAxes(idx)
	keys = np.arange(len(id_map_keys))
	hz.visualizeArrayTransform(voxels, keys, id_map_values, axes['x'], axes['y'], axes['z'])

	return {'k': id_map_keys, 'v': id_map_values}
	
def getSimulationProperties(idx):
	sim = hz.getSimulation(idx)

	freq_in_Hz = sim.SetupSettings.FrequencyProp.ValueAs(s4l.Unit.Hertz)
	periods = sim.SetupSettings.Periods
	
	print "f:", freq_in_Hz/1000000., "MHz, periods:", periods 
	d_min = 99999999999999999999 #H Fix
	wavelength_min = 99999999999999999999 #H Fix
	for settings in sim.AllSettings:
		if isinstance(settings, acoustic.MaterialSettings):
			if settings.IsReflector == False:
				c_in_ms = settings.SpeedOfSoundProp.ValueAs(s4l.Unit.Velocity)
				wavelength = 1.0 * c_in_ms / freq_in_Hz * 1000 # in mm
				total_distance = wavelength * periods # in mm
				att = settings.AttenuationCoefficientProp.Value
				att_per_mm = att / 1000.
				print settings.Name, "c:", c_in_ms, "m/s, lambda:", wavelength, "mm, d:", total_distance, "mm, a: (1,5,10mm)", \
					round(100*np.exp(-att_per_mm*1),1), round(100*np.exp(-att_per_mm*5),1), round(100*np.exp(-att_per_mm*10),1), "% Pressure"
				d_min = np.min( (d_min, total_distance) )
				wavelength_min = np.min( (wavelength_min, wavelength) )
	print "Min distance covered:", d_min, "mm"
	print "Min wavelength:", wavelength_min, "mm"