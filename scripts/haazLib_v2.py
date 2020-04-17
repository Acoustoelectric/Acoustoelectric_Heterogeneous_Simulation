# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
import s4l_v1 as s4l
import s4l_v1.document as document
import s4l_v1.analysis as analysis
import s4l_v1.model as model
import h5py
import shutil
import numpy as np
import subprocess
import vtk
from vtk.util import numpy_support
import os
import time
import matplotlib.pyplot as plt
import s4l_v1.simulation.acoustic as acoustic
import time
import struct 
import uuid
import scipy.io

def interpn(*args, **kw):
	"""
	Interpolation on N-D. 
	ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
	where the arrays x, y, z, ... define a rectangular grid
	and a.shape == (len(x), len(y), len(z), ...) are the values
	interpolate at xi, yi, zi, ...

	https://github.com/JohannesBuchner/regulargrid/blob/master/regulargrid/interpn.py
	"""
	
	# method = kw.pop('method', 'cubic')
	# method = kw.pop('method', 'nearest')
	method = kw.pop('method', 'linear')
	# method = kw.pop('method', 'slinear')
	if kw:
		raise ValueError("Unknown arguments: " % kw.keys())
	nd = (len(args)-1)//2
	if len(args) != 2*nd+1:
		raise ValueError("Wrong number of arguments")
	q = args[:nd]
	qi = args[nd+1:]
	a = args[nd]
	for j in range(nd):
		# a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
		# a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False, fill_value=np.nan)(qi[j])
		a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False)(qi[j])
		
	return a
	
def interpn_nearest(*args, **kw):
	"""
	Interpolation on N-D. 
	ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
	where the arrays x, y, z, ... define a rectangular grid
	and a.shape == (len(x), len(y), len(z), ...) are the values
	interpolate at xi, yi, zi, ...

	https://github.com/JohannesBuchner/regulargrid/blob/master/regulargrid/interpn.py
	"""
	
	# method = kw.pop('method', 'cubic')
	method = kw.pop('method', 'nearest')
	# method = kw.pop('method', 'linear')
	# method = kw.pop('method', 'slinear')
	if kw:
		raise ValueError("Unknown arguments: " % kw.keys())
	nd = (len(args)-1)//2
	if len(args) != 2*nd+1:
		raise ValueError("Wrong number of arguments")
	q = args[:nd]
	qi = args[nd+1:]
	a = args[nd]
	for j in range(nd):
		# a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
		# a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False, fill_value=np.nan)(qi[j])
		a = interp1d(q[j], a, axis=j, kind=method, bounds_error=False)(qi[j])
		
	return a
	
def getSimulation(idx):
	"""
	Get Simulation from S4L file. 
	out = getSimulation(idx)
	idx - integer / string - simulation index or name
	return - simulation object
	"""
	return document.AllSimulations[idx]


def getInputFileName(s):
	"""
	Get input file from simulation.
	out = getInputFileName(s)
	s - simulation object / integer / string - simulation object or index or name
	return - string - input file
	"""
	if isinstance(s, int):
		return getSimulation(s).raw.InputFileName(0)
	elif isinstance(s, str):
		return getSimulation(s).raw.InputFileName(0)
	else:
		return s.raw.InputFileName(0)


def getOutputFileName(s):
	"""
	Get input file from simulation. 
	out = getOutputFileName(s)
	s - simulation object / integer / string - simulation object or index or name
	return - string - input file
	"""
	if isinstance(s, int):
		return getSimulation(s).raw.OutputFileName(0)
	elif isinstance(s, str):
		return getSimulation(s).raw.OutputFileName(0)
	else:
		return s.raw.OutputFileName(0)

	
def writeInputFile(idx):
	"""
	Write input file. 
	writeInputFile(idx)
	idx - simulation object / integer / string - simulation object or index or name
	return - boolean - success
	"""
	return getSimulation(idx).raw.WriteInputFile()

	
def getImage(image_name):
	"""
	Get image from analysis. 
	out = getImage(image_name)
	image_name - string - image data name (or part of name)
	return - image data
	"""
	a = analysis.FindAlgorithm(image_name)[0]
	a.Update(0)
	return a

	
def getImageDimensions(image_name):
	"""
	Get image grid from analysis. 
	out = getImageDimensions(image_name)
	image_name - string - image data name (or part of name)
	return - image grid size
	"""
	return getImage(image_name).Outputs[0].Data.Grid.Dimensions

	
def getImageAxes_uniform(image_name, scaling=0.001):
	"""
	Get image axes assuming uniform grid. 
	{'x','y','z'} = getImageAxes_uniform(image_name, scaling=0.001)
	image_name - string - image data name (or part of name)
	scaling - float - scaling factor axes are multiplied by
	return - image data
	"""
	# Assume uniform grid
	image = getImage(image_name)
	dims = getImageDimensions(image_name)
	pt_last = image.Outputs[0].Data.Grid.GetPoint( image.Outputs[0].Data.Grid.NumberOfPoints - 1 )
	pt_first = image.Outputs[0].Data.Grid.GetPoint(0)

	# Apply grid transform
	pt_last = image.Outputs[0].Data.Grid.Transform.Transform( pt_last )
	pt_first = image.Outputs[0].Data.Grid.Transform.Transform( pt_first )

	temp = pt_last-pt_first
	x_step = temp[0] / (dims[0]-1) # Divide by number of voxels
	y_step = temp[1] / (dims[1]-1)
	z_step = temp[2] / (dims[2]-1)
	axis_x = np.linspace(pt_first[0], pt_last[0], dims[0]) * scaling
	axis_y = np.linspace(pt_first[1], pt_last[1], dims[1]) * scaling
	axis_z = np.linspace(pt_first[2], pt_last[2], dims[2]) * scaling
	return {'x': axis_x, 'y': axis_y, 'z': axis_z}
	
	
def getImageData_flat(image_name):
	"""
	Get flattened image data from analysis. 
	out = getImageData_flat(image_name)
	image_name - string - image data name (or part of name)
	return - 1D array - flattened image data
	"""
	return getImage(image_name).Outputs[0].Data.Field(0)

	
def getImageData_ndarray(image_name):
	"""
	Get image data from analysis reshaped according to its grid. 
	out = getImageData_ndarray(image_name)
	image_name - string - image data name (or part of name)
	return - ndarray - image data
	"""
	data = getImageData_flat(image_name)
	grid = getImageDimensions(image_name)
	if data.size == grid[0]*grid[1]*grid[2]:
		return data.reshape(grid, order='F')
	elif data.size == (grid[0]-1)*(grid[1]-1)*(grid[2]-1):
		return data.reshape(np.array(grid)-1, order='F')

	
def getImageDataAndAxis(image_name, scaling=0.001):
	"""
	Get image data and axis from analysis reshaped according to its grid
	and the grid transform. 
	{'out','x','y','z'} = getImageDataAndAxis(image_name, scaling=0.001)
	image_name - string - image data name (or part of name)
	scaling - float - scaling factor axes are multiplied by
	return - dict - {'out', 'x', 'y', 'z'}
	"""

	# Assume uniform grid
	image = getImage(image_name)
	data = getImageData_ndarray(image_name)
	dims = getImageDimensions(image_name)
	pt_last = image.Outputs[0].Data.Grid.GetPoint( image.Outputs[0].Data.Grid.NumberOfPoints - 1 )
	pt_first = image.Outputs[0].Data.Grid.GetPoint(0)
	# Apply grid transform
	pt_last = image.Outputs[0].Data.Grid.Transform.Transform( pt_last )
	pt_first = image.Outputs[0].Data.Grid.Transform.Transform( pt_first )

	# Check if any axes are reversed / mirrored. Correct axes and data accorrdingly
	temp = pt_last-pt_first
	for i,t in enumerate(temp):
		if t < 0:
			data = np.flip(data,i)
			pt_last[i], pt_first[i] = pt_first[i], pt_last[i]
			
	axis_x = np.linspace(pt_first[0], pt_last[0], dims[0]) * scaling
	axis_y = np.linspace(pt_first[1], pt_last[1], dims[1]) * scaling
	axis_z = np.linspace(pt_first[2], pt_last[2], dims[2]) * scaling

	return {'out': data, 'x': axis_x, 'y': axis_y, 'z': axis_z}
	
	
def geth5Dataset(h5f, group_name, dataset_name):
	"""
	Find and get dataset from h5 file. 
	out = geth5Dataset(h5f, group_name, dataset_name)
	h5f - string - h5 file path and name
	group_name - string - where to initiate search, '/' for root
	dataset_name - string - dataset to be found
	return - numpy array 
	"""
	
	def find_dataset(name):
		""" Find first object with dataset_name anywhere in the name """
		if dataset_name in name:
			return name
	
	with h5py.File(h5f, 'r') as f:
		k = f[group_name].visit(find_dataset)
		return f[group_name + '/' + k][()]


def seth5Dataset(h5f, group_name, dataset_name, data):
	"""
	Set numpy array dataset from h5 file. 
	bool = seth5Dataset(h5f, group_name, dataset_name, data)
	h5f - string - h5 file path and name
	group_name - string - where to initiate search, '/' for root
	dataset_name - string - dataset to be found
	data - numpy array - data to be copied over
	return - bool - success
	"""
	
	def find_dataset(name):
		""" 
		Find first object with dataset_name anywhere in the name 
		"""
		if dataset_name in name:
			return name
	
	with h5py.File(h5f, 'r+') as f:
		k = f[group_name].visit(find_dataset)
		if k == None:
			return False
		f[group_name + '/' + k][()] = data
	return True

	
def setInputFileForVSDebug(idx, f=r"D:\_haza\tmp\tmp.h5", write_input_file=True):
	"""
	Copy simulation's input file to location for Visual Studio debugging
	setInputFileForVSDebug(idx, f=r"D:\_haza\tmp\tmp.h5", write_input_file=True)
	idx - int - index of simulation in S4L
	f - string - h5 file path and name where to copy input file
	write_input_file - bool - select whether to (re)write input file before copying over
	"""
	sim = getSimulation(idx)
	if write_input_file:
		writeInputFile(idx)
	
	if sim:
		orig_input_file = sim.InputFileName(0)
		print "Simulation:", sim.SimulationName
		print "Copying ", orig_input_file
		print "to ", f
		shutil.copyfile(orig_input_file,f)

	
def compare(a,b):
	"""
	Compare two arrays for equality within tolerance
	bool = compare(a,b)
	a - array - array
	b - array - array
	return - boolean - 
	"""
	if (a==b).all():
		print "Identical equality"
		return True
	else:
		rt = 1
		while True:
			if np.allclose(o0, o1, rt):
				rt = rt / 1.01
			else:
				rt = rt*1.001
				if rt > 0.5:
					print "Not equal"
				else:
					print "Equality with relative tolerance of: ", rt
				break
		return False

	
def compareh5Dataset(h5f1, h5f2, dataset_name, group_name='/'):
	"""
	Compare two h5 file's dataset for equality within tolerance
	bool = compareh5Dataset(h5f1, h5f2, dataset_name, group_name='/')
	a - array - array
	b - array - array
	return - boolean - equality within tolerance
	"""
	a = geth5Dataset(h5f1, group_name, dataset_name)
	b = geth5Dataset(h5f2, group_name, dataset_name)
	return compare(a,b)

	
def compareh5Voxels(idx1, idx2):
	"""
	Compare two simulation's voxels for equality
	bool = compareh5Voxels(idx1, idx2)
	idx1 - int - index of simulation 1
	idx2 - int - index of simulation 2
	return - boolean - equality within tolerance
	"""
	h5f1 = getInputFileName(idx1)
	h5f2 = getInputFileName(idx2)
	return compareh5Dataset(h5f1, h5f2, 'voxels', 'Meshes')

	
def compareh5OutputField(idx1, idx2):
	"""
	Compare two simulation's output fields for equality within a tolerance
	bool = compareh5OutputField(idx1, idx2)
	idx1 - int - index of simulation 1
	idx2 - int - index of simulation 2
	return - boolean - equality within tolerance
	"""
	h5f1 = getOutputFileName(idx1)
	h5f2 = getOutputFileName(idx2)
	dataset_name = '/AllFields/p(x,y,z,f)/_Object/Snapshots/0/comp0'
	return compareh5Dataset(h5f1, h5f2, dataset_name, 'FieldGroups')

	
def getVoxels(i):
	"""
	Get Simulation's Voxel array
	out = getVoxels(i)
	i - int / str - index of simulation / name of simulation / input file
	return - np array - voxel array
	"""
	if isinstance(i, str) and i.endswith("put.h5"):
		h5f = i
	else:
		h5f = getInputFileName(i)
	return geth5Dataset(h5f, 'Meshes', 'voxels')

	
def getAxes(i):
	"""
	Get Simulation's x, y, and z axes
	{'x','y','z'} = getAxes(i)
	i - int / str - index of simulation / name of simulation / input file
	return - dict - 'x', 'y', 'z' axes
	"""
	if isinstance(i, str) and i.endswith("put.h5"):
		h5f = i
	else:
		h5f = getInputFileName(i)
	x = geth5Dataset(h5f, 'Meshes', 'axis_x')
	y = geth5Dataset(h5f, 'Meshes', 'axis_y')
	z = geth5Dataset(h5f, 'Meshes', 'axis_z')
	return {'x': x, 'y': y, 'z': z}

	
def openInputFile(idx):
	"""
	Open Input File in default hdf viewer
	openInputFile(idx)
	idx - int - index of simulation
	"""
	h5f = getInputFileName(idx)
	os.startfile(h5f)

	
def openOutputFile(idx):
	"""
	Open Output File in default hdf viewer
	openOutputFile(idx)
	idx - int - index of simulation
	"""
	h5f = getOutputFileName(idx)
	os.startfile(h5f)

	
def visualizeArray(data, x=None,y=None,z=None, unit_name="", unit="", name="", scaling = 1.0, visualize_max = False, scalarrange = None, visualize_isosurface = False, debug = False):
	"""
	Create field in postpro to visualize data, if no axis provided
	then they are automatically assigned based on number of elements in data
	visualizeArray(data, x=None,y=None,z=None,unit_name="", unit="", name="", scaling = 0.001, goToMax = False)
	data - numpy array - 3D numpy array
	x - numpy array - 1D numpy array of x axis coordinates
	y - numpy array - 1D numpy array of y axis coordinates
	z - numpy array - 1D numpy array of z axis coordinates
	unit_name - string - unit name to be displayed
	unit - string - unit to be displayed
	name - string - name of field in postpro
	scaling - float - scaling factor for axes (default 0.001)
	goToMax - boolean - visualize slice field viewer of maximum in array
	visualize_slice
	scalarrange
	visualize_isosurface
	isosurface_value
	debug
	"""
	
	if x is None:
		x = np.arange(data.shape[0]+1)*.001
	if y is None:
		y = np.arange(data.shape[1]+1)*.001
	if z is None:
		z = np.arange(data.shape[2]+1)*.001

	grid = analysis.core.RectilinearGrid()
	grid.XAxis = np.array(x)*scaling
	grid.YAxis = np.array(y)*scaling
	grid.ZAxis = np.array(z)*scaling

	field = analysis.core.DoubleFieldData()
	field.Grid = grid
	if data.size == (len(x) * len(y) * len(z)):
		field.ValueLocation = analysis.core.eValueLocation.kNode
	elif data.size == (len(x)-1) * (len(y)-1) * (len(z)-1):
		field.ValueLocation = analysis.core.eValueLocation.kCellCenter
	else:
		print "ERROR: Grid and Data don't match"
		return False
	field.NumberOfComponents = 1
	field.NumberOfSnapshots = 1
	field.Quantity.Unit = s4l.Unit(unit)
	field.Quantity.Name = unit_name

	# Note: memory layout is such that x is fastest, z slowest dimension
	values = data.ravel('F')	
	values = values.astype(np.float64)
	field.SetField( 0, values )
	assert field.Check()

	producer = analysis.core.TrivialProducer()
	if name != "":
		producer.Description = name
	producer.SetDataObject(field)

	if visualize_max:
		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.XY
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)

		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.YZ
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)
		
		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.XZ
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)
	
	if visualize_isosurface:
		iso_surface_viewer = analysis.viewers.IsoSurfaceViewer()
		iso_surface_viewer.Inputs[0].Connect( producer.Outputs[0] )
		iso_surface_viewer.Data.Mode = iso_surface_viewer.Data.Mode.enum.QuantityRealModulus #H CHECK - maybe should just use default (delete this line)
		iso_surface_viewer.Visualization.ScalarBarVisible = False
		iso_surface_viewer.UpdateAttributes()
		iso_surface_viewer.Update(0)
		document.AllAlgorithms.Add(iso_surface_viewer)

	document.AllAlgorithms.Add(producer)
	
	return producer

def visualizeComplexArray(data, x=None,y=None,z=None, unit_name="", unit="", name="", scaling = 1.0, visualize_max = False, scalarrange = None, visualize_isosurface = False, debug = False):
	"""
	Create field in postpro to visualize data, if no axis provided
	then they are automatically assigned based on number of elements in data
	visualizeArray(data, x=None,y=None,z=None,unit_name="", unit="", name="", scaling = 0.001, goToMax = False)
	data - numpy array - 3D numpy array
	x - numpy array - 1D numpy array of x axis coordinates
	y - numpy array - 1D numpy array of y axis coordinates
	z - numpy array - 1D numpy array of z axis coordinates
	unit_name - string - unit name to be displayed
	unit - string - unit to be displayed
	name - string - name of field in postpro
	scaling - float - scaling factor for axes (default 0.001)
	goToMax - boolean - visualize slice field viewer of maximum in array
	visualize_slice
	scalarrange
	visualize_isosurface
	isosurface_value
	debug
	"""
	
	if x is None:
		x = np.arange(data.shape[0]+1)*.001
	if y is None:
		y = np.arange(data.shape[1]+1)*.001
	if z is None:
		z = np.arange(data.shape[2]+1)*.001

	grid = analysis.core.RectilinearGrid()
	grid.XAxis = np.array(x)*scaling
	grid.YAxis = np.array(y)*scaling
	grid.ZAxis = np.array(z)*scaling

	field = analysis.core.ComplexDoubleFieldData()
	field.Grid = grid
	if data.size == (len(x) * len(y) * len(z)):
		field.ValueLocation = analysis.core.eValueLocation.kNode
	elif data.size == (len(x)-1) * (len(y)-1) * (len(z)-1):
		field.ValueLocation = analysis.core.eValueLocation.kCellCenter
	else:
		print "ERROR: Grid and Data don't match"
		return False
	field.NumberOfComponents = 1
	field.NumberOfSnapshots = 1
	field.Quantity.Unit = s4l.Unit(unit)
	field.Quantity.Name = unit_name

	# Note: memory layout is such that x is fastest, z slowest dimension
	values = data.ravel('F')	
	#values = values.astype(np.complex64)
	field.SetField( 0, values )
	assert field.Check()

	producer = analysis.core.TrivialProducer()
	if name != "":
		producer.Description = name
	producer.SetDataObject(field)

	if visualize_max:
		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.XY
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)

		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.YZ
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)
		
		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		sfv.Slice.Plane = sfv.Slice.Plane.enum.XZ
		sfv.Update(0)
		sfv.GotoMaxSlice()
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)
	
	if visualize_isosurface:
		iso_surface_viewer = analysis.viewers.IsoSurfaceViewer()
		iso_surface_viewer.Inputs[0].Connect( producer.Outputs[0] )
		iso_surface_viewer.Data.Mode = iso_surface_viewer.Data.Mode.enum.QuantityRealModulus #H CHECK - maybe should just use default (delete this line)
		iso_surface_viewer.Visualization.ScalarBarVisible = False
		iso_surface_viewer.UpdateAttributes()
		iso_surface_viewer.Update(0)
		document.AllAlgorithms.Add(iso_surface_viewer)

	document.AllAlgorithms.Add(producer)
	
	return producer

def visualize2DArray(data, x=None,y=None,z=None, unit_name="", unit="", name="", scaling = 1., visualize_slice = False, scalarrange = None, visualize_isosurface = False, isosurface_value = None, debug = False):
	"""
	Create field in postpro to visualize data, if no axis provided
	then they are automatically assigned based on number of elements in data
	visualize2DArray(data, x=None,y=None,z=None, unit_name="", unit="", name="", scaling = 0.001, visualize = False, scalarrange = None)
	data - numpy array - 3D numpy array
	x - numpy array - 1D numpy array of x axis coordinates
	y - numpy array - 1D numpy array of y axis coordinates
	z - numpy array - 1D numpy array of z axis coordinates
	unit_name - string - unit name to be displayed
	unit - string - unit to be displayed
	name - string - name of field in postpro
	scaling - float - scaling factor for axes (default 0.001)
	visualize - bool - automatically extract field
	"""
	
	from numpy import newaxis

	# Deal with dimension of data
	if data.ndim == 2:
		data = data[:,:,newaxis]
	elif data.ndim < 2 or data.ndim > 3:
		print "Data Dimensions Error"
		return
		
	# Deal with scalar axis by turning into array
	if np.isscalar(x):
		x = [x]
	elif np.isscalar(y):
		y = [y]
	elif np.isscalar(z):
		z = [z]
	
	#If axes are not set, then make axes
	if x is None:
		x = np.arange(data.shape[0]+1)*.001 
	if y is None:
		y = np.arange(data.shape[1]+1)*.001
	if z is None:
		z = np.arange(data.shape[2]+1)*.001

	# Deal with monotonically decreasing axes
	if np.all(np.diff(x) < 0) and np.size(x) > 1:
		x = x[::-1]
		data = data[::-1]
		if debug == True:
			print "Warning: Monotonically decreasing x axes"
	if np.all(np.diff(y) < 0) and np.size(y) > 1:
		y = y[::-1]
		data = data[:,::-1]
		if debug == True:
			print "Warning: Monotonically decreasing y axes"
	if np.all(np.diff(z) < 0) and np.size(z) > 1:
		z = z[::-1]
		data = data[:,:,::-1]
		if debug == True:
			print "Warning: Monotonically decreasing z axes"
	
	grid = analysis.core.RectilinearGrid()
	grid.XAxis = np.array(x)*scaling
	grid.YAxis = np.array(y)*scaling
	grid.ZAxis = np.array(z)*scaling

	field = analysis.core.DoubleFieldData()
	field.Grid = grid
	if data.size == (len(x) * len(y) * len(z)):
		field.ValueLocation = analysis.core.eValueLocation.kNode
	elif data.size == (len(x)-1) * (len(y)-1) * (len(z)-1):
		field.ValueLocation = analysis.core.eValueLocation.kCellCenter
	else:
		print "ERROR: Grid and Data don't match"
		return
	field.NumberOfComponents = 1
	field.NumberOfSnapshots = 1
	field.Quantity.Unit = s4l.Unit(unit)
	field.Quantity.Name = unit_name

	# Note: memory layout is such that x is fastest, z slowest dimension
	values = data.ravel('F')	
	values = values.astype(np.float64)
	field.SetField( 0, values )
	assert field.Check()

	producer = analysis.core.TrivialProducer()
	if name != "":
		producer.Description = name
	producer.SetDataObject(field)
	
	# Adding a SliceSurfaceViewer
	if visualize_slice:
		sfv = analysis.viewers.SliceFieldViewer()
		sfv.Inputs[0].Connect( producer.Outputs[0] )
		if len(x) == 1:
			sfv.Slice.Plane = sfv.Slice.Plane.enum.YZ
		elif len(y) == 1:
			sfv.Slice.Plane = sfv.Slice.Plane.enum.XZ
		elif len(z) == 1:
			sfv.Slice.Plane = sfv.Slice.Plane.enum.XY
		sfv.Visualization.ScalarBarVisible = False
		if scalarrange != None:
			sfv.ScalarRange = scalarrange
		sfv.Update(0)
		document.AllAlgorithms.Add(sfv)
		
	# Adding a IsoSurfaceViewer
	if visualize_isosurface: 
		iso_surface_viewer = analysis.viewers.IsoSurfaceViewer()
		iso_surface_viewer.Inputs[0].Connect( producer.Outputs[0] )
		if isosurface_value is not None:
			iso_surface_viewer.IsoValues = (isosurface_value,)
		iso_surface_viewer.Data.Mode = iso_surface_viewer.Data.Mode.enum.QuantityRealModulus #H CHECK - maybe should just use default (delete this line)
		iso_surface_viewer.Visualization.ScalarBarVisible = False
		iso_surface_viewer.UpdateAttributes()
		iso_surface_viewer.Update()
		if scalarrange != None:
			iso_surface_viewer.ScalarRange = scalarrange
		iso_surface_viewer.UpdateAttributes()
		iso_surface_viewer.Update()
		document.AllAlgorithms.Add(iso_surface_viewer)

	document.AllAlgorithms.Add(producer)

	return producer

	
def visualizeVoxels(idx, mat_idxs=None):
	"""
	Create field in postpro to visualize voxel area and data stored in voxel array
	visualizeVoxels(idx, mat_idxs=None)
	idx - int - index of simulation
	mat_idxs - list - list of material indices to visualize
	"""
	out = getVoxels(idx)
	if mat_idxs != None:
		out = mapKeystoValues(out, mat_idxs, np.ones_like(mat_idxs))
	axes = getAxes(idx)
	visualizeArray(out, axes['x'],axes['y'],axes['z'])
	
def mapFunctiontoArrayElements(data, f):
	"""
	Map a given function to all array elements 
	out = mapFunctiontoArrayElements(data, f)
	data - numpy array - 3D numpy array
	f - function - function to be applied element wise to array
	return - array - array with transformed elements
	"""	
	f = np.vectorize(f)
	return f(data)
	
def runSolver(idx, bin_path=r'C:\Users\montanaro\Documents\SourceCode\Solvers_build\_bin\Debug/'):
	"""
	Map a given function to all array elements 
	runSolver(idx, mode='Debug', bin_path=r'F:\Users\ioannis\Documents\SourceCode\Solvers_build\_bin\Debug/')
	idx - int - index of simulation
	bin_path - str - path where iSolve is located
	"""	
	prog_name = r'iSolve.exe'
	h5f = getInputFileName(idx)
	subprocess.Popen(['cmd',  '/k', bin_path+prog_name, h5f])	
	
def toVTI(imgArray,name,spacing = [0.1,0.1,0.1], origin = [0,0,0], sigma = 0, imp=True):
	""" 
	Convert an image array into VTI data and import it into S4L
	writeVTI(imgArray,name,spacing = [0.1,0.1,0.1], origin = [0,0,0], sigma = 0)
	"""
	dim = imgArray.shape
	
	image = vtk.vtkImageData()
	#image.SetScalarTypeToUnsignedInt()
	image.SetDimensions(dim[0],dim[1],dim[2])
	image.SetOrigin(origin[0],origin[1],origin[2])
	image.SetSpacing(spacing[0],spacing[1],spacing[2])


	array = np.ravel(imgArray,'F')#.ravel()
	vtkarray = numpy_support.numpy_to_vtk(array,True,vtk.VTK_FLOAT)
	vtkarray.SetName("ImageScalars")
	#image.SetScalarTypeToDouble()
	image.GetPointData().SetScalars(vtkarray)

	writer = vtk.vtkXMLImageDataWriter()
	
	if sigma != 0:
		smoother = vtk.vtkImageGaussianSmooth()
		smoother.SetInputData(image)
		smoother.SetStandardDeviation(1)
		smoother.Update()
		writer.SetInputConnection(smoother.GetOutputPort(0))
	else:
		writer.SetInputData(image)

	writer.SetFileName(name)
	writer.Update()

	if imp==True:
		XCoreModeling.Import(name)
		print "Update this function before it becomes deprecated" 

def maskArray(data, mask, id):
	"""
	Takes two arrays, one of values, one of indices, and returns data array 
	filtered for given index id, with 'np.nan' outside the masked region. 
	out = maskArray(data, mask, id)
	data - numpy array - 3D numpy array of data values
	mask - numpy array - 3D numpy array of indices
	id - int - index from indices to be masked
	return - numpy array - filtered data array values 
	"""
	# out = np.zeros_like(data)
	out = np.full_like(data, np.nan)
	np.copyto(out, data, where= mask == id)
	return out

def mapKeystoValues(data, keys, values):
	"""
	Takes an array filled with keys and converts them to their corresponding array of values
	mapKeystoValues(data, keys, values)
	data - numpy array - 3D numpy array of indices
	keys - 1D numpy array - list of indices
	values - 1D numpy array - list of values 
	return - numpy array - transformed data array
	"""
	# out = np.zeros_like(data)
	out = np.full_like(data, np.nan, dtype=np.double) #H
	# for i in np.arange(len(keys)):
	for i in np.arange(len(keys)):
		k = keys[i]
		v = values[i]
		out[data==k] = v
	return out

def visualizeArrayTransform(data, keys, values, x=None, y=None, z=None, unit_name="", unit="", name=""):
	"""
	Create field in postpro to visualize data transformed according to a key value mapping
	visualizeArrayTransform(data, keys, values, x=None, y=None, z=None, unit_name="", unit="", name="")
	data - numpy array - 3D numpy array of indices / keys 
	keys - list - indeces / keys 
	values - list - values indices / keys will be replaced by
	x - numpy array - 1D numpy array of x axis coordinates
	y - numpy array - 1D numpy array of y axis coordinates
	z - numpy array - 1D numpy array of z axis coordinates
	unit_name - string - unit name to be displayed
	unit - string - unit to be displayed
	name - string - field name in postpro
	"""
	out = mapKeystoValues(data, keys, values)
	visualizeArray(out,x,y,z, unit_name, unit, name)
	
def doVoxelandWrite(idx):
	"""
	Create voxels and write input file
	bool = doVoxelandWrite(idx)
	idx - int - index of simulation
	"""
	sim = getSimulation(idx)
	sim.CreateVoxels()
	return writeInputFile(sim)

def getModelEntity(name):
	"""
	Get first model entity matching name 
	out = getModelEntity(name)
	name - str - model entity name
	"""
	entities = model.AllEntities()
	return entities[name]

def getModelEntityInGroup(group_name, name):
	"""
	Get first model entity matching name in group 
	out = getModelEntityInGroup(group_name, name)
	group_name - str - group name
	name - str - model entity name
	"""
	entities = model.AllEntities()[group_name]
	for e in entities.Entities:
		if e.Name == name:
			out = e
	return out

class Timer(object):
	"""
	Timer class to print time elapsed, to use:
		with Timer('foo_stuff'):
			# code here
	""" 
	def __init__(self, name=None):
		self.name = name

	def __enter__(self):
		self.tstart = time.time()

	def __exit__(self, type, value, traceback):
		if self.name:
			print '[%s]' % self.name,
		print 'Time Elapsed: %s' % (time.time() - self.tstart)

def getBoundingBox(element_name, visualize=False):
	"""
	out = getBoundingBox(element_name, visualize=False)
	"""
	print "Update this function before it becomes deprecated" 
	e = getModelEntity(element_name)
	bb = XCoreModeling.GetBoundingBox([e])
	if visualize:
		wire_block = XCoreModeling.CreateWireBlock(bb[0], bb[1])
		wire_block.Name = element_name + " bounding box"
	return bb

	
def visualizeVectorArray(data, x=None,y=None,z=None, unit_name="", unit=""):
	"""
	Create field in postpro to visualize vector data, if no axis provided
	then they are automatically assigned based on number of elements in data
	visualizeArray(data, x=None,y=None,z=None,unit_name="", unit="")
	data - numpy array - 3D numpy array
	x - numpy array - 1D numpy array of x axis coordinates
	y - numpy array - 1D numpy array of y axis coordinates
	z - numpy array - 1D numpy array of z axis coordinates
	unit_name - string - unit name to be displayed
	unit - string - unit to be displayed
	"""
	import numpy as np
	
	scaling = 0.001	
	if x is None:
		x = np.arange(data.shape[0]+1)*scaling 
	if y is None:
		y = np.arange(data.shape[1]+1)*scaling 
	if z is None:
		z = np.arange(data.shape[2]+1)*scaling 

	grid = analysis.core.RectilinearGrid()
	grid.XAxis = x
	grid.YAxis = y
	grid.ZAxis = z

	field = analysis.core.DoubleFieldData()
	field.Grid = grid
	if data.size == (len(x) * len(y) * len(z) *3):
		field.ValueLocation = analysis.core.eValueLocation.kNode
	elif data.size == (len(x)-1) * (len(y)-1) * (len(z)-1) * 3:
		field.ValueLocation = analysis.core.eValueLocation.kCellCenter
	else:
		print "ERROR: Grid and Data don't match"
		return False
	field.NumberOfComponents = 3
	field.NumberOfSnapshots = 1
	field.Quantity.Unit = s4l.Unit(unit)
	field.Quantity.Name = unit_name

	# Note: memory layout is such that x is fastest, z slowest dimension
	# values = data.ravel('F')
	values = data.reshape(-1, data.shape[-1], order='F')	
	values = values.astype(np.float64)
	field.SetField( 0, values )
	assert field.Check()

	producer = analysis.core.TrivialProducer()
	producer.Description = "Vector Data Name " + unit_name
	producer.SetDataObject(field)

	document.AllAlgorithms.Add(producer)

	
def visualizeVectorList(points,vectors, unit_name="", unit="", name="", visualize=False):
	"""
	visualizeVectorList(points,vectors, unit_name="", unit="", name="", visualize=False)
	"""

	# Visualize Active Neurons
	# Create target grid for where you need interpolated values
	target_grid = analysis.core.PointSet()
	target_grid.Unit = s4l.model.LengthUnits()

	for point in points:
		target_grid.InsertNextPoint(1000*Vec3(point))

	field = analysis.core.FloatFieldData()
	field.Grid = target_grid
	field.Allocate(1,target_grid.NumberOfPoints, 3)
	if unit_name != "":
		field.Quantity.Name = unit_name
	if unit != "":
		field.Quantity.Unit = s4l.Unit(unit)

	# set values of field
	data = np.array(vectors).astype('f')
	field.SetField(0,data)

	source = analysis.core.TrivialProducer()
	source.AddDataObject(field)
	if name != "":
		source.Description = name
	
	if visualize:
		vv = analysis.viewers.VectorFieldViewer() # Create a vector viewer
		vv.Inputs[0].Connect( source.Outputs[0] )
		vv.ArrowSize = 0.1 # Change arrow size
		vv.Update(0) # Update the vector viewer
		document.AllAlgorithms.Add(vv)

	document.AllAlgorithms.Add(source)

	
def visualizePoint(vec3_coord, scaling = 1000, name=""):
	"""
	visualizePoint(vec3_coord, scaling = 1000, name="")
	"""
	vec3_coord[0] = vec3_coord[0]*scaling
	vec3_coord[1] = vec3_coord[1]*scaling
	vec3_coord[2] = vec3_coord[2]*scaling
	pt = s4l.model.CreatePoint(vec3_coord)
	if name != "":
		pt.Name = name
	return pt
	
def getImageStatistics(image_name, image_label = None, scaling = 0.001):
	"""
	getImageStatistics(image_name, image_label = None, scaling = 0.001)
	"""
	print image_name
	data = getImageDataAndAxis(image_name, scaling)
	out = data['out']
	printStatistics(out)
	
	if image_label is not None:
		data_label = getImageDataAndAxis(image_label, scaling)['out']
		if out.shape == data_label.shape:
			for lab in np.unique(data_label):
				print "Label " + str(lab)
				masked_data = maskArray(out, data_label, lab)
				printStatistics(masked_data)
		else:
			print "Image Data and Label don't match"

def getStatistics(data):
	"""
	[min_val, max_val, median_val, mean_val, std_val, num_val] = getStatistics(data)
	"""
	max_val = np.nanmax(data)
	min_val = np.nanmin(data)
	median_val = np.nanmedian(data)
	mean_val = np.nanmean(data)
	std_val = np.nanstd(data)
	num_val = np.isnan(data)[np.isnan(data) == False].size
	return [min_val, max_val, median_val, mean_val, std_val, num_val]
	
def printStatistics(data):
	"""
	[min_val, max_val, median_val, mean_val, std_val, num_val] = printStatistics(data)
	"""
	stats = getStatistics(data)
	template1 = "{0:10} {1:10} {2:10} {3:10} {4:10} {5:10}" # column widths: 10, 10, 10, 10, 10
	template2 = "{0:10.2f} {1:10.2f} {2:10.2f} {3:10.2f} {4:10.2f} {5:10}" # column widths: 10, 10, 10, 10, 10
	print template1.format("Min", "Max", "Median", "Mean", "Std", "N")
	print template2.format(stats[0], stats[1], stats[2], stats[3], stats[4], stats[5])
	return stats
	
def visualizeImageHistogram(image_name, image_label = None, scaling = 0.001, bins = None):
	"""
	visualizeImageHistogram(image_name, image_label = None, scaling = 0.001, bins = None)
	"""
	print image_name
	data = getImageDataAndAxis(image_name, scaling)
	out = data['out']
	visualizeHistogram(out, bins, False)
	
	if image_label is not None:
		data_label = getImageDataAndAxis(image_label, scaling)['out']
		if out.shape == data_label.shape:
			for lab in np.unique(data_label):
				print "Label " + str(lab)
				masked_data = maskArray(out, data_label, lab)
				visualizeHistogram(masked_data, show = False)
		else:
			print "Image Data and Label don't match"			
	plt.show()
			
def visualizeHistogram(data, bins = None, show = True):
	"""
	visualizeHistogram(data, bins = None, show = True)
	"""
	if bins is None:
		N = np.isnan(data)[np.isnan(data) == False].size
		bins = N/1000
		if bins<10:
			bins = 10
	# hist, x = np.histogram(data[~np.isnan(data)], bins)
	hist, x = np.histogram(data[~np.isnan(data)], bins, density = True)
	center = (x[:-1] + x[1:]) / 2
	
	plot(center,hist, name = "Histogram")

	# plt.bar(center, hist, align='center')
	# plt.plot(center, hist)
	# if show == True:
		# plt.show()
	return {'hist' : hist, 'centers' : center, 'x' : x}
	
def getNumberOfSimulations():
	"""
	out = getNumberOfSimulations()
	"""
	return len(s4l.document.AllSimulations)
	
def getResult_field(sim_idx, field_name, sensor_name = "Overall Field"):
	"""
	{'out','x','y','z'} = getResult_field(sim_idx, field_name, sensor_name = "Overall Field")
	"""
	simulation = getSimulation(sim_idx)
	simulation_extractor = simulation.Results()
	sensor_extractor = simulation_extractor[sensor_name]
	output = sensor_extractor.Outputs[field_name]
	output.Update()
	field_flat = output.Data.Field(0)
	dims = np.array(output.Data.Grid.Dimensions)
	if np.prod(dims) != len(field_flat):
		dims = dims-1
	field = field_flat.reshape(dims, order='F')
	
	xaxis = output.Data.Grid.XAxis
	yaxis = output.Data.Grid.YAxis
	zaxis = output.Data.Grid.ZAxis
	return {'out' : field, 'x' : xaxis, 'y' : yaxis, 'z' : zaxis}

def copyToClipboard(string1):
	"""
	copyToClipboard(string1)
	"""
	from Tkinter import Tk
	r = Tk()
	r.withdraw()
	r.clipboard_clear()
	r.clipboard_append(string1)
	r.update() # now it stays on the clipboard after the window is closed
	r.destroy()
	
def savenpArrayToMatlab(data, f_name = r'tmp_matlab', f_path=r'C:\Users\montanaro\Desktop/'):
	scipy.io.savemat(f_path + f_name, dict(data=data))
	
def plot(x,y,xname='x',xunit='m',yname='y',yunit='', name = ""):
	import XPostProcessor as xp
	y=np.asarray(y).astype(np.float32)

	xy_data = xp.FloatXYData()
	xy_data.Allocate(y.size,[1])
	xy_data.Axis = x
	xy_data.AxisQuantity.Name=xname
	xy_data.AxisQuantity.Unit=xp.Unit(xunit)

	# Set attributes
	xy_data.Quantity.Name = yname
	xy_data.Quantity.Unit = xp.Unit(yunit)

	# Accessing component
	xy_data.SetComponent(0,y)

	# Create simple pipeline with 2D plot
	source = xp.TrivialProducer()
	source.SetDataObject(xy_data)
	source.Description=name

	import XPostProcessorUI as xpui
	plot = xpui.PlotViewer(1)
	plot.SetInputConnection(source.GetOutputPort(0))

	# Add to GUI
	reg = xp.AlgorithmRegistry()
	reg.AddAlgorithm(source)
	reg.AddAlgorithm(plot)

	return source

def getHistogramPeak(data, start_from = None, bins = None):
	out = visualizeHistogram(data, bins = bins, show = False)
	hist = out['hist']
	centers = out['centers']
	if start_from is not None:
		filtered_idx = np.where(centers >= start_from)
		hist = hist[filtered_idx]
		centers = centers[filtered_idx]
	if hist.size == 0:
		peak = np.nan
	else:
		peak = centers[np.argmax(hist)]
	return peak
	
def getImageHistogramPeak(image_name, image_label = None, start_from = None, bins = None, scaling = 0.001):
	"""
	visualizeImageHistogram(image_name, image_label = None, scaling = 0.001, bins = None)
	"""
	peaks = np.array([])
	
	print image_name + " peaks"
	print "Overall",
	data = getImageDataAndAxis(image_name, scaling)
	out = data['out']
	
	peaks = np.append(peaks, getHistogramPeak(out, start_from, bins))
	print peaks[-1]
	
	if image_label is not None:
		data_label = getImageDataAndAxis(image_label, scaling)['out']
		if out.shape == data_label.shape:
			for lab in np.unique(data_label):
				print "Label " + str(lab),
				masked_data = maskArray(out, data_label, lab)
				peaks = np.append(peaks, getHistogramPeak(masked_data, start_from, bins))
				print peaks[-1]
		else:
			print "Image Data and Label don't match"			
	return peaks
	
def getImageLabelAndHistogramPeak(image_name, image_label = None, start_from = None, bins = None, scaling = 0.001):
	"""
	visualizeImageHistogram(image_name, image_label = None, scaling = 0.001, bins = None)
	TODO
	-999 assigned to 'overall' 
	"""

	# Create empty array to hold label indices and their peaks
	peaks = np.array([]).reshape(0,2)
	
	print image_name + " peaks"
	print "Overall",
	data = getImageDataAndAxis(image_name, scaling)
	out = data['out']
	
	# Arbitrarily assign the label -999 to the whole CT and its peak 
	peaks = np.vstack( (peaks, [-999,getHistogramPeak(out, start_from, bins)]) )
	print peaks[-1,1]
	
	if image_label is not None:
		data_label = getImageDataAndAxis(image_label, scaling)['out']
		if out.shape == data_label.shape:
			for lab in np.unique(data_label):
				print "Label " + str(lab),
				masked_data = maskArray(out, data_label, lab)
				peaks = np.vstack( (peaks, [lab, getHistogramPeak(masked_data, start_from, bins)]) )
				print peaks[-1,1]
		else:
			print "Image Data and Label don't match"			
	return peaks