import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def extract_elevations(filename):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    
    structured_grid = reader.GetOutput()
    points = structured_grid.GetPoints()
    dimensions = structured_grid.GetDimensions()
    
    locs = vtk_to_numpy(points.GetData())[:,:2]
    x_length = dimensions[0]
    
    surf_elev = locs[:x_length]
    base_elev = locs[len(locs)-x_length:len(locs)]
    
    return surf_elev,base_elev

def vtk_to_xarray(filename):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    
    structured_grid = reader.GetOutput()
    points = structured_grid.GetPoints()
    dimensions = (np.array(structured_grid.GetDimensions()[:2]))[::-1]
    
    locs = vtk_to_numpy(points.GetData())[:,:2]
    x_length = dimensions[1]
    
    xads = xr.Dataset()
    
    point_data = structured_grid.GetPointData()
    num_point_arrays = point_data.GetNumberOfArrays()
    for i in range(num_point_arrays):
        array_name = point_data.GetArrayName(i)
        array = vtk_to_numpy(point_data.GetArray(i))
        if(len(np.shape(array))>1):
            dimensions_temp = np.append(dimensions,np.shape(array)[1])
            dims=('point_z','point_x','direction')
            array = array.reshape(dimensions_temp)
        else:
            dims=('point_z','point_x')
            array = array.reshape(dimensions)
        xada = xr.DataArray(array,dims=dims)
        xads[f'{array_name}'] = xada
    
    cell_data = structured_grid.GetCellData()
    num_cell_arrays = cell_data.GetNumberOfArrays()
    cell_dims = dimensions - 1
    for i in range(num_cell_arrays):
        array_name = cell_data.GetArrayName(i)
        array = vtk_to_numpy(cell_data.GetArray(i))
        array = array.reshape(cell_dims)
        xada = xr.DataArray(array,dims=('cell_z','cell_x'))
        xads[f'{array_name}'] = xada
                             
    return xads

def get_grid_locs_and_dimensions(filename):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update() 
    structured_grid = reader.GetOutput()
    dimensions = structured_grid.GetDimensions()   
    points = structured_grid.GetPoints()
    locs = vtk_to_numpy(points.GetData())[:,:2]*1000
    return locs,dimensions

def extract_front_profile(filename):
    locs,dimensions = get_grid_locs_and_dimensions(filename)
    num_x_cells = dimensions[0]
    front = locs[::num_x_cells]
    return front
    