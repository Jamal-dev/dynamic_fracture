#!/usr/bin/env python
# coding: utf-8

import numpy as np
#import matplotlib.pyplot as plt
import meshio
import glob

# folder path
folder_path = '*.vtk'
# query point data
# exclude regions which are not needed
width = 10
height = 10
num_points = 100
grid_x, grid_y = np.mgrid[0:width:num_points*1j, 0:height:num_points*1j]
print(grid_x.shape,grid_x.min(),grid_x.max())
print(grid_y.shape,grid_y.min(),grid_y.max())




# read all vtk files in the folder
# specify the path to the folder
files = glob.glob(folder_path)
files.sort()
# print(files)


# In[17]:


def path_crack(phi):
    # left to right
    x0=0.5 * width
    x = x0
    # change delta_x to change the resolution in x
    delta_x = 1e-3
    # change 100 num to change the resolution in y
    y = np.linspace(0, height, 100)
    pos_x = []
    pos_y = []
    phi_values = []
    i = 0
    while (1==1):
        x +=  delta_x
        if   x > width:
            break
        phi_list = phi(x,y)
        min_phi = float(min(phi_list))
        if  min_phi > 0.3 or min_phi < 0:
            continue
        index_y = np.argmin(phi_list)
        # print(x, y[index_y], phi_list[index_y])
        pos_x.append(x)
        pos_y.append(float(y[index_y]))
        phi_values.append(min_phi)
    # pos_x = np.array(pos_x)
    # pos_y = np.array(pos_y)
    # phi_values = np.array(phi_values)
    # print(pos_x.shape, pos_y.shape, phi_values.shape)
    # print(f'max(phi_values) = {max(phi_values)}')
    return pos_x, pos_y, phi_values


# In[23]:


from scipy.interpolate import NearestNDInterpolator
import json
# import yaml
crack_path_data ={'time_steps':[], 'pos_x':{}, 'pos_y':{}, 'phi_values':{}}
for timestep, filename in enumerate(files):
    data = meshio.read(filename)
    phi = NearestNDInterpolator(data.points[:,0:2], data.point_data['phase_field'])
    pos_x, pos_y, phi_values = path_crack(phi)
    # phi_values is not empty
    if len(phi_values) > 0:
        crack_path_data['time_steps'].append(timestep)
        crack_path_data['pos_x'][timestep] = pos_x
        crack_path_data['pos_y'][timestep] = pos_y
        crack_path_data['phi_values'][timestep] = phi_values
with open('crack_path_data.json', 'w') as outfile:
    json.dump(crack_path_data, outfile)





