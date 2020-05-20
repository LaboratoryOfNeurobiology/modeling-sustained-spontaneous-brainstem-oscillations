# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 21:46:02 2020

@author: David
"""

# PRELIMINARIES ===============================================================

# Setting working directory
import os
import sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
# Importing function that carries out interval halving
from IntHalving import IHalving
# Importing numpy library
import numpy as np

# FUNCTIONS WHOSE ROOTS WE SEEK ===============================================

def sphere_p3(x):    # n-dimensonal unit sphere with p=3 metric
    '''Evaluates the function at a particular x list or array of coordinates 
    whose roots give a unit circle with a p=3 metric.
    
    Args.: single dimensional array of numbers, has a length of 2
    Returns: a single float number
    '''
    X = np.array(x)
    f = np.sum(np.abs(X)**(1/3))**3 - 1
    return f

def sphere_p2(x):    # n-dimensonal unit sphere with p=2 metric
    '''Evaluates the function at a particular x list or array of coordinates 
    whose roots give a unit circle with a p=2 metric.
    
    Args.: single dimensional array of numbers, has a length of 2
    Returns: a single float number
    '''
    X = np.array(x)
    f = np.sum(np.abs(X)**(1/2))**2 - 1
    return f

def sphere(x):    # n-dimensonal unit sphere with Eucledian metric
    '''Evaluates the function at a particular x list or array of coordinates 
    whose roots give a unit circle with Eucledian metric.
    
    Args.: single dimensional array of numbers, has a length of 2
    Returns: a single float number
    '''
    X = np.array(x)
    f = np.sum(np.abs(X)**2)**(1/2) - 1
    return f

# PARAMETERS ==================================================================

# Chosen function
fun = sphere_p2
# Name given to iteration (must be a string)
name = fun.__name__
# Range of solution
#sol_range = [[0.1, 1], [-0.1, 1], [0.01, 1]]
sol_range = [[-0.1, 1], [-0.1, 1]]
# Number of steps along the associated coordinates, must be integers
Num_div = [2, 2]
# Number of iteration steps
N_it = 8

# INTERVAL-HALVING ITERATIONS =================================================

# Computing results
sol, It = IHalving(fun, sol_range, Num_div, N_it, name=name, save='end', tol='eps')

# Printing details of iteration
print(It)

# PLOTTING ====================================================================

# Would you like to save figure?
FigSav = True

# marker type
mark = ","
# point size for boundary points
size_pt_br = 0.001
# point size for grid points
size_pt_gr = 0.0001
# color of boundary points
col_pt_br = [1, 0, 0]
# color of grid points
col_pt_gr = [0, 0, 0]
# margin of plot in [%]
marg = 0

# For camera position control of 3D plot:

# azimuthal angle (rotation about vertical axis) in [deg]
azim = -45
# altitude angle (rotation around horizon) in [deg]
alti = 45

import matplotlib.pyplot as plt

if len(Num_div) == 2: # 2D
    
    fig, ax = plt.subplots()
    # Plotting gridpoints for all iteration steps:
    grid = It.get_nodes_scaled(It.get_last_step())
    for i in range(0, np.shape(grid)[0]):
        ax.scatter(grid[i,0], grid[i,1], s=size_pt_gr, marker=mark, c=[col_pt_gr])
    # Plotting boundary points (solution):
    for i in range(0, np.shape(sol)[0]):
        ax.scatter(sol[i,0], sol[i,1], s=size_pt_br, marker=mark, c=[col_pt_br])
    # Customizing plotrange
    sol_range = np.array(sol_range)
    dim_x = sol_range[0, 1] - sol_range[0, 0]
    dim_y = sol_range[1, 1] - sol_range[1, 0]
    plt.xlim(sol_range[0, 0] - dim_x*marg/100, sol_range[0, 1] + dim_x*marg/100)
    plt.ylim(sol_range[1, 0] - dim_y*marg/100, sol_range[1, 1] + dim_y*marg/100)
    # Saving figure:
    if FigSav:
        fig.savefig(fun.__name__ + '_Nx_' + str(Num_div[0]) + '_Ny_' + str(Num_div[1]) + '_Nit_' + str(N_it) + '.pdf', bbox_inches='tight')

elif len(Num_div) == 3: # 3D

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
#    # Plotting gridpoints for all iteration steps:
#    grid = It.get_nodes_scaled(It.get_last_step())
#    for i in range(0, np.shape(grid)[0]):
#        ax.scatter(grid[i,0], grid[i,1], grid[i,2], marker=mark, s=size_pt_gr, c=[col_pt_gr])
    # Plotting boundary points (solution):
    for i in range(0, np.shape(sol)[0]):
        ax.scatter(sol[i,0], sol[i,1], sol[i,2], marker=mark, s=size_pt_br, c=[col_pt_br])
    # Customizing plotrange
    sol_range = np.array(sol_range)
    dim_x = sol_range[0, 1] - sol_range[0, 0]
    dim_y = sol_range[1, 1] - sol_range[1, 0]
    dim_z = sol_range[2, 1] - sol_range[2, 0]
    xlim = (sol_range[0, 0] - dim_x*marg/100, sol_range[0, 1] + dim_x*marg/100)
    ylim = (sol_range[1, 0] - dim_y*marg/100, sol_range[1, 1] + dim_y*marg/100)
    zlim = (sol_range[2, 0] - dim_z*marg/100, sol_range[2, 1] + dim_z*marg/100)    
    ax.set_xlim(left=xlim[0], right=xlim[1])
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    ax.set_zlim(bottom=zlim[0], top=zlim[1])
    
    ax.view_init(elev=alti, azim=azim)
    
    # Saving figure:
    if FigSav:
        fig.savefig(fun.__name__ + '_Nx_' + str(Num_div[0]) + '_Ny_' + str(Num_div[1]) + '_Nz_' + str(Num_div[2]) + '_Nit_' + str(N_it) + '.pdf', bbox_inches='tight')
    
else:    
    print('This scipt cannot plot solution boundaries in ' + str(len(Num_div)) + ' dimensions.')