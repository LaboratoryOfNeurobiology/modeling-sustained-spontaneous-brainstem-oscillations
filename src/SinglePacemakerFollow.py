# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 21:46:02 2020

@author: David
"""

# PRELIMINARIES ===============================================================

# Setting working directory
from __future__ import division
import os
import sys
import numpy as np
import random as rdm
from math import sin, cos, pi
from neuron import h, nrn, gui
import matplotlib.pyplot as plt
import faulthandler
import PN_Modeling as pnm
import pickle as pkl
from IntHalving_follow import IHalving_follow


# FUNCTIONS WHOSE ROOTS WE SEEK ===============================================
def pacemaker_cell_func(LoParameters): #node):
    """ Evaluates the simulation of the pacemaker nucleus
    with the parameters defined in LoParameters and returns
    the frequency of oscillations.

    @LoParameters: [EK, gKp]
    @returns: frequency of oscillations either 1e-15 or freq > 0

    """
    parconfig = "[EK, gKs, gKa]"
    X = np.array(LoParameters)
    freq = None
    # NEURON utilities
    h.cvode_active(1)
    h.finitialize(-65)
    h.celsius = 27
    # Biophysical parameters
    ek = LoParameters[0]
    ena = 50
    # Pacemaker cell soma specific
    ps_EL = -70
    ps_gNa = 1.0
    ps_gK = LoParameters[1]
    ps_gL = 0.0001
    J = 30
    # Pacemaker cell axon specific
    pa_EL = -70
    pa_gNa = 0.5
    pa_gK = LoParameters[2]
    pa_gL = 0.001
    M = 45

    # Duration Parameters
    T_STOP = 100  # (ms)
    # Object storing
    pacemaker_cells = []
    pacemaker_cells.append(pnm.PacemakerCell([ek, ena, ps_EL, ps_gNa, ps_gK, ps_gL, J],
                                             [ek, ena, pa_EL, pa_gNa, pa_gK, pa_gL, M],
                                             0))
    """
    Orient objects in 3D-space with polar coordinates (position, rotation)
    where the center of the coordinate system corresponds to
    the center of the pacemaker nucleus cell network.

    Default neuronal orientation before repositioning and rotation
    y    z         y
    ^  ^>          ^
    | /            | _______
    |/             |(       )_______________________________________
    |------>  (0,0)+(-So>ma-)______________Axon_____________________----> x
    |              |(_______)
    |              |
    V              v
    """
    t_pace = pacemaker_cells[0]
    len_pace = t_pace.give_len("soma") + t_pace.give_len("axon")
    dt_pace = 2 * pi / 1

    # Position for pacemaker cells (position, rotation)
    for pace, cell in enumerate(pacemaker_cells):
        cell.set_position((len_pace + 10.001) * cos(pi + (dt_pace * pace)),
                          (len_pace + 10.001) * sin(pi + (dt_pace * pace)),
                          0)
        cell.rotateZ(pace * (2 * pi / 1))

    # Begin simulation of model
    print(f"Starting simulation {list(X)} on pc={pc.id()}")
    h.tstop = T_STOP
    h.run()

    # Simulation analysis
    single_pace_s = []
    single_pace_a = []
    # File setup
    try:
        #os.mkdir(f"/Users/daniel/Desktop/Development/simulations/{list(X)}")
        os.mkdir(f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/ekgkgk_pacemakercell_sims/{list(X)}")
    except OSError:
        pass
    spikes_file = open(
        f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/ekgkgk_pacemakercell_sims/{list(X)}/{list(X)}_spiketimes.txt",
        "w")
    #spikes_file = open(
    #    f"/Users/daniel/Desktop/Development/simulations/{list(X)}/{list(X)}_spiketimes.txt",
    #    "w")
    # Determine if oscillating spontaneously
    last_s = None
    last_a = None
    for i, cell in enumerate(pacemaker_cells):
        soma_v, axon_v, t_v = cell.give_spikes()
        s_v, a_v = cell.give_vec_halves()
        if (len(soma_v) > 0 and len(axon_v) > 0):
            last_s = soma_v[-1]
            last_a = axon_v[-1]
        spikes_file.write(f"{['volt', i, list(s_v), list(a_v)]}\n")
        if i == 0:
            spikes_file.write(f"{['time', list(t_v)]}\n")
            single_pace_s = list(soma_v)
            single_pace_a = list(axon_v)
        # Store raw cellular spike data for this simulation.
        spikes_file.write(f"{['soma', i, list(soma_v)]}\n")
        spikes_file.write(f"{['axon', i, list(axon_v)]}\n")
    p_s_f = len(single_pace_s) / (T_STOP * 0.001)
    p_a_f = len(single_pace_a) / (T_STOP * 0.001)
    spikes_file.write(f"{['freqs', p_s_f, p_a_f]}\n")
    spikes_file.close()

    if (len(single_pace_s) > 3
        and last_s is not None
        and last_a is not None
        and 100 - last_s < 25
        and 100 - last_a < 25):
        if p_s_f == p_a_f:
            freq = p_s_f
        else:
            freq = p_a_f
    else:
        freq = -1e-15

    print("pacemaker_soma freq = " + str(p_s_f) + " Hz.")
    print("pacemaker_axon freq = " + str(p_a_f) + " Hz.")
    # pnm.cellular_potentials(all_cells, f"{X}={parconfig}", X, p_s_f, r_s_f)
    # pnm.raster(all_cells, True, f"{X}={parconfig}", X, p_s_f, p_a_f, r_s_f, r_a_f)
    #simulations_file = open(f"/Users/daniel/Desktop/Development/simulations/SimulationRawOutput.txt", "a")
    simulations_file = open(
        f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/ekgkgk_pacemakercell_sims/SimulationRawOutput3d.txt", "a")
    simulations_file.write(f"{[list(X), p_s_f, p_a_f]}\n")
    simulations_file.close()
    print(f"Final = {freq}\n")
    print(f"Simulation {list(X)} on pc={pc.id()} is done.\n")
    #return node, freq
    return freq

if __name__ == '__main__':
    h.nrnmpi_init()
    faulthandler.enable()
    pc = h.ParallelContext()
    pc.runworker()
    file_dir = os.path.dirname(__file__)
    sys.path.append(file_dir)
    h.load_file("nrngui.hoc")
    # RE-LOADING DATA =============================================================
    # file from which data will be retrieved
    file_name = '/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/It_pacemaker_cell_func_2020_May_14_16h_11m'

    file_handle = open(file_name + '.pkl', 'rb')
    print('loading pkl file')
    It = pkl.load(file_handle)
    file_handle.close()
    print('load complete')
    # funtion that is used for solution patch following
    fun = pacemaker_cell_func
    # number of iterations
    print('getting last step')
    N_it = It.get_last_step()
    # solution range
    print('getting sol rnage')
    sol_range = np.array(It.get_sol_range())
    # number of divisions
    print('getting n_div')
    Num_div = np.array(It.get_n_div())

    # INTERVAL-HALVING ITERATIONS =================================================

    # Computing results
    print('starting int halving')
    sol, Sol = IHalving_follow(pc, fun, It)

    # Getting all nodes
    gr = Sol.get_nodes_scaled()

    # Printing details of iteration
    print(Sol)

    # PLOTTING ====================================================================
    """
    # Would you like to save figure?
    FigSav = True

    # marker type
    mark = ","
    # point size for boundary points
    size_pt_br = 0.001
    # point size for grid points
    size_pt_gr = 0.1
    # color of boundary points
    col_pt_br = [1, 0, 0]
    # color of grid points
    col_pt_gr = [0, 0, 0]
    # margin of plot in [%]
    marg = 5

    # For camera position control of 3D plot:

    # azimuthal angle (rotation about vertical axis) in [deg]
    azim = -45
    # altitude engle (rotation around horizon) in [deg]
    alti = 45

    import matplotlib.pyplot as plt

    if len(Num_div) == 2: # 2D

        fig, ax = plt.subplots()

        # plotting gridpoints
        for i in range(np.shape(gr)[0]):
            ax.scatter(gr[i, 0], gr[i, 1] , s=size_pt_br, marker=mark, c=[col_pt_gr])

        # plotting points along solution
        for i in range(np.shape(sol)[0]):
            ax.scatter(sol[i, 0], sol[i, 1] , s=size_pt_br, marker=mark, c=[col_pt_br])

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
    #
    #    # plotting gridpoints
    #    for i in range(np.shape(gr)[0]):
    #        ax.scatter(gr[i, 0], gr[i, 1], gr[i, 2] , s=size_pt_br, marker=mark, c=[col_pt_gr])

        # plotting points along solution
        for i in range(np.shape(sol)[0]):
            ax.scatter(sol[i, 0], sol[i, 1], sol[i, 2] , s=size_pt_br, marker=mark, c=[col_pt_br])

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
        """
    h.quit()
