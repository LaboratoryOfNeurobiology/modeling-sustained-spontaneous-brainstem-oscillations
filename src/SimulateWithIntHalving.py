#!/home/zupanc-lab/anaconda3/bin/python
from __future__ import division
import os
import numpy as np
import random as rdm
import sys
from math import sin, cos, pi
import PN_Modeling as pnm
from neuron import h, nrn, gui
from IntHalving import IHalving
import matplotlib.pyplot as plt
import faulthandler


def network_func(LoParameters): #node):
    """ Evaluates the simulation of the pacemaker nucleus
    with the parameters defined in LoParameters and returns
    the frequency of oscillations.

    @LoParameters: [EK, ENa, gNar, gKr, gNap, gKp]
    @returns: frequency of oscillations either 0 or freq > 0

    """
    parconfig = "[EK, ENa, gNar, gKr, gNap, gKp]"
    X = np.array(LoParameters)
    freq = None
    # NEURON utilities
    h.cvode_active(1)
    h.finitialize(-65)
    h.celsius = 27
    # Biophysical parameters
    ek = LoParameters[0]
    ena = LoParameters[1]
    # Pacemaker cell soma specific
    ps_EL = - 70
    ps_gNa = LoParameters[4]
    ps_gK = LoParameters[5]
    ps_gL = 0.0001
    J = 30
    # Pacemaker cell axon specific
    pa_EL = - 70
    pa_gNa = 0.5
    pa_gK = 0.02
    pa_gL = 0.001
    M = 45
    # Relay cell soma specific
    rs_EL = -70
    rs_gNa = LoParameters[2]
    rs_gK = LoParameters[3]
    rs_gL = 0.0003
    K = 60
    # Relay cell axon specific
    ra_EL = -70
    ra_gNa = 0.5
    ra_gK = 0.05
    ra_gL = 0.001
    N = 40
    # Synapse Parameters
    conduct_rng = [0.5, 10]  # nS
    # Topology Parameters
    n_pacemakers = 87
    n_relays = 20
    n_p2p_projections = 6
    n_p2r_projections = 7
    # Duration Parameters
    T_STOP = 100  # (ms)
    # Object storing
    pacemaker_cells = []
    relay_cells = []
    pace_network_graph = pnm.Graph()

    # Build Synapse Adjacency List and initialize cell objects.
    for i in range(n_pacemakers + n_relays):
        # num vertices = num cells in network
        pace_network_graph.add_vertex(i)
        if i < n_pacemakers:  # LOOKING AT A PACEMAKER CELL
            p_p_projections = rdm.sample(range(0, n_pacemakers),
                                         n_p2p_projections)
            p_r_projections = rdm.sample(range(n_pacemakers,
                                               n_relays + n_pacemakers), n_p2r_projections)

            while i in p_p_projections:  # DON'T WANT TO SYNAPSE ONTO SELF
                p_p_projections = rdm.sample(range(0, n_pacemakers),
                                             n_p2p_projections)

            # Now have 2 unique lists of all cells that cell #i synapses to.
            # Add cell i's projections to graph
            all_projections = p_p_projections + p_r_projections
            for proj in all_projections:
                pace_network_graph.add_edge([i, proj])
            pacemaker_cells.append(pnm.PacemakerCell([ek, ena, ps_EL, ps_gNa, ps_gK, ps_gL, J],
                                                     [ek, ena, pa_EL, pa_gNa, pa_gK, pa_gL, M],
                                                     i))

        else:  # LOOKING AT A RELAY CELL
            # Relay cells don't synapse to anything.
            relay_cells.append(pnm.RelayCell([ek, ena, rs_EL, rs_gNa, rs_gK, rs_gL, K],
                                             [ek, ena, ra_EL, ra_gNa, ra_gK, ra_gL, N],
                                             i))
    all_cells = pacemaker_cells + relay_cells
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
    t_relay = relay_cells[0]
    len_pace = t_pace.give_len("soma") + t_pace.give_len("axon")
    len_relay = t_relay.give_len("soma") + t_relay.give_len("axon")
    dt_pace = 2 * pi / n_pacemakers
    dt_relay = 2 * pi / n_relays
    # First for relay cells (position, rotation)
    for relay, cell in enumerate(relay_cells):
        cell.set_position((len_relay + 10) * cos(pi + (dt_relay * relay)),
                          (len_relay + 10) * sin(pi + (dt_relay * relay)),
                          0)
        cell.rotateZ(relay * (2 * pi / n_relays))
    # Second for pacemaker cells (position, rotation)
    for pace, cell in enumerate(pacemaker_cells):
        cell.set_position((len_pace + len_relay + 10.001) * cos(pi + (dt_pace * pace)),
                          (len_pace + len_relay + 10.001) * sin(pi + (dt_pace * pace)),
                          0)
        cell.rotateZ(pace * (2 * pi / n_pacemakers))
    # For dict entries, create a Synapse b/t the key cell & each payload cell.
    for key in pace_network_graph.vertices():
        pre_syn_cell = all_cells[key]
        for ident in pace_network_graph.edges_of_vertex(key):
            post_syn_cell = all_cells[ident]
            pre_syn_cell.add_synapse(post_syn_cell, conduct_rng)

    # Begin simulation of model
    print(f"Starting simulation {list(X)} on pc={pc.id()}")
    h.tstop = T_STOP
    h.run()

    # Simulation analysis
    single_pace_s = []
    single_pace_a = []
    single_relay_s = []
    single_relay_a = []
    x = rdm.randint(0, 86)
    y = rdm.randint(87, 106)
    # File setup
    try:
        #os.mkdir(f"/Users/daniel/Desktop/Development/simulations/{list(X)}")
        os.mkdir(f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/simulations/{list(X)}")
    except OSError:
        pass
    spikes_file = open(
        f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/simulations/{list(X)}/{list(X)}_spiketimes.txt",
        "w")
    #spikes_file = open(
    #    f"/Users/daniel/Desktop/Development/simulations/{list(X)}/{list(X)}_spiketimes.txt",
    #    "w")
    # Determine if oscillating spontaneously
    for i, cell in enumerate(all_cells):
        soma_v, axon_v, t_v = cell.give_spikes()
        s_v, a_v = cell.give_vec_halves()
        spikes_file.write(f"{['volt', i, list(s_v), list(a_v)]}\n")
        if i == 0:
            spikes_file.write(f"{['time', list(t_v)]}\n")
        if i == x:
            single_pace_s = list(soma_v)
            single_pace_a = list(axon_v)
        if i == y:
            single_relay_s = list(soma_v)
            single_relay_a = list(axon_v)
        # Store raw cellular spike data for this simulation.
        spikes_file.write(f"{['soma', i, list(soma_v)]}\n")
        spikes_file.write(f"{['axon', i, list(axon_v)]}\n")
    p_s_f = len(single_pace_s) / (T_STOP * 0.001)
    p_a_f = len(single_pace_a) / (T_STOP * 0.001)
    r_s_f = len(single_relay_s) / (T_STOP * 0.001)
    r_a_f = len(single_relay_a) / (T_STOP * 0.001)
    spikes_file.write(f"{['freqs', p_s_f, p_a_f, r_s_f, r_a_f]}\n")
    spikes_file.close()

    if ((len(single_relay_s) - (len(single_pace_s)) <= 1)
            and (len(single_pace_s) - (len(single_relay_s)) <= 1)
            and (len(single_pace_s) > 3)
            and (len(single_relay_s) > 3)):
        if p_s_f == r_s_f:
            freq = p_s_f
        else:
            freq = r_s_f
    else:
        freq = 1e-15

    print("pacemaker_soma freq = " + str(p_s_f) + " Hz.")
    print("pacemaker_axon freq = " + str(p_a_f) + " Hz.")
    print("relay_soma freq = " + str(r_s_f) + " Hz.")
    print("relay_axon freq = " + str(r_a_f) + " Hz.")
    # pnm.cellular_potentials(all_cells, f"{X}={parconfig}", X, p_s_f, r_s_f)
    # pnm.raster(all_cells, True, f"{X}={parconfig}", X, p_s_f, p_a_f, r_s_f, r_a_f)
    #simulations_file = open(f"/Users/daniel/Desktop/Development/simulations/SimulationRawOutput.txt", "a")
    simulations_file = open(
        f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/simulations/SimulationRawOutput.txt", "a")
    simulations_file.write(f"{[list(X), p_s_f, p_a_f, r_s_f, r_a_f]}\n")
    simulations_file.close()
    print(f"Final = {freq}\n")
    print(f"Simulation {list(X)} on pc={pc.id()} is done.\n")
    #return node, freq
    return freq


def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno) + " on PCID:%d" % (pc.id()))
    return trace


if __name__ == '__main__':
    # FILE STRUCTURE SETUP ======================================================
    h.nrnmpi_init()
    faulthandler.enable()
    #sys.settrace(trace)
    pc = h.ParallelContext()
    pc.runworker()
    file_dir = os.path.dirname(__file__)
    sys.path.append(file_dir)
    # Neuron setup
    h.load_file("nrngui.hoc")
    fun = network_func
    name = fun.__name__
    # Range of solution
    sol_range = [[-150, -50], [50, 150], [0, 1], [0, 0.1], [0, 0.75], [0, 0.09]]
    # Number of steps along the corresponding coordinates
    Num_div = [4, 4, 3, 3, 3, 3]
    # Number of iterations
    N_it = 6
    FigSav = True
    # Generate parameter boundaries with IntervalHalving
    sol, It = IHalving(pc, fun, sol_range, Num_div, N_it, name=name, save='step', tol='eps')
    pc.done()
    dir_path = "/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/simulations"
    pnm.post_sim_analysis(dir_path)
    # print(It)

    # PLOTTING ====================================================================
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
    marg = 20

    if len(Num_div) == 2:  # 2D

        fig, ax = plt.subplots()
        # Plotting gridpoints for all iteration steps:
        grid = It.get_nodes_scaled(It.get_last_step())
        for i in range(0, np.shape(grid)[0]):
            ax.scatter(grid[i, 0], grid[i, 1], s=size_pt_gr, marker=mark, c=[col_pt_gr])
        # Plotting boundary points (solution):
        for i in range(0, np.shape(sol)[0]):
            ax.scatter(sol[i, 0], sol[i, 1], s=size_pt_br, marker=mark, c=[col_pt_br])
        # Customizing plotrange
        sol_range = np.array(sol_range)
        dim_x = sol_range[0, 1] - sol_range[0, 0]
        dim_y = sol_range[1, 1] - sol_range[1, 0]
        plt.xlim(sol_range[0, 0] - dim_x * marg / 100, sol_range[0, 1] + dim_x * marg / 100)
        plt.ylim(sol_range[1, 0] - dim_y * marg / 100, sol_range[1, 1] + dim_y * marg / 100)
        # Saving figure:
        if FigSav:
            fig.savefig(
                fun.__name__ + '_Nx_' + str(Num_div[0]) + '_Ny_' + str(Num_div[1]) + '_Nit_' + str(N_it) + '.pdf',
                bbox_inches='tight')

    elif len(Num_div) == 3:  # 3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Plotting gridpoints for all iteration steps:
        grid = It.get_nodes_scaled(It.get_last_step())
        for i in range(0, np.shape(grid)[0]):
            ax.scatter(grid[i, 0], grid[i, 1], grid[i, 2], marker=mark, s=size_pt_gr, c=[col_pt_gr])
        # Plotting boundary points (solution):
        for i in range(0, np.shape(sol)[0]):
            ax.scatter(sol[i, 0], sol[i, 1], sol[i, 2], marker=mark, s=size_pt_br, c=[col_pt_br])
        # Customizing plotrange
        sol_range = np.array(sol_range)
        dim_x = sol_range[0, 1] - sol_range[0, 0]
        dim_y = sol_range[1, 1] - sol_range[1, 0]
        dim_z = sol_range[2, 1] - sol_range[2, 0]
        xlim = (sol_range[0, 0] - dim_x * marg / 100, sol_range[0, 1] + dim_x * marg / 100)
        ylim = (sol_range[1, 0] - dim_y * marg / 100, sol_range[1, 1] + dim_y * marg / 100)
        zlim = (sol_range[2, 0] - dim_z * marg / 100, sol_range[2, 1] + dim_z * marg / 100)
        ax.set_xlim(left=xlim[0], right=xlim[1])
        ax.set_ylim(bottom=ylim[0], top=ylim[1])
        ax.set_zlim(bottom=zlim[0], top=zlim[1])
        # Saving figure:
        if FigSav:
            fig.savefig(fun.__name__ + '_Nx_' + str(Num_div[0]) + '_Ny_' + str(Num_div[1]) + '_Nz_' + str(
                Num_div[2]) + '_Nit_' + str(N_it) + '.pdf', bbox_inches='tight')

    else:
        print('This scipt cannot plot solution boundaries in ' + str(len(Num_div)) + ' dimensions.')
    h.quit()
