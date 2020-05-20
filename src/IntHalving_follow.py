# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:01:13 2020

@author: David
"""
# PRELIMINARIES ===========================================================
# Importing numpy library for computations
import numpy as np
# Importing pickle for saving data structures as objects
import pickle as pkl
# Importing module for getting actual time
from datetime import date
# Importing classes used for data saving
from Iteration import Iteration
from Solution import Solution

def fib(lis):
    bif = lis
    n = bif[0]
    d = bif[1]
    s_range = bif[2]
    fun = bif[3]
    # scaling node using stepsize
    node_scaled = np.array(n)*d + s_range[:, 0]
    f_val = fun(node_scaled)
    return (n, f_val)

def IHalving_follow(pc, fun, It_or_Sol, save='end'):
    """
    Using the cube size of last finished interval halving step, it calculates
    the solution of unfinished solution surface patches of the interval halving
    problem associated with function fun and with Iteration class instance It.

    Args.: fun -        A function whose argument x is an array of coordinates
                        along which we look for solutions.
                        Its output must be a scalar number.
                        Its argument is an array of numbers.

           It_or_Sol -  An Iteration or Solution class instance containing data
                        that is used when following unfinished patches of
                        solution surface.

           save -       It must be one of the following strings:
                        'end'   - data is saved when iteration ended
                        'node'  - data is saved each time a new node is
                                  computed
    """

    # UTILITY FUNCTIONS =======================================================

    def save_Sol(Sol, file_name):
        '''
        Args.:      - Sol:          Solution class instance to be saved
                    - file_name:    name used when saving
        Returns:    None
        '''
        filehandler  = open(file_name, 'wb')
        pkl.dump(Sol, filehandler)
        filehandler.close()

    # CHECKING INPUT DATA =====================================================
    PC = pc

    msg = 'Argument \'It_or_Sol\' must be an instance of either the Iteration'\
    + ' or the Solution class!'
    if isinstance(It_or_Sol, Iteration):
        # Creating new Solution instance:
        Sol = Solution(It_or_Sol)
    elif isinstance(It_or_Sol, Solution):
        # Loading Solution instance:
        Sol = It_or_Sol
    else:
        raise ValueError(msg)

    save_allowed = ['end', 'node']
    if save not in save_allowed:
        msg = '\', \''.join(save_allowed)
        raise ValueError('Argument \'save\' must be one of the following strings: \'' + msg + '\' !')

    # PRE-COMPUTATIONS ========================================================

    # Checking whether function name is the same as in previous iteration,
    # changing function name if necessary:
    print('sol cahange fun')
    Sol.change_fun(fun)

    # Adding log on the continuation of solution patches
    print('log follow')
    Sol.log_follow()

    # Epsilon value is inherited from the final iteration step
    print('get eps')
    eps = Sol.get_tol()

    # Retriving parameters from Iteration instance:
    s_range = np.array(Sol.get_sol_range())  # solution range
    N_div = Sol.get_n_div()                  # initial number of divisions
    name = Sol.get_name()                    # name of iteration

    # Generating file name used for saving
    date_string = Sol.get_time().strftime("%Y_%B_%d_%Hh_%Mm")
    name_used = name + (len(name) > 0)*'_'
    file_name = 'Sol_' + name_used + date_string + '.pkl'

    print('computing other parameters')
    # Computing other parameters:
    N_dim = len(N_div)              # number of dimensions
    N_node = 2**N_dim               # number of nodes per cube
    N_node_ini, N_cube_ini = 1, 1   # number of all initial nodes & cubes
    for i in range(N_dim):
        N_node_ini *= (N_div[i]+1)
        N_cube_ini *= N_div[i]
    D_0 = (s_range[:, 1] - s_range[:, 0])/N_div     # initial stepsizes

    print('Cube increment matrix computation:')
    D_cube = np.zeros((N_node, N_dim), dtype= int)  # cube increment matrix
    for i in range(N_node):
        D_cube[i, :] = i // (2**(N_dim - 1 - np.arange(0, N_dim))) % 2

    print('Cube neighbor matrix computation (used for generating all neighbors of a cube')
    D_neig = np.zeros((2*N_dim, N_dim), dtype= int) # cube neighbor matrix
    for i in range(N_dim):
        D_neig[2*i, i] = -1
        D_neig[2*i + 1, i] = 1

    # Getting number of completed interval-halving steps
    load_step = Sol.get_last_step()

    # Getting nodes of last iteration step
    No = Sol.get_nodes_raw()

    # Stepsizes of last iteration step:
    D = D_0/(2**load_step)

    print('Maximum discrete (integer-number) boundaries of coordinates (min are 0):')
    Coord_max = np.array(N_div)*(2**load_step)

    try:
        print('start try')
        if isinstance(It_or_Sol, Solution):
            No_new = {}
            for node in No:
                print(f"{node} feelsbadman")
                if np.isnan(No[node]):
                    No_new[node] = No[node]

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZED
            print('starting sim submission')
            sims_completed = 0
            for node in No_new.keys():
                PC.submit(fib, [node, D, s_range, fun])
            while(PC.working()):
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No[n] = val
                Sol.add_node(n, val)
                if sims_completed % 2000 == 0:
                    save_Sol(Sol, file_name)
                    print("sim saved")
            print("This iteration's sims are done. Proceeding with saving data.")
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##

        # Computing bracketing cubes of last completed iteration step:
        print('doing cube stuff')
        Cu_pre = Sol.get_cubes()
        Cu_br = {}
        key = 0
        for cube in Cu_pre:
            nodes = Cu_pre[cube]
            sum_sign = 0
            for node in nodes:
                gr = No[node] > eps
                le = No[node] < -eps
                sum_sign += int(gr) - int(le)
            if (sum_sign != N_node) and (sum_sign != -N_node):
                key += 1
                Cu_br[key] = nodes

        print('Generating neighbors of all bracketing cubes:')

        # dictionary for storing nodes that are generated by shifting bracketing
        # cubes in all coordinate directions
        No_new = {}
        # dictionary of newly generated cubes
        Cu_new = {}
        # list containing all cubes
        Cu_list = list(Cu_pre.values())
        key = 0
        print('triple forloop')
        print(len(Cu_br))
        print(2*N_dim)
        ct = 0
        for cube in Cu_br: # checking each bracketing cube
            nodes = Cu_br[cube] # getting node coordinates of a particular br. cube
            for i in range(2*N_dim): # shifting br. cube in all coordinate diections
                D_i = D_neig[i, :] # shift in a particular direction
                Cu_i = () # tuple collecting nodes of new cube
                count = 0 # counter of nodes that lie outside the solution range
                for node in nodes: # for all nodes of br. cube
                    node_ij = tuple(D_i + np.array(node)) # compute shifted value
                    Cu_i += (node_ij,) # add shifted node to new cube
                    # boolean marking whether node is outside of solution range
                    outside = np.any(np.array(node_ij) < 0) or np.any(np.array(node_ij) > Coord_max)
                    count += (outside)
                    # if new node is non-existent and does not fall outside of solution range
                    if (node_ij not in No) and not outside:
                        No_new[node_ij] = np.nan # add new node
                        Sol.add_node(node_ij, np.nan) # add new node to Sol
                if Cu_i not in Cu_list and count == 0: # if any of the nodes of new cube is new then store it as a new cube
                    key += 1
                    Cu_new[key] = Cu_i
                ct +=1
                print(ct)

        print('Adding new cubes to cube dictionary:')
        key = len(Cu_pre)
        for cube in Cu_new:
            key += 1
            Cu_pre[key] = Cu_new[cube]

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZED
        # Computing function values at new nodes
        sims_completed = 0
        for node in No_new.keys():
            PC.submit(fib, [node, D, s_range, fun])
        while(PC.working()):
            sims_completed += 1
            n, val = PC.pyret()
            No[n] = val
            Sol.add_node(n, val)
            if sims_completed % 2000 == 0:
                save_Sol(Sol, file_name)
                print("sim saved")
        print("This iteration's sims are done. Proceeding with saving data.")
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##

        # Finding new bracketing cubes:
        Cu_new_br = {}
        key = 0
        for cube in Cu_new:
            nodes = Cu_new[cube]
            sum_sign = 0
            for node in nodes:
                gr = No[node] > eps
                le = No[node] < -eps
                sum_sign += int(gr) - int(le)
            if (sum_sign != N_node) and (sum_sign != -N_node):
                key += 1
                Cu_new_br[key] = nodes

        while len(Cu_new_br) > 0:

            # add new bracketing cubes to bracketing cube dictionary
            key = len(Cu_br)
            for cube in Cu_new_br:
                key += 1
                Cu_br[key] = Cu_new_br[cube]

            No_new = {}
            Cu_new = {}
            Cu_list = list(Cu_pre.values())
            key = 0
            for cube in Cu_new_br:
                nodes = Cu_new_br[cube]
                for i in range(2*N_dim):
                    D_i = D_neig[i, :]
                    Cu_i = ()
                    count = 0
                    for node in nodes:
                        node_ij = tuple(D_i + np.array(node))
                        Cu_i += (node_ij,)
                        outside = np.any(np.array(node_ij) < 0) or np.any(np.array(node_ij) > Coord_max)
                        count += (outside)
                        if (node_ij not in No) and not outside:
                            No_new[node_ij] = np.nan
                            Sol.add_node(node_ij, np.nan)
                    if Cu_i not in Cu_list and count == 0:
                        key += 1
                        Cu_new[key] = Cu_i

            # Adding new cubes to cube dictionary:
            key = len(Cu_pre)
            for cube in Cu_new:
                key += 1
                Cu_pre[key] = Cu_new[cube]

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZED
            # Computing function values at new nodes
            sims_completed = 0
            for node in No_new.keys():
                PC.submit(fib, [node, D, s_range, fun])
            while(PC.working()):
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No[n] = val
                Sol.add_node(n, val)
                if sims_completed % 2000 == 0:
                    save_Sol(Sol, file_name)
                    print("sim saved")
            print("This iteration's sims are done. Proceeding with saving data.")
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##

            # Finding new bracketing cubes:
            Cu_new_br = {}
            key = 0
            for cube in Cu_new:
                nodes = Cu_new[cube]
                sum_sign = 0
                for node in nodes:
                    gr = No[node] > eps
                    le = No[node] < -eps
                    sum_sign += int(gr) - int(le)
                if (sum_sign != N_node) and (sum_sign != -N_node):
                    key += 1
                    Cu_new_br[key] = nodes

        # Iteration ended, save solution
        if save == 'end':
            save_Sol(Sol, file_name)

        # Printing summary:
        N_all = len(Sol.get_nodes_raw())
        N_new = N_all - len(Sol.get_iterations().get_nodes_raw(load_step))
        msg = '\nOverall number of gridpoints: ' + str(N_all)  + '\n'
        msg += 'Number of gridpoints added by continuation of solution patches: ' +\
        str(N_new) + ' (%.2f' % float(N_new*100/N_all) + '%)'
        print(msg)

        return Sol.get_sol(), Sol

    except (KeyboardInterrupt, SystemExit):
        raise Exception('Execution terminated by user.\nData is saved as: ' + file_name)

    finally:
        # Saving iteration instance
        save_Sol(Sol, file_name)
