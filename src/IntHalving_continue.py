# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:01:13 2020

@author: lehotzky
"""
# PRELIMINARIES ===========================================================
# Importing numpy library for computations
import numpy as np
# Importing pickle for saving data structures as objects
import pickle as pkl
# Importing module for getting actual time
from datetime import date
# Importing function that defines class for data saving
from Iteration import Iteration

def fib(lis):
    bif = lis
    n = bif[0]
    d = bif[1]
    s_range = bif[2]
    fun = bif[3]
    # scaling node using stepsize
    node_scaled = np.array(n)*d + s_range[:, 0]
    f_val = fun(node_scaled)
    return n, f_val

def IHalving_cont(pc, fun, It, dN_it, save='end', tol='same'):
    '''
    - Continues an existing interval-halving iteration for equation fun(x)=0
      whose data is stored in the Iteration class instance It.
    - After completing the last iteration step of It that might be unfinished,
      it computes dN_it number of additional iteration steps.
    - All iteration data is saved according to the 'save' argument and when a
      KeyboardInterrupt Exception is raised.
    - Data is saved as an instance of the Iteration class using pickle.

    Args.: fun -        A function whose argument x is an array of coordinates
                        along which sol_range is defined.
                        Its output must be a scalar number.
                        Its argument is an array of numbers.

           It -         an Iteration class instance containing data that is
                        used when continuing iteration

           dN_it        a non-negative integer that marks the number of
                        iteration steps to be added

           save -       It must be one of the following strings:
                        'end'   - data is saved when iteration ended
                        'step'  - data is saved per iteration step
                        'node'  - data is saved each time a new node is
                                  computed

           tol -        Tolerance for zeros: |fun(x)| < tol will be treated as
                        fun(x) having a zero at x.
                        If not 'same' then the updated tol value will be
                        applied only after the last iteration step of It has
                        been completed.pc.done()
                        It must be a positive float, or one of the following
                        strings:
                        'eps'   - indicates machine epsilon
                        'same'  - tolerance is inherited from previous
                                  simulation

    Returns: A tuple whose elements are:

           sol -        A 2D numpy array whose rows correspond to points along
                        the solution surface defined by fun(x) = 0.
                        Its columns are associated with the respective
                        coordinates of x.

           It -         An Iteration object instance corresponding to the
                        executed interval-halving iteration.
    '''


    # UTILITY FUNCTIONS =======================================================

    def save_It(It, file_name):
        '''
        Args.:      - It: Iteration class instance to be saved
                    - file_name: name used when saving
        Returns:    None
        '''
        filehandler = open(f"/media/zupanc-lab/Seagate Backup Plus Drive/2020_Pn_Oscillation/{file_name}", 'wb')
        pkl.dump(It, filehandler)
        filehandler.close()

    # CHECKING INPUT DATA =====================================================
    PC = pc

    if not isinstance(It, Iteration):
        raise ValueError('Argument \'It\' must be an instance of the Iteration class!')

    if not isinstance(dN_it, int):
        raise ValueError('Argument \'dN_it\' must be a non-negative integer!')
    elif dN_it < 0:
        raise ValueError('Argument \'dN_it\' must non-negative!')

    save_allowed = ['end', 'step', 'node']
    if save not in save_allowed:
        msg = '\', \''.join(save_allowed)
        raise ValueError('Argument \'save\' must be one of the following strings: \'' + msg + '\' !')

    if tol != 'eps' and tol != 'same' and not isinstance(tol, float):
        raise ValueError('Argument \'tol\' must be either a positive float or one of the strings \'eps\', \'same\'!')
    else:
        if tol != 'eps' and tol != 'same':
            if tol < 0:
                raise ValueError('Argument \'tol\' must be positive!')

    # Checking whether function name is the same as in previous iteration,
    # changing function name if necessary:
    #It.change_fun(fun)

    # PRE-COMPUTATIONS ========================================================

    # Retriving parameters from Iteration instance:
    s_range = np.array(It.get_sol_range())  # solution range
    N_div = It.get_n_div()                  # number of divisions
    last_step = It.get_last_step()          # last iteration step
    name = It.get_name()                    # name of iteration

    # Computing other parameters:
    N_dim = len(N_div)              # number of dimensions
    N_node = 2**N_dim               # number of nodes per cube
    N_node_ini, N_cube_ini = 1, 1   # number of all initial nodes & cubes
    for i in range(N_dim):
        N_node_ini *= (N_div[i]+1)
        N_cube_ini *= N_div[i]
    D_0 = (s_range[:, 1] - s_range[:, 0])/N_div     # initial stepsizes

    # Cube increment matrix computation:
    D_cube = np.zeros((N_node, N_dim), dtype=int)  # cube increment matrix
    for i in range(N_node):
        D_cube[i, :] = i // (2**(N_dim - 1 - np.arange(0, N_dim))) % 2

    # Cube shift matrix computation (applies only to initial grid generation):
    P = np.zeros(N_dim)     # denominators
    for k in range(N_dim):
        Pj = 1
        if k + 1 < N_dim:
            for j in range(k + 1, N_dim):
                Pj *= N_div[j]
        P[k] = Pj
    D_shift = np.zeros((N_cube_ini, N_dim), dtype=int)     # cube shift matrix
    for i in range(1, N_cube_ini):
        D_shift[i, :] = i // P % N_div

    # Printing massage on resuming iteration:
    print('\nIteration continued from iteration step #' + str(last_step) + '.')

    # Addign log about continuing iteration:
    It.log_continue()

    # Epsilon value used for identifying zeros while completing (possibly)
    # unfinished iteration step:
    eps = It.get_tol()

    # Generating file name used for saving
    date_string = It.get_time().strftime("%Y_%B_%d_%Hh_%Mm")
    name_used = name + (len(name) > 0)*'_'
    file_name = 'It_' + name_used + date_string + '.pkl'

    try:
        # COMPLETING LAST STEP ================================================

        if last_step == 0: # Resume from initialization

            # We recreate all node coordinates and bracketing cubes because it
            # is uncertain where exactly the iteration was interrupted.

            # Reading node dictionary:
            No_read = It.get_nodes_raw(last_step)
            # Cube dictionary is regenerated (as it does not require much time
            # and provided dictionary might be incomplete):
            Cu_pre = {}
            # Generating nodes of 1st cube and adding/rading associated
            # function values:
            No_it = {}
            for i in range(N_node):
                node_i = tuple(D_cube[i, :])
                if node_i not in No_read:
                    No_it[node_i] = np.nan
                    It.add_node(node_i, np.nan)
                else:
                    No_it[node_i] = No_read[node_i]
            # Generating 1st cube:
            Cu1 = ()
            for key in No_it:
                Cu1 += (key,)
            Cu_pre[1] = Cu1
            It.add_cube(1, Cu1)
            # Generating rest of the nodes and cubes and reading associated
            # function values:
            for i in range(1, N_cube_ini):
                D_i = D_shift[i, :]
                Cui = ()
                for node in Cu1:
                    node_ij = tuple(D_i + np.array(node))
                    Cui += (node_ij,)
                    if node_ij not in No_read:
                        No_it[node_ij] = np.nan
                        It.add_node(node_ij, np.nan)
                    else:
                        No_it[node_ij] = No_read[node_ij]
                Cu_pre[i+1] = Cui
                It.add_cube(i+1, Cui)
            # Separating nodes to which no function value has been assigned
            No_new = {}
            count = 0
            for node in No_it:
                if np.isnan(No_it[node]):
                    count += 1
                    No_new[node] = No_it[node]

            print('\nInitialization has been started, ' + str(len(No_new)) + ' number of gridpoints are being added.')

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZE
            # Computing new function values:
            sims_completed = 0
            save_It(It, file_name)
            for node in No_new.keys():
                PC.submit(fib, [node, D_0, s_range, fun])
            while PC.working():
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No_it[n] = val
                It.add_node(n, val)
                if sims_completed % 20 == 0:
                    save_It(It, file_name)
                    print("sim saved")
            print("This iteration's sims are done. Proceeding with saving data.")

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##

            # Collecting bracketing cubes:
            Cu_br = {}
            key = 0
            for cube in Cu_pre:
                nodes = Cu_pre[cube]
                sum_sign = 0
                for node in nodes:
                    gr = No_it[node] > eps
                    le = No_it[node] < -eps
                    sum_sign += int(gr) - int(le)
                if (sum_sign != N_node) and (sum_sign != -N_node):
                    key += 1
                    Cu_br[key] = nodes
            # Saving step data:
            if save == 'step':
                save_It(It, file_name)

        else: # Resume from last interval-halving step

            # We recreate all node coordinates and bracketing cubes because it
            # is undertain where exactly the iteration was interrupted.

            # Reading nodes of last step:
            No_it = It.get_nodes_raw(last_step)
            # Reading nodes of prior step:
            No_prev = It.get_nodes_raw(last_step - 1)

            # Getting cubes of previous iteration:
            Cu_pre = It.get_cubes(last_step - 1)
            # Computing bracketing cubes:
            Cu_br = {}
            key = 0
            for cube in Cu_pre:
                nodes = Cu_pre[cube]
                sum_sign = 0
                for node in nodes:
                    gr = No_prev[node] > eps
                    le = No_prev[node] < -eps
                    sum_sign += int(gr) - int(le)
                if (sum_sign != N_node) and (sum_sign != -N_node):
                    key += 1
                    Cu_br[key] = nodes
            # Rescaling discrete node coordinates:
            No_it_new = {}
            for key in No_prev:
                key_new = tuple(np.array(key)*2)
                f_val = No_prev[key]
                No_it_new[key_new] = f_val
                It.add_node(key_new, f_val)
            # Halving bracketing cubes:
            Cu_pre = {}
            key = 0
            for cube in Cu_br:
                nodes = Cu_br[cube]
                node_min = nodes[0]
                node_min_scaled = np.array(node_min)*2
                Cu1 = (tuple(node_min_scaled),)
                # generating 1st cube and its nodes
                for i in range(1, N_node):
                    node_ij = tuple(node_min_scaled + D_cube[i, :])
                    Cu1 += (node_ij,)
                    if node_ij not in No_it:
                        No_it_new[node_ij] = np.nan
                        It.add_node(node_ij, np.nan)
                    else:
                        No_it_new[node_ij] = No_it[node_ij]
                key += 1
                Cu_pre[key] = Cu1
                It.add_cube(key, Cu1)
                # generating rest of the cubes and their nodes
                for i in range(1, N_node):
                    D_i = D_cube[i, :]
                    key += 1
                    Cui = ()
                    for node in Cu1:
                        node_ij = tuple(D_i + np.array(node))
                        Cui += (node_ij,)
                        if node_ij not in No_it:
                            No_it_new[node_ij] = np.nan
                            It.add_node(node_ij, np.nan)
                        else:
                            No_it_new[node_ij] = No_it[node_ij]
                    Cu_pre[key] = Cui
                    It.add_cube(key, Cui)
            # Halving stepsize:
            D = D_0/(2**last_step)

            # Separating nodes to which no function value has been assigned
            No_new = {}
            count = 0
            for node in No_it_new:
                if np.isnan(No_it_new[node]):
                    count += 1
                    No_new[node] = No_it_new[node]

            print('\nIteration step #' + str(last_step) + ' has been started, ' + str(len(No_new)) + ' number of gridpoints are being added.')
            save_It(It, file_name)
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZED
            # Computing new function values:
            sims_completed = 0
            for node in No_new.keys():
                PC.submit(fib, [node, D, s_range, fun])
            while PC.working():
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No_it_new[n] = val
                It.add_node(n, val)
                if sims_completed % 200000 == 0:
                    save_It(It, file_name)
                    print("sim saved")
            print("This iteration's sims are done. Proceeding with saving data.")
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
            # Collecting bracketing cubes:
            Cu_br = {}
            key = 0
            for cube in Cu_pre:
                nodes = Cu_pre[cube]
                sum_sign = 0
                for node in nodes:
                    sum_sign += np.sign(No_it_new[node])
                if abs(sum_sign) < N_node:
                    key += 1
                    Cu_br[key] = nodes

            Cu_br = {}
            key = 0
            for cube in Cu_pre:
                nodes = Cu_pre[cube]
                sum_sign = 0
                for node in nodes:
                    gr = No_it_new[node] > eps
                    le = No_it_new[node] < -eps
                    sum_sign += int(gr) - int(le)
                if (sum_sign != N_node) and (sum_sign != -N_node):
                    key += 1
                    Cu_br[key] = nodes
            No_it = No_it_new.copy()
            # Saving step data:
            if save == 'step':
                save_It(It, file_name)

        # ADDING NEW INTERVAL-HALVING ITERATION STEPS =========================

        # Updating epsilon if requested by user:
        if tol == 'same':
            eps = It.get_tol()
        elif tol == 'eps':
            eps = np.finfo(float).eps
            It.set_tol(eps)
        else:
            eps = tol
            It.set_tol(eps)

        # Total number of gridpoints we have so far:
        count_total = len(It.get_nodes_raw(last_step))

        for it in range(last_step + 1, last_step + dN_it + 1):

            # Adding iteration step to instance:
            It.add_step()
            # Rescaling discrete node coordinates (considering halved stepsize):
            No_it_new = {}
            for key in No_it:
                key_new = tuple(np.array(key)*2)
                f_val = No_it[key]
                No_it_new[key_new] = f_val
                It.add_node(key_new, f_val)
            # Halving bracketing cubes:
            Cu_pre = {}
            key = 0
            for cube in Cu_br:
                nodes = Cu_br[cube]
                node_min = nodes[0]
                node_min_scaled = np.array(node_min)*2
                # generating 1st cube and its nodes
                Cu1 = (tuple(node_min_scaled),)
                for i in range(1, N_node):
                    node_ij = tuple(node_min_scaled + D_cube[i, :])
                    Cu1 += (node_ij,)
                    if node_ij not in No_it_new:
                        No_it_new[node_ij] = np.nan
                        It.add_node(node_ij, np.nan)
                key += 1
                Cu_pre[key] = Cu1
                It.add_cube(key, Cu1)
                # generating rest of the nodes and cubes
                for i in range(1, N_node):
                    D_i = D_cube[i, :]
                    key += 1
                    Cui = ()
                    for node in Cu1:
                        node_ij = tuple(D_i + np.array(node))
                        Cui += (node_ij,)
                        if node_ij not in No_it_new:
                            No_it_new[node_ij] = np.nan
                            It.add_node(node_ij, np.nan)
                    Cu_pre[key] = Cui
                    It.add_cube(key, Cui)
            # Halved step:
            D = D_0/(2**it)

            # Separating nodes to which no function value has been assigned
            No_new = {}
            count = 0
            for node in No_it_new:
                if np.isnan(No_it_new[node]):
                    count += 1
                    No_new[node] = No_it_new[node]
            count_total += count

            print('\nIteration step #' + str(it) + ' has been started, ' + str(len(No_new)) + ' number of gridpoints are being added.')
            save_It(It, file_name)

## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## PARALLELIZED
            # Computing new function values:
            sims_completed = 0
            for node in No_new.keys():
                PC.submit(fib, [node, D, s_range, fun])
            while PC.working():
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No_it_new[n] = val
                It.add_node(n, val)
                if sims_completed % 400000 == 0:
                    save_It(It, file_name)
                    print("sim saved")
            print("This iteration's sims are done. Proceeding with saving data.")
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##

            # Collecting bracketing cubes:
            Cu_br = {}
            key = 0
            for cube in Cu_pre:
                nodes = Cu_pre[cube]
                sum_sign = 0
                for node in nodes:
                    gr = No_it_new[node] > eps
                    le = No_it_new[node] < -eps
                    sum_sign += int(gr) - int(le)
                if (sum_sign != N_node) and (sum_sign != -N_node):
                    key += 1
                    Cu_br[key] = nodes
            No_it = No_it_new.copy()
            # Saving per step:
            if save == 'step':
                save_It(It, file_name)

        # Collecting nodes that belong to bracketing cubes:
        No_br = {}
        Sol = []
        if It.get_last_step() == 0:
            D_end = D_0
        else:
            D_end = D
        for cube in Cu_br:
            nodes = Cu_br[cube]
            for node in nodes:
                if node not in No_br:
                    No_br[node] = No_it[node]
                    Sol_i = np.zeros(N_dim + 1)
                    Sol_i[:N_dim] = np.array(node)*D_end + s_range[:, 0]
                    Sol_i[N_dim] = No_it[node]
                    Sol.append(Sol_i)
        Sol = np.array(Sol)
        N_SolPt = np.shape(Sol)[0]

        print('\nIteration completed.\nNumber of gridpoints along the final boundary: ' + str(N_SolPt) + '\nEfficiency: %.2f' % float(N_SolPt*100/count_total) + '%')

        # Saving only at end:
        if save == 'end':
            save_It(It, file_name)

        print('\nData is saved as: ' + file_name)

        return (Sol, It)

    except (KeyboardInterrupt, SystemExit):
        raise Exception('Execution terminated by user.\nData is saved as: ' + file_name)

    finally:
        # Saving iteration instance
        save_It(It, file_name)
