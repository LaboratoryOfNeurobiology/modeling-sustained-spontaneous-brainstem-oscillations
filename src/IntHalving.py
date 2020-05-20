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
    node_scaled = np.array(n) * d + s_range[:, 0]
    f_val = fun(node_scaled)
    return n, f_val


def IHalving(pc, fun, sol_range, Num_div, N_it, name='', save='end', tol='eps'):
    '''
    - Carries out N_it number of interval halving iterations in order to find
      the solution to equation fun(x) = 0 along the range of x defined by
      sol_range.
    - Iteration starts from an initial equidistant grid of sol_range, defined
      by Num_div.
    - All iteration data is saved according to the 'save' argument and when a
      KeyboardInterrupt Exception is raised.
    - Data is saved as an instance of the Iteration class using pickle.

    Args.: fun -        A function whose argument x is an array of coordinates
                        along which sol_range is defined.
                        Its output must be a scalar number.
                        Its argument is an array of numbers.

           sol_range -  Range of the argument x of fun(x) along which the zeros
                        of fun(x) = 0 are sought.
                        It must be a list whose elements are lists of two
                        numbers that mark the lower and upper bounds for the
                        corresponding coordinate of x.
                        Its length must be equal to that of x.

           Num_div -    Number of divisions used along coordinates of x.
                        It defines the granularity of the initial grid on
                        sol_range.
                        It must be a list of positive integers.
                        Its length must be equal to that of x.

           N_it -       Number of interval halving iterations to be completed.
                        The final accuracy (final cube dimension) along
                        coordinate i of x is:

                        (sol_range[i][1] - sol_range[i][0])/(Num_div[i]*2**N_it)

           name -       Arbitrary name, that is used when saving data.
                        It must be a string.

           save -       It must be one of the following strings:
                        'end'   - data is saved when iteration ended
                        'step'  - data is saved per iteration step
                        'node'  - data is saved each time a new node is
                                  computed

           tol -        Tolerance for zeros: |fun(x)| < tol will be treated as
                        fun(x) having a zero at x.
                        It must be a positive float, or the string 'eps',
                        indicating machine epsilon.

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
        return

    # CHECKING INPUT DATA =====================================================
    s_range = np.array(sol_range);
    N_div = np.array(Num_div, dtype=int);
    # Number of dimensions
    N_dim = len(N_div)
    PC = pc
    # Checking whether function inputs are consistent
    if len(s_range) != len(N_div):
        raise ValueError('Length of arguments \'sol_range\' and \'N_div\' must be equal!')

    if not isinstance(name, str):
        raise ValueError('Argument \'name\' must be string!')

    save_allowed = ['end', 'step', 'node']
    if save not in save_allowed:
        msg = '\', \''.join(save_allowed)
        raise ValueError('Argument \'save\' must be one of the following strings: \'' + msg + '\' !')

    if tol != 'eps' and not isinstance(tol, float):
        raise ValueError('Argument \'tol\' must be either a positive float or the string \'eps\'!')
    else:
        if tol != 'eps':
            if tol < 0:
                raise ValueError('Argument \'tol\' must be positive!')

    # PRE-COMPUTATIONS ========================================================
    # epsilon used for identifying zeros
    if tol == 'eps':
        eps = np.finfo(float).eps
    else:
        eps = tol

    # Number of nodes per cube
    N_node = 2 ** N_dim
    print(f"{N_node} nodes per cube\n")

    # Computing number of all initial nodes & cubes
    N_node_ini, N_cube_ini = 1, 1
    for i in range(N_dim):
        N_node_ini *= (N_div[i] + 1)
        N_cube_ini *= N_div[i]

    print(f"{N_node_ini} initial nodes")
    print(f"{N_cube_ini} initial cubes")
    # Stepsize vector:
    # - numpy array of (initial) stepsizes along coordinates
    D_0 = (s_range[:, 1] - s_range[:, 0]) / N_div

    # Cube increment matrix:
    # - 2D numpy array of shape (N_node, N_dim)
    # - its rows contain the non-scaled (integer) coordinates of each node within
    #   a cube
    # - coordinates are given relative to the cube's node with minimal coordinates
    #   (minimal node)
    # - coordinates of this vector are given in stepsize units
    D_cube = np.zeros((N_node, N_dim), dtype=int)
    for i in range(N_node):
        D_cube[i, :] = i // (2 ** (N_dim - 1 - np.arange(0, N_dim))) % 2

    # Denominators in cube shift matrix calculations:
    P = np.zeros(N_dim)
    for k in range(N_dim):
        Pj = 1
        if k + 1 < N_dim:
            for j in range(k + 1, N_dim):
                Pj *= N_div[j]
        P[k] = Pj

    # Cube shift matrix for initial grid:
    # - 2D numpy array of shape (N_cube_ini, N_dim)
    # - we generate new cubes by shifting the cube with minimal coordinates
    #   (minimal cube)
    # - each row of this matrix corresponds to an N_dim vector that defines the
    #   vector of shifting from the minimal cube
    # - coordinates of this vector are given in stepsize units
    D_shift = np.zeros((N_cube_ini, N_dim), dtype=int)
    for i in range(1, N_cube_ini):
        D_shift[i, :] = i // P % N_div

    # Creation of Iteration class instance and file name
    It = Iteration(name, fun, sol_range, Num_div)

    # recording tolerance for zero
    It.set_tol(eps)

    # file name for saving
    date_string = It.get_time().strftime("%Y_%B_%d_%Hh_%Mm")
    name_used = name + (len(name) > 0) * '_'
    file_name = 'It_' + name_used + date_string + '.pkl'

    print('\nIteration started.\n')

    try:

        # INITIALIZATION ======================================================
        # Here we generate the initial grid and determine its bracketing cubes.

        # Node dictionary:
        # - dictionary of non-scaled (integer) nodes at (initial) iteration
        #   step
        # - keys:   tuples of integers corresponding to the coordiates of each
        #           node (coordinates are in stepsize units)
        # - vals:   floats or numpy NaN-s corresponding to function values at
        #           the particular node
        No_it = {}

        # Dictionary of cubes corresponding to a particular (here 0th)
        # iteration step:
        # - contains cubes and node keys of their nodes from No_it
        # - keys:   integer numbers (cube IDs)
        # - vals:   a tuple of tuples that are keys from No_it
        Cu_pre = {}

        # Adding iteration step to instance:
        It.add_step()

        # Nodes of 1st cube are the rows of the cube increment matrix:
        for i in range(N_node):
            node_i = tuple(D_cube[i, :])
            No_it[node_i] = np.nan
            # adding nodes to step in instance
            It.add_node(node_i, np.nan)

        # Generating 1st cube:
        Cu1 = ()
        for key in No_it:
            Cu1 += (key,)
        Cu_pre[1] = Cu1
        # adding cube to iteration step in instance
        It.add_cube(1, Cu1)

        # Generating rest of the cubes:
        for i in range(1, N_cube_ini):
            D_i = D_shift[i, :]
            Cui = ()
            for node in Cu1:
                node_ij = tuple(D_i + np.array(node))
                Cui += (node_ij,)
                if node_ij not in No_it:
                    No_it[node_ij] = np.nan
                    # adding nodes to step in instance
                    It.add_node(node_ij, np.nan)
            Cu_pre[i + 1] = Cui
            # adding cube to iteration step in instance
            It.add_cube(i + 1, Cui)

        print('Initialization has been started, ' + str(len(No_it)) + ' number of gridpoints are being added.')

        ## ------------------------------------------------------------------------- ##
        ## ------------------------------------------------------------------------- ##
        ## PARALLELIZED
        # Evaluating the function at all nodes with scaled coordinates
        save_It(It, file_name)
        sims_completed = 0
        for node in No_it.keys():
            PC.submit(fib, [node, D_0, s_range, fun])
            #node_scaled = np.array(node) * D_0 + s_range[:, 0]
            #PC.submit(fun, node_scaled, node)
        while PC.working():
            print(f"sim result {sims_completed} appended")
            sims_completed += 1
            n, val = PC.pyret()
            No_it[n] = val
            It.add_node(n, val)
            if sims_completed % 2000 == 0:
                save_It(It, file_name)
                print("sim saved")
        print(f"This iteration's sims are done. = {sims_completed}")

        ## ------------------------------------------------------------------------- ##
        ## ------------------------------------------------------------------------- ##
        print("All saved.")
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

        # Saving per step:
        if save == 'step':
            save_It(It, file_name)

        # ITERATION (INTERVAL-HALVING) ============================================
        # Here we refine (halve along each dimension) bracketing cubes and
        # determine new (smaller) bracketing cubes.

        count_total = N_node_ini
        for it in range(N_it):
            # Adding iteration step to instance:
            It.add_step()

            # Rescaling discrete node cpyplotoordinates (considering halved stepsize):
            No_it_new = {}
            for key in No_it:
                key_new = tuple(np.array(key) * 2)
                f_val = No_it[key]
                No_it_new[key_new] = f_val
                # adding nodes to step in instance
                It.add_node(key_new, f_val)
            # Halving bracketing cubes:
            Cu_pre = {}
            key = 0
            for cube in Cu_br:
                nodes = Cu_br[cube]
                # Minimal node is used (1st node):
                node_min = nodes[0]
                # Scaling coordinates:
                node_min_scaled = np.array(node_min) * 2
                Cu1 = (tuple(node_min_scaled),)
                # Generating 1st (minimal) cube and its nodes:
                for i in range(1, N_node):
                    # Adding elements of the cube increment matrix:
                    node_ij = tuple(node_min_scaled + D_cube[i, :])
                    Cu1 += (node_ij,)
                    if node_ij not in No_it_new:
                        No_it_new[node_ij] = np.nan
                        # adding nodes to step in instance
                        It.add_node(node_ij, np.nan)
                key += 1
                Cu_pre[key] = Cu1
                It.add_cube(key, Cu1)
                # Shifting 1st cube (and thereby generating new ones):
                for i in range(1, N_node):
                    D_i = D_cube[i, :]
                    key += 1
                    Cui = ()
                    for node in Cu1:
                        node_ij = tuple(D_i + np.array(node))
                        Cui += (node_ij,)
                        if node_ij not in No_it_new:
                            No_it_new[node_ij] = np.nan
                            # adding nodes to step in instance
                            It.add_node(node_ij, np.nan)
                    Cu_pre[key] = Cui
                    It.add_cube(key, Cui)
            # Halved step:
            D = D_0 / (2 ** (it + 1))
            # Separating nodes to which no function value has been assigned
            No_new = {}
            count = 0
            for node in No_it_new:
                if np.isnan(No_it_new[node]):
                    count += 1
                    No_new[node] = No_it_new[node]
            count_total += count

            print('\nIteration step #' + str(it + 1) + ' has been started, ' + str(
                len(No_new)) + ' number of gridpoints are being added.')

            ## ------------------------------------------------------------------------- ##
            ## ------------------------------------------------------------------------- ##
            ## PARALLELIZED
            # Evaluating the function at all new nodes using scaled
            # (non-integer) coordinates:
            sims_completed = 0
            for node in No_new.keys():
                PC.submit(fib, [node, D, s_range, fun])
                #node_scaled = np.array(node) * D + s_range[:, 0]
                #PC.submit(fun, node_scaled, node)
            while PC.working():
                print(f"sim result {sims_completed} appended")
                sims_completed += 1
                n, val = PC.pyret()
                No_it_new[n] = val
                It.add_node(n, val)
                if sims_completed % 2000 == 0:
                    save_It(It, file_name)
                    print("sim saved")
            print("This iteration's sims are done.")
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

        print("step 0")
        # Collecting nodes that belong to bracketing cubes:
        No_br = {}
        Sol = []
        if It.get_last_step() == 0:
            D_end = D_0
        else:
            D_end = D

        print("Step 1")
        for cube in Cu_br:
            nodes = Cu_br[cube]
            for node in nodes:
                print("step 3")
                if node not in No_br:
                    No_br[node] = No_it[node]
                    Sol_i = np.zeros(N_dim + 1)
                    Sol_i[:N_dim] = np.array(node) * D_end + s_range[:, 0]
                    Sol_i[N_dim] = No_it[node]
                    Sol.append(Sol_i)
                    print("step 2")
        Sol = np.array(Sol)
        N_SolPt = np.shape(Sol)[0]

        print('\nIteration completed.\nNumber of gridpoints along the final boundary: ' + str(
            N_SolPt) + '\nEfficiency: %.2f' % float(N_SolPt * 100 / count_total) + '%')

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
