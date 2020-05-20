# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 19:48:51 2020

@author: David
"""

from datetime import datetime
import numpy as np


class Iteration(object):
    '''
    Creates an instance that corresponds to a specific interval halving
    iteration.
    It saves all data to a dictionary, with keys corresponding to the seq.
    number of the particular iteration step.
    '''

    def __init__(self, name, fun, sol_range, n_div):
        # name of the iteration
        self.name = name
        # keep track of function name
        self.fun_name = fun.__name__
        # keep track of when the iteration step was started
        self.time = datetime.now()
        # log string which keeps track of start, continuation and associated
        # tolerance values
        self.log = ''
        # keep track of the parameters of the interval halving iteration
        self.sol_range = sol_range
        self.n_div = n_div
        # keep track of tolerance for zero
        self.tol = None
        # dictionary that stores all data related to particular interval
        # halving
        self.steps = {}
        # keeping track of iteration steps
        self.step = -1

    def change_fun(self, fun):
        """
        Changes function name associated with the iteration instance.
        This function change occurs only after the approval of user, in order
        to make sure that interval-halving iteration of a particular function
        is not continued using a different function.
        """
        new_name = fun.__name__
        old_name = self.fun_name
        if new_name != old_name:
            print('The function name used in current iteration (' + new_name \
                  + ') does not match with that of used in previous iteration (' \
                  + old_name + ').')
            while True:
                usr_inp = str(input('Do you want to continue? (y/n)'))
                inp = usr_inp.lower()
                if inp != 'y' and inp != 'n':
                    print('Input \'' + usr_inp + '\' is not valid.' + \
                          ' Please enter \'y\' for yes, or \'n\' for no!')
                elif inp == 'y':
                    self.log += '\n' \
                                + str(datetime.now().strftime("%B %d, %Y, %H:%M:%S")) + \
                                ': function name was changed from ' + old_name + ' to ' + \
                                new_name + '.'
                    break
                elif inp == 'n':
                    raise Exception('Interval-halving iteration has been aborted by user.' + \
                                    '\nIteration instance is left unchanged.')

    def add_step(self):
        '''
        Adds a new, empty iteration step to the dictionary of iteration steps.
        Updates node and cube dictionaries accordingly. Adds log about starting
        new iteration step.

        Returns:    - None
        '''

        # keep track of iteration step
        self.step += 1
        # dictionary of nodes corresponding to the particular iteration step
        self.nodes = {}
        # dictionary of cubes prior to interval halving
        self.cubes = {}
        # dictionary of iteration steps containing node cube dictionaries
        self.steps[self.step] = [self.nodes, self.cubes]
        # record the start of iteration step in log
        self.log += '\n' + str(datetime.now().strftime("%B %d, %Y, %H:%M:%S")) + ': iteration step #' + str(
            self.step) + ' was started.'

    def add_node(self, node, val):
        '''
        Adds a new node to the node dictionary or overrrides an existing one
        together with its associated function value corresponding to the last
        iteration step.

        Args.:      - node: a tuple of raw (unscaled, integer) node coordinates
                      val: value associated with the node
        Returns:    - None
        '''
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        self.nodes[node] = val
        self.steps[self.step][0] = self.nodes

    def add_cube(self, cube, nodes):
        '''
        Adds a new cube to the cube dictionary or overrides an existing one
        both associated with the last iteration step.

        Args.:      - cube: an integer, the sequence number of the cube
                      nodes: a tuple of tuples containing the raw (unscaled,
                      integer) node coordinates associated with the cube
        Returns:    - None
        '''
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        self.cubes[cube] = nodes
        self.steps[self.step][1] = self.cubes

    # adding log about continuing an iteration
    def log_continue(self):
        '''
        Adds log about continuing iteration at a particular iteration step.

        Returns:    - None
        '''
        self.log += '\n' + str(datetime.now().strftime(
            "%B %d, %Y, %H:%M:%S")) + ': interval halving was continued from iteration step #' + str(self.step) + '.'

    # setters
    def set_tol(self, tol):
        '''
        Sets tolerance for zero. Adds a log about this change in tolerance.

        Arg.:       - tol: a positive number
        Returns:    - None
        '''
        if tol != 'eps' and not isinstance(tol, (int, float)):
            raise Exception('Tolerance must be a positive number or the string \'eps\'!')
        elif tol != 'eps' and tol < 0:
            raise Exception('Tolerance must be a positive number or the string \'eps\'!')
        self.tol = tol
        self.log += '\n' + str(
            datetime.now().strftime("%B %d, %Y, %H:%M:%S")) + ': tolerance for zero was set to tol=' + str(tol) + '.'

    # getter functions

    def get_name(self):
        return self.name

    def get_tol(self):
        return self.tol

    def get_sol_range(self):
        return self.sol_range

    def get_n_div(self):
        return self.n_div

    def get_time(self):
        return self.time

    def get_last_step(self):
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        return self.step

    def get_cubes(self, step):
        '''
        Arg.:       iteration step for which nodes are requested
        Returns:    dictionary of cubes where keys are cube identifiers and
                    values are tuples of the raw coordinates of associated
                    nodes
        '''
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        elif step > self.step:
            msg = 'Available iteration steps: '
            for i in range(self.step + 1):
                if i > 0:
                    msg += ', '
                msg += str(i)
            raise Exception('Iteration step does not exist!\n' + msg + '\n')
        else:
            if len(self.nodes) == 0:
                raise Exception('No nodes have been added to this iteration step!\n')
        return self.steps[step][1]

    def get_nodes_raw(self, step):
        '''
        Arg.:       iteration step for which nodes are requested
        Returns:    dictionary of nodes where keys are a tuple of raw (integer,
                    without scaling to solution range) coordinates and values
                    are the corresponding function value evauated at that node
        '''
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        elif step > self.step:
            msg = 'Available iteration steps: '
            for i in range(self.step + 1):
                if i > 0:
                    msg += ', '
                msg += str(i)
            raise Exception('Iteration step does not exist!\n' + msg + '\n')
        else:
            if len(self.nodes) == 0:
                raise Exception('No nodes have been added to this iteration step!\n')
        return self.steps[step][0]

    def get_nodes_scaled(self, step):
        '''
        Arg.:       iteration step for which nodes are requested
        Returns:    2D numpy array with rows corresponding to different nodes.
                    Columns are associated with the node's coordinates and the
                    resulting function values in the last column
        '''
        if len(self.steps) == 0:
            raise Exception('No iteration steps have been added yet!')
        elif step > self.step:
            msg = 'Available iteration steps: '
            for i in range(self.step + 1):
                if i > 0:
                    msg += ', '
                msg += str(i)
            raise Exception('Iteration step does not exist!\n' + msg + '\n')
        else:
            if len(self.nodes) == 0:
                raise Exception('No nodes have been added to this iteration step!\n')

        s_range = np.array(self.sol_range)
        nodes = self.steps[step][0]

        D = (s_range[:, 1] - s_range[:, 0]) / self.n_div
        for i in range(step):
            D /= 2
        nodes_scaled = np.zeros((len(nodes), len(self.n_div) + 1))
        i = 0
        for node in nodes:
            node_scaled = np.array(node) * D + s_range[:, 0]
            nodes_scaled[i, :-1] = node_scaled
            nodes_scaled[i, -1] = nodes[node]
            i += 1

        return nodes_scaled

    def get_sol(self, step):
        '''
        Arg.:       iteration step for which nodes defining solution are
                    requested
        Returns:    2D numpy array with rows corresponding to different nodes.
                    Columns are associated with the node's coordinates and the
                    resulting function values in the last column
        '''
        eps = self.tol
        N_it = self.step + (sum(np.isnan(np.array(list(self.nodes.values())))) == 0)
        if N_it == 0:
            raise Exception('No iteration steps have been added yet!')
        elif step > (N_it - 1):
            msg = 'Available iteration steps: '
            for i in range(N_it):
                if i > 0:
                    msg += ', '
                msg += str(i)
            if (N_it - 1) == self.step:
                raise Exception('Iteration step does not exist!\n' + msg + '\n')
            else:
                raise Exception('Iteration step is incomplete!\n' + msg + '\n')
        else:
            if len(self.nodes) == 0:
                raise Exception('No nodes have been added to this iteration step!\n')
        # computing bracketing cubes
        cubes_br = {}
        nodes = self.steps[step][0]
        cubes = self.steps[step][1]
        N_node = len(cubes[1])
        key = 0
        for cube in cubes:
            cube_nodes = cubes[cube]
            sum_sign = 0
            for node in cube_nodes:
                gr = nodes[node] > eps
                le = nodes[node] < -eps
                sum_sign += int(gr) - int(le)
            if (sum_sign != N_node) and (sum_sign != -N_node):
                key += 1
                cubes_br[key] = cube_nodes
        # collecting nodes of bracketing cubes
        No_br = {}
        Sol = []
        N_dim = len(self.n_div)
        s_range = np.array(self.sol_range)
        D = (s_range[:, 1] - s_range[:, 0]) / self.n_div
        for i in range(step):
            D /= 2
        for cube in cubes_br:
            cube_nodes = cubes_br[cube]
            for node in cube_nodes:
                if node not in No_br:
                    No_br[node] = nodes[node]
                    Sol_i = np.zeros(N_dim + 1)
                    Sol_i[:N_dim] = np.array(node) * D + s_range[:, 0]
                    Sol_i[N_dim] = nodes[node]
                    Sol.append(Sol_i)
        Sol = np.array(Sol)
        return Sol

    # printing

    def __str__(self):
        msg = '\nName of interval-halving iteration: ' + self.name + '\n'
        msg += 'Time of start: ' + str(self.time.strftime("%B %d, %Y, %H:%M:%S")) + '\n'
        msg += 'Range of solution: ' + str(self.sol_range) + '\n'
        msg += 'Initial number of divisions: ' + str(self.n_div) + '\n'
        if len(self.steps) > 0:
            msg += 'Sequence number of last iteration step: ' + str(self.step) + '\n'
            if len(self.nodes) > 0:
                n_incomplete = sum(np.isnan(np.array(list(self.nodes.values()))))
                if n_incomplete > 0:
                    it_complete = 'NO\nNumber of incomplete gridpoints: ' + str(n_incomplete)
                else:
                    it_complete = 'YES'
                msg += 'Last iteration step completed? ' + it_complete + '\n'
            else:
                msg += 'No nodes have been associated with last iteration step.\n'
        else:
            msg += 'Iteration has not been started.\n'
        msg += 'Logs:' + self.log + '\n'

        return (msg)
