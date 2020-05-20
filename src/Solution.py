# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 19:48:51 2020

@author: David
"""

class Solution(object):
    '''
    Creates an instance that corresponds to an interval halving Iteration 
    instance.
    After following unfinished solution branches until boundary of solution 
    range is reached, it stores bracketing cubes and their nodes. It also 
    stores interpolated surface parameters fitted on these nodes.
    '''
    
    def __init__(self, It):
        """
        Arg.:   It -    an Iteration instance that contains all data 
                        corresponding to interval-halving iteration
        """
        import numpy as np
        from Iteration import Iteration
        # Checking whether It is an Iteration instance:
        if not isinstance(It, Iteration):
            raise ValueError('Argument \'It\' must be an instance of the Iteration class!')
        # Checking whether last iteration step is completed:
        last_step = It.get_last_step()
        step = last_step
        No = It.get_nodes_raw(step).copy()
        if np.any(np.isnan(np.array(list(No.values())))):
            print('Last iteration step (' + str(last_step) + ') is not completed.')     
            msg = 'Do you want to continue by using data from iteration step #' \
            + str(last_step-1) + '? (y/n): '
            while True:
                usr_inp = str(input(msg))
                inp = usr_inp.lower()
                if inp != 'y' and inp != 'n':
                    print('Input \'' + usr_inp + '\' is not valid.' + \
                          ' Please enter \'y\' for yes, or \'n\' for no!')
                elif inp == 'y':
                    step = last_step - 1
                    No = It.get_nodes_raw(step).copy()
                    break
                elif inp == 'n':
                    raise Exception('Solution patch continuation has been aborted by user.')
        # asociated Iteration class instance
        self.iterations = It
        # name of the iteration
        self.name = It.get_name()
        # time when the iteration step was started
        self.time = It.get_time()
        # logs
        self.log = It.log[:]
        # solution range
        self.sol_range = It.get_sol_range()
        # initial number of divisions along coordinates
        self.n_div = It.get_n_div()
        # tolerance for zero (in last iteration step)
        self.tol = It.get_tol()
        # name of function that was used (name corresponds to last iteration 
        # step)
        self.fun_name = It.fun_name[:]
        # seqence number of last iteration step
        self.step = step
        # dictionary of nodes
        self.nodes = No
        # dictionary of cubes
        self.cubes = It.get_cubes(step).copy()
        
    def change_fun(self, fun):
        """
        Changes function name from what is associated with the Iteration 
        instance.
        This function change occurs only after the approval of user, in order 
        to make sure that solution curve continuation is not executed for a 
        function different from the one recorded in Iteration instance.
        """
        from datetime import datetime
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
                    raise Exception('Interval-halving iteration has been aborted by user.'+\
                                    '\nIteration instance is left unchanged.')
                
    def add_node(self, node, val):
        '''
        Adds a new node to the node dictionary or overrrides an existing one 
        together with its associated function value
        
        Args.:      - node: a tuple of raw (unscaled, integer) node coordinates
                      val: value associated with the node
        Returns:    - None
        '''
        self.nodes[node] = val
    
    def add_cube(self, cube, nodes):
        '''
        Adds a new cube to the cube dictionary or overrides an existing one 
        both associated with the last iteration step.
        
        Args.:      - cube: an integer, the sequence number of the cube
                      nodes: a tuple of tuples containing the raw (unscaled, 
                      integer) node coordinates associated with the cube
        Returns:    - None
        '''
        self.cubes[cube] = nodes
    
    # adding log about following unfinished solution patches
    def log_follow(self):
        '''
        Adds log about continuation of unfinished solution patches produced by 
        final iteration step.
        
        Returns:    - None
        '''
        from datetime import datetime
        self.log += '\n' + str(datetime.now().strftime("%B %d, %Y, %H:%M:%S")) + \
        ': continuation of solution patches was started from iteration step #' + str(self.step) + '.'
        
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
        return self.step
    
    def get_iterations(self):
        '''
        Returns:    an Iteration instance that contains data of all 
                    interval-halving iteration steps
        '''
        return self.iterations
        
    def get_cubes(self):
        '''
        Returns:    dictionary of cubes where keys are cube identifiers and 
                    values are tuples of the raw coordinates of associated 
                    nodes
        '''
        return self.cubes
    
    def get_nodes_raw(self):
        '''
        Returns:    dictionary of nodes where keys are tuples of raw (integer,
                    without scaling to solution range) coordinates and values 
                    are corresponding function values evauated at the nodes
        '''
        return self.nodes
    
    def get_nodes_scaled(self):
        '''
        Returns:    2D numpy array with rows corresponding to different nodes. 
                    Columns are associated with the node's coordinates and the 
                    associated function values in the last column.
        '''          
        import numpy as np
        
        s_range = np.array(self.sol_range)
        nodes = self.nodes
        
        D = (s_range[:, 1] - s_range[:, 0])/self.n_div
        for i in range(self.step):
            D /= 2
        nodes_scaled = np.zeros((len(nodes), len(self.n_div) + 1))
        i = 0
        for node in nodes:
            node_scaled = np.array(node)*D + s_range[:, 0]
            nodes_scaled[i, :-1] = node_scaled
            nodes_scaled[i, -1] = nodes[node]
            i += 1
            
        return nodes_scaled
    
    def get_sol(self):
        '''
        Returns:    2D numpy array with rows corresponding to different nodes. 
                    Columns are associated with the node's coordinates and
                    corresponding function values in the last column
        '''
        import numpy as np
        # tolerance for zero value
        eps = self.tol
        # computing bracketing cubes
        nodes = self.nodes
        cubes = self.cubes
        N_node = len(cubes[1])
        key = 0
        cubes_br = {}
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
        D = (s_range[:, 1] - s_range[:, 0])/self.n_div
        for i in range(self.step):
            D /= 2
        for cube in cubes_br:
            cube_nodes = cubes_br[cube]
            for node in cube_nodes:
                if node not in No_br:
                    No_br[node] = nodes[node]
                    Sol_i = np.zeros(N_dim + 1)
                    Sol_i[:N_dim] = np.array(node)*D + s_range[:, 0]
                    Sol_i[N_dim] = nodes[node]
                    Sol.append(Sol_i)
        Sol = np.array(Sol)
        return Sol
    
    def get_precision(self):
        '''
        Returns:    a list of floats that contains the final size of bracketing 
                    cubes along the respective coordinates
        '''
        import numpy as np
        s_range = np.array(self.sol_range)
        n_div = np.array(self.n_div)
        return list((s_range[:, 1] - s_range[:, 0])/(n_div*(2**self.step)))
    
    def get_efficiency(self):
        '''
        Returns:    ratio of nodes along solution surface vs. all nodes, 
                    characterizes the efficiency of the interval-halving method
        '''
        import numpy as np
        N_all = len(self.get_nodes_raw())
        N_sol = np.shape(self.get_sol())[0]
        return N_sol/N_all        
    # printing    
    def __str__(self):
        msg = '\nName of interval-halving iteration: ' + self.name + '\n'
        msg += 'Name of analyzed function: ' + self.fun_name + '\n'
        msg += 'Time of start: ' + str(self.time.strftime("%B %d, %Y, %H:%M:%S")) + '\n'
        msg += 'Range of solution: ' + str(self.sol_range) + '\n'
        msg += 'Initial number of divisions: ' + str(self.n_div) + '\n'
        msg += 'Number of interval-halving iterations: ' + str(self.step) + '\n'
        msg += 'Final step sizes: ' + str(self.get_precision()) + '\n'
        msg += 'Final efficiency: %.2f' % float(self.get_efficiency()*100) + '%\n'
        msg += 'Logs:' + self.log + '\n'
        return(msg)