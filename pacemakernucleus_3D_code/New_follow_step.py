import numpy as np
import time


def new_follow_step(No, Cu_add, Cu_br, sol_range, num_div, it_step):
    """
    Carries out an iteration step for the solution surface following algorithm by looking for bracketing cubes among
    Cu_add and then generating their neighbors for further check.

    :param No:  dictionary that contains all discrete nodes as keys, and corresponding frequencies as values
    :param Cu_add:  a tuple of tuples that contain (discrete) nodes of newly-generated n-cubes which might be bracketing
    :param Cu_br:   a tuple of tuples that contain (discrete) nodes of all bracketing n-cubes
    :param sol_range:   2D array where first dimension corresponds to coordinates and the second to lower and upper
                        limits
    :param num_div:     1D array that corresponds to the initial number of divisions along each coordinate
    :param it_step:     integer that marks the final iteration step of interval-halving
    :return:    nodes_scaled:   2D numpy array containing scaled coordinates (columns) of all newly-added nodes (rows)
                nodes_discrete: a tuple of tuples that contain (discrete) coordinates of all newly-added nodes
                Cu_add: a tuple of tuples that contain (discrete) nodes of newly-found cubes which might be bracketing
                Cu_br:  a tuple of tuples that contain (discrete) all bracketing cubes including newly-found ones
    """

    # PRE-COMPUTATIONS -------------------------------------------------------------------------------------------------

    eps = np.finfo(float).eps               # tolerance for zero detection
    s_range = np.array(sol_range)           # solution range
    N_div = np.array(num_div, dtype=int)    # number of divisions
    N_dim = len(N_div)                      # number of dimensions
    N_node = 2 ** N_dim                     # number of nodes per cube
    N_node_ini, N_cube_ini = 1, 1           # number of all initial nodes & cubes
    for i in range(N_dim):
        N_node_ini *= (N_div[i] + 1)
        N_cube_ini *= N_div[i]
    D_0 = (s_range[:, 1] - s_range[:, 0]) / N_div  # initial step sizes

    # Cube increment matrix computation:
    D_cube = np.zeros((N_node, N_dim), dtype=int)  # cube increment matrix
    for i in range(N_node):
        D_cube[i, :] = i // (2 ** (N_dim - 1 - np.arange(0, N_dim))) % 2

    # Converting cube increment matrix into a tuple
    D_c = ()
    for i in range(N_node):
        D_c += (tuple(D_cube[i, :]),)

    # Neighboring cubes relative to a cube
    Cu_neig = ()
    for i in range(N_dim):
        Cu_pi = ()
        Cu_ni = ()
        for node in D_c:
            D_p, D_n = np.array(node), np.array(node)
            D_p[i] += 1
            D_n[i] -= 1
            Cu_pi += (tuple(D_p),)
            Cu_ni += (tuple(D_n),)
        Cu_neig += (Cu_pi,)
        Cu_neig += (Cu_ni,)

    # Step sizes of last iteration step:
    D = D_0 / (2 ** it_step)

    # Maximum discrete (integer-number) boundaries of coordinates (min are 0):
    Coord_max = np.array(N_div, dtype=int) * (2 ** it_step)

    # FINDING NEW BRACKETING CUBES AMONG ADDED CUBES -------------------------------------------------------------------

    t0_start = time.perf_counter()

    # Collecting new bracketing cubes:
    Cu_new_br = ()
    for cube in Cu_add:
        sum_sign = 0
        for node in cube:
            gr = No[node] > eps
            le = No[node] < -eps
            sum_sign += int(gr) - int(le)
        if (sum_sign != N_node) and (sum_sign != -N_node):
            Cu_new_br += (cube,)

    # add new bracketing cubes to bracketing cube dictionary
    for cube in Cu_new_br:
        if cube not in Cu_br:
            Cu_br += (cube,)

    t0_end = time.perf_counter()
    print('    New bracketing cubes have been identified. | Time consumption: ', t0_end - t0_start, '[second]')

    # FINDING ALL NEIGHBORING CUBES OF NEW BRACKETING CUBES ------------------------------------------------------------

    t1_start = time.perf_counter()

    # This dictionary stores all nodes that make up the neighboring cubes of new bracketing cubes (these nodes
    # may be outside the search range)
    # keys: discrete nodes
    # values: keys of associated cubes in Cu_new
    No_new = {}

    # This dictionary stores all neighboring cubes of new bracketing cubes
    # keys: sequence numbers used to identify unique cubes
    # values: tuple of discrete nodes
    Cu_new = {}

    # Collecting all neighboring cubes of new bracketing cubes
    key = 0
    for cube in Cu_new_br:  # checking for each bracketing cube
        node_minimal = np.array(cube[0], dtype=int)  # minimal node coordinates of a particular br. cube
        for Cu_n in Cu_neig:
            key += 1
            Cu_i = ()
            for node in Cu_n:
                node_ij = tuple(node_minimal + np.array(node, dtype=int))
                if node_ij in No_new:
                    No_new[node_ij] += (key,)
                else:
                    No_new[node_ij] = (key,)
                Cu_i += (node_ij,)
            Cu_new[key] = Cu_i

    t1_end = time.perf_counter()
    print('    Neighbors of new bracketing cubes and their nodes have been collected. | Time consumption: ',
          t1_end - t1_start, '[second]')

    # FINDING NEW NODES TO BE ADDED ------------------------------------------------------------------------------------

    # Identifying nodes in No_new that
    # - are not yet inside the node dictionary,
    # - and do not fall outside the search range

    # Dictionary that contains only nodes that have not been computed yet and do not fall outside the search
    # range
    # keys: discrete nodes
    # values: keys associated with cubes in Cu_new

    t2_start = time.perf_counter()

    No_add = {}
    for node in No_new:
        node_np = np.array(node, dtype=int)
        diff = Coord_max - node_np
        if np.all(node_np >= 0) and np.all(diff >= 0):
            if node not in No:
                No_add[node] = No_new[node]

    t2_end = time.perf_counter()
    print('    Nodes to be added have been collected. | Time consumption: ', t2_end - t2_start, '[second]')

    # FINDING NEW CUBES TO BE ADDED ------------------------------------------------------------------------------------

    t3_start = time.perf_counter()

    # Tuple that contains keys inside Cu_new, that are associated with cubes incorporating added nodes
    keys_add = ()
    for node in No_add:
        keys = No_add[node]
        for key in keys:
            if key not in keys_add:
                keys_add += (key,)

    # Tuple containing all cubes that incorporate added nodes
    Cu_add = ()
    for key in keys_add:
        if Cu_new[key] not in Cu_add:
            Cu_add += (Cu_new[key],)

    t3_end = time.perf_counter()
    print('    Cubes to be added have been collected. | Time consumption: ', t3_end - t3_start, '[second]')

    # SCALING ADDED NODES ----------------------------------------------------------------------------------------------

    nodes_discrete = tuple(No_add.keys())           # tuple of discrete nodes to be added
    n_node_new = len(No_add)                        # number of new nodes to be added
    nodes_scaled = np.zeros((n_node_new, N_dim))    # scaled coordinates of nodes to be added
    for i in range(n_node_new):
        nodes_scaled[i, :] = np.array(nodes_discrete[i]) * D + s_range[:, 0]

    return nodes_scaled, nodes_discrete, Cu_add, Cu_br

