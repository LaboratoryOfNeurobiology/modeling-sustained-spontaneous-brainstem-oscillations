import numpy as np


def initialize(sol_range, num_div):
    """
    Generates initial grid (nodes and cubes) for interval halving.

    :param sol_range: 2D array where first dimension corresponds to coordinates and the second to lower and upper limits
    :param num_div: 1D array that corresponds to the initial number of divisions along each coordinate
    :return: (nodes_compute, No, Cu), where
            -   nodes_compute: a numpy array where rows correspond to grid points at which simulations have to be run,
                and columns are associated with respective coordinate values
            -   nodes_discrete: a tuple of tuples that contain coordinates for all nodes particular to the
                initialization step
            -   Cu: a tuple of tuples that contain (discrete) nodes of n-cubes, Cu is particular to the initialization
                step
    """
    s_range = np.array(sol_range)
    n_div = np.array(num_div)
    n_dim = len(num_div)            # number of coordinates
    n_node_per_cube = 2**n_dim      # number of nodes per cube

    # number of initial nodes & cubes
    n_node_ini, n_cube_ini = 1, 1
    for i in range(n_dim):
        n_node_ini *= (n_div[i] + 1)
        n_cube_ini *= n_div[i]

    # step sizes along coordinates
    D_0 = (s_range[:, 1] - s_range[:, 0])/n_div

    # cube increment matrix
    D_cube = np.zeros((n_node_per_cube, n_dim), dtype=int)
    for i in range(n_node_per_cube):
        D_cube[i, :] = i // (2 ** (n_dim - 1 - np.arange(0, n_dim))) % 2

    # denominators in cube shift matrix
    P = np.zeros(n_dim)
    for k in range(n_dim):
        Pj = 1
        if k + 1 < n_dim:
            for j in range(k + 1, n_dim):
                Pj *= n_div[j]
        P[k] = Pj

    # cube shift matrix
    D_shift = np.zeros((n_cube_ini, n_dim), dtype=int)
    for i in range(1, n_cube_ini):
        D_shift[i, :] = i // P % n_div

    nodes_discrete = ()     # nodes
    Cu = ()     # cubes

    # nodes of the 1-st cube
    for i in range(n_node_per_cube):
        node_i = tuple(D_cube[i, :])
        nodes_discrete += (node_i,)

    # generating 1-st cube
    Cu_1 = ()
    for key in nodes_discrete:
        Cu_1 += (key,)
    Cu += (Cu_1,)

    # generating the rest of the cubes
    for i in range(1, n_cube_ini):
        D_i = D_shift[i, :]
        Cu_i = ()
        for node in Cu_1:
            node_ij = tuple(D_i + np.array(node))
            Cu_i += (node_ij,)
            if node_ij not in nodes_discrete:
                nodes_discrete += (node_ij,)
        Cu += (Cu_i,)

    # nodes
    nodes_compute = np.zeros((n_node_ini, n_dim))
    for (i, node) in zip(range(len(nodes_discrete)), nodes_discrete):
        nodes_compute[i, :] = np.array(node)*D_0 + s_range[:, 0]

    return nodes_compute, nodes_discrete, Cu
