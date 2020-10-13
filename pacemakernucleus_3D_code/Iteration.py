import math
import pickle as pkl
import numpy as np
import FileTree as ft
from datetime import datetime


class Iteration(object):
    """
    class Iteration
        Contains all relevant data for saving and reconstruction of an
        n-dimensional parameter space grid.
        :param func_name | name of the set of simulations that this iteration
        belongs to
        :param sol_range | (List of Lists of Integers) = Ranges of each
        dimension in the parameter space
        :param num_div | (List of Integers) = Values corresponding to dimension
        granularity
        :param it_step | (Integer) = Current iteration step
        :param nodes_scaled | (2D np Array) = Rows correspond to new grid points
        at which simulations have to be
            run, and columns are associated with respective coordinate values
        :param nodes_discrete | (Tuple of Tuples) = Contains discrete
        coordinates for all nodes at which frequency
            values have
                to be computed if self.status == "pre" or "inter"
                just been computed if self.status == "post"
            and is associated with iteration step it_step
        :param cu | (Tuple of Tuples) = Contains (discrete) nodes of n-cubes,
        Cu is particular to it_step
        :param no | (Dictionary) = Contains all discrete nodes that have been
         computed as keys, and corresponding
            frequencies as values, both up to iteration step it_step
    """

    def __init__(self, func_name, sol_range, num_div, it_step, nodes_scaled,
                 nodes_discrete, cu, no, status=None, next_node_idx=0):
        self.func_name = func_name  # Name of the simulation
        self.sol_range = sol_range
        self.num_div = num_div
        self.it_step = it_step
        self.nodes_scaled = nodes_scaled  # nodes to to be simulated
        self.nodes_discrete = nodes_discrete  # nodes already simulated
        self.Cu = cu
        self.No = no  # Node dictionary contains nodes whose freqs have been
        # acquired
        self.file_tree = ft.FileTree()
        # Identifiers
        self.status = status  # the current date&time + one of "_pre",
        # "_inter", "_post", updated at each save after initialization
        self.file_name = f"{str(self.status)}_sim_iteration_" + str(self.it_step)
        self.next_node_idx = next_node_idx  # The next node to be simulated
        self.results = None

    def save(self):
        """Save the iteration file based on how many sims are left"""
        # prefix is one of "pre", "inter", "post"

        if self.results is None:
            new_prefix = "_pre"
        elif len(self.results) < len(self.nodes_scaled):
            new_prefix = "_inter"
        elif len(self.results) == len(self.nodes_scaled):
            new_prefix = "_post"

        now = str(datetime.now())
        file_name = f"{now + new_prefix}_sim_iteration_{self.it_step}"
        self.status = new_prefix
        with open(f"{self.file_tree.get_its_num_dir(0, self.it_step)}/{file_name}.pkl", 'wb') as f:
            pkl.dump(self, f)
        print(f"Iteration saved as {file_name}")
        return

    def give_shapes(self):
        return self.results, self.nodes_discrete, self.No, self.Cu, \
                self.num_div, self.sol_range, self.it_step

    def give_it_step(self):
        return self.it_step

    def give_file_name(self):
        return self.file_name

    def do_split(self, split_dirs, num_splits, max_num_sims=None):
        """
        Given a max number of simulations to allocate, a number of different
        jobs to split them between,and two directories representing the
        different file system storage locations, split the nodes
        among the jobs and save them in their respective directories.
        Then save update the iteration file and save it.
        """
        first_sim_idx = self.next_node_idx
        scaled_length = len(self.nodes_scaled)
        if (max_num_sims is None) \
                or (scaled_length - first_sim_idx <= max_num_sims):
            n_nodes_to_split = scaled_length - first_sim_idx
        else:
            n_nodes_to_split = max_num_sims


        print(f"{n_nodes_to_split} nodes to distribute among {num_splits} "
              f"jobs.\n")

        n_left = n_nodes_to_split
        for split in range(num_splits):
            d_nodes = math.ceil(n_nodes_to_split / num_splits)
            split_dir = split_dirs[split]

            if n_left <= 0:
                # No nodes left to distribute
                break
            elif n_left <= d_nodes:
                sub_scaled = self.nodes_scaled[first_sim_idx:, :]
            else:  #
                sub_scaled = self.nodes_scaled[first_sim_idx:first_sim_idx
                                                             + d_nodes, :]
            first_sim_idx += sub_scaled.shape[0]
            n_left -= sub_scaled.shape[0]

            # Pickle the subarray
            file_name = f"{split_dir}/split_{split}.pkl"
            with open(file_name, "wb") as f:
                pkl.dump(sub_scaled, f)
            print(f"{file_name} pickled array saved\n")
        print(f"split from {self.next_node_idx}")
        self.next_node_idx += n_nodes_to_split
        print(f"up to {self.next_node_idx}")
        self.save()
        return

    def update_results(self, np_arr):
        if self.results is None:
            self.results = np_arr
        else:
            self.results = np.vstack(self.results, np_arr)

