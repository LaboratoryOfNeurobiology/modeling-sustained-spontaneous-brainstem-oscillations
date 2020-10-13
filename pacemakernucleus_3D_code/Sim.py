

class Sim(object):
    """
    Represents a neuron simulation that has already been run.
        :param node_arr = (np.Array) gives parameters passed to simulation
        :param frequency = (real) the frequency of spontaneous oscillations
        :param all_cells = (List of PacemakerCell and RelayCell) Contains all cell objects which can be used
            for data analysis of network activity
        :param run_time = (float) total time the simulation took to execute.
        :param results_index = (int) index of this node's position in the array for this split.
    """
    def __init__(self, node_arr, frequency, all_cells, run_time, results_index):
        self.node = node_arr
        self.freq = frequency
        self.all_cells = all_cells
        self.run_time = run_time
        self.results_index = results_index

    def give_node(self):
        return self.node

    def give_freq(self):
        return self.freq

    def give_cells(self):
        return self.all_cells

    def give_run_time(self):
        return self.run_time

    def give_results_idx(self):
        return self.results_index