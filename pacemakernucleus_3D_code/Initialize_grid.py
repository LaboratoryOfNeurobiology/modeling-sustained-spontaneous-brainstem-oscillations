
from Initialize import initialize
import Iteration as it
import FileTree as ft

if __name__ == '__main__':
    current_step = 0
    sim_name = "network_simulations"
    # Range of solution
    sol_range = [[-150, -50], [0, 0.1], [0, 0.1]]
    # Number of steps along the associated coordinates, must be integers
    num_div = [8, 6, 6]
    # File structure
    file_tree = ft.FileTree()
    file_tree.create_tree()
    # Initialize grid
    nodes_scaled, nodes_discrete, Cu = initialize(sol_range, num_div)
    No = {}
    # Create and save iteration object
    initial_it = it.Iteration(sim_name, sol_range, num_div, current_step,
                              nodes_scaled, nodes_discrete, Cu, No)
    initial_it.save()
    # Exit
    print("initial iteration successfully saved, proceed with Split_iteration")
    # Done
