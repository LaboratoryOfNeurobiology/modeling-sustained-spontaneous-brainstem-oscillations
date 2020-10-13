import Iteration_Follow as it
import FileTree as ft
import sys
import pickle as pkl
from os import listdir
from os.path import isfile, join
from Initialize_follow import initialize_follow

if __name__ == '__main__':
    "=================================="
    sim_name = "network_simulations_3D"
    it_step = 6
    time_5_prefix = sys.argv[1]
    time_6_prefix = sys.argv[2]
    isFollowing = True
    "=================================="
    if isFollowing:
        root = 0
    else:
        root = 0
    # Range of solution
    file_tree = ft.FileTree()
    it_dir_5 = file_tree.get_its_num_dir(0, it_step - 1)
    it_dir_6 = file_tree.get_its_num_dir(0, it_step)
    only_files_5 = [f for f in listdir(it_dir_5) if isfile(join(it_dir_5, f))]
    only_files_6 = [f for f in listdir(it_dir_6) if isfile(join(it_dir_6, f))]
    it_file_5 = ""
    for i in only_files_5:
        if time_5_prefix in i:
            it_file_5 = i
            break
    # Find the correct save based on `time`
    it_file_6 = ""
    for i in only_files_6:
        if time_6_prefix in i:
            it_file_6 = i
            break

    with open(f"{it_dir_5}/{it_file_5}", "rb") as f:
        It_obj_5 = pkl.load(f)

    with open(f"{it_dir_6}/{it_file_6}", "rb") as f:
        It_obj_6 = pkl.load(f)

    results_5, nodes_discrete_5, No_5, Cu_5, num_div_5, sol_range_5, \
    it_step_5 = It_obj_5.give_shapes()
    results_6, nodes_discrete_6, No_6, Cu_6, num_div_6, sol_range_6, \
    it_step_6 = It_obj_6.give_shapes()

    nodes_scaled, nodes_discrete, Cu_add, Cu_br \
        = initialize_follow(No_6, Cu_5, sol_range_5, num_div_5, it_step_5)
    new_it = it.Iteration(sim_name, sol_range_5, num_div_5, it_step_5,
                          nodes_scaled, nodes_discrete, Cu_br, Cu_add, No_6, isFollowing=False)
    print(f"{len(nodes_scaled)} nodes to simulate in this follow step")
    new_it.save()
    # Exit
    print("initial iteration successfully saved, proceed with Split_iteration")
    # Done

