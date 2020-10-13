import sys
import pickle as pkl
import FileTree as ft
from os import listdir
from os.path import isfile, join


def split(t, iteration, max_sims, splits):
    sim_name = "network_simulations_3D"
    it_to_split = iteration
    prefix = "_pre"
    prefix_to_split = f"{t}{prefix}_sim_iteration_{it_to_split}"
    n_splits = splits
    max_sims_to_split = max_sims

    file_tree = ft.FileTree()
    # Get path to iteration save
    it_dir_0 = file_tree.get_its_num_dir(0, it_to_split)
    only_files = [f for f in listdir(it_dir_0) if isfile(join(it_dir_0, f))]
    # Find the correct save based on `time`
    it_file = ""
    for i in only_files:
        if prefix_to_split in i:
            it_file = i
            break
    print(it_file)
    # Load iteration
    split_paths = [file_tree.get_splits_num_dir(0, it_to_split, j)
                   for j in range(n_splits)]
    """split_paths += [file_tree.get_splits_num_dir((1, it_to_split, j)
                    for j in range(64, 83))]"""
    with open(f"{it_dir_0}/{it_file}", 'rb') as file_handle:
        It = pkl.load(file_handle)
        # Split Iteration object into different files
        It.do_split(split_paths, n_splits, max_sims_to_split)
    print("It split completed successfully.")


if __name__ == '__main__':
    # Split setup
    time = sys.argv[1]
    iteration = int(sys.argv[2])
    n_splits = int(sys.argv[3])
    max_sims = None
    split(time, iteration, max_sims, n_splits)
