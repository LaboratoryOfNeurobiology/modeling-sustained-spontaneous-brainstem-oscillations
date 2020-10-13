import sys
import FileTree as ft
import numpy as np
import pickle as pkl
import os
from os import listdir
from os.path import isfile, join

if __name__ == "__main__":
    it_file_time = str(sys.argv[1])
    iteration = int(sys.argv[2])
    n_splits = int(sys.argv[3])
    tree = ft.FileTree()
    it_dir = tree.get_its_num_dir(0, iteration)
    it_dir_files = [f for f in listdir(it_dir) if isfile(join(it_dir, f))]
    results_dir_list = [f"{tree.get_split_results_dir(0, iteration, j)}" for j in range(n_splits)]
    size_total = 0
    ArrList = []
    it_file = None

    # Find correct iteration file
    for i in it_dir_files:
        if it_file_time in i:
            it_file = i
            break

    # Load iteration object
    if it_file is not None:
        with open(f"{it_dir}/{it_file}", "rb") as f:
            It = pkl.load(f)
    else:
        exit("Iteration file not found, exiting")

    print(It.results)

    # Gather results arrays into ArrList
    for direc in results_dir_list:
        results_files = [f for f in listdir(direc) if isfile(join(direc, f))]
        for file in results_files:
            with open(f"{direc}/{file}", "rb") as f:
                arr = pkl.load(f)
                size_total += arr.shape[0]
                ArrList.append(arr)
                #os.remove(f"{direc}/{file}")

    # Concatenate all results arrays and append them to the iteration's results
    result_arr = np.vstack(ArrList)
    It.update_results(result_arr)
    print(It.results)
    It.save()
    print(len(It.results))
