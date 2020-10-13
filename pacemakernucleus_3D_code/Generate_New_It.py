import sys
import pickle as pkl
from FileTree import FileTree
from Iteration import Iteration
from New_iteration import new_iteration
from os import listdir
from os.path import isfile, join
import time


def gen_new_it(pre, it):
    # Get correct file path
    ft = FileTree()
    path = f"{ft.get_its_num_dir(0, it)}"
    only_files = [f for f in listdir(path) if isfile(join(path, f))]
    it_file = ""
    for i in only_files:
        if pre in i:
            it_file = i
            break
    # Load file buffer
    with open(f"{path}/{it_file}", "rb") as f:
        It = pkl.load(f)
    # Get shapes
    results, nodes_discrete, No, Cu, num_div, sol_range, it_step = \
        It.give_shapes()
    nodes_scaled, nodes_discrete, No, Cu = new_iteration(results,
                                                         nodes_discrete, No, Cu,
                                                         num_div, sol_range,
                                                         it_step)
    New_It = Iteration(It.sol_range, It.num_div, it_step + 1, nodes_scaled,
                       nodes_discrete, Cu, No)
    New_It.save()


if __name__ == "__main__":
    pre = sys.argv[1]
    iteration = 0
    start = time.time()
    gen_new_it(pre, iteration)
    end = time.time()
    print(
        f"{((end - start) / 60)} minutes to generate a new iteration"
        f" {int(iteration) + 1}")
