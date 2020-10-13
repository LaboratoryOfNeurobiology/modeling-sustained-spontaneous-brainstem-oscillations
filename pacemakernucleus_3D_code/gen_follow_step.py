import sys
import pickle as pkl
from FileTree import FileTree
from Iteration_Follow import Iteration
from os import listdir
from os.path import isfile, join
from New_follow_step import new_follow_step
import time


def gen_new_it(pre, it, isFollowing):
    if isFollowing:
        root = 0
    else:
        root = 0
    ft = FileTree()
    path = f"{ft.get_its_num_dir(root, it)}"
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
    results, nodes_discrete, No, Cu, num_div, sol_range, it_step = It.give_shapes()
    # Updating node dictionary before generating new follow_step
    #print(Cu)
    for i in range(len(nodes_discrete)):
        tup = nodes_discrete[i]
        No[tup] = results[i, -1]

    Cu_add = It.give_cu_add()
    nodes_scaled, nodes_discrete, Cu_add, Cu_br = \
        new_follow_step(No, Cu_add, Cu, sol_range, num_div, it_step)
    if len(nodes_discrete) > 0:
        New_It = Iteration(sim_name, sol_range, num_div,
                           it_step, nodes_scaled, nodes_discrete, Cu_br, Cu_add,
                           No)
    else:
        New_It = Iteration(sim_name, sol_range, num_div,
                           it_step, nodes_scaled, nodes_discrete, Cu, Cu_add,
                           No)
        print("All done final following step iteration saved."
              " Proceed to data analysis")
    New_It.save()
    print(f"number of nodes in next follow step = {len(nodes_discrete)}\n")


if __name__ == "__main__":
    "========================================"
    sim_name = "network_simulations_3D"
    isFollowing = True
    "========================================="
    pre = sys.argv[1]
    iteration = int(sys.argv[2])
    start = time.time()

    gen_new_it(pre, iteration, isFollowing)

    end = time.time()
    print(
        f"{((end - start) / 60)} minutes to generate a new iteration"
        f" {int(iteration) + 1}")

