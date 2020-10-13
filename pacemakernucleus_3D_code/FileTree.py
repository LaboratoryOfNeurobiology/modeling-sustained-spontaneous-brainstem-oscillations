import os


class FileTree(object):
    """
    This class contains the paths and methods for structuring storage of data
    on the Discovery cluster.
    -------------------------------------------------------------------------
    ├── root1
│   └── iterations
│       ├── iteration_0
│       │   ├── 2020-06-09\ 14:32:37.496718_pre_sim_iteration_0.pkl
│       │   └── splits
│       │       ├── split_0
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_1
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_2
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_3
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_4
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_5
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_6
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       └── split_7
│       │           ├── pickled_simulation_objects
│       │           ├── results_dictionaries
│       │           └── tarred_simulation_objects
│       ├── iteration_1
│       │   └── splits
│       │       ├── split_0
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_1
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_2
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_3
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_4
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_5
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_6
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       └── split_7
│       │           ├── pickled_simulation_objects
│       │           ├── results_dictionaries
│       │           └── tarred_simulation_objects
│       ├── iteration_2
│       │   └── splits
│       │       ├── split_0
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_1
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_2
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_3
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_4
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_5
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_6
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       └── split_7
│       │           ├── pickled_simulation_objects
│       │           ├── results_dictionaries
│       │           └── tarred_simulation_objects
│       ├── iteration_3
│       │   └── splits
│       │       ├── split_0
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_1
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_2
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_3
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_4
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_5
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_6
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       └── split_7
│       │           ├── pickled_simulation_objects
│       │           ├── results_dictionaries
│       │           └── tarred_simulation_objects
│       ├── iteration_4
│       │   └── splits
│       │       ├── split_0
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_1
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_2
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_3
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_4
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_5
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       ├── split_6
│       │       │   ├── pickled_simulation_objects
│       │       │   ├── results_dictionaries
│       │       │   └── tarred_simulation_objects
│       │       └── split_7
│       │           ├── pickled_simulation_objects
│       │           ├── results_dictionaries
│       │           └── tarred_simulation_objects
│       └── iteration_5
│           └── splits
│               ├── split_0
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_1
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_2
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_3
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_4
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_5
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               ├── split_6
│               │   ├── pickled_simulation_objects
│               │   ├── results_dictionaries
│               │   └── tarred_simulation_objects
│               └── split_7
│                   ├── pickled_simulation_objects
│                   ├── results_dictionaries
│                   └── tarred_simulation_objects
└── root2
    └── iterations
        ├── iteration_0
        │   └── splits
        │       ├── split_10
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_11
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_12
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_13
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_14
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_15
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_8
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       └── split_9
        │           ├── pickled_simulation_objects
        │           ├── results_dictionaries
        │           └── tarred_simulation_objects
        ├── iteration_1
        │   └── splits
        │       ├── split_10
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_11
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_12
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_13
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_14
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_15
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_8
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       └── split_9
        │           ├── pickled_simulation_objects
        │           ├── results_dictionaries
        │           └── tarred_simulation_objects
        ├── iteration_2
        │   └── splits
        │       ├── split_10
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_11
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_12
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_13
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_14
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_15
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_8
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       └── split_9
        │           ├── pickled_simulation_objects
        │           ├── results_dictionaries
        │           └── tarred_simulation_objects
        ├── iteration_3
        │   └── splits
        │       ├── split_10
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_11
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_12
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_13
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_14
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_15
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_8
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       └── split_9
        │           ├── pickled_simulation_objects
        │           ├── results_dictionaries
        │           └── tarred_simulation_objects
        ├── iteration_4
        │   └── splits
        │       ├── split_10
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_11
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_12
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_13
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_14
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_15
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       ├── split_8
        │       │   ├── pickled_simulation_objects
        │       │   ├── results_dictionaries
        │       │   └── tarred_simulation_objects
        │       └── split_9
        │           ├── pickled_simulation_objects
        │           ├── results_dictionaries
        │           └── tarred_simulation_objects
        └── iteration_5
            └── splits
                ├── split_10
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_11
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_12
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_13
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_14
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_15
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                ├── split_8
                │   ├── pickled_simulation_objects
                │   ├── results_dictionaries
                │   └── tarred_simulation_objects
                └── split_9
                    ├── pickled_simulation_objects
                    ├── results_dictionaries
                    └── tarred_simulation_objects
    """

    def __init__(self):
        # Root directories
        self.root_0 = "/scratch/hartman.da/network_simulations_3D"
        self.root_1 = "/scratch/hartman.da/network_simulations_3D"
        # One dir below root
        self.its_dir = "/iterations"
        self.its_prefix = "/iteration_"
        self.frequencies_text = "/aggregated_frequencies"
        # Two dirs below root
        self.pickled_dicts_dir = "/pickled_simulation_objects"
        self.tarred_simulations_dir = "/tarred_simulation_objects"
        self.splits_dir = "/splits"
        self.splits_prefix = "/split_"
        self.results_dir = "/results_array"
        self.results_prefix = "/result_"

    # Utility Functions
    def create_tree(self):
        try:
            for i in range(1):  # For each root
                os.mkdir(self.get_root_dir(i))
                os.mkdir(self.get_its_dir(i))
                for j in range(7):  # For each iteration
                    if i == 0:
                        k_low = 0
                        k_rng = 128
                    os.mkdir(self.get_its_num_dir(i, j))
                    os.mkdir(self.get_splits_dir(i, j))
                    for k in range(k_low, k_rng): # For each split
                        os.mkdir(self.get_splits_num_dir(i, j, k))
                        os.mkdir(self.get_split_results_dir(i, j, k))
                        os.mkdir(self.get_pickled_sims_dir(i, j, k))
                        os.mkdir(self.get_tarred_sim_dir(i, j, k))

        except OSError:
            print("File Structure has already been created, and should not be overridden. Exiting.")
            pass

    def create_follow_tree(self, root):
        try:
            os.mkdir(self.get_root_dir(root))
            os.mkdir(self.get_follow_dir(root))
            for j in range(7,27):  # For each iteration
                os.mkdir(self.get_its_num_dir(root, j, True))
                os.mkdir(self.get_splits_dir(root, j, True))
                for k in range(128):  # For each split
                    os.mkdir(self.get_splits_num_dir(root, j, k, True))
                    os.mkdir(self.get_split_results_dir(root, j, k, True))
                    os.mkdir(self.get_pickled_sims_dir(root, j, k, True))
                    os.mkdir(self.get_tarred_sim_dir(root, j, k, True))
        except OSError:
            print(
                "File Structure has already been created, and should not be "
                "overridden. Exiting.")
            pass

    def get_root_dir(self, root):
        if root == 0:
            return self.root_0
        else:
            return self.root_1

    # In network_simulations/
    def get_freqs_file(self, root):
        return self.get_root_dir(root) + "/aggregated_frequencies.txt"

    def get_its_dir(self, root):
        return self.get_root_dir(root) + self.its_dir

    def get_follow_dir(self, root):
        return self.get_root_dir(root) + self.follow_dir

    # In Iterations/
    def get_its_num_dir(self, root, it, follow=False):
        if follow:
            return self.get_follow_dir(root) + self.its_prefix + str(it)
        else:
            return self.get_its_dir(root) + self.its_prefix + str(it)

    # In Iteration_n/
    def get_splits_dir(self, root, it, follow=False):
        return self.get_its_num_dir(root, it, follow) + self.splits_dir

    def get_splits_num_dir(self, root, it, split, follow=False):
        return self.get_splits_dir(root, it, follow) \
               + self.splits_prefix + str(split)

    # In splits/
    def get_split_results_dir(self, root, it, split, follow=False):
        return self.get_splits_num_dir(root, it, split, follow) + \
               self.results_dir

    def get_pickled_sims_dir(self, root, it, split, follow=False):
        return self.get_splits_num_dir(root, it, split, follow) \
               + self.pickled_dicts_dir

    def get_tarred_sim_dir(self, root, it, split, follow=False):
        return self.get_splits_num_dir(root, it, split, follow) \
               + self.tarred_simulations_dir
