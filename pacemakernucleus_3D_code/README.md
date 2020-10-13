Northeastern's Discovery Supercluster is being used to run the simulations in the 6-dimensional parameter space 
[EK, ENa, gNar, gKr, gNap, gKp]. Exploration of this parameter space is achieved with the interval-halving method.
This method can take in an arbitrary number of dimensions, their ranges, and an initial granularity, to more efficiently
explore the manifolds representing the equilibrium of the system.

Exploration of the pacemaker nucleus's parameter space, whose cells are modeled based on the Hodgkin-Huxley
equations for neuronal signaling, requires 6 iterations of Interval Halving in order to interpolate the surface
of the manifold. 
We explore this system across parameter space ranges that encompass - to our knowledge - biologically plausible values. 
The initial iteration requires simulation of several thousand points in the 6-dimensional parameter space. Successive 
iterations require orders or magnitude more simulations to be completed than the previous. 
Completion of all simulations across all iterations requires parallelizing as many simulations as possible. All
simulations are independant of each other which means that a batch style parallelization schema, employing 
Neuron's ParallelContext() module, and OpenMPI, is sufficient to complete the task. 

We plan to allocate 2048 cpu cores split across 40 different nodes on the Discovery Supercluster, where each node
runs an equal portion of the simulations for that iteration. This will enable us to run 2048 simulations in parallel 
where the current worst-case run-time for a single simulation is 9 minutes. 

We must create a modular workflow, employing several different python and bash scripts to do the following:
    1. Create the initial iteration's parameter-space grid
    2. Store data required for simulation of this grid in an iteration object and pickle it
    3. Load pickled iteration object and distribute the nodes to be simulated into 40 different "splits" pickle files
    4. Generic script takes in a single split pickle file and simulates in parallel, that file's nodes,
    then stores the frequency data into a new pickle file where nodes are matched with frequencies.
    5. Bash script employs slurm to allocate resources on Discovery and submit array jobs with
        mpirun -n 52 python genericscript.py file.pkl 
       where all pkl files are "splits" stored in an array
    6. Parse all simulated "splits", reconstruct the entire grid with frequencies and 
        generate the next iteration's grid. Store the grid into a new Iteration object.
    7. Pass new Iteration to step 3 and continue the cycle until all 6 iterations are done.