# PacemakerNucleus
Code storage for the models and simulations of the pacemaker nucleus and its sub-structures.

Files:
PN_Modeling.py contains all of the model specifications and functions for construction of networks using pacemaker and relay cells.

HalfGap.mod is an NMODL file containing specifications for the rectifying electrical gap junctions connecting cells in the network of the pacemaker nucleus. This must be in the working directory of the simulation script, or in a parent directory, and nrnivmodl must be run prior to simulation of any model employing this mod file.


Folders:
examples: contains example scripts for the use of the IntervalHalving method for parameter selection

src: contains the most recent up-to-date versions of the simulation scripts for models of the pacemaker nucleus and its substructures.

solutions: contains pickled solution files representing n-dimensional data for an interval halved simulation.

![Closed boundary in 2 dimensions [EK, gK] for sustained, spontaneous oscillations of an isolated pacemaker cell.](images/pacemaker_cell_func_Nx_15_Ny_20_Nit6.pdf)
