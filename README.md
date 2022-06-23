# WaterClusters_PIMD
This code is designed to perform both classical (MD) or Path Integral Molecular Dynamics (PIMD) simulations at constant temperature of finite water clusters. For the water interactions, subroutines implementing the flexibles SPC/F and qTIP4P/F potentials are included. The canonical ensemble is efficienlty sampled through the use of a Langevin thermostat coupled to each degree of freedom. In addition, the energy and forces calculation of each replica (or bead) is parallelized with OMP, which speeds up the code significantly. 

The structure of the program in Fortran modules was designed and implemented by Eduardo Hernández, in a similar way to the Trocadero software (R. Rurali, E. Hernández, Computational Materials Science 28 (2003) 85–106). The main program and the subroutines needed to implement the path integral method were written by Alfonso Gijón.

The code is not currently maintained, but any question can be asked to Alfonso Gijón (alfonso.gijon2@gmail.com). If you download and use the code for your work, please cite our [paper](https://doi.org/10.1039/D2CP01088G) 
---
