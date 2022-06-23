# WaterClusters_PIMD
This code is designed to perform both classical (MD) or Path Integral Molecular Dynamics (PIMD) simulations at constant temperature of finite water clusters. For the water interactions, subroutines implementing the flexibles [SPC/F](https://doi.org/10.1103/PhysRevB.31.2643) and [qTIP4P/F](https://doi.org/10.1063/1.3167790) potentials are included. The canonical ensemble is efficiently sampled through the use of a [Langevin thermostat coupled to each degree of freedom](https://doi.org/10.1103/PhysRevE.75.056707). In addition, the energy and forces calculation of each replica (or bead) of the system is parallelized with OMP, which speeds up the code significantly. 

The structure of the program in Fortran modules was designed and implemented by Eduardo Hernández, in a similar way to the [Trocadero software](https://doi.org/10.1016/S0927-0256(03)00100-9). The main program and the subroutines needed to implement the path integral method were written and validated by Alfonso Gijón.

The code is not currently maintained, but any question can be asked to Alfonso Gijón (alfonso.gijon2@gmail.com). If you download and use the code for your work, please cite this repository and our [paper](https://doi.org/10.1039/D2CP01088G). 
---

```
@Article{Gijon2021,
  author ="Gijón, A. and Hernández, E. R.",
  title  ="Quantum simulations of neutral water clusters and singly-charged water cluster anions",
  journal  ="Phys. Chem. Chem. Phys.",
  year  ="2022",
  volume  ="24",
  issue  ="23",
  pages  ="14440-14451",
  publisher  ="The Royal Society of Chemistry",
  doi  ="10.1039/D2CP01088G",
  url  ="http://dx.doi.org/10.1039/D2CP01088G",
  abstract  ="We report a computational study of the structural and energetic properties of water clusters and singly-charged water cluster anions containing from 20 to 573 water molecules. We have used both a classical and a quantum description of the molecular degrees of freedom. Water intra and inter-molecular interactions have been modelled through the SPC/F model{,} while the water-excess electron interaction has been described via the well-known Turi–Borgis potential. We find that in general the quantum effects of the water degrees of freedom are small{,} but they do influence the cluster-size at which the excess electron stabilises inside the cluster{,} which occurs at smaller cluster sizes when quantum effects are taken into consideration."
}
```
