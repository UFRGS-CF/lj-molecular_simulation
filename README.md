# lj-molecular_simulation
Molecular Dynamics Simulation in C based on "Understanding Molecular Simulation From Algorithms to Applications" by Frenkel and Smit.


## MolecularSimulation.c

This is a molecular dynamics simulations based on Lennard-Jones potential. The paramaters for this simulation are number of particles (N), density (rho), initial temperature (t0), time-step (dt), end time (tmax) and equilibrium time (teq). The program use a Verlet Algorithm for the integration. It also computes the radial distribution and mean squared displacement (MSD). The output are 4 files:

 - energias.txt: contains the time (T), potential energy (U), kinetic energy (K) and total energy (E) per particle;
 - g_r.txt: contains the radius of the sphere analysed (r) and the radial distribution (g(r));
 - msd.txt: contains the time (T) and the mean squared displacemente (dr);
 - posi.txt: are a group of files each one containing the position of each particle in i 100 steps in simulation. This files are used in the program Jmol (http://jmol.sourceforge.net/) to visualize the particles evolution;
 
**Verify before running the program that in the same branch you have a folder called "podicoes". In this folder will be saved all the posi.txt files.**

## Graficos
This is a directory with python programs to plot simulation results. You can choose to plot them separately or in a single run (use "Graficos_UNI.py"). The graphs are: of the energy of the system, the MSD, the radial distribution, the temperature, the velocity autocorrelation (**still in progress**) and the velocity distribution. According to Frenkel's book, the energy graph must have the following behaviour: 

- the total energy must remain constant; 
- the kinetic and potential energy may vary initially, but they must oscillate around their equilibrium value near the end;

Any other behaviour means there is an error in the simulation.
