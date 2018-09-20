# MonteCarloSimulatedAnnealing
Monte Carlo code of Simulated Annealing written in Fortran, which can be submitted to a cluster, which uses SLURM.
It is also possible to use the code without the use of a cluster by just deleting the sbatch command in thermo_script. The operating system assumed is a Linux distribution. The Fortran compiler used was GCC 4.9.3 and the OpenMPI version 1.8.8.

OBJECT OF SIMULATION:

The system under consideration consists of a 2D lattice (hexagonal or square), on which classical particles can reside. They can only be located on the lattice sites. The energy or Hamiltonian of the system is written as
\mathcal{H}=\sum_{i,j}V(i,j)n_in_j,
where n_i represents the particle occupation number and is either 0 or 1 and V(i,j) is -J for nearest neighbors and C*exp(-|r_i-r_j|r_s), where J and C are the parameters of the interaction. J represents the attraction between nearest neighbors and C the repulsion between all other pairs. The number of particles in the system is held constant. The boundary conditions are periodic. The Markov chain process consists of the destruction of one random particle and generation of another on one of the nearest neighbor sites of the destroyed one.
The code performs thermalization from a high temperature to a low temperature and is capable of measuring the energy of the system as well as the heat capacity. The initial configuration is  One can then use the Mathematica file to check results of the simulation.

SOFTWARE REQUIREMENTS:

- Linux operating system (Rocks 6.1)
- GNU Fortran compiler (GCC 4.9.3)
- OpenMPI (OpenMPI 1.8.8)
- SLURM (optional use on a cluster)

HOW TO RUN THE CODE:

First of all it is important to make a few modifications:

- In the file thermodynamics.f90 one needs to define the working directory of the program in the varible pwd (for example: pwd="/mydirectory/").
- In the file thermodynamics change << ENTER THE WORKING DIRECTORY >> to your working directory (for example: /mydirectory/)
- In the file thermo_script change << ENTER THE WORKING DIRECTORY >> to your working directory (for example: /mydirectory/)
- In the file CodeCheck.nb change << ENTER THE WORKING DIRECTORY >> to your working directory (keep in mind the Mathematica standar of directory names)

Afterwards, input the desired parameters in thermo_script and run it. The output should consist of four files starting with all_konf_, konf_, en_ and cv_. Then run the necessary code in CodeCheck.nb to view the results. A fair warning though, calculations for system sizes up to 30x30 sites are quite fast but get exponentially slower with increasing the size.

The uploaded files already contain an example of a frustrated system, where there is no crystallization. There is a glassy transition present.

PARAMETER DEFINITIONS:

answer : Set this to 0 unless you want to use a predefined particle configuration.

lattice : Options are "tri" or "sq" for triagular or square lattice.

tempz : The starting temperature of the system (high temperature).

tempk : The final temperature of the system (low temperature).

tempi : Temperature interval (should be small).

intname : Name your interaction type.

jay : The attraction magnitude. It ca also be repulsion if jay is negative.

lc : how many cs you want to simulate.

c : array of all the cs you want to simulate. one c is the magnitude of the repulsion between particles.

radij : The cutoff of the interaction (maximum is half the system size).

mrad : The range of the hop of a particle in Monte Carlo dynamics.

ls : r_s.

le : Number of system sizes to simulate.

side1 : One side of the system in number of lattice sites along the side.

side2 : Other side of the system in number of lattice sites along the side.

frac : The ration between the number of particles in the system and the number of lattice sites.

ni : The number of Monte Carlo sweeps (ni will get multiplied by the number of particles in the system) to perform at each temperature.

intv : Number of iterations (will also get multiplied by the number of particles in the system) to let pass between sampling the energy and specific heat in the system.

intt : Number of sweeps (-||-) before sampling begins.
