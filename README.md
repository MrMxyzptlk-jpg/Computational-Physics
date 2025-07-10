# Computational-Physics

## Guia 1: Chaotic systems
Used the following tools for analysis:
- Discrete Fourier transforms
- Poincare sections
- Fractals' visualization
The systems considered where the Logistic Equation ([ex2](Guia1/ex2)), Double Pendulum ([ex3](Guia1/ex3)), Pullen-Edmonds Hamiltonian ([ex4](Guia1/ex4)).

## Guia 2: Partial differential equations.
Solved the heat equation with three different methods:
- Forward Euler (explicit method)
- Backward Euler (implicit method)
- Crank-Nicolson
Considered null Dirichlet boundary condition in [ex1](Guia2/ex1/) and null Neumann boundary condition in [ex2](Guia2/ex2/).

## Guia 3: Random number generators (RNGs)
Introductory level concepts such as:
- Random walks ([ex2](Guia3/ex2/))
- Simple Monte Carlo integration ([ex3](Guia3/ex3/))
- N-dimensional ball Monte Carlo integration with gaussian sampling ([ex4](Guia3/ex4/))
Implemented three different RNG's: ran2, MZRAN and "Mersenne Twister" (all referenced can be found in the used modules).


## Guia 4: 2D Ising model
Calculated different observables and visualized the evolution, using the RNG's from the previous section.

## Guia 5: Molecular Dynamic
Simulated a many body system interacting through a two particle potential, chosen to be the Lennard-Jones potential. Implemented Molecular Dynamics (velocity-Verlet algorithm only) and also Monte Carlo (Metropolis algorithm only) to compare. Various observables and the evolution can be saved.

## Guia 6: Brownian Dynamics
Added to the program in the previous section a new formalism to simulate particles in a colloidal suspension (deterministic and stochastic contributions).

## Particle Dynamics
Added a non-local potential: Coulomb potential. Exploited the periodicity of the bulk simulation to calculate long range interaction using Ewald summation (see [details](Particle_Dynamics/README.md)). This code is the full modularization of the code in the previous two sections. It also allows for non-equivalent charges in the system.

