# Dragonball-style Heller Dynamics

A small Fortran code to test the feasability of Heller TGWP dynamics on an ab initio trajectory in the Dragonball format. 

## Workflow

The expected workflow:

1. read equilibrium geometry and Hessian
2. compute normal modes and normal mode rotation matrix
3. read the trajectory and compute the normal mode
4. evolve width using Hessian matrix
5. compute instanteneous correlation function
6. Fourier-transform

## TO DO

- [x] input reading
- [x] equilibrium reading
   - [x] Mass database
   - [x] Symbol matching
   - [x] actual reading
- [x] normal modes
- [x] trajectory into normal modes
- [x] evolve width
- [x] evolve $\Delta \gamma$
- [x] correlation function w/ MAPLE
- [x] compute initial correlation function
- [x] code correlation function
- [x] fourier
