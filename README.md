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
- [ ] equilibrium reading
   - [ ] Mass database
   - [ ] Symbol matching
   - [ ] actual reading
- [ ] normal modes
- [ ] trajectory into normal modes
- [ ] evolve width
- [ ] correlation function w/ MAPLE
- [ ] code correlation function
- [ ] fourier
