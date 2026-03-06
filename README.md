# FEM Beam Analysis (MATLAB + ANSYS)

## Overview
Finite element analysis of a multi-section beam using Euler–Bernoulli beam elements.  
The project includes a MATLAB implementation of the FEM formulation and validation using ANSYS simulation.

## Objectives
- Develop a beam finite element model
- Compute nodal displacements and rotations
- Calculate reaction forces and bending moments
- Perform modal analysis to determine natural frequencies
- Validate results using ANSYS

## Engineering Approach
- Derivation of the beam element stiffness matrix
- Assembly of the global stiffness matrix
- Application of boundary conditions
- Solution of the displacement system
- Eigenvalue analysis for natural frequencies

## Tools
- MATLAB
- ANSYS
- Finite Element Method (FEM)

## Results
- Deflection distribution along the beam
- Rotation and bending moment diagrams
- Natural frequencies of the structure

## Repository Contents
- MATLAB code for FEM beam analysis
- ANSYS simulation results
- Plots of structural response
- Project documentation

## MATLAB Code

The finite element formulation was implemented in MATLAB using Euler-Bernoulli beam elements.

The code includes:
- global stiffness and mass matrix assembly
- static solution for nodal displacements
- reaction force calculation
- deflection, rotation, and bending moment plotting
- midpoint deflection evaluation
- natural frequency extraction by solving the generalized eigenvalue problem

Main script:
- `code/fem_beam_analysis.m`
