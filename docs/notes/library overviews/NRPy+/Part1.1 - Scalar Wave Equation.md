# Overview
- In this notebook we want to solve the scalar wave equation.
  $$
  \partial_t^2 u (x^\mu) = c^2 \nabla^2 u(x^\mu)
  $$
- We make two coupled first order partial differential equations from it
  $$
  \begin{align}
  \partial_t u &= v\\
  \partial_t v &= c^2\nabla^2 u
  \end{align}
  $$
- Once we have set our initial conditions we evolve the system forward in time.
# Method of lines
 The method of lines would be used because:
 1. Enables us to handle the spatial derivatives of an initial value problem PDE using standard finite difference methods
 2. Enables us to handle the temporal derivative of an initial value problem PDE using standard strategies for solving ordinary differential equations
 ***All is true as long as the initial value problem PDE can be written in the form:***
 $$
 \partial_t \vec f = M\cdot \vec f
 $$
 So for temporal derivatives we implement Runge Kutta 5 as our standard approach.

# Basic Algorithm
1. Set up the numerical grid and free parameters
2. Using the exact solutions we allocate memory for grid functions including temporary storage for the RK4 steps
3. Set the grid-function equal to the initial data
4. Evolve:
	1. At each step evaluate the RHS of the equations (spatial parts)
	2. Apply the boundary conditions

# Discussion
This notebook provides very small detail about the reason behind this kind of implementation. To the eye of a computational physicist it might look fine but for learning, the reason behind grids, the exact solution implementation, etc. might look vague.
 