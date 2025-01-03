<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Generating C Code for the Scalar Wave Equation in Cartesian Coordinates

## Authors: Zach Etienne & Thiago Assumpção
### Formatting improvements courtesy Brandon Clark

## This module generates the C Code for the Scalarwave in Cartesian coordinates and sets up either a monochromatic plane wave or spherical Gaussian [Initial Data](https://en.wikipedia.org/wiki/Initial_value_problem). The module emphasizes the use of the Method of Lines, thus transforming the partial differential equation problem into a more manageable ordinary differential equation problem.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented below ([right-hand-side expressions](#code_validation1); [initial data expressions](#code_validation2)). In addition, all expressions have been validated against a trusted code (the [original SENR/NRPy+ code](https://bitbucket.org/zach_etienne/nrpy)).

### NRPy+ Source Code for this module: 
* [ScalarWave/ScalarWave_RHSs.py](../edit/ScalarWave/ScalarWave_RHSs.py)
* [ScalarWave/InitialData.py](../edit/ScalarWave/InitialData.py)

## Introduction: 
### Problem Statement

We wish to numerically solve the scalar wave equation as an [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem) in Cartesian coordinates:
$$\partial_t^2 u = c^2 \nabla^2 u \text{,}$$
where $u$ (the amplitude of the wave) is a function of time and space: $u = u(t,x,y,...)$ (spatial dimension as-yet unspecified) and $c$ is the wave speed, subject to some initial condition

$$u(0,x,y,...) = f(x,y,...)$$

and suitable spatial boundary conditions.

As described in the next section, we will find it quite useful to define
$$v(t,x,y,...) = \partial_t u(t,x,y,...).$$

In this way, the second-order PDE is reduced to a set of two coupled first-order PDEs

\begin{align}
\partial_t u &= v \\
\partial_t v &= c^2 \nabla^2 u.
\end{align}

We will use NRPy+ to generate efficient C codes capable of generating both initial data $u(0,x,y,...) = f(x,y,...)$; $v(0,x,y,...)=g(x,y,...)$, as well as finite-difference expressions for the right-hand sides of the above expressions. These expressions are needed within the *Method of Lines* to "integrate" the solution forward in time.

### The Method of Lines

Once we have initial data, we "evolve it forward in time", using the [Method of Lines](https://reference.wolfram.com/language/tutorial/NDSolveMethodOfLines.html). In short, the Method of Lines enables us to handle 
1. the **spatial derivatives** of an initial value problem PDE using **standard finite difference approaches**, and
2. the **temporal derivatives** of an initial value problem PDE using **standard strategies for solving ordinary differential equations (ODEs)**, so long as the initial value problem PDE can be written in the form
$$\partial_t \vec{f} = \mathbf{M}\ \vec{f},$$
where $\mathbf{M}$ is an $N\times N$ matrix filled with differential operators that act on the $N$-element column vector $\vec{f}$. $\mathbf{M}$ may not contain $t$ or time derivatives explicitly; only *spatial* partial derivatives are allowed to appear inside $\mathbf{M}$. The scalar wave equation as written in the [previous module](Tutorial-ScalarWave.ipynb),
\begin{equation}
\partial_t 
\begin{bmatrix}
u \\
v 
\end{bmatrix}=
\begin{bmatrix}
0 & 1 \\
c^2 \nabla^2  & 0 
\end{bmatrix}
\begin{bmatrix}
u \\
v 
\end{bmatrix},
\end{equation}
satisfies this requirement. 

Thus we can treat the spatial derivatives $\nabla^2 u$ of the scalar wave equation  using **standard finite-difference approaches**, and the temporal derivatives $\partial_t u$ and $\partial_t v$ using **standard approaches for solving ODEs**. In [the next module](Tutorial-Start_to_Finish-ScalarWave.ipynb), we will apply the highly robust [explicit Runge-Kutta fourth-order scheme](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) (RK4), used widely for numerically solving ODEs, to "march" (integrate) the solution vector $\vec{f}$ forward in time from its initial value ("initial data").

### Basic Algorithm

The basic algorithm for solving the scalar wave equation [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem), based on the Method of Lines (see section above), is outlined below, with NRPy+-based components highlighted in <font color='green'>green</font>. We will review how NRPy+ generates these core components in this module.

1. Allocate memory for gridfunctions, including temporary storage for the RK4 time integration.
1. <font color='green'>Set gridfunction values to initial data.</font>
1. Evolve the system forward in time using RK4 time integration. At each RK4 substep, do the following:
    1. <font color='green'>Evaluate scalar wave RHS expressions.</font>
    1. Apply boundary conditions.

**We refer to the right-hand side of the equation $\partial_t \vec{f} = \mathbf{M}\ \vec{f}$ as the RHS. In this case, we refer to the $\mathbf{M}\ \vec{f}$ as the "scalar wave RHSs".** In the following sections we will 

1. use NRPy+ to cast the scalar wave RHS expressions -- in finite difference form -- into highly efficient C code, 
    1. first in one spatial dimension with fourth-order finite differences,
    1. and then in three spatial dimensions with tenth-order finite differences; we will
1. use NRPy+ to generate monochromatic plane-wave initial data for the scalar wave equation, where the wave propagates in an arbitrary direction.

As for the $\nabla^2 u$ term, spatial derivatives are handled in NRPy+ via [finite differencing](https://en.wikipedia.org/wiki/Finite_difference).

We will sample the solution $\{u,v\}$ at discrete, uniformly-sampled points in space and time. For simplicity, let's assume that we consider the wave equation in one spatial dimension. Then the solution at any sampled point in space and time is given by
$$u^n_i = u(t_n,x_i) = u(t_0 + n \Delta t, x_0 + i \Delta x),$$
where $\Delta t$ and $\Delta x$ represent the temporal and spatial resolution, respectively. $v^n_i$ is sampled at the same points in space and time.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

1. [Step 1](#initializenrpy): Initialize core NRPy+ modules
1. [Step 2](#rhss1d): Scalar Wave RHSs in One Spatial Dimension, Fourth-Order Finite Differencing
    1. [Step 2.a](#ccode1d): C-code output example: Scalar wave RHSs with 4th order finite difference stencils
1. [Step 3](#rhss3d): Scalar Wave RHSs in Three Spatial Dimensions
    1. [Step 3.a](#code_validation1): Code Validation against `ScalarWave.ScalarWave_RHSs` NRPy+ module
    1. [Step 3.b](#ccode3d): C-code output example: Scalar wave RHSs with 10th order finite difference stencils and SIMD enabled
1. [Step 4](#id): Setting up Initial Data for the Scalar Wave Equation
    1. [Step 4.a](#planewave): The Monochromatic Plane-Wave Solution
    1. [Step 4.b](#sphericalgaussian): The Spherical Gaussian Solution (*Courtesy Thiago Assumpção*)
1. [Step 5](#code_validation2): Code Validation against `ScalarWave.InitialData` NRPy+ module
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from NRPy+:


```python
# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
from outputC import lhrh         # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
```

<a id='rhss1d'></a>

# Step 2: Scalar Wave RHSs in One Spatial Dimension \[Back to [top](#toc)\]
$$\label{rhss1d}$$

To minimize complication, we will first restrict ourselves to solving the wave equation in one spatial dimension, so
\begin{align}
\partial_t u &= v \\
\partial_t v &= c^2 \nabla^2 u \\
 &= c^2 \partial_x^2 u.
\end{align}

We will construct SymPy expressions of the right-hand sides of $u$ and $v$ using [NRPy+ finite-difference notation](Tutorial-Finite_Difference_Derivatives.ipynb) to represent the derivative, so that finite-difference C-code kernels can be easily constructed.

Extension of this operator to higher spatial dimensions when using NRPy+ is straightforward, as we will see below.


```python
# Step P2: Define the C parameter wavespeed. The `wavespeed`
#          variable is a proper SymPy variable, so it can be
#          used in below expressions. In the C code, it acts
#          just like a usual parameter, whose value is
#          specified in the parameter file.
thismodule = "ScalarWave"
wavespeed = par.Cparameters("REAL",thismodule,"wavespeed", 1.0)

# Step 1: Set the spatial dimension parameter, and then read
#         the parameter as DIM.
par.set_parval_from_str("grid::DIM", 1)
DIM = par.parval_from_str("grid::DIM")

# Step 2: Register gridfunctions that are needed as input
#         to the scalar wave RHS expressions.
uu, vv = gri.register_gridfunctions("EVOL", ["uu", "vv"])

# Step 3: Declare the rank-2 indexed expression \partial_{ij} u,
#         which is symmetric about interchange of indices i and j.
#         Derivative variables like these must have an underscore
#         in them, so the finite difference module can parse the
#         variable name properly.
uu_dDD = ixp.declarerank2("uu_dDD", "sym01")

# Step 4: Define right-hand sides for the evolution.
uu_rhs = vv
vv_rhs = 0
for i in range(DIM):
    vv_rhs += wavespeed*wavespeed*uu_dDD[i][i]

vv_rhs = sp.simplify(vv_rhs)
```

<a id='ccode1d'></a>

## Step 2.a: C-code output example: Scalar wave RHSs with 4th order finite difference stencils \[Back to [top](#toc)\]
$$\label{ccode1d}$$

As was discussed in [the finite difference section of the tutorial](Tutorial-Finite_Difference_Derivatives.ipynb), NRPy+ approximates derivatives using [finite-difference methods](https://en.wikipedia.org/wiki/Finite_difference),  the second-order derivative $\partial_x^2$ accurate to fourth-order in uniform grid spacing $\Delta x$ (from fitting the unique 4th-degree polynomial to 5 sample points of $u$) is given by
\begin{equation}
\left[\partial_x^2 u(t,x)\right]_j = \frac{1}{(\Delta x)^2}
\left(
-\frac{1}{12} \left(u_{j+2} + u_{j-2}\right) 
+ \frac{4}{3}  \left(u_{j+1} + u_{j-1}\right)
- \frac{5}{2} u_j \right)
+ \mathcal{O}\left((\Delta x)^4\right).
\end{equation}


```python
# Step 5: Set the finite differencing order to 4.
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

# Step 6: Generate C code for scalarwave evolution equations,
#         print output to the screen (standard out, or stdout).
fin.FD_outputC("stdout",
               [lhrh(lhs=gri.gfaccess("rhs_gfs", "uu"), rhs=uu_rhs),
                lhrh(lhs=gri.gfaccess("rhs_gfs", "vv"), rhs=vv_rhs)])
```

    {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
       */
      /*
       *  Original SymPy expression:
       *  "const double uu_dDD00 = invdx0**2*(-5*uu/2 + 4*uu_i0m1/3 - uu_i0m2/12 + 4*uu_i0p1/3 - uu_i0p2/12)"
       */
      const double uu_i0m2 = in_gfs[IDX2S(UUGF, i0-2)];
      const double uu_i0m1 = in_gfs[IDX2S(UUGF, i0-1)];
      const double uu = in_gfs[IDX2S(UUGF, i0)];
      const double uu_i0p1 = in_gfs[IDX2S(UUGF, i0+1)];
      const double uu_i0p2 = in_gfs[IDX2S(UUGF, i0+2)];
      const double vv = in_gfs[IDX2S(VVGF, i0)];
      const double FDPart1_Rational_5_2 = 5.0/2.0;
      const double FDPart1_Rational_1_12 = 1.0/12.0;
      const double FDPart1_Rational_4_3 = 4.0/3.0;
      const double uu_dDD00 = ((invdx0)*(invdx0))*(FDPart1_Rational_1_12*(-uu_i0m2 - uu_i0p2) + FDPart1_Rational_4_3*(uu_i0m1 + uu_i0p1) - FDPart1_Rational_5_2*uu);
      /*
       * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
       */
      /*
       *  Original SymPy expressions:
       *  "[rhs_gfs[IDX2S(UUGF, i0)] = vv,
       *    rhs_gfs[IDX2S(VVGF, i0)] = uu_dDD00*wavespeed**2]"
       */
      rhs_gfs[IDX2S(UUGF, i0)] = vv;
      rhs_gfs[IDX2S(VVGF, i0)] = uu_dDD00*((wavespeed)*(wavespeed));
    }


**Success!** Notice that indeed NRPy+ was able to compute the spatial derivative operator,
\begin{equation}
\left[\partial_x^2 u(t,x)\right]_j \approx \frac{1}{(\Delta x)^2}
\left(
-\frac{1}{12} \left(u_{j+2} + u_{j-2}\right) 
+ \frac{4}{3}  \left(u_{j+1} + u_{j-1}\right)
- \frac{5}{2} u_j \right),
\end{equation}
correctly (it's easier to read in the "Original SymPy expressions" comment block at the top of the C output).

As NRPy+ is designed to generate codes in arbitrary coordinate systems, instead of sticking with Cartesian notation for 3D coordinates, $x,y,z$, we instead adopt $x_0,x_1,x_2$ for our coordinate labels. Thus you will notice the appearance of `invdx0`$=1/\Delta x_0$, where $\Delta x_0$ is the (uniform) grid spacing in the zeroth or $x_0$ direction. In this case, $x_0$ represents the $x$ direction.

<a id='rhss3d'></a>

# Step 3: Scalar Wave RHSs in Three Spatial Dimensions \[Back to [top](#toc)\]
$$\label{rhss3d}$$

Let's next repeat the same process, only this time for the scalar wave equation in **3 spatial dimensions** (3D).


```python
# Step 1: Define the C parameter wavespeed. The `wavespeed`
#         variable is a proper SymPy variable, so it can be
#         used in below expressions. In the C code, it acts
#         just like a usual parameter, whose value is
#         specified in the parameter file.
wavespeed = par.Cparameters("REAL", thismodule, "wavespeed", 1.0)

# Step 2: Set the spatial dimension parameter
#         to *FOUR* this time, and then read
#         the parameter as DIM.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

# Step 3a: Reset gridfunctions registered in 1D case above,
#          to avoid NRPy+ throwing an error about double-
#          registering gridfunctions, which is not allowed.
gri.glb_gridfcs_list = []

# Step 3b: Register gridfunctions that are needed as input
#          to the scalar wave RHS expressions.
uu, vv = gri.register_gridfunctions("EVOL", ["uu", "vv"])

# Step 4: Declare the rank-2 indexed expression \partial_{ij} u,
#         which is symmetric about interchange of indices i and j
#         Derivative variables like these must have an underscore
#         in them, so the finite difference module can parse the
#         variable name properly.
uu_dDD = ixp.declarerank2("uu_dDD", "sym01")

# Step 5: Define right-hand sides for the evolution.
uu_rhs = vv
vv_rhs = 0
for i in range(DIM):
    vv_rhs += wavespeed*wavespeed*uu_dDD[i][i]

# Step 6: Simplify the expression for c^2 \nabla^2 u (a.k.a., vv_rhs):
vv_rhs = sp.simplify(vv_rhs)
```

<a id='code_validation1'></a>

## Step 3.a:  Validate SymPy expressions against `ScalarWave.ScalarWave_RHSs` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation1}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the three-spatial-dimension Scalar Wave equation (i.e., `uu_rhs` and `vv_rhs`) between

1. this tutorial and 
2. the [NRPy+ ScalarWave.ScalarWave_RHSs](../edit/ScalarWave/ScalarWave_RHSs.py) module.


```python
# Step 10: We already have SymPy expressions for uu_rhs and vv_rhs in
#          terms of other SymPy variables. Even if we reset the list
#          of NRPy+ gridfunctions, these *SymPy* expressions for
#          uu_rhs and vv_rhs *will remain unaffected*.
#
#          Here, we will use the above-defined uu_rhs and vv_rhs to
#          validate against the same expressions in the
#          ScalarWave/ScalarWave_RHSs.py module,
#          to ensure consistency between this tutorial
#          (historically speaking, the tutorial was written first)
#          and the ScalarWave_RHSs.py module itself.
#
# Reset the list of gridfunctions, as registering a gridfunction
#   twice will spawn an error.
gri.glb_gridfcs_list = []

# Step 11: Call the ScalarWave_RHSs() function from within the
#         ScalarWave/ScalarWave_RHSs.py module,
#         which should do exactly the same as in Steps 1-10 above.
import ScalarWave.ScalarWave_RHSs as swrhs
swrhs.ScalarWave_RHSs()

# Step 12: Consistency check between the tutorial notebook above
#         and the ScalarWave_RHSs() function from within the
#         ScalarWave/ScalarWave_RHSs.py module.
print("Consistency check between ScalarWave tutorial and NRPy+ module:")
print("uu_rhs - swrhs.uu_rhs = "+str(sp.simplify(uu_rhs - swrhs.uu_rhs))+"\t\t (should be zero)")
print("vv_rhs - swrhs.vv_rhs = "+str(sp.simplify(vv_rhs - swrhs.vv_rhs))+"\t\t (should be zero)")
```

    Consistency check between ScalarWave tutorial and NRPy+ module:
    uu_rhs - swrhs.uu_rhs = 0		 (should be zero)
    vv_rhs - swrhs.vv_rhs = 0		 (should be zero)


<a id='ccode3d'></a>

## Step 3.b: C-code output example: Scalar wave RHSs with 10th order finite difference stencils and SIMD enabled \[Back to [top](#toc)\]
$$\label{ccode3d}$$

Next, we'll output the above expressions as Ccode, using the [NRPy+ finite-differencing C code kernel generation infrastructure](Tutorial-Finite_Difference_Derivatives.ipynb). This code will represent spatial derivatives as 10th-order finite differences and output the C code with [SIMD](https://en.wikipedia.org/wiki/SIMD) enabled. ([Common-subexpression elimination](https://en.wikipedia.org/wiki/Common_subexpression_elimination) is enabled by default.)


```python
# Step 7: Set the finite differencing order to 10.
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 10)

# Step 8: Generate C code for scalarwave evolution equations,
#         print output to the screen (standard out, or stdout).
fin.FD_outputC("stdout",
               [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
                lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)], params="enable_SIMD=True")
```

    {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
       */
      /*
       *  Original SymPy expressions:
       *  "[const REAL_SIMD_ARRAY uu_dDD00 = invdx0**2*(-5269*uu/1800 + 5*uu_i0m1_i1_i2/3 - 5*uu_i0m2_i1_i2/21 + 5*uu_i0m3_i1_i2/126 - 5*uu_i0m4_i1_i2/1008 + uu_i0m5_i1_i2/3150 + 5*uu_i0p1_i1_i2/3 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/126 - 5*uu_i0p4_i1_i2/1008 + uu_i0p5_i1_i2/3150),
       *    const REAL_SIMD_ARRAY uu_dDD11 = invdx1**2*(-5269*uu/1800 + 5*uu_i0_i1m1_i2/3 - 5*uu_i0_i1m2_i2/21 + 5*uu_i0_i1m3_i2/126 - 5*uu_i0_i1m4_i2/1008 + uu_i0_i1m5_i2/3150 + 5*uu_i0_i1p1_i2/3 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/126 - 5*uu_i0_i1p4_i2/1008 + uu_i0_i1p5_i2/3150),
       *    const REAL_SIMD_ARRAY uu_dDD22 = invdx2**2*(-5269*uu/1800 + 5*uu_i0_i1_i2m1/3 - 5*uu_i0_i1_i2m2/21 + 5*uu_i0_i1_i2m3/126 - 5*uu_i0_i1_i2m4/1008 + uu_i0_i1_i2m5/3150 + 5*uu_i0_i1_i2p1/3 - 5*uu_i0_i1_i2p2/21 + 5*uu_i0_i1_i2p3/126 - 5*uu_i0_i1_i2p4/1008 + uu_i0_i1_i2p5/3150)]"
       */
      const REAL_SIMD_ARRAY uu_i0_i1_i2m5 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-5)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2m4 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-4)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2m3 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-3)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2m2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-2)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2m1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2-1)]);
      const REAL_SIMD_ARRAY uu_i0_i1m5_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-5,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1m4_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-4,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1m3_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-3,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1m2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-2,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1m1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1-1,i2)]);
      const REAL_SIMD_ARRAY uu_i0m5_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-5,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0m4_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-4,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0m3_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-3,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0m2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-2,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0m1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0-1,i1,i2)]);
      const REAL_SIMD_ARRAY uu = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0p1_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+1,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0p2_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+2,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0p3_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+3,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0p4_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+4,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0p5_i1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0+5,i1,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1p1_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+1,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1p2_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+2,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1p3_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+3,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1p4_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+4,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1p5_i2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1+5,i2)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2p1 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+1)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2p2 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+2)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2p3 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+3)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2p4 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+4)]);
      const REAL_SIMD_ARRAY uu_i0_i1_i2p5 = ReadSIMD(&in_gfs[IDX4S(UUGF, i0,i1,i2+5)]);
      const REAL_SIMD_ARRAY vv = ReadSIMD(&in_gfs[IDX4S(VVGF, i0,i1,i2)]);
      const double tmpFDPart1_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(tmpFDPart1_NegativeOne_);
      const double tmpFDPart1_Rational_1_3150 = 1.0/3150.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_1_3150 = ConstSIMD(tmpFDPart1_Rational_1_3150);
      const double tmpFDPart1_Rational_5269_1800 = 5269.0/1800.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_5269_1800 = ConstSIMD(tmpFDPart1_Rational_5269_1800);
      const double tmpFDPart1_Rational_5_1008 = 5.0/1008.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_5_1008 = ConstSIMD(tmpFDPart1_Rational_5_1008);
      const double tmpFDPart1_Rational_5_126 = 5.0/126.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_5_126 = ConstSIMD(tmpFDPart1_Rational_5_126);
      const double tmpFDPart1_Rational_5_21 = 5.0/21.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_5_21 = ConstSIMD(tmpFDPart1_Rational_5_21);
      const double tmpFDPart1_Rational_5_3 = 5.0/3.0;
      const REAL_SIMD_ARRAY FDPart1_Rational_5_3 = ConstSIMD(tmpFDPart1_Rational_5_3);
      const REAL_SIMD_ARRAY FDPart1_0 = MulSIMD(FDPart1_Rational_5269_1800, uu);
      const REAL_SIMD_ARRAY uu_dDD00 = MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0m3_i1_i2, uu_i0p3_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0m1_i1_i2, uu_i0p1_i1_i2), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0m5_i1_i2, uu_i0p5_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0m4_i1_i2, uu_i0p4_i1_i2), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0m2_i1_i2, uu_i0p2_i1_i2), FDPart1_0))))));
      const REAL_SIMD_ARRAY uu_dDD11 = MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0_i1m3_i2, uu_i0_i1p3_i2), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0_i1m1_i2, uu_i0_i1p1_i2), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0_i1m5_i2, uu_i0_i1p5_i2), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0_i1m4_i2, uu_i0_i1p4_i2), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0_i1m2_i2, uu_i0_i1p2_i2), FDPart1_0))))));
      const REAL_SIMD_ARRAY uu_dDD22 = MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(FDPart1_Rational_5_126, AddSIMD(uu_i0_i1_i2m3, uu_i0_i1_i2p3), FusedMulAddSIMD(FDPart1_Rational_5_3, AddSIMD(uu_i0_i1_i2m1, uu_i0_i1_i2p1), FusedMulSubSIMD(FDPart1_Rational_1_3150, AddSIMD(uu_i0_i1_i2m5, uu_i0_i1_i2p5), FusedMulAddSIMD(FDPart1_Rational_5_1008, AddSIMD(uu_i0_i1_i2m4, uu_i0_i1_i2p4), FusedMulAddSIMD(FDPart1_Rational_5_21, AddSIMD(uu_i0_i1_i2m2, uu_i0_i1_i2p2), FDPart1_0))))));
      /*
       * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
       */
      /*
       *  Original SymPy expressions:
       *  "[const REAL_SIMD_ARRAY __RHS_exp_0 = vv,
       *    const REAL_SIMD_ARRAY __RHS_exp_1 = wavespeed**2*(uu_dDD00 + uu_dDD11 + uu_dDD22)]"
       */
      const REAL_SIMD_ARRAY __RHS_exp_0 = vv;
      const REAL_SIMD_ARRAY __RHS_exp_1 = MulSIMD(MulSIMD(wavespeed, wavespeed), AddSIMD(uu_dDD00, AddSIMD(uu_dDD11, uu_dDD22)));
      WriteSIMD(&rhs_gfs[IDX4S(UUGF, i0, i1, i2)], __RHS_exp_0);
      WriteSIMD(&rhs_gfs[IDX4S(VVGF, i0, i1, i2)], __RHS_exp_1);
    }


<a id='id'></a>

# Step 4: Setting up Initial Data for the Scalar Wave Equation \[Back to [top](#toc)\]
$$\label{id}$$

<a id='planewave'></a>

## Step 4.a: The Monochromatic Plane-Wave Solution \[Back to [top](#toc)\]
$$\label{planewave}$$

The solution to the scalar wave equation for a monochromatic (single-wavelength) wave traveling in the $\hat{k}$ direction is
$$u(\vec{x},t) = f(\hat{k}\cdot\vec{x} - c t),$$
where $\hat{k}$ is a unit vector. We choose $f(\hat{k}\cdot\vec{x} - c t)$ to take the form
$$
f(\hat{k}\cdot\vec{x} - c t) = \sin\left(\hat{k}\cdot\vec{x} - c t\right) + 2,
$$
where we add the $+2$ to ensure that the exact solution never crosses through zero. In places where the exact solution passes through zero, the relative error (i.e., the most common error measure used to check that the numerical solution converges to the exact solution) is undefined. Also, $f(\hat{k}\cdot\vec{x} - c t)$ plus a constant is still a solution to the wave equation.


```python
# Step 1: Set parameters defined in other modules
xx = gri.xx # Sets the Cartesian coordinates xx[0]=x; xx[1]=y; xx[2]=z

# Step 2: Declare free parameters intrinsic to these initial data
time = par.Cparameters("REAL", thismodule, "time",0.0)
kk = par.Cparameters("REAL", thismodule, ["kk0", "kk1", "kk2"],[1.0,1.0,1.0])

# Step 3: Normalize the k vector
kk_norm = sp.sqrt(kk[0]**2 + kk[1]**2 + kk[2]**2)

# Step 4: Compute k.x
dot_product = sp.sympify(0)
for i in range(DIM):
    dot_product += xx[i]*kk[i]
dot_product /= kk_norm

# Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
uu_ID_PlaneWave = sp.sin(dot_product - wavespeed*time)+2
vv_ID_PlaneWave = sp.diff(uu_ID_PlaneWave, time)
```

Next we verify that $f(\hat{k}\cdot\vec{x} - c t)$ satisfies the wave equation, by computing
$$\left(c^2 \nabla^2 - \partial_t^2 \right)\ f\left(\hat{k}\cdot\vec{x} - c t\right),$$
and confirming the result is exactly zero.


```python
sp.simplify(wavespeed**2*(sp.diff(uu_ID_PlaneWave,xx[0],2) +
                          sp.diff(uu_ID_PlaneWave,xx[1],2) +
                          sp.diff(uu_ID_PlaneWave,xx[2],2))
            - sp.diff(uu_ID_PlaneWave,time,2))
```




$\displaystyle 0$



<a id='sphericalgaussian'></a>

## Step 4.b: The Spherical Gaussian Solution \[Back to [top](#toc)\]
$$\label{sphericalgaussian}$$

Here we will implement the spherical Gaussian solution, which consists of ingoing and outgoing wavefronts:
\begin{align}
u(r,t) &= u_{\rm out}(r,t) + u_{\rm in}(r,t) + 1,\ \ \text{where}\\
u_{\rm out}(r,t) &=\frac{r-ct}{r} \exp\left[\frac{-(r-ct)^2}{2 \sigma^2}\right] \\
u_{\rm in}(r,t) &=\frac{r+ct}{r} \exp\left[\frac{-(r+ct)^2}{2 \sigma^2}\right] \\
\end{align}
where $c$ is the wavespeed, and $\sigma$ is the width of the Gaussian (i.e., the "standard deviation").


```python
# Step 1: Set parameters defined in other modules
xx = gri.xx # Sets the Cartesian coordinates xx[0]=x; xx[1]=y; xx[2]=z

# Step 2: Declare free parameters intrinsic to these initial data
time  = par.Cparameters("REAL", thismodule, "time", 0.0)
sigma = par.Cparameters("REAL", thismodule, "sigma", 2.0)

# Step 4: Compute r
r = sp.sympify(0)
for i in range(DIM):
    r += xx[i]**2
r = sp.sqrt(r)

# Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
uu_ID_SphericalGaussianOUT = +(r - wavespeed*time)/r * sp.exp( -(r - wavespeed*time)**2 / (2*sigma**2) )
uu_ID_SphericalGaussianIN  = +(r + wavespeed*time)/r * sp.exp( -(r + wavespeed*time)**2 / (2*sigma**2) )
uu_ID_SphericalGaussian = uu_ID_SphericalGaussianOUT + uu_ID_SphericalGaussianIN + sp.sympify(2)
vv_ID_SphericalGaussian = sp.diff(uu_ID_SphericalGaussian, time)
```

Since the wave equation is linear, both the leftgoing and rightgoing waves must satisfy the wave equation, which implies that their sum also satisfies the wave equation. 

Next we verify that $u(r,t)$ satisfies the wave equation, by computing
$$\left(c^2 \nabla^2 - \partial_t^2 \right)\left\{u_{\rm R}(r,t)\right\},$$

and

$$\left(c^2 \nabla^2 - \partial_t^2 \right)\left\{u_{\rm L}(r,t)\right\},$$

are separately zero. We do this because SymPy has difficulty simplifying the combined expression.


```python
print(sp.simplify(wavespeed**2*(sp.diff(uu_ID_SphericalGaussianOUT,xx[0],2) +
                                sp.diff(uu_ID_SphericalGaussianOUT,xx[1],2) +
                                sp.diff(uu_ID_SphericalGaussianOUT,xx[2],2))
                - sp.diff(uu_ID_SphericalGaussianOUT,time,2)) )

print(sp.simplify(wavespeed**2*(sp.diff(uu_ID_SphericalGaussianIN,xx[0],2) +
                                sp.diff(uu_ID_SphericalGaussianIN,xx[1],2) +
                                sp.diff(uu_ID_SphericalGaussianIN,xx[2],2))
                    - sp.diff(uu_ID_SphericalGaussianIN,time,2)))
```

    0
    0


<a id='code_validation2'></a>

# Step 5: Code Validation against `ScalarWave.InitialData` NRPy+ module  \[Back to [top](#toc)\]
$$\label{code_validation2}$$

As a code validation check, we will verify agreement in the SymPy expressions for plane-wave initial data for the Scalar Wave equation between
1. this tutorial and 
2. the NRPy+ [ScalarWave.InitialData](../edit/ScalarWave/InitialData.py) module.


```python
# We just defined SymPy expressions for uu_ID and vv_ID in
# terms of other SymPy variables. Here, we will use the
# above-defined uu_ID and vv_ID to validate against the
# same expressions in the ScalarWave/InitialData.py
# module, to ensure consistency between this tutorial
# (historically speaking, the tutorial was written first)
# and the PlaneWave ID module itself.
#
# Step 6: Call the InitialData(Type="PlaneWave") function from within the
#         ScalarWave/InitialData.py module,
#         which should do exactly the same as in Steps 1-5 above.
import ScalarWave.InitialData as swid
swid.InitialData(Type="PlaneWave")

# Step 7: Consistency check between the tutorial notebook above
#         and the PlaneWave option from within the
#         ScalarWave/InitialData.py module.
print("Consistency check between ScalarWave tutorial and NRPy+ module: PlaneWave Case")
if sp.simplify(uu_ID_PlaneWave - swid.uu_ID) != 0:
    print("TEST FAILED: uu_ID_PlaneWave - swid.uu_ID = "+str(sp.simplify(uu_ID_PlaneWave - swid.uu_ID))+"\t\t (should be zero)")
    sys.exit(1)
if sp.simplify(vv_ID_PlaneWave - swid.vv_ID) != 0:
    print("TEST FAILED: vv_ID_PlaneWave - swid.vv_ID = "+str(sp.simplify(vv_ID_PlaneWave - swid.vv_ID))+"\t\t (should be zero)")
    sys.exit(1)
print("TESTS PASSED!")


# Step 8: Consistency check between the tutorial notebook above
#         and the SphericalGaussian option from within the
#         ScalarWave/InitialData.py module.
import sys
swid.InitialData(Type="SphericalGaussian")
print("Consistency check between ScalarWave tutorial and NRPy+ module: SphericalGaussian Case")
if sp.simplify(uu_ID_SphericalGaussian - swid.uu_ID) != 0:
    print("TEST FAILED: uu_ID_SphericalGaussian - swid.uu_ID = "+str(sp.simplify(uu_ID_SphericalGaussian - swid.uu_ID))+"\t\t (should be zero)")
    sys.exit(1)
if sp.simplify(vv_ID_SphericalGaussian - swid.vv_ID) != 0:
    print("TEST FAILED: vv_ID_SphericalGaussian - swid.vv_ID = "+str(sp.simplify(vv_ID_SphericalGaussian - swid.vv_ID))+"\t\t (should be zero)")
    sys.exit(1)
print("TESTS PASSED!")
```

    Consistency check between ScalarWave tutorial and NRPy+ module: PlaneWave Case
    TESTS PASSED!
    Consistency check between ScalarWave tutorial and NRPy+ module: SphericalGaussian Case
    TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ScalarWave.pdf](Tutorial-ScalarWave.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ScalarWave")
```

    Created Tutorial-ScalarWave.tex, and compiled LaTeX file to PDF file
        Tutorial-ScalarWave.pdf

