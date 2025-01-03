<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Generating C code for the right-hand-side of the scalar wave equation, in ***curvilinear*** coordinates, using a reference metric formalism

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook outlines a method for solving the right-hand side of the scalar wave equation in curvilinear coordinates using a reference metric formalism, with an emphasis on spherical coordinates. The tutorial utilizes NRPy+ for deriving the necessary contracted Christoffel symbols and handling the scalar wave equation in these curvilinear coordinates.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). In addition, all expressions have been validated against a trusted code (the [original SENR/NRPy+ code](https://bitbucket.org/zach_etienne/nrpy)).

### NRPy+ Source Code for this module: [ScalarWaveCurvilinear/ScalarWaveCurvilinear_RHSs.py](../edit/ScalarWaveCurvilinear/ScalarWaveCurvilinear_RHSs.py)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

0. [Preliminaries](#prelim): Reference Metrics and Picking Best Coordinate System to Solve the PDE
1. [Example](#example): The scalar wave equation in spherical coordinates
1. [Step 1](#contracted_christoffel): Contracted Christoffel symbols $\hat{\Gamma}^i = \hat{g}^{ij}\hat{\Gamma}^k_{ij}$ in spherical coordinates, using NRPy+
1. [Step 2](#rhs_scalarwave_spherical): The right-hand side of the scalar wave equation in spherical coordinates, using NRPy+
1. [Step 3](#code_validation): Code Validation against `ScalarWave.ScalarWaveCurvilinear_RHSs`  NRPy+ Module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='prelim'></a>

# Preliminaries: Reference Metrics and Picking Best Coordinate System to Solve the PDE \[Back to [top](#toc)\]
$$\label{prelim}$$

Recall from [NRPy+ tutorial notebook on the Cartesian scalar wave equation](Tutorial-ScalarWave.ipynb) that the scalar wave equation in 3D Cartesian coordinates is given by

$$\partial_t^2 u = c^2 \nabla^2 u \text{,}$$
where $u$ (the amplitude of the wave) is a function of time and Cartesian coordinates in space: $u = u(t,x,y,z)$ (spatial dimension as-yet unspecified), and subject to some initial condition
$$u(0,x,y,z) = f(x,y,z),$$

with suitable (sometimes approximate) spatial boundary conditions.

To simplify this equation, let's first choose units such that $c=1$. Alternative wave speeds can be constructed
by simply rescaling the time coordinate, with the net effect being that the time $t$ is replaced with time in dimensions of space; i.e., $t\to c t$.

$$\partial_t^2 u = \nabla^2 u$$

As we learned in the [NRPy+ tutorial notebook on reference metrics](Tutorial-Reference_Metric.ipynb), reference metrics are a means to pick the best coordinate system for the PDE we wish to solve. However, to take advantage of reference metrics requires first that we generalize the PDE. In the case of the scalar wave equation, this involves first rewriting in [Einstein notation](https://en.wikipedia.org/wiki/Einstein_notation) (with implied summation over repeated indices) via

$$(-\partial_t^2 + \nabla^2) u = \eta^{\mu\nu} u_{,\ \mu\nu} = 0,$$

where $u_{,\mu\nu} = \partial_\mu \partial_\nu u$, and $\eta^{\mu\nu}$ is the contravariant flat-space metric tensor with components $\text{diag}(-1,1,1,1)$.

Next, we apply the "comma-goes-to-semicolon rule" and replace $\eta^{\mu\nu}$ with $\hat{g}^{\mu\nu}$ to generalize the scalar wave equation to an arbitrary reference metric $\hat{g}^{\mu\nu}$:

$$\hat{g}^{\mu\nu} u_{;\ \mu\nu} = \hat{g}^{\mu\nu} \hat{\nabla}_{\mu} \hat{\nabla}_{\nu} u = 0,$$

where $\hat{\nabla}_{\mu}$ denotes the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) with respect to the reference metric basis vectors $\hat{x}^{\mu}$, and $\hat{g}^{\mu \nu} \hat{\nabla}_{\mu} \hat{\nabla}_{\nu} u$ is the covariant
[D'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) of $u$.

For example, suppose we wish to model a short-wavelength wave that is nearly spherical. In this case, if we were to solve the wave equation PDE in Cartesian coordinates, we would in principle need high resolution in all three cardinal directions. If instead, we chose spherical coordinates centered at the center of the wave, we might need high resolution only in the radial direction, with only a few points required in the angular directions. Thus choosing spherical coordinates would be far more computationally efficient than modeling the wave in Cartesian coordinates.

Let's now expand the covariant scalar wave equation in arbitrary coordinates. Since the covariant derivative of a scalar is equivalent to its partial derivative, we have
\begin{align}
0 &= \hat{g}^{\mu \nu} \hat{\nabla}_{\mu} \hat{\nabla}_{\nu} u \\
&= \hat{g}^{\mu \nu} \hat{\nabla}_{\mu} \partial_{\nu} u.
\end{align}

$\partial_{\nu} u$ transforms as a one-form under covariant differentiation, so we have
$$\hat{\nabla}_{\mu} \partial_{\nu} u = \partial_{\mu} \partial_{\nu} u - \hat{\Gamma}^\tau_{\mu\nu} \partial_\tau u,$$
where 

$$\hat{\Gamma}^\tau_{\mu\nu} = \frac{1}{2} \hat{g}^{\tau\alpha} \left(\partial_\nu \hat{g}_{\alpha\mu} + \partial_\mu \hat{g}_{\alpha\nu} - \partial_\alpha \hat{g}_{\mu\nu} \right)$$
are the [Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols) associated with the reference metric $\hat{g}_{\mu\nu}$.

Then the scalar wave equation is written as
$$0 = \hat{g}^{\mu \nu} \left( \partial_{\mu} \partial_{\nu} u - \hat{\Gamma}^\tau_{\mu\nu} \partial_\tau u\right).$$

Define the contracted Christoffel symbols as
$$\hat{\Gamma}^\tau = \hat{g}^{\mu\nu} \hat{\Gamma}^\tau_{\mu\nu}.$$

Then the scalar wave equation is given by
$$0 = \hat{g}^{\mu \nu} \partial_{\mu} \partial_{\nu} u - \hat{\Gamma}^\tau \partial_\tau u.$$

The reference metrics we adopt satisfy
$$\hat{g}^{t \nu} = -\delta^{t \nu},$$
where $\delta^{t \nu}$ is the [Kronecker delta](https://en.wikipedia.org/wiki/Kronecker_delta). Therefore the scalar wave equation in curvilinear coordinates can be written
\begin{align}
0 &= \hat{g}^{\mu \nu} \partial_{\mu} \partial_{\nu} u - \hat{\Gamma}^\tau \partial_\tau u \\
&= -\partial_t^2 u + \hat{g}^{i j} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u \\
\implies \partial_t^2 u &= \hat{g}^{i j} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u,
\end{align}
where repeated Latin indices denote implied summation over *spatial* components only. This module implements the bottom equation for arbitrary reference metrics satisfying $\hat{g}^{t \nu} = -\delta^{t \nu}$. To gain an appreciation for what NRPy+ accomplishes automatically, let's first work out the scalar wave equation in spherical coordinates by hand.

<a id='example'></a>

# Example: The scalar wave equation in spherical coordinates \[Back to [top](#toc)\]
$$\label{example}$$

For example, the spherical reference metric is written

$$\hat{g}_{\mu\nu} = \begin{pmatrix}
-1 & 0 & 0 & 0 \\
 0 & 1 & 0 & 0 \\
 0 & 0 & r^2 & 0 \\
 0 & 0 & 0 & r^2 \sin^2 \theta \\
\end{pmatrix}.
$$

Since the inverse of a diagonal matrix is simply the inverse of the diagonal elements, we can write 
$$\hat{g}^{\mu\nu} = \begin{pmatrix}
-1 & 0 & 0 & 0 \\
 0 & 1 & 0 & 0 \\
 0 & 0 & \frac{1}{r^2} & 0 \\
 0 & 0 & 0 & \frac{1}{r^2 \sin^2 \theta} \\
\end{pmatrix}.$$

The scalar wave equation in these coordinates can thus be written
\begin{align}
0 &= \hat{g}^{\mu \nu} \partial_{\mu} \partial_{\nu} u - \hat{\Gamma}^\tau \partial_\tau u \\
&= \hat{g}^{tt} \partial_t^2 u + \hat{g}^{rr} \partial_r^2 u + \hat{g}^{\theta\theta} \partial_\theta^2 u  + \hat{g}^{\phi\phi} \partial_\phi^2 u - \hat{\Gamma}^\tau \partial_\tau u \\
&= -\partial_t^2 u + \partial_r^2 u + \frac{1}{r^2} \partial_\theta^2
u + \frac{1}{r^2 \sin^2 \theta} \partial_\phi^2 u - \hat{\Gamma}^\tau \partial_\tau u\\
\implies \partial_t^2 u &= \partial_r^2 u + \frac{1}{r^2} \partial_\theta^2
u + \frac{1}{r^2 \sin^2 \theta} \partial_\phi^2 u - \hat{\Gamma}^\tau \partial_\tau u.
\end{align}

The contracted Christoffel symbols 
$\hat{\Gamma}^\tau$ can then be computed directly from the metric $\hat{g}_{\mu\nu}$.

It can be shown (exercise to the reader) that the only nonzero
components of $\hat{\Gamma}^\tau$ in static spherical polar coordinates are
given by
\begin{align}
\hat{\Gamma}^r &= -\frac{2}{r} \\
\hat{\Gamma}^\theta &= -\frac{\cos\theta}{r^2 \sin\theta}.
\end{align}

Thus we have found the Laplacian in spherical coordinates is simply:

\begin{align}
\nabla^2 u &= 
\partial_r^2 u + \frac{1}{r^2} \partial_\theta^2 u + \frac{1}{r^2 \sin^2 \theta} \partial_\phi^2 u - \hat{\Gamma}^\tau \partial_\tau u\\
&= \partial_r^2 u + \frac{1}{r^2} \partial_\theta^2 u + \frac{1}{r^2 \sin^2 \theta} \partial_\phi^2 u +  \frac{2}{r} \partial_r u + \frac{\cos\theta}{r^2 \sin\theta} \partial_\theta u
\end{align}
(cf. http://mathworld.wolfram.com/SphericalCoordinates.html; though note that they defined the angle $\phi$ as $\theta$ and $\theta$ as $\phi$.)

<a id='contracted_christoffel'></a>

# Step 1: Contracted Christoffel symbols $\hat{\Gamma}^i = \hat{g}^{ij}\hat{\Gamma}^k_{ij}$ in spherical coordinates, using NRPy+ \[Back to [top](#toc)\]
$$\label{contracted_christoffel}$$

Let's next use NRPy+ to derive the contracted Christoffel symbols
$$\hat{g}^{ij} \hat{\Gamma}^k_{ij}$$
in spherical coordinates, where $i\in\{1,2,3\}$ and $j\in\{1,2,3\}$ are spatial indices.

As discussed in the [NRPy+ tutorial notebook on reference metrics](Tutorial-Reference_Metric.ipynb), several reference-metric-related quantities in spherical coordinates are computed in NRPy+ (provided the parameter **`reference_metric::CoordSystem`** is set to **`"Spherical"`**), including the inverse spatial spherical reference metric $\hat{g}^{ij}$ and the Christoffel symbols from this reference metric $\hat{\Gamma}^{i}_{jk}$. 


```python
import sympy as sp              # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par  # NRPy+: Parameter interface
import grid as gri              # NRPy+: Functionality for handling numerical grids
import indexedexp as ixp        # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm  # NRPy+: Reference metric support

# reference_metric::CoordSystem can be set to Spherical, SinhSpherical, SinhSphericalv2,
#                           Cylindrical, SinhCylindrical, SinhCylindricalv2, etc.
#                           See reference_metric.py and NRPy+ tutorial notebook on
#                           reference metrics for full list and description of how
#                           to extend.
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
par.set_parval_from_str("grid::DIM",3)

rfm.reference_metric()

contractedGammahatU = ixp.zerorank1()
for k in range(3):
    for i in range(3):
        for j in range(3):
            contractedGammahatU[k] += rfm.ghatUU[i][j] * rfm.GammahatUDD[k][i][j]

for k in range(3):
    print("contracted GammahatU["+str(k)+"]:")
    print(sp.simplify(contractedGammahatU[k]))
    # Sadly pretty_print results in garbage output in the generated PDF at the bottom of this notebook.
#     sp.pretty_print(sp.simplify(contractedGammahatU[k]))
    if k<2:
        print("\n\n")
```

    contracted GammahatU[0]:
    -2/xx0
    
    
    
    contracted GammahatU[1]:
    -1/(xx0**2*tan(xx1))
    
    
    
    contracted GammahatU[2]:
    0


<a id='rhs_scalarwave_spherical'></a>

# Step 2: The right-hand side of the scalar wave equation in spherical coordinates, using NRPy+ \[Back to [top](#toc)\]
$$\label{rhs_scalarwave_spherical}$$

Following our [implementation of the scalar wave equation in Cartesian coordinates](Tutorial-ScalarWave.ipynb), we will introduce a new variable $v=\partial_t u$ that will enable us to split the second time derivative into two first-order time derivatives.

\begin{align}
\partial_t u &= v \\
\partial_t v &= \hat{g}^{ij} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u
\end{align}

Adding back the sound speed $c$, we have a choice of a single factor of $c$ multiplying both right-hand sides, or a factor of $c^2$ multiplying the second equation only. We'll choose the latter.

\begin{align}
\partial_t u &= v \\
\partial_t v &= c^2 \left(\hat{g}^{ij} \partial_{i} \partial_{j} u - \hat{\Gamma}^i \partial_i u\right)
\end{align}

Now let's generate the C code for the finite-difference representations of the right-hand sides of the above "time evolution" equations for $u$ and $v$. Since the right-hand side of $\partial_t v$ contains implied sums over $i$ and $j$ in the first term, and an implied sum over $k$ in the second term, we'll find it useful to split the right-hand side into two parts,

\begin{equation}
\partial_t v = c^2 \left(
{\underbrace {\textstyle \hat{g}^{ij} \partial_{i} \partial_{j} u}_{\text{Part 1}}} 
{\underbrace {\textstyle -\hat{\Gamma}^i \partial_i u}_{\text{Part 2}}}\right),
\end{equation}

and perform the implied sums in two pieces.


```python
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import reference_metric as rfm   # NRPy+: Reference metric support
from outputC import lhrh         # NRPy+: Core C code output module
```


```python
# The name of this module ("scalarwave") is given by __name__:
thismodule = __name__

# Step 0: Read the spatial dimension parameter as DIM.
DIM = par.parval_from_str("grid::DIM")

# Step 1: Set the finite differencing order to 4.
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",4)

# Step 2a: Reset the gridfunctions list; below we define the
#          full complement of gridfunctions needed by this
#          tutorial. This line of code enables us to re-run this
#          tutorial without resetting the running Python kernel.
gri.glb_gridfcs_list = []
# Step 2b: Register gridfunctions that are needed as input
#          to the scalar wave RHS expressions.
uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])

# Step 3a: Declare the rank-1 indexed expression \partial_{i} u,
#          Derivative variables like these must have an underscore
#          in them, so the finite difference module can parse the
#          variable name properly.
uu_dD = ixp.declarerank1("uu_dD")

# Step 3b: Declare the rank-2 indexed expression \partial_{ij} u,
#          which is symmetric about interchange of indices i and j
#          Derivative variables like these must have an underscore
#          in them, so the finite difference module can parse the
#          variable name properly.
uu_dDD = ixp.declarerank2("uu_dDD","sym01")

# Step 4: Define the C parameter wavespeed. The `wavespeed`
#         variable is a proper SymPy variable, so it can be
#         used in below expressions. In the C code, it acts
#         just like a usual parameter, whose value is
#         specified in the parameter file.
wavespeed = par.Cparameters("REAL",thismodule,"wavespeed", 1.0)

# Step 5: Define right-hand sides for the evolution.
uu_rhs = vv
# Step 5b: The right-hand side of the \partial_t v equation
#          is given by:
#          \hat{g}^{ij} \partial_i \partial_j u - \hat{\Gamma}^i \partial_i u.
#          ^^^^^^^^^^^^ PART 1 ^^^^^^^^^^^^^^^^ ^^^^^^^^^^ PART 2 ^^^^^^^^^^^
vv_rhs = 0
for i in range(DIM):
    # PART 2:
    vv_rhs -= contractedGammahatU[i]*uu_dD[i]
    for j in range(DIM):
        # PART 1:
        vv_rhs += rfm.ghatUU[i][j]*uu_dDD[i][j]

vv_rhs *= wavespeed*wavespeed

# Step 6: Generate C code for scalarwave evolution equations,
#         print output to the screen (standard out, or stdout).
fin.FD_outputC("stdout",
               [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
                lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])
```

    {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
       */
      /*
       *  Original SymPy expressions:
       *  "[const double uu_dD0 = invdx0*(-2*uu_i0m1_i1_i2/3 + uu_i0m2_i1_i2/12 + 2*uu_i0p1_i1_i2/3 - uu_i0p2_i1_i2/12),
       *    const double uu_dD1 = invdx1*(-2*uu_i0_i1m1_i2/3 + uu_i0_i1m2_i2/12 + 2*uu_i0_i1p1_i2/3 - uu_i0_i1p2_i2/12),
       *    const double uu_dDD00 = invdx0**2*(-5*uu/2 + 4*uu_i0m1_i1_i2/3 - uu_i0m2_i1_i2/12 + 4*uu_i0p1_i1_i2/3 - uu_i0p2_i1_i2/12),
       *    const double uu_dDD11 = invdx1**2*(-5*uu/2 + 4*uu_i0_i1m1_i2/3 - uu_i0_i1m2_i2/12 + 4*uu_i0_i1p1_i2/3 - uu_i0_i1p2_i2/12),
       *    const double uu_dDD22 = invdx2**2*(-5*uu/2 + 4*uu_i0_i1_i2m1/3 - uu_i0_i1_i2m2/12 + 4*uu_i0_i1_i2p1/3 - uu_i0_i1_i2p2/12)]"
       */
      const double uu_i0_i1_i2m2 = in_gfs[IDX4S(UUGF, i0,i1,i2-2)];
      const double uu_i0_i1_i2m1 = in_gfs[IDX4S(UUGF, i0,i1,i2-1)];
      const double uu_i0_i1m2_i2 = in_gfs[IDX4S(UUGF, i0,i1-2,i2)];
      const double uu_i0_i1m1_i2 = in_gfs[IDX4S(UUGF, i0,i1-1,i2)];
      const double uu_i0m2_i1_i2 = in_gfs[IDX4S(UUGF, i0-2,i1,i2)];
      const double uu_i0m1_i1_i2 = in_gfs[IDX4S(UUGF, i0-1,i1,i2)];
      const double uu = in_gfs[IDX4S(UUGF, i0,i1,i2)];
      const double uu_i0p1_i1_i2 = in_gfs[IDX4S(UUGF, i0+1,i1,i2)];
      const double uu_i0p2_i1_i2 = in_gfs[IDX4S(UUGF, i0+2,i1,i2)];
      const double uu_i0_i1p1_i2 = in_gfs[IDX4S(UUGF, i0,i1+1,i2)];
      const double uu_i0_i1p2_i2 = in_gfs[IDX4S(UUGF, i0,i1+2,i2)];
      const double uu_i0_i1_i2p1 = in_gfs[IDX4S(UUGF, i0,i1,i2+1)];
      const double uu_i0_i1_i2p2 = in_gfs[IDX4S(UUGF, i0,i1,i2+2)];
      const double vv = in_gfs[IDX4S(VVGF, i0,i1,i2)];
      const double FDPart1_Rational_2_3 = 2.0/3.0;
      const double FDPart1_Rational_1_12 = 1.0/12.0;
      const double FDPart1_Rational_5_2 = 5.0/2.0;
      const double FDPart1_Rational_4_3 = 4.0/3.0;
      const double FDPart1_0 = -FDPart1_Rational_5_2*uu;
      const double uu_dD0 = invdx0*(FDPart1_Rational_1_12*(uu_i0m2_i1_i2 - uu_i0p2_i1_i2) + FDPart1_Rational_2_3*(-uu_i0m1_i1_i2 + uu_i0p1_i1_i2));
      const double uu_dD1 = invdx1*(FDPart1_Rational_1_12*(uu_i0_i1m2_i2 - uu_i0_i1p2_i2) + FDPart1_Rational_2_3*(-uu_i0_i1m1_i2 + uu_i0_i1p1_i2));
      const double uu_dDD00 = ((invdx0)*(invdx0))*(FDPart1_0 + FDPart1_Rational_1_12*(-uu_i0m2_i1_i2 - uu_i0p2_i1_i2) + FDPart1_Rational_4_3*(uu_i0m1_i1_i2 + uu_i0p1_i1_i2));
      const double uu_dDD11 = ((invdx1)*(invdx1))*(FDPart1_0 + FDPart1_Rational_1_12*(-uu_i0_i1m2_i2 - uu_i0_i1p2_i2) + FDPart1_Rational_4_3*(uu_i0_i1m1_i2 + uu_i0_i1p1_i2));
      const double uu_dDD22 = ((invdx2)*(invdx2))*(FDPart1_0 + FDPart1_Rational_1_12*(-uu_i0_i1_i2m2 - uu_i0_i1_i2p2) + FDPart1_Rational_4_3*(uu_i0_i1_i2m1 + uu_i0_i1_i2p1));
      /*
       * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
       */
      /*
       *  Original SymPy expressions:
       *  "[rhs_gfs[IDX4S(UUGF, i0, i1, i2)] = vv,
       *    rhs_gfs[IDX4S(VVGF, i0, i1, i2)] = wavespeed**2*(2*uu_dD0/xx0 + uu_dD1*sin(2*xx1)/(2*xx0**2*sin(xx1)**2) + uu_dDD00 + uu_dDD11/xx0**2 + uu_dDD22/(xx0**2*sin(xx1)**2))]"
       */
      const double FDPart3_0 = (1.0/((xx0)*(xx0)));
      const double FDPart3_1 = FDPart3_0/((sin(xx1))*(sin(xx1)));
      rhs_gfs[IDX4S(UUGF, i0, i1, i2)] = vv;
      rhs_gfs[IDX4S(VVGF, i0, i1, i2)] = ((wavespeed)*(wavespeed))*(FDPart3_0*uu_dDD11 + (1.0/2.0)*FDPart3_1*uu_dD1*sin(2*xx1) + FDPart3_1*uu_dDD22 + 2*uu_dD0/xx0 + uu_dDD00);
    }


<a id='code_validation'></a>

# Step 3: Code Validation against `ScalarWave.ScalarWaveCurvilinear_RHSs`  NRPy+ Module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the Curvilinear Scalar Wave equation (i.e., uu_rhs and vv_rhs) between

1. this tutorial and 
2. the NRPy+ [ScalarWave.ScalarWaveCurvilinear_RHSs](../edit/ScalarWaveCurvilinear/ScalarWaveCurvilinear_RHSs.py) module.

By default, we analyze the RHSs in Spherical coordinates, though other coordinate systems may be chosen.


```python
# Step 7: We already have SymPy expressions for uu_rhs and vv_rhs in
#         terms of other SymPy variables. Even if we reset the list
#         of NRPy+this gridfunctions, these *SymPy* expressions for
#         uu_rhs and vv_rhs *will remain unaffected*.
#
#         Here, we will use the above-defined uu_rhs and vv_rhs to
#         validate against the same expressions in the
#         ScalarWaveCurvilinear/ScalarWaveCurvilinear module,
#         to ensure consistency between the tutorial and the
#         module itself.
#
# Reset the list of gridfunctions, as registering a gridfunction
#   twice will spawn an error.
gri.glb_gridfcs_list = []

# Step 8: Call the ScalarWaveCurvilinear_RHSs() function from within the
#         ScalarWaveCurvilinear/ScalarWaveCurvilinear_RHSs.py module,
#         which should do exactly the same as in Steps 1-6 above.
import ScalarWave.ScalarWaveCurvilinear_RHSs as swcrhs
swcrhs.ScalarWaveCurvilinear_RHSs()

# Step 9: Consistency check between the tutorial notebook above
#         and the ScalarWaveCurvilinear_RHSs() function from within the
#         ScalarWaveCurvilinear/ScalarWaveCurvilinear_RHSs.py module.
print("Consistency check between ScalarWaveCurvilinear tutorial and NRPy+ module:")
if sp.simplify(uu_rhs - swcrhs.uu_rhs) != 0:
    print("TEST FAILED: uu_ID_SphericalGaussian - swid.uu_ID = "+str(sp.simplify(uu_rhs - swcrhs.uu_rhs))+"\t\t (should be zero)")
    sys.exit(1)
if sp.simplify(vv_rhs - swcrhs.vv_rhs) != 0:
    print("TEST FAILED: vv_ID_SphericalGaussian - swid.vv_ID = "+str(sp.simplify(vv_rhs - swcrhs.vv_rhs))+"\t\t (should be zero)")
    sys.exit(1)
print("TESTS PASSED!")
```

    Consistency check between ScalarWaveCurvilinear tutorial and NRPy+ module:
    TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ScalarWaveCurvilinear.pdf](Tutorial-ScalarWaveCurvilinear.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ScalarWaveCurvilinear")
```

    Created Tutorial-ScalarWaveCurvilinear.tex, and compiled LaTeX file to PDF
        file Tutorial-ScalarWaveCurvilinear.pdf

