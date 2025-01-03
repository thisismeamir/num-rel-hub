<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Start-to-Finish Example: Applying Boundary Conditions in Curvilinear Coordinates in Three Dimensions

## Authors: Zach Etienne & Terrence Pierre Jacques

## This paper validates boundary condition algorithms for curvilinear coordinates like Spherical and Cylindrical, as mentioned in the [SENR/NRPy+ paper](https://arxiv.org/abs/1712.07658). Using an Eigen-Coordinate system, it differentiates between inner and outer boundaries, applying parity conditions to arbitrary-rank tensors. The approach is validated by successfully navigating coordinate singularities through a combination of cell-centered grids, tensor rescaling, and a stable Method of Lines time stepping algorithm.

<a id='intro'></a>

## Introduction: Applying boundary conditions in curvilinear coordinates
$$\label{intro}$$

We face four challenges when solving PDEs in curvilinear coordinates:

1. [Challenge 1](#innerouterbcs): Unlike ordinary Cartesian coordinate grids, not all boundary points on uniform curvilinear coordinate grids are outer boundary points.
1. [Challenge 2](#coordinversion): Figuring out the locations to which boundary points map requires a coordinate inversion, but some coordinate systems are not easily invertible.
1. [Challenge 3](#parity): Tensors and vectors in curvilinear coordinates can change directions across "inner" boundaries, changing sign as a result. 
1. [Challenge 4](#coordsingularities): Coordinate singularities can appear, causing a stiffening or divergence of terms in PDEs.

<a id='innerouterbcs'></a>

### Challenge 1: Inner versus outer boundary points
$$\label{innerouterbcs}$$

Unlike Cartesian coordinates, the boundaries of our grid in generic curvilinear coordinates are not all outer boundaries. 

Consider first a computational grid in Cartesian coordinates, with $(x_0,x_1,x_2)=(x,y,z)$, that is *uniform* (i.e., with $\Delta x$, $\Delta y$, and $\Delta z$ all set to constant values). This grid may extend over arbitrary coordinate ranges $x_i \in [x_{i, \rm min},x_{i, \rm max}]$.

By contrast, consider now a uniform grid in spherical coordinates $(x_0,x_1,x_2)=(r,\theta,\phi)$ with constant spacing $\Delta r$, $\Delta \theta$, and $\Delta \phi$ between grid points in $r$, $\theta$, and $\phi$, respectively. Further, let's assume that these grids span all possible values of $\theta$ and $\phi$, with $r=0$ included in the domain. Then our numerical domain must satisfy the following relations

+ $x_0 = r \in [0,{\rm RMAX}]$,
+ $x_1 = \theta \in [0,\pi]$, and
+ $x_2 = \phi \in [-\pi,\pi]$. (Notice how we do not choose $x_2= \phi \in [0,2\pi]$ so that our conversion from Cartesian to spherical coordinates is compatible with the output range from the ${\rm atan2}(y,x)$ function: $\phi={\rm atan2}(y,x)\in[-\pi,\pi]$.)

Notice that, unlike Cartesian coordinates, the boundaries of this numerical grid in spherical coordinates are not all outer boundaries. For example, data stored at $\phi=-\pi$ will be identical to data at $\phi=\pi$, *regardless of $r$ or $\theta$*. I.e., $\phi$ satisfies *periodic* boundary conditions only. Further, $\theta=0$ presents a more complicated boundary condition, in which points with negative $\theta$ map to points with $|\theta|$ but at an angle of $\phi\to \phi+\pi$. Finally, negative $r$ points will map to positive $r$ points on the other side of the origin. We call these boundaries *inner* boundaries, as they generally map to other points in the interior (as opposed to the outer boundaries) of the grid. 

As we generally cannot apply an outer boundary condition to the inner boundaries, these boundaries will need to be treated differently.

On our numerical grids, this poses some difficulty, as finite difference derivatives we compute within the numerical domain require that the grid be extended beyond the domain boundaries. In spherical coordinates, this means that we need additional grid points at, e.g.,  $r<0$, $\theta<0$, and $\phi>\pi$, just to name a few. Whether they be on outer or inner boundaries, we call grid points in the extended region *ghost zones*.

Numerical grids of $N$th order accuracy generally possess $N/2$ ghost zone points in the boundary regions (i.e., $x_i < x_{i,\rm min}$ and $x_i > x_{i, \rm max}$). While in Cartesian coordinates, these ghost zone points map to regions outside the grid domain $x_i \in [x_{i, \rm min},x_{i, \rm max}]$, in spherical coordinates, *most* ghost zone points map to regions *inside* the grid domain. For example, for some $\tilde{r}\in [0,{\rm RMAX}]$ and $\tilde{\theta}\in[0,\pi]$, the ghost zone point $(\tilde{r},\tilde{\theta},2\pi+\Delta \phi/2)$ would map to the interior point $(\tilde{r},\tilde{\theta},\Delta \phi/2)$ because the $\phi$ coordinate is periodic. Thus when given a ghost zone point in some arbitrary curvilinear coordinate system, we are faced with the problem of addressing the following two questions:
1. Does a given ghost point map to an interior point, or is it an outer boundary point (i.e., a point exterior to the domain)?
1. If the ghost zone point maps to an interior point, to which interior point does it map?

<a id='coordinversion'></a>

### Challenge 2: Inverting coordinates
$$\label{coordinversion}$$

Coordinate systems within NRPy+ are generally Spherical-like, Cylindrical-like, SymTP-like (where SymTP is a prolate-spheroidal coordinate system), or Cartesian-like. For example, SinhSphericalv2 coordinates are exactly the same as Spherical coordinates, except we choose an odd function for the radial coordinate $r$ as a function of $x_0$:

$$
r(x_0) = {\rm AMPL} \left[ {\rm const\_dr} x_0 + \sinh\left(\frac{x_0}{\rm SINHW}\right) / \sinh\left(\frac{1}{\rm SINHW}\right) \right].
$$

While this coordinate choice exhibits nice properties for certain cases, the function $x_0(r)$ is not a closed-form expression. Thus finding the mapping of ghost zone points in the radial direction would require a root finder. 

*Is there an easier way of dealing with this problem than with a root finder?*

<a id='parity'></a>

### Challenge 3: Parity: changes of direction in vectors and tensors across inner boundaries
$$\label{parity}$$

When applying inner boundary conditions to vectors and tensors, we must consider how the direction or *parity* of vector and tensor components change across the inner boundary.

Suppose we have a vector $v^\rho$ defined at ghost zone $(-\rho,\phi,z)$ ($\rho>0$) in cylindrical coordinates. This will map to an interior point at $(\rho,\phi+\pi,z)$. At this point, the direction of the $\hat{\rho}$ unit vector flips sign. Thus we cannot simply set the value of $v^\rho$ to the value it possesses at interior point $(\rho,\phi+\pi,z)$; that would result in a sign error. Instead, we have
\begin{align}
v^\rho(-\rho,\phi,z)&=-v^\rho(\rho,\phi+\pi,z) \\
&= \mathbf{e}^\rho\left(-\rho,\phi,z\right) \cdot \mathbf{e}^\rho\left(\rho,\phi+\pi,z\right)v^\rho(\rho,\phi+\pi,z),
\end{align}
where $\mathbf{e}^\rho\left(\rho,\phi,z\right)$ is the $\rho$ unit vector evaluated at point $(\rho,\phi,z)$, and $\mathbf{e}^\rho\left(-\rho,\phi,z\right) \cdot \mathbf{e}^\rho\left(\rho,\phi+\pi,z\right)$ is the dot product of the two unit vectors, which must evaluate to $\pm 1$ (i.e., the **parity**). Contrast this with scalars, which do not possess a sense of direction/parity.

<a id='coordsingularities'></a>

### Challenge 4: Coordinate singularities
$$\label{coordsingularities}$$

Most non-Cartesian, orthogonal coordinate systems (like spherical coordinates) possess *coordinate singularities*. 

For example, coordinate singularities in spherical coordinates lie along $\theta=0$ and $\theta=\pi$; these are points where the coordinate system focuses to a single point. For example, the coordinate singularity at the North Pole is the reason why all directions are south there. Critically, these singularities manifest as points where the reference metric or its inverse crosses through zero or diverges to $\infty$. As we derived in a [previous module](Tutorial-ScalarWaveCurvilinear.ipynb), the Laplacian in spherical polar coordinates takes the form
$$
\nabla^2 u = \partial_r^2 u + \frac{1}{r^2} \partial_\theta^2 u + \frac{1}{r^2 \sin^2 \theta} \partial_\phi^2 u +  \frac{2}{r} \partial_r u + \frac{\cos\theta}{r^2 \sin\theta} \partial_\theta u,
$$
which diverges at $r=0$ and $\sin\theta=0-$precisesly at the $\theta=0$ and $\theta=\pi$ coordinate singularity. 

To avoid this divergence, we simply choose that our numerical grids be **cell-centered**. 

I.e., given the desired bounds of the grid interior to be 

\begin{align}
x_0 &\in [x_{0,\ \rm min},x_{0,\ \rm max}]\\
x_1 &\in [x_{1,\ \rm min},x_{1,\ \rm max}]\\
x_2 &\in [x_{2,\ \rm min},x_{2,\ \rm max}],
\end{align}

${\rm NGHOSTS}$ to be the number of ghost zones (assumed the same in all directions), and $\{N_0,N_1,N_2\}$ to be the desired number of points in the grid interior in the $\{x_0,x_1,x_2\}$ directions, respectively, then the numerical grid spacing in each respective direction will be given by

\begin{align}
dx_0 &= \frac{x_{0,\ \rm max} - x_{0,\ \rm min}}{N_0} \\
dx_1 &= \frac{x_{1,\ \rm max} - x_{1,\ \rm min}}{N_1} \\
dx_2 &= \frac{x_{2,\ \rm max} - x_{2,\ \rm min}}{N_2}.
\end{align}

Given the above definitions, the complete set of indices $\{{\rm i0},{\rm i1},{\rm i2}\}$ located at $\{x_{0,{\rm i0}},x_{1,{\rm i1}},x_{2,{\rm i2}}\}$ as follows:

\begin{align}
x_{0,{\rm i0}} &= x_{0,\ \rm min} + \left[({\rm i0}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_0 \\
x_{1,{\rm i1}} &= x_{1,\ \rm min} + \left[({\rm i1}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_1 \\
x_{2,{\rm i2}} &= x_{2,\ \rm min} + \left[({\rm i2}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_2 \\
\end{align}
where

* ${\rm i0}\in [0,N_0+2\cdot{\rm NGHOSTS})$
* ${\rm i1}\in [0,N_1+2\cdot{\rm NGHOSTS})$
* ${\rm i2}\in [0,N_2+2\cdot{\rm NGHOSTS})$,

which guarantees the interior is covered by exactly $\{N_0,N_1,N_2\}$ grid points, the boundaries are covered by ${\rm NGHOSTS}$ ghost zones, and we maintain cell centering.

So for example, if we choose a numerical grid in *spherical* coordinates $\{r,\theta,\phi\}$, with 3 ghost zone points (needed for e.g., 6th-order-accurate centered finite differencing), and we want the grid interior to be sampled with $\{N_r,N_\theta,N_\phi\}$ grid points, then we have

\begin{align}
{\rm NGHOSTS} &= 3 \\
dr &= \frac{r_{\rm max} - 0}{N_r} \\
d\theta &= \frac{\pi - 0}{N_\theta} \\
d\phi &= \frac{\pi - (-\pi)}{N_\phi} \\
r_{{\rm i0}} &= 0 + \left[({\rm i0}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_0 \\
&= \left[({\rm i0}-3) + \frac{1}{2}\right] dx_0 \\
\theta_{{\rm i1}} &= 0 + \left[({\rm i1}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_1 \\
&= \left[({\rm i1}-3) + \frac{1}{2}\right] dx_1 \\
\phi_{{\rm i2}} &= -\pi + \left[({\rm i2}-{\rm NGHOSTS}) + \frac{1}{2}\right] dx_2 \\
&= -\pi + \left[({\rm i2}-3) + \frac{1}{2}\right] dx_2, \\
\end{align}

where again
* ${\rm i0}\in [0,N_r+2\cdot{\rm NGHOSTS})$
* ${\rm i1}\in [0,N_\theta+2\cdot{\rm NGHOSTS})$
* ${\rm i2}\in [0,N_\phi+2\cdot{\rm NGHOSTS})$,

which guarantees the interior is covered by exactly $\{N_r,N_\theta,N_\phi\}$ grid points, the boundaries are covered by ${\rm NGHOSTS}$ ghost zones, and we maintain cell centering.

Notice that in NRPy+, we use the [physics](https://en.wikipedia.org/wiki/Spherical_coordinate_system) notation for spherical coordinates, where $\theta$ is the polar and $\phi$ is the azimuthal angle. Also, we choose $\phi$ to range from $-\pi$ to $+\pi$, which is most useful since it is compatible with output from [`atan2`](https://en.wikipedia.org/wiki/Atan2).

**Exercise to student**: Given the prescription above, why do the integers $N_\theta$ and $N_\phi$ need to be even?

As Laplacians like these appear on the right-hand sides of, e.g., the scalar wave equation in curvilinear coordinates, we still have a problem of some terms becoming quite large as the coordinate singularity is approached. This issue manifests as a stiffening of the PDE, requiring that we be very careful about the precise [Method of Lines](Tutorial-Method_of_Lines-C_Code_Generation.ipynb) timestepping algorithm used. See [Cordero-Carrión & Cerdá-Durán](https://arxiv.org/abs/1211.5930) for information on dealing with this subtlety in a second-order Runge-Kutta Method of Lines context; it was later found that the standard RK4 method maintains stable solutions to PDEs affected by this sort of stiffening.

The above discussion focuses primarily on scalar fields. However, when solving PDEs involving vectors and tensors, the vectors and tensors themselves can exhibit divergent behavior at coordinate singularities. The good news is, this singular behavior is well-understood in terms of the scale factors of the reference metric, enabling us to define rescaled versions of these quantities that are well-behaved (so that, e.g., they can be finite differenced).

For example, given a smooth vector *in a 3D Cartesian basis* $\bar{\Lambda}^{i}$, all components $\bar{\Lambda}^{x}$, $\bar{\Lambda}^{y}$, and $\bar{\Lambda}^{z}$ will be smooth (by assumption). When changing the basis to spherical coordinates (applying the appropriate Jacobian matrix transformation), we will find that since $\phi = \arctan(y/x)$, $\bar{\Lambda}^{\phi}$ is given by

\begin{align}
\bar{\Lambda}^{\phi} &= \frac{\partial \phi}{\partial x} \bar{\Lambda}^{x} + 
\frac{\partial \phi}{\partial y} \bar{\Lambda}^{y} + 
\frac{\partial \phi}{\partial z} \bar{\Lambda}^{z} \\
&= -\frac{y}{x^2+y^2} \bar{\Lambda}^{x} + 
\frac{x}{x^2+y^2} \bar{\Lambda}^{y} \\
&= -\frac{y}{(r \sin\theta)^2} \bar{\Lambda}^{x} + 
\frac{x}{(r \sin\theta)^2} \bar{\Lambda}^{y} \\
&= -\frac{r \sin\theta \sin\phi}{(r \sin\theta)^2} \bar{\Lambda}^{x} + 
\frac{r \sin\theta \cos\phi}{(r \sin\theta)^2} \bar{\Lambda}^{y}\\
&= -\frac{\sin\phi}{r \sin\theta} \bar{\Lambda}^{x} + 
\frac{\cos\phi}{r \sin\theta} \bar{\Lambda}^{y}\\
\end{align}

Thus $\bar{\Lambda}^{\phi}$ diverges at all points where $r\sin\theta=0$ (or equivalently where $x=y=0$; i.e., the $z$-axis) due to the $\frac{1}{r\sin\theta}$ that appear in the Jacobian transformation.

This divergence might pose no problem on cell-centered grids that avoid $r \sin\theta=0$, except that the BSSN equations require that *first and second derivatives* of quantities like $\bar{\Lambda}^{\phi}$ be computed. Usual strategies for numerical approximation of these derivatives (e.g., finite difference methods) will "see" these divergences and errors generally will not drop to zero with increased numerical sampling of the functions at points near where the functions diverge.

However, notice that if we define $\lambda^{\phi}$ such that

$$\bar{\Lambda}^{\phi} = \frac{1}{r\sin\theta} \lambda^{\phi},$$

then $\lambda^{\phi}$ will be smooth and non-divergent as well. The strategy when computing derivatives of $\bar{\Lambda}^{\phi}$, therefore, is to perform the product rule on the above expression, computing derivatives of the scale factors *analytically* (i.e., exactly using a computer algebra system like SymPy), and smooth terms like $\lambda^{\phi}$ with finite-difference derivatives.

Avoiding such singularities can be generalized to arbitrary coordinate systems, so long as $\lambda^i$ is defined as:

$$\bar{\Lambda}^{i} = \frac{\lambda^i}{\text{scalefactor[i]}} ,$$

where scalefactor\[i\] is the $i$th scale factor in the given coordinate system. This idea can be extended to covariant (lowered-index) vectors and arbitrary tensors, as described in [the BSSN quantities tutorial notebook](Tutorial-BSSN_quantities.ipynb#rescaling_tensors).

**In summary, Challenge 4 is addressed via a combination of cell-centered grids, tensor rescaling, and a stable Method of Lines time-stepping algorithm. This tutorial notebook will therefore focus on addressing Challenges 1 through 3, which, coincidentally, are addressed via an appropriate boundary condition algorithm.**

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#basic_algorithm): Overview of boundary condition algorithm in curvilinear coordinates
    1. [Step 1.a](#challenge1): Addressing Challenge 1: Distinguishing inner from outer boundary points
    1. [Step 1.b](#challenge2): Addressing Challenge 2: Eigen-Coordinate systems
    1. [Step 1.c](#challenge3): Addressing Challenge 3: Applying parity conditions to arbitrary-rank tensors
1. [Step 2](#ccode_bc_struct): `bc_struct` & friends: Data structures for storing boundary condition information
1. [Step 3](#nrpycodegen): NRPy+-based C code generation for parity conditions
    1. [Step 3.a](#dotproducts): Python function `parity_conditions_symbolic_dot_products()`: Set up C code for computing unit-vector dot products (=parity) for each of the 10 parity condition types
    1. [Step 3.b](#set_parity_type): Populate `NRPy_basic_defines.h` with gridfunction parity types, which are set based on gridfunction name
1. [Step 4](#bcstruct_set_up): `bcstruct_set_up()`: C function for setting up `bcstruct`
    1. [Step 4.a](#bcstruct_set_up_eigencoord_mapping): `EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()`: Perform $(x,y,z)\to (x_0,x_1,x_2), \left(i_0,i_1,i_2\right)$ mappings using Eigen-Coordinate approach
    1. [Step 4.b](#bcstruct_set_up_inner_bc_parity): `set_parity_for_inner_boundary_single_pt()`: Perform $(x,y,z)\to \{(x_0,x_1,x_2), \left(i_0,i_1,i_2\right)\}$ mappings using Eigen-Coordinate approach
    1. [Step 4.c](#bcstruct_set_up_add_to_cfunc): Add `bcstruct_set_up()` to NRPy+ C function dictionary
1. [Step 5](#apply_bcs_inner_only): `apply_bcs_inner_only()`: Inner boundary condition C function
1. [Step 6](#outer_bcs): Outer boundary condition C functions
    1. [Step 6.a](#extrapolation_bcs) `"extrapolation"` outer boundary conditions: apply quadratic polynomial extrapolation
    1. [Step 6.b](#radiation_bcs): `"radiation"` outer boundary conditions
        1. [Step 6.b.i](#radiation_bcs_theory): Core *ansatz*, and implications
        1. [Step 6.b.ii](#radiation_bcs_numerical): Numerical implementation
        1. [Step 6.b.iii](#radiation_bcs_numerical_inv_jacobian): The $\partial_r f$ term: Computing $\frac{\partial x^i}{\partial r}$
        1. [Step 6.b.iv](#radiation_bcs_numerical_partial_i_f): The $\partial_r f$ term: Computing $\partial_i f$ with arbitrary-offset finite-difference derivatives
        1. [Step 6.b.v](#radiation_bcs_numerical_compute_partial_r_f): The $\partial_r f$ term: Putting it all together in `compute_partial_r_f()`
        1. [Step 6.b.vi](#radiation_bcs_evaluating_k): Evaluating the $k$ in the $k/r^3$ term
        1. [Step 6.b.vii](#radiation_bcs_apply_bcs_radiation): `radiation_bcs_single_pt()`: Apply radiation BCs to a single point
        1. [Step 6.b.viii](#apply_bcs_outerradiation_and_inner): `apply_bcs_outerradiation_and_inner()`: Apply radiation BCs at all outer boundary points and inner BCs at all inner boundary points
1. [Step 7](#start2finish): `CurviBC_Playground.c`: Start-to-Finish C code module for testing & validating curvilinear boundary conditions
    1. [Step 7.a](#register_gfs): Register gridfunctions of all 10 parity types in NRPy+
    1. [Step 7.b](#validate): Set up test data for Curvilinear Boundary Conditions code validation
    1. [Step 7.c](#mainc): `CurviBC_Playground`'s `main.c` Code
    1. [Step 7.d](#curvibc_setupall): Add all CurviBC C codes to C function dictionary, and add CurviBC definitions to `NRPy_basic_defines.h`
1. [Step 8](#add_cfunction_dicts): Add all C functions to `Cfunction_dict`, also output `NRPy_basic_defines.h` and `NRPy_function_prototypes.h`
1. [Step 9](#senr_compare): Validation: Compare with original SENR results
1. [Step 10](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='basic_algorithm'></a>

# Step 1: Overview of boundary condition algorithm in curvilinear coordinates \[Back to [top](#toc)\]
$$\label{basic_algorithm}$$

Here we **review** the basic algorithm for addressing Challenges [1](#innerouterbcs) [2](#coordinversion), and [3](#parity) discussed in the [Introduction](#intro) above. 

The algorithm itself is **implemented** as C code in Steps [**2**](#bc_struct) (data structures), [**3**](#set_up__bc_gz_map_and_parity_condns) (searching entire grid for inner and outer boundary points and setting parities), [**4**](#set_bcstruct) (setting data structures for quick and efficient implementation of outer boundaries), and [**5**](#apply_bcs_curvilinear) (function to apply inner & outer boundary conditions).

<a id='challenge1'></a>

## Step 1.a: Addressing Challenge 1: Distinguishing inner from outer boundary points \[Back to [top](#toc)\]
$$\label{challenge1}$$

At each ghost zone grid point $\mathbf{d}_{\rm gz}=(x_0,x_1,x_2)$, we will do the following:

1. Evaluate the Cartesian coordinate $\left(x(x_0,x_1,x_2),y(x_0,x_1,x_2),z(x_0,x_1,x_2)\right)$, corresponding to this grid point. Then evaluate the inverse $\mathbf{d}_{\rm new}=\left(x_0(x,y,z),x_1(x,y,z),x_2(x,y,z)\right)$. 
    1. If $\mathbf{d}_{\rm new} \ne \mathbf{d}_{\rm gz}$, then the ghost zone grid point maps to a point in the grid interior, *which is exactly the case described in the above section*. To distinguish this case from an "outer boundary condition", we shall henceforth refer to it variously as an application of an "interior", "inner", or "parity" boundary condition.
    1. Ghost zone points for which $\mathbf{d}_{\rm new} \equiv \mathbf{d}_{\rm gz}$ are on the outer boundary of the grid, and standard outer boundary conditions should be applied.

In detail, the algorithm is as follows:

1. Convert the coordinate $(x_0,x_1,x_2)$ for the ghost zone point to Cartesian coordinates $\left(x(x_0,x_1,x_2),y(x_0,x_1,x_2),z(x_0,x_1,x_2)\right)$. For example, if we choose ordinary spherical coordinates $(x_0,x_1,x_2)=(r,\theta,\phi)$, then
    + $x(r,\theta,\phi) = r \sin(\theta) \cos(\phi) = x_0 \sin(x_1) \cos(x_2)$
    + $y(r,\theta,\phi) = r \sin(\theta) \sin(\phi) = x_0 \sin(x_1) \sin(x_2)$
    + $z(r,\theta,\phi) = r \cos(\theta)            = x_0 \cos(x_1)$
1. Once we have $(x,y,z)$, we then find the corresponding value $(x_0,x_1,x_2)_{\rm in/OB}=(r,\theta,\phi)_{\rm in/OB}$ *in the grid interior or outer boundary*, via the simple inverse formula:
    + $r_{\rm in/OB}      = x_{0, \rm in/OB} = \sqrt{x^2+y^2+z^2} \in [0,\infty)$
    + $\theta_{\rm in/OB} = x_{1, \rm in/OB} = {\rm acos}\left(\frac{z}{\sqrt{x^2+y^2+z^2}}\right) \in [0,\pi]$
    + $\phi_{\rm in/OB}   = x_{2, \rm in/OB} = {\rm atan2}(y,x) \in [-\pi,\pi]$ [Wikipedia article on atan2()](https://en.wikipedia.org/w/index.php?title=Atan2&oldid=859313982)

1. If $(x_0,x_1,x_2)_{\rm in/OB}$ is the same as the original $(x_0,x_1,x_2)$, then we know $(x_0,x_1,x_2)$ is an outer boundary point (in spherical coordinates, at $r>{\rm RMAX}$), and we store `(i0,i1,i2)`$_{\rm in/OB} = (-1,-1,-1)$. Otherwise, we know that $(x_0,x_1,x_2)$ maps to some interior point at index `(i0,i1,i2)`, which we store:
    + $\rm{i0}_{\rm in/OB}=\frac{r_{\rm in/OB}      -      r_{\rm min}}{\Delta r}      - \frac{1}{2}$
    + $\rm{i1}_{\rm in/OB}=\frac{\theta_{\rm in/OB} - \theta_{\rm min}}{\Delta \theta} - \frac{1}{2}$
    + $\rm{i2}_{\rm in/OB}=\frac{\phi_{\rm in/OB}   -   \phi_{\rm min}}{\Delta \phi}   - \frac{1}{2}$

1. When updating a ghost zone point `(i0,i1,i2)` in the domain exterior, if the corresponding `(i0,i1,i2)`$_{\rm in/OB}$ was set to $(-1,-1,-1)$, then we apply outer boundary conditions. Otherwise, we simply copy the data from the interior point at `(i0,i1,i2)`$_{\rm in/OB}$ to `(i0,i1,i2)`.

Following the prescription in the [SENR/NRPy+ paper](https://arxiv.org/abs/1712.07658), we will implement curvilinear boundary conditions for rank-0, rank-1, and symmetric rank-2 tensors in three dimensions; as this is the same dimension and highest rank needed for BSSN.

<a id='challenge2'></a>

## Step 1.b: Addressing Challenge 2: Eigen-Coordinate Systems \[Back to [top](#toc)\]
$$\label{challenge2}$$

Suppose we were to rewrite the spherical coordinate $r$ as an arbitrary odd function of $x_0$ instead of $x_0$ itself. In that case, $r(-x_0)=-r(x_0)$, and all parity conditions remain unchanged. However, the inverse function, $x_0(r)$, may not be writable as a closed-form expression, requiring a Newton-Raphson root finder to find the appropriate boundary mappings. 

To greatly simplify the algorithm in the case of arbitrary $r(x_0)$ in Spherical-like coordinates, or $\rho(x_0)$ or $z(x_2)$ in Cylindrical-like coordinates, we note that the coordinate mappings *and* parities for all Spherical-like coordinate systems are identical to the mappings and parities for ordinary Spherical coordinates. The same holds true for Cylindrical-like and SymTP-like coordinate systems. Thus so long as we know the correct "Eigen-Coordinate system" (i.e., Spherical in the case of SinhSpherical or SinhSphericalv2; Cylindrical in the case of SinhCylindrical; SymTP in the case of SinhSymTP; etc.) there is no need for a Newton-Raphson root finder to set up the boundary conditions.

<a id='challenge3'></a>

## Step 1.c: Addressing Challenge 3: Applying parity conditions to arbitrary-rank tensors \[Back to [top](#toc)\]
$$\label{challenge3}$$

Above we presented the strategy for applying parity boundary conditions to a single component of a vector. Here we outline the generic algorithm for arbitrary-rank tensors.

Continuing the discussion from the previous section, we assume $\mathbf{d}_{\rm new} \ne \mathbf{d}_{\rm gz}$ (otherwise we would apply the *outer* boundary condition algorithm). Next suppose we are given a generic rank-$N$ tensor ($N>0$).

1. The first component of the rank-$N$ tensor corresponds to some direction with unit vector $\mathbf{e}^i$; e.g., $v^r$ corresponds to the $\mathbf{e}^r$ direction. Compute the dot product of the unit vector $\mathbf{e}^i$ evaluated at points $\mathbf{d}_{\rm gz}$ and $\mathbf{d}_{\rm new}$. Define this dot product as $P_1$ ("$P$" for "parity"):
$$
P_1 = \mathbf{e}^i\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^i\left(\mathbf{d}_{\rm new}\right).
$$
1. $P_1$ will take the value of $\pm 1$, depending on the unit-vector direction and the points $\mathbf{d}_{\rm gz}$ and $\mathbf{d}_{\rm new}$
1. Repeat the above for the remaining components of the rank-$N$ tensor $j\in \{2,3,...,N\}$, storing each $P_j$.
1. The tensor mapping from $\mathbf{d}_{\rm gz}$ to $\mathbf{d}_{\rm new}$ for this tensor $T^{ijk...}_{mnp...}$ will be given by
$$
T^{ijk...}_{lmn...}(x_0,x_1,x_2)_{\rm gz} = \prod_{\ell=1}^N P_\ell T^{ijk...}_{mnp...}(x_0,x_1,x_2)_{\rm new}.
$$

In this formulation of BSSN, we only need to deal with rank-0, rank-1, and *symmetric* rank-2 tensors. Further, our basis consists of 3 directions, so there are a total of 
+ 1 parity condition (the trivial +1) for scalars (rank-0 tensors)
+ 3 parity conditions for all rank-1 tensors (corresponding to each direction)
+ 6 parity conditions for all *symmetric* rank-2 tensors (corresponding to the number of elements in the lower or upper triangle of a $3\times3$ matrix, including the diagonal)

Thus we must keep track of the behavior of **10 separate parity conditions**, which can be evaluated once the numerical grid has been set up, for all time. The following Table outlines the correct conditions for each:

The appropriate dot products determining parity condition are assigned to each gridfunction based on the following numbering:

Tensor type | Parity type | Dot product(s) determining parity condition (see equation above)
--- | --- | ---
Scalar (Rank-0 tensor) | 0 | (*none*)
Rank-1 tensor in **i0** direction | 1 | $\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)$
Rank-1 tensor in **i1** direction | 2 | $\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)$
Rank-1 tensor in **i2** direction | 3 | $\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)$
Rank-2 tensor in **i0-i0** direction | 4 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$
Rank-2 tensor in **i0-i1** direction | 5 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i0-i2** direction | 6 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i1-i1** direction | 7 | $\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$
Rank-2 tensor in **i1-i2** direction | 8 | $\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i2-i2** direction | 9 | $\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$


In the following few steps, we will document the data structures for and implementation of this boundary condition algorithm.

<a id='ccode_bc_struct'></a>

# Step 2:  `bc_struct` & friends: Data structures for storing boundary condition information \[Back to [top](#toc)\]
$$\label{ccode_bc_struct}$$

Here we define `bc_struct`, a C data structure that stores all information needed to apply boundary conditions at all boundary points. 

Information needed to fill each ghost zone is based on whether the ghost zone is an inner or an outer boundary point.

**`innerpt_bc_struct`: C struct for inner boundary points**

If the ghost zone point is an inner boundary point, it maps to *a different point* either on the grid outer boundary or the grid interior. Further, vectors or tensors will flip sign based on their parity across this boundary (see section above). To fill in each inner boundary point we store data to the C struct `inner_bc_struct`:

```C
typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;
```

**`outerpt_bc_struct`: C struct for outer boundary points**

Unlike inner boundary points, ghost zone points that are outer boundary points map to *themselves*, and are filled using outer boundary conditions. Further, outer boundary ghost zones must be filled in from the inside outward, as e.g., the outermost ghost zones may depend on ghost zones closer to the grid interior being set. Thus filling in outer boundary points generally requires data at neighboring gridpoints *in the direction of the grid interior*. For example in Cartesian coordinates, filling in the $x=x_{\rm max}$ face of the grid depends on data at $x<x_{\rm max}$. Similarly, filling in the $y=y_{\rm min}$ face will require data at $y>y_{\rm min}$.

So some indication of not only the gridpoint `(i0,i1,i2)` but also the face on which it exists, must be stored. We set the latter to the 1-byte integers `FACEX0`, `FACEX1`, and `FACEX2`, such that if we are on the `x0=x0_max` face, `FACEX0=-1`, `FACEX1=FACEX2=0`, to indicate that `x0=x0_max` points will be filled in using data at `x0<x0_max`. Thus to fill in each outer boundary point we store data to the C struct `outer_bc_struct`:

```C
typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;
```

**`bc_info_struct`: C struct for other boundary information**

The `bc_info_struct` simply stores the total number of inner and outer boundary points. To be scalable when updating boundary conditions on multiple cores, as many boundary points are updated simultaneously as possible. Inner boundary points can all be updated simultaneously, while only one ghostzone level of min and max faces of outer boundary points can be updated simultaneously.

```C
typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;
```

**`bc_struct`: A struct of structs for storing all boundary condition information**

As described above, 

```C
typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;  // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
```


```python
# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC      # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict  # NRPy+: Core C code output module
from outputC import add_to_Cfunction_dict  # NRPy+: Core C code output module
from outputC import Cfunction    # NRPy+: Core C code output module
from outputC import indent_Ccode # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import finite_difference as fin  # NRPy+: Finite-difference module
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions
from UnitTesting.assert_equal import check_zero  # NRPy+: Checks whether an expression evaluates to zero.

# Step P2: Create C code output directory:
Ccodesrootdir = os.path.join("CurviBoundaryConditions_Ccodes/")
# P2.a: First remove C code output directory if it exists
# Courtesy https://stackoverflow.com/questions/303200/how-do-i-remove-delete-a-folder-that-is-not-empty
shutil.rmtree(Ccodesrootdir, ignore_errors=True)
# P2.b: Then create a fresh directory
# Then create a fresh directory
cmd.mkdir(Ccodesrootdir)

# CoordSystem = "Cylindrical"
CoordSystem = "SinhSpherical"
outer_bcs_type = "extrapolation"
# outer_bcs_type = "radiation"
RADIATION_BC_FD_ORDER = 4  # Set to -1 to adopt the same value as finite_difference::FD_CENTDERIVS_ORDER

# Set the finite differencing order to 4; although this module doesn't compute
#   finite difference derivatives, it does need to set the number of ghost zone
#   cells, which is generally based on NGHOSTS, which itself depends
#   on finite_difference::FD_CENTDERIVS_ORDER.
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

_unused = par.Cparameters("int", __name__, "outer_bc_type", "EXTRAPOLATION_OUTER_BCS")
```

Define/declare core data structures for curvilinear BCs within `NRPy_basic_defines.h`:


```python
# First register basic C data structures/macros inside NRPy_basic_defines.h
def NRPy_basic_defines_CurviBC_data_structures():
    return r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

#define EXTRAPOLATION_OUTER_BCS 0  // used to identify/specify params.outer_bc_type
#define RADIATION_OUTER_BCS     1  // used to identify/specify params.outer_bc_type

typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;  // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""
```

<a id='nrpycodegen'></a>

#  Step 3: NRPy+-based C code generation for parity conditions \[Back to [top](#toc)\]
$$\label{nrpycodegen}$$

Much of the algorithm needed for setting up `bcstruct` requires a loop over all gridpoints on the numerical grid. As the precise numerical grids are chosen at C runtime, that part of the algorithm must be run entirely within a static C code.

However, there are two parts to the overall algorithm that must be generated by NRPy+, namely

1. [Step 3.a](#dotproducts): `parity_conditions_symbolic_dot_products()`: Based on the chosen reference metric, sets up the needed unit-vector dot products for each of the 10 parity condition types.
1. [Step 3.b](#set_parity_type): Set parity type for each gridfunction registered within NRPy+, based on the digits at the end of each gridfunction name, append result to `dirname+gridfunction_defines.h`

<a id='dotproducts'></a>

## Step 3.a: Python function `parity_conditions_symbolic_dot_products()`: Set up C code for computing unit-vector dot products (=parity) for each of the 10 parity condition types \[Back to [top](#toc)\]
$$\label{dotproducts}$$

Next, we generate the C code necessary to perform needed dot products for filling in the parity condition arrays inside `bcstruct`. 

Using the unit vectors defined in `rfm.UnitVectors[][]` (in `reference_metric.py`), each unit vector takes as input either $\mathbf{d}_{\rm gz} = (x_0,x_1,x_2)_{\rm IB}$=`(xx0,xx1,xx2)` or $\mathbf{d}_{\rm new} = (x_0,x_1,x_2)_{\rm in}$=`(xx0_inbounds,xx1_inbounds,xx2_inbounds)` as summarized in the table above in [Step 1](#challenge2). We paste the table here again, for quick reference. Pay special attention to the parity type numbering, as this is the convention adopted in the code:

Tensor type | Parity type | Dot product(s) determining parity condition (
--- | --- | ---
Scalar (Rank-0 tensor) | 0 | (*none*)
Rank-1 tensor in **i0** direction | 1 | $\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)$
Rank-1 tensor in **i1** direction | 2 | $\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)$
Rank-1 tensor in **i2** direction | 3 | $\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)$
Rank-2 tensor in **i0-i0** direction | 4 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$
Rank-2 tensor in **i0-i1** direction | 5 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i0-i2** direction | 6 | $\left[\mathbf{e}^0\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^0\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i1-i1** direction | 7 | $\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$
Rank-2 tensor in **i1-i2** direction | 8 | $\left[\mathbf{e}^1\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^1\left(\mathbf{d}_{\rm new}\right)\right]\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]$
Rank-2 tensor in **i2-i2** direction | 9 | $\left[\mathbf{e}^2\left(\mathbf{d}_{\rm gz}\right) \cdot \mathbf{e}^2\left(\mathbf{d}_{\rm new}\right)\right]^2 = 1$

Looping over all 10 parity types, the corresponding symbolic expressions for dot product(s) is output to C code. For example, in Spherical coordinates, parity type 1's dot product output as C code is given by: 

```
parity[1] = sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds);
```

To wit (as described above), there are 10 parity types for BSSN evolved variables, which include scalars, vectors, and symmetric rank-2 tensors.


```python
# Set unit-vector dot products (=parity) for each of the 10 parity condition types
def parity_conditions_symbolic_dot_products():
    parity = ixp.zerorank1(DIM=10)
    UnitVectors_inner = ixp.zerorank2()
    xx0_inbounds,xx1_inbounds,xx2_inbounds = sp.symbols("xx0_inbounds xx1_inbounds xx2_inbounds", real=True)
    for i in range(3):
        for j in range(3):
            UnitVectors_inner[i][j] = rfm.UnitVectors[i][j].subs(rfm.xx[0],xx0_inbounds).subs(rfm.xx[1],xx1_inbounds).subs(rfm.xx[2],xx2_inbounds)
    # Type 0: scalar
    parity[0] = sp.sympify(1)
    # Type 1: i0-direction vector or one-form
    # Type 2: i1-direction vector or one-form
    # Type 3: i2-direction vector or one-form
    for i in range(3):
        for Type in range(1,4):
            parity[Type] += rfm.UnitVectors[Type-1][i]*UnitVectors_inner[Type-1][i]
    # Type 4: i0i0-direction rank-2 tensor
    # parity[4] = parity[1]*parity[1]
    # Type 5: i0i1-direction rank-2 tensor
    # Type 6: i0i2-direction rank-2 tensor
    # Type 7: i1i1-direction rank-2 tensor
    # Type 8: i1i2-direction rank-2 tensor
    # Type 9: i2i2-direction rank-2 tensor
    count = 4
    for i in range(3):
        for j in range(i,3):
            parity[count] = parity[i+1]*parity[j+1]
            count = count + 1

    lhs_strings = []
    for i in range(10):
        lhs_strings.append("REAL_parity_array["+str(i)+"]")
    outstr = """
    // NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
    //      ten parity conditions, in given coordinate system.
    //      Needed for automatically determining sign of tensor across coordinate boundary.
    // Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
"""
    return outstr + outputC(parity, lhs_strings, filename="returnstring", params="preindent=2")

#     print("\n\nExample: parity type 1's dot product is given by: \n"+lhs_strings[1]+" = "+str(parity[1]))
```

<a id='set_parity_type'></a>

## Step 3.b: Populate `NRPy_basic_defines.h` with gridfunction parity types, which are set based on gridfunction name \[Back to [top](#toc)\]
$$\label{set_parity_type}$$

For example, if the gridfunction name ends with "01", then (based on the table above) the `set_parity_types()` function below will set the parity_type of that gridfunction to 5. We can be assured this is a rather robust algorithm, because `gri.register_gridfunctions()` in [grid.py](../edit/grid.py) will throw an error if a gridfunction's **base** name ends in an integer. This strict syntax was added with the express purpose of making it easier to set parity types based solely on the gridfunction name.

After each parity type is found, we store the parity type of each gridfunction to `const int8_t arrays` `evol_gf_parity` and `aux_gf_parity`, appended to the end of `NRPy_basic_defines.h`.


```python
def NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=True):
    # First add human-readable gridfunction aliases (grid.py) to NRPy_basic_defines dictionary,
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = gri.gridfunction_lists()[0:3]

    # Step 3.b: set the parity conditions on all gridfunctions in gf_list,
    #       based on how many digits are at the end of their names
    def set_parity_types(list_of_gf_names):
        parity_type = []
        for name in list_of_gf_names:
            for gf in gri.glb_gridfcs_list:
                if gf.name == name:
                    parity_type__orig_len = len(parity_type)
                    if gf.DIM < 3 or gf.DIM > 4:
                        print("Error: Cannot currently specify parity conditions on gridfunctions with DIM<3 or >4.")
                        sys.exit(1)
                    if gf.rank == 0:
                        parity_type.append(0)
                    elif gf.rank == 1:
                        if gf.DIM == 3:
                            parity_type.append(int(gf.name[-1]) + 1)  # = 1 for e.g., beta^0; = 2 for e.g., beta^1, etc.
                        elif gf.DIM == 4:
                            parity_type.append(int(gf.name[-1]))  # = 0 for e.g., b4^0; = 1 for e.g., beta^1, etc.
                    elif gf.rank == 2:
                        if gf.DIM == 3:
                            # element of a list; a[-2] the
                            # second-to-last element, etc.
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            if idx0 == "0" and idx1 == "0":
                                parity_type.append(4)
                            elif (idx0 == "0" and idx1 == "1") or (idx0 == "1" and idx1 == "0"):
                                parity_type.append(5)
                            elif (idx0 == "0" and idx1 == "2") or (idx0 == "2" and idx1 == "0"):
                                parity_type.append(6)
                            elif idx0 == "1" and idx1 == "1":
                                parity_type.append(7)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(8)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(9)
                        elif gf.DIM == 4:
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            # g4DD00 = g_{tt} : parity type = 0
                            # g4DD01 = g_{tx} : parity type = 1
                            # g4DD02 = g_{ty} : parity type = 2
                            # g4DD0a = g_{ta} : parity type = a
                            if idx0 == "0":
                                parity_type.append(int(idx1))
                            elif idx1 == "0":
                                parity_type.append(int(idx0))
                            if idx0 == "1" and idx1 == "1":
                                parity_type.append(4)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(5)
                            elif (idx0 == "1" and idx1 == "3") or (idx0 == "3" and idx1 == "1"):
                                parity_type.append(6)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(7)
                            elif (idx0 == "2" and idx1 == "3") or (idx0 == "3" and idx1 == "2"):
                                parity_type.append(8)
                            elif idx0 == "3" and idx1 == "3":
                                parity_type.append(9)
                    if len(parity_type) == parity_type__orig_len:
                        print("Error: Could not figure out parity type for "+gf.gftype+" gridfunction: " + gf.name,gf.DIM,gf.name[-2],gf.name[-1],gf.rank)
                        sys.exit(1)
        if len(parity_type) != len(list_of_gf_names):
            print("Error: For some reason the length of the parity types list did not match the length of the gf list.")
            sys.exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)
    auxevol_parity_type = set_parity_types(auxevol_variables_list)

    # Output all gridfunctions to Ccodesrootdir/gridfunction_defines.h
    # ... then append to the file the parity type for each gridfunction.
    outstr = """
/* PARITY TYPES FOR ALL GRIDFUNCTIONS.
 * SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */
"""
    if len(evolved_variables_list) > 0:
        outstr += "static const int8_t evol_gf_parity[" + str(len(evolved_variables_list)) + "] = { "
        for i in range(len(evolved_variables_list) - 1):
            outstr += str(evol_parity_type[i]) + ", "
        outstr += str(evol_parity_type[len(evolved_variables_list) - 1]) + " };\n"

    if len(auxiliary_variables_list) > 0:
        outstr += "static const int8_t aux_gf_parity[" + str(len(auxiliary_variables_list)) + "] = { "
        for i in range(len(auxiliary_variables_list) - 1):
            outstr += str(aux_parity_type[i]) + ", "
        outstr += str(aux_parity_type[len(auxiliary_variables_list) - 1]) + " };\n"

    if len(auxevol_variables_list) > 0:
        outstr += "static const int8_t auxevol_gf_parity[" + str(len(auxevol_variables_list)) + "] = { "
        for i in range(len(auxevol_variables_list) - 1):
            outstr += str(auxevol_parity_type[i]) + ", "
        outstr += str(auxevol_parity_type[len(auxevol_variables_list) - 1]) + " };\n"

    if verbose == True:
        for i in range(len(evolved_variables_list)):
            print("Evolved gridfunction \"" + evolved_variables_list[i] + "\" has parity type " + str(
                evol_parity_type[i]) + ".")
        for i in range(len(auxiliary_variables_list)):
            print("Auxiliary gridfunction \"" + auxiliary_variables_list[i] + "\" has parity type " + str(
                aux_parity_type[i]) + ".")
        for i in range(len(auxevol_variables_list)):
            print("AuxEvol gridfunction \"" + auxevol_variables_list[i] + "\" has parity type " + str(
                auxevol_parity_type[i]) + ".")
    return outstr
```

<a id='bcstruct_set_up'></a>

# Step 4: `bcstruct_set_up()`: C function for setting up `bcstruct` \[Back to [top](#toc)\]
$$\label{bcstruct_set_up}$$

`bcstruct_set_up()` sets up `bcstruct`, the C struct that contains all information needed to fill in grid boundaries. It depends on two C helper functions, which are included in `bcstruct_set_up.c`:

1. `EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds()` uses the [Eigen-Coordinate](#challenge2) approach described above to find the mapping of boundary point `(x0,x1,x2)` $\to$ `(x,y,z)` $\to$ `(x0',x1',x2')`. If `(x0,x1,x2)` is an inner boundary point, then `(x0',x1',x2')` $\equiv$ `(x0,x1,x2)_inbounds` will map to a different gridpoint-either in the grid interior or at an outer boundary point. I.e., when `(x0,x1,x2)`$\ne$`(x0',x1',x2')`, `(x0,x1,x2)` is an *inner* boundary gridpoint. Otherwise the boundary point `(x0,x1,x2)` is an *outer* boundary gridpoint.
1. `set_parity_for_inner_boundary_single_pt()` sets the parity condition for inner boundary points.

<a id='bcstruct_set_up_eigencoord_mapping'></a>

## Step 4.a:  `EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()`: Perform $(x,y,z)\to (x_0,x_1,x_2), \left(i_0,i_1,i_2\right)$ mappings using Eigen-Coordinate approach \[Back to [top](#toc)\]
$$\label{bcstruct_set_up_eigencoord_mapping}$$


First we generate the C code needed for applying boundary conditions in generic coordinate systems, using the [Eigen-Coordinate](#challenge2) approach described above:
* $\left(x(x_0,x_1,x_2),y(x_0,x_1,x_2),z(x_0,x_1,x_2)\right) \to \left(x_0(x,y,z),x_1(x,y,z),x_2(x,y,z)\right), \left(i_0,i_1,i_2\right)$:


```python
def Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
    desc = """EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
  A coordinate system's "eigencoordinate" is the simplest member
  of its family; all spherical-like coordinate systems have
  Spherical as their eigencoordinate. The same is true for
  cylindrical-like (Cylindrical is eigencoordinate),
  Cartesian-like (Cartesian is the eigencoordinate), and
  SymTP-like (SymTP is the eigencoordinate) coordinates.

  For a given gridpoint (i0,i1,i2) and corresponding coordinate
  (x0,x1,x2), this function performs the dual mapping
  (x0,x1,x2) -> (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  Note that (x0,x1,x2) IS NOT ALWAYS equal to (x0,x1,x2)';
  For example consider in Spherical coordinates
  (x0,x1,x2)=(r,theta,phi)=(-0.1,pi/4,pi/4).
  This point will map to (x0,x1,x2)', in which x0>0,
  because the inversion r=sqrt(Cartx^2+Carty^2+Cartz^2)
  is always positive. In this case, (x0,x1,x2) is considered
  an *inner* boundary point, and on a cell-centered grid
  is guaranteed to map to a grid point in the grid interior;
  filling in this point requires copying data, and possibly
  multiplying by a +/- 1 if the data is from a gridfunction
  storing tensors/vectors.
"""
    c_type = "static void"
    name = "EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt"
    params = """const paramstruct *restrict params, REAL *restrict xx[3],
                                                               const int i0, const int i1, const int i2,
                                                               REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]"""
    body = r"""
  // This is a 3-step algorithm:
  // Step 1: (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  // Step 2: (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  //         Find the interior eigencoordinate point
  //         (x0,x1,x2)' to which (Cartx,Carty,Cartz)
  //         maps, as well as the corresponding
  //         gridpoint integer index (i0,i1,i2). For
  //         cell-centered grids, (x0,x1,x2) will always
  //         overlap exactly (to roundoff error) a point
  //         on the numerical grid.
  // Step 3: Sanity check
  //         Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!
"""
    # Load up the EigenCoordinate corresponding to reference_metric::CoordSystem
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem",rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 1: Output C code for the Eigen-Coordinate mapping from xx->Cartesian':
    body += r"""
  // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
  REAL xCart[3];  // where (x,y,z) is output
  {
    // xx_to_Cart for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart[0]","xCart[1]","xCart[2]"],
            "returnstring", params="preindent=2")+"  }\n"
    body += r"""
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];"""

    # Step 2: Output C code for the Eigen-Coordinate mapping from Cartesian->xx':
    body += r"""
  // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
  //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
  //      the point is an outer boundary point.
  //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
  //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
  //      appropriate parity condition (+/- 1).
  REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
"""
    # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        print("ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for")
        print("       reference_metric::CoordSystem = "+par.parval_from_str("reference_metric::CoordSystem"))
        print("       Boundary conditions in curvilinear coordinates REQUiRE this be set.")
        sys.exit(1)
    # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    body += """  // Cart_to_xx for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+");\n"
    body += outputC([rfm.Cart_to_xx[0],rfm.Cart_to_xx[1],rfm.Cart_to_xx[2]],
                    ["Cart_to_xx0_inbounds","Cart_to_xx1_inbounds","Cart_to_xx2_inbounds"],
                    filename="returnstring", params="preindent=1")
    body += r"""
  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const REAL xxmin[3] = {
    xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
    xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
    xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };

  // Finally compute i{0,1,2}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );
"""

    # Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
    rfm.reference_metric()

    # Step 3:
    body += """
  // Step 3: Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!

  // Step 3.a: Compute x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds):
  const REAL x0_inbounds = xx[0][i0_inbounds];
  const REAL x1_inbounds = xx[1][i1_inbounds];
  const REAL x2_inbounds = xx[2][i2_inbounds];

  // Step 3.b: Compute {x,y,z}Cart_from_xx, as a
  //           function of i0,i1,i2
  REAL xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {
    // xx_to_Cart for Coordinate """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart_from_xx","yCart_from_xx","zCart_from_xx"],
            "returnstring", params="preindent=2,includebraces=False")
    body += r"""  }

  // Step 3.c: Compute {x,y,z}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  REAL xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {
    // xx_to_Cart_inbounds for Coordinate """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0_inbounds];
    REAL xx1 = xx[1][i1_inbounds];
    REAL xx2 = xx[2][i2_inbounds];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart_from_xx_inbounds","yCart_from_xx_inbounds","zCart_from_xx_inbounds"],
            "returnstring", params="preindent=2,includebraces=False")
    body += r"""  }

  // Step 3.d: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
#define EPS_REL 1e-8
  const REAL norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (double)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {
    fprintf(stderr,"Error in """+CoordSystem_orig+r""" coordinate system: Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (double)xCart_from_xx,(double)yCart_from_xx,(double)zCart_from_xx,
            (double)xCart_from_xx_inbounds,(double)yCart_from_xx_inbounds,(double)zCart_from_xx_inbounds,
            xx[0][i0],xx[1][i1],xx[2][i2],
            xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
    exit(1);
  }

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
"""
    __function_prototype_ignore, Cfunc = Cfunction(
        desc=desc, c_type=c_type, name=name, params=params,
        body=body)
    return Cfunc
```

<a id='bcstruct_set_up_inner_bc_parity'></a>

## Step 4.b:  `set_parity_for_inner_boundary_single_pt()`: Perform $(x,y,z)\to \{(x_0,x_1,x_2), \left(i_0,i_1,i_2\right)\}$ mappings using Eigen-Coordinate approach \[Back to [top](#toc)\]
$$\label{bcstruct_set_up_inner_bc_parity}$$

Next, we define the C function `set_parity_for_inner_boundary_single_pt()`, which implements the algorithm described above in [Challenge 3](#challenge3).

Given the inputs: 
* inner boundary point $(x_0,x_1,x_2)$ situated at grid index `(i0,i1,i2)`, and
* the interior (non-boundary) point to which it maps $(x0,x1,x2)'=$ `(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])`

this function evaluates dot products of the unit vectors evaluated at `(i0_inbounds,i1_inbounds,i2_inbounds)` and `(i0,i1,i2)`, for all 10 gridfunction parities currently supported in NRPy+. The C code for computing the needed symbolic dot products is generated above in [Step 3](#dotproducts).


```python
def Cfunction__set_parity_for_inner_boundary_single_pt():
    desc = """set_parity_for_inner_boundary_single_pt():
  Given (x0,x1,x2)=(xx0,xx1,xx2) and
  (x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
  (see description of
  EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
  above for more details), here we compute the parity conditions
  for all 10 tensor types supported by NRPy+.
"""
    c_type = "static void"
    name = "set_parity_for_inner_boundary_single_pt"
    params = """const paramstruct *restrict params, const REAL xx0,const REAL xx1,const REAL xx2,
                                             const REAL x0x1x2_inbounds[3], const int idx,
                                             innerpt_bc_struct *restrict innerpt_bc_arr"""
    body = r"""
  const REAL xx0_inbounds = x0x1x2_inbounds[0];
  const REAL xx1_inbounds = x0x1x2_inbounds[1];
  const REAL xx2_inbounds = x0x1x2_inbounds[2];

  REAL REAL_parity_array[10];
  {
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+
""" + parity_conditions_symbolic_dot_products() + r"""
  }

  // Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
  for(int whichparity=0;whichparity<10;whichparity++) {
    //printf("Good? Parity %d evaluated to %e\n",whichparity,(double)REAL_parity_array[whichparity]);
    if( fabs(REAL_parity_array[whichparity]) < 1 - 1e-8 || fabs(REAL_parity_array[whichparity]) > 1 + 1e-8 ) {
      fprintf(stderr,"Error at point (%e %e %e), which maps to (%e %e %e).\n",
              xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
      fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
              REAL_parity_array[whichparity]);
      exit(1);
    }
    // The typecast (int8_t)REAL_parity_array[whichparity] *does not work*.
    //  Thankfully, we've already checked whether REAL_parity_array[whichparity]
    //  is within 1e-8 of +/- 1, so here we just check the sign of the
    //  REAL_parity_array to find the correct value of innerpt_bc_arr[idx].parity[parity].
    for(int parity=0;parity<10;parity++) {
      innerpt_bc_arr[idx].parity[parity] = 1;
      if(REAL_parity_array[parity] < 0) innerpt_bc_arr[idx].parity[parity] = -1;
    }
  } // END for(int whichparity=0;whichparity<10;whichparity++)
"""
    __function_prototype_ignore, Cfunc = Cfunction(
        desc=desc, c_type=c_type, name=name, params=params,
        body=body)
    return Cfunc
```

<a id='bcstruct_set_up_add_to_cfunc'></a>

## Step 4.c: Add `bcstruct_set_up()` to NRPy+ C function dictionary \[Back to [top](#toc)\]
$$\label{bcstruct_set_up_add_to_cfunc}$$

With the core routines `EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt` and `set_parity_for_inner_boundary_single_pt()` now defined, we can now add `bcstruct_set_up()` to the NRPy+ C function dictionary.

At each inner and outer boundary point, this function sets the `innerpt_bc_struct` and `outerpt_bc_struct`, respectively. Note that

1. Inner boundary points always map to a different grid point. An inner boundary point falls into one of the following two categories
    1. pure inner boundary points (that map to a grid point entirely within the grid interior), or
    1. inner-maps-to-outer boundary points (inner boundary points that map to a grid point on the outer boundary).
1. Outer boundary points always map to themselves.


```python
def add_to_Cfunction_dict_bcstruct_set_up(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    prefunc  = Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
    prefunc += Cfunction__set_parity_for_inner_boundary_single_pt()
    desc = r"""At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):

Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
  Recall that at each inner boundary point we must set innerpt_bc_struct:
    typedef struct __innerpt_bc_struct__ {
      int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
      int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
      int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
    } innerpt_bc_struct;
  At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
    Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
        This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
        Cartesian coordinate (x,y,z), then finds the grid point
        (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
        corresponding to this Cartesian coordinate (x,y,z).
    If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an inner boundary point. We must set
        Set bcstruct->inner_bc_array for this point, which requires we specify
        both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
        conditions for this gridpoint. The latter is found & specified within the
        function set_parity_for_inner_boundary_single_pt().
    If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an outer boundary point. Take care of outer BCs in Step 2.
Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
  Recall that at each inner boundary point we must set outerpt_bc_struct:
    typedef struct __outerpt_bc_struct__ {
      short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
      int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
      //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
      //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
      //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
      //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
      //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
      //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
    } outerpt_bc_struct;
  Outer boundary points are filled from the inside out, two faces at a time.
    E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
    including the ghostzones, with NGHOSTS=2.
    We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
    points in first, since they will in general (at least in the case of extrapolation
    outer BCs) depend on e.g., i0=2 and i0=3 points.
    Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
    since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
    Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
    depend on i0 min and max faces being filled. The remaining pattern goes like this:
    Upper x1 face: (i0={1,12},i1=12,i2={2,11})
    Lower x2 face: (i0={1,12},i1={1,12},i2=1)
    Upper x2 face: (i0={1,12},i1={1,12},i2=12)
    Lower x0 face: (i0=0,i1={1,12},i2={1,12})
    Upper x0 face: (i0=13,i1={1,12},i2={1,12})
    Lower x1 face: (i0={0,13},i1=0,i2={2,11})
    Upper x1 face: (i0={0,13},i1=13,i2={2,11})
    Lower x2 face: (i0={0,13},i1={0,13},i2=0)
    Upper x2 face: (i0={0,13},i1={0,13},i2=13)
  Note that we allocate a outerpt_bc_struct at *all* boundary points,
    regardless of whether the point is an outer or inner point. However
    the struct is set only at outer boundary points. This is slightly
    wasteful, but only in memory, not in CPU.
"""
    c_type = "void"
    name = "bcstruct_set_up"
    params = "const paramstruct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct"
    body = r"""
  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)",
             i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        }
      }
    }
    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;

    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc( sizeof(innerpt_bc_struct)*num_inner );
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3S(i0,i1,i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3S(i0i1i2_inbounds[0],i0i1i2_inbounds[1],i0i1i2_inbounds[2]);
          //printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(params, xx[0][i0],xx[1][i1],xx[2][i2],
                                                  x0x1x2_inbounds, which_inner, bcstruct->inner_bc_array);

          which_inner++;
        }
      }
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  // First set up loop bounds for outer boundary condition updates,
  //   store to bc_info->bc_loop_bounds[which_gz][face][]. Also
  //   allocate memory for outer_bc_array[which_gz][face][]:
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) {
    const int x0min_face_range[6] = { imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2] };  imin[0]--;
    const int x0max_face_range[6] = { imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2] };  imax[0]++;
    const int x1min_face_range[6] = { imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2] };  imin[1]--;
    const int x1max_face_range[6] = { imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2] };  imax[1]++;
    const int x2min_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2] };  imin[2]--;
    const int x2max_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1 };  imax[2]++;

    int face=0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x0min_face_range[1]-x0min_face_range[0]) *
                                                                                              (x0min_face_range[3]-x0min_face_range[2]) *
                                                                                              (x0min_face_range[5]-x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i]; }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i]; }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x1min_face_range[1]-x1min_face_range[0]) *
                                                                                              (x1min_face_range[3]-x1min_face_range[2]) *
                                                                                              (x1min_face_range[5]-x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i]; }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i]; }
    face++;
    ////////////////////////


    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3*which_gz + face/2] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                             ((x2min_face_range[1]-x2min_face_range[0]) *
                                                                                              (x2min_face_range[3]-x2min_face_range[2]) *
                                                                                              (x2min_face_range[5]-x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i]; }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i]; }
    face++;
    ////////////////////////
  }

  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn*2;
#define IDX2D_BCS(i0,i0min,i0max, i1,i1min,i1max ,i2,i2min,i2max)       \
        ( ((i0)-(i0min)) + ((i0max)-(i0min)) * ( ((i1)-(i1min)) + ((i1max)-(i1min)) * ((i2)-(i2min)) ) )
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn*2+1;
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    }
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
```

<a id='apply_bcs_inner_only'></a>

# Step 5: `apply_bcs_inner_only()`: Inner boundary condition C function \[Back to [top](#toc)\]
$$\label{apply_bcs_inner_only}$$

After `bcstruct_set_up()`, `bcstruct` contains all information needed to fill inner boundary points. This includes

1. The interior or outer boundary point to which the inner point maps.
1. The tensor parity multiplier (+1 or -1) that the value of the gridfunction will need to be multiplied.

This function will generally need to be called after outer boundary conditions have been applied, as inner boundary points can depend on outer boundary points having been filled.


```python
def add_to_Cfunction_dict_apply_bcs_inner_only(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h")]
    desc = r"""
Apply BCs to inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    c_type = "void"
    name = "apply_bcs_inner_only"
    params = "const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  // collapse(2) results in a nice speedup here, esp in 2D. Two_BHs_collide goes from
  //    5550 M/hr to 7264 M/hr on a Ryzen 9 5950X running on all 16 cores with core affinity.
#pragma omp parallel for collapse(2)  // spawn threads and distribute across them
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      gfs[IDX4ptS(which_gf, dstpt)] = bcstruct->inner_bc_array[pt].parity[evol_gf_parity[which_gf]] * gfs[IDX4ptS(which_gf, srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++)
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
```

<a id='outer_bcs'></a>

# Step 6: Outer boundary condition C functions \[Back to [top](#toc)\]
$$\label{outer_bcs}$$

<a id='extrapolation_bcs'></a>

## Step 6.a: `"extrapolation"` outer boundary conditions: apply quadratic polynomial extrapolation \[Back to [top](#toc)\]
$$\label{extrapolation_bcs}$$

As an option, quadratic extrapolation may be applied to each outer boundary point `(i0,i1,i2)`, as follows.

Suppose the outer boundary point is at the `i0=max(i0)` face. Then we fit *known* data at `i0-3`, `i0-2`, and `i0-1` [i.e., $f_{-3}=f(x_{i0-3}=x_{-3})$, $f_{-2}=f(x_{i0-2}=x_{-2})$, and $f_{-1}=f(x_{i0-1}=x_{-1})$] to the unique quadratic polynomial:

\begin{align}
f_{-3} &= c_2 x_{-3}^2 + c_1 x_{-3} + c_0 \\
f_{-2} &= c_2 x_{-2}^2 + c_1 x_{-2} + c_0 \\
f_{-1} &= c_2 x_{-1}^2 + c_1 x_{-1} + c_0 \\
\end{align}

We wish to extrapolate to $f_0=f(x_0)$. Since our grid has uniform spacing, 

* $x_{-3}=x_0-3\Delta x$, 
* $x_{-2}=x_0-2\Delta x$, and
* $x_{-1}=x_0-\Delta x$.

The extrapolated value $f_0$ cannot depend on the choice of the fiducial $x_0$ (i.e., it will hold for *any* choice of $x_0$), so without loss of generality we will set $x_0=0$:

$$
\mathbf{A c} =
\left[
\begin{array}{ccc}
 1 & x_{-3}  & x_{-3}^2 \\
 1 & x_{-2}  & x_{-2}^2 \\
 1 & x_{-1}  & x_{-1}^2 \\
\end{array}
\right]
\left[
\begin{array}{c}
 c_0 \\
 c_1 \\
 c_2 \\
\end{array}
\right]
=
\left[
\begin{array}{ccc}
 1 & -3 \Delta x  & 9 \Delta x^2 \\
 1 & -2 \Delta x  & 4 \Delta x^2 \\
 1 & -1 \Delta x  &   \Delta x^2 \\
\end{array}
\right]
\left[
\begin{array}{c}
 c_0 \\
 c_1 \\
 c_2 \\
\end{array}
\right]
=
\left[
\begin{array}{ccc}
 1 & -3  & 9 \\
 1 & -2  & 4 \\
 1 & -1  & 1  \\
\end{array}
\right]
\left[
\begin{array}{c}
 c_0 \\
 c_1 \Delta x \\
 c_2 \Delta x^2 \\
\end{array}
\right]
=
\left[
\begin{array}{c}
 f_{-3} \\
 f_{-2} \\
 f_{-1} \\
\end{array}
\right]
= \mathbf{f}.
$$

This is known as the [Vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) for the quadratic polynomial, and its solution for $c_0$, $c_1$, and $c_2$ will be unique. But before we invert the matrix, we note that as we wish to solve for $f(x_0)=f(0)$ (again, the extrapolated value $f(x_0)$ cannot depend on the choice of the fiducial $x_0$; i.e., it will hold for *any* choice of $x_0$ so without loss of generality we will set $x_0=0$), our quadratic polynomial simplifies to:

$$
f_{0} = c_2 x_{0}^2 + c_1 x_{0} + c_0 = c_0.
$$

Thus we need only extract the value of $c_0$, which is done in the next cell.


```python
from sympy import symbols,Matrix,factor,pretty_print,simplify,latex
MM = Matrix([[1,-3,9],
             [1,-2,4],
             [1,-1,1]])
# print(latex(factor(MM.inv())))
pretty_print(factor(MM.inv()))
```

    ⎡ 1   -3   3 ⎤
    ⎢            ⎥
    ⎢3/2  -4  5/2⎥
    ⎢            ⎥
    ⎣1/2  -1  1/2⎦


Thus we get 

$$
f_0 = c_0 = f_{-3} - 3 f_{-2} + 3 f_{-1}.
$$
To determine the coefficient at the `i0=min(i0)` face, the above analysis can be repeated:


```python
from sympy import symbols,Matrix,factor,pretty_print,simplify,latex
MM = Matrix([[1,+3,9],
             [1,+2,4],
             [1,+1,1]])
# print(latex(factor(MM.inv())))
pretty_print(factor(MM.inv()))
```

    ⎡ 1    -3   3  ⎤
    ⎢              ⎥
    ⎢-3/2  4   -5/2⎥
    ⎢              ⎥
    ⎣1/2   -1  1/2 ⎦


... and we find the $c_0$ coefficient is basically the same as on the `i0=max(i0)` face; just replace $f_{-3}\to f_{+3}$, etc:

$$
f_0 = c_0 = f_{+3} - 3 f_{+2} + 3 f_{+1}.
$$

The resulting extrapolation algorithm appears in the `apply_bcs_outerextrap_and_inner()` function below. I.e.,

```c
// *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
gfs[IDX4ptS(which_gf, idx_offset0)] =
  +3.0*gfs[IDX4ptS(which_gf, idx_offset1)]
  -3.0*gfs[IDX4ptS(which_gf, idx_offset2)]
  +1.0*gfs[IDX4ptS(which_gf, idx_offset3)];
```

This function is meant to be called within a [Method of Lines timestepping algorithm](Tutorial-Method_of_Lines-C_Code_Generation.ipynb), at the very end of each MoL substep. As its name implies, it also updates the inner boundary points, which must be updated after the outer boundary points as inner points can map to outer boundary points.


```python
def add_to_Cfunction_dict_apply_bcs_outerextrap_and_inner(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    desc = r"""
Suppose the outer boundary point is at the i0=max(i0) face. Then we fit known data at i0-3, i0-2, and i0-1
  to the unique quadratic polynomial that passes through those points, and fill the data at
  i0 with the value implied from the polynomial.
As derived in nrpytutorial's Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
  the coefficients must be f_{i0} = f_{i0-3} - 3 f_{i0-2} + 3 f_{i0-1}.
  To check these coefficients are correct, consider
  * f(x0 = constant. Then f_{i0} = f_{i0-3} <- CHECK!
  * f(x) = x. WOLOG suppose x0=0. Then f_{i0} = (-3dx) - 3(-2dx) + 3(-dx) = + dx(-3+6-3) = 0 <- CHECK!
  * f(x) = x^2. WOLOG suppose x0=0. Then f_{i0} = (-3dx)^2 - 3(-2dx)^2 + 3(-dx)^2 = + dx^2(9-12+3) = 0 <- CHECK!
"""
    c_type = "void"
    name = "apply_bcs_outerextrap_and_inner"
    params = "const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx_offset0 = IDX3S(i0,i1,i2);
            const int idx_offset1 = IDX3S(i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2);
            const int idx_offset2 = IDX3S(i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2);
            const int idx_offset3 = IDX3S(i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
              gfs[IDX4ptS(which_gf, idx_offset0)] =
                +3.0*gfs[IDX4ptS(which_gf, idx_offset1)]
                -3.0*gfs[IDX4ptS(which_gf, idx_offset2)]
                +1.0*gfs[IDX4ptS(which_gf, idx_offset3)];
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(params, bcstruct, gfs);
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
```

<a id='radiation_bcs'></a>

## Step 6.b: `"radiation"` outer boundary conditions \[Back to [top](#toc)\]
$$\label{radiation_bcs}$$

<a id='radiation_bcs_theory'></a>

### Step 6.b.i: Core *ansatz*, and implications \[Back to [top](#toc)\]
$$\label{radiation_bcs_theory}$$

Here we describe outgoing radiation (a.k.a., Sommerfeld) boundary conditions.

Our implementation follows that of the [Einstein Toolkit](https://einsteintoolkit.org/)'s [`NewRad` thorn](https://einsteintoolkit.org/thornguide/EinsteinEvolve/NewRad/documentation.html), but extends to arbitrary finite-difference order.

The basic *ansatz* is that at **each outer boundary point**, an arbitrary field $f$ behaves locally as an outgoing spherical wave:

$$
f(r,t) = f_{r\to\infty} + \frac{w(r - ct)}{r} + \frac{K}{r^2},
$$
where

* $f_{r\to\infty}$ is the (constant) value the field $f$ tends to as $r\to\infty$.
* $w(r - ct)/r$ is a solution to the spherically symmetric wave equation ($\partial_t^2 (r w) - c^2 \partial_r^2 (r w) = 0$) for an *outgoing* spherical wave.
* $K/r^2$ corrects for the next-leading-order radial falloff, where $K$ is a constant determined at each boundary point $(r,\theta,\phi)$ at each time $t$ by analyzing the behavior of $f(r)$ just inside the outer boundary of the numerical domain. The full numerical implementation of this term is described below.

This boundary condition approach thus contains two free parameters: $c$ (the wavespeed for field $f$) and $f_{r\to\infty}$. Both of these should be known for any given field prior to attempting a solution.

For convenience, we apply this boundary condition not to $f$ itself but to $\partial_t f$, which will always be computed in a [Method of Lines](Tutorial-Method_of_Lines-C_Code_Generation.ipynb) approach within `NRPy+`.

Taking the first time derivative we get:

$$
\partial_t f = -c \frac{w'(r - ct)}{r}.
$$

As $w$ represents an outgoing wave, in which temporal and spatial derivatives are directly related, we will find it convenient to also compute the radial derivative $\partial_r f$ as well:

$$
\partial_r f = \frac{w'(r - ct)}{r} - \frac{w(r - ct)}{r^2} - 2\frac{K}{r^3}.
$$

Combining these two equations we get:

$$
\partial_t f = -c \left(\partial_r f + \frac{w(r - ct)}{r^2} + 2\frac{K}{r^3}\right).
$$

Next we use the *ansatz* to compute $w(r - ct)/r^2$:

\begin{align}
f(r,t) &= f_{r\to\infty} + \frac{w(r - ct)}{r} + \frac{K}{r^2} \\
\implies \frac{w(r - ct)}{r^2} &= \frac{f - f_{r\to\infty}}{r} - \frac{K}{r^3}.
\end{align}

This enables us to rewrite $\partial_t f$ as

\begin{align}
\partial_t f &= -c \left(\partial_r f + \frac{f - f_{r\to\infty}}{r} - \frac{K}{r^3} + 2\frac{K}{r^3}\right) \\
&= -c \left(\partial_r f + \frac{f - f_{r\to\infty}}{r} + \frac{K}{r^3}\right) \\
&= -\frac{c}{r} \left[r \partial_r f + \left(f - f_{r\to\infty}\right)\right] + \frac{k}{r^3},
\end{align}
where $k=-Kc$ just re-expresses the unknown function.

<a id='radiation_bcs_numerical'></a>

### Step 6.b.ii: Numerical implementation \[Back to [top](#toc)\]
$$\label{radiation_bcs_numerical}$$


We start with the equation derived in the previous section:
$$
\partial_t f = -\frac{c}{r} \left[r \partial_r f + \left(f - f_{r\to\infty}\right)\right] + \frac{k}{r^3}.
$$

First note that the right-hand side of this equation is applied *at the current time* (i.e., the same time at which $\partial_t f$ is evaluated in the Method of Lines timestepping). I.e., it is applied prior to the MoL update. Further, at the current time, $f$ is known at *all* gridpoints (interior and boundary).

As this boundary condition must be applicable to *any* curvilinear coordinate system $x^i$ (and not just spherical coordinates), $\partial_r f$ must be computed via

$$
\partial_r f = \frac{\partial x^i}{\partial r} \partial_{i} f = \frac{\partial x^i}{\partial r} \frac{\partial f}{\partial x^{i}},
$$

where $\partial_i f$ is evaluated using finite-difference derivatives, suitably upwinded to avoid extending beyond the domain of the grid.

<a id='radiation_bcs_numerical_inv_jacobian'></a>

### Step 6.b.iii: The $\partial_r f$ term: Computing $\frac{\partial x^i}{\partial r}$ \[Back to [top](#toc)\]
$$\label{radiation_bcs_numerical_inv_jacobian}$$

The term $\partial x^i/\partial r$ cannot be immediately computed, as $x^i$ is never written as a function of $r$. However, in all curvilinear coordinates supported by NRPy+, $r$, $\theta$, and $\phi$ must be specified as functions of $x^i$. Thus we have exact expressions for the Jacobian:

$$
J_i^j = \frac{\partial x^j_{\rm Sph}}{\partial x^i_{\rm Curv}},
$$

and $\partial x^i/\partial r$ can be computed from the inverse of this matrix (computed within `NRPy+`):

$$
(J^{-1})_j^i = \frac{\partial x^j_{\rm Curv}}{\partial x^i_{\rm Sph}},
$$
and only looking at the $i=0$ component (as $x^i_{\rm Sph}=r$).

For example let's set the coordinate system to Cylindrical, as exact expressions do exist for cylindrical coordinates in terms of spherical coordinates for this case:

* $\rho = r \sin\theta$
* $\phi = \phi$
* $z = r \cos\theta$

Thus

$$
\partial_r \rho = \sin\theta = \frac{\rho}{\sqrt{\rho^2 + z^2}},
$$

where the second equality can be computed directly from the inverse Jacobian. Note also that the negative root is not considered, as $\theta \in [0,\pi]$.

Again expressions for cylindrical coordinates in terms of spherical coordinates *do not exist* within `NRPy+` (as they are not generally easy to describe).

Next, let's check that NRPy+ yields the correct result for this case, choosing Cylindrical coordinates, where $x^i=$ `(xx0,xx1,xx2)` =$(\rho,\phi,z)$:


```python
par.set_parval_from_str("reference_metric::CoordSystem", "Cylindrical")
rfm.reference_metric()

# Jac_dUrfm_dDSphUD is the inverse Jacobian
Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Spherical()

print("partial rho / partial r  =  " + str(sp.simplify(Jac_dUrfm_dDSphUD[0][0])))
print("partial phi / partial r  =  " + str(sp.simplify(Jac_dUrfm_dDSphUD[1][0])))
print("partial z   / partial r  =  " + str(sp.simplify(Jac_dUrfm_dDSphUD[2][0])))
```

    partial rho / partial r  =  xx0/sqrt(xx0**2 + xx2**2)
    partial phi / partial r  =  0
    partial z   / partial r  =  xx2/sqrt(xx0**2 + xx2**2)


Next we define a NRPy+-generated inlined C function, `r_and_partial_xi_partial_r_derivs()` for computing $\frac{\partial x^i}{\partial r}$. It turns out to be convenient for this function to also compute $r$, since both are needed:


```python
def setup_Cfunction_r_and_partial_xi_partial_r_derivs():
    desc = "Compute r(xx0,xx1,xx2) and partial_r x^i."
    c_type = "static inline void"
    name = "r_and_partial_xi_partial_r_derivs"
    params = """const paramstruct *restrict params,const REAL xx0,const REAL xx1,const REAL xx2,
                                                     REAL *r,
                                                     REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    body = outputC([rfm.xxSph[0],
                    Jac_dUrfm_dDSphUD[0][0], # sp.simplify(expr) is too slow here for SinhCylindrical
                    Jac_dUrfm_dDSphUD[1][0], # sp.simplify(expr) is too slow here for SinhCylindrical
                    Jac_dUrfm_dDSphUD[2][0]],# sp.simplify(expr) is too slow here for SinhCylindrical
                   ["*r", "*partial_x0_partial_r", "*partial_x1_partial_r", "*partial_x2_partial_r"], filename="returnstring",
                   params="preindent=1,outCverbose=False,includebraces=False")
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func
```

<a id='radiation_bcs_numerical_partial_i_f'></a>

### Step 6.b.iv: The $\partial_r f$ term: Computing $\partial_i f$ with arbitrary-offset finite-difference derivatives \[Back to [top](#toc)\]
$$\label{radiation_bcs_numerical_partial_i_f}$$

Coefficients computed here can be and have been validated against the tables in this [Wikipedia article](https://en.wikipedia.org/wiki/Finite_difference_coefficient).


```python
import finite_difference as fin
# import sympy as sp

def get_arb_offset_FD_coeffs_indices(FDORDER, offset, deriv):
    # deriv = 1 <-- 1st derivative
    Minv = fin.setup_FD_matrix__return_inverse(FDORDER+1, offset)
    indices = []
    coeffs = []
    for i in range(FDORDER+1):
        indices.append(i-int(FDORDER/2) + offset)
        coeffs.append(Minv[i, deriv])
    return coeffs, indices

FDORDER=4
for offset in range(-int(FDORDER/2), int(FDORDER/2)+1):
    coeffs, indices = get_arb_offset_FD_coeffs_indices(FDORDER, offset, 1)
    print("At FD order=" + str(FDORDER) + " & stencil " + str(indices) + ", we get coeffs = " + str(coeffs))
```

    At FD order=4 & stencil [-4, -3, -2, -1, 0], we get coeffs = [1/4, -4/3, 3, -4, 25/12]
    At FD order=4 & stencil [-3, -2, -1, 0, 1], we get coeffs = [-1/12, 1/2, -3/2, 5/6, 1/4]
    At FD order=4 & stencil [-2, -1, 0, 1, 2], we get coeffs = [1/12, -2/3, 0, 2/3, -1/12]
    At FD order=4 & stencil [-1, 0, 1, 2, 3], we get coeffs = [-1/4, -5/6, 3/2, -1/2, 1/12]
    At FD order=4 & stencil [0, 1, 2, 3, 4], we get coeffs = [-25/12, 4, -3, 4/3, -1/4]


Next we define a NRPy+-generated inlined C function, `FD1_arbitrary_upwind_xN_dirn()` for computing $\partial_i f$ with an arbitrary upwind/downwind:


```python
def setup_Cfunction_FD1_arbitrary_upwind(dirn, radiation_BC_FD_order=-1):
    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if radiation_BC_FD_order == -1:
        radiation_BC_FD_order = default_FDORDER

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", radiation_BC_FD_order)

    includes = []
    desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
    c_type = "static inline REAL"
    name = "FD1_arbitrary_upwind_x"+str(dirn)+"_dirn"
    params = """const paramstruct *restrict params, const REAL *restrict gf,
                                                const int i0,const int i1,const int i2, const int offset"""
    body = r"""switch(offset) {
"""
    tmp_list = []
    for offset in range(0, int(radiation_BC_FD_order / 2) + 1):
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += "case " + str(offset) + ":\n"
        body += "  return ("
        coeffs, indices = get_arb_offset_FD_coeffs_indices(radiation_BC_FD_order, offset, 1)
        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue  # skip this iteration if coeff=0
            offset = str(indices[i])
            if i > 0:
                body += "          "
            if offset == "0":
                body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0+"+offset+",i1,i2)]\n"
                elif dirn == 1:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1+"+offset+",i2)]\n"
                elif dirn == 2:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2+"+offset+")]\n"
        body = body[:-1].replace("+-", "-") + ") * invdx"+str(dirn)+";\n"
    body += """}
return 0.0 / 0.0;  // poison output if offset computed incorrectly
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", default_FDORDER)

    return func
```

<a id='radiation_bcs_numerical_compute_partial_r_f'></a>

### Step 6.b.v: The $\partial_r f$ term: Putting it all together in `compute_partial_r_f()` \[Back to [top](#toc)\]
$$\label{radiation_bcs_numerical_compute_partial_r_f}$$

The NRPy+ generated inlined C function `compute_partial_r_f()` evaluates $\partial_r f$ via

$$
\partial_r f = \frac{\partial x^i}{\partial r}\partial_i f,
$$

where

* $\frac{\partial x^i}{\partial r}$ is computed with C function `r_and_partial_xi_partial_r_derivs()`, and
* $\partial_i f$ is computed with C function `FD1_arbitrary_upwind_xN_dirn()`;

both of these C functions are constructed above.


```python
def setup_Cfunction_compute_partial_r_f(radiation_BC_FD_order=-1):
    desc = "Compute \partial_r f"
    c_type = "static inline REAL"
    name = "compute_partial_r_f"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs,
                                       const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
                                       const int FACEi0,const int FACEi1,const int FACEi2,
                                       const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Spherical()

    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if radiation_BC_FD_order == -1:
        radiation_BC_FD_order = default_FDORDER

    body = r"""  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_FD_order/2 = """ + str(int(radiation_BC_FD_order/2)) + r"""
  const int FD1_stencil_radius = """ + str(int(radiation_BC_FD_order/2)) + r""";

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {i0,i1,i2}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
"""
    for i in range(3):
        si = str(i)
        if check_zero(Jac_dUrfm_dDSphUD[i][0]):
            body += "  const REAL partial_x"+si+"_f=0.0;\n"
        else:
            body += "  int i"+si+"_offset = FACEi"+si+";  // Shift stencil away from the face we're updating.\n"
            body += "  // Next adjust i"+si+"_offset so that FD stencil never goes out of bounds.\n"
            body += "  if(dest_i"+si+" < FD1_stencil_radius) i"+si+"_offset = FD1_stencil_radius-dest_i"+si+";\n"
            body += "  else if(dest_i"+si+" > (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1)) i"+si+"_offset = (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1) - dest_i"+si+";\n"
            body += "  const REAL partial_x"+si+"_f=FD1_arbitrary_upwind_x"+si+"_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i"+si+"_offset);\n"
    body += "  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;\n"
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func
```

<a id='radiation_bcs_evaluating_k'></a>

### Step 6.b.vi: Evaluating the $k$ in the $k/r^3$ term \[Back to [top](#toc)\]
$$\label{radiation_bcs_evaluating_k}$$

First we note that in the *ansatz*, if $f(r,t)$ perfectly satisfied the outgoing wave equation, then $k=0$ and

$$
\left[\partial_t f\right]_{\rm Outgoing\ Wave} = -\frac{c}{r} \left[r \partial_r f + \left(f - f_{r\to\infty}\right)\right].
$$

However, numerical error, errors associated with outer boundaries being placed at finite (as opposed to infinite) radii, and higher-order radial falloffs will invariably lead to $k$ being *nonzero*. 

Given that $\partial_t f$ is evaluated *at all points in the interior* of the grid, the difference between perfect satisfaction of the radially outgoing wave equation and the actual solution is *known* at these points. Let's call this difference $\xi$:

$$
\xi = \partial_t f - \left[\partial_t f\right]_{\rm Outgoing\ Wave} \equiv \frac{k}{r^3}.
$$

We compute $\xi_{\rm int} = \xi(r_{\rm int})$ at a neighboring interior point at radius $r=r_{\rm int}$, so that

$$
\xi_{\rm int} = [\partial_t f]_{\rm int} - \left[\partial_t f\right]_{\rm Outgoing\ Wave,\ int} \equiv \frac{k}{r_{\rm int}^3}.
$$

In this way, we obtain $k$:

$$
k = r_{\rm int}^3 \left([\partial_t f]_{\rm int} - \left[\partial_t f\right]_{\rm Outgoing\ Wave,\ int}\right)
$$

To determine the appropriate interior point, we simply keep track of the current boundary face being updated and choose the nearest neighbor in the direction of the grid center (assumed to be at $r=0$).

As computing $k$ depends on all other functions, we incorporate its computation into the core `apply_bcs_outerradiation_and_inner()` routine - the subject of the next subsection.

<a id='radiation_bcs_apply_bcs_radiation'></a>

### Step 6.b.vii: Putting it all together: `radiation_bcs_single_pt()` \[Back to [top](#toc)\]
$$\label{radiation_bcs_apply_bcs_radiation}$$

`radiation_bcs_single_pt()` combines all the above algorithms to evaluate radiation boundary conditions at a single point.

Recall that at each point on the boundary, we have
$$
\partial_t f = \left[\partial_t f\right]_{\rm Outgoing\ Wave} + \frac{k}{r^3},
$$
where
$$
\left[\partial_t f\right]_{\rm Outgoing\ Wave} = -\frac{c}{r} \left[r \partial_r f + \left(f - f_{r\to\infty}\right)\right],
$$
and
$$
k = r_{\rm int}^3 \left([\partial_t f]_{\rm int} - \left[\partial_t f\right]_{\rm Outgoing\ Wave,\ int}\right).
$$


```python
def setup_Cfunction_radiation_bcs(radiation_BC_FD_order=-1):
    includes = []
    prefunc = ""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    for i in range(3):
        # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
        if not check_zero(Jac_dUrfm_dDSphUD[i][0]):
            prefunc += setup_Cfunction_FD1_arbitrary_upwind(dirn=i, radiation_BC_FD_order=radiation_BC_FD_order)
    prefunc += setup_Cfunction_r_and_partial_xi_partial_r_derivs()
    prefunc += setup_Cfunction_compute_partial_r_f(radiation_BC_FD_order=radiation_BC_FD_order)
    desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
    c_type = "static inline REAL"
    name = "radiation_bcs"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct,REAL *restrict xx[3],
                                 const REAL *restrict gfs, REAL *restrict gfs_rhss,
                                 const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
                                 const int dest_i0,const int dest_i1,const int dest_i2,
                                 const short FACEi0,const short FACEi1,const short FACEi2"""
    body = r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
REAL r, partial_x0_partial_r,partial_x1_partial_r,partial_x2_partial_r;
REAL r_int, partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int;
r_and_partial_xi_partial_r_derivs(params,xx[0][dest_i0],xx[1][dest_i1],xx[2][dest_i2],
                                  &r, &partial_x0_partial_r, &partial_x1_partial_r,  &partial_x2_partial_r);
r_and_partial_xi_partial_r_derivs(params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int],
                                  &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
const REAL partial_r_f     = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
const REAL partial_r_f_int = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);

const int idx3 = IDX3S(dest_i0,dest_i1,dest_i2);
const int idx3_int = IDX3S(dest_i0_int,dest_i1_int,dest_i2_int);

const REAL partial_t_f_int = gfs_rhss[IDX4ptS(which_gf, idx3_int)];

const REAL c = gf_wavespeed;
const REAL f_infinity = gf_f_infinity;
const REAL f     = gfs[IDX4ptS(which_gf, idx3)];
const REAL f_int = gfs[IDX4ptS(which_gf, idx3_int)];
const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

const REAL rinv = 1.0 / r;
const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

return partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=includes, prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func
```

<a id='apply_bcs_outerradiation_and_inner'></a>

### Step 6.b.viii: `apply_bcs_outerradiation_and_inner()`: Apply radiation BCs at all outer boundary points and inner BCs at all inner boundary points \[Back to [top](#toc)\]
$$\label{apply_bcs_outerradiation_and_inner}$$


```python
def add_to_Cfunction_dict_apply_bcs_outerradiation_and_inner(rel_path_to_Cparams=os.path.join("."), radiation_BC_FD_order=2):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    prefunc = setup_Cfunction_radiation_bcs(radiation_BC_FD_order=radiation_BC_FD_order)
    desc = ""
    c_type = "void"
    name = "apply_bcs_outerradiation_and_inner"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                        const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS],
                                        REAL *restrict gfs, REAL *restrict rhs_gfs"""
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx3 = IDX3S(i0,i1,i2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply radiation BCs to all outer boundary points. ***
              rhs_gfs[IDX4ptS(which_gf, idx3)] = radiation_bcs(params, bcstruct, xx, gfs, rhs_gfs, which_gf,
                                                               custom_wavespeed[which_gf], custom_f_infinity[which_gf],
                                                               i0,i1,i2, FACEX0,FACEX1,FACEX2);
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
```

<a id='start2finish'></a>

# Step 7: `CurviBC_Playground.c`: Start-to-Finish C code module for testing & validating curvilinear boundary conditions \[Back to [top](#toc)\]
$$\label{start2finish}$$

<a id='register_gfs'></a>

## Step 7.a: Register gridfunctions of all 10 parity types \[Back to [top](#toc)\]
$$\label{register_gfs}$$ 

Here we register within NRPy+ one gridfunction per each of the 10 parity conditions. These will be output automatically to `NRPy_basic_defines.h` below.


```python
# Step 7.a: Register gridfunctions of all 10 parity types

# 6 gridfunctions, corresponding to all unique rank-2 tensor components:
ranktwosymmDD = ixp.register_gridfunctions_for_single_rank2("EVOL","ranktwosymmDD", "sym01", f_infinity=0.0, wavespeed=1.0)
# 3 gridfunctions, corresponding to all unique rank-1 tensor components:
rankoneU = ixp.register_gridfunctions_for_single_rank1("EVOL","rankoneU", f_infinity=0.0, wavespeed=1.0)
# 1 rank-0 (scalar) gridfunction
# WAS: rankzero = ixp.gri.register_gridfunctions("EVOL","rankzero", f_infinity=1.0, wavespeed=sp.sqrt(2.0))
# NOW: consistent with ScalarWaveCurvilinear:
rankzero = ixp.gri.register_gridfunctions("EVOL","rankzero", f_infinity=1.0, wavespeed=1.0)
```

<a id='validate'></a>

## Step 7.b: Set up test data for Curvilinear Boundary Conditions code validation \[Back to [top](#toc)\]
$$\label{validate}$$

We will validate this curvilinear boundary condition module by comparing its results with the original (trusted) SENR code, as follows:

* **Discrete data test**:
    1. Fill all 10 gridfunctions at each gridpoint with the unique gridpoint integer index `IDX3S(i0,i1,i2)`
    1. Apply curvilinear boundary conditions
    1. Compare output data at all gridpoints with those from the original SENR code. Agreement should be perfect.

Another (future, to-be-implemented) test, which will enable us to validate coordinate systems that do not exist within the original SENR code, is described below:

* **Smooth data test** (TODO):
    1. Fill all 10 gridfunctions with data that are smooth in the Cartesian basis.
    1. Apply Jacobian transformation to all data points, to convert to curvilinear basis
    1. Apply curvilinear boundary conditions
    1. Apply Jacobian transformation to all data points, to convert back to Cartesian basis
    1. Compute difference between original Cartesian data and transformed data. Difference should be zero (to within roundoff) at all points except those that are influenced by outer boundary conditions.


```python
def add_to_Cfunction_dict_CurviBC_Discrete_initial_data():
    add_to_Cfunction_dict(
        includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"],
        desc     ="Populate in_gfs[] with discrete, regular data",
        c_type     ="void",
        name     ="CurviBC_Discrete_initial_data",
        params   ="const paramstruct *restrict params,REAL *restrict in_gfs",
        body     ="""
#pragma omp parallel for
    for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i++) {
        in_gfs[i] = (REAL)i;
    }
""", rel_path_to_Cparams = os.path.join("."))
```


```python
import loop as lp
import ScalarWave.InitialData as swid
def add_to_Cfunction_dict_CurviBC_Smooth_ScalarWave_initial_data():
    swid.SphericalGaussian(CoordSystem=CoordSystem)
    add_to_Cfunction_dict(
        includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"],
        desc     ="Populate in_gfs[] and in_gfs_rhss[] with smooth, regular data",
        c_type     ="void",
        name     ="CurviBC_Smooth_ScalarWave_initial_data",
        params   ="""const paramstruct *restrict params,REAL *restrict xx[3],
                                            REAL *restrict in_gfs,REAL *restrict in_gfs_rhss""",
        body    = """
  const REAL sigma = 3.0;
  const REAL time = 0.0;
""" + lp.simple_loop('AllPoints, Read_xxs',
                                  outputC([swid.uu_ID,swid.vv_ID],
              ["in_gfs[IDX4S(RANKZEROGF,i0,i1,i2)]","in_gfs_rhss[IDX4S(RANKZEROGF,i0,i1,i2)]"],
               filename="returnstring", params="outCverbose=False,includebraces=False")),
        rel_path_to_Cparams = os.path.join("."))
```

<a id='mainc'></a>

## Step 7.c: `CurviBC_Playground`'s `main.c` Code \[Back to [top](#toc)\]
$$\label{mainc}$$


```python
def add_to_Cfunction_dict_main__CurviBC_Playground():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Write test data to gridfunctions
// Step 2: Overwrite all data in ghost zones with NaNs
// Step 3: Apply curvilinear boundary conditions
// Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
// Step 5: Free all allocated memory
"""
    c_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"
    body  = r"""  // Step 0a: Read command-line input, error out if nonconformant
  if(argc != 2) {
    fprintf(stderr,"Error: Expected one command-line argument: ./CurviBC_Playground [test type: Smooth or Discrete],\n");
    exit(1);
  }

  griddata_struct griddata;
  set_Cparameters_to_default(&griddata.params);

  char CoordSystem_name[50];
  snprintf(CoordSystem_name, 50, """
    body += "\""+CoordSystem+"\""
    body += r""");

  // Step 0b: Set number of gridpoints...
  const int Nxx[3] = { 4, 4, 4 };

  // Step 0c: Set test type to Smooth or Discrete
  char test_type[100];
  snprintf(test_type, 100, "%s", argv[1]);
  if(strncmp("Smooth",  test_type, 100) != 0 &&
     strncmp("Discrete",test_type, 100) != 0) {
    fprintf(stderr,"Error: test type = %s not supported. Choose Smooth or Discrete (CASE SENSITIVE).\n",test_type);
    exit(1);
  }

  // Step 0d: Uniform coordinate grids are stored to *xx[3]
  // Step 0d.i: Set bcstruct
  {
    int EigenCoord;
    EigenCoord = 1;
    // Step 0d.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0e: Find ghostzone mappings; set up bcstruct
    bcstruct_set_up(&griddata.params, griddata.xx, &griddata.bcstruct);
    // Step 0e.i: Free allocated space for xx[][] array
    for(int i=0;i<3;i++) free(griddata.xx[i]);

    // Step 0f: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //          params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //          chosen (non-Eigen) CoordSystem.
    EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
  }

  // Step 0g: Set all C parameters "blah" for params.blah, including
  //          Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0, etc.
  const int Nxx_plus_2NGHOSTS0 = griddata.params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata.params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata.params.Nxx_plus_2NGHOSTS2;

  // Step 0h: Allocate memory for gridfunctions and store to local pointers test_gfs & test_gfs_rhss
  griddata.test_gfs      = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2);
  griddata.test_gfs_rhss = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2);
  REAL *restrict test_gfs      = griddata.test_gfs;
  REAL *restrict test_gfs_rhss = griddata.test_gfs_rhss;

  // Step 1: Write test data to gridfunctions
  if(strncmp("Discrete", test_type, 100)==0) {
    CurviBC_Discrete_initial_data(&griddata.params, test_gfs);
    CurviBC_Discrete_initial_data(&griddata.params, test_gfs_rhss);
  } else {
    fprintf(stderr, "Sorry, curvilinear boundary conditions test = %s not yet supported. Feel free to contribute!\n",test_type);
  }

  // Step 2: Overwrite all data in ghost zones with NaNs
  LOOP_REGION(0,Nxx_plus_2NGHOSTS0, 0,Nxx_plus_2NGHOSTS1, 0,Nxx_plus_2NGHOSTS2) {
    for(int gf=0;gf<NUM_EVOL_GFS;gf++) {
      const int idx4 = IDX4S(gf,i0,i1,i2);
      if(i0 < NGHOSTS || i0 >= Nxx_plus_2NGHOSTS0-NGHOSTS) test_gfs_rhss[idx4] = +(0.0 / 0.0);
      if(i1 < NGHOSTS || i1 >= Nxx_plus_2NGHOSTS1-NGHOSTS) test_gfs_rhss[idx4] = +(0.0 / 0.0);
      if(i2 < NGHOSTS || i2 >= Nxx_plus_2NGHOSTS2-NGHOSTS) test_gfs_rhss[idx4] = +(0.0 / 0.0);
    }
  }

  // Step 3: Apply outer and inner boundary conditions
  if(griddata.params.outer_bc_type == EXTRAPOLATION_OUTER_BCS) {
    // Warning: Extrapolation BCs should generally be applied to test_gfs, not test_gfs_rhss.
    //          We apply them to rhss here just to make the coding simpler (gf_rhss have the
    //          same data in the Discrete test as gfs).
    apply_bcs_outerextrap_and_inner(&griddata.params, &griddata.bcstruct, test_gfs_rhss);
  } else if(griddata.params.outer_bc_type == RADIATION_OUTER_BCS) {
    apply_bcs_outerradiation_and_inner(&griddata.params, &griddata.bcstruct, griddata.xx,
                                       gridfunctions_wavespeed,gridfunctions_f_infinity,
                                       test_gfs,test_gfs_rhss);
  }

  // Step 4: Print gridfunction data after curvilinear boundary conditions have been applied:
  char filename[120];  sprintf(filename,"out4x4x4-%s-NGHOSTS4oFD.txt", CoordSystem_name);
  FILE *outfile = fopen(filename, "w");

  LOOP_REGION(0,Nxx_plus_2NGHOSTS0, 0,Nxx_plus_2NGHOSTS1, 0,Nxx_plus_2NGHOSTS2) {
    fprintf(outfile, "%d %d %d | ", i0,i1,i2);
    for(int gf=0;gf<NUM_EVOL_GFS;gf++) {
      const int idx4 = IDX4S(gf,i0,i1,i2);
      if(!isnan(test_gfs_rhss[idx4])) {
        fprintf(outfile, "%d ", (int)test_gfs_rhss[idx4]);
      } else {
        fprintf(stderr, "ERROR: found NaN %d %d %d %d %d\n", gf, i0,i1,i2, NUM_EVOL_GFS);
        //exit(1);
      }
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);

  // Step 5: Free all allocated memory
  free(griddata.bcstruct.inner_bc_array);
  for(int ng=0;ng<NGHOSTS*3;ng++) free(griddata.bcstruct.pure_outer_bc_array[ng]);
  free(griddata.test_gfs);  free(griddata.test_gfs_rhss);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body, enableCparameters=False)
```

<a id='curvibc_setupall'></a>

## Step 7.d: `CurviBoundaryConditions_register_C_functions_and_NRPy_basic_defines()`: Add all CurviBC C codes to C function dictionary, and add CurviBC definitions to `NRPy_basic_defines.h` \[Back to [top](#toc)\]
$$\label{curvibc_setupall}$$

Populate `NRPy_basic_defines.h` with `bcstruct` and friends, as well as parity types for registered gridfunctions.


```python
# Only call this after ALL gridfunctions have been registered!
def CurviBoundaryConditions_register_NRPy_basic_defines(verbose=True):
    # Then set up the dictionary entry for CurviBC in NRPy_basic_defines
    Nbd_str  = NRPy_basic_defines_CurviBC_data_structures()
    Nbd_str += NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=verbose)
    outC_NRPy_basic_defines_h_dict["CurviBoundaryConditions"] = Nbd_str

    # Register griddata_struct variables for this module,
    #   where griddata_struct is declared in NRPy_basic_defines.h
    gri.glb_griddata_struct_list += [gri.glb_griddata(__name__, "bc_struct bcstruct;")]
```


```python
def CurviBoundaryConditions_register_C_functions(rel_path_to_Cparams=os.path.join("./"),
                                                 radiation_BC_FD_order=4):
    add_to_Cfunction_dict_bcstruct_set_up(rel_path_to_Cparams=rel_path_to_Cparams)
    add_to_Cfunction_dict_apply_bcs_outerradiation_and_inner(rel_path_to_Cparams=rel_path_to_Cparams,
                                                             radiation_BC_FD_order=radiation_BC_FD_order)
    add_to_Cfunction_dict_apply_bcs_inner_only(rel_path_to_Cparams=rel_path_to_Cparams)
    add_to_Cfunction_dict_apply_bcs_outerextrap_and_inner(rel_path_to_Cparams=rel_path_to_Cparams)
```

<a id='add_cfunction_dicts'></a>

# Step 8: Add all C functions to `Cfunction_dict`, also output `NRPy_basic_defines.h` and `NRPy_function_prototypes.h` \[Back to [top](#toc)\]
$$\label{add_cfunction_dicts}$$


```python
import outputC as outC
# Then we set the coordinate system for the numerical grid
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()
# Step 5: Generate & register set_Nxx_dxx_invdx_params__and__xx(), which sets
#         params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
#         chosen (not necessarily Eigen-) CoordSystem.
rfm.register_NRPy_basic_defines(enable_rfm_precompute=False)
rfm.add_to_Cfunc_dict_set_Nxx_dxx_invdx_params__and__xx()
rfm.register_NRPy_basic_defines()

CurviBoundaryConditions_register_NRPy_basic_defines()
CurviBoundaryConditions_register_C_functions()

add_to_Cfunction_dict_apply_bcs_inner_only()
add_to_Cfunction_dict_apply_bcs_outerextrap_and_inner()
add_to_Cfunction_dict_apply_bcs_outerradiation_and_inner(radiation_BC_FD_order=RADIATION_BC_FD_ORDER)


outC.outputC_register_C_functions_and_NRPy_basic_defines()  # #define M_PI,  etc.
# Declare paramstruct, register set_Cparameters_to_default(),
#   and output declare_Cparameters_struct.h and set_Cparameters[].h:
outC.NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(os.path.join(Ccodesrootdir))
par.register_NRPy_basic_defines()  # add `paramstruct params` to griddata struct.

gri.register_C_functions_and_NRPy_basic_defines(list_of_extras_in_griddata_struct=
                                                ["REAL *restrict test_gfs","REAL *restrict test_gfs_rhss"])  # #define IDX3S(),  etc.
fin.register_C_functions_and_NRPy_basic_defines(NGHOSTS_account_for_onezone_upwind=True)  # #define NGHOSTS, etc.
# Comment out above line and uncomment below line to confirm independent Python module agrees.
# import CurviBoundaryConditions.CurviBoundaryConditions as CBC
# CBC.CurviBoundaryConditions_register_C_functions_and_NRPy_basic_defines()

# initial data function:
add_to_Cfunction_dict_CurviBC_Discrete_initial_data()
add_to_Cfunction_dict_CurviBC_Smooth_ScalarWave_initial_data()

# main function:
add_to_Cfunction_dict_main__CurviBC_Playground()

outC.construct_NRPy_basic_defines_h(Ccodesrootdir)
outC.construct_NRPy_function_prototypes_h(Ccodesrootdir)
# The following is run from inside cmd.new_C_compile() in the next code cell.
# outC.construct_Makefile_from_outC_function_dict(Ccodesrootdir, "CurviBC_Playground")
# print(outC.outC_function_dict["set_Nxx_dxx_invdx_params__and__xx"])
```

    Evolved gridfunction "rankoneU0" has parity type 1.
    Evolved gridfunction "rankoneU1" has parity type 2.
    Evolved gridfunction "rankoneU2" has parity type 3.
    Evolved gridfunction "ranktwosymmDD00" has parity type 4.
    Evolved gridfunction "ranktwosymmDD01" has parity type 5.
    Evolved gridfunction "ranktwosymmDD02" has parity type 6.
    Evolved gridfunction "ranktwosymmDD11" has parity type 7.
    Evolved gridfunction "ranktwosymmDD12" has parity type 8.
    Evolved gridfunction "ranktwosymmDD22" has parity type 9.
    Evolved gridfunction "rankzero" has parity type 0.


<a id='senr_compare'></a>

# Step 9: Validation: Compile & compare with original (trusted) SENR results \[Back to [top](#toc)\]
$$\label{senr_compare}$$


```python
import cmdline_helper as cmd
# from outputC import construct_Makefile_from_outC_function_dict
# construct_Makefile_from_outC_function_dict(Ccodesrootdir, "CurviBC_Playground", compiler_opt_option="fast")
cmd.new_C_compile(Ccodesrootdir, "CurviBC_Playground", compiler_opt_option="fast") # fastdebug or debug also supported
os.chdir(Ccodesrootdir)
cmd.Execute("CurviBC_Playground", "Discrete")
#cmd.Execute("CurviBC_Playground", "4 4 4 Discrete", "out4x4x4-Spherical-NGHOSTS4oFD.txt")
os.chdir("..")
```

    (EXEC): Executing `make -j18`...
    (BENCH): Finished executing in 0.41 seconds.
    Finished compilation.
    (EXEC): Executing `taskset -c 1,3,5,7,9,11,13,15 ./CurviBC_Playground Discrete`...
    (BENCH): Finished executing in 0.20 seconds.



```python
import filecmp
if outer_bcs_type == "extrapolation":
    if "Cylindrical" in CoordSystem:
        if filecmp.cmp(os.path.join(Ccodesrootdir, "out4x4x4-"+CoordSystem+"-NGHOSTS4oFD.txt"),
                       os.path.join("CurviBoundaryConditions", "SENRout4x4x4-Cylindrical_NGHOSTS4oFD.txt")) == False:
            print("ERROR: "+CoordSystem+" boundary conditions malfunction!")
            sys.exit(1)
    elif "Spherical" in CoordSystem:
        if filecmp.cmp(os.path.join(Ccodesrootdir, "out4x4x4-"+CoordSystem+"-NGHOSTS4oFD.txt"),
                       os.path.join("CurviBoundaryConditions", "SENRout4x4x4-Spherical_NGHOSTS4oFD.txt")) == False:
            print("ERROR: "+CoordSystem+" boundary conditions malfunction!")
            sys.exit(1)
    else:
        print("ERROR: "+CoordSystem+" coordinate system comparison unavailable!")
        sys.exit(1)
    print(CoordSystem + " boundary condition comparison test between this tutorial notebook & trusted original SENR code: PASSED")
```

    SinhSpherical boundary condition comparison test between this tutorial notebook & trusted original SENR code: PASSED


<a id='latex_pdf_output'></a>

# Step 10: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Start_to_Finish-Curvilinear_BCs.pdf](Tutorial-Start_to_Finish-Curvilinear_BCs.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Start_to_Finish-Curvilinear_BCs")
```

    Created Tutorial-Start_to_Finish-Curvilinear_BCs.tex, and compiled LaTeX
        file to PDF file Tutorial-Start_to_Finish-Curvilinear_BCs.pdf

