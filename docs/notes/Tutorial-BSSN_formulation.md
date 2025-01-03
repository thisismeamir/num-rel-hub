<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# The BSSN Formulation of General Relativity in Generic Curvilinear Coordinates: An Overview

## Author: Zach Etienne

## This notebook outlines the BSSN formulation of Einstein's equations in generic curvilinear coordinates via NRPy+. It discusses the derivation, numerical implementation, and time evolution of BSSN variables and gauge quantities. It also explores the application of Hamiltonian and momentum constraints in constructing spacetime from initial data.

\*As Einstein's equations in this formalism take the form of highly nonlinear, coupled *wave equations*, the [tutorial notebook on the scalar wave equation in curvilinear coordinates](Tutorial-ScalarWaveCurvilinear.ipynb) is *required* reading before beginning this module. That module, as well as its own prerequisite [module on reference metrics within NRPy+](Tutorial-Reference_Metric.ipynb) provides the needed overview of how NRPy+ handles reference metrics.

## Introduction:
NRPy+'s original purpose was to be an easy-to-use code capable of generating Einstein's equations in a broad class of [singular](https://en.wikipedia.org/wiki/Coordinate_singularity), curvilinear coordinate systems, where the user need only input the scale factors of the underlying reference metric. Upon generating these equations, NRPy+ would then leverage SymPy's [common-expression-elimination (CSE)](https://en.wikipedia.org/wiki/Common_subexpression_elimination) and C code generation routines, coupled to its own [single-instruction, multiple-data (SIMD)](https://en.wikipedia.org/wiki/SIMD) functions, to generate highly-optimized C code.

### Background Reading/Lectures:

* Mathematical foundations of BSSN and 3+1 initial value problem decompositions of Einstein's equations:
    * [Thomas Baumgarte's lectures on mathematical formulation of numerical relativity](https://www.youtube.com/watch?v=t3uo2R-yu4o&list=PLRVOWML3TL_djTd_nsTlq5aJjJET42Qke)
    * [Yuichiro Sekiguchi's introduction to BSSN](http://www2.yukawa.kyoto-u.ac.jp/~yuichiro.sekiguchi/3+1.pdf) 
* Extensions to the standard BSSN approach used in NRPy+
    * [Brown's covariant "Lagrangian" formalism of BSSN](https://arxiv.org/abs/0902.3652)
    * [BSSN in spherical coordinates, using the reference-metric approach of Baumgarte, Montero, Cordero-Carrión, and Müller (2012)](https://arxiv.org/abs/1211.6632)
    * [BSSN in generic curvilinear coordinates, using the extended reference-metric approach of Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658)

### A Note on Notation:

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This module lays out the mathematical foundation for the BSSN formulation of Einstein's equations, as detailed in the references in the above Background Reading/Lectures section. It is meant to provide an overview of the basic equations and point of reference for **full tutorial notebooks** linked below:

1. [Step 1](#brownslagrangebssn): [Brown](https://arxiv.org/abs/0902.3652)'s covariant formulation of the BSSN time-evolution equations (see next section for gauge conditions)
    1. [Step 1.a](#fullequations): Numerical implementation of BSSN time-evolution equations
        1. [Step 1.a.i](#liederivs) ([**BSSN quantities module [start here]**](Tutorial-BSSN_quantities.ipynb); [**BSSN time-evolution module**](Tutorial-BSSN_time_evolution-BSSN_RHSs.ipynb)): Expanding the Lie derivatives; the BSSN time-evolution equations in their final form
1. [Step 2](#gaugeconditions) ([**full tutorial notebook**](Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb)): Time-evolution equations for the BSSN |gauge quantities $\alpha$ and $\beta^i$
1. [Step 3](#constraintequations) ([**full tutorial notebook**](Tutorial-BSSN_constraints.ipynb)): The BSSN constraint equations
    1. [Step 3.a](#hamiltonianconstraint): The Hamiltonian constraint
    1. [Step 3.b](#momentumconstraint): The momentum constraint
1. [Step 4](#gammaconstraint) ([**full tutorial notebook**](Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.ipynb)): The BSSN algebraic constraint $\hat{\gamma}=\bar{\gamma}$
1. [Step 5](#latex_pdf_output) Output this notebook to $\LaTeX$-formatted PDF file

<a id='brownslagrangebssn'></a>

# Step 1: [Brown](https://arxiv.org/abs/0902.3652)'s covariant formulation of BSSN \[Back to [top](#toc)\]
$$\label{brownslagrangebssn}$$

The covariant "Lagrangian" BSSN formulation of [Brown (2009)](https://arxiv.org/abs/0902.3652), which requires

$$
\partial_t \bar{\gamma} = 0,
$$

results in the BSSN equations taking the following form (Eqs. 11 and 12 in [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658)):

\begin{align}
  \partial_{\perp} \bar{\gamma}_{i j} {} = {} & \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right ) - 2 \alpha \bar{A}_{i j} \; , \\
  \partial_{\perp} \bar{A}_{i j} {} = {} & -\frac{2}{3} \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K \nonumber \\
  & + e^{-4 \phi} \left \{-2 \alpha \bar{D}_{i} \bar{D}_{j} \phi + 4 \alpha \bar{D}_{i} \phi \bar{D}_{j} \phi \right . \nonumber \\
    & \left . + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_{i} \bar{D}_{j} \alpha + \alpha \bar{R}_{i j} \right \}^{\text{TF}} \; , \\
  \partial_{\perp} \phi {} = {} & \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) \; , \\
  \partial_{\perp} K {} = {} & \frac{1}{3} \alpha K^{2} + \alpha \bar{A}_{i j} \bar{A}^{i j} \nonumber \\
  & - e^{-4 \phi} \left (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi \right ) \; , \\
  \partial_{\perp} \bar{\Lambda}^{i} {} = {} & \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} \nonumber \\
  & - 2 \bar{A}^{i j} \left (\partial_{j} \alpha - 6 \partial_{j} \phi \right ) + 2 \bar{A}^{j k} \Delta_{j k}^{i} \nonumber \\
  & -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K \\
\end{align}
where 
* the $\text{TF}$ superscript denotes the trace-free part.
* $\bar{\gamma}_{ij} = \varepsilon_{i j} + \hat{\gamma}_{ij}$, where $\bar{\gamma}_{ij} = e^{-4\phi} \gamma_{ij}$ is the conformal metric, $\gamma_{ij}$ is the physical metric (see below), and $\varepsilon_{i j}$ encodes information about the non-hatted metric.
* $\gamma_{ij}$, $\beta^i$, and $\alpha$ are the physical (as opposed to conformal) spatial 3-metric, shift vector, and lapse, respectively, which may be defined via the 3+1 decomposition line element (in [$G=c=1$ units](https://en.wikipedia.org/wiki/Planck_units)):
$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij}\left(dx^i + \beta^i dt\right)\left(dx^j + \beta^j dt\right).$$
* $\bar{R}_{ij}$ is the conformal Ricci tensor, computed via
\begin{align}
  \bar{R}_{i j} {} = {} & - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} + \bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k} + \Delta^{k} \Delta_{(i j) k} \nonumber \\
  & + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right ) \; .
\end{align}
* $\partial_{\perp} = \partial_t - \mathcal{L}_\beta$; $\mathcal{L}_\beta$ is the [Lie derivative](https://en.wikipedia.org/wiki/Lie_derivative) along the shift vector $\beta^i$.
* $\partial_0 = \partial_t - \beta^i \partial_i$ is an advective time derivative.
* $\hat{D}_j$ is the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) with respect to the reference metric $\hat{\gamma}_{ij}$.
* $\bar{D}_j$ is the [covariant derivative](https://en.wikipedia.org/wiki/Covariant_derivative) with respect to the barred spatial 3-metric $\bar{\gamma}_{ij}$
* $\Delta^i_{jk}$ is the tensor constructed from the difference of barred and hatted Christoffel symbols:
$$\Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}$$
    * The related quantity $\Delta^i$ is defined $\Delta^i \equiv \bar{\gamma}^{jk} \Delta^i_{jk}$.
* $\bar{A}_{ij}$ is the conformal, trace-free extrinsic curvature: 
$$\bar{A}_{ij} = e^{-4\phi} \left(K_{ij} - \frac{1}{3}\gamma_{ij} K\right),$$
where $K$ is the trace of the extrinsic curvature $K_{ij}$.

<a id='fullequations'></a>

## Step 1.a: Numerical implementation of BSSN time-evolution equations \[Back to [top](#toc)\]
$$\label{fullequations}$$

Regarding the numerical implementation of the above equations, first, notice the left-hand sides of the equations include the time derivatives. Numerically, these equations are solved using an [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem), where data are specified at a given time $t$, so that the solution at any later time can be obtained using the [Method of Lines (MoL)](https://en.wikipedia.org/wiki/Method_of_lines). MoL requires that the equations be written in the form:

$$\partial_t \vec{U} = \vec{f}\left(\vec{U},\partial_i \vec{U}, \partial_i \partial_j \vec{U},...\right),$$

for the vector of "evolved quantities" $\vec{U}$, where the right-hand side vector $\vec{f}$ *does not* contain *explicit* time derivatives of $\vec{U}$.

Thus we must first rewrite the above equations so that *only* partial derivatives of time appear on the left-hand sides of the equations, meaning that the Lie derivative terms must be moved to the right-hand sides of the equations.

<a id='liederivs'></a>

### Step 1.a.i: Expanding the Lie derivatives; BSSN equations in their final form \[Back to [top](#toc)\]
$$\label{liederivs}$$

In this Step, we provide explicit expressions for the [Lie derivatives](https://en.wikipedia.org/wiki/Lie_derivative) $\mathcal{L}_\beta$ appearing inside the $\partial_\perp = \partial_t - \mathcal{L}_\beta$ operators for $\left\{\bar{\gamma}_{i j},\bar{A}_{i j},W, K, \bar{\Lambda}^{i}\right\}$.

In short, the Lie derivative of tensor weight $w$ is given by (from [the wikipedia article on Lie derivatives](https://en.wikipedia.org/wiki/Lie_derivative))
\begin{align}
(\mathcal {L}_X T) ^{a_1 \ldots a_r}{}_{b_1 \ldots b_s} &= X^c(\partial_c T^{a_1 \ldots a_r}{}_{b_1 \ldots b_s}) \\
&\quad - (\partial_c X ^{a_1}) T ^{c a_2 \ldots a_r}{}_{b_1 \ldots b_s} - \ldots - (\partial_c X^{a_r}) T ^{a_1 \ldots a_{r-1}c}{}_{b_1 \ldots b_s} \\
 &\quad  +  (\partial_{b_1} X^c) T ^{a_1 \ldots a_r}{}_{c b_2 \ldots b_s} + \ldots + (\partial_{b_s} X^c) T ^{a_1 \ldots a_r}{}_{b_1 \ldots b_{s-1} c} + w (\partial_{c} X^c) T ^{a_1 \ldots a_r}{}_{b_1 \ldots b_{s}}
\end{align}

Thus to evaluate the Lie derivative, one must first know the tensor density weight $w$ for each tensor. In this formulation of Einstein's equations, **all evolved quantities have density weight $w=0$**, so according to the definition of Lie derivative above,
\begin{align}
\mathcal{L}_\beta \bar{\gamma}_{ij} &= \beta^k \partial_k \bar{\gamma}_{ij} + \partial_i \beta^k \bar{\gamma}_{kj} + \partial_j \beta^k \bar{\gamma}_{ik}, \\
\mathcal{L}_\beta \bar{A}_{ij} &= \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}, \\
\mathcal{L}_\beta \phi &= \beta^k \partial_k \phi, \\
\mathcal{L}_\beta K &= \beta^k \partial_k K, \\
\mathcal{L}_\beta \bar{\Lambda}^i &= \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
\end{align}

With these definitions, the BSSN equations for the un-rescaled evolved variables in the form $\partial_t \vec{U} = f\left(\vec{U},\partial_i \vec{U}, \partial_i \partial_j \vec{U},...\right)$ become

\begin{align}
  \partial_t \bar{\gamma}_{i j} {} = {} & \left[\beta^k \partial_k \bar{\gamma}_{ij} + \partial_i \beta^k \bar{\gamma}_{kj} + \partial_j \beta^k \bar{\gamma}_{ik} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right ) - 2 \alpha \bar{A}_{i j} \; , \\
  \partial_t \bar{A}_{i j} {} = {} & \left[\beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K \nonumber \\
  & + e^{-4 \phi} \left \{-2 \alpha \bar{D}_{i} \bar{D}_{j} \phi + 4 \alpha \bar{D}_{i} \phi \bar{D}_{j} \phi  + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_{i} \bar{D}_{j} \alpha + \alpha \bar{R}_{i j} \right \}^{\text{TF}} \; , \\
  \partial_t \phi {} = {} & \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) \; , \\
  \partial_{t} K {} = {} & \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^{2} + \alpha \bar{A}_{i j} \bar{A}^{i j} - e^{-4 \phi} \left (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi \right ) \; , \\
  \partial_t \bar{\Lambda}^{i} {} = {} & \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} \nonumber \\
  & - 2 \bar{A}^{i j} \left (\partial_{j} \alpha - 6 \partial_{j} \phi \right ) + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i}  -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
\end{align}

where the terms moved from the right-hand sides to the left-hand sides are enclosed in square braces. 

Notice that the shift advection operator $\beta^k \partial_k \left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\}$ appears on the right-hand side of *every* expression. As the shift determines how the spatial coordinates $x^i$ move on the next 3D slice of our 4D manifold, we find that representing $\partial_k$ in these shift advection terms via an *upwinded* finite difference stencil results in far lower numerical errors. This trick is implemented below in all shift advection terms.

As discussed in the [NRPy+ tutorial notebook on BSSN quantities](Tutorial-BSSN_quantities.ipynb), tensorial expressions can diverge at coordinate singularities, so each tensor in the set of BSSN variables

$$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$

is written in terms of the corresponding rescaled quantity in the set

$$\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\},$$ 

respectively, as defined in the [BSSN quantities tutorial](Tutorial-BSSN_quantities.ipynb). 

<a id='gaugeconditions'></a>

# Step 2: Time-evolution equations for the BSSN gauge quantities $\alpha$ and $\beta^i$ \[Back to [top](#toc)\]
$$\label{gaugeconditions}$$

As described in the **Background Reading/Lectures** linked to above, the gauge quantities $\alpha$ and $\beta^i$ specify how coordinate time and spatial points adjust from one spatial hypersurface to the next, in our 3+1 decomposition of Einstein's equations.

As choosing $\alpha$ and $\beta^i$ is equivalent to choosing coordinates for where we sample our solution to Einstein's equations, we are completely free to choose $\alpha$ and $\beta^i$ on any given spatial hypersurface. It has been found that fixing $\alpha$ and $\beta^i$ to constant values in the context of dynamical spacetimes results in instabilities, so we generally need to define expressions for $\partial_t \alpha$ and $\partial_t \beta^i$ and couple these equations to the rest of the BSSN time-evolution equations.

Though we are free to choose the form of the right-hand sides of the gauge time evolution equations, very few have been found robust in the presence of (puncture) black holes.

The most commonly adopted gauge conditions for BSSN (i.e., time-evolution equations for the BSSN gauge quantities $\alpha$ and $\beta^i$) are the 

* $1+\log$ lapse condition:
$$
\partial_0 \alpha = -2 \alpha K
$$
* Second-order Gamma-driving shift condition:
\begin{align}
\partial_0 \beta^i &= B^{i} \\
\partial_0 B^i &= \frac{3}{4} \partial_{0} \bar{\Lambda}^{i} - \eta B^{i},
\end{align}

where $\partial_0$ is the advection operator; i.e., $\partial_0 A^i = \partial_t A^i - \beta^j \partial_j A^i$. Note that $\partial_{0} \bar{\Lambda}^{i}$ in the right-hand side of the $\partial_{0} B^{i}$ equation is computed by adding $\beta^j \partial_j \bar{\Lambda}^i$ to the right-hand side expression given for $\partial_t \bar{\Lambda}^i$, so no explicit time dependence occurs in the right-hand sides of the BSSN evolution equations and the Method of Lines can be applied directly.

While it is incredibly robust in Cartesian coordinates, [Brown](https://arxiv.org/abs/0902.3652) pointed out that the above time-evolution equation for the shift is not covariant. In fact, we have found this non-covariant version to result in very poor results when solving Einstein's equations in spherical coordinates for a spinning black hole with the spin axis pointed in the $\hat{x}$ direction. Therefore we adopt Brown's covariant version as described in the [**full time-evolution equations for the BSSN gauge quantities $\alpha$ and $\beta^i$ tutorial notebook**](Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb).

<a id='constraintequations'></a>

# Step 3: The BSSN constraint equations \[Back to [top](#toc)\]
$$\label{constraintequations}$$

In a way analogous to Maxwell's equations, the BSSN decomposition of Einstein's equations is written as a set of time-evolution equations and a set of constraint equations. In this step, we present the BSSN constraints

\begin{align}
\mathcal{H} &= 0 \\
\mathcal{M^i} &= 0,
\end{align}

where $\mathcal{H}=0$ is the **Hamiltonian constraint**, and $\mathcal{M^i} = 0$ is the **momentum constraint**. When constructing our spacetime from the initial data, one spatial hypersurface at a time, to confirm that at a given time, the Hamiltonian and momentum constraint violations converge to zero as expected with increased numerical resolution.

<a id='hamiltonianconstraint'></a>

## Step 3.a: The Hamiltonian constraint $\mathcal{H}$ \[Back to [top](#toc)\]
$$\label{hamiltonianconstraint}$$

The Hamiltonian constraint is written (Eq. 13 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
$$
\mathcal{H} = \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right)
$$

<a id='momentumconstraint'></a>

## Step 3.b: The momentum constraint $\mathcal{M}^i$ \[Back to [top](#toc)\]
$$\label{momentumconstraint}$$

The momentum constraint is written (Eq. 47 of [Ruchlin, Etienne, & Baumgarte](https://arxiv.org/pdf/1712.07658.pdf)):

$$    \mathcal{M}^i = e^{-4\phi} \left(
    \frac{1}{\sqrt{\bar{\gamma}}} \hat{D}_j\left(\sqrt{\bar{\gamma}}\bar{A}^{ij}\right) +
    6 \bar{A}^{ij}\partial_j \phi -
    \frac{2}{3} \bar{\gamma}^{ij}\partial_j K +
    \bar{A}^{jk} \Delta\Gamma^i_{jk} + \bar{A}^{ik} \Delta\Gamma^j_{jk}\right)
$$

Notice the momentum constraint as written in [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf) is missing a term, as described in [Ruchlin, Etienne, & Baumgarte](https://arxiv.org/pdf/1712.07658.pdf).

<a id='gammaconstraint'></a>

# Step 4: The BSSN algebraic constraint: $\hat{\gamma}=\bar{\gamma}$ \[Back to [top](#toc)\]
$$\label{gammaconstraint}$$

[Brown](https://arxiv.org/abs/0902.3652)'s covariant Lagrangian formulation of BSSN, which we adopt, requires that $\partial_t \bar{\gamma} = 0$, where $\bar{\gamma}=\det \bar{\gamma}_{ij}$. We generally choose to set $\bar{\gamma}=\hat{\gamma}$ in our initial data.

Numerical errors will cause $\bar{\gamma}$ to deviate from a constant in time. This actually disrupts the hyperbolicity of the PDEs (causing crashes), so to cure this, we adjust $\bar{\gamma}_{ij}$ at the end of each Runge-Kutta timestep, so that its determinant satisfies $\bar{\gamma}=\hat{\gamma}$ at all times. We adopt the following, rather standard prescription (Eq. 53 of [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658)):

$$
\bar{\gamma}_{ij} \to \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{ij}.
$$
Notice the expression on the right is guaranteed to have a determinant equal to $\hat{\gamma}$.

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_formulation.pdf](Tutorial-BSSN_formulation.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_formulation")
```

    Created Tutorial-BSSN_formulation.tex, and compiled LaTeX file to PDF file
        Tutorial-BSSN_formulation.pdf

