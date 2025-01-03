<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# The BSSN Time-Evolution Equations

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook focuses on a detailed construction and explanation of the BSSN formulation's time evolution equations, given in [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658), within the NRPy+ environment. It incorporates an upwinded finite difference stencil approach in handling the shift advection terms, enhancing the numerical precision.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** All expressions generated in this module have been validated against a trusted code (the original NRPy+/SENR code, which itself was validated against [Baumgarte's code](https://arxiv.org/abs/1211.6632)).

### NRPy+ Source Code for this module: [BSSN/BSSN_RHS.py](../edit/BSSN/BSSN_RHSs.py)

## Introduction:
This module documents and constructs the time evolution equations of the BSSN formulation of Einstein's equations, as defined  in [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658) (see also [Baumgarte, Montero, Cordero-Carrión, and Müller (2012)](https://arxiv.org/abs/1211.6632)). 

**This module is part of the following set of NRPy+ tutorial notebooks on the BSSN formulation of general relativity:**

* An overview of the BSSN formulation of Einstein's equations, as well as links for background reading/lectures, are provided in [the NRPy+ tutorial notebook on the BSSN formulation](Tutorial-BSSN_formulation.ipynb).
* Basic BSSN quantities are defined in the [BSSN quantities NRPy+ tutorial notebook](Tutorial-BSSN_quantities.ipynb).
* Other BSSN equation tutorial notebooks:
    * [Time-evolution equations the BSSN gauge quantities $\alpha$, $\beta^i$, and $B^i$](Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb).
    * [BSSN Hamiltonian and momentum constraints](Tutorial-BSSN_constraints.ipynb)
    * [Enforcing the $\bar{\gamma} = \hat{\gamma}$ constraint](Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.ipynb)

### A Note on Notation

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

0. [Preliminaries](#bssntimeevolequations): BSSN time-evolution equations, as described in the [BSSN formulation NRPy+ tutorial notebook](Tutorial-BSSN_formulation.ipynb)
1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#gammabar): Right-hand side of $\partial_t \bar{\gamma}_{ij}$
    1. [Step 2.a](#term1_partial_gamma): Term 1 of $\partial_t \bar{\gamma}_{i j}$
    1. [Step 2.b](#term2_partial_gamma): Term 2 of $\partial_t \bar{\gamma}_{i j}$
    1. [Step 2.c](#term3_partial_gamma): Term 3 of $\partial_t \bar{\gamma}_{i j}$
1. [Step 3](#abar): Right-hand side of $\partial_t \bar{A}_{ij}$
    1. [Step 3.a](#term1_partial_upper_a): Term 1 of $\partial_t \bar{A}_{i j}$
    1. [Step 3.c](#term2_partial_upper_a): Term 2 of $\partial_t \bar{A}_{i j}$
    1. [Step 3.c](#term3_partial_upper_a): Term 3 of $\partial_t \bar{A}_{i j}$
1. [Step 4](#cf): Right-hand side of $\partial_t \phi \to \partial_t (\text{cf})$
1. [Step 5](#trk): Right-hand side of $\partial_t \text{tr} K$
1. [Step 6](#lambdabar): Right-hand side of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.a](#term1_partial_lambda): Term 1 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.b](#term2_partial_lambda): Term 2 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.c](#term3_partial_lambda): Term 3 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.d](#term4_partial_lambda): Term 4 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.e](#term5_partial_lambda): Term 5 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.f](#term6_partial_lambda): Term 6 of $\partial_t \bar{\Lambda}^i$
    1. [Step 6.g](#term7_partial_lambda): Term 7 of $\partial_t \bar{\Lambda}^i$
1. [Step 7](#rescalingrhss): Rescaling the BSSN right-hand sides; rewriting them in terms of the rescaled quantities $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$
1. [Step 8](#code_validation): Code Validation against `BSSN.BSSN_RHSs` NRPy+ module
1. [Step 9](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='bssntimeevolequations'></a>

# Preliminaries: BSSN time-evolution equations \[Back to [top](#toc)\]
$$\label{bssntimeevolequations}$$

As described in the [BSSN formulation NRPy+ tutorial notebook](Tutorial-BSSN_formulation.ipynb), the BSSN time-evolution equations are given by

\begin{align}
  \partial_t \bar{\gamma}_{i j} {} = {} & \left[\beta^k \partial_k \bar{\gamma}_{ij} + \partial_i \beta^k \bar{\gamma}_{kj} + \partial_j \beta^k \bar{\gamma}_{ik} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right ) - 2 \alpha \bar{A}_{i j} \; , \\
  \partial_t \bar{A}_{i j} {} = {} & \left[\beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K \nonumber \\
  & + e^{-4 \phi} \left \{-2 \alpha \bar{D}_{i} \bar{D}_{j} \phi + 4 \alpha \bar{D}_{i} \phi \bar{D}_{j} \phi  + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_{i} \bar{D}_{j} \alpha + \alpha \bar{R}_{i j} \right \}^{\text{TF}} \; , \\
  \partial_t \phi {} = {} & \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) \; , \\
  \partial_{t} K {} = {} & \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^{2} + \alpha \bar{A}_{i j} \bar{A}^{i j} - e^{-4 \phi} \left (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi \right ) \; , \\
  \partial_t \bar{\Lambda}^{i} {} = {} & \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} \nonumber \\
  & - 2 \bar{A}^{i j} \left (\partial_{j} \alpha - 6 \partial_{j} \phi \right ) + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i}  -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
\end{align}

where the Lie derivative terms (often seen on the left-hand side of these equations) are enclosed in square braces.

Notice that the shift advection operator $\beta^k \partial_k \left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}\right\}$ appears on the right-hand side of *every* expression. As the shift determines how the spatial coordinates $x^i$ move on the next 3D slice of our 4D manifold, we find that representing $\partial_k$ in these shift advection terms via an *upwinded* finite difference stencil results in far lower numerical errors. This trick is implemented below in all shift advection terms. Upwinded derivatives are indicated in NRPy+ by the `_dupD` variable suffix.


As discussed in the [NRPy+ tutorial notebook on BSSN quantities](Tutorial-BSSN_quantities.ipynb), tensorial expressions can diverge at coordinate singularities, so each tensor in the set of BSSN variables

$$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$

is written in terms of the corresponding rescaled quantity in the set

$$\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\},$$ 

respectively, as defined in the [BSSN quantities tutorial](Tutorial-BSSN_quantities.ipynb). 

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from NRPy+:


```python
# Step 1.a: import all needed modules from Python/NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python modules for multiplatform OS-level functions
import BSSN.BSSN_quantities as Bq # NRPy+: Basic BSSN quantities

# Step 1.b: Set the coordinate system for the numerical grid
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")

# Step 1.c: Given the chosen coordinate system, set up
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
rfm.reference_metric()

# Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
#           a 3+1-dimensional decomposition of the general
#           relativistic field equations)
DIM = 3

# Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
Bq.BSSN_basic_tensors()
gammabarDD = Bq.gammabarDD
AbarDD     = Bq.AbarDD
LambdabarU = Bq.LambdabarU
trK   = Bq.trK
alpha = Bq.alpha
betaU = Bq.betaU

# Step 1.f: Import all needed rescaled BSSN tensors:
cf  = Bq.cf
lambdaU = Bq.lambdaU
```

<a id='gammabar'></a>

# Step 2: Right-hand side of $\partial_t \bar{\gamma}_{ij}$ \[Back to [top](#toc)\]
$$\label{gammabar}$$

Let's start with

$$
\partial_t \bar{\gamma}_{i j} = 
{\underbrace {\textstyle \left[\beta^k \partial_k \bar{\gamma}_{ij} + \partial_i \beta^k \bar{\gamma}_{kj} + \partial_j \beta^k \bar{\gamma}_{ik} \right]}_{\text{Term 1}}} + 
{\underbrace {\textstyle \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )}_{\text{Term 2}}}
{\underbrace {\textstyle -2 \alpha \bar{A}_{i j}}_{\text{Term 3}}}.
$$


First, note that the term containing $\bar{A}_{k}^{k}$ is there to drive $\bar{A}_{k}^{k}$ to zero. It was implemented in this form in T. Baumgarte's BSSN-in-spherical-coordinates code, against which NRPy+ and SENR were originally validated. You will find this term in [Brown](https://arxiv.org/abs/0902.3652)'s covariant formulation of the BSSN time-evolution equations (see third term in Eq 21a).

<a id='term1_partial_gamma'></a>

## Step 2.a: Term 1 of $\partial_t \bar{\gamma}_{i j}$  \[Back to [top](#toc)\] 
$$\label{term1_partial_gamma}$$

Term 1 of $\partial_t \bar{\gamma}_{i j} =$ `gammabar_rhsDD[i][j]`: $\beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}$


First we import derivative expressions for betaU defined in the [NRPy+ BSSN quantities tutorial notebook](Tutorial-BSSN_quantities.ipynb)


```python
# Step 2.a.i: Import derivative expressions for betaU defined in the BSSN.BSSN_quantities module:
Bq.betaU_derivs()
betaU_dD = Bq.betaU_dD
betaU_dupD = Bq.betaU_dupD
betaU_dDD = Bq.betaU_dDD
# Step 2.a.ii: Import derivative expression for gammabarDD
Bq.gammabar__inverse_and_derivs()
gammabarDD_dupD = Bq.gammabarDD_dupD

# Step 2.a.iii: First term of \partial_t \bar{\gamma}_{i j} right-hand side:
# \beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}
gammabar_rhsDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            gammabar_rhsDD[i][j] += betaU[k]*gammabarDD_dupD[i][j][k] + betaU_dD[k][i]*gammabarDD[k][j] \
                                                                      + betaU_dD[k][j]*gammabarDD[i][k]
```

<a id='term2_partial_gamma'></a>

## Step 2.b: Term 2 of $\partial_t \bar{\gamma}_{i j}$ \[Back to [top](#toc)\]
$$\label{term2_partial_gamma}$$

Term 2 of $\partial_t \bar{\gamma}_{i j} =$ `gammabar_rhsDD[i][j]`: $\frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )$

Let's first convert this expression to be in terms of the evolved variables $a_{ij}$ and $\mathcal{B}^i$, starting with $\bar{A}_{ij} = a_{ij} \text{ReDD[i][j]}$. Then $\bar{A}^k_{k} = \bar{\gamma}^{ij} \bar{A}_{ij}$, and we have already defined $\bar{\gamma}^{ij}$ in terms of the evolved quantity $h_{ij}$.

Next, we wish to compute 

$$\bar{D}_{k} \beta^{k} = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}},$$

where $\bar{\gamma}$ is the determinant of the conformal metric $\bar{\gamma}_{ij}$. ***Exercise to student: Prove the above relation.***
[Solution.](https://physics.stackexchange.com/questions/81453/general-relativity-christoffel-symbol-identity)

Usually (i.e., so long as we make the parameter choice `detgbarOverdetghat_equals_one = False` ) we will choose $\bar{\gamma}=\hat{\gamma}$, so $\bar{\gamma}$ will in general possess coordinate singularities. Thus we would prefer to rewrite derivatives of $\bar{\gamma}$ in terms of derivatives of $\bar{\gamma}/\hat{\gamma} = 1$.


```python
# Step 2.b.i: First import \bar{A}_{ij} = AbarDD[i][j], and its contraction trAbar = \bar{A}^k_k
#           from BSSN.BSSN_quantities
Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()
trAbar = Bq.trAbar

# Step 2.b.ii: Import detgammabar quantities from BSSN.BSSN_quantities:
Bq.detgammabar_and_derivs()
detgammabar = Bq.detgammabar
detgammabar_dD = Bq.detgammabar_dD

# Step 2.b.ii: Compute the contraction \bar{D}_k \beta^k = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}
Dbarbetacontraction = sp.sympify(0)
for k in range(DIM):
    Dbarbetacontraction += betaU_dD[k][k] + betaU[k]*detgammabar_dD[k]/(2*detgammabar)

# Step 2.b.iii: Second term of \partial_t \bar{\gamma}_{i j} right-hand side:
# \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )
for i in range(DIM):
    for j in range(DIM):
        gammabar_rhsDD[i][j] += sp.Rational(2,3)*gammabarDD[i][j]*(alpha*trAbar - Dbarbetacontraction)
```

<a id='term3_partial_gamma'></a>

## Step 2.c: Term 3 of $\partial_t \bar{\gamma}_{i j}$ \[Back to [top](#toc)\] 
$$\label{term3_partial_gamma}$$

Term 3 of $\partial_t \bar{\gamma}_{i j}$ = `gammabar_rhsDD[i][j]`: $-2 \alpha \bar{A}_{ij}$



```python
# Step 2.c: Third term of \partial_t \bar{\gamma}_{i j} right-hand side:
# -2 \alpha \bar{A}_{ij}
for i in range(DIM):
    for j in range(DIM):
        gammabar_rhsDD[i][j] += -2*alpha*AbarDD[i][j]
```

<a id='abar'></a>

# Step 3: Right-hand side of $\partial_t \bar{A}_{ij}$ \[Back to [top](#toc)\]
$$\label{abar}$$

$$\partial_t \bar{A}_{i j} = 
{\underbrace {\textstyle \left[\beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik} \right]}_{\text{Term 1}}}
{\underbrace {\textstyle - \frac{2}{3} \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K}_{\text{Term 2}}} + 
{\underbrace {\textstyle e^{-4 \phi} \left \{-2 \alpha \bar{D}_{i} \bar{D}_{j} \phi + 4 \alpha \bar{D}_{i} \phi \bar{D}_{j} \phi  + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_{i} \bar{D}_{j} \alpha + \alpha \bar{R}_{i j} \right \}^{\text{TF}}}_{\text{Term 3}}}$$

<a id='term1_partial_upper_a'></a>

## Step 3.a: Term 1 of $\partial_t \bar{A}_{i j}$ \[Back to [top](#toc)\]
$$\label{term1_partial_upper_a}$$

Term 1 of $\partial_t \bar{A}_{i j}$ = `Abar_rhsDD[i][j]`: $\left[\beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik} \right]$


Notice the first subexpression has a $\beta^k \partial_k A_{ij}$ advection term, which will be upwinded.


```python
# Step 3.a: First term of \partial_t \bar{A}_{i j}:
# \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}

AbarDD_dupD = Bq.AbarDD_dupD # From Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()

Abar_rhsDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            Abar_rhsDD[i][j] += betaU[k]*AbarDD_dupD[i][j][k] + betaU_dD[k][i]*AbarDD[k][j] \
                                                              + betaU_dD[k][j]*AbarDD[i][k]
```

<a id='term2_partial_upper_a'></a>

## Step 3.b: Term 2 of $\partial_t \bar{A}_{i j}$ \[Back to [top](#toc)\]
$$\label{term2_partial_upper_a}$$

Term 2 of $\partial_t \bar{A}_{i j}$ = `Abar_rhsDD[i][j]`: $- \frac{2}{3} \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} \bar{A}^{k}_{j} + \alpha \bar{A}_{i j} K$


Note that $\bar{D}_{k} \beta^{k}$ was already defined as `Dbarbetacontraction`.


```python
# Step 3.b: Second term of \partial_t \bar{A}_{i j}:
# - (2/3) \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K
gammabarUU = Bq.gammabarUU # From Bq.gammabar__inverse_and_derivs()
AbarUD     = Bq.AbarUD     # From Bq.AbarUU_AbarUD_trAbar()
for i in range(DIM):
    for j in range(DIM):
        Abar_rhsDD[i][j] += -sp.Rational(2,3)*AbarDD[i][j]*Dbarbetacontraction + alpha*AbarDD[i][j]*trK
        for k in range(DIM):
            Abar_rhsDD[i][j] += -2*alpha * AbarDD[i][k]*AbarUD[k][j]
```

<a id='term3_partial_upper_a'></a>

## Step 3.c: Term 3 of $\partial_t \bar{A}_{i j}$ \[Back to [top](#toc)\]
$$\label{term3_partial_upper_a}$$


Term 3 of $\partial_t \bar{A}_{i j}$ = `Abar_rhsDD[i][j]`: $e^{-4 \phi} \left \{-2 \alpha \bar{D}_{i} \bar{D}_{j} \phi + 4 \alpha \bar{D}_{i} \phi \bar{D}_{j} \phi  + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_{i} \bar{D}_{j} \alpha + \alpha \bar{R}_{i j} \right \}^{\text{TF}}$ 

The first covariant derivatives of $\phi$ and $\alpha$ are simply partial derivatives. However, $\phi$ is not a gridfunction; `cf` is. cf = $W$ (default value) denotes that the evolved variable is $W=e^{-2 \phi}$, which results in smoother spacetime fields around puncture black holes (desirable).


```python
# Step 3.c.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
Bq.phi_and_derivs()
phi_dD   = Bq.phi_dD
phi_dupD = Bq.phi_dupD
exp_m4phi  = Bq.exp_m4phi
phi_dBarD  = Bq.phi_dBarD  # phi_dBarD = Dbar_i phi = phi_dD (since phi is a scalar)
phi_dBarDD = Bq.phi_dBarDD # phi_dBarDD = Dbar_i Dbar_j phi (covariant derivative)

# Step 3.c.ii: Define RbarDD
Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
RbarDD = Bq.RbarDD

# Step 3.c.iii: Define first and second derivatives of \alpha, as well as
#         \bar{D}_i \bar{D}_j \alpha, which is defined just like phi
alpha_dD = ixp.declarerank1("alpha_dD")
alpha_dDD = ixp.declarerank2("alpha_dDD","sym01")
alpha_dBarD = alpha_dD
alpha_dBarDD = ixp.zerorank2()
GammabarUDD = Bq.GammabarUDD # Defined in Bq.gammabar__inverse_and_derivs()
for i in range(DIM):
    for j in range(DIM):
        alpha_dBarDD[i][j] = alpha_dDD[i][j]
        for k in range(DIM):
            alpha_dBarDD[i][j] += - GammabarUDD[k][i][j]*alpha_dD[k]

# Step 3.c.iv: Define the terms in curly braces:
curlybrackettermsDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        curlybrackettermsDD[i][j] = -2*alpha*phi_dBarDD[i][j] + 4*alpha*phi_dBarD[i]*phi_dBarD[j] \
                                    +2*alpha_dBarD[i]*phi_dBarD[j] \
                                    +2*alpha_dBarD[j]*phi_dBarD[i] \
                                    -alpha_dBarDD[i][j] + alpha*RbarDD[i][j]

# Step 3.c.v: Compute the trace:
curlybracketterms_trace = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        curlybracketterms_trace += gammabarUU[i][j]*curlybrackettermsDD[i][j]

# Step 3.c.vi: Third and final term of Abar_rhsDD[i][j]:
for i in range(DIM):
    for j in range(DIM):
        Abar_rhsDD[i][j] += exp_m4phi*(curlybrackettermsDD[i][j] - \
                                       sp.Rational(1,3)*gammabarDD[i][j]*curlybracketterms_trace)
```

<a id='cf'></a>

# Step 4: Right-hand side of $\partial_t \phi \to \partial_t (\text{cf})$ \[Back to [top](#toc)\]
$$\label{cf}$$

$$\partial_t \phi = 
{\underbrace {\textstyle \left[\beta^k \partial_k \phi \right]}_{\text{Term 1}}} + 
{\underbrace {\textstyle \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right)}_{\text{Term 2}}}$$

The right-hand side of $\partial_t \phi$ is trivial except for the fact that the actual evolved variable is `cf` (short for conformal factor), which could represent
* cf = $\phi$
* cf = $W = e^{-2 \phi}$ (default)
* cf = $\chi = e^{-4 \phi}$

Thus we are actually computing the right-hand side of the equation $\partial_t $cf, which is related to $\partial_t \phi$ via simple relations:
* cf = $\phi$: $\partial_t $cf$ = \partial_t \phi$ (unchanged)
* cf = $W$: $\partial_t $cf$ = \partial_t (e^{-2 \phi}) = -2 e^{-2\phi}\partial_t \phi = -2 W \partial_t \phi$. Thus we need to multiply the right-hand side by $-2 W = -2$cf when cf = $W$.
* cf = $\chi$: Same argument as for $W$, except the right-hand side must be multiplied by $-4 \chi=-4$cf.


```python
# Step 4: Right-hand side of conformal factor variable "cf". Supported
#          options include: cf=phi, cf=W=e^(-2*phi) (default), and cf=chi=e^(-4*phi)
# \partial_t phi = \left[\beta^k \partial_k \phi \right] <- TERM 1
#                  + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) <- TERM 2

cf_rhs = sp.Rational(1,6) * (Dbarbetacontraction - alpha*trK) # Term 2
for k in range(DIM):
    cf_rhs += betaU[k]*phi_dupD[k] # Term 1

# Next multiply to convert phi_rhs to cf_rhs.
if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "phi":
    pass # do nothing; cf_rhs = phi_rhs
elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "W":
    cf_rhs *= -2*cf # cf_rhs = -2*cf*phi_rhs
elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "chi":
    cf_rhs *= -4*cf # cf_rhs = -4*cf*phi_rhs
else:
    print("Error: EvolvedConformalFactor_cf == "+
          par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf")+" unsupported!")
    sys.exit(1)
```

<a id='trk'></a>

# Step 5: Right-hand side of $\partial_t K$ \[Back to [top](#toc)\]
$$\label{trk}$$

$$
\partial_{t} K = 
{\underbrace {\textstyle \left[\beta^i \partial_i K \right]}_{\text{Term 1}}} + 
{\underbrace {\textstyle \frac{1}{3} \alpha K^{2}}_{\text{Term 2}}} +
{\underbrace {\textstyle \alpha \bar{A}_{i j} \bar{A}^{i j}}_{\text{Term 3}}}
{\underbrace {\textstyle - e^{-4 \phi} \left (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi \right )}_{\text{Term 4}}}
$$


```python
# Step 5: right-hand side of trK (trace of extrinsic curvature):
# \partial_t K = \beta^k \partial_k K <- TERM 1
#           + \frac{1}{3} \alpha K^{2} <- TERM 2
#           + \alpha \bar{A}_{i j} \bar{A}^{i j} <- TERM 3
#           - - e^{-4 \phi} (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi ) <- TERM 4
# TERM 2:
trK_rhs = sp.Rational(1,3)*alpha*trK*trK
trK_dupD = ixp.declarerank1("trK_dupD")
for i in range(DIM):
    # TERM 1:
    trK_rhs += betaU[i]*trK_dupD[i]
for i in range(DIM):
    for j in range(DIM):
        # TERM 4:
        trK_rhs += -exp_m4phi*gammabarUU[i][j]*(alpha_dBarDD[i][j] + 2*alpha_dBarD[j]*phi_dBarD[i])
AbarUU = Bq.AbarUU # From Bq.AbarUU_AbarUD_trAbar()
for i in range(DIM):
    for j in range(DIM):
        # TERM 3:
        trK_rhs += alpha*AbarDD[i][j]*AbarUU[i][j]
```

<a id='lambdabar'></a>

# Step 6: Right-hand side of $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{lambdabar}$$

\begin{align}
\partial_t \bar{\Lambda}^{i} &= 
{\underbrace {\textstyle \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right]}_{\text{Term 1}}} + 
{\underbrace {\textstyle \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i}}_{\text{Term 2}}} + 
{\underbrace {\textstyle \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}}_{\text{Term 3}}} + 
{\underbrace {\textstyle \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}}_{\text{Term 4}}} \nonumber \\
  & 
{\underbrace {\textstyle - 2 \bar{A}^{i j} \left (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi \right )}_{\text{Term 5}}} + 
{\underbrace {\textstyle 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i}}_{\text{Term 6}}} 
{\underbrace {\textstyle -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K}_{\text{Term 7}}}
\end{align}

<a id='term1_partial_lambda'></a>

## Step 6.a: Term 1 of  $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term1_partial_lambda}$$

Term 1 of  $\partial_t \bar{\Lambda}^{i}$: $\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k$

Computing this term requires that we define $\bar{\Lambda}^i$ and $\bar{\Lambda}^i_{,j}$ in terms of the rescaled (i.e., actual evolved) variable $\lambda^i$ and derivatives: 
\begin{align}
\bar{\Lambda}^i &= \lambda^i \text{ReU[i]} \\
\bar{\Lambda}^i_{,\ j} &= \lambda^i_{,\ j} \text{ReU[i]} + \lambda^i \text{ReUdD[i][j]} 
\end{align}


```python
# Step 6: right-hand side of \partial_t \bar{\Lambda}^i:
# \partial_t \bar{\Lambda}^i = \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k <- TERM 1
#                            + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} <- TERM 2
#                            + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} <- TERM 3
#                            + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} <- TERM 4
#                            - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \partial_{j} \phi) <- TERM 5
#                            + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i} <- TERM 6
#                            - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K <- TERM 7

# Step 6.a: Term 1 of \partial_t \bar{\Lambda}^i: \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
# First we declare \bar{\Lambda}^i and \bar{\Lambda}^i_{,j} in terms of \lambda^i and \lambda^i_{,j}
LambdabarU_dupD = ixp.zerorank2()
lambdaU_dupD = ixp.declarerank2("lambdaU_dupD","nosym")
for i in range(DIM):
    for j in range(DIM):
        LambdabarU_dupD[i][j] = lambdaU_dupD[i][j]*rfm.ReU[i] + lambdaU[i]*rfm.ReUdD[i][j]

Lambdabar_rhsU = ixp.zerorank1()
for i in range(DIM):
    for k in range(DIM):
        Lambdabar_rhsU[i] += betaU[k]*LambdabarU_dupD[i][k] - betaU_dD[i][k]*LambdabarU[k] # Term 1
```

<a id='term2_partial_lambda'></a>

## Step 6.b: Term 2 of $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term2_partial_lambda}$$

Term 2 of $\partial_t \bar{\Lambda}^{i}$: $\bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i}$

This is a relatively difficult term to compute, as it requires we evaluate the second covariant derivative of the shift vector, with respect to the hatted (i.e., reference) metric.

Based on the definition of the covariant derivative, we have
$$
\hat{D}_{k} \beta^{i} = \beta^i_{,k} + \hat{\Gamma}^i_{mk} \beta^m
$$

Since $\hat{D}_{k} \beta^{i}$ is a tensor, the covariant derivative of this will have the same indexing as a tensor $T_k^i$:

$$
\hat{D}_{j} T^i_k = T^i_{k,j} + \hat{\Gamma}^i_{dj} T^d_k - \hat{\Gamma}^d_{kj} T^i_d.
$$

Therefore,
\begin{align}
\hat{D}_{j} \left(\hat{D}_{k} \beta^{i}\right) &= \left(\beta^i_{,k} + \hat{\Gamma}^i_{mk} \beta^m\right)_{,j} + \hat{\Gamma}^i_{dj} \left(\beta^d_{,k} + \hat{\Gamma}^d_{mk} \beta^m\right) - \hat{\Gamma}^d_{kj} \left(\beta^i_{,d} + \hat{\Gamma}^i_{md} \beta^m\right) \\
&= \beta^i_{,kj} + \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j} + \hat{\Gamma}^i_{dj}\beta^d_{,k} + \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \beta^i_{,d} - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m \\
&= {\underbrace {\textstyle \beta^i_{,kj}}_{\text{Term 2a}}} +
{\underbrace {\textstyle \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j} + \hat{\Gamma}^i_{dj}\beta^d_{,k} - \hat{\Gamma}^d_{kj} \beta^i_{,d}}_{\text{Term 2b}}} +
{\underbrace {\textstyle \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m}_{\text{Term 2c}}},
\end{align}

where 
$$
\text{Term 2} = \bar{\gamma}^{jk} \left(\text{Term 2a} + \text{Term 2b} + \text{Term 2c}\right)
$$


```python
# Step 6.b: Term 2 of \partial_t \bar{\Lambda}^i = \bar{\gamma}^{jk} (Term 2a + Term 2b + Term 2c)
# Term 2a: \bar{\gamma}^{jk} \beta^i_{,kj}
Term2aUDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            Term2aUDD[i][j][k] += betaU_dDD[i][k][j]
# Term 2b: \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j}
#          + \hat{\Gamma}^i_{dj}\beta^d_{,k} - \hat{\Gamma}^d_{kj} \beta^i_{,d}
Term2bUDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for m in range(DIM):
                Term2bUDD[i][j][k] += rfm.GammahatUDDdD[i][m][k][j]*betaU[m]    \
                                      + rfm.GammahatUDD[i][m][k]*betaU_dD[m][j] \
                                      + rfm.GammahatUDD[i][m][j]*betaU_dD[m][k] \
                                      - rfm.GammahatUDD[m][k][j]*betaU_dD[i][m]
# Term 2c: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m
Term2cUDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for m in range(DIM):
                for d in range(DIM):
                    Term2cUDD[i][j][k] += ( rfm.GammahatUDD[i][d][j]*rfm.GammahatUDD[d][m][k] \
                                           -rfm.GammahatUDD[d][k][j]*rfm.GammahatUDD[i][m][d])*betaU[m]

Lambdabar_rhsUpieceU = ixp.zerorank1()

# Put it all together to get Term 2:
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            Lambdabar_rhsU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
            Lambdabar_rhsUpieceU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
```

<a id='term3_partial_lambda'></a>

## Step 6.c: Term 3 of  $\partial_t \bar{\Lambda}^{i}$: $\frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}$ \[Back to [top](#toc)\]
$$\label{term3_partial_lambda}$$

Term 3 of $\partial_t \bar{\Lambda}^{i}$: $\frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}$
 
This term is the simplest to implement, as $\bar{D}_{j} \beta^{j}$ and $\Delta^i$ have already been defined, as `Dbarbetacontraction` and `DGammaU[i]`, respectively:


```python
# Step 6.c: Term 3 of \partial_t \bar{\Lambda}^i:
#    \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}
DGammaU = Bq.DGammaU # From Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
for i in range(DIM):
    Lambdabar_rhsU[i] += sp.Rational(2,3)*DGammaU[i]*Dbarbetacontraction # Term 3
```

<a id='term4_partial_lambda'></a>

## Step 6.d: Term 4 of  $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term4_partial_lambda}$$

$\partial_t \bar{\Lambda}^{i}$: $\frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}$

Recall first that

$$\bar{D}_{k} \beta^{k} = \beta^k_{,\ k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}},$$
which is a scalar, so

\begin{align}
\bar{D}_m \bar{D}_{j} \beta^{j} &= \left(\beta^k_{,\ k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}\right)_{,m} \\
&= \beta^k_{\ ,km} + \frac{\beta^k_{\ ,m} \bar{\gamma}_{,k} + \beta^k \bar{\gamma}_{\ ,km}}{2 \bar{\gamma}} - \frac{\beta^k \bar{\gamma}_{,k} \bar{\gamma}_{,m}}{2 \bar{\gamma}^2}
\end{align}

Thus, 
\begin{align}
\bar{D}^i \bar{D}_{j} \beta^{j} 
&= \bar{\gamma}^{im} \bar{D}_m \bar{D}_{j} \beta^{j} \\
&= \bar{\gamma}^{im} \left(\beta^k_{\ ,km} + \frac{\beta^k_{\ ,m} \bar{\gamma}_{,k} + \beta^k \bar{\gamma}_{\ ,km}}{2 \bar{\gamma}} - \frac{\beta^k \bar{\gamma}_{,k} \bar{\gamma}_{,m}}{2 \bar{\gamma}^2} \right)
\end{align}


```python
# Step 6.d: Term 4 of \partial_t \bar{\Lambda}^i:
#           \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}
detgammabar_dDD = Bq.detgammabar_dDD # From Bq.detgammabar_and_derivs()
Dbarbetacontraction_dBarD = ixp.zerorank1()
for k in range(DIM):
    for m in range(DIM):
        Dbarbetacontraction_dBarD[m] += betaU_dDD[k][k][m] + \
                                        (betaU_dD[k][m]*detgammabar_dD[k] +
                                         betaU[k]*detgammabar_dDD[k][m])/(2*detgammabar) \
                                        -betaU[k]*detgammabar_dD[k]*detgammabar_dD[m]/(2*detgammabar*detgammabar)
for i in range(DIM):
    for m in range(DIM):
        Lambdabar_rhsU[i] += sp.Rational(1,3)*gammabarUU[i][m]*Dbarbetacontraction_dBarD[m]
```

<a id='term5_partial_lambda'></a>

## Step 6.e: Term 5 of  $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term5_partial_lambda}$$

Term 5 of  $\partial_t \bar{\Lambda}^{i}$: $- 2 \bar{A}^{i j} \left (\partial_{j} \alpha - 6\alpha \partial_{j} \phi\right)$


```python
# Step 6.e: Term 5 of \partial_t \bar{\Lambda}^i:
#           - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi)
for i in range(DIM):
    for j in range(DIM):
        Lambdabar_rhsU[i] += -2*AbarUU[i][j]*(alpha_dD[j] - 6*alpha*phi_dD[j])
```

<a id='term6_partial_lambda'></a>

## Step 6.f: Term 6 of $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term6_partial_lambda}$$

Term 6 of $\partial_t \bar{\Lambda}^{i}$: $2\alpha \bar{A}^{j k} \Delta_{j k}^{i}$


```python
# Step 6.f: Term 6 of \partial_t \bar{\Lambda}^i:
#           2 \alpha \bar{A}^{j k} \Delta^{i}_{j k}
DGammaUDD = Bq.DGammaUDD # From RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            Lambdabar_rhsU[i] += 2*alpha*AbarUU[j][k]*DGammaUDD[i][j][k]
```

<a id='term7_partial_lambda'></a>

## Step 6.g: Term 7 of $\partial_t \bar{\Lambda}^{i}$ \[Back to [top](#toc)\]
$$\label{term7_partial_lambda}$$

$\partial_t \bar{\Lambda}^{i}$: $-\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K$


```python
# Step 6.g: Term 7 of \partial_t \bar{\Lambda}^i:
#           -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
trK_dD = ixp.declarerank1("trK_dD")
for i in range(DIM):
    for j in range(DIM):
        Lambdabar_rhsU[i] += -sp.Rational(4,3)*alpha*gammabarUU[i][j]*trK_dD[j]
```

<a id='rescalingrhss'></a>

# Step 7: Rescaling the BSSN right-hand sides; rewriting them in terms of the rescaled quantities $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$ \[Back to [top](#toc)\]
$$\label{rescalingrhss}$$

Next we rescale the right-hand sides of the BSSN equations so that the evolved variables are $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}\right\}$


```python
# Step 7: Rescale the RHS quantities so that the evolved
#          variables are smooth across coord singularities
h_rhsDD     = ixp.zerorank2()
a_rhsDD     = ixp.zerorank2()
lambda_rhsU = ixp.zerorank1()
for i in range(DIM):
    lambda_rhsU[i] = Lambdabar_rhsU[i] / rfm.ReU[i]
    for j in range(DIM):
        h_rhsDD[i][j] = gammabar_rhsDD[i][j] / rfm.ReDD[i][j]
        a_rhsDD[i][j] =     Abar_rhsDD[i][j] / rfm.ReDD[i][j]
#print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
#print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
#print(betaU_dD)
#print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
#print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
```

<a id='code_validation'></a>

# Step 8: Code Validation against `BSSN.BSSN_RHSs` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_RHSs](../edit/BSSN/BSSN_RHSs.py) module.

By default, we analyze the RHSs in Spherical coordinates, though other coordinate systems may be chosen.


```python
# Step 8: We already have SymPy expressions for BSSN RHS expressions
#         in terms of other SymPy variables. Even if we reset the
#         list of NRPy+ gridfunctions, these *SymPy* expressions for
#         BSSN RHS variables *will remain unaffected*.
#
#         Here, we will use the above-defined BSSN RHS expressions
#         to validate against the same expressions in the
#         BSSN/BSSN_RHSs.py file, to ensure consistency between
#         this tutorial and the module itself.
#
# Reset the list of gridfunctions, as registering a gridfunction
#   twice will spawn an error.
gri.glb_gridfcs_list = []


# Step 9.a: Call the BSSN_RHSs() function from within the
#           BSSN/BSSN_RHSs.py module,
#           which should do exactly the same as in Steps 1-16 above.
import BSSN.BSSN_RHSs as bssnrhs
bssnrhs.BSSN_RHSs()

print("Consistency check between BSSN_RHSs tutorial and NRPy+ module: ALL SHOULD BE ZERO.")

print("trK_rhs - bssnrhs.trK_rhs = " + str(trK_rhs - bssnrhs.trK_rhs))
print("cf_rhs - bssnrhs.cf_rhs = " + str(cf_rhs - bssnrhs.cf_rhs))

for i in range(DIM):
    print("lambda_rhsU["+str(i)+"] - bssnrhs.lambda_rhsU["+str(i)+"] = " +
          str(lambda_rhsU[i] - bssnrhs.lambda_rhsU[i]))
    for j in range(DIM):
        print("h_rhsDD["+str(i)+"]["+str(j)+"] - bssnrhs.h_rhsDD["+str(i)+"]["+str(j)+"] = "
              + str(h_rhsDD[i][j] - bssnrhs.h_rhsDD[i][j]))
        print("a_rhsDD["+str(i)+"]["+str(j)+"] - bssnrhs.a_rhsDD["+str(i)+"]["+str(j)+"] = "
              + str(a_rhsDD[i][j] - bssnrhs.a_rhsDD[i][j]))
```

    Consistency check between BSSN_RHSs tutorial and NRPy+ module: ALL SHOULD BE ZERO.
    trK_rhs - bssnrhs.trK_rhs = 0
    cf_rhs - bssnrhs.cf_rhs = 0
    lambda_rhsU[0] - bssnrhs.lambda_rhsU[0] = 0
    h_rhsDD[0][0] - bssnrhs.h_rhsDD[0][0] = 0
    a_rhsDD[0][0] - bssnrhs.a_rhsDD[0][0] = 0
    h_rhsDD[0][1] - bssnrhs.h_rhsDD[0][1] = 0
    a_rhsDD[0][1] - bssnrhs.a_rhsDD[0][1] = 0
    h_rhsDD[0][2] - bssnrhs.h_rhsDD[0][2] = 0
    a_rhsDD[0][2] - bssnrhs.a_rhsDD[0][2] = 0
    lambda_rhsU[1] - bssnrhs.lambda_rhsU[1] = 0
    h_rhsDD[1][0] - bssnrhs.h_rhsDD[1][0] = 0
    a_rhsDD[1][0] - bssnrhs.a_rhsDD[1][0] = 0
    h_rhsDD[1][1] - bssnrhs.h_rhsDD[1][1] = 0
    a_rhsDD[1][1] - bssnrhs.a_rhsDD[1][1] = 0
    h_rhsDD[1][2] - bssnrhs.h_rhsDD[1][2] = 0
    a_rhsDD[1][2] - bssnrhs.a_rhsDD[1][2] = 0
    lambda_rhsU[2] - bssnrhs.lambda_rhsU[2] = 0
    h_rhsDD[2][0] - bssnrhs.h_rhsDD[2][0] = 0
    a_rhsDD[2][0] - bssnrhs.a_rhsDD[2][0] = 0
    h_rhsDD[2][1] - bssnrhs.h_rhsDD[2][1] = 0
    a_rhsDD[2][1] - bssnrhs.a_rhsDD[2][1] = 0
    h_rhsDD[2][2] - bssnrhs.h_rhsDD[2][2] = 0
    a_rhsDD[2][2] - bssnrhs.a_rhsDD[2][2] = 0


<a id='latex_pdf_output'></a>

# Step 9: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_time_evolution-BSSN_RHSs.pdf](Tutorial-BSSN_time_evolution-BSSN_RHSs.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_time_evolution-BSSN_RHSs")
```

    Created Tutorial-BSSN_time_evolution-BSSN_RHSs.tex, and compiled LaTeX file
        to PDF file Tutorial-BSSN_time_evolution-BSSN_RHSs.pdf

