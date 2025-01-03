<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# BSSN Quantities

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This module documents and constructs a number of quantities useful for building symbolic (SymPy) expressions in terms of the core BSSN quantities $\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$, as defined  in [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658) (see also [Baumgarte, Montero, Cordero-Carrión, and Müller (2012)](https://arxiv.org/abs/1211.6632)).

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

[comment]: <> (Introduction: TODO)

### A Note on Notation:

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

Each family of quantities is constructed within a given function (**boldfaced** below). This notebook is organized as follows


1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
1. [Step 2](#declare_bssn_gfs): **`declare_BSSN_gridfunctions_if_not_declared_already()`**: Declare all of the core BSSN variables $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$ and register them as gridfunctions
1. [Step 3](#rescaling_tensors) Rescaling tensors to avoid coordinate singularities
    1. [Step 3.a](#bssn_basic_tensors) **`BSSN_basic_tensors()`**: Define all basic conformal BSSN tensors $\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\bar{\Lambda}^{i}, \beta^i, B^i\right\}$ in terms of BSSN gridfunctions
1. [Step 4](#bssn_barred_metric__inverse_and_derivs): **`gammabar__inverse_and_derivs()`**: $\bar{\gamma}^{ij}$, and spatial derivatives of $\bar{\gamma}_{ij}$ including $\bar{\Gamma}^{i}_{jk}$
    1. [Step 4.a](#bssn_barred_metric__inverse): Inverse conformal 3-metric: $\bar{\gamma^{ij}}$
    1. [Step 4.b](#bssn_barred_metric__derivs): Derivatives of the conformal 3-metric $\bar{\gamma}_{ij,k}$ and $\bar{\gamma}_{ij,kl}$, and associated "barred" Christoffel symbols $\bar{\Gamma}^{i}_{jk}$
1. [Step 5](#detgammabar_and_derivs): **`detgammabar_and_derivs()`**: $\det \bar{\gamma}_{ij}$ and its derivatives
1. [Step 6](#abar_quantities): **`AbarUU_AbarUD_trAbar()`**: Quantities related to conformal traceless extrinsic curvature $\bar{A}_{ij}$: $\bar{A}^{ij}$, $\bar{A}^i_j$, and $\bar{A}^k_k$
1. [Step 7](#rbar): **`RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()`**: The conformal ("barred") Ricci tensor $\bar{R}_{ij}$ and associated quantities
    1. [Step 7.a](#rbar_part1): Conformal Ricci tensor, part 1: The $\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}$ term
    1. [Step 7.b](#rbar_part2): Conformal Ricci tensor, part 2: The $\bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k}$ term
    1. [Step 7.c](#rbar_part3): Conformal Ricci tensor, part 3: The $\Delta^{k} \Delta_{(i j) k}  + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right )$ terms
    1. [Step 7.d](#summing_rbar_terms): Summing the terms and defining $\bar{R}_{ij}$
1. [Step 8](#beta_derivs): **`betaU_derivs()`**: Unrescaled shift vector $\beta^i$ and spatial derivatives $\beta^i_{,j}$ and $\beta^i_{,jk}$
1. [Step 9](#phi_and_derivs): **`phi_and_derivs()`**: Standard BSSN conformal factor $\phi$, and its derivatives $\phi_{,i}$, $\phi_{,ij}$, $\bar{D}_j \phi$, and $\bar{D}_j\bar{D}_k \phi$
    1. [Step 9.a](#phi_ito_cf): $\phi$ in terms of the chosen (possibly non-standard) conformal factor variable `cf` (e.g., `cf`$=W=e^{-4\phi}$)
    1. [Step 9.b](#phi_covariant_derivs): Partial and covariant derivatives of $\phi$
1. [Step 10](#code_validation): Code Validation against `BSSN.BSSN_quantities` NRPy+ module
1. [Step 11](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules \[Back to [top](#toc)\]

$$\label{initializenrpy}$$


```python
# Step 1: Import all needed modules from NRPy+:
import NRPy_param_funcs as par
import sympy as sp
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import sys

# Step 1.a: Set the coordinate system for the numerical grid
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")

# Step 1.b: Given the chosen coordinate system, set up
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
rfm.reference_metric()

# Step 1.c: Set spatial dimension (must be 3 for BSSN, as BSSN is
#           a 3+1-dimensional decomposition of the general
#           relativistic field equations)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Step 1.d: Declare/initialize parameters for this module
thismodule = "BSSN_quantities"
par.initialize_param(par.glb_param("char", thismodule, "EvolvedConformalFactor_cf", "W"))
par.initialize_param(par.glb_param("bool", thismodule, "detgbarOverdetghat_equals_one", "True"))
```

<a id='declare_bssn_gfs'></a>

# Step 2: `declare_BSSN_gridfunctions_if_not_declared_already()`: Declare all of the core BSSN variables $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$ and register them as gridfunctions \[Back to [top](#toc)\]
$$\label{declare_bssn_gfs}$$


```python
# Step 2: Register all needed BSSN gridfunctions.
# Step 2.a: Register indexed quantities, using ixp.register_... functions
hDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "hDD", "sym01")
aDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "aDD", "sym01")
lambdaU = ixp.register_gridfunctions_for_single_rank1("EVOL", "lambdaU")
vetU = ixp.register_gridfunctions_for_single_rank1("EVOL", "vetU")
betU = ixp.register_gridfunctions_for_single_rank1("EVOL", "betU")
# Step 2.b: Register scalar quantities, using gri.register_gridfunctions()
trK, cf, alpha = gri.register_gridfunctions("EVOL",["trK", "cf", "alpha"])
```

<a id='rescaling_tensors'></a>

# Step 3: Rescaling tensors to avoid coordinate singularities \[Back to [top](#toc)\]
$$\label{rescaling_tensors}$$

While the [covariant form of the BSSN evolution equations](Tutorial-BSSNCurvilinear.ipynb) are properly covariant (with the potential exception of the shift evolution equation, since the shift is a [freely specifiable gauge quantity](https://en.wikipedia.org/wiki/Gauge_fixing)), components of the rank-1 and rank-2 tensors $\varepsilon_{i j}$, $\bar{A}_{i j}$, and $\bar{\Lambda}^{i}$ will drop to zero (destroying information) or diverge (to $\infty$) at coordinate singularities. 

The good news is, this singular behavior is well-understood in terms of the scale factors of the reference metric, enabling us to define rescaled version of these quantities that are well behaved (so that, e.g., they can be finite differenced).

For example, given a smooth vector *in a 3D Cartesian basis* $\bar{\Lambda}^{i}$, all components $\bar{\Lambda}^{x}$, $\bar{\Lambda}^{y}$, and $\bar{\Lambda}^{z}$ will be smooth (by assumption). When changing the basis to spherical coordinates (applying the appropriate Jacobian matrix transformation), we will find that since $\phi = \arctan(y/x)$, $\bar{\Lambda}^{\phi}$ is given by

\begin{align}
\bar{\Lambda}^{\phi} &= \frac{\partial \phi}{\partial x} \bar{\Lambda}^{x} + 
\frac{\partial \phi}{\partial y} \bar{\Lambda}^{y} + 
\frac{\partial \phi}{\partial z} \bar{\Lambda}^{z} \\
&= -\frac{y}{x^2+y^2} \bar{\Lambda}^{x} + 
\frac{x}{x^2+y^2} \bar{\Lambda}^{y} \\
&= -\frac{y}{(r \sin\theta)^2} \bar{\Lambda}^{x} + 
\frac{x}{(r \sin\theta)^2} \bar{\Lambda}^{y}.
\end{align}

Thus $\bar{\Lambda}^{\phi}$ diverges at all points where $r\sin\theta=0$ (or equivalently where $x=y=0$; i.e., the $z$-axis) due to the $\frac{1}{(r\sin\theta)^2}$ that appear in the Jacobian transformation. 

This divergence might pose no problem on cell-centered grids that avoid $r \sin\theta=0$, except that the BSSN equations require that *first and second derivatives* of these quantities be taken. Usual strategies for numerical approximation of these derivatives (e.g., finite difference methods) will "see" these divergences and errors generally will not drop to zero with increased numerical sampling of the functions at points near where the functions diverge.

However, notice that if we define $\lambda^{\phi}$ such that

$$\bar{\Lambda}^{\phi} = \frac{1}{r\sin\theta} \lambda^{\phi},$$

then $\lambda^{\phi}$ will be smooth as well. 

Avoiding such singularities can be generalized to other coordinate systems, so long as $\lambda^i$ is defined as:

$$\bar{\Lambda}^{i} = \frac{\lambda^i}{\text{scalefactor[i]}} ,$$

where scalefactor\[i\] is the $i$th scale factor in the given coordinate system. In an identical fashion, we define the smooth versions of $\beta^i$ and $B^i$ to be $\mathcal{V}^i$ and $\mathcal{B}^i$, respectively. We refer to $\mathcal{V}^i$ and $\mathcal{B}^i$ as vet\[i\] and bet\[i\] respectively in the code after the Hebrew letters that bear some resemblance. 

Similarly, we define the smooth versions of $\bar{A}_{ij}$ and $\varepsilon_{ij}$ ($a_{ij}$ and $h_{ij}$, respectively) via

\begin{align}
\bar{A}_{ij} &= \text{scalefactor[i]}\ \text{scalefactor[j]}\  a_{ij} \\
\varepsilon_{ij} &= \text{scalefactor[i]}\ \text{scalefactor[j]}\  h_{ij},
\end{align}

where in this case we *multiply* due to the fact that these tensors are purely covariant (as opposed to contravariant). To slightly simplify the notation, in NRPy+ we define the *rescaling matrices* `ReU[i]` and `ReDD[i][j]`, such that

\begin{align}
\text{ReU[i]} &= 1 / \text{scalefactor[i]} \\
\text{ReDD[i][j]} &= \text{scalefactor[i] scalefactor[j]}.
\end{align}

Thus, for example, $\bar{A}_{ij}$ and $\bar{\Lambda}^i$ can be expressed as the [Hadamard product](https://en.wikipedia.org/w/index.php?title=Hadamard_product_(matrices)&oldid=852272177) of matrices :

\begin{align}
\bar{A}_{ij} &= \mathbf{ReDD}\circ\mathbf{a} = \text{ReDD[i][j]} a_{ij} \\
\bar{\Lambda}^{i} &= \mathbf{ReU}\circ\mathbf{\lambda} = \text{ReU[i]} \lambda^i,
\end{align}
where no sums are implied by the repeated indices.

Further, since the scale factors are *time-independent*, 

\begin{align}
\partial_t \bar{A}_{ij} &= \text{ReDD[i][j]}\  \partial_t a_{ij} \\
\partial_t \bar{\gamma}_{ij} &= \partial_t \left(\varepsilon_{ij} + \hat{\gamma}_{ij}\right)\\
&= \partial_t \varepsilon_{ij} \\
&= \text{scalefactor[i]}\ \text{scalefactor[j]}\ \partial_t h_{ij}.
\end{align}

Thus instead of taking space or time derivatives of BSSN quantities

$$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$ 

across coordinate singularities, we instead factor out the singular scale factors according to this prescription so that space or time derivatives of BSSN quantities are written in terms of finite-difference derivatives of the *rescaled* variables 

$$\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\},$$ 

and *exact* expressions for (spatial) derivatives of scale factors. Note that `cf` is the chosen conformal factor (supported choices for `cf` are discussed in [Step 6.a](#phi_ito_cf)). 

As an example, let's evaluate $\bar{\Lambda}^{i}_{\, ,\, j}$ according to this prescription:

\begin{align}
\bar{\Lambda}^{i}_{\, ,\, j} &= -\frac{\lambda^i}{(\text{ReU[i]})^2}\ \partial_j \left(\text{ReU[i]}\right) + \frac{\partial_j \lambda^i}{\text{ReU[i]}} \\
&= -\frac{\lambda^i}{(\text{ReU[i]})^2}\ \text{ReUdD[i][j]} + \frac{\partial_j \lambda^i}{\text{ReU[i]}}.
\end{align}

Here, the derivative `ReUdD[i][j]` **is computed symbolically and exactly** using SymPy, and the derivative $\partial_j \lambda^i$ represents a derivative of a *smooth* quantity (so long as $\bar{\Lambda}^{i}$ is smooth in the Cartesian basis).

<a id='bssn_basic_tensors'></a>

## Step 3.a: `BSSN_basic_tensors()`: Define all basic conformal BSSN tensors $\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\bar{\Lambda}^{i}, \beta^i, B^i\right\}$ in terms of BSSN gridfunctions \[Back to [top](#toc)\]
$$\label{bssn_basic_tensors}$$

The `BSSN_vars__tensors()` function defines the tensorial BSSN quantities $\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\bar{\Lambda}^{i}, \beta^i, B^i\right\}$, in terms of the rescaled "base" tensorial quantities $\left\{h_{i j},a_{i j}, \lambda^{i}, \mathcal{V}^i, \mathcal{B}^i\right\},$ respectively:

\begin{align}
\bar{\gamma}_{i j} &= \hat{\gamma}_{ij} + \varepsilon_{ij}, \text{ where } \varepsilon_{ij} = h_{ij} \circ \text{ReDD[i][j]} \\
\bar{A}_{i j} &= a_{ij} \circ \text{ReDD[i][j]} \\
\bar{\Lambda}^{i} &= \lambda^i \circ \text{ReU[i]} \\
\beta^{i} &= \mathcal{V}^i \circ \text{ReU[i]} \\
B^{i} &= \mathcal{B}^i \circ \text{ReU[i]}
\end{align}

Rescaling vectors and tensors are built upon the scale factors for the chosen (in a general, singular) coordinate system, which is defined in NRPy+'s [reference_metric.py](../edit/reference_metric.py) ([Tutorial](Tutorial-Reference_Metric.ipynb)), and the rescaled variables are defined in the stub function [BSSN/BSSN_rescaled_vars.py](../edit/BSSN/BSSN_rescaled_vars.py). 

Here we implement `BSSN_vars__tensors()`:


```python
# Step 3.a: Define all basic conformal BSSN tensors in terms of BSSN gridfunctions

# Step 3.a.i: gammabarDD and AbarDD:
gammabarDD = ixp.zerorank2()
AbarDD     = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
        gammabarDD[i][j] = hDD[i][j]*rfm.ReDD[i][j] + rfm.ghatDD[i][j]
        # Abar_{ij}      = a_{ij}*ReDD[i][j]
        AbarDD[i][j]     = aDD[i][j]*rfm.ReDD[i][j]

# Step 3.a.ii: LambdabarU, betaU, and BU:
LambdabarU = ixp.zerorank1()
betaU      = ixp.zerorank1()
BU         = ixp.zerorank1()
for i in range(DIM):
    LambdabarU[i] = lambdaU[i]*rfm.ReU[i]
    betaU[i]      = vetU[i]   *rfm.ReU[i]
    BU[i]         = betU[i]   *rfm.ReU[i]
```

<a id='bssn_barred_metric__inverse_and_derivs'></a>

# Step 4: `gammabar__inverse_and_derivs()`: $\bar{\gamma}^{ij}$, and spatial derivatives of $\bar{\gamma}_{ij}$ including $\bar{\Gamma}^{i}_{jk}$ \[Back to [top](#toc)\]
$$\label{bssn_barred_metric__inverse_and_derivs}$$

<a id='bssn_barred_metric__inverse'></a>

## Step 4.a: Inverse conformal 3-metric: $\bar{\gamma^{ij}}$ \[Back to [top](#toc)\]
$$\label{bssn_barred_metric__inverse}$$

Since $\bar{\gamma}^{ij}$ is the inverse of $\bar{\gamma}_{ij}$, we apply a $3\times 3$ symmetric matrix inversion to compute $\bar{\gamma}^{ij}$. 


```python
# Step 4.a: Inverse conformal 3-metric gammabarUU:
# Step 4.a.i: gammabarUU:
gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(gammabarDD)
```

<a id='bssn_barred_metric__derivs'></a>

## Step 4.b: Derivatives of the conformal 3-metric $\bar{\gamma}_{ij,k}$ and $\bar{\gamma}_{ij,kl}$, and associated "barred" Christoffel symbols $\bar{\Gamma}^{i}_{jk}$ \[Back to [top](#toc)\]
$$\label{bssn_barred_metric__derivs}$$

In the BSSN-in-curvilinear coordinates formulation, all quantities must be defined in terms of rescaled quantities $h_{ij}$ and their derivatives (evaluated using finite differences), as well as reference-metric quantities and their derivatives (evaluated exactly using SymPy). 

For example, $\bar{\gamma}_{ij,k}$ is given by:
\begin{align}
\bar{\gamma}_{ij,k} &= \partial_k \bar{\gamma}_{ij} \\
&= \partial_k \left(\hat{\gamma}_{ij} + \varepsilon_{ij}\right) \\
&= \partial_k \left(\hat{\gamma}_{ij} + h_{ij} \text{ReDD[i][j]}\right) \\
&= \hat{\gamma}_{ij,k} + h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]},
\end{align}
where `ReDDdD[i][j][k]` is computed within `rfm.reference_metric()`.


```python
# Step 4.b.i gammabarDDdD[i][j][k]
#            = \hat{\gamma}_{ij,k} + h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]}.

gammabarDD_dD = ixp.zerorank3()
hDD_dD = ixp.declarerank3("hDD_dD","sym01")
hDD_dupD = ixp.declarerank3("hDD_dupD","sym01")
gammabarDD_dupD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            gammabarDD_dD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                     hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k]

            # Compute associated upwinded derivative, needed for the \bar{\gamma}_{ij} RHS
            gammabarDD_dupD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                       hDD_dupD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k]
```

By extension, the second derivative $\bar{\gamma}_{ij,kl}$ is given by
\begin{align}
\bar{\gamma}_{ij,kl} &= \partial_l \left(\hat{\gamma}_{ij,k} + h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]}\right)\\
&= \hat{\gamma}_{ij,kl} + h_{ij,kl} \text{ReDD[i][j]} + h_{ij,k} \text{ReDDdD[i][j][l]} + h_{ij,l} \text{ReDDdD[i][j][k]} + h_{ij} \text{ReDDdDD[i][j][k][l]}
\end{align}


```python
# Step 4.b.ii: Compute gammabarDD_dDD in terms of the rescaled BSSN quantity hDD
#      and its derivatives, as well as the reference metric and rescaling
#      matrix, and its derivatives (expression given below):
hDD_dDD = ixp.declarerank4("hDD_dDD","sym01_sym23")
gammabarDD_dDD = ixp.zerorank4()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                # gammabar_{ij,kl} = gammahat_{ij,kl}
                #                  + h_{ij,kl} ReDD[i][j]
                #                  + h_{ij,k} ReDDdD[i][j][l] + h_{ij,l} ReDDdD[i][j][k]
                #                  + h_{ij} ReDDdDD[i][j][k][l]
                gammabarDD_dDD[i][j][k][l]  = rfm.ghatDDdDD[i][j][k][l]
                gammabarDD_dDD[i][j][k][l] += hDD_dDD[i][j][k][l]*rfm.ReDD[i][j]
                gammabarDD_dDD[i][j][k][l] += hDD_dD[i][j][k]*rfm.ReDDdD[i][j][l] + \
                                              hDD_dD[i][j][l]*rfm.ReDDdD[i][j][k]
                gammabarDD_dDD[i][j][k][l] += hDD[i][j]*rfm.ReDDdDD[i][j][k][l]

```

Finally, we compute the Christoffel symbol associated with the barred 3-metric: $\bar{\Gamma}^{i}_{kl}$:
$$
\bar{\Gamma}^{i}_{kl} = \frac{1}{2} \bar{\gamma}^{im} \left(\bar{\gamma}_{mk,l} + \bar{\gamma}_{ml,k} - \bar{\gamma}_{kl,m} \right)
$$


```python
# Step 4.b.iii: Define barred Christoffel symbol \bar{\Gamma}^{i}_{kl} = GammabarUDD[i][k][l] (see expression below)
GammabarUDD = ixp.zerorank3()
for i in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            for m in range(DIM):
                # Gammabar^i_{kl} = 1/2 * gammabar^{im} ( gammabar_{mk,l} + gammabar_{ml,k} - gammabar_{kl,m}):
                GammabarUDD[i][k][l] += sp.Rational(1,2)*gammabarUU[i][m]* \
                                        (gammabarDD_dD[m][k][l] + gammabarDD_dD[m][l][k] - gammabarDD_dD[k][l][m])
```

<a id='detgammabar_and_derivs'></a>

# Step 5: `detgammabar_and_derivs()`: $\det \bar{\gamma}_{ij}$ and its derivatives \[Back to [top](#toc)\]
$$\label{detgammabar_and_derivs}$$


As described just before Section III of [Baumgarte *et al* (2012)](https://arxiv.org/pdf/1211.6632.pdf), we are free to choose $\det \bar{\gamma}_{ij}$, which should remain fixed in time.

As in [Baumgarte *et al* (2012)](https://arxiv.org/pdf/1211.6632.pdf) generally we make the choice $\det \bar{\gamma}_{ij} = \det \hat{\gamma}_{ij}$, but *this need not be the case; we could choose to set $\det \bar{\gamma}_{ij}$ to another expression.*

In case we do not choose to set $\det \bar{\gamma}_{ij}/\det \hat{\gamma}_{ij}=1$, below we begin the implementation of a gridfunction, `detgbarOverdetghat`, which defines an alternative expression in its place. 

$\det \bar{\gamma}_{ij}/\det \hat{\gamma}_{ij}$=`detgbarOverdetghat`$\ne 1$ is not yet implemented. However, we can define `detgammabar` and its derivatives in terms of a generic `detgbarOverdetghat` and $\det \hat{\gamma}_{ij}$ and their derivatives:

\begin{align}
\text{detgammabar} &= \det \bar{\gamma}_{ij} = \text{detgbarOverdetghat} \cdot \left(\det \hat{\gamma}_{ij}\right) \\
\text{detgammabar}\_\text{dD[k]} &= \left(\det \bar{\gamma}_{ij}\right)_{,k} = \text{detgbarOverdetghat}\_\text{dD[k]} \det \hat{\gamma}_{ij} +  \text{detgbarOverdetghat} \left(\det \hat{\gamma}_{ij}\right)_{,k} \\
\end{align}
https://en.wikipedia.org/wiki/Determinant#Properties_of_the_determinant


```python
# Step 5: det(gammabarDD) and its derivatives
detgbarOverdetghat = sp.sympify(1)
detgbarOverdetghat_dD = ixp.zerorank1()
detgbarOverdetghat_dDD = ixp.zerorank2()

if par.parval_from_str(thismodule+"::detgbarOverdetghat_equals_one") == "False":
    print("Error: detgbarOverdetghat_equals_one=\"False\" is not fully implemented yet.")
    sys.exit(1)
## Approach for implementing detgbarOverdetghat_equals_one=False:
#     detgbarOverdetghat = gri.register_gridfunctions("AUX", ["detgbarOverdetghat"])
#     detgbarOverdetghatInitial = gri.register_gridfunctions("AUX", ["detgbarOverdetghatInitial"])
#     detgbarOverdetghat_dD = ixp.declarerank1("detgbarOverdetghat_dD")
#     detgbarOverdetghat_dDD = ixp.declarerank2("detgbarOverdetghat_dDD", "sym01")

# Step 5.b: Define detgammabar, detgammabar_dD, and detgammabar_dDD (needed for
#           \partial_t \bar{\Lambda}^i below)detgammabar = detgbarOverdetghat * rfm.detgammahat
detgammabar = detgbarOverdetghat * rfm.detgammahat

detgammabar_dD = ixp.zerorank1()
for i in range(DIM):
    detgammabar_dD[i] = detgbarOverdetghat_dD[i] * rfm.detgammahat + detgbarOverdetghat * rfm.detgammahatdD[i]

detgammabar_dDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        detgammabar_dDD[i][j] = detgbarOverdetghat_dDD[i][j] * rfm.detgammahat + \
                                detgbarOverdetghat_dD[i] * rfm.detgammahatdD[j] + \
                                detgbarOverdetghat_dD[j] * rfm.detgammahatdD[i] + \
                                detgbarOverdetghat * rfm.detgammahatdDD[i][j]
```

<a id='abar_quantities'></a>

# Step 6: `AbarUU_AbarUD_trAbar_AbarDD_dD()`: Quantities related to conformal traceless extrinsic curvature $\bar{A}_{ij}$: $\bar{A}^{ij}$, $\bar{A}^i_j$, and $\bar{A}^k_k$ \[Back to [top](#toc)\]
$$\label{abar_quantities}$$

$\bar{A}^{ij}$ is given by application of the raising operators (a.k.a., the inverse 3-metric) $\bar{\gamma}^{jk}$ on both of the covariant ("down") components:
$$
\bar{A}^{ij} = \bar{\gamma}^{ik}\bar{\gamma}^{jl} \bar{A}_{kl}.
$$

$\bar{A}^i_j$ is given by a single application of the raising operator (a.k.a., the inverse 3-metric) $\bar{\gamma}^{ik}$ on $\bar{A}_{kj}$:
$$
\bar{A}^i_j = \bar{\gamma}^{ik}\bar{A}_{kj}.
$$

The trace of $\bar{A}_{ij}$, $\bar{A}^k_k$, is given by a contraction with the barred 3-metric:
$$
\text{Tr}(\bar{A}_{ij}) = \bar{A}^k_k = \bar{\gamma}^{kj}\bar{A}_{jk}.
$$

Note that while $\bar{A}_{ij}$ is defined as the *traceless* conformal extrinsic curvature, it may acquire a nonzero trace (assuming the initial data impose tracelessness) due to numerical error. $\text{Tr}(\bar{A}_{ij})$ is included in the BSSN equations to drive $\text{Tr}(\bar{A}_{ij})$ to zero.

In terms of rescaled BSSN quantities, $\bar{A}_{ij}$ is given by
$$
\bar{A}_{ij} = \text{ReDD[i][j]} a_{ij},
$$
so in terms of the same quantities, $\bar{A}_{ij,k}$ is given by
$$
\bar{A}_{ij,k} = \text{ReDDdD[i][j][k]} a_{ij} + \text{ReDD[i][j]} a_{ij,k}.
$$


```python
# Step 6: Quantities related to conformal traceless extrinsic curvature

# Step 6.a.i: Compute Abar^{ij} in terms of Abar_{ij} and gammabar^{ij}
AbarUU = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                # Abar^{ij} = gammabar^{ik} gammabar^{jl} Abar_{kl}
                AbarUU[i][j] += gammabarUU[i][k]*gammabarUU[j][l]*AbarDD[k][l]

# Step 6.a.ii: Compute Abar^i_j in terms of Abar_{ij} and gammabar^{ij}
AbarUD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Abar^i_j = gammabar^{ik} Abar_{kj}
            AbarUD[i][j] += gammabarUU[i][k]*AbarDD[k][j]

# Step 6.a.iii: Compute Abar^k_k = trace of Abar:
trAbar = sp.sympify(0)
for k in range(DIM):
    for j in range(DIM):
        # Abar^k_k = gammabar^{kj} Abar_{jk}
        trAbar += gammabarUU[k][j]*AbarDD[j][k]

# Step 6.a.iv: Compute Abar_{ij,k}
AbarDD_dD   = ixp.zerorank3()
AbarDD_dupD = ixp.zerorank3()
aDD_dD   = ixp.declarerank3("aDD_dD"  ,"sym01")
aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            AbarDD_dupD[i][j][k] = rfm.ReDDdD[i][j][k]*aDD[i][j] + rfm.ReDD[i][j]*aDD_dupD[i][j][k]
            AbarDD_dD[i][j][k]   = rfm.ReDDdD[i][j][k]*aDD[i][j] + rfm.ReDD[i][j]*aDD_dD[  i][j][k]
```

<a id='rbar'></a>

# Step 7: `RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()`: The conformal ("barred") Ricci tensor $\bar{R}_{ij}$ and associated quantities \[Back to [top](#toc)\]
$$\label{rbar}$$

Let's compute perhaps the most complicated expression in the BSSN evolution equations, the conformal Ricci tensor:

\begin{align}
  \bar{R}_{i j} {} = {} & - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} + \bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k} + \Delta^{k} \Delta_{(i j) k} \nonumber \\
  & + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right ) \; .
\end{align}

Let's tackle the $\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}$ term first:

<a id='rbar_part1'></a>

## Step 7.a: Conformal Ricci tensor, part 1: The $\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}$ term \[Back to [top](#toc)\]
$$\label{rbar_part1}$$

First note that the covariant derivative of a metric with respect to itself is zero
$$\hat{D}_{l} \hat{\gamma}_{ij} = 0,$$
so 
$$\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} = \hat{D}_{k} \hat{D}_{l} \left(\hat{\gamma}_{i j} + \varepsilon_{ij}\right) = \hat{D}_{k} \hat{D}_{l} \varepsilon_{ij}.$$

Next, the covariant derivative of a tensor is given by (from the [wikipedia article on covariant differentiation](https://en.wikipedia.org/wiki/Covariant_derivative)):
\begin{align}
  {(\nabla_{e_c} T)^{a_1 \ldots a_r}}_{b_1 \ldots b_s} = {}
    &\frac{\partial}{\partial x^c}{T^{a_1 \ldots a_r}}_{b_1 \ldots b_s} \\
    &+ \,{\Gamma ^{a_1}}_{dc} {T^{d a_2 \ldots a_r}}_{b_1 \ldots b_s} + \cdots + {\Gamma^{a_r}}_{dc} {T^{a_1 \ldots a_{r-1}d}}_{b_1 \ldots b_s} \\
    &-\,{\Gamma^d}_{b_1 c} {T^{a_1 \ldots a_r}}_{d b_2 \ldots b_s} - \cdots - {\Gamma^d}_{b_s c} {T^{a_1 \ldots a_r}}_{b_1 \ldots b_{s-1} d}.
\end{align}

Therefore, 
$$\hat{D}_{l} \bar{\gamma}_{i j} = \hat{D}_{l} \varepsilon_{i j} = \varepsilon_{i j,l} - \hat{\Gamma}^m_{i l} \varepsilon_{m j} -\hat{\Gamma}^m_{j l} \varepsilon_{i m}.$$

Since the covariant first derivative is a tensor, the covariant second derivative is given by (same as [Eq. 27 in Baumgarte et al (2012)](https://arxiv.org/pdf/1211.6632.pdf))

\begin{align}
\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} &= \hat{D}_{k} \hat{D}_{l} \varepsilon_{i j}  \\
&= \partial_k \hat{D}_{l} \varepsilon_{i j}
 - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right) 
 - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right)
 - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right),
\end{align}

where the first term is the partial derivative of the expression already derived for $\hat{D}_{l} \varepsilon_{i j}$:

\begin{align}
\partial_k \hat{D}_{l} \varepsilon_{i j} &= \partial_k \left(\varepsilon_{ij,l} - \hat{\Gamma}^m_{i l} \varepsilon_{m j} -\hat{\Gamma}^m_{j l} \varepsilon_{i m} \right) \\
&= \varepsilon_{ij,lk} - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j} - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k} - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m} - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}.
\end{align}

In terms of the evolved quantity $h_{ij}$, the derivatives of $\varepsilon_{ij}$ are given by:
\begin{align}
\varepsilon_{ij,k} &= \partial_k \left(h_{ij} \text{ReDD[i][j]}\right) \\
&= h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]},
\end{align}
and
\begin{align}
\varepsilon_{ij,kl} &= \partial_l \left(h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]} \right)\\
&= h_{ij,kl} \text{ReDD[i][j]} + h_{ij,k} \text{ReDDdD[i][j][l]} + h_{ij,l} \text{ReDDdD[i][j][k]} + h_{ij} \text{ReDDdDD[i][j][k][l]}.
\end{align}


```python
# Step 7: Conformal Ricci tensor, part 1: The \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} term

# Step 7.a.i: Define \varepsilon_{ij} = epsDD[i][j]
epsDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        epsDD[i][j] = hDD[i][j]*rfm.ReDD[i][j]

# Step 7.a.ii: Define epsDD_dD[i][j][k]
hDD_dD = ixp.declarerank3("hDD_dD","sym01")
epsDD_dD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            epsDD_dD[i][j][k] = hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k]

# Step 7.a.iii: Define epsDD_dDD[i][j][k][l]
hDD_dDD = ixp.declarerank4("hDD_dDD","sym01_sym23")
epsDD_dDD = ixp.zerorank4()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                epsDD_dDD[i][j][k][l] = hDD_dDD[i][j][k][l]*rfm.ReDD[i][j] + \
                                        hDD_dD[i][j][k]*rfm.ReDDdD[i][j][l] + \
                                        hDD_dD[i][j][l]*rfm.ReDDdD[i][j][k] + \
                                        hDD[i][j]*rfm.ReDDdDD[i][j][k][l]
```

We next compute three quantities derived above:

* `gammabarDD_DhatD[i][j][l]` = $\hat{D}_{l} \bar{\gamma}_{i j} = \hat{D}_{l} \varepsilon_{i j} = \varepsilon_{i j,l} - \hat{\Gamma}^m_{i l} \varepsilon_{m j} -\hat{\Gamma}^m_{j l} \varepsilon_{i m}$,
* `gammabarDD_DhatD\_dD[i][j][l][k]` = $\partial_k \hat{D}_{l} \bar{\gamma}_{i j} = \partial_k \hat{D}_{l} \varepsilon_{i j} = \varepsilon_{ij,lk} - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j} - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k} - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m} - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}$, and
* `gammabarDD_DhatDD[i][j][l][k]` = $\hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} = \partial_k \hat{D}_{l} \varepsilon_{i j} - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right)  - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right) - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right)$.


```python
# Step 7.a.iv: DhatgammabarDDdD[i][j][l] = \bar{\gamma}_{ij;\hat{l}}
# \bar{\gamma}_{ij;\hat{l}} = \varepsilon_{i j,l}
#                           - \hat{\Gamma}^m_{i l} \varepsilon_{m j}
#                           - \hat{\Gamma}^m_{j l} \varepsilon_{i m}
gammabarDD_dHatD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for l in range(DIM):
            gammabarDD_dHatD[i][j][l] = epsDD_dD[i][j][l]
            for m in range(DIM):
                gammabarDD_dHatD[i][j][l] += - rfm.GammahatUDD[m][i][l]*epsDD[m][j] \
                                             - rfm.GammahatUDD[m][j][l]*epsDD[i][m]

# Step 7.a.v: \bar{\gamma}_{ij;\hat{l},k} = DhatgammabarDD_dHatD_dD[i][j][l][k]:
#        \bar{\gamma}_{ij;\hat{l},k} = \varepsilon_{ij,lk}
#                                      - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j}
#                                      - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k}
#                                      - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m}
#                                      - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}
gammabarDD_dHatD_dD = ixp.zerorank4()
for i in range(DIM):
    for j in range(DIM):
        for l in range(DIM):
            for k in range(DIM):
                gammabarDD_dHatD_dD[i][j][l][k] = epsDD_dDD[i][j][l][k]
                for m in range(DIM):
                    gammabarDD_dHatD_dD[i][j][l][k] += -rfm.GammahatUDDdD[m][i][l][k]*epsDD[m][j]  \
                                                       -rfm.GammahatUDD[m][i][l]*epsDD_dD[m][j][k] \
                                                       -rfm.GammahatUDDdD[m][j][l][k]*epsDD[i][m]  \
                                                       -rfm.GammahatUDD[m][j][l]*epsDD_dD[i][m][k]

# Step 7.a.vi: \bar{\gamma}_{ij;\hat{l}\hat{k}} = DhatgammabarDD_dHatDD[i][j][l][k]
#          \bar{\gamma}_{ij;\hat{l}\hat{k}} = \partial_k \hat{D}_{l} \varepsilon_{i j}
#                                           - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right)
#                                           - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right)
#                                           - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right)
gammabarDD_dHatDD = ixp.zerorank4()
for i in range(DIM):
    for j in range(DIM):
        for l in range(DIM):
            for k in range(DIM):
                gammabarDD_dHatDD[i][j][l][k] = gammabarDD_dHatD_dD[i][j][l][k]
                for m in range(DIM):
                    gammabarDD_dHatDD[i][j][l][k] += - rfm.GammahatUDD[m][l][k]*gammabarDD_dHatD[i][j][m] \
                                                     - rfm.GammahatUDD[m][i][k]*gammabarDD_dHatD[m][j][l] \
                                                     - rfm.GammahatUDD[m][j][k]*gammabarDD_dHatD[i][m][l]
```

<a id='rbar_part2'></a>

## Step 7.b: Conformal Ricci tensor, part 2: The $\bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k}$ term \[Back to [top](#toc)\]
$$\label{rbar_part2}$$

By definition, the index symmetrization operation is given by:
$$\bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k} = \frac{1}{2} \left( \bar{\gamma}_{ki} \hat{D}_{j} \bar{\Lambda}^{k} + \bar{\gamma}_{kj} \hat{D}_{i} \bar{\Lambda}^{k} \right),$$

and $\bar{\gamma}_{ij}$ is trivially computed ($=\varepsilon_{ij} + \hat{\gamma}_{ij}$) so the only nontrival part to computing this term is in evaluating $\hat{D}_{j} \bar{\Lambda}^{k}$.

The covariant derivative is with respect to the hatted metric (i.e. the reference metric), so
$$\hat{D}_{j} \bar{\Lambda}^{k} = \partial_j \bar{\Lambda}^{k} + \hat{\Gamma}^{k}_{mj} \bar{\Lambda}^m,$$
except we cannot take derivatives of $\bar{\Lambda}^{k}$ directly due to potential issues with coordinate singularities. Instead we write it in terms of the rescaled quantity $\lambda^k$ via
$$\bar{\Lambda}^{k} = \lambda^k \text{ReU[k]}.$$

Then the expression for $\hat{D}_{j} \bar{\Lambda}^{k}$ becomes
$$
\hat{D}_{j} \bar{\Lambda}^{k} = \lambda^{k}_{,j} \text{ReU[k]} + \lambda^{k} \text{ReUdD[k][j]} + \hat{\Gamma}^{k}_{mj} \lambda^{m} \text{ReU[m]},
$$
and the NRPy+ code for this expression is written


```python
# Step 7.b: Second term of RhatDD: compute \hat{D}_{j} \bar{\Lambda}^{k} = LambarU_dHatD[k][j]
lambdaU_dD = ixp.declarerank2("lambdaU_dD","nosym")
LambarU_dHatD = ixp.zerorank2()
for j in range(DIM):
    for k in range(DIM):
        LambarU_dHatD[k][j] = lambdaU_dD[k][j]*rfm.ReU[k] + lambdaU[k]*rfm.ReUdD[k][j]
        for m in range(DIM):
            LambarU_dHatD[k][j] += rfm.GammahatUDD[k][m][j]*lambdaU[m]*rfm.ReU[m]
```

<a id='rbar_part3'></a>

## Step 7.c: Conformal Ricci tensor, part 3: The $\Delta^{k} \Delta_{(i j) k}  + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right )$ terms \[Back to [top](#toc)\]
$$\label{rbar_part3}$$

Our goal here is to compute the quantities appearing as the final terms of the conformal Ricci tensor:
$$
\Delta^{k} \Delta_{(i j) k}  + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right).
$$

* `DGammaUDD[k][i][j]`$= \Delta^k_{ij}$ is simply the difference in Christoffel symbols: $\Delta^{k}_{ij} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}$, and 
* `DGammaU[k]`$= \Delta^k$ is the contraction: $\bar{\gamma}^{ij} \Delta^{k}_{ij}$

Adding these expressions to Ricci is straightforward, since $\bar{\Gamma}^i_{jk}$ and $\bar{\gamma}^{ij}$ were defined above in [Step 4](#bssn_barred_metric__inverse_and_derivs), and $\hat{\Gamma}^i_{jk}$ was computed within NRPy+'s `reference_metric()` function:


```python
# Step 7.c: Conformal Ricci tensor, part 3: The \Delta^{k} \Delta_{(i j) k}
#           + \bar{\gamma}^{k l}*(2 \Delta_{k(i}^{m} \Delta_{j) m l}
#           + \Delta_{i k}^{m} \Delta_{m j l}) terms

# Step 7.c.i: Define \Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk} = DGammaUDD[i][j][k]
DGammaUDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            DGammaUDD[i][j][k] = GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]

# Step 7.c.ii: Define \Delta^i = \bar{\gamma}^{jk} \Delta^i_{jk}
DGammaU = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            DGammaU[i] += gammabarUU[j][k] * DGammaUDD[i][j][k]
```

Next we define $\Delta_{ijk}=\bar{\gamma}_{im}\Delta^m_{jk}$:


```python
# Step 7.c.iii: Define \Delta_{ijk} = \bar{\gamma}_{im} \Delta^m_{jk}
DGammaDDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for m in range(DIM):
                DGammaDDD[i][j][k] += gammabarDD[i][m] * DGammaUDD[m][j][k]
```

<a id='summing_rbar_terms'></a>

## Step 7.d: Summing the terms and defining $\bar{R}_{ij}$ \[Back to [top](#toc)\]
$$\label{summing_rbar_terms}$$

We have now constructed all of the terms going into $\bar{R}_{ij}$:

\begin{align}
  \bar{R}_{i j} {} = {} & - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j} + \bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^{k} + \Delta^{k} \Delta_{(i j) k} \nonumber \\
  & + \bar{\gamma}^{k l} \left (2 \Delta_{k(i}^{m} \Delta_{j) m l} + \Delta_{i k}^{m} \Delta_{m j l} \right ) \; .
\end{align}


```python
# Step 7.d: Summing the terms and defining \bar{R}_{ij}

# Step 7.d.i: Add the first term to RbarDD:
#         Rbar_{ij} += - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}
RbarDD = ixp.zerorank2()
RbarDDpiece = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                RbarDD[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]
                RbarDDpiece[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]

# Step 7.d.ii: Add the second term to RbarDD:
#         Rbar_{ij} += (1/2) * (gammabar_{ki} Lambar^k_{;\hat{j}} + gammabar_{kj} Lambar^k_{;\hat{i}})
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            RbarDD[i][j] += sp.Rational(1,2) * (gammabarDD[k][i]*LambarU_dHatD[k][j] + \
                                                gammabarDD[k][j]*LambarU_dHatD[k][i])

# Step 7.d.iii: Add the remaining term to RbarDD:
#      Rbar_{ij} += \Delta^{k} \Delta_{(i j) k} = 1/2 \Delta^{k} (\Delta_{i j k} + \Delta_{j i k})
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            RbarDD[i][j] += sp.Rational(1,2) * DGammaU[k] * (DGammaDDD[i][j][k] + DGammaDDD[j][i][k])

# Step 7.d.iv: Add the final term to RbarDD:
#      Rbar_{ij} += \bar{\gamma}^{k l} (\Delta^{m}_{k i} \Delta_{j m l}
#                   + \Delta^{m}_{k j} \Delta_{i m l}
#                   + \Delta^{m}_{i k} \Delta_{m j l})
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    RbarDD[i][j] += gammabarUU[k][l] * (DGammaUDD[m][k][i]*DGammaDDD[j][m][l] +
                                                        DGammaUDD[m][k][j]*DGammaDDD[i][m][l] +
                                                        DGammaUDD[m][i][k]*DGammaDDD[m][j][l])
```

<a id='beta_derivs'></a>

# Step 8: **`betaU_derivs()`**: The unrescaled shift vector $\beta^i$ spatial derivatives: $\beta^i_{,j}$ & $\beta^i_{,jk}$, written in terms of the rescaled shift vector $\mathcal{V}^i$ \[Back to [top](#toc)\]
$$\label{beta_derivs}$$

This step, which documents the function `betaUbar_and_derivs()` inside the [BSSN.BSSN_unrescaled_and_barred_vars](../edit/BSSN/BSSN_unrescaled_and_barred_vars) module, defines three quantities:

[comment]: <> (Fix Link Above: TODO)

* `betaU_dD[i][j]`$=\beta^i_{,j} = \left(\mathcal{V}^i \circ \text{ReU[i]}\right)_{,j} = \mathcal{V}^i_{,j} \circ \text{ReU[i]} + \mathcal{V}^i \circ \text{ReUdD[i][j]}$
* `betaU_dupD[i][j]`: the same as above, except using *upwinded* finite-difference derivatives to compute $\mathcal{V}^i_{,j}$ instead of *centered* finite-difference derivatives.
* `betaU_dDD[i][j][k]`$=\beta^i_{,jk} = \mathcal{V}^i_{,jk} \circ \text{ReU[i]} + \mathcal{V}^i_{,j} \circ \text{ReUdD[i][k]} + \mathcal{V}^i_{,k} \circ \text{ReUdD[i][j]}+\mathcal{V}^i \circ \text{ReUdDD[i][j][k]}$


```python
# Step 8: The unrescaled shift vector betaU spatial derivatives:
#         betaUdD & betaUdDD, written in terms of the
#         rescaled shift vector vetU
vetU_dD = ixp.declarerank2("vetU_dD","nosym")
vetU_dupD = ixp.declarerank2("vetU_dupD","nosym") # Needed for upwinded \beta^i_{,j}
vetU_dDD = ixp.declarerank3("vetU_dDD","sym12")   # Needed for \beta^i_{,j}
betaU_dD = ixp.zerorank2()
betaU_dupD = ixp.zerorank2() # Needed for, e.g., \beta^i RHS
betaU_dDD = ixp.zerorank3() # Needed for, e.g., \bar{\Lambda}^i RHS
for i in range(DIM):
    for j in range(DIM):
        betaU_dD[i][j] = vetU_dD[i][j]*rfm.ReU[i] + vetU[i]*rfm.ReUdD[i][j]
        betaU_dupD[i][j] = vetU_dupD[i][j]*rfm.ReU[i] + vetU[i]*rfm.ReUdD[i][j] # Needed for \beta^i RHS
        for k in range(DIM):
            # Needed for, e.g., \bar{\Lambda}^i RHS:
            betaU_dDD[i][j][k] = vetU_dDD[i][j][k]*rfm.ReU[i]  + vetU_dD[i][j]*rfm.ReUdD[i][k] + \
                                 vetU_dD[i][k]*rfm.ReUdD[i][j] + vetU[i]*rfm.ReUdDD[i][j][k]
```

<a id='phi_and_derivs'></a>

# Step 9: **`phi_and_derivs()`**: Standard BSSN conformal factor $\phi$, and its derivatives $\phi_{,i}$, $\phi_{,ij}$, $\bar{D}_j \phi$, and $\bar{D}_j\bar{D}_k \phi$, all written in terms of BSSN gridfunctions like $\text{cf}$ \[Back to [top](#toc)\]
$$\label{phi_and_derivs}$$

<a id='phi_ito_cf'></a>

## Step 9.a: $\phi$ in terms of the chosen (possibly non-standard) conformal factor variable $\text{cf}$ (e.g., $\text{cf}=\chi=e^{-4\phi}$) \[Back to [top](#toc)\]
$$\label{phi_ito_cf}$$

When solving the BSSN time evolution equations across the coordinate singularity (i.e., the "puncture") inside puncture black holes for example, the standard conformal factor $\phi$ becomes very sharp, whereas $\chi=e^{-4\phi}$ is far smoother (see, e.g., [Campanelli, Lousto, Marronetti, and Zlochower (2006)](https://arxiv.org/abs/gr-qc/0511048) for additional discussion). Thus if we choose to rewrite derivatives of $\phi$ in the BSSN equations in terms of finite-difference derivatives `cf`$=\chi$, numerical errors will be far smaller near the puncture.

The BSSN modules in NRPy+ support three options for the conformal factor variable `cf`:

1. `cf`$=\phi$,
1. `cf`$=\chi=e^{-4\phi}$, and
1. `cf`$=W = e^{-2\phi}$.

The BSSN equations are written in terms of $\phi$ (actually only $e^{-4\phi}$ appears) and derivatives of $\phi$, we now define $e^{-4\phi}$ and derivatives of $\phi$ in terms of the chosen `cf`.

First, we define the base variables needed within the BSSN equations:


```python
# Step 9: Standard BSSN conformal factor phi,
#         and its partial and covariant derivatives,
#         all in terms of BSSN gridfunctions like cf

# Step 9.a.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
cf_dD = ixp.declarerank1("cf_dD")
cf_dupD = ixp.declarerank1("cf_dupD") # Needed for \partial_t \phi next.
cf_dDD = ixp.declarerank2("cf_dDD","sym01")
phi_dD = ixp.zerorank1()
phi_dupD = ixp.zerorank1()
phi_dDD = ixp.zerorank2()
exp_m4phi = sp.sympify(0)
```

Then we define $\phi_{,i}$, $\phi_{,ij}$, and $e^{-4\phi}$ for each of the choices of `cf`.

For `cf`$=\phi$, this is trivial:


```python
# Step 9.a.ii: Assuming cf=phi, define exp_m4phi, phi_dD,
#              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
if par.parval_from_str(thismodule+"::EvolvedConformalFactor_cf") == "phi":
    for i in range(DIM):
        phi_dD[i] = cf_dD[i]
        phi_dupD[i] = cf_dupD[i]
        for j in range(DIM):
            phi_dDD[i][j] = cf_dDD[i][j]
    exp_m4phi = sp.exp(-4*cf)
```

For `cf`$=W=e^{-2\phi}$, we have

* $\phi_{,i} = -\text{cf}_{,i} / (2 \text{cf})$
* $\phi_{,ij} = (-\text{cf}_{,ij} + \text{cf}_{,i}\text{cf}_{,j}/\text{cf}) / (2 \text{cf})$
* $e^{-4\phi} = \text{cf}^2$

***Exercise to student: Prove the above relations***


```python
# Step 9.a.iii: Assuming cf=W=e^{-2 phi}, define exp_m4phi, phi_dD,
#               phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
if par.parval_from_str(thismodule+"::EvolvedConformalFactor_cf") == "W":
    # \partial_i W = \partial_i (e^{-2 phi}) = -2 e^{-2 phi} \partial_i phi
    # -> \partial_i phi = -\partial_i cf / (2 cf)
    for i in range(DIM):
        phi_dD[i] = - cf_dD[i] / (2*cf)
        phi_dupD[i] = - cf_dupD[i] / (2*cf)
        for j in range(DIM):
            # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (2 cf)]
            #                           = - cf_{,ij} / (2 cf) + \partial_i cf \partial_j cf / (2 cf^2)
            phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i]*cf_dD[j] / cf) / (2*cf)
    exp_m4phi = cf*cf
```

For `cf`$=\chi=e^{-4\phi}$, we have

* $\phi_{,i} = -\text{cf}_{,i} / (4 \text{cf})$
* $\phi_{,ij} = (-\text{cf}_{,ij} + \text{cf}_{,i}\text{cf}_{,j}/\text{cf}) / (4 \text{cf})$
* $e^{-4\phi} = \text{cf}$

***Exercise to student: Prove the above relations***


```python
# Step 9.a.iv: Assuming cf=chi=e^{-4 phi}, define exp_m4phi, phi_dD,
#              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
if par.parval_from_str(thismodule+"::EvolvedConformalFactor_cf") == "chi":
    # \partial_i chi = \partial_i (e^{-4 phi}) = -4 e^{-4 phi} \partial_i phi
    # -> \partial_i phi = -\partial_i cf / (4 cf)
    for i in range(DIM):
        phi_dD[i] = - cf_dD[i] / (4*cf)
        phi_dupD[i] = - cf_dupD[i] / (4*cf)
        for j in range(DIM):
            # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (4 cf)]
            #                           = - cf_{,ij} / (4 cf) + \partial_i cf \partial_j cf / (4 cf^2)
            phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i]*cf_dD[j] / cf) / (4*cf)
    exp_m4phi = cf

# Step 9.a.v: Error out if unsupported EvolvedConformalFactor_cf choice is made:
cf_choice = par.parval_from_str(thismodule+"::EvolvedConformalFactor_cf")
if cf_choice not in ('phi', 'W', 'chi'):
    print("Error: EvolvedConformalFactor_cf == "+par.parval_from_str(thismodule+"::EvolvedConformalFactor_cf")+" unsupported!")
    sys.exit(1)
```

<a id='phi_covariant_derivs'></a>

## Step 9.b: Covariant derivatives of $\phi$ \[Back to [top](#toc)\]
$$\label{phi_covariant_derivs}$$

Since $\phi$ is a scalar, $\bar{D}_i \phi = \partial_i \phi$.

Thus the second covariant derivative is given by
\begin{align}
\bar{D}_i \bar{D}_j \phi &= \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}\\
                         &= \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}.
\end{align}


```python
# Step 9.b: Define phi_dBarD = phi_dD (since phi is a scalar) and phi_dBarDD (covariant derivative)
#          \bar{D}_i \bar{D}_j \phi = \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}
#                                   = \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}
phi_dBarD = phi_dD
phi_dBarDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        phi_dBarDD[i][j] = phi_dDD[i][j]
        for k in range(DIM):
            phi_dBarDD[i][j] += - GammabarUDD[k][i][j]*phi_dD[k]
```

<a id='code_validation'></a>

# Step 10: Code validation against `BSSN.BSSN_quantities` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_quantities](../edit/BSSN/BSSN_quantities.py) module.

By default, we analyze the RHSs in Spherical coordinates, though other coordinate systems may be chosen.


```python
def comp_func(expr1,expr2,basename,prefixname2="Bq."):
    if str(expr1-expr2)!="0":
        print(basename+" - "+prefixname2+basename+" = "+ str(expr1-expr2))
        return 1
    return 0

def gfnm(basename,idx1,idx2=None,idx3=None):
    if idx2 is None:
        return basename+"["+str(idx1)+"]"
    if idx3 is None:
        return basename+"["+str(idx1)+"]["+str(idx2)+"]"
    return basename+"["+str(idx1)+"]["+str(idx2)+"]["+str(idx3)+"]"

expr_list = []
exprcheck_list = []
namecheck_list = []

# Step 3:
import BSSN.BSSN_quantities as Bq
Bq.BSSN_basic_tensors()
for i in range(DIM):
    namecheck_list.extend([gfnm("LambdabarU",i),gfnm("betaU",i),gfnm("BU",i)])
    exprcheck_list.extend([Bq.LambdabarU[i],Bq.betaU[i],Bq.BU[i]])
    expr_list.extend([LambdabarU[i],betaU[i],BU[i]])
    for j in range(DIM):
        namecheck_list.extend([gfnm("gammabarDD",i,j),gfnm("AbarDD",i,j)])
        exprcheck_list.extend([Bq.gammabarDD[i][j],Bq.AbarDD[i][j]])
        expr_list.extend([gammabarDD[i][j],AbarDD[i][j]])

# Step 4:
Bq.gammabar__inverse_and_derivs()
for i in range(DIM):
    for j in range(DIM):
        namecheck_list.extend([gfnm("gammabarUU",i,j)])
        exprcheck_list.extend([Bq.gammabarUU[i][j]])
        expr_list.extend([gammabarUU[i][j]])
        for k in range(DIM):
            namecheck_list.extend([gfnm("gammabarDD_dD",i,j,k),
                                   gfnm("gammabarDD_dupD",i,j,k),
                                   gfnm("GammabarUDD",i,j,k)])
            exprcheck_list.extend([Bq.gammabarDD_dD[i][j][k],Bq.gammabarDD_dupD[i][j][k],Bq.GammabarUDD[i][j][k]])
            expr_list.extend(        [gammabarDD_dD[i][j][k],gammabarDD_dupD[i][j][k],GammabarUDD[i][j][k]])

# Step 5:
Bq.detgammabar_and_derivs()
namecheck_list.extend(["detgammabar"])
exprcheck_list.extend([Bq.detgammabar])
expr_list.extend([detgammabar])
for i in range(DIM):
    namecheck_list.extend([gfnm("detgammabar_dD",i)])
    exprcheck_list.extend([Bq.detgammabar_dD[i]])
    expr_list.extend([detgammabar_dD[i]])
    for j in range(DIM):
        namecheck_list.extend([gfnm("detgammabar_dDD",i,j)])
        exprcheck_list.extend([Bq.detgammabar_dDD[i][j]])
        expr_list.extend([detgammabar_dDD[i][j]])

# Step 6:
Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()
namecheck_list.extend(["trAbar"])
exprcheck_list.extend([Bq.trAbar])
expr_list.extend([trAbar])
for i in range(DIM):
    for j in range(DIM):
        namecheck_list.extend([gfnm("AbarUU",i,j),gfnm("AbarUD",i,j)])
        exprcheck_list.extend([Bq.AbarUU[i][j],Bq.AbarUD[i][j]])
        expr_list.extend([AbarUU[i][j],AbarUD[i][j]])
        for k in range(DIM):
            namecheck_list.extend([gfnm("AbarDD_dD",i,j,k)])
            exprcheck_list.extend([Bq.AbarDD_dD[i][j][k]])
            expr_list.extend([AbarDD_dD[i][j][k]])

# Step 7:
Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
for i in range(DIM):
    namecheck_list.extend([gfnm("DGammaU",i)])
    exprcheck_list.extend([Bq.DGammaU[i]])
    expr_list.extend([DGammaU[i]])

    for j in range(DIM):
        namecheck_list.extend([gfnm("RbarDD",i,j)])
        exprcheck_list.extend([Bq.RbarDD[i][j]])
        expr_list.extend([RbarDD[i][j]])
        for k in range(DIM):
            namecheck_list.extend([gfnm("DGammaUDD",i,j,k),gfnm("gammabarDD_dHatD",i,j,k)])
            exprcheck_list.extend([Bq.DGammaUDD[i][j][k],Bq.gammabarDD_dHatD[i][j][k]])
            expr_list.extend([DGammaUDD[i][j][k],gammabarDD_dHatD[i][j][k]])

# Step 8:
Bq.betaU_derivs()
for i in range(DIM):
    for j in range(DIM):
        namecheck_list.extend([gfnm("betaU_dD",i,j),gfnm("betaU_dupD",i,j)])
        exprcheck_list.extend([Bq.betaU_dD[i][j],Bq.betaU_dupD[i][j]])
        expr_list.extend([betaU_dD[i][j],betaU_dupD[i][j]])
        for k in range(DIM):
            namecheck_list.extend([gfnm("betaU_dDD",i,j,k)])
            exprcheck_list.extend([Bq.betaU_dDD[i][j][k]])
            expr_list.extend([betaU_dDD[i][j][k]])

# Step 9:
Bq.phi_and_derivs()
#phi_dD,phi_dupD,phi_dDD,exp_m4phi,phi_dBarD,phi_dBarDD
namecheck_list.extend(["exp_m4phi"])
exprcheck_list.extend([Bq.exp_m4phi])
expr_list.extend([exp_m4phi])
for i in range(DIM):
    namecheck_list.extend([gfnm("phi_dD",i),gfnm("phi_dupD",i),gfnm("phi_dBarD",i)])
    exprcheck_list.extend([Bq.phi_dD[i],Bq.phi_dupD[i],Bq.phi_dBarD[i]])
    expr_list.extend(        [phi_dD[i],phi_dupD[i],phi_dBarD[i]])
    for j in range(DIM):
        namecheck_list.extend([gfnm("phi_dDD",i,j),gfnm("phi_dBarDD",i,j)])
        exprcheck_list.extend([Bq.phi_dDD[i][j],Bq.phi_dBarDD[i][j]])
        expr_list.extend([phi_dDD[i][j],phi_dBarDD[i][j]])

num_failures = 0
for i in range(len(expr_list)):
    num_failures += comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])

if num_failures == 0:
    print("ALL TESTS PASSED!")
else:
    print(str(num_failures) + " TESTS FAILED")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 11: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_quantities.pdf](Tutorial-BSSN_quantities.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_quantities")
```

    Created Tutorial-BSSN_quantities.tex, and compiled LaTeX file to PDF file
        Tutorial-BSSN_quantities.pdf

