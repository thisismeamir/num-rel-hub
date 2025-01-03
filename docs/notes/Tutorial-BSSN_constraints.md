<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [BSSN](http://www2.yukawa.kyoto-u.ac.jp/~yuichiro.sekiguchi/3+1.pdf) Hamiltonian and momentum constraint equations, in ***curvilinear*** coordinates, using a covariant reference metric approach: C code generation

## Authors: Ian Ruchlin & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook demonstrates the construction of the BSSN Hamiltonian and momentum constraint equations as symbolic (SymPy) expressions, in terms of the core BSSN quantities $\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$, as defined in [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658). The module implements a generic curvilinear coordinate reference metric approach which is an extension of the spherical coordinate reference metric approach of [Baumgarte, Montero, Cordero-Carrión, and Müller (2012)](https://arxiv.org/abs/1211.6632), building upon the covariant Lagrangian BSSN formalism of [Brown (2009)](https://arxiv.org/abs/0902.3652). 

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** All expressions generated in this module have been validated against a trusted code where applicable (the original NRPy+/SENR code, which itself was validated against [Baumgarte's code](https://arxiv.org/abs/1211.6632)).

### NRPy+ Source Code for this module: [BSSN/BSSN_constraints.py](../edit/BSSN/BSSN_constraints.py)


[comment]: <> (Introduction: TODO)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
1. [Step 2](#hamiltonianconstraint): Construct the Hamiltonian constraint $\mathcal{H}$.
1. [Step 3](#momentumconstraint): Construct the momentum constraint $\mathcal{M}^i$.
1. [Step 4](#code_validation): Code Validation against `BSSN.BSSN_constraints` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

We start by loading the needed modules. Notably, this module depends on several quantities defined in the [BSSN/BSSN_quantities.py](../edit/BSSN/BSSN_quantities.py) Python code, documented in the NRPy+ [BSSN quantities](Tutorial-BSSN_quantities.ipynb). In [Step 2](#hamiltonianconstraint) we call functions within [BSSN.BSSN_quantities](../edit/BSSN/BSSN_quantities.py) to define quantities needed in this module.


```python
# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp               # SymPy, Python's core symbolic algebra package on which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import reference_metric as rfm   # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq

# Step 1.a: Set spatial dimension (must be 3 for BSSN, as BSSN is
#           a 3+1-dimensional decomposition of the general
#           relativistic field equations)
DIM = 3

# Step 1.b: Given the chosen coordinate system, set up
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
rfm.reference_metric()
```

<a id='hamiltonianconstraint'></a>

# Step 2: $\mathcal{H}$, the Hamiltonian constraint \[Back to [top](#toc)\]
$$\label{hamiltonianconstraint}$$

Next we define the Hamiltonian constraint. Eq. 13 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf) yields:
$$
\mathcal{H} = {\underbrace {\textstyle \frac{2}{3} K^2}_{\rm Term\ 1}} - 
{\underbrace {\textstyle \bar{A}_{ij} \bar{A}^{ij}}_{\rm Term\ 2}} + 
{\underbrace {\textstyle e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right)}_{\rm Term\ 3}}
$$


```python
# Step 2: The Hamiltonian constraint.
# First declare all needed variables
Bq.declare_BSSN_gridfunctions_if_not_declared_already() # Sets trK
Bq.BSSN_basic_tensors()   # Sets AbarDD
Bq.gammabar__inverse_and_derivs() # Sets gammabarUU
Bq.AbarUU_AbarUD_trAbar_AbarDD_dD() # Sets AbarUU and AbarDD_dD
Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU() # Sets RbarDD
Bq.phi_and_derivs() # Sets phi_dBarD & phi_dBarDD

# Term 1: 2/3 K^2
H = sp.Rational(2,3)*Bq.trK**2

# Term 2: -A_{ij} A^{ij}
for i in range(DIM):
    for j in range(DIM):
        H += -Bq.AbarDD[i][j]*Bq.AbarUU[i][j]

# Term 3a: trace(Rbar)
Rbartrace = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        Rbartrace += Bq.gammabarUU[i][j]*Bq.RbarDD[i][j]

# Term 3b: -8 \bar{\gamma}^{ij} \bar{D}_i \phi \bar{D}_j \phi = -8*phi_dBar_times_phi_dBar
# Term 3c: -8 \bar{\gamma}^{ij} \bar{D}_i \bar{D}_j \phi      = -8*phi_dBarDD_contraction
phi_dBar_times_phi_dBar = sp.sympify(0) # Term 3b
phi_dBarDD_contraction  = sp.sympify(0) # Term 3c
for i in range(DIM):
    for j in range(DIM):
        phi_dBar_times_phi_dBar += Bq.gammabarUU[i][j]*Bq.phi_dBarD[i]*Bq.phi_dBarD[j]
        phi_dBarDD_contraction  += Bq.gammabarUU[i][j]*Bq.phi_dBarDD[i][j]

# Add Term 3:
H += Bq.exp_m4phi*(Rbartrace - 8*(phi_dBar_times_phi_dBar + phi_dBarDD_contraction))
```

<a id='momentumconstraint'></a>

# Step 3: $\mathcal{M}^i$, the momentum constraint \[Back to [top](#toc)\]
$$\label{momentumconstraint}$$

***Courtesy Ian Ruchlin***

The following definition of the momentum constraint is a simplification of Eq. 47 or [Ruchlin, Etienne, & Baumgarte (2018)](https://arxiv.org/pdf/1712.07658.pdf), which itself was a corrected version of the momentum constraint presented in Eq. 14 of [Baumgarte *et al*](https://arxiv.org/pdf/1211.6632.pdf).

Start with the physical momentum constraint
$$
\mathcal{M}^{i} \equiv D_{j} \left ( K^{i j} - \gamma^{i j} K \right ) = 0 \; .
$$
Expanding and using metric compatibility with the physical covariant derivative $D_{i}$ yields
$$
\mathcal{M}^{i} = D_{j} K^{i j} - \gamma^{i j} \partial_{j} K \; .
$$
The physical extrinsic curvature $K_{i j}$ is related to the trace-free extrinsic curvature $A_{i j}$ by
$$
K_{i j} = A_{i j} + \frac{1}{3} \gamma_{i j} K \; .
$$
Thus,
$$
\mathcal{M}^{i} = D_{j} A^{i j} - \frac{2}{3} \gamma^{i j} \partial_{j} K \; .
$$
The physical metric $\gamma_{i j}$ is related to the conformal metric $\bar{\gamma}_{i j}$ by the conformal rescaling
$$
\gamma_{i j} = e^{4 \phi} \bar{\gamma}_{i j} \; ,
$$
and similarly for the trace-free extrinsic curvature
$$
A_{i j} = e^{4 \phi} \bar{A}_{i j} \; .
$$
It can be shown (Eq. (3.34) in Baumgarte & Shapiro (2010) with $\alpha = -4$ and $\psi = e^{\phi}$) that the physical and conformal covariant derivatives obey
$$
D_{j} A^{i j} = e^{-10 \phi} \bar{D}_{j} \left (e^{6 \phi} \bar{A}^{i j} \right ) \; .
$$
Then, the constraint becomes
$$
\mathcal{M}^i = e^{-4\phi} \left(
{\underbrace {\textstyle \bar{D}_j \bar{A}^{ij}}_{\rm Term\ 1}} + 
{\underbrace {\textstyle 6 \bar{A}^{ij}\partial_j \phi}_{\rm Term\ 2}} - 
{\underbrace {\textstyle \frac{2}{3} \bar{\gamma}^{ij}\partial_j K}_{\rm Term\ 3}}\right) \; .
$$

Let's first implement Terms 2 and 3:


```python
# Step 3: M^i, the momentum constraint

MU = ixp.zerorank1()

# Term 2: 6 A^{ij} \partial_j \phi:
for i in range(DIM):
    for j in range(DIM):
        MU[i] += 6*Bq.AbarUU[i][j]*Bq.phi_dD[j]

# Term 3: -2/3 \bar{\gamma}^{ij} K_{,j}
trK_dD = ixp.declarerank1("trK_dD") # Not defined in BSSN_RHSs; only trK_dupD is defined there.
for i in range(DIM):
    for j in range(DIM):
        MU[i] += -sp.Rational(2,3)*Bq.gammabarUU[i][j]*trK_dD[j]
```

Now, we turn our attention to Term 1. The covariant divergence involves upper indices in $\bar{A}^{i j}$, but it would be easier for us to finite difference the rescaled $\bar{A}_{i j}$. A simple application of the inverse conformal metric yields
$$
\bar{D}_{j} \bar{A}^{i j} = \bar{\gamma}^{i k} \bar{\gamma}^{j l} \bar{D}_{j} \bar{A}_{k l} \; .
$$
As usual, the covariant derivative is related to the ordinary derivative using the conformal Christoffel symbols
$$
\bar{D}_{k} \bar{A}_{i j} = \partial_{k} \bar{A}_{i j} - \bar{\Gamma}^{l}_{k i} \bar{A}_{l j} - \bar{\Gamma}^{l}_{k j} \bar{A}_{i l} \; .
$$
It is the ordinary derivative above that is approximated by finite difference. The BSSN formulation used here does not rely on spatial derivatives $\partial_{k} \bar{A}_{i j}$ in any of the right-hand-sides (except for the advection term, which uses the upwinded derivative), and so we must declare a new ordinary, centered stencil derivative field of rank 3.


```python
# First define aDD_dD:
aDD_dD = ixp.declarerank3("aDD_dD","sym01")

# Then evaluate the conformal covariant derivative \bar{D}_j \bar{A}_{lm}
AbarDD_dBarD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            AbarDD_dBarD[i][j][k] = Bq.AbarDD_dD[i][j][k]
            for l in range(DIM):
                AbarDD_dBarD[i][j][k] += -Bq.GammabarUDD[l][k][i]*Bq.AbarDD[l][j]
                AbarDD_dBarD[i][j][k] += -Bq.GammabarUDD[l][k][j]*Bq.AbarDD[i][l]

# Term 1: Contract twice with the metric to make \bar{D}_{j} \bar{A}^{ij}
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                MU[i] += Bq.gammabarUU[i][k]*Bq.gammabarUU[j][l]*AbarDD_dBarD[k][l][j]

# Finally, we multiply by e^{-4 phi} and rescale the momentum constraint:
for i in range(DIM):
    MU[i] *= Bq.exp_m4phi / rfm.ReU[i]
```

<a id='code_validation'></a>

# Step 4: Code Validation against `BSSN.BSSN_constraints` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_constraints](../edit/BSSN/BSSN_constraints.py) module.

By default, we analyze these expressions in Spherical coordinates, though other coordinate systems may be chosen.


```python
# Step 4: Code Validation against BSSN.BSSN_constraints NRPy+ module

# We already have SymPy expressions for BSSN constraints
#         in terms of other SymPy variables. Even if we reset the
#         list of NRPy+ gridfunctions, these *SymPy* expressions for
#         BSSN constraint variables *will remain unaffected*.
#
#         Here, we will use the above-defined BSSN constraint expressions
#         to validate against the same expressions in the
#         BSSN/BSSN_constraints.py file, to ensure consistency between
#         this tutorial and the module itself.
#
# Reset the list of gridfunctions, as registering a gridfunction
#   twice (in the bssnrhs.BSSN_RHSs() call) will spawn an error.
gri.glb_gridfcs_list = []

# Call the BSSN_RHSs() function from within the
#          BSSN/BSSN_RHSs.py module,
#          which should do exactly the same as in Steps 1-16 above.
import BSSN.BSSN_constraints as bssncon
bssncon.BSSN_constraints()

print("Consistency check between BSSN_constraints tutorial and NRPy+ module: ALL SHOULD BE ZERO.")

print("H - bssncon.H = " + str(H - bssncon.H))
for i in range(DIM):
    print("MU["+str(i)+"] - bssncon.MU["+str(i)+"] = " + str(MU[i] - bssncon.MU[i]))
```

    Consistency check between BSSN_constraints tutorial and NRPy+ module: ALL SHOULD BE ZERO.
    H - bssncon.H = 0
    MU[0] - bssncon.MU[0] = 0
    MU[1] - bssncon.MU[1] = 0
    MU[2] - bssncon.MU[2] = 0


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_constraints.pdf](Tutorial-BSSN_constraints.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_constraints")
```

    Created Tutorial-BSSN_constraints.tex, and compiled LaTeX file to PDF file
        Tutorial-BSSN_constraints.pdf

