<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# BSSN Time-Evolution Equations for the Gauge Fields $\alpha$ and $\beta^i$

## Authors: Zach Etienne & Terrence Pierre Jacques
### Formatting improvements courtesy Brandon Clark

[comment]: <> (Abstract: TODO, or make the introduction an abstract and additional notes section, and write a new Introduction)

## This notebook constructs SymPy expressions for the right-hand sides of the time-evolution equations for the gauge fields $\alpha$ (the lapse) and $\beta^i$ (the shift) in the BSSN formalism. The tutorial explores various gauge conditions and their robustness in the presence of black holes. Various lapse and shift conditions, including the $1+\log$ lapse, harmonic slicing, frozen lapse, and gamma-driving shift conditions, are elaborated.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** All expressions generated in this module have been validated against a trusted code (the original NRPy+/SENR code, which itself was validated against [Baumgarte's code](https://arxiv.org/abs/1211.6632)).

### NRPy+ Source Code for this module: [BSSN/BSSN_gauge_RHSs.py](../edit/BSSN/BSSN_gauge_RHSs.py)


## Introduction:
This tutorial demonstrates the generation of SymPy expressions for the time-evolution equations concerning the gauge fields $\alpha$ and $\beta^i$, fundamental components in the 3+1 solution to Einstein’s equations. $\alpha$, or the lapse, dictates the amount of proper time elapsing between one timestep to the next at each point, while $\beta^i$, known as the shift, determines the proper distance numerical grid points travel from one timestep to the next. Although there are no strict requirements for the selection of gauge conditions, few have shown robustness in scenarios with black holes, hence the focus of this notebook is limited to the most stable selections.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
1. [Step 2](#lapseconditions): Right-hand side of $\partial_t \alpha$
    1. [Step 2.a](#onepluslog): $1+\log$ lapse
    1. [Step 2.b](#harmonicslicing): Harmonic slicing
    1. [Step 2.c](#frozen): Frozen lapse
    1. [Step 2.d](#statictrumpet_onepluslog): Alternative 1+log condition for Static Trumpet initial data
1. [Step 3](#shiftconditions): Right-hand side of $\partial_t \beta^i$: Second-order Gamma-driving shift conditions
    1. [Step 3.a](#origgammadriving): Original, non-covariant Gamma-driving shift condition
    1. [Step 3.b](#covgammadriving): [Brown](https://arxiv.org/abs/0902.3652)'s suggested covariant Gamma-driving shift condition
        1. [Step 3.b.i](#partial_beta): The right-hand side of the $\partial_t \beta^i$ equation
        1. [Step 3.b.ii](#partial_upper_b): The right-hand side of the $\partial_t B^i$ equation
    1. [Step 3.c](#statictrumpet_nonadvecgammadriving): Non-advecting Gamma-driving shift condition (used for evolving "Static Trumpet" initial data)
1. [Step 4](#rescale): Rescale right-hand sides of BSSN gauge equations
1. [Step 5](#code_validation): Code Validation against `BSSN.BSSN_gauge_RHSs` NRPy+ module
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from Python/NRPy+:


```python
# Step 1: Import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities
import BSSN.BSSN_RHSs as Brhs     # NRPy+: Constructs BSSN right-hand-side expressions
import sys                        # Standard Python modules for multiplatform OS-level functions

# Step 1.c: Declare/initialize parameters for this module
thismodule = "BSSN_gauge_RHSs"
par.initialize_param(par.glb_param("char", thismodule, "LapseEvolutionOption", "OnePlusLog"))
par.initialize_param(par.glb_param("char", thismodule, "ShiftEvolutionOption", "GammaDriving2ndOrder_Covariant"))

# Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
#           a 3+1-dimensional decomposition of the general
#           relativistic field equations)
DIM = 3

# Step 1.e: Given the chosen coordinate system, set up
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
rfm.reference_metric()

# Step 1.f: Define BSSN scalars & tensors (in terms of rescaled BSSN quantities)
import BSSN.BSSN_quantities as Bq
Bq.BSSN_basic_tensors()
Bq.betaU_derivs()

import BSSN.BSSN_RHSs as Brhs
Brhs.BSSN_RHSs()
```

<a id='lapseconditions'></a>

# Step 2: Right-hand side of $\partial_t \alpha$ \[Back to [top](#toc)\]
$$\label{lapseconditions}$$

<a id='onepluslog'></a>

## Step 2.a: $1+\log$ lapse \[Back to [top](#toc)\]
$$\label{onepluslog}$$

The [$1=\log$ lapse condition](https://arxiv.org/abs/gr-qc/0206072) is a member of the [Bona-Masso family of lapse choices](https://arxiv.org/abs/gr-qc/9412071), which has the desirable property of singularity avoidance. As is common (e.g., see [Campanelli *et al* (2005)](https://arxiv.org/abs/gr-qc/0511048)), we make the replacement $\partial_t \to \partial_0 = \partial_t - \beta^i \partial_i$ to ensure lapse characteristics advect with the shift. The bracketed term in the $1+\log$ lapse condition below encodes the shift advection term:

\begin{align}
\partial_0 \alpha &= -2 \alpha K \\
\implies \partial_t \alpha &= \left[\beta^i \partial_i \alpha\right] - 2 \alpha K
\end{align}


```python
# Step 2.a: The 1+log lapse condition:
#   \partial_t \alpha = \beta^i \alpha_{,i} - 2*\alpha*K
# First import expressions from BSSN_quantities
cf    = Bq.cf
trK   = Bq.trK
alpha = Bq.alpha
betaU = Bq.betaU

# Implement the 1+log lapse condition
if par.parval_from_str(thismodule+"::LapseEvolutionOption") == "OnePlusLog":
    alpha_rhs = -2*alpha*trK
    alpha_dupD = ixp.declarerank1("alpha_dupD")
    for i in range(DIM):
        alpha_rhs += betaU[i]*alpha_dupD[i]
```

<a id='harmonicslicing'></a>

## Step 2.b: Harmonic slicing \[Back to [top](#toc)\]
$$\label{harmonicslicing}$$

As defined on Pg 2 of https://arxiv.org/pdf/gr-qc/9902024.pdf , this is given by 

$$
\partial_t \alpha = \partial_t e^{6 \phi} = 6 e^{6 \phi} \partial_t \phi
$$

If 

$$\text{cf} = W = e^{-2 \phi},$$ 

then

$$
6 e^{6 \phi} \partial_t \phi = 6 W^{-3} \partial_t \phi.
$$

However,
$$
\partial_t \phi = -\partial_t \text{cf} / (2 \text{cf})$$

(as described above), so if `cf`$=W$, then
\begin{align}
\partial_t \alpha &= 6 e^{6 \phi} \partial_t \phi \\
&= 6 W^{-3} \left(-\frac{\partial_t W}{2 W}\right) \\
&= -3 \text{cf}^{-4} \text{cf}\_\text{rhs}
\end{align}

**Exercise to students: Implement Harmonic slicing for `cf`$=\chi$** 


```python
# Step 2.b: Implement the harmonic slicing lapse condition
if par.parval_from_str(thismodule+"::LapseEvolutionOption") == "HarmonicSlicing":
    if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "W":
        alpha_rhs = -3*cf**(-4)*Brhs.cf_rhs
    elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "phi":
        alpha_rhs = 6*sp.exp(6*cf)*Brhs.cf_rhs
    else:
        print("Error LapseEvolutionOption==HarmonicSlicing unsupported for EvolvedConformalFactor_cf!=(W or phi)")
        sys.exit(1)
```

<a id='frozen'></a>

## Step 2.c: Frozen lapse \[Back to [top](#toc)\]
$$\label{frozen}$$

This slicing condition is given by
$$\partial_t \alpha = 0,$$

which is rarely a stable lapse condition.


```python
# Step 2.c: Frozen lapse
#    \partial_t \alpha = 0
if par.parval_from_str(thismodule+"::LapseEvolutionOption") == "Frozen":
    alpha_rhs = sp.sympify(0)
```

<a id='statictrumpet_onepluslog'></a>

## Step 2.d: Alternative $1+\log$ lapse for Static Trumpet initial data \[Back to [top](#toc)\]
$$\label{statictrumpet_onepluslog}$$

An alternative to the standard 1+log condition to be used with Static Trumpet initial data, the lapse is evolved according to a
condition consistent with staticity, given by equation 67 in [Ruchlin, Etienne, & Baumgarte (2018)](https://arxiv.org/pdf/1712.07658.pdf)

\begin{align}
\partial_0 \alpha &= -\alpha(1 - \alpha) K \\
\implies \partial_t \alpha &= \left[\beta^i \partial_i \alpha\right] -\alpha(1 - \alpha) K
\end{align}


```python
# Step 2.d: Alternative 1+log lapse condition:
#   \partial_t \alpha = \beta^i \alpha_{,i} -\alpha*(1 - \alpha)*K

# Implement the alternative 1+log lapse condition
if par.parval_from_str(thismodule+"::LapseEvolutionOption") == "OnePlusLogAlt":
    alpha_rhs = -alpha*(1 - alpha)*trK
    alpha_dupD = ixp.declarerank1("alpha_dupD")
    for i in range(DIM):
        alpha_rhs += betaU[i]*alpha_dupD[i]
```

<a id='shiftconditions'></a>

# Step 3: Right-hand side of $\partial_t \beta^i$: Second-order Gamma-driving shift conditions \[Back to [top](#toc)\]
$$\label{shiftconditions}$$

The motivation behind Gamma-driving shift conditions is well documented in the book [*Numerical Relativity* by Baumgarte & Shapiro](https://www.amazon.com/Numerical-Relativity-Einsteins-Equations-Computer/dp/052151407X/).

<a id='origgammadriving'></a>

## Step 3.a: Original, non-covariant Gamma-driving shift condition \[Back to [top](#toc)\]
$$\label{origgammadriving}$$

**Option 1: Non-Covariant, Second-Order Shift**

We adopt the [*shifting (i.e., advecting) shift*](https://arxiv.org/abs/gr-qc/0605030) non-covariant, second-order shift condition:
\begin{align}
\partial_0 \beta^i &= B^{i} \\
\partial_0 B^i &= \frac{3}{4} \partial_{0} \bar{\Lambda}^{i} - \eta B^{i} \\
\implies \partial_t \beta^i &= \left[\beta^j \partial_j \beta^i\right] + B^{i} \\
\partial_t B^i &= \left[\beta^j \partial_j B^i\right] + \frac{3}{4} \partial_{0} \bar{\Lambda}^{i} - \eta B^{i},
\end{align}
where $\eta$ is the shift damping parameter, and $\partial_{0} \bar{\Lambda}^{i}$ in the right-hand side of the $\partial_{0} B^{i}$ equation is computed by adding $\beta^j \partial_j \bar{\Lambda}^i$ to the right-hand side expression given for $\partial_t \bar{\Lambda}^i$ in the BSSN time-evolution equations as listed [here](Tutorial-BSSN_formulation.ipynb), so no explicit time dependence occurs in the right-hand sides of the BSSN evolution equations and the Method of Lines can be applied directly.


```python
# Step 3.a: Set \partial_t \beta^i
# First import expressions from BSSN_quantities
BU         = Bq.BU
betU       = Bq.betU
betaU_dupD = Bq.betaU_dupD
# Define needed quantities
beta_rhsU = ixp.zerorank1()
B_rhsU = ixp.zerorank1()
if par.parval_from_str(thismodule+"::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
    # Step 3.a.i: Compute right-hand side of beta^i
    # *  \partial_t \beta^i = \beta^j \beta^i_{,j} + B^i
    for i in range(DIM):
        beta_rhsU[i] += BU[i]
        for j in range(DIM):
            beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]
    # Compute right-hand side of B^i:
    eta = par.Cparameters("REAL", thismodule, ["eta"],2.0)

    # Step 3.a.ii: Compute right-hand side of B^i
    # *  \partial_t B^i     = \beta^j B^i_{,j} + 3/4 * \partial_0 \Lambda^i - eta B^i
    # Step 3.a.iii: Define BU_dupD, in terms of derivative of rescaled variable \bet^i
    BU_dupD = ixp.zerorank2()
    betU_dupD = ixp.declarerank2("betU_dupD","nosym")
    for i in range(DIM):
        for j in range(DIM):
            BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]

    # Step 3.a.iv: Compute \partial_0 \bar{\Lambda}^i = (\partial_t - \beta^i \partial_i) \bar{\Lambda}^j
    Lambdabar_partial0 = ixp.zerorank1()
    for i in range(DIM):
        Lambdabar_partial0[i] = Brhs.Lambdabar_rhsU[i]
    for i in range(DIM):
        for j in range(DIM):
            Lambdabar_partial0[j] += -betaU[i]*Brhs.LambdabarU_dupD[j][i]

    # Step 3.a.v: Evaluate RHS of B^i:
    for i in range(DIM):
        B_rhsU[i] += sp.Rational(3,4)*Lambdabar_partial0[i] - eta*BU[i]
        for j in range(DIM):
            B_rhsU[i] += betaU[j]*BU_dupD[i][j]
```

<a id='covgammadriving'></a>

## Step 3.b: [Brown](https://arxiv.org/abs/0902.3652)'s suggested covariant Gamma-driving shift condition \[Back to [top](#toc)\]
$$\label{covgammadriving}$$

<a id='partial_beta'></a>

### Step 3.b.i: The right-hand side of the $\partial_t \beta^i$ equation \[Back to [top](#toc)\]
$$\label{partial_beta}$$

This is [Brown's](https://arxiv.org/abs/0902.3652) suggested formulation (Eq. 20b; note that Eq. 20a is the same as our lapse condition, as $\bar{D}_j \alpha = \partial_j \alpha$ for scalar $\alpha$):
$$\partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}$$
Based on the definition of the covariant derivative, we have
$$
\bar{D}_{j} \beta^{i} = \beta^i_{,j} + \bar{\Gamma}^i_{mj} \beta^m,
$$
so the above equation becomes
\begin{align}
\partial_t \beta^i &= \left[\beta^j \left(\beta^i_{,j} + \bar{\Gamma}^i_{mj} \beta^m\right)\right] + B^{i}\\
&= {\underbrace {\textstyle \beta^j \beta^i_{,j}}_{\text{Term 1}}} + 
{\underbrace {\textstyle \beta^j \bar{\Gamma}^i_{mj} \beta^m}_{\text{Term 2}}} + 
{\underbrace {\textstyle B^i}_{\text{Term 3}}} 
\end{align}


```python
# Step 3.b: The right-hand side of the \partial_t \beta^i equation
if par.parval_from_str(thismodule+"::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
    # Step 3.b Option 2: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}
    # First we need GammabarUDD, defined in Bq.gammabar__inverse_and_derivs()
    Bq.gammabar__inverse_and_derivs()
    GammabarUDD = Bq.GammabarUDD
    # Then compute right-hand side:
    # Term 1: \beta^j \beta^i_{,j}
    for i in range(DIM):
        for j in range(DIM):
            beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]

    # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
    for i in range(DIM):
        for j in range(DIM):
            for m in range(DIM):
                beta_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*betaU[m]
    # Term 3: B^i
    for i in range(DIM):
        beta_rhsU[i] += BU[i]
```

<a id='partial_upper_b'></a>

### Step 3.b.ii: The right-hand side of the $\partial_t B^i$ equation \[Back to [top](#toc)\]
$$\label{partial_upper_b}$$

$$\partial_t B^i = \left[\beta^j \bar{D}_j B^i\right] + \frac{3}{4}\left( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} \right) - \eta B^{i}$$

Based on the definition of the covariant derivative, we have for vector $V^i$
$$
\bar{D}_{j} V^{i} = V^i_{,j} + \bar{\Gamma}^i_{mj} V^m,
$$
so the above equation becomes
\begin{align}
\partial_t B^i &= \left[\beta^j \left(B^i_{,j} + \bar{\Gamma}^i_{mj} B^m\right)\right] + \frac{3}{4}\left[ \partial_t \bar{\Lambda}^{i} - \beta^j \left(\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m\right) \right] - \eta B^{i} \\
&= {\underbrace {\textstyle \beta^j B^i_{,j}}_{\text{Term 1}}} + 
{\underbrace {\textstyle \beta^j \bar{\Gamma}^i_{mj} B^m}_{\text{Term 2}}} + 
{\underbrace {\textstyle \frac{3}{4}\partial_t \bar{\Lambda}^{i}}_{\text{Term 3}}} -
{\underbrace {\textstyle \frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}}_{\text{Term 4}}} -
{\underbrace {\textstyle \frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m}_{\text{Term 5}}} -
{\underbrace {\textstyle \eta B^i}_{\text{Term 6}}}
\end{align}


```python
if par.parval_from_str(thismodule+"::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
    # Step 3.c: Covariant option:
    #  \partial_t B^i = \beta^j \bar{D}_j B^i
    #               + \frac{3}{4} ( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} )
    #               - \eta B^{i}
    #                 = \beta^j B^i_{,j} + \beta^j \bar{\Gamma}^i_{mj} B^m
    #               + \frac{3}{4}[ \partial_t \bar{\Lambda}^{i}
    #                            - \beta^j (\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m)]
    #               - \eta B^{i}
    # Term 1, part a: First compute B^i_{,j} using upwinded derivative
    BU_dupD = ixp.zerorank2()
    betU_dupD = ixp.declarerank2("betU_dupD","nosym")
    for i in range(DIM):
        for j in range(DIM):
            BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]
    # Term 1: \beta^j B^i_{,j}
    for i in range(DIM):
        for j in range(DIM):
            B_rhsU[i] += betaU[j]*BU_dupD[i][j]
    # Term 2: \beta^j \bar{\Gamma}^i_{mj} B^m
    for i in range(DIM):
        for j in range(DIM):
            for m in range(DIM):
                B_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*BU[m]
    # Term 3: \frac{3}{4}\partial_t \bar{\Lambda}^{i}
    for i in range(DIM):
        B_rhsU[i] += sp.Rational(3,4)*Brhs.Lambdabar_rhsU[i]
    # Term 4: -\frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}
    for i in range(DIM):
        for j in range(DIM):
            B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*Brhs.LambdabarU_dupD[i][j]
    # Term 5: -\frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m
    for i in range(DIM):
        for j in range(DIM):
            for m in range(DIM):
                B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*GammabarUDD[i][m][j]*Bq.LambdabarU[m]
    # Term 6: - \eta B^i
    # eta is a free parameter; we declare it here:
    eta = par.Cparameters("REAL", thismodule, ["eta"],2.0)
    for i in range(DIM):
        B_rhsU[i] += -eta*BU[i]
```

<a id='statictrumpet_nonadvecgammadriving'></a>

## Step 3.c: Non-advecting Gamma-driving shift condition (used for evolving "Static Trumpet" initial data) \[Back to [top](#toc)\]
$$\label{statictrumpet_nonadvecgammadriving}$$


For the shift vector evolution equation, we desire only that the right-hand sides vanish analytically (although numerical error is expected to result in specious evolution). To this end, we adopt the nonadvecting Gamma-driver condition, given by equations 68a and 68b in [Ruchlin, Etienne, & Baumgarte (2018)](https://arxiv.org/pdf/1712.07658.pdf)

\begin{align}
\partial_t \beta^i &= B^{i} \\
\partial_t B^i &= \frac{3}{4} \partial_{t} \bar{\Lambda}^{i} - \eta B^{i},
\end{align}


```python
# Step 3.c: Set \partial_t \beta^i

if par.parval_from_str(thismodule+"::ShiftEvolutionOption") == "NonAdvectingGammaDriving":
    # Step 3.c.i: Compute right-hand side of beta^i
    # *  \partial_t \beta^i = B^i
    for i in range(DIM):
        beta_rhsU[i] += BU[i]

    # Compute right-hand side of B^i:
    eta = par.Cparameters("REAL", thismodule, ["eta"],2.0)

    # Step 3.c.ii: Compute right-hand side of B^i
    # *  \partial_t B^i     = 3/4 * \partial_t \Lambda^i - eta B^i
    # Step 3.c.iii: Evaluate RHS of B^i:
    for i in range(DIM):
        B_rhsU[i] += sp.Rational(3,4)*Brhs.Lambdabar_rhsU[i] - eta*BU[i]

```

<a id='rescale'></a>

# Step 4: Rescale right-hand sides of BSSN gauge equations \[Back to [top](#toc)\]
$$\label{rescale}$$

Next we rescale the right-hand sides of the BSSN equations so that the evolved variables are $\left\{h_{i j},a_{i j},\text{cf}, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}$


```python
# Step 4: Rescale the BSSN gauge RHS quantities so that the evolved
#         variables may remain smooth across coord singularities
vet_rhsU    = ixp.zerorank1()
bet_rhsU    = ixp.zerorank1()
for i in range(DIM):
    vet_rhsU[i]    =   beta_rhsU[i] / rfm.ReU[i]
    bet_rhsU[i]    =      B_rhsU[i] / rfm.ReU[i]
#print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
#print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
#print(betaU_dD)
#print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
#print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
```

<a id='code_validation'></a>

# Step 5: Code Validation against `BSSN.BSSN_gauge_RHSs` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN gauge equations between

1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_gauge_RHSs](../edit/BSSN/BSSN_gauge_RHSs.py) module.

By default, we analyze the RHSs in Spherical coordinates and with the covariant Gamma-driving second-order shift condition, though other coordinate systems & gauge conditions may be chosen.


```python
# Step 5: We already have SymPy expressions for BSSN gauge RHS expressions
#         in terms of other SymPy variables. Even if we reset the
#         list of NRPy+ gridfunctions, these *SymPy* expressions for
#         BSSN RHS variables *will remain unaffected*.
#
#         Here, we will use the above-defined BSSN gauge RHS expressions
#         to validate against the same expressions in the
#         BSSN/BSSN_gauge_RHSs.py file, to ensure consistency between
#         this tutorial and the module itself.
#
# Reset the list of gridfunctions, as registering a gridfunction
#   twice will spawn an error.
gri.glb_gridfcs_list = []


# Step 5.a: Call the BSSN_gauge_RHSs() function from within the
#           BSSN/BSSN_gauge_RHSs.py module,
#           which should generate exactly the same expressions as above.
import BSSN.BSSN_gauge_RHSs as Bgrhs
par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption","GammaDriving2ndOrder_Covariant")
Bgrhs.BSSN_gauge_RHSs()

print("Consistency check between BSSN.BSSN_gauge_RHSs tutorial and NRPy+ module: ALL SHOULD BE ZERO.")

print("alpha_rhs - bssnrhs.alpha_rhs = " + str(alpha_rhs - Bgrhs.alpha_rhs))

for i in range(DIM):
    print("vet_rhsU["+str(i)+"] - bssnrhs.vet_rhsU["+str(i)+"] = " + str(vet_rhsU[i] - Bgrhs.vet_rhsU[i]))
    print("bet_rhsU["+str(i)+"] - bssnrhs.bet_rhsU["+str(i)+"] = " + str(bet_rhsU[i] - Bgrhs.bet_rhsU[i]))
```

    Consistency check between BSSN.BSSN_gauge_RHSs tutorial and NRPy+ module: ALL SHOULD BE ZERO.
    alpha_rhs - bssnrhs.alpha_rhs = 0
    vet_rhsU[0] - bssnrhs.vet_rhsU[0] = 0
    bet_rhsU[0] - bssnrhs.bet_rhsU[0] = 0
    vet_rhsU[1] - bssnrhs.vet_rhsU[1] = 0
    bet_rhsU[1] - bssnrhs.bet_rhsU[1] = 0
    vet_rhsU[2] - bssnrhs.vet_rhsU[2] = 0
    bet_rhsU[2] - bssnrhs.bet_rhsU[2] = 0


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.pdf](Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs")
```

    Created Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.tex, and compiled
        LaTeX file to PDF file Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.pdf

