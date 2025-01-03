<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Equations of General Relativistic Magnetohydrodynamics (GRMHD)

## Author: Zach Etienne

## This notebook documents and constructs a number of quantities useful for building symbolic (SymPy) expressions for the equations of general relativistic magnetohydrodynamics (GRMHD), using the same (Valencia-like) formalism as `IllinoisGRMHD`.

**Notebook Status:** <font color='orange'><b> Self-Validated; induction equation not yet implemented </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

## Introduction

We write the equations of general relativistic magnetohydrodynamics in conservative form as follows (Eqs. 41-44 of [Duez *et al*](https://arxiv.org/pdf/astro-ph/0503420.pdf):

\begin{eqnarray}
\ \partial_t \rho_* &+& \partial_j \left(\rho_* v^j\right) = 0 \\
\partial_t \tilde{\tau} &+& \partial_j \left(\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j \right) = s \\
\partial_t \tilde{S}_i &+& \partial_j \left(\alpha \sqrt{\gamma} T^j{}_i \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i} \\
\partial_t \tilde{B}^i &+& \partial_j \left(v^j \tilde{B}^i - v^i \tilde{B}^j\right) = 0,
\end{eqnarray}
where 

$$
s = \alpha \sqrt{\gamma}\left[\left(T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \right)K_{ij}
- \left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha \right].
$$

We represent $T^{\mu\nu}$ as the sum of the stress-energy tensor of a perfect fluid $T^{\mu\nu}_{\rm GRHD}$, plus the stress-energy associated with the electromagnetic fields in the force-free electrodynamics approximation $T^{\mu\nu}_{\rm GRFFE}$ (equivalently, $T^{\mu\nu}_{\rm em}$ in Duez *et al*):

$$
T^{\mu\nu} = T^{\mu\nu}_{\rm GRHD} + T^{\mu\nu}_{\rm GRFFE},
$$

where 

* $T^{\mu\nu}_{\rm GRHD}$ is constructed from rest-mass density $\rho_0$, pressure $P$, internal energy $\epsilon$, 4-velocity $u^\mu$, and ADM metric quantities as described in the [NRPy+ GRHD equations tutorial notebook](Tutorial-GRHD_Equations-Cartesian.ipynb); and 
* $T^{\mu\nu}_{\rm GRFFE}$ is constructed from the magnetic field vector $B^i$ and ADM metric quantities as described in the [NRPy+ GRFFE equations tutorial notebook](Tutorial-GRFFE_Equations-Cartesian.ipynb).

All quantities can be written in terms of the full GRMHD stress-energy tensor $T^{\mu\nu}$ in precisely the same way they are defined in the GRHD equations. ***Therefore, we will not define special functions for generating these quantities, and instead refer the user to the appropriate functions in the [GRHD module](../edit/GRHD/equations.py)*** Namely,

* The GRMHD conservative variables:
    * $\rho_* = \alpha\sqrt{\gamma} \rho_0 u^0$, via `GRHD.compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U)`
    * $\tilde{\tau} = \alpha^2\sqrt{\gamma} T^{00} - \rho_*$, via `GRHD.compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star)`
    * $\tilde{S}_i  = \alpha \sqrt{\gamma} T^0{}_i$, via `GRHD.compute_S_tildeD(alpha, sqrtgammaDET, T4UD)`
* The GRMHD fluxes:
    * $\rho_*$ flux: $\left(\rho_* v^j\right)$, via `GRHD.compute_rho_star_fluxU(vU, rho_star)`
    * $\tilde{\tau}$ flux: $\left(\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j \right)$, via `GRHD.compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU,T4UU,rho_star)`
    * $\tilde{S}_i$ flux: $\left(\alpha \sqrt{\gamma} T^j{}_i \right)$, via `GRHD.compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD)`
* GRMHD source terms:
    * $\tilde{\tau}$ source term $s$: defined above, via `GRHD.compute_s_source_term(KDD,betaU,alpha, sqrtgammaDET,alpha_dD, T4UU)` 
    * $\tilde{S}_i$ source term: $\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}$, via `GRHD.compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DD_zerotimederiv_dD, T4UU)`

In summary, all terms in the GRMHD equations can be constructed once the full GRMHD stress-energy tensor $T^{\mu\nu} = T^{\mu\nu}_{\rm GRHD} + T^{\mu\nu}_{\rm GRFFE}$ is constructed. For completeness, the full set of input variables include:
* Spacetime quantities:
    * ADM quantities $\alpha$, $\beta^i$, $\gamma_{ij}$, $K_{ij}$
* Hydrodynamical quantities:
    * Rest-mass density $\rho_0$
    * Pressure $P$
    * Internal energy $\epsilon$
    * 4-velocity $u^\mu$
* Electrodynamical quantities
    * Magnetic field $B^i= \tilde{B}^i / \gamma$

### A Note on Notation

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

For instance, in calculating the first term of $b^2 u^\mu u^\nu$, we use Greek indices:

```python
T4EMUU = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        # Term 1: b^2 u^{\mu} u^{\nu}
        T4EMUU[mu][nu] = smallb2*u4U[mu]*u4U[nu]
```

When we calculate $\beta_i = \gamma_{ij} \beta^j$, we use Latin indices:
```python
betaD = ixp.zerorank1(DIM=3)
for i in range(3):
    for j in range(3):
        betaD[i] += gammaDD[i][j] * betaU[j]
```

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook). This can be seen when we handle $\frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}$:
```python
# \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2
for i in range(3):
    for mu in range(4):
        for nu in range(4):
            S_tilde_rhsD[i] += alpsqrtgam * T4EMUU[mu][nu] * g4DD_zerotimederiv_dD[mu][nu][i+1] / 2
```

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

Each family of quantities is constructed within a given function (**boldfaced** below). This notebook is organized as follows


1. [Step 1](#importmodules): Import needed NRPy+ & Python modules
1. [Step 2](#stressenergy): Define the GRMHD stress-energy tensor $T^{\mu\nu}$ and $T^\mu{}_\nu$:
    * **compute_T4UU()**, **compute_T4UD()**: 
1. [Step 3](#declarevarsconstructgrhdeqs): Construct $T^{\mu\nu}$ from GRHD & GRFFE modules with ADM and GRMHD input variables, and construct GRMHD equations from the full GRMHD stress-energy tensor.
1. [Step 4](#code_validation): Code Validation against `GRMHD.equations` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='importmodules'></a>

# Step 1: Import needed NRPy+ & Python modules \[Back to [top](#toc)\]
$$\label{importmodules}$$


```python
# Step 1: Import needed core NRPy+ modules
import sympy as sp        # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
```

<a id='stressenergy'></a>

# Step 2: Define the GRMHD stress-energy tensor $T^{\mu\nu}$ and $T^\mu{}_\nu$ \[Back to [top](#toc)\]
$$\label{stressenergy}$$

Recall from above that

$$
T^{\mu\nu} = T^{\mu\nu}_{\rm GRHD} + T^{\mu\nu}_{\rm GRFFE},
$$
where

* $T^{\mu\nu}_{\rm GRHD}$ is constructed from the `GRHD.compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)` [GRHD](../edit/GRHD/equations.py) [(**tutorial**)](Tutorial-GRHD_Equations-Cartesian.ipynb) function, and
* $T^{\mu\nu}_{\rm GRFFE}$ is constructed from the `GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallbsquared,u4U)` [GRFFE](../edit/GRFFE/equations.py) [(**tutorial**)](Tutorial-GRFFE_Equations-Cartesian.ipynb) function

Since a lowering operation on a sum of tensors is equivalent to the lowering operation applied to the individual tensors in the sum,

$$
T^\mu{}_{\nu} = T^\mu{}_{\nu}{}^{\rm GRHD} + T^\mu{}_{\nu}{}^{\rm GRFFE},
$$

where

* $T^\mu{}_{\nu}{}^{\rm GRHD}$ is constructed from the `GRHD.compute_T4UD(gammaDD,betaU,alpha, T4UU)` [GRHD](../edit/GRHD/equations.py) [(**tutorial**)](Tutorial-GRHD_Equations-Cartesian.ipynb) function, and
* $T^{\mu\nu}_{\rm GRFFE}$ is constructed from the `GRFFE.compute_TEM4UD(gammaDD,betaU,alpha, TEM4UU)` [GRFFE](../edit/GRFFE/equations.py) [(**tutorial**)](Tutorial-GRFFE_Equations-Cartesian.ipynb) function.


```python
import GRHD.equations as GRHD
import GRFFE.equations as GRFFE

# Step 2.a: Define the GRMHD T^{mu nu} (a 4-dimensional tensor)
def compute_GRMHD_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U, smallb4U, smallbsquared):
    global GRHDT4UU
    global GRFFET4UU
    global T4UU

    GRHD.compute_T4UU(   gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)
    GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallbsquared,u4U)

    GRHDT4UU  = ixp.zerorank2(DIM=4)
    GRFFET4UU = ixp.zerorank2(DIM=4)
    T4UU      = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            GRHDT4UU[mu][nu]  = GRHD.T4UU[mu][nu]
            GRFFET4UU[mu][nu] =                     GRFFE.TEM4UU[mu][nu]
            T4UU[mu][nu]      = GRHD.T4UU[mu][nu] + GRFFE.TEM4UU[mu][nu]

# Step 2.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_GRMHD_T4UD(gammaDD,betaU,alpha, GRHDT4UU,GRFFET4UU):
    global T4UD

    GRHD.compute_T4UD(   gammaDD,betaU,alpha,  GRHDT4UU)
    GRFFE.compute_TEM4UD(gammaDD,betaU,alpha, GRFFET4UU)

    T4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UD[mu][nu] = GRHD.T4UD[mu][nu] + GRFFE.TEM4UD[mu][nu]
```

<a id='declarevarsconstructgrhdeqs'></a>

# Step 3: Declare ADM and hydrodynamical input variables, and construct all terms in GRMHD equations \[Back to [top](#toc)\]
$$\label{declarevarsconstructgrhdeqs}$$


```python
# First define hydrodynamical quantities
u4U = ixp.declarerank1("u4U", DIM=4)
rho_b,P,epsilon = sp.symbols('rho_b P epsilon',real=True)
B_tildeU = ixp.declarerank1("B_tildeU", DIM=3)

# Then ADM quantities
gammaDD = ixp.declarerank2("gammaDD","sym01",DIM=3)
KDD     = ixp.declarerank2("KDD"    ,"sym01",DIM=3)
betaU   = ixp.declarerank1("betaU", DIM=3)
alpha   = sp.symbols('alpha', real=True)

# Then numerical constant
sqrt4pi = sp.symbols('sqrt4pi', real=True)

# First compute smallb4U & smallbsquared from BtildeU, which are needed
#      for GRMHD stress-energy tensor T4UU and T4UD:
GRHD.compute_sqrtgammaDET(gammaDD)
GRFFE.compute_B_notildeU(GRHD.sqrtgammaDET,      B_tildeU)
GRFFE.compute_smallb4U(     gammaDD,betaU,alpha, u4U,GRFFE.B_notildeU, sqrt4pi)
GRFFE.compute_smallbsquared(gammaDD,betaU,alpha, GRFFE.smallb4U)

# Then compute the GRMHD stress-energy tensor:
compute_GRMHD_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U, GRFFE.smallb4U, GRFFE.smallbsquared)
compute_GRMHD_T4UD(gammaDD,betaU,alpha, GRHDT4UU,GRFFET4UU)

# Compute conservative variables in terms of primitive variables
GRHD.compute_rho_star( alpha, GRHD.sqrtgammaDET, rho_b,u4U)
GRHD.compute_tau_tilde(alpha, GRHD.sqrtgammaDET, T4UU,GRHD.rho_star)
GRHD.compute_S_tildeD( alpha, GRHD.sqrtgammaDET, T4UD)

# Then compute v^i from u^mu
GRHD.compute_vU_from_u4U__no_speed_limit(u4U)

# Next compute fluxes of conservative variables
GRHD.compute_rho_star_fluxU(                           GRHD.vU,     GRHD.rho_star)
GRHD.compute_tau_tilde_fluxU(alpha, GRHD.sqrtgammaDET, GRHD.vU,T4UU,GRHD.rho_star)
GRHD.compute_S_tilde_fluxUD( alpha, GRHD.sqrtgammaDET,         T4UD)

# Then declare derivatives & compute g4DD_zerotimederiv_dD
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Then compute source terms on tau_tilde and S_tilde equations
GRHD.compute_s_source_term(KDD,betaU,alpha, GRHD.sqrtgammaDET,alpha_dD,                   T4UU)
GRHD.compute_S_tilde_source_termD(   alpha, GRHD.sqrtgammaDET,GRHD.g4DD_zerotimederiv_dD, T4UU)
```

<a id='code_validation'></a>

# Step 4: Code Validation against `GRMHD.equations` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the GRHD equations generated in
1. this tutorial versus
2. the NRPy+ [GRMHD.equations](../edit/GRMHD/equations.py) module.


```python
import GRMHD.equations as GRMHD

# Compute stress-energy tensor T4UU and T4UD:
GRMHD.compute_GRMHD_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U, GRFFE.smallb4U, GRFFE.smallbsquared)
GRMHD.compute_GRMHD_T4UD(gammaDD,betaU,alpha, GRMHD.GRHDT4UU,GRMHD.GRFFET4UU)
```


```python
def comp_func(expr1,expr2,basename,prefixname2="Ge."):
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

for mu in range(4):
    for nu in range(4):
        namecheck_list.extend([gfnm("GRMHD.GRHDT4UU",mu,nu),gfnm("GRMHD.GRFFET4UU",mu,nu),
                               gfnm("GRMHD.T4UU",    mu,nu),gfnm("GRMHD.T4UD",     mu,nu)])
        exprcheck_list.extend([GRMHD.GRHDT4UU[mu][nu],GRMHD.GRFFET4UU[mu][nu],
                                   GRMHD.T4UU[mu][nu],     GRMHD.T4UD[mu][nu]])
        expr_list.extend([GRHDT4UU[mu][nu],GRFFET4UU[mu][nu],
                          T4UU[mu][nu],         T4UD[mu][nu]])

num_failures = 0
for i in range(len(expr_list)):
    num_failures += comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])

import sys
if num_failures == 0:
    print("ALL TESTS PASSED!")
else:
    print("ERROR: "+str(num_failures)+" TESTS DID NOT PASS")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-GRMHD_Equations-Cartesian.pdf](Tutorial-GRMHD_Equations-Cartesian.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-GRMHD_Equations-Cartesian")
```

    Created Tutorial-GRMHD_Equations-Cartesian.tex, and compiled LaTeX file to
        PDF file Tutorial-GRMHD_Equations-Cartesian.pdf

