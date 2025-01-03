<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Equations of General Relativistic Force-Free Electrodynamics (GRFFE)

## Authors: Patrick Nelson, Zach Etienne, and Leo Werneck

## This notebook documents and constructs a number of quantities useful for building symbolic (SymPy) expressions for the equations of general relativistic force-free electrodynamics (GRFFE), using the same (Valencia) formalism as `IllinoisGRMHD`. The formulation further incorporates the computation of flux and source terms for the induction equation, with detailed code validation against the existing GRFFE equations NRPy+ module.

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

## Introduction

We write the equations of general relativistic hydrodynamics in conservative form as follows (adapted from Eqs. 41-44 of [Duez et al](https://arxiv.org/pdf/astro-ph/0503420.pdf)):

\begin{eqnarray}
\partial_t \tilde{S}_i &+& \partial_j \left(\alpha \sqrt{\gamma} T^j_{{\rm EM}i} \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu}_{\rm EM} g_{\mu\nu,i},
\end{eqnarray}
where we assume $T^{\mu\nu}_{\rm EM}$ is the electromagnetic stress-energy tensor:
$$
T^{\mu\nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu\nu} - b^\mu b^\nu,
$$
and 
$$
v^j = \frac{u^j}{u^0} \\
$$

Some quantities can be defined in precisely the same way they are defined in the GRHD equations. ***Therefore, we will not define special functions for generating these quantities, and instead refer the user to the appropriate functions in the [GRHD module](../edit/GRHD/equations.py)*** Namely,

* The GRFFE conservative variables:
    * $\tilde{S}_i  = \alpha \sqrt{\gamma} T^0{}_i$, via `GRHD.compute_S_tildeD(alpha, sqrtgammaDET, T4UD)`
* The GRFFE fluxes:
    * $\tilde{S}_i$ flux: $\left(\alpha \sqrt{\gamma} T^j{}_i \right)$, via `GRHD.compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD)`
* GRFFE source terms:
    * $\tilde{S}_i$ source term: $\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}$, via `GRHD.compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DD_zerotimederiv_dD, T4UU)`

Also, we will write the 4-metric in terms of the ADM 3-metric, lapse, and shift using standard equations.

Thus the full set of input variables includes:
* Spacetime quantities:
    * ADM quantities $\alpha$, $\beta^i$, $\gamma_{ij}$
* Hydrodynamical quantities:
    * 4-velocity $u^\mu$
* Electrodynamical quantities
    * Magnetic field $B^i$

Additionally, we need to evolve the vector $A_i$ according to
$$
\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j)
$$
and the scalar potential $\Phi$ according to
$$
\partial_t [\sqrt{\gamma} \Phi] + \partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) = - \xi \alpha [\sqrt{\gamma} \Phi],
$$
where $\xi$ determines the strength of the damping term. This is typically set to $0.1$.

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
1. [Step 2](#u0bu): Define needed quantities for $T^{\mu\nu}_{\rm EM}$, the EM part of the stress-energy tensor
1. [Step 3](#stressenergy): **compute_TEM4UU()**, **compute_TEM4UD()**: Define the stress-energy tensor $T^{\mu\nu}_{\rm EM}$ and $T^\mu_{{\rm EM}\nu}$
1. [Step 4](#inductioneq): Vector potential induction equation, assuming generalized Lorenz gauge
    1. [Step 4.a](#inductionterms) Compute the flux term and the source term for the induction equation
        1. **compute_AD_flux_term()**,**compute_AD_source_term_operand_for_FD()**
    1. [Step 4.b](#gaugeeq) Compute the damping term and flux term for the gauge equation
        1. **compute_psi6Phi_rhs_flux_term_operand()**,**compute_psi6Phi_rhs_damping_term()**
1. [Step 5](#declarevarsconstructgrffeeqs): Declare ADM and hydrodynamical input variables, and construct GRFFE equations
1. [Step 6](#code_validation): Code Validation against `GRFFE.equations` NRPy+ module for GRMHD
1. [Step 7](#code_validation_2): Code Validation against `GRFFE.equations` NRPy+ module for GRFFE
1. [Step 8](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='importmodules'></a>

# Step 1: Import needed NRPy+ & Python modules \[Back to [top](#toc)\]
$$\label{importmodules}$$


```python
# Step 1: Import needed core NRPy+ modules
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
```

<a id='u0bu'></a>

# Step 2: Define needed quantities for $T^{\mu\nu}_{\rm EM}$, the EM part of the stress-energy tensor \[Back to [top](#toc)\]
$$\label{u0bu}$$

We are given $B^i$, the magnetic field as measured by a *normal* observer, yet $T^{\mu\nu}_{\rm EM}$ depends on $b^{\mu}$, the magnetic field as measured by an observer comoving with the plasma $B^{\mu}_{\rm (u)}$, divided by $\sqrt{4\pi}$.

In the ideal MHD limit, $B^{\mu}_{\rm (u)}$ is orthogonal to the plasma 4-velocity $u^\mu$, which sets the $\mu=0$ component. 

$B^{\mu}_{\rm (u)}$ is related to the magnetic field as measured by a *normal* observer $B^i$ via a simple projection (Eq 21 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)), which results in the expressions (Eqs 23 and 24 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)):

\begin{align}
\sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
\sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}\\
\end{align}

Further, $B^i$ is related to the actual magnetic field evaluated in IllinoisGRMHD, $\tilde{B}^i$ via

$$B^i = \frac{\tilde{B}^i}{\gamma},$$

where $\gamma$ is the determinant of the spatial 3-metric:


```python
# Step 2.a: Define B^i = Btilde^i / sqrt(gamma)
def compute_B_notildeU(sqrtgammaDET, B_tildeU):
    global B_notildeU
    B_notildeU = ixp.zerorank1(DIM=3)
    for i in range(3):
        B_notildeU[i] = B_tildeU[i]/sqrtgammaDET
```

Next we compute Eqs 23 and 24 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf):

\begin{align}
\sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
\sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}.
\end{align}

In doing so, we will store the scalar $u_j B^j$ to `u4_dot_B_notilde`:


```python
# Step 2.b.i: Define b^mu.
def compute_smallb4U(gammaDD,betaU,alpha, u4U,B_notildeU, sqrt4pi):
    global smallb4U
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    u4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            u4D[mu] += AB4m.g4DD[mu][nu]*u4U[nu]
    smallb4U = ixp.zerorank1(DIM=4)
    u4_dot_B_notilde = sp.sympify(0)
    for i in range(3):
        u4_dot_B_notilde += u4D[i+1]*B_notildeU[i]

    # b^0 = (u_j B^j)/[alpha * sqrt(4 pi)]
    smallb4U[0] = u4_dot_B_notilde / (alpha*sqrt4pi)
    # b^i = [B^i + (u_j B^j) u^i]/[alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4U[i+1] = (B_notildeU[i] + u4_dot_B_notilde*u4U[i+1]) / (alpha*u4U[0]*sqrt4pi)
```

However, in the case of pure GRFFE (that is, if we assume the pressure and density are zero), we can make an additional simplifying assumption. When deriving the equations of GRFFE, one has a choice as to how to define the velocity. Following the example set by [McKinney (2006)](https://arxiv.org/pdf/astro-ph/0601410.pdf), section 2.1, and [Paschalidis (2013)](https://arxiv.org/pdf/1310.3274v2.pdf), Appendix A, we choose a referencing frame comoving with the plasma that fulfills the ideal MHD condition (i.e., the electric field is zero). However, in GRFFE, this choice yields a one-parameter family of timelike vectors. By choosing the drift velocity $v^i = u^i/u^0$ that minimizes the Lorentz factor, we find that the four-velocity is orthogonal to the magnetic field, or $u_j B^j = 0$. With this assumption, 

\begin{align}
\sqrt{4\pi} b^0 &= 0 \\
\sqrt{4\pi} b^i &= \frac{B^i}{\alpha u^0}.
\end{align}

This simplification also gives the inversion from $\tilde{S}_i$ to the drift velocity $v^i$ a closed form, greatly simplifying the conservative-to-primitive solver and removing the need to 


```python
# Step 2.b.ii: Define b^mu when u4 and B are orthogonal
def compute_smallb4U_with_driftvU_for_FFE(gammaDD,betaU,alpha, u4U,B_notildeU, sqrt4pi):
    global smallb4_with_driftv_for_FFE_U
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    u4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            u4D[mu] += AB4m.g4DD[mu][nu]*u4U[nu]
    smallb4_with_driftv_for_FFE_U = ixp.zerorank1(DIM=4)

    # b^0 = 0
    smallb4_with_driftv_for_FFE_U[0] = 0
    # b^i = B^i / [alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4_with_driftv_for_FFE_U[i+1] = B_notildeU[i] / (alpha*u4U[0]*sqrt4pi)
```

Finally we compute `smallbsquared`=$b^2 = b_{\mu} b^{\mu} = g_{\mu \nu} b^{\nu}b^{\mu}$:


```python
# Step 2.c: Define b^2.
def compute_smallbsquared(gammaDD,betaU,alpha, smallb4U):
    global smallbsquared
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    smallbsquared = sp.sympify(0)
    for mu in range(4):
        for nu in range(4):
            smallbsquared += AB4m.g4DD[mu][nu]*smallb4U[mu]*smallb4U[nu]
```

<a id='stressenergy'></a>

# Step 3: Define the electromagnetic stress-energy tensor $T^{\mu\nu}_{\rm EM}$ and $T^\mu_{{\rm EM}\nu}$ \[Back to [top](#toc)\]
$$\label{stressenergy}$$

Recall from above that

$$
T^{\mu\nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{1}{2} b^2 g^{\mu\nu} - b^\mu b^\nu.
$$
Also 

$$
T^\mu_{{\rm EM}\nu} = T^{\mu\delta}_{\rm EM} g_{\delta \nu}
$$


```python
# Step 3.a: Define T_{EM}^{mu nu} (a 4-dimensional tensor)
def compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallbsquared,u4U):
    global TEM4UU

    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    # Finally compute T^{mu nu}
    TEM4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            TEM4UU[mu][nu] = smallbsquared*u4U[mu]*u4U[nu] \
                             + sp.Rational(1,2)*smallbsquared*AB4m.g4UU[mu][nu] \
                             - smallb4U[mu]*smallb4U[nu]

# Step 3.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_TEM4UD(gammaDD,betaU,alpha, TEM4UU):
    global TEM4UD
    # Next compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
    # First we'll need g_{alpha nu} in terms of ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    TEM4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                TEM4UD[mu][nu] += TEM4UU[mu][delta]*AB4m.g4DD[delta][nu]
```

<a id='inductioneq'></a>

# Step 4: Vector potential induction equation, assuming generalized Lorenz gauge \[Back to [top](#toc)\]
$$\label{inductioneq}$$

Now, we will turn our attention to the induction equation, which evolves $A_i$, as well as an additional electromagnetic gauge evolution equation, which evolves $\Phi$. For a cell-centered formulation, they are as follows:
$$
\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j),
$$
from Eq. 17 of the [original `GiRaFFE` paper](https://arxiv.org/pdf/1704.00599.pdf), and 
$$
\partial_t [\sqrt{\gamma} \Phi] + \partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) = - \xi \alpha [\sqrt{\gamma} \Phi],
$$
from Eq. 19. When it comes to taking the derivatives in these equations, we will write these functions with the intention of computing the operand first, storing it as a gridfunction, and then finite differencing that in a later step. 

<a id='inductionterms'></a>

## Step 4.a: Compute the flux term and the source term for the induction equation \[Back to [top](#toc)\]
$$\label{inductionterms}$$

We'll now take a closer look at the induction equation, the right-hand side of which has two terms:
$$
\partial_t A_i = \underbrace{\epsilon_{ijk} v^j B^k}_{\rm "Flux"\ term} - \partial_i (\underbrace{\alpha \Phi - \beta^j A_j}_{\rm "Source"\ term}),
$$
The flux term here is a simple cross-product between the drift velocity and the magnetic field. The source term is the gradient of a gauge-dependent combination of the lapse $\alpha$, the shift $\beta^i$, the vector potential $A_i$, and the scalar potential $\Phi$; we will write a function to compute that operand, and save the finite-difference derivative until later. 


```python
def compute_AD_flux_term(sqrtgammaDET,driftvU,BU):
    # Levi-Civita tensor for cross products
    LeviCivitaDDD = ixp.LeviCivitaTensorDDD_dim3_rank3(sqrtgammaDET)
    global A_fluxD
    A_fluxD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # \epsilon_{ijk} v^j B^k
                A_fluxD[i] += LeviCivitaDDD[i][j][k]*driftvU[j]*BU[k]

def compute_AD_source_term_operand_for_FD(sqrtgammaDET,betaU,alpha,psi6Phi,AD):
    Phi = psi6Phi/sqrtgammaDET
    global AevolParen
    # \alpha \Phi - \beta^j A_j
    AevolParen = alpha * Phi
    for j in range(3):
        AevolParen += -betaU[j] * AD[j]

```

<a id='gaugeeq'></a>

## Step 4.b: Compute the damping term and flux term for the gauge equation \[Back to [top](#toc)\]
$$\label{gaugeeq}$$

Next, we will build the expressions for the RHS of the evolution equation 
$$
\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\underbrace{\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]}_{\rm Gauge\ evolution\ term}) - \underbrace{\xi \alpha [\sqrt{\gamma} \Phi]}_{\rm Damping\ term}.
$$
Once again, we will do this in two parts. First, we will compute the operand of the divergence in the flux term, leaving the finite-difference derivative for later; then, we will compute the damping term.


```python
def compute_psi6Phi_rhs_flux_term_operand(gammaDD,sqrtgammaDET,betaU,alpha,AD,psi6Phi):
    gammaUU,_gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    AU = ixp.zerorank1()
    # Raise the index on A in the usual way:
    for i in range(3):
        for j in range(3):
            AU[i] += gammaUU[i][j] * AD[j]

    global PhievolParenU
    PhievolParenU = ixp.zerorank1(DIM=3)

    for j in range(3):
        # \alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]
        PhievolParenU[j] += alpha*sqrtgammaDET*AU[j] - betaU[j]*psi6Phi

def compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping):
    # - \xi \alpha [\sqrt{\gamma} \Phi]
    # Combine the divergence and the damping term
    global psi6Phi_damping
    psi6Phi_damping = - xi_damping * alpha * psi6Phi

```

<a id='declarevarsconstructgrffeeqs'></a>

# Step 5: Declare ADM and hydrodynamical input variables, and construct GRFFE equations \[Back to [top](#toc)\]
$$\label{declarevarsconstructgrffeeqs}$$

The GRFFE equations are given by the induction equation (handled later) and the evolution equation for $\tilde{S}_i$, the Poynting one-form:

\begin{eqnarray}
\partial_t \tilde{S}_i &+& \partial_j \left(\alpha \sqrt{\gamma} T^j_{{\rm EM}i} \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu}_{\rm EM} g_{\mu\nu,i}.
\end{eqnarray}

Notice that all terms in this equation ($\tilde{S}_i$, $\left(\alpha \sqrt{\gamma} T^j_{{\rm EM}i} \right)$, and the source term) are *identical* to those in the evolution equation for $\tilde{S}_i$ in [general relativistic hydrodynamics (GRHD)](Tutorial-GRHD_Equations-Cartesian.ipynb); one need only replace $T^{\mu\nu}$ of GRHD with the $T^{\mu\nu}_{\rm EM}$ defined above. 

Thus we will reuse expressions from the [general relativistic hydrodynamics (GRHD)](Tutorial-GRHD_Equations-Cartesian.ipynb) module:


```python
# First define hydrodynamical quantities
u4U = ixp.declarerank1("u4U", DIM=4)
B_tildeU = ixp.declarerank1("B_tildeU", DIM=3)
AD       = ixp.declarerank1("D",        DIM=3)
B_tildeU = ixp.declarerank1("B_tildeU", DIM=3)
psi6Phi = sp.symbols('psi6Phi', real=True)

# Then ADM quantities
gammaDD = ixp.declarerank2("gammaDD","sym01",DIM=3)
betaU   = ixp.declarerank1("betaU", DIM=3)
alpha   = sp.symbols('alpha', real=True)

# Then numerical constant
sqrt4pi = sp.symbols('sqrt4pi', real=True)
xi_damping = sp.symbols('xi_damping', real=True)

# First compute stress-energy tensor T4UU and T4UD:
import GRHD.equations as GRHD
GRHD.compute_sqrtgammaDET(gammaDD)
compute_B_notildeU(GRHD.sqrtgammaDET, B_tildeU)
compute_smallb4U(gammaDD,betaU,alpha, u4U, B_notildeU, sqrt4pi)
compute_smallbsquared(gammaDD,betaU,alpha, smallb4U)

compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallbsquared,u4U)
compute_TEM4UD(gammaDD,betaU,alpha, TEM4UU)

# Compute conservative variables in terms of primitive variables
GRHD.compute_S_tildeD( alpha, GRHD.sqrtgammaDET, TEM4UD)
S_tildeD = GRHD.S_tildeD

# Next compute fluxes of conservative variables
GRHD.compute_S_tilde_fluxUD(alpha, GRHD.sqrtgammaDET,    TEM4UD)
S_tilde_fluxUD = GRHD.S_tilde_fluxUD

# Then declare derivatives & compute g4DDdD
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Finally compute source terms on tau_tilde and S_tilde equations
GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET, GRHD.g4DD_zerotimederiv_dD, TEM4UU)
S_tilde_source_termD = GRHD.S_tilde_source_termD

# We must also compute the terms for the induction equation and gauge evolution equation
GRHD.compute_vU_from_u4U__no_speed_limit(u4U) # We need the drift velocity
compute_AD_flux_term(GRHD.sqrtgammaDET,GRHD.vU,B_notildeU)
compute_AD_source_term_operand_for_FD(GRHD.sqrtgammaDET,betaU,alpha,psi6Phi,AD)
compute_psi6Phi_rhs_flux_term_operand(gammaDD,GRHD.sqrtgammaDET,betaU,alpha,AD,psi6Phi)
compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping)
```

<a id='code_validation'></a>

# Step 6: Code Validation against `GRFFE.equations` NRPy+ module for GRMHD \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the GRFFE equations generated in
1. this tutorial notebook versus
2. the NRPy+ [GRFFE.equations](../edit/GRFFE/equations.py) module.


```python
import GRFFE.equations as GRFFE

# First compute B^i from Btilde^i:
GRFFE.compute_B_notildeU(GRHD.sqrtgammaDET, B_tildeU)

# Then compute b^mu and b^2:
GRFFE.compute_smallb4U(gammaDD, betaU, alpha, u4U, GRFFE.B_notildeU, sqrt4pi)
GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)

# Next construct stress-energy tensor T4UU and T4UD:
GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, GRFFE.smallb4U, GRFFE.smallbsquared,u4U)
GRFFE.compute_TEM4UD(gammaDD,betaU,alpha, GRFFE.TEM4UU)

# Compute conservative variables in terms of primitive variables
GRHD.compute_S_tildeD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)
Ge_S_tildeD = GRHD.S_tildeD

# Next compute fluxes of conservative variables
GRHD.compute_S_tilde_fluxUD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)
Ge_S_tilde_fluxUD = GRHD.S_tilde_fluxUD

# Then declare derivatives & compute g4DDdD
# gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
# betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
# alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Finally compute source terms on S_tilde equations
GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET,GRHD.g4DD_zerotimederiv_dD,GRFFE.TEM4UU)
Ge_S_tilde_source_termD = GRHD.S_tilde_source_termD

# We must also compute the terms for the induction equation and gauge evolution equation
# GRHD.compute_vU_from_u4U__no_speed_limit(u4U) # We need the drift velocity
GRFFE.compute_AD_flux_term(GRHD.sqrtgammaDET,GRHD.vU,B_notildeU)
GRFFE.compute_AD_source_term_operand_for_FD(GRHD.sqrtgammaDET,betaU,alpha,psi6Phi,AD)
GRFFE.compute_psi6Phi_rhs_flux_term_operand(gammaDD,GRHD.sqrtgammaDET,betaU,alpha,AD,psi6Phi)
GRFFE.compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping)
```


```python
def comp_func(expr1,expr2,basename,prefixname2="GRFFE."):
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
# PhievolParenU
# psi6Phi_damping
# A_fluxD
# AevolParen
for mu in range(4):
    for nu in range(4):
        namecheck_list.extend([gfnm("TEM4UU",mu,nu),gfnm("TEM4UD",mu,nu)])
        exprcheck_list.extend([GRFFE.TEM4UU[mu][nu],GRFFE.TEM4UD[mu][nu]])
        expr_list.extend([TEM4UU[mu][nu],TEM4UD[mu][nu]])

for i in range(3):
    namecheck_list.extend([gfnm("S_tildeD",i),gfnm("S_tilde_source_termD",i),gfnm("A_fluxD",i),gfnm("PhievolParenU",i)])
    exprcheck_list.extend([Ge_S_tildeD[i],Ge_S_tilde_source_termD[i],GRFFE.A_fluxD[i],GRFFE.PhievolParenU[i]])
    expr_list.extend([S_tildeD[i],S_tilde_source_termD[i],A_fluxD[i],PhievolParenU[i]])
    for j in range(3):
        namecheck_list.extend([gfnm("S_tilde_fluxUD",i,j)])
        exprcheck_list.extend([Ge_S_tilde_fluxUD[i][j]])
        expr_list.extend([S_tilde_fluxUD[i][j]])

namecheck_list.extend(["AevolParen","psi6Phi_damping"])
exprcheck_list.extend([GRFFE.AevolParen,GRFFE.psi6Phi_damping])
expr_list.extend([AevolParen,psi6Phi_damping])

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


<a id='code_validation_2'></a>

# Step 7: Code Validation against `GRFFE.equations` NRPy+ module for GRFFE \[Back to [top](#toc)\]
$$\label{code_validation_2}$$

Additionally, we verify agreement in the SymPy expressions for the GRFFE equations generated in
1. this tutorial notebook versus
2. the NRPy+ [GRFFE.equations](../edit/GRFFE/equations.py) module
for the case $u_j B^j = 0$. 

We will only recompute $b^\mu$ and the expressions that depend on it.


```python
# Generate the expressions within the tutorial, starting with:
compute_smallb4U_with_driftvU_for_FFE(gammaDD,betaU,alpha, u4U, B_notildeU, sqrt4pi)
compute_smallbsquared(gammaDD,betaU,alpha, smallb4_with_driftv_for_FFE_U)

compute_TEM4UU(gammaDD,betaU,alpha, smallb4_with_driftv_for_FFE_U, smallbsquared,u4U)
compute_TEM4UD(gammaDD,betaU,alpha, TEM4UU)

# Compute conservative variables in terms of primitive variables
GRHD.compute_S_tildeD( alpha, GRHD.sqrtgammaDET, TEM4UD)
S_tildeD = GRHD.S_tildeD

# Next compute fluxes of conservative variables
GRHD.compute_S_tilde_fluxUD(alpha, GRHD.sqrtgammaDET,    TEM4UD)
S_tilde_fluxUD = GRHD.S_tilde_fluxUD

# Then declare derivatives & compute g4DDdD
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Finally compute source terms on tau_tilde and S_tilde equations
GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET, GRHD.g4DD_zerotimederiv_dD, TEM4UU)
S_tilde_source_termD = GRHD.S_tilde_source_termD

# Now compute the expressions from the module
# Then compute b^mu and b^2:
GRFFE.compute_smallb4U_with_driftvU_for_FFE(gammaDD, betaU, alpha, u4U, GRFFE.B_notildeU, sqrt4pi)
GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4_with_driftv_for_FFE_U)

# Next construct stress-energy tensor T4UU and T4UD:
GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, GRFFE.smallb4_with_driftv_for_FFE_U, GRFFE.smallbsquared,u4U)
GRFFE.compute_TEM4UD(gammaDD,betaU,alpha, GRFFE.TEM4UU)

# Compute conservative variables in terms of primitive variables
GRHD.compute_S_tildeD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)
Ge_S_tildeD = GRHD.S_tildeD

# Next compute fluxes of conservative variables
GRHD.compute_S_tilde_fluxUD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)
Ge_S_tilde_fluxUD = GRHD.S_tilde_fluxUD

# Then declare derivatives & compute g4DDdD
# gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
# betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
# alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Finally compute source terms on S_tilde equations
GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET,GRHD.g4DD_zerotimederiv_dD,GRFFE.TEM4UU)
Ge_S_tilde_source_termD = GRHD.S_tilde_source_termD
```


```python
all_passed=True

expr_list = []
exprcheck_list = []
namecheck_list = []

for mu in range(4):
    for nu in range(4):
        namecheck_list.extend([gfnm("TEM4UU",mu,nu),gfnm("TEM4UD",mu,nu)])
        exprcheck_list.extend([GRFFE.TEM4UU[mu][nu],GRFFE.TEM4UD[mu][nu]])
        expr_list.extend([TEM4UU[mu][nu],TEM4UD[mu][nu]])

for i in range(3):
    namecheck_list.extend([gfnm("S_tildeD",i),gfnm("S_tilde_source_termD",i)])
    exprcheck_list.extend([Ge_S_tildeD[i],Ge_S_tilde_source_termD[i]])
    expr_list.extend([S_tildeD[i],S_tilde_source_termD[i]])
    for j in range(3):
        namecheck_list.extend([gfnm("S_tilde_fluxUD",i,j)])
        exprcheck_list.extend([Ge_S_tilde_fluxUD[i][j]])
        expr_list.extend([S_tilde_fluxUD[i][j]])

for i in range(len(expr_list)):
    comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])

import sys
if all_passed:
    print("ALL TESTS PASSED!")
else:
    print("ERROR: AT LEAST ONE TEST DID NOT PASS")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 8: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-GRFFE_Equations-Cartesian.pdf](Tutorial-GRFFE_Equations-Cartesian.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-GRFFE_Equations-Cartesian")
```

    Created Tutorial-GRFFE_Equations-Cartesian.tex, and compiled LaTeX file to
        PDF file Tutorial-GRFFE_Equations-Cartesian.pdf

