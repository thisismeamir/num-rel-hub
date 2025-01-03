<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Equations of General Relativistic Hydrodynamics (GRHD)

## Authors: Zach Etienne & Patrick Nelson

## This notebook documents and constructs a number of quantities useful for building symbolic (SymPy) expressions for the equations of general relativistic hydrodynamics (GRHD), using the same (Valencia) formalism as `IllinoisGRMHD`.

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

## Introduction

We write the equations of general relativistic hydrodynamics in conservative form as follows (adapted from Eqs. 41-44 of [Duez et al](https://arxiv.org/pdf/astro-ph/0503420.pdf):

\begin{eqnarray}
\ \partial_t \rho_* &+& \partial_j \left(\rho_* v^j\right) = 0 \\
\partial_t \tilde{\tau} &+& \partial_j \left(\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j \right) = s \\
\partial_t \tilde{S}_i &+& \partial_j \left(\alpha \sqrt{\gamma} T^j{}_i \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i},
\end{eqnarray}
where we assume $T^{\mu\nu}$ is the stress-energy tensor of a perfect fluid:
$$
T^{\mu\nu} = \rho_0 h u^{\mu} u^{\nu} + P g^{\mu\nu},
$$
the $s$ source term is given in terms of ADM quantities via

$$
s = \alpha \sqrt{\gamma}\left[\left(T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \right)K_{ij}
- \left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha \right],
$$

and 
\begin{align}
v^j &= \frac{u^j}{u^0} \\
\rho_* &= \alpha\sqrt{\gamma} \rho_0 u^0 \\
h &= 1 + \epsilon + \frac{P}{\rho_0}.
\end{align}

Also, we will write the 4-metric in terms of the ADM 3-metric, lapse, and shift using standard equations.

Thus the full set of input variables includes:
* Spacetime quantities:
    * ADM quantities $\alpha$, $\beta^i$, $\gamma_{ij}$, $K_{ij}$
* Hydrodynamical quantities:
    * Rest-mass density $\rho_0$
    * Pressure $P$
    * Internal energy $\epsilon$
    * 4-velocity $u^\mu$

For completeness, the rest of the conservative variables are given by
\begin{align}
\tilde{\tau} &= \alpha^2\sqrt{\gamma} T^{00} - \rho_* \\
\tilde{S}_i  &= \alpha \sqrt{\gamma} T^0{}_i
\end{align}

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
1. [Step 2](#stressenergy): Define the stress-energy tensor $T^{\mu\nu}$ and $T^\mu{}_\nu$:
    * **compute_enthalpy()**, **compute_T4UU()**, **compute_T4UD()**: 
1. [Step 3](#primtoconserv): Writing the conservative variables in terms of the primitive variables:
    * **compute_sqrtgammaDET()**, **compute_rho_star()**, **compute_tau_tilde()**, **compute_S_tildeD()**
1. [Step 4](#grhdfluxes): Define the fluxes for the GRHD equations
    1. [Step 4.a](#rhostarfluxterm): Define $\rho_*$ flux term for GRHD equations:
        * **compute_vU_from_u4U__no_speed_limit()**, **compute_rho_star_fluxU()**:
    1. [Step 4.b](#taustildesourceterms) Define $\tilde{\tau}$ and $\tilde{S}_i$ flux terms for GRHD equations:
        * **compute_tau_tilde_fluxU()**, **compute_S_tilde_fluxUD()**
1. [Step 5](#grhdsourceterms): Define source terms on RHSs of GRHD equations
    1. [Step 5.a](#ssourceterm): Define $s$ source term on RHS of $\tilde{\tau}$ equation:
        * **compute_s_source_term()**
    1. [Step 5.b](#stildeisourceterm): Define source term on RHS of $\tilde{S}_i$ equation
        1. [Step 5.b.i](#fourmetricderivs): Compute $g_{\mu\nu,i}$ in terms of ADM quantities and their derivatives:
            * **compute_g4DD_zerotimederiv_dD()**
        1. [Step 5.b.ii](#stildeisource): Compute source term of the $\tilde{S}_i$ equation: $\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}$:
            * **compute_S_tilde_source_termD()**
1. [Step 6](#convertvtou): Conversion of $v^i$ to $u^\mu$ (Courtesy Patrick Nelson):
    * **u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit()**, **u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit()**
1. [Step 7](#declarevarsconstructgrhdeqs): Declare ADM and hydrodynamical input variables, and construct GRHD equations
1. [Step 8](#code_validation): Code Validation against `GRHD.equations` NRPy+ module
1. [Step 9](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='importmodules'></a>

# Step 1: Import needed NRPy+ & Python modules \[Back to [top](#toc)\]
$$\label{importmodules}$$


```python
# Step 1: Import needed core NRPy+ modules
from outputC import nrpyAbs      # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
```

<a id='stressenergy'></a>

# Step 2: Define the stress-energy tensor $T^{\mu\nu}$ and $T^\mu{}_\nu$ \[Back to [top](#toc)\]
$$\label{stressenergy}$$

Recall from above that

$$
T^{\mu\nu} = \rho_0 h u^{\mu} u^{\nu} + P g^{\mu\nu},
$$
where $h = 1 + \epsilon + \frac{P}{\rho_0}$. Also 

$$
T^\mu{}_{\nu} = T^{\mu\delta} g_{\delta \nu}
$$


```python
# Step 2.a: First define h, the enthalpy:
def compute_enthalpy(rho_b,P,epsilon):
    global h
    h = 1 + epsilon + P/rho_b

# Step 2.b: Define T^{mu nu} (a 4-dimensional tensor)
def compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U):
    global T4UU

    compute_enthalpy(rho_b,P,epsilon)
    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    # Finally compute T^{mu nu}
    T4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = rho_b * h * u4U[mu]*u4U[nu] + P*AB4m.g4UU[mu][nu]

# Step 2.c: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_T4UD(gammaDD,betaU,alpha, T4UU):
    global T4UD
    # Next compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
    # First we'll need g_{alpha nu} in terms of ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    T4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                T4UD[mu][nu] += T4UU[mu][delta]*AB4m.g4DD[delta][nu]
```

<a id='primtoconserv'></a>

# Step 3: Writing the conservative variables in terms of the primitive variables \[Back to [top](#toc)\]
$$\label{primtoconserv}$$

Recall from above that the conservative variables may be written as
\begin{align}
\rho_* &= \alpha\sqrt{\gamma} \rho_0 u^0 \\
\tilde{\tau} &= \alpha^2\sqrt{\gamma} T^{00} - \rho_* \\
\tilde{S}_i  &= \alpha \sqrt{\gamma} T^0{}_i
\end{align}

$T^{\mu\nu}$ and $T^\mu{}_\nu$ have already been defined $-$ all in terms of primitive variables. Thus we'll just need $\sqrt{\gamma}=$`gammaDET`, and all conservatives can then be written in terms of other defined quantities, which themselves are written in terms of primitive variables and the ADM metric.


```python
# Step 3: Writing the conservative variables in terms of the primitive variables
def compute_sqrtgammaDET(gammaDD):
    global sqrtgammaDET
    _gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    sqrtgammaDET = sp.sqrt(gammaDET)

def compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U):
    global rho_star
    # Compute rho_star:
    rho_star = alpha*sqrtgammaDET*rho_b*u4U[0]

def compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star):
    global tau_tilde
    tau_tilde = alpha**2*sqrtgammaDET*T4UU[0][0] - rho_star

def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    global S_tildeD
    S_tildeD = ixp.zerorank1(DIM=3)
    for i in range(3):
        S_tildeD[i] = alpha*sqrtgammaDET*T4UD[0][i+1]
```

<a id='grhdfluxes'></a>

# Step 4: Define the fluxes for the GRHD equations \[Back to [top](#toc)\]
$$\label{grhdfluxes}$$

<a id='rhostarfluxterm'></a>

## Step 4.a: Define $\rho_*$ flux term for GRHD equations \[Back to [top](#toc)\]
$$\label{rhostarfluxterm}$$

Recall from above that
\begin{array}
\ \partial_t \rho_* &+ \partial_j \left(\rho_* v^j\right) = 0.
\end{array}

Here we will define the $\rho_* v^j$ that constitutes the flux of $\rho_*$, first defining $v^j=u^j/u^0$:


```python
# Step 4: Define the fluxes for the GRHD equations
# Step 4.a: vU from u4U may be needed for computing rho_star_flux from u4U
def compute_vU_from_u4U__no_speed_limit(u4U):
    global vU
    # Now compute v^i = u^i/u^0:
    vU = ixp.zerorank1(DIM=3)
    for j in range(3):
        vU[j] = u4U[j+1]/u4U[0]

# Step 4.b: rho_star flux
def compute_rho_star_fluxU(vU, rho_star):
    global rho_star_fluxU
    rho_star_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        rho_star_fluxU[j] = rho_star*vU[j]
```

<a id='taustildesourceterms'></a>

## Step 4.b: Define $\tilde{\tau}$ and $\tilde{S}_i$ flux terms for GRHD equations \[Back to [top](#toc)\]
$$\label{taustildesourceterms}$$

Recall from above that
\begin{array}
\ \partial_t \tilde{\tau} &+ \partial_j \underbrace{\left(\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j \right)} &= s \\
\partial_t \tilde{S}_i &+ \partial_j \underbrace{\left(\alpha \sqrt{\gamma} T^j{}_i \right)} &= \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}.
\end{array}

Here we will define all terms that go inside the $\partial_j$'s on the left-hand side of the above equations (i.e., the underbraced expressions):


```python
# Step 4.c: tau_tilde flux
def compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU,T4UU, rho_star):
    global tau_tilde_fluxU
    tau_tilde_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        tau_tilde_fluxU[j] = alpha**2*sqrtgammaDET*T4UU[0][j+1] - rho_star*vU[j]

# Step 4.d: S_tilde flux
def compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD):
    global S_tilde_fluxUD
    S_tilde_fluxUD = ixp.zerorank2(DIM=3)
    for j in range(3):
        for i in range(3):
            S_tilde_fluxUD[j][i] = alpha*sqrtgammaDET*T4UD[j+1][i+1]
```

<a id='grhdsourceterms'></a>

# Step 5: Define source terms on RHSs of GRHD equations \[Back to [top](#toc)\]
$$\label{grhdsourceterms}$$

<a id='ssourceterm'></a>

## Step 5.a: Define $s$ source term on RHS of $\tilde{\tau}$ equation \[Back to [top](#toc)\]
$$\label{ssourceterm}$$


Recall again from above the $s$ source term on the right-hand side of the $\tilde{\tau}$ evolution equation is given in terms of ADM quantities and the stress-energy tensor via
$$
s = \underbrace{\alpha \sqrt{\gamma}}_{\text{Term 3}}\left[\underbrace{\left(T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \right)K_{ij}}_{\text{Term 1}}
\underbrace{- \left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha}_{\text{Term 2}} \right],
$$


```python
def compute_s_source_term(KDD,betaU,alpha, sqrtgammaDET,alpha_dD, T4UU):
    global s_source_term
    s_source_term = sp.sympify(0)
    # Term 1:
    for i in range(3):
        for j in range(3):
            s_source_term += (T4UU[0][0]*betaU[i]*betaU[j] + 2*T4UU[0][i+1]*betaU[j] + T4UU[i+1][j+1])*KDD[i][j]
    # Term 2:
    for i in range(3):
        s_source_term += -(T4UU[0][0]*betaU[i] + T4UU[0][i+1])*alpha_dD[i]
    # Term 3:
    s_source_term *= alpha*sqrtgammaDET
```

<a id='stildeisourceterm'></a>

## Step 5.b: Define source term on RHS of $\tilde{S}_i$ equation \[Back to [top](#toc)\]
$$\label{stildeisourceterm}$$

Recall from above
$$
\partial_t \tilde{S}_i + \partial_j \left(\alpha \sqrt{\gamma} T^j{}_i \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}.
$$
Our goal here will be to compute
$$
\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}.
$$

<a id='fourmetricderivs'></a>

### Step 5.b.i: Compute $g_{\mu\nu,i}$ in terms of ADM quantities and their derivatives \[Back to [top](#toc)\]
$$\label{fourmetricderivs}$$


To compute $g_{\mu\nu,i}$ we need to evaluate the first derivative of $g_{\mu\nu}$ in terms of ADM variables.

We are given $\gamma_{ij}$, $\alpha$, and $\beta^i$, and the 4-metric is given in terms of these quantities via
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$

Thus 
$$
g_{\mu\nu,k} = \begin{pmatrix} 
-2 \alpha\alpha_{,i} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k} & \beta_{i,k} \\
\beta_{j,k} & \gamma_{ij,k}
\end{pmatrix},
$$
where $\beta_{i} = \gamma_{ij} \beta^j$, so
$$
\beta_{i,k} = \gamma_{ij,k} \beta^j + \gamma_{ij} \beta^j_{,k}
$$


```python
def compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD):
    global g4DD_zerotimederiv_dD
    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]

    betaDdD = ixp.zerorank2(DIM=3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
                betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

    # Eq. 2.122 in B&S
    g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)
    for k in range(3):
        # Recall that g4DD[0][0] = -alpha^2 + betaU[j]*betaD[j]
        g4DD_zerotimederiv_dD[0][0][k+1] += -2*alpha*alpha_dD[k]
        for j in range(3):
            g4DD_zerotimederiv_dD[0][0][k+1] += betaU_dD[j][k]*betaD[j] + betaU[j]*betaDdD[j][k]

    for i in range(3):
        for k in range(3):
            # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
            g4DD_zerotimederiv_dD[i+1][0][k+1] = g4DD_zerotimederiv_dD[0][i+1][k+1] = betaDdD[i][k]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that g4DD[i][j] = gammaDD[i][j]
                g4DD_zerotimederiv_dD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]
```

<a id='stildeisource'></a>

### Step 5.b.ii: Compute source term of the $\tilde{S}_i$ equation: $\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}$ \[Back to [top](#toc)\]
$$\label{stildeisource}$$

Now that we've computed `g4DD_zerotimederiv_dD`$=g_{\mu\nu,i}$, the $\tilde{S}_i$ evolution equation source term may be quickly constructed.


```python
# Step 5.b.ii: Compute S_tilde source term
def compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DD_zerotimederiv_dD, T4UU):
    global S_tilde_source_termD
    S_tilde_source_termD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                S_tilde_source_termD[i] += sp.Rational(1,2)*alpha*sqrtgammaDET*T4UU[mu][nu]*g4DD_zerotimederiv_dD[mu][nu][i+1]
```

<a id='convertvtou'></a>

# Step 6: Conversion of $v^i$ to $u^\mu$ (Courtesy Patrick Nelson) \[Back to [top](#toc)\]
$$\label{convertvtou}$$

According to Eqs. 9-11 of [the IllinoisGRMHD paper](https://arxiv.org/pdf/1501.07276.pdf), the Valencia 3-velocity $v^i_{(n)}$ is related to the 4-velocity $u^\mu$ via

\begin{align}
\alpha v^i_{(n)} &= \frac{u^i}{u^0} + \beta^i \\
\implies u^i &= u^0 \left(\alpha v^i_{(n)} - \beta^i\right)
\end{align}

Defining $v^i = \frac{u^i}{u^0}$, we get

$$v^i = \alpha v^i_{(n)} - \beta^i,$$

and in terms of this variable we get

\begin{align}
g_{00} \left(u^0\right)^2 + 2 g_{0i} u^0 u^i + g_{ij} u^i u^j &= \left(u^0\right)^2 \left(g_{00} + 2 g_{0i} v^i + g_{ij} v^i v^j\right)\\
\implies u^0 &= \pm \sqrt{\frac{-1}{g_{00} + 2 g_{0i} v^i + g_{ij} v^i v^j}} \\
&= \pm \sqrt{\frac{-1}{(-\alpha^2 + \beta^2) + 2 \beta_i v^i + \gamma_{ij} v^i v^j}} \\
&= \pm \sqrt{\frac{1}{\alpha^2 - \gamma_{ij}\left(\beta^i + v^i\right)\left(\beta^j + v^j\right)}}\\
&= \pm \sqrt{\frac{1}{\alpha^2 - \alpha^2 \gamma_{ij}v^i_{(n)}v^j_{(n)}}}\\
&= \pm \frac{1}{\alpha}\sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}
\end{align}

Generally speaking, numerical errors will occasionally drive expressions under the radical to either negative values or potentially enormous values (corresponding to enormous Lorentz factors). Thus a reliable approach for computing $u^0$ requires that we first rewrite the above expression in terms of the Lorentz factor squared: $\Gamma^2=\left(\alpha u^0\right)^2$:
\begin{align}
u^0 &= \pm \frac{1}{\alpha}\sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}\\
\implies \left(\alpha u^0\right)^2 &= \frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}} \\
\implies \gamma_{ij}v^i_{(n)}v^j_{(n)} &= 1 - \frac{1}{\left(\alpha u^0\right)^2} \\
&= 1 - \frac{1}{\Gamma^2}
\end{align}

In order for the bottom expression to hold true, the left-hand side must be between 0 and 1. Again, this is not guaranteed due to the appearance of numerical errors. In fact, a robust algorithm will not allow $\Gamma^2$ to become too large (which might contribute greatly to the stress-energy of a given gridpoint), so let's define the largest allowed Lorentz factor as $\Gamma_{\rm max}$.

Then our algorithm for computing $u^0$ is as follows:

If
$$R=\gamma_{ij}v^i_{(n)}v^j_{(n)}>1 - \frac{1}{\Gamma_{\rm max}^2},$$ 
then adjust the 3-velocity $v^i$ as follows:

$$v^i_{(n)} \to \sqrt{\frac{1 - \frac{1}{\Gamma_{\rm max}^2}}{R}}v^i_{(n)}.$$

After this rescaling, we are then guaranteed that if $R$ is recomputed, it will be set to its ceiling value $R=R_{\rm max} = 1 - \frac{1}{\Gamma_{\rm max}^2}$.

Then, regardless of whether the ceiling on $R$ was applied, $u^0$ can be safely computed via

$$
u^0 = \frac{1}{\alpha \sqrt{1-R}},
$$
and the remaining components $u^i$ via
$$
u^i = u^0 v^i.
$$

In summary our algorithm for computing $u^{\mu}$ from $v^i = \frac{u^i}{u^0}$ is as follows:

1. Choose a maximum Lorentz factor $\Gamma_{\rm max}$=`GAMMA_SPEED_LIMIT`, and define $v^i_{(n)} = \frac{1}{\alpha}\left( \frac{u^i}{u^0} + \beta^i\right)$.
1. Compute $R=\gamma_{ij}v^i_{(n)}v^j_{(n)}=1 - \frac{1}{\Gamma^2}$
1. If $R \le 1 - \frac{1}{\Gamma_{\rm max}^2}$, then skip the next step.
1. Otherwise if $R > 1 - \frac{1}{\Gamma_{\rm max}^2}$ then adjust $v^i_{(n)}\to \sqrt{\frac{1 - \frac{1}{\Gamma_{\rm max}^2}}{R}}v^i_{(n)}$, which will force $R=R_{\rm max}$.
1. Given the $R$ computed in the above step, $u^0 = \frac{1}{\alpha \sqrt{1-R}}$, and $u^i=u^0 v^i$.

While the above algorithm is quite robust, its `if()` statement in the fourth step is not very friendly to NRPy+ or an optimizing C compiler, as it would require NRPy+ to generate separate C kernels for each branch of the `if()`. Let's instead try the following trick, which Roland Haas taught us. 

Define $R^*$ as

$$
R^* = \frac{1}{2} \left(R_{\rm max} + R - |R_{\rm max} - R| \right).
$$

If $R>R_{\rm max}$, then $|R_{\rm max} - R|=R - R_{\rm max}$, and we get:

$$
R^* = \frac{1}{2} \left(R_{\rm max} + R - (R - R_{\rm max}) \right) = \frac{1}{2} \left(2 R_{\rm max}\right) = R_{\rm max}
$$

If $R\le R_{\rm max}$, then $|R_{\rm max} - R|=R_{\rm max} - R$, and we get:

$$
R^* = \frac{1}{2} \left(R_{\rm max} + R - (R_{\rm max} - R) \right) = \frac{1}{2} \left(2 R\right) = R
$$

Then we can rescale *all* $v^i_{(n)}$ via

$$
v^i_{(n)} \to v^i_{(n)} \sqrt{\frac{R^*}{R}},
$$

though we must be very careful to carefully handle the case in which $R=0$. To avoid any problems in this case, we simply adjust the above rescaling by adding a tiny number [`TINYDOUBLE`](https://en.wikipedia.org/wiki/Tiny_Bubbles) to $R$ in the denominator, typically `1e-100`:

$$
v^i_{(n)} \to v^i_{(n)} \sqrt{\frac{R^*}{R + {\rm TINYDOUBLE}}}.
$$

Finally, $u^0$ can be immediately and safely computed, via:
$$
u^0 = \frac{1}{\alpha \sqrt{1-R^*}},
$$
and $u^i$ via 
$$
u^i = u^0 v^i = u^0 \left(\alpha v^i_{(n)} - \beta^i\right).
$$


```python
# Step 6.a: Convert Valencia 3-velocity v_{(n)}^i into u^\mu, and apply a speed limiter
#           Speed-limited ValenciavU is output to rescaledValenciavU global.
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, ValenciavU):
    # Inputs:  Metric lapse alpha, shift betaU, 3-metric gammaDD, Valencia 3-velocity ValenciavU
    # Outputs (as globals): u4U_ito_ValenciavU, rescaledValenciavU

    # R = gamma_{ij} v^i v^j
    R = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    thismodule = "GRHD"
    # The default value isn't terribly important here, since we can overwrite in the main C code
    GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)  # Default value based on
                                                                                        # IllinoisGRMHD.
    # GiRaFFE default = 2000.0
    Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
    # Now, we set Rstar = min(Rmax,R):
    # If R <  Rmax, then Rstar = 0.5*(Rmax+R-Rmax+R) = R
    # If R >= Rmax, then Rstar = 0.5*(Rmax+R+Rmax-R) = Rmax
    Rstar = sp.Rational(1, 2) * (Rmax + R - nrpyAbs(Rmax - R))

    # We add TINYDOUBLE to R below to avoid a 0/0, which occurs when
    #    ValenciavU == 0 for all Valencia 3-velocity components.
    # "Those tiny *doubles* make me warm all over
    #  with a feeling that I'm gonna love you till the end of time."
    #    - Adapted from Connie Francis' "Tiny Bubbles"
    TINYDOUBLE = par.Cparameters("#define",thismodule,"TINYDOUBLE",1e-100)

    # The rescaled (speed-limited) Valencia 3-velocity
    #     is given by, v_{(n)}^i = sqrt{Rstar/R} v^i
    global rescaledValenciavU
    rescaledValenciavU = ixp.zerorank1(DIM=3)
    for i in range(3):
        # If R == 0, then Rstar == 0, so sqrt( Rstar/(R+TINYDOUBLE) )=sqrt(0/1e-100) = 0
        #   If your velocities are of order 1e-100 and this is physically
        #   meaningful, there must be something wrong with your unit conversion.
        rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rstar/(R + TINYDOUBLE))

    # Finally compute u^mu in terms of Valenciav^i
    # u^0 = 1/(alpha-sqrt(1-R^*))
    global u4U_ito_ValenciavU
    u4U_ito_ValenciavU = ixp.zerorank1(DIM=4)
    u4U_ito_ValenciavU[0] = 1/(alpha*sp.sqrt(1-Rstar))
    # u^i = u^0 ( alpha v^i_{(n)} - beta^i ), where v^i_{(n)} is the Valencia 3-velocity
    for i in range(3):
        u4U_ito_ValenciavU[i+1] = u4U_ito_ValenciavU[0] * (alpha * rescaledValenciavU[i] - betaU[i])

# Step 6.b: Convert v^i into u^\mu, and apply a speed limiter.
#           Speed-limited vU is output to rescaledvU global.
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha,betaU,gammaDD, vU):
    ValenciavU = ixp.zerorank1(DIM=3)
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i])/alpha
    u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, ValenciavU)
    # Since ValenciavU is written in terms of vU,
    #   u4U_ito_ValenciavU is actually u4U_ito_vU
    global u4U_ito_vU
    u4U_ito_vU = ixp.zerorank1(DIM=4)
    for mu in range(4):
        u4U_ito_vU[mu] = u4U_ito_ValenciavU[mu]
    # Finally compute the rescaled (speed-limited) vU
    global rescaledvU
    rescaledvU = ixp.zerorank1(DIM=3)
    for i in range(3):
        rescaledvU[i] = alpha * rescaledValenciavU[i] - betaU[i]
```

<a id='declarevarsconstructgrhdeqs'></a>

# Step 7: Declare ADM and hydrodynamical input variables, and construct GRHD equations \[Back to [top](#toc)\]
$$\label{declarevarsconstructgrhdeqs}$$


```python
# First define hydrodynamical quantities
u4U = ixp.declarerank1("u4U", DIM=4)
rho_b,P,epsilon = sp.symbols('rho_b P epsilon',real=True)

# Then ADM quantities
gammaDD = ixp.declarerank2("gammaDD","sym01",DIM=3)
KDD     = ixp.declarerank2("KDD"    ,"sym01",DIM=3)
betaU   = ixp.declarerank1("betaU", DIM=3)
alpha   = sp.symbols('alpha', real=True)

# First compute stress-energy tensor T4UU and T4UD:
compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)
compute_T4UD(gammaDD,betaU,alpha, T4UU)

# Next sqrt(gamma)
compute_sqrtgammaDET(gammaDD)

# Compute conservative variables in terms of primitive variables
compute_rho_star( alpha, sqrtgammaDET, rho_b,u4U)
compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star)
compute_S_tildeD( alpha, sqrtgammaDET, T4UD)

# Then compute v^i from u^mu
compute_vU_from_u4U__no_speed_limit(u4U)

# Next compute fluxes of conservative variables
compute_rho_star_fluxU(                                   vU,      rho_star)
compute_tau_tilde_fluxU(alpha,              sqrtgammaDET, vU,T4UU, rho_star)
compute_S_tilde_fluxUD( alpha,              sqrtgammaDET,    T4UD)

# Then declare derivatives & compute g4DD_zerotimederiv_dD
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Then compute source terms on tau_tilde and S_tilde equations
compute_s_source_term(KDD,betaU,alpha, sqrtgammaDET,alpha_dD, T4UU)
compute_S_tilde_source_termD(    alpha, sqrtgammaDET,g4DD_zerotimederiv_dD,   T4UU)

# Then compute the 4-velocities in terms of an input Valencia 3-velocity testValenciavU[i]
testValenciavU = ixp.declarerank1("testValenciavU",DIM=3)
u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, testValenciavU)

# Finally compute the 4-velocities in terms of an input 3-velocity testvU[i] = u^i/u^0
testvU = ixp.declarerank1("testvU",DIM=3)
u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha,betaU,gammaDD, testvU)
```

<a id='code_validation'></a>

# Step 8: Code Validation against `GRHD.equations` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the GRHD equations generated in
1. this tutorial versus
2. the NRPy+ [GRHD.equations](../edit/GRHD/equations.py) module.


```python
import GRHD.equations as Ge

# First compute stress-energy tensor T4UU and T4UD:
Ge.compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)
Ge.compute_T4UD(gammaDD,betaU,alpha, Ge.T4UU)

# Next sqrt(gamma)
Ge.compute_sqrtgammaDET(gammaDD)

# Compute conservative variables in terms of primitive variables
Ge.compute_rho_star( alpha, Ge.sqrtgammaDET, rho_b,u4U)
Ge.compute_tau_tilde(alpha, Ge.sqrtgammaDET, Ge.T4UU,Ge.rho_star)
Ge.compute_S_tildeD( alpha, Ge.sqrtgammaDET, Ge.T4UD)

# Then compute v^i from u^mu
Ge.compute_vU_from_u4U__no_speed_limit(u4U)

# Next compute fluxes of conservative variables
Ge.compute_rho_star_fluxU (                                     Ge.vU,         Ge.rho_star)
Ge.compute_tau_tilde_fluxU(alpha,              Ge.sqrtgammaDET, Ge.vU,Ge.T4UU, Ge.rho_star)
Ge.compute_S_tilde_fluxUD (alpha,              Ge.sqrtgammaDET,       Ge.T4UD)

# Then declare derivatives & compute g4DD_zerotimederiv_dD
# gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
# betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
# alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
Ge.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Finally compute source terms on tau_tilde and S_tilde equations
Ge.compute_s_source_term(KDD,betaU,alpha, Ge.sqrtgammaDET,alpha_dD, Ge.T4UU)
Ge.compute_S_tilde_source_termD(   alpha, Ge.sqrtgammaDET,Ge.g4DD_zerotimederiv_dD,Ge.T4UU)

GetestValenciavU = ixp.declarerank1("testValenciavU")
Ge.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, GetestValenciavU)

GetestvU = ixp.declarerank1("testvU")
Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(                alpha, betaU, gammaDD, GetestvU)
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

namecheck_list.extend(["sqrtgammaDET","rho_star","tau_tilde","s_source_term"])
exprcheck_list.extend([Ge.sqrtgammaDET,Ge.rho_star,Ge.tau_tilde,Ge.s_source_term])
expr_list.extend([sqrtgammaDET,rho_star,tau_tilde,s_source_term])
for mu in range(4):
    namecheck_list.extend([gfnm("u4_ito_ValenciavU",mu),gfnm("u4U_ito_vU",mu)])
    exprcheck_list.extend([   Ge.u4U_ito_ValenciavU[mu],   Ge.u4U_ito_vU[mu]])
    expr_list.extend(        [   u4U_ito_ValenciavU[mu],      u4U_ito_vU[mu]])
    for nu in range(4):
        namecheck_list.extend([gfnm("T4UU",mu,nu),gfnm("T4UD",mu,nu)])
        exprcheck_list.extend([Ge.T4UU[mu][nu],Ge.T4UD[mu][nu]])
        expr_list.extend([T4UU[mu][nu],T4UD[mu][nu]])
        for delta in range(4):
            namecheck_list.extend([gfnm("g4DD_zerotimederiv_dD",mu,nu,delta)])
            exprcheck_list.extend([Ge.g4DD_zerotimederiv_dD[mu][nu][delta]])
            expr_list.extend([g4DD_zerotimederiv_dD[mu][nu][delta]])


for i in range(3):
    namecheck_list.extend([gfnm("S_tildeD",i),gfnm("vU",i),gfnm("rho_star_fluxU",i),
                           gfnm("tau_tilde_fluxU",i),gfnm("S_tilde_source_termD",i),
                           gfnm("rescaledValenciavU",i), gfnm("rescaledvU",i)])
    exprcheck_list.extend([Ge.S_tildeD[i],Ge.vU[i],Ge.rho_star_fluxU[i],
                           Ge.tau_tilde_fluxU[i],Ge.S_tilde_source_termD[i],
                           Ge.rescaledValenciavU[i],Ge.rescaledvU[i]])
    expr_list.extend([S_tildeD[i],vU[i],rho_star_fluxU[i],
                      tau_tilde_fluxU[i],S_tilde_source_termD[i],
                      rescaledValenciavU[i],rescaledvU[i]])
    for j in range(3):
        namecheck_list.extend([gfnm("S_tilde_fluxUD",i,j)])
        exprcheck_list.extend([Ge.S_tilde_fluxUD[i][j]])
        expr_list.extend([S_tilde_fluxUD[i][j]])

num_failures = 0
for i in range(len(expr_list)):
    num_failures += comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])

import sys
if num_failures == 0:
    print("ALL TESTS PASSED!")
else:
    print("ERROR: " + str(num_failures) + " TESTS DID NOT PASS")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 9: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-GRHD_Equations-Cartesian.pdf](Tutorial-GRHD_Equations-Cartesian.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-GRHD_Equations-Cartesian")
```

    Created Tutorial-GRHD_Equations-Cartesian.tex, and compiled LaTeX file to
        PDF file Tutorial-GRHD_Equations-Cartesian.pdf

