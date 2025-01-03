<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Fishbone-Moncrief Initial Data

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook outlines the implementation of Fishbone-Moncrief initial data for GRMHD simulations in the Einstein Toolkit (ETK) using the NRPy+ module. It transforms the spherical coordinate data into Cartesian coordinates and constructs key physics quantities like four-velocity, Kerr-Schild metric, magnetic field, and Lorentz factor. Finally, the module outputs these calculations to C code.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** The expressions in this notebook have been validated against trusted versions, as part of the [Event Horizon Telescope GRMHD code comparison project](https://arxiv.org/abs/1904.04923) ([see this tutorial notebook for the analysis](Tutorial-Start_to_Finish-FishboneMoncriefID_standalone.ipynb)), and research performed as part of the TCAN project [80NSSC18K1488](https://compact-binaries.org/research/area/tcan). Also, this tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**


### NRPy+ Source Code for this module: [FishboneMoncriefID/FishboneMoncriefID.py](../edit/FishboneMoncriefID/FishboneMoncriefID.py)

## Introduction:
The goal of this module will be to construct Fishbone-Moncrief initial data for GRMHD simulations in a format suitable for the Einstein Toolkit (ETK). We will be using the equations as derived in [the original paper](http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1976ApJ...207..962F&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf), which will hereafter be called "***the FM paper***". Since we want to use this with the ETK, our final result will be in Cartesian coordinates. The natural coordinate system for these data is spherical, however, so we will use [reference_metric.py](../edit/reference_metric.py) ([**Tutorial**](Tutorial-Reference_Metric.ipynb)) to help with the coordinate transformation.

This notebook documents the equations in the NRPy+ module [FishboneMoncrief.py](../edit/FishboneMoncriefID/FishboneMoncriefID.py).  Then, we will build an Einstein Toolkit [thorn](Tutorial-ETK_thorn-FishboneMoncriefID.ipynb) to set this initial data.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#fishbonemoncrief): Implementing Fishbone-Moncrief initial data within NRPy+
    1. [Step 2.a](#registergridfunctions): Register within NRPy+ needed gridfunctions and initial parameters
    1. [Step 2.b](#l_of_r): Specific angular momentum $l(r)$
    1. [Step 2.c](#enthalpy): Specific enthalpy $h$
    1. [Step 2.d](#pressure_density): Pressure and density, from the specific enthalpy
    1. [Step 2.e](#covariant_velocity): Nonzero covariant velocity components $u_\mu$
    1. [Step 2.f](#inverse_bl_metric): Inverse metric $g^{\mu\nu}$ for the black hole in Boyer-Lindquist coordinates
    1. [Step 2.g](#xform_to_ks): Transform components of four-veloicty $u^\mu$ to Kerr-Schild
    1. [Step 2.h](#ks_metric): Define Kerr-Schild metric $g_{\mu\nu}$ and extrinsic curvature $K_{ij}$
    1. [Step 2.i](#magnetic_field): Seed poloidal magnetic field $B^i$
    1. [Step 2.j](#adm_metric): Set the ADM quantities $\alpha$, $\beta^i$, and $\gamma_{ij}$ from the spacetime metric $g_{\mu\nu}$
    1. [Step 2.k](#magnetic_field_comoving_frame): Set the magnetic field components in the comoving frame $b^\mu$, and $b^2$, which is twice the magnetic pressure
    1. [Step 2.l](#lorentz_fac_valencia): Lorentz factor $\Gamma = \alpha u^0$ and Valencia 3-velocity $v^i_{(n)}$ 
1. [Step 3](#output_to_c): Output SymPy expressions to C code, using NRPy+

1. [Step 4](#code_validation): Code Validation against Code Validation against `FishboneMoncriefID.FishboneMoncriefID` NRPy+ module NRPy+ module 
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$
 
We begin by importing the packages and NRPy+ modules that we will need. We will also set some of the most commonly used parameters. 


```python
# Step 1a: Import needed NRPy+ core modules:
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm   # NRPy+: Reference metric support

par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()
#Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

thismodule = "FishboneMoncriefID"
```

<a id='fishbonemoncrief'></a>

# Step 2: The Fishbone-Moncrief Initial Data Prescription \[Back to [top](#toc)\]
$$\label{fishbonemoncrief}$$

With NRPy's most important functions now available to us, we can start to set up the rest of the tools we will need to build the initial data. 

<a id='registergridfunctions'></a>

## Step 2.a: Register within NRPy+ needed gridfunctions and initial parameters \[Back to [top](#toc)\]
$$\label{registergridfunctions}$$

We will now register the gridfunctions we expect to use. Critically, we register the physical metric and extrinsic curvature tensors. 


```python
gPhys4UU = ixp.register_gridfunctions_for_single_rank2("AUX","gPhys4UU", "sym01", DIM=4)
KDD = ixp.register_gridfunctions_for_single_rank2("EVOL","KDD", "sym01")

# Variables needed for initial data given in spherical basis
r, th, ph = gri.register_gridfunctions("AUX",["r","th","ph"])

r_in,r_at_max_density,a,M = par.Cparameters("REAL",thismodule,
                                            ["r_in","r_at_max_density",    "a","M"],
                                            [   6.0,              12.0, 0.9375,1.0])

kappa,gamma = par.Cparameters("REAL",thismodule,["kappa","gamma"], [1.0e-3, 4.0/3.0])

# The return value from gri.register_gridfunctions("AUX","LorentzFactor") is unused, so we ignore it here:
gri.register_gridfunctions("AUX","LorentzFactor")
```




$\displaystyle LorentzFactor$



<a id='l_of_r'></a>

## Step 2.b: Specific angular momentum $l(r)$ \[Back to [top](#toc)\]
$$\label{l_of_r}$$

Now, we can begin actually building the ID equations. We will start with the value of the angular momentum $l$ at the position $r \equiv$`r_at_max_density` where the density is at a maximum, as in equation 3.8 of the FM paper:
\begin{align}
l(r) &= \pm \left( \frac{M}{r^3} \right) ^{1/2}
    \left[ \frac{r^4+r^2a^2-2Mra^2 \mp a(Mr)^{1/2}(r^2-a^2)}
    {r^2 -3Mr \pm 2a(Mr)^{1/2}} \right].
\end{align}


```python
def calculate_l_at_r(r):
    l  = sp.sqrt(M/r**3) * (r**4 + r**2*a**2 - 2*M*r*a**2 - a*sp.sqrt(M*r)*(r**2-a**2))
    l /= r**2 - 3*M*r + 2*a*sp.sqrt(M*r)
    return l

# First compute angular momentum at r_at_max_density, TAKING POSITIVE ROOT. This way disk is co-rotating with black hole
# Eq 3.8:
l = calculate_l_at_r(r_at_max_density)
```

<a id='enthalpy'></a>

## Step 2.c: Specific enthalpy $h$ \[Back to [top](#toc)\]
$$\label{enthalpy}$$


Next, we will follow equation 3.6 of the FM paper to compute the enthalpy $h$ by first finding its logarithm $\ln h$. Fortunately, we can make this process quite a bit simpler by first identifying the common subexpressions. Let 
\begin{align}
\Delta &= r^2 - 2Mr + a^2 \\
\Sigma &= r^2 + a^2 \cos^2 (\theta) \\
A &= (r^2+a^2)^2 - \Delta a^2 \sin^2(\theta);
\end{align}
furthermore, let 
\begin{align}
\text{tmp3} &= \sqrt{\frac{1 + 4 l^2 \Sigma^2 \Delta}{A \sin^2 (\theta)}}. \\
\end{align}
(These terms reflect the radially-independent part of the log of the enthalpy, `ln_h_const`.)
So, 
$$
{\rm ln\_h\_const} = \frac{1}{2} * \log \left( \frac{1+\text{tmp3}}{\Sigma \Delta/A} \right) - \frac{1}{2} \text{tmp3} - \frac{2aMrl}{A}
$$


```python
# Eq 3.6:
# First compute the radially-independent part of the log of the enthalpy, ln_h_const
Delta = r**2 - 2*M*r + a**2
Sigma = r**2 + a**2*sp.cos(th)**2
A = (r**2 + a**2)**2 - Delta*a**2*sp.sin(th)**2

# Next compute the radially-dependent part of log(enthalpy), ln_h
tmp3 = sp.sqrt(1 + 4*l**2*Sigma**2*Delta/(A*sp.sin(th))**2)
# Term 1 of Eq 3.6
ln_h  = sp.Rational(1,2)*sp.log( ( 1 + tmp3) / (Sigma*Delta/A))
# Term 2 of Eq 3.6
ln_h -= sp.Rational(1,2)*tmp3
# Term 3 of Eq 3.6
ln_h -= 2*a*M*r*l/A

```

Additionally, let
\begin{align}
\Delta_{\rm in} &= r_{\rm in}^2 - 2Mr_{\rm in} + a^2 \\
\Sigma_{\rm in} &= r_{\rm in}^2 + a^2 \cos^2 (\pi/2) \\
A_{\rm in} &= (r_{\rm in}^2+a^2)^2 - \Delta_{\rm in} a^2 \sin^2(\pi/2)
\end{align}
and 
\begin{align}
\text{tmp3in} &= \sqrt{\frac{1 + 4 l^2 \Sigma_{\rm in}^2 \Delta_{\rm in}}{A_{\rm in} \sin^2 (\theta)}}, \\
\end{align}
corresponding to the radially Independent part of log(enthalpy), $\ln h$:
\begin{align}
{\rm mln\_h\_in} = -\frac{1}{2} * \log \left( \frac{1+\text{tmp3in}}{\Sigma_{\rm in} \Delta_{\rm in}/A_{\rm in}} \right) + \frac{1}{2} \text{tmp3in} + \frac{2aMr_{\rm in}l}{A_{\rm in}}. \\
\end{align}
(Note that there is some typo in the expression for these terms given in Eq 3.6, so we opt to just evaluate the negative of the first three terms at r=`r_in` and th=pi/2 (the integration constant), as described in the text below Eq. 3.6.)

So, then, we exponentiate:
\begin{align}
\text{hm1} \equiv h-1 &= e^{{\rm ln\_h}+{\rm mln\_h\_in}}-1. \\
\end{align}



```python
# Next compute the radially-INdependent part of log(enthalpy), ln_h
# Note that there is some typo in the expression for these terms given in Eq 3.6, so we opt to just evaluate
#   negative of the first three terms at r=r_in and th=pi/2 (the integration constant), as described in
#   the text below Eq. 3.6, basically just copying the above lines of code.
# Delin = Delta_in ; Sigin = Sigma_in ; Ain = A_in .
Delin = r_in**2 - 2*M*r_in + a**2
Sigin = r_in**2 + a**2*sp.cos(sp.pi/2)**2
Ain   = (r_in**2 + a**2)**2 - Delin*a**2*sp.sin(sp.pi/2)**2

tmp3in = sp.sqrt(1 + 4*l**2*Sigin**2*Delin/(Ain*sp.sin(sp.pi/2))**2)
# Term 4 of Eq 3.6
mln_h_in  = -sp.Rational(1,2)*sp.log( ( 1 + tmp3in) / (Sigin*Delin/Ain))
# Term 5 of Eq 3.6
mln_h_in += sp.Rational(1,2)*tmp3in
# Term 6 of Eq 3.6
mln_h_in += 2*a*M*r_in*l/Ain

hm1 = sp.exp(ln_h + mln_h_in) - 1
```

<a id='pressure_density'></a>

## Step 2.d: Pressure and density, from the specific enthalpy \[Back to [top](#toc)\]
$$\label{pressure_density}$$

Python 3.4 + SymPy 1.0.0 has a serious problem taking the power here; it hangs forever, so instead we use the identity $x^{1/y} = \exp(\frac{1}{y} * \log(x))$. Thus, our expression for density becomes (in Python 2.7 + SymPy 0.7.4.1):

\begin{align}
\rho_0 &= \left( \frac{(h-1)(\gamma-1)}{\kappa \gamma} \right)^{1/(\gamma-1)} \\
&= \exp \left[ {\frac{1}{\gamma-1} \log \left( \frac{(h-1)(\gamma-1)}{\kappa \gamma}\right)} \right]
\end{align}

Additionally, the pressure $P_0 = \kappa \rho_0^\gamma$


```python
rho_initial,Pressure_initial = gri.register_gridfunctions("AUX",["rho_initial","Pressure_initial"])

# Python 3.4 + sympy 1.0.0 has a serious problem taking the power here, hangs forever.
# so instead we use the identity x^{1/y} = exp( [1/y] * log(x) )
# Original expression (works with Python 2.7 + sympy 0.7.4.1):
# rho_initial = ( hm1*(gamma-1)/(kappa*gamma) )**(1/(gamma - 1))
# New expression (workaround):
rho_initial = sp.exp( (1/(gamma-1)) * sp.log( hm1*(gamma-1)/(kappa*gamma) ))
Pressure_initial = kappa * rho_initial**gamma

```

<a id='covariant_velocity'></a>

## Step 2.e: Nonzero covariant velocity components $u_\mu$ \[Back to [top](#toc)\]
$$\label{covariant_velocity}$$

We now want to compute eq 3.3; we will start by finding $e^{-2 \chi}$ in Boyer-Lindquist (BL) coordinates. By eq 2.16, $\chi = \psi - \nu$, so, by eqs. 3.5,
\begin{align}
e^{2 \nu} &= \frac{\Sigma \Delta}{A} \\
e^{2 \psi} &= \frac{A \sin^2 \theta}{\Sigma} \\
e^{-2 \chi} &= e^{2 \nu} / e^{2 \psi} = e^{2(\nu - \psi)}.
\end{align}

Next, we will calculate the 4-velocity $u_i$ of the fluid disk in BL coordinates. We start with eqs. 3.3 and 2.13, finding
\begin{align}
u_{(r)} = u_{(\theta)} &= 0 \\
u_{(\phi)} &= \sqrt{-1+ \frac{1}{2}\sqrt{1 + 4l^2e^{-2 \chi}}} \\
u_{(t)} &= - \sqrt{1 + u_{(\phi)}^2}.
\end{align}

Given that $\omega = 2aMr/A$, we then find that, in BL coordinates,
\begin{align}
u_r = u_{\theta} &= 0 \\
u_{\phi} &= u_{(\phi)} \sqrt{e^{2 \psi}} \\
u_t &= u_{(t)} \sqrt{e^{2 \nu}} - \omega u_{\phi},
\end{align}
using eq. 2.13 to get the last relation.


```python
# Eq 3.3: First compute exp(-2 chi), assuming Boyer-Lindquist coordinates
#    Eq 2.16: chi = psi - nu, so
#    Eq 3.5 -> exp(-2 chi) = exp(-2 (psi - nu)) = exp(2 nu)/exp(2 psi)
exp2nu  = Sigma*Delta / A
exp2psi = A*sp.sin(th)**2 / Sigma
expm2chi = exp2nu / exp2psi

# Eq 3.3: Next compute u_(phi).
u_pphip = sp.sqrt((-1 + sp.sqrt(1 + 4*l**2*expm2chi))/2)
# Eq 2.13: Compute u_(t)
u_ptp   = -sp.sqrt(1 + u_pphip**2)

# Next compute spatial components of 4-velocity in Boyer-Lindquist coordinates:
uBL4D = ixp.zerorank1(DIM=4) # Components 1 and 2: u_r = u_theta = 0
# Eq 2.12 (typo): u_(phi) = e^(-psi) u_phi -> u_phi = e^(psi) u_(phi)
uBL4D[3] = sp.sqrt(exp2psi)*u_pphip

# Assumes Boyer-Lindquist coordinates:
omega = 2*a*M*r/A
# Eq 2.13: u_(t) = 1/sqrt(exp2nu) * ( u_t + omega*u_phi )
#     -->  u_t = u_(t) * sqrt(exp2nu) - omega*u_phi
#     -->  u_t = u_ptp * sqrt(exp2nu) - omega*uBL4D[3]
uBL4D[0] = u_ptp*sp.sqrt(exp2nu) - omega*uBL4D[3]
```

<a id='inverse_bl_metric'></a>

## Step 2.f: Inverse metric $g^{\mu\nu}$ for the black hole in Boyer-Lindquist coordinates \[Back to [top](#toc)\]
$$\label{inverse_bl_metric}$$

Next, we will use eq. 2.1 to find the inverse physical (as opposed to conformal) metric in BL coordinates, using the shorthands defined in eq. 3.5:
\begin{align}
g_{tt}                  &= - \frac{\Sigma \Delta}{A} + \omega^2 \sin^2 \theta \frac{A}{\Sigma} \\
g_{t \phi} = g_{\phi t} &= - \omega \sin^2 \theta \frac{A}{\Sigma} \\
g_{\phi \phi}           &= \sin^2 \theta \frac{A}{\Sigma},
\end{align}
which can be inverted to show that 
\begin{align}
g^{tt}                  &= - \frac{A}{\Delta \Sigma} \\
g^{t \phi} = g^{\phi t} &= \frac{2aMr}{\Delta \Sigma} \\
g^{\phi \phi}           &= - \frac{4a^2M^2r^2}{\Delta A \Sigma} + \frac{\Sigma^2}{A \Sigma \sin^2 \theta}.
\end{align}

With this, we will now be able to raise the index on the BL $u_i$: $u^i = g^{ij} u_j$


```python
# Eq. 3.5:
# w = 2*a*M*r/A;
# Eqs. 3.5 & 2.1:
# gtt = -Sig*Del/A + w^2*Sin[th]^2*A/Sig;
# gtp = w*Sin[th]^2*A/Sig;
# gpp = Sin[th]^2*A/Sig;
# FullSimplify[Inverse[{{gtt,gtp},{gtp,gpp}}]]
gPhys4BLUU = ixp.zerorank2(DIM=4)
gPhys4BLUU[0][0] = -A/(Delta*Sigma)
# DO NOT NEED TO SET gPhys4BLUU[1][1] or gPhys4BLUU[2][2]!
gPhys4BLUU[0][3] = gPhys4BLUU[3][0] = -2*a*M*r/(Delta*Sigma)
gPhys4BLUU[3][3] = -4*a**2*M**2*r**2/(Delta*A*Sigma) + Sigma**2/(A*Sigma*sp.sin(th)**2)

uBL4U = ixp.zerorank1(DIM=4)
for i in range(4):
    for j in range(4):
        uBL4U[i] += gPhys4BLUU[i][j]*uBL4D[j]
```

<a id='xform_to_ks'></a>

## Step 2.g: Transform components of four-velocity $u^\mu$ to Kerr-Schild \[Back to [top](#toc)\]
$$\label{xform_to_ks}$$

Now, we will transform the 4-velocity from the Boyer-Lindquist to the Kerr-Schild basis. This algorithm is adapted from [HARM](https://github.com/atchekho/harmpi/blob/master/init.c). This definees the tensor `transformBLtoKS`, where the diagonal elements are $1$, and the non-zero off-diagonal elements are 
\begin{align}
\text{transformBLtoKS}_{tr} &= \frac{2r}{r^2-2r+a^2} \\
\text{transformBLtoKS}_{\phi r} &=  \frac{a}{r^2-2r+a^2} \\
\end{align}


```python
# https://github.com/atchekho/harmpi/blob/master/init.c
# Next transform Boyer-Lindquist velocity to Kerr-Schild basis:
transformBLtoKS = ixp.zerorank2(DIM=4)
for i in range(4):
    transformBLtoKS[i][i] = 1
transformBLtoKS[0][1] = 2*r/(r**2 - 2*r + a*a)
transformBLtoKS[3][1] =   a/(r**2 - 2*r + a*a)
#uBL4U = ixp.declarerank1("UBL4U",DIM=4)
# After the xform below, print(uKS4U) outputs:
# [UBL4U0 + 2*UBL4U1*r/(a**2 + r**2 - 2*r), UBL4U1, UBL4U2, UBL4U1*a/(a**2 + r**2 - 2*r) + UBL4U3]
uKS4U = ixp.zerorank1(DIM=4)
for i in range(4):
    for j in range(4):
        uKS4U[i] += transformBLtoKS[i][j]*uBL4U[j]
```

<a id='ks_metric'></a>

## Step 2.h: Define Kerr-Schild metric $g_{\mu\nu}$ and extrinsic curvature $K_{ij}$ \[Back to [top](#toc)\]
$$\label{ks_metric}$$

We will also adopt the Kerr-Schild metric for Fishbone-Moncrief disks. Further details can be found in [Cook's Living Review](http://gravity.psu.edu/numrel/jclub/jc/Cook___LivRev_2000-5.pdf) article on initial data, or in the appendix of [this](https://arxiv.org/pdf/1704.00599.pdf) article. So, in KS coordinates,
\begin{align}
\rho^2 &= r^2 + a^2 \cos^2 \theta \\
\Delta &= r^2 - 2Mr + a^2 \\
\alpha &= \left(1 + \frac{2Mr}{\rho^2}\right)^{-1/2} \\
\beta^0 &= \frac{2 \alpha^2 Mr}{\rho^2} \\
\gamma_{00} &= 1 + \frac{2Mr}{\rho^2} \\
\gamma_{02} = \gamma_{20} &= -\left(1+\frac{2Mr}{\rho^2}\right) a \sin^2 \theta \\
\gamma_{11} &= \rho^2 \\
\gamma_{22} &= \left(r^2+a^2+\frac{2Mr}{\rho^2} a^2 \sin^2 \theta\right) \sin^2 \theta.
\end{align}
(Note that only the non-zero components of $\beta^i$ and $\gamma_{ij}$ are defined here.)


```python
# Adopt the Kerr-Schild metric for Fishbone-Moncrief disks
# http://gravity.psu.edu/numrel/jclub/jc/Cook___LivRev_2000-5.pdf
# Alternatively, Appendix of https://arxiv.org/pdf/1704.00599.pdf
rhoKS2  = r**2 + a**2*sp.cos(th)**2 # Eq 79 of Cook's Living Review article
DeltaKS = r**2 - 2*M*r + a**2    # Eq 79 of Cook's Living Review article
alphaKS = 1/sp.sqrt(1 + 2*M*r/rhoKS2)
betaKSU = ixp.zerorank1()
betaKSU[0] = alphaKS**2*2*M*r/rhoKS2
gammaKSDD = ixp.zerorank2()
gammaKSDD[0][0] = 1 + 2*M*r/rhoKS2
gammaKSDD[0][2] = gammaKSDD[2][0] = -(1 + 2*M*r/rhoKS2)*a*sp.sin(th)**2
gammaKSDD[1][1] = rhoKS2
gammaKSDD[2][2] = (r**2 + a**2 + 2*M*r/rhoKS2 * a**2*sp.sin(th)**2) * sp.sin(th)**2
```

We can also define the following useful quantities, continuing in KS coordinates:
\begin{align}
A &= a^2 \cos (2 \theta) + a^2 +2r^2 \\
B &= A + 4Mr \\
D &= \sqrt{\frac{2Mr}{a^2 \cos^2 \theta +r^2}+1};
\end{align}
we will also define the extrinsic curvature:
\begin{align}
K_{00}          &= D\frac{A+2Mr}{A^2 B} (4M(a^2 \cos(2 \theta)+a^s-2r^2)) \\
K_{01} = K_{10} &= \frac{D}{AB} (8a^2Mr\sin \theta \cos \theta) \\
K_{02} = K_{20} &= \frac{D}{A^2} (-2aM \sin^2 \theta (a^2\cos(2 \theta)+a^2-2r^2)) \\
K_{11}          &= \frac{D}{B} (4Mr^2) \\
K_{12} = K_{21} &= \frac{D}{AB} (-8a^3Mr \sin^3 \theta \cos \theta) \\
K_{22}          &= \frac{D}{A^2 B} (2Mr \sin^2 \theta (a^4(r-M) \cos(4 \theta) + a^4 (M+3r) + 4a^2 r^2 (2r-M) + 4a^2 r \cos(2 \theta) (a^2 + r(M+2r)) + 8r^5)). \\
\end{align}
Note that the indexing for extrinsic curvature only runs from 0 to 2, since there are no time components to the tensor. 


```python
AA = a**2 * sp.cos(2*th) + a**2 + 2*r**2
BB = AA + 4*M*r
DD = sp.sqrt(2*M*r / (a**2 * sp.cos(th)**2 + r**2) + 1)
KDD[0][0] = DD*(AA + 2*M*r)/(AA**2*BB) * (4*M*(a**2 * sp.cos(2*th) + a**2 - 2*r**2))
KDD[0][1] = KDD[1][0] = DD/(AA*BB) * 8*a**2*M*r*sp.sin(th)*sp.cos(th)
KDD[0][2] = KDD[2][0] = DD/AA**2 * (-2*a*M*sp.sin(th)**2 * (a**2 * sp.cos(2*th) + a**2 - 2*r**2))
KDD[1][1] = DD/BB * 4*M*r**2
KDD[1][2] = KDD[2][1] = DD/(AA*BB) * (-8*a**3*M*r*sp.sin(th)**3*sp.cos(th))
KDD[2][2] = DD/(AA**2*BB) * \
            (2*M*r*sp.sin(th)**2 * (a**4*(r-M)*sp.cos(4*th) + a**4*(M+3*r) +
             4*a**2*r**2*(2*r-M) + 4*a**2*r*sp.cos(2*th)*(a**2 + r*(M+2*r)) + 8*r**5))
```

We must also compute the inverse and determinant of the KS metric. We can use the NRPy+ [indexedexp.py](../edit/indexedexp.py) function to do this easily for the inverse physical 3-metric $\gamma^{ij}$, and then use the lapse $\alpha$ and the shift $\beta^i$ to find the full, inverse 4-dimensional metric, $g^{ij}$. We use the general form relating the 3- and 4- metric from (B&S 2.122)
\begin{equation}
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta\cdot\beta & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix},
\end{equation}
and invert it. That is,
\begin{align}
g^{00}          &= -\frac{1}{\alpha^2} \\
g^{0i} = g^{i0} &= \frac{\beta^{i-1}}{\alpha^2} \\
g^{ij} = g^{ji} &= \gamma^{(i-1) (j-1)} - \frac{\beta^{i-1} \beta^{j-1}}{\alpha^2},
\end{align}
keeping careful track of the differences in the indexing conventions for 3-dimensional quantities and 4-dimensional quantities (Python always indexes lists from 0, but in four dimensions, the 0 direction corresponds to time, while in 3+1, the connection to time is handled by other variables).


```python
# For compatibility, we must compute gPhys4UU
gammaKSUU,gammaKSDET = ixp.symm_matrix_inverter3x3(gammaKSDD)
# See, e.g., Eq. 4.49 of https://arxiv.org/pdf/gr-qc/0703035.pdf , where N = alpha
gPhys4UU[0][0] = -1 / alphaKS**2
for i in range(1,4):
    if i>0:
        # if the quantity does not have a "4", then it is assumed to be a 3D quantity.
        #  E.g., betaKSU[] is a spatial vector, with indices ranging from 0 to 2:
        gPhys4UU[0][i] = gPhys4UU[i][0] = betaKSU[i-1]/alphaKS**2
for i in range(1,4):
    for j in range(1,4):
        # if the quantity does not have a "4", then it is assumed to be a 3D quantity.
        #  E.g., betaKSU[] is a spatial vector, with indices ranging from 0 to 2,
        #    and gammaKSUU[][] is a spatial tensor, with indices again ranging from 0 to 2.
        gPhys4UU[i][j] = gPhys4UU[j][i] = gammaKSUU[i-1][j-1] - betaKSU[i-1]*betaKSU[j-1]/alphaKS**2
```

<a id='magnetic_field'></a>

## Step 2.i: Seed poloidal magnetic field $B^i$ \[Back to [top](#toc)\]
$$\label{magnetic_field}$$

The original Fishbone-Moncrief initial data prescription describes a non-self-gravitating accretion disk in hydrodynamical equilibrium about a black hole. The following assumes that a very weak magnetic field seeded into this disk will not significantly disturb this equilibrium, at least on a dynamical (free-fall) timescale.

Now, we will set up the magnetic field that, when simulated with a GRMHD code, will give us insight into the electromagnetic emission from the disk. We define the vector potential $A_i$ to be proportional to $\rho_0$, and, as usual, let the magnetic field $B^i$ be the curl of the vector potential.


```python
A_b = par.Cparameters("REAL",thismodule,"A_b",1.0)

A_3vecpotentialD = ixp.zerorank1()
# Set A_phi = A_b*rho_initial FIXME: why is there a sign error?
A_3vecpotentialD[2] = -A_b * rho_initial

BtildeU = ixp.register_gridfunctions_for_single_rank1("EVOL","BtildeU")
# Eq 15 of https://arxiv.org/pdf/1501.07276.pdf:
# B = curl A -> B^r = d_th A_ph - d_ph A_th
BtildeU[0] = sp.diff(A_3vecpotentialD[2],th) - sp.diff(A_3vecpotentialD[1],ph)
# B = curl A -> B^th = d_ph A_r - d_r A_ph
BtildeU[1] = sp.diff(A_3vecpotentialD[0],ph) - sp.diff(A_3vecpotentialD[2],r)
# B = curl A -> B^ph = d_r A_th - d_th A_r
BtildeU[2] = sp.diff(A_3vecpotentialD[1],r)  - sp.diff(A_3vecpotentialD[0],th)
```

<a id='adm_metric'></a>

## Step 2.j: Set the ADM quantities $\alpha$, $\beta^i$, and $\gamma_{ij}$ from the spacetime metric $g_{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{adm_metric}$$

Now, we wish to build the 3+1-dimensional variables in terms of the inverse 4-dimensional spacetime metric $g^{ij},$ as demonstrated in eq. 4.49 of [Gourgoulhon's lecture notes on 3+1 formalisms](https://arxiv.org/pdf/gr-qc/0703035.pdf) (letting $N=\alpha$). So,
\begin{align}
\alpha &= \sqrt{-\frac{1}{g^{00}}} \\
\beta^i &= \alpha^2 g^{0 (i+1)} \\
\gamma^{ij} &= g^{(i+1) (j+1)} + \frac{\beta^i \beta^j}{\alpha^2},
\end{align}
again keeping careful track of the differences in the indexing conventions for 3-dimensional quantities and 4-dimensional quantities. We will also take the inverse of $\gamma^{ij}$, obtaining (naturally) $\gamma_{ij}$ and its determinant $|\gamma|$. (Note that the function we use gives the determinant of $\gamma^{ij}$, which is the reciprocal of $|\gamma|$.)


```python
# Construct spacetime metric in 3+1 form:
# See, e.g., Eq. 4.49 of https://arxiv.org/pdf/gr-qc/0703035.pdf , where N = alpha
# The return values from gri.register_gridfunctions() & ixp.register_gridfunctions_for_single_rank1() are
#   unused, so we ignore them below:
gri.register_gridfunctions("EVOL",["alpha"])
ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")

alpha = sp.sqrt(1/(-gPhys4UU[0][0]))
betaU = ixp.zerorank1()
for i in range(3):
    betaU[i] = alpha**2 * gPhys4UU[0][i+1]
gammaUU = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        gammaUU[i][j] = gPhys4UU[i+1][j+1] + betaU[i]*betaU[j]/alpha**2

# The return value from ixp.register_gridfunctions_for_single_rank2() is unused so we ignore it below:
ixp.register_gridfunctions_for_single_rank2("EVOL","gammaDD","sym01")
gammaDD,igammaDET = ixp.symm_matrix_inverter3x3(gammaUU)
gammaDET = 1/igammaDET
```

Now, we will lower the index on the shift vector $\beta_j = \gamma_{ij} \beta^i$ and use that to calculate the 4-dimensional metric tensor, $g_{ij}$. So, we have 
\begin{align}
g_{00}                    &= -\alpha^2 + \beta^2 \\
g_{0 (i+1)} = g_{(i+1) 0} &= \beta_i \\
g_{(i+1) (j+1)}           &= \gamma_{ij},
\end{align}
where $\beta^2 \equiv \beta^i \beta_i$.


```python
###############
# Next compute g_{\alpha \beta} from lower 3-metric, using
# Eq 4.47 of https://arxiv.org/pdf/gr-qc/0703035.pdf
betaD = ixp.zerorank1()
for i in range(3):
    for j in range(3):
        betaD[i] += gammaDD[i][j]*betaU[j]

beta2 = sp.sympify(0)
for i in range(3):
    beta2 += betaU[i]*betaD[i]

gPhys4DD = ixp.zerorank2(DIM=4)
gPhys4DD[0][0] = -alpha**2 + beta2
for i in range(3):
    gPhys4DD[0][i+1] = gPhys4DD[i+1][0] = betaD[i]
    for j in range(3):
        gPhys4DD[i+1][j+1] = gammaDD[i][j]
```

<a id='magnetic_field_comoving_frame'></a>

## Step 2.k: Set the magnetic field components in the comoving frame $b^\mu$, and $b^2$, which is twice the magnetic pressure \[Back to [top](#toc)\]
$$\label{magnetic_field_comoving_frame}$$

Next compute $b^{\mu}$ using Eqs 23, 24, 27 and 31 of [this paper](https://arxiv.org/pdf/astro-ph/0503420.pdf):
\begin{align}
B^i &= \frac{\tilde{B}}{\sqrt{|\gamma|}} \\
B^0_{(u)} &= \frac{u_{i+1} B^i}{\alpha} \\ 
b^0 &= \frac{B^0_{(u)}}{\sqrt{4 \pi}} \\
b^{i+1} &= \frac{\frac{B^i}{\alpha} + B^0_{(u)} u^{i+1}}{u^0 \sqrt{4 \pi}}
\end{align}


```python
###############
# Next compute b^{\mu} using Eqs 23 and 31 of https://arxiv.org/pdf/astro-ph/0503420.pdf
uKS4D = ixp.zerorank1(DIM=4)
for i in range(4):
    for j in range(4):
        uKS4D[i] += gPhys4DD[i][j] * uKS4U[j]

# Eq 27 of https://arxiv.org/pdf/astro-ph/0503420.pdf
BU = ixp.zerorank1()
for i in range(3):
    BU[i] = BtildeU[i]/sp.sqrt(gammaDET)

# Eq 23 of https://arxiv.org/pdf/astro-ph/0503420.pdf
BU0_u = sp.sympify(0)
for i in range(3):
    BU0_u += uKS4D[i+1]*BU[i]/alpha

smallbU = ixp.zerorank1(DIM=4)
smallbU[0] = BU0_u   / sp.sqrt(4 * sp.pi)
# Eqs 24 and 31 of https://arxiv.org/pdf/astro-ph/0503420.pdf
for i in range(3):
    smallbU[i+1] = (BU[i]/alpha + BU0_u*uKS4U[i+1])/(sp.sqrt(4*sp.pi)*uKS4U[0])

smallbD = ixp.zerorank1(DIM=4)
for i in range(4):
    for j in range(4):
        smallbD[i] += gPhys4DD[i][j]*smallbU[j]

smallb2 = sp.sympify(0)
for i in range(4):
    smallb2 += smallbU[i]*smallbD[i]
```

<a id='lorentz_fac_valencia'></a>

## Step 2.l: Lorentz factor $\Gamma = \alpha u^0$ and Valencia 3-velocity $v^i_{(n)}$ \[Back to [top](#toc)\]
$$\label{lorentz_fac_valencia}$$

Now, we will define the Lorentz factor ($= \alpha u^0$) and the Valencia 3-velocity $v^i_{(n)}$, which sets the 3-velocity as measured by normal observers to the spatial slice: 
\begin{align}
v^i_{(n)} &= \frac{u^i}{u^0 \alpha} + \frac{\beta^i}{\alpha}, \\
\end{align}
as shown in eq 11 of [this](https://arxiv.org/pdf/1501.07276.pdf) paper. We will also compute the product of the square root of the determinant of the 3-metric with the lapse.


```python
###############
LorentzFactor = alpha * uKS4U[0]
# Define Valencia 3-velocity v^i_(n), which sets the 3-velocity as measured by normal observers to the spatial slice:
#  v^i_(n) = u^i/(u^0*alpha) + beta^i/alpha. See eq 11 of https://arxiv.org/pdf/1501.07276.pdf
Valencia3velocityU = ixp.zerorank1()
for i in range(3):
    Valencia3velocityU[i] = uKS4U[i + 1] / (alpha * uKS4U[0]) + betaU[i] / alpha

# sqrtgamma4DET = sp.sqrt(gammaDET)*alpha
```

<a id='output_to_c'></a>

## Step 3: Output above-generated expressions to C code, using NRPy+ \[Back to [top](#toc)\]
$$\label{output_to_c}$$

Finally, we have constructed the underlying expressions necessary for the Fishbone-Moncrief initial data. By means of demonstration, we will use NRPy+'s `FD_outputC()` to print the expressions. (The actual output statements are commented out right now, to save time in testing.)


```python
from outputC import lhrh   # NRPy+: Core C code output module
KerrSchild_CKernel = [\
                     lhrh(lhs=gri.gfaccess("out_gfs","alpha"),rhs=alpha),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU0"),rhs=betaU[0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU1"),rhs=betaU[1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU2"),rhs=betaU[2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD00"),rhs=gammaDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD01"),rhs=gammaDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD02"),rhs=gammaDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD11"),rhs=gammaDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD12"),rhs=gammaDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD22"),rhs=gammaDD[2][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD00"),rhs=KDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD01"),rhs=KDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD02"),rhs=KDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD11"),rhs=KDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD12"),rhs=KDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD22"),rhs=KDD[2][2]),\
                     ]

#fin.FD_outputC("stdout",KerrSchild_CKernel)

FMdisk_Lorentz_uUs_CKernel = [\
                             lhrh(lhs=gri.gfaccess("out_gfs","LorentzFactor"),rhs=LorentzFactor),\
#                              lhrh(lhs=gri.gfaccess("out_gfs","uKS4U1"),rhs=uKS4U[1]),\
#                              lhrh(lhs=gri.gfaccess("out_gfs","uKS4U2"),rhs=uKS4U[2]),\
#                              lhrh(lhs=gri.gfaccess("out_gfs","uKS4U3"),rhs=uKS4U[3]),\
                             ]

#fin.FD_outputC("stdout",FMdisk_Lorentz_uUs_CKernel)

FMdisk_hm1_rho_P_CKernel = [\
#                           lhrh(lhs=gri.gfaccess("out_gfs","hm1"),rhs=hm1),\
                           lhrh(lhs=gri.gfaccess("out_gfs","rho_initial"),rhs=rho_initial),\
                           lhrh(lhs=gri.gfaccess("out_gfs","Pressure_initial"),rhs=Pressure_initial),\
                           ]

#fin.FD_outputC("stdout",FMdisk_hm1_rho_P_CKernel)

udotu = sp.sympify(0)
for i in range(4):
    udotu += uKS4U[i]*uKS4D[i]
#NRPy_file_output(OUTDIR+"/standalone-spherical_coords/NRPy_codegen/FMdisk_Btildes.h", [],[],[],
#                 ID_protected_variables + ["r","th","ph"],
#                 [],[uKS4U[0], "uKS4Ut", uKS4U[1],"uKS4Ur", uKS4U[2],"uKS4Uth", uKS4U[3],"uKS4Uph",
#                     uKS4D[0], "uKS4Dt", uKS4D[1],"uKS4Dr", uKS4D[2],"uKS4Dth", uKS4D[3],"uKS4Dph",
#                     uKS4D[1] * BU[0] / alpha, "Bur", uKS4D[2] * BU[1] / alpha, "Buth", uKS4D[3] * BU[2] / alpha, "Buph",
#                     gPhys4DD[0][0], "g4DD00", gPhys4DD[0][1], "g4DD01",gPhys4DD[0][2], "g4DD02",gPhys4DD[0][3], "g4DD03",
#                     BtildeU[0], "BtildeUr", BtildeU[1], "BtildeUth",BtildeU[2], "BtildeUph",
#                     smallbU[0], "smallbUt", smallbU[1], "smallbUr", smallbU[2], "smallbUth",smallbU[3], "smallbUph",
#                     smallb2,"smallb2",udotu,"udotu"])

FMdisk_Btildes_CKernel = [\
                         lhrh(lhs=gri.gfaccess("out_gfs","BtildeU0"),rhs=BtildeU[0]),\
                         lhrh(lhs=gri.gfaccess("out_gfs","BtildeU1"),rhs=BtildeU[1]),\
                         lhrh(lhs=gri.gfaccess("out_gfs","BtildeU2"),rhs=BtildeU[2]),\
                         ]

#fin.FD_outputC("stdout",FMdisk_Btildes_CKernel)
```

We will now use the relationships between coordinate systems provided by [reference_metric.py](../edit/reference_metric.py) to convert our expressions to Cartesian coordinates. See [Tutorial-Reference_Metric](Tutorial-Reference_Metric.ipynb) for more detail.


```python
# Now that all derivatives of ghat and gbar have been computed,
# we may now substitute the definitions r = rfm.xxSph[0], th=rfm.xxSph[1],...
# WARNING: Substitution only works when the variable is not an integer. Hence the if not isinstance(...,...) stuff.
# If the variable isn't an integer, we revert transcendental functions inside to normal variables. E.g., sin(x2) -> sinx2
#  Reverting to normal variables in this way makes expressions simpler in NRPy, and enables transcendental functions
#  to be pre-computed in SENR.

alpha = alpha.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
for i in range(DIM):
    betaU[i] = betaU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    for j in range(DIM):
        gammaDD[i][j] = gammaDD[i][j].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
        KDD[i][j]     = KDD[i][j].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

# GRMHD variables:
# Density and pressure:
hm1           = hm1.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
rho_initial          = rho_initial.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
Pressure_initial     = Pressure_initial.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
LorentzFactor = LorentzFactor.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

# "Valencia" three-velocity
for i in range(DIM):
    BtildeU[i] = BtildeU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    uKS4U[i+1] = uKS4U[i+1].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    uBL4U[i+1] = uBL4U[i+1].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    Valencia3velocityU[i] = Valencia3velocityU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

```

At last, we will use our reference metric formalism and the Jacobian associated with the two coordinate systems to convert the spherical initial data to Cartesian coordinates. The module reference_metric.py provides us with the definition of $r, \theta, \phi$ in Cartesian coordinates. To find the Jacobian to then transform from spherical to Cartesian, we must find the tensor \begin{equation} \frac{\partial x_i}{\partial y_j}, \end{equation} where $x_i \in \{r,\theta,\phi\}$ and $y_i \in \{x,y,z\}$. We will also compute its inverse.


```python

# uUphi = uKS4U[3]
# uUphi = sympify_integers__replace_rthph(uUphi,r,th,ph,rfm.xxSph[0],rfm.xxSph[1],rfm.xxSph[2])
# uUt = uKS4U[0]
# uUt = sympify_integers__replace_rthph(uUt,r,th,ph,rfm.xxSph[0],rfm.xxSph[1],rfm.xxSph[2])

# Transform initial data to our coordinate system:
# First compute Jacobian and its inverse
drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff( rfm.xxSph[0],rfm.xx[1]), sp.diff( rfm.xxSph[0],rfm.xx[2])],
                                       [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                       [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv()

# Declare as gridfunctions the final quantities we will output for the initial data
IDalpha = gri.register_gridfunctions("EVOL","IDalpha")
IDgammaDD = ixp.register_gridfunctions_for_single_rank2("EVOL","IDgammaDD","sym01")
IDKDD = ixp.register_gridfunctions_for_single_rank2("EVOL","IDKDD","sym01")
IDbetaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","IDbetaU")
IDValencia3velocityU = ixp.register_gridfunctions_for_single_rank1("EVOL","IDValencia3velocityU")

IDalpha = alpha
for i in range(3):
    IDbetaU[i] = 0
    IDValencia3velocityU[i] = 0
    for j in range(3):
        # Matrices are stored in row, column format, so (i,j) <-> (row,column)
        IDbetaU[i]   += dx__drrefmetric_0UDmatrix[(i,j)]*betaU[j]
        IDValencia3velocityU[i]   += dx__drrefmetric_0UDmatrix[(i,j)]*Valencia3velocityU[j]
        IDgammaDD[i][j] = 0
        IDKDD[i][j] = 0
        for k in range(3):
            for l in range(3):
                IDgammaDD[i][j] += drrefmetric__dx_0UDmatrix[(k,i)]*drrefmetric__dx_0UDmatrix[(l,j)]*gammaDD[k][l]
                IDKDD[i][j]     += drrefmetric__dx_0UDmatrix[(k,i)]*drrefmetric__dx_0UDmatrix[(l,j)]*    KDD[k][l]


# -={ Spacetime quantities: Generate C code from expressions and output to file }=-
KerrSchild_to_print = [\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDalpha"),rhs=IDalpha),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDbetaU0"),rhs=IDbetaU[0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDbetaU1"),rhs=IDbetaU[1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDbetaU2"),rhs=IDbetaU[2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD00"),rhs=IDgammaDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD01"),rhs=IDgammaDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD02"),rhs=IDgammaDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD11"),rhs=IDgammaDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD12"),rhs=IDgammaDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDgammaDD22"),rhs=IDgammaDD[2][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD00"),rhs=IDKDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD01"),rhs=IDKDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD02"),rhs=IDKDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD11"),rhs=IDKDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD12"),rhs=IDKDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","IDKDD22"),rhs=IDKDD[2][2]),\
                     ]

# -={ GRMHD quantities: Generate C code from expressions and output to file }=-
FMdisk_GRHD_hm1_to_print = [lhrh(lhs=gri.gfaccess("out_gfs","rho_initial"),rhs=rho_initial)]

FMdisk_GRHD_velocities_to_print = [\
                                 lhrh(lhs=gri.gfaccess("out_gfs","IDValencia3velocityU0"),rhs=IDValencia3velocityU[0]),\
                                 lhrh(lhs=gri.gfaccess("out_gfs","IDValencia3velocityU1"),rhs=IDValencia3velocityU[1]),\
                                 lhrh(lhs=gri.gfaccess("out_gfs","IDValencia3velocityU2"),rhs=IDValencia3velocityU[2]),\
                                 ]
```

To verify this against the old version of FishboneMoncriefID from the old version of NRPy, we use the `mathematica_code()` output function. 


```python
# Comment out debug code for now, to reduce this file's size.

#from mathematica_output import *

# print("ID1alpha = " + sp.mathematica_code(IDalpha) + ";")
# print("ID1beta0 = " + sp.mathematica_code(IDbetaU[0]) + ";")
# print("ID1beta1 = " + sp.mathematica_code(IDbetaU[1]) + ";")
# print("ID1beta2 = " + sp.mathematica_code(IDbetaU[2]) + ";")
# print("ID1gamma00 = " + sp.mathematica_code(IDgammaDD[0][0]) + ";")
# print("ID1gamma01 = " + sp.mathematica_code(IDgammaDD[0][1]) + ";")
# print("ID1gamma02 = " + sp.mathematica_code(IDgammaDD[0][2]) + ";")
# print("ID1gamma11 = " + sp.mathematica_code(IDgammaDD[1][1]) + ";")
# print("ID1gamma12 = " + sp.mathematica_code(IDgammaDD[1][2]) + ";")
# print("ID1gamma22 = " + sp.mathematica_code(IDgammaDD[2][2]) + ";")
# print("ID1K00 = " + sp.mathematica_code(IDKDD[0][0]) + ";")
# print("ID1K01 = " + sp.mathematica_code(IDKDD[0][1]) + ";")
# print("ID1K02 = " + sp.mathematica_code(IDKDD[0][2]) + ";")
# print("ID1K11 = " + sp.mathematica_code(IDKDD[1][1]) + ";")
# print("ID1K12 = " + sp.mathematica_code(IDKDD[1][2]) + ";")
# print("ID1K22 = " + sp.mathematica_code(IDKDD[2][2]) + ";")

# print("hm11 = " + sp.mathematica_code(hm1) + ";")

# print("ID1Valencia3velocityU0 = " + sp.mathematica_code(IDValencia3velocityU[0]) + ";")
# print("ID1Valencia3velocityU1 = " + sp.mathematica_code(IDValencia3velocityU[1]) + ";")
# print("ID1Valencia3velocityU2 = " + sp.mathematica_code(IDValencia3velocityU[2]) + ";")

```

<a id='code_validation'></a>

# Step 4: Code Validation against `FishboneMoncriefID.FishboneMoncriefID` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for these Fishbone-Moncrief initial data between

1. this tutorial and 
2. the NRPy+ [FishboneMoncriefID.FishboneMoncriefID](../edit/FishboneMoncriefID/FishboneMoncriefID.py) module.


```python
gri.glb_gridfcs_list = []

import FishboneMoncriefID.FishboneMoncriefID as fmid
fmid.FishboneMoncriefID()

print("IDalpha - fmid.IDalpha = " + str(IDalpha - fmid.IDalpha))
print("rho_initial - fmid.rho_initial = " + str(rho_initial - fmid.rho_initial))
print("hm1 - fmid.hm1 = " + str(hm1 - fmid.hm1))

for i in range(DIM):
    print("IDbetaU["+str(i)+"] - fmid.IDbetaU["+str(i)+"] = " + str(IDbetaU[i] - fmid.IDbetaU[i]))
    print("IDValencia3velocityU["+str(i)+"] - fmid.IDValencia3velocityU["+str(i)+"] = "\
          + str(IDValencia3velocityU[i] - fmid.IDValencia3velocityU[i]))
    for j in range(DIM):
        print("IDgammaDD["+str(i)+"]["+str(j)+"] - fmid.IDgammaDD["+str(i)+"]["+str(j)+"] = "
              + str(IDgammaDD[i][j] - fmid.IDgammaDD[i][j]))
        print("IDKDD["+str(i)+"]["+str(j)+"] - fmid.IDKDD["+str(i)+"]["+str(j)+"] = "
              + str(IDKDD[i][j] - fmid.IDKDD[i][j]))
```

    IDalpha - fmid.IDalpha = 0
    rho_initial - fmid.rho_initial = 0
    hm1 - fmid.hm1 = 0
    IDbetaU[0] - fmid.IDbetaU[0] = 0
    IDValencia3velocityU[0] - fmid.IDValencia3velocityU[0] = 0
    IDgammaDD[0][0] - fmid.IDgammaDD[0][0] = 0
    IDKDD[0][0] - fmid.IDKDD[0][0] = 0
    IDgammaDD[0][1] - fmid.IDgammaDD[0][1] = 0
    IDKDD[0][1] - fmid.IDKDD[0][1] = 0
    IDgammaDD[0][2] - fmid.IDgammaDD[0][2] = 0
    IDKDD[0][2] - fmid.IDKDD[0][2] = 0
    IDbetaU[1] - fmid.IDbetaU[1] = 0
    IDValencia3velocityU[1] - fmid.IDValencia3velocityU[1] = 0
    IDgammaDD[1][0] - fmid.IDgammaDD[1][0] = 0
    IDKDD[1][0] - fmid.IDKDD[1][0] = 0
    IDgammaDD[1][1] - fmid.IDgammaDD[1][1] = 0
    IDKDD[1][1] - fmid.IDKDD[1][1] = 0
    IDgammaDD[1][2] - fmid.IDgammaDD[1][2] = 0
    IDKDD[1][2] - fmid.IDKDD[1][2] = 0
    IDbetaU[2] - fmid.IDbetaU[2] = 0
    IDValencia3velocityU[2] - fmid.IDValencia3velocityU[2] = 0
    IDgammaDD[2][0] - fmid.IDgammaDD[2][0] = 0
    IDKDD[2][0] - fmid.IDKDD[2][0] = 0
    IDgammaDD[2][1] - fmid.IDgammaDD[2][1] = 0
    IDKDD[2][1] - fmid.IDKDD[2][1] = 0
    IDgammaDD[2][2] - fmid.IDgammaDD[2][2] = 0
    IDKDD[2][2] - fmid.IDKDD[2][2] = 0


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename [Tutorial-FishboneMoncriefID.pdf](Tutorial-FishboneMoncriefID.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-FishboneMoncriefID")
```

    Created Tutorial-FishboneMoncriefID.tex, and compiled LaTeX file to PDF
        file Tutorial-FishboneMoncriefID.pdf

