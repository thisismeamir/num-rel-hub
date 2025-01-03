<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# The Spinning Effective One-Body Hamiltonian: "v4P"

## Author: Tyler Knowles

## This module documents the reduced spinning effective one-body Hamiltonian as numerically implemented in LALSuite's SEOBNRv4P gravitational waveform approximant. Through a canonical transformation, $H_{\text{real}}$ translates to $H_{\text{eff}}$ representing a test particle's trajectory in a deformed Kerr background.  The notebook presents Python code for $H_{\text{real}}$ computation, utilizing inputs such as black hole masses, the Euler-Mascheroni constant, a tortoise coordinate, and twelve dynamic variables.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated against the LALSuite [SEOBNRv4P code]( https://git.ligo.org/lscsoft/lalsuite.) that was reviewed and approved for LIGO parameter estimation by the LIGO Scientific Collaboration.  That is, the value $H_{\rm real}$ output from this notebook agrees to roundoff error with the value of $H_{\rm real}$ computed by the LALSuite function XLALSimIMRSpinPrecEOBHamiltonian().

### NRPy+ Source Code for this module: [SEOBNR_v4P_Hamiltonian.py](../edit/SEOBNR/SEOBNR_v4P_Hamiltonian.py)

<a id='intro'></a>

## Introduction
$$\label{intro}$$

### The Physical System of Interest

Consider two black holes with masses $m_{1}$, $m_{2}$ and spins ${\bf S}_{1}$, ${\bf S}_{2}$ in a binary system.  The spinning effective one-body ("SEOB") Hamiltonian $H_{\rm real}$ (defined in [this cell](#hreal)) describes the dynamics of this system; we will define $H_{\rm real}$ as in [Barausse and Buonanno (2010)](https://arxiv.org/abs/0912.3517) Section VE.  There, $H_{\rm real}$ is canonically transformed and mapped to an effective Hamiltonian $H_{\rm eff}$ (defined in [this cell](#heff)) describing the motion of a test particle of mass $\mu$ (defined in [this cell](#mu)) and spin ${\bf S}^{*}$ (defined in [this cell](#sstar)) moving in a deformed Kerr background.  Here we seek to break up $H_{\rm real}$ and document the terms in such a way that the resulting Python code can be used to numerically evaluate $H_{\rm real}$.

We write $H_{\rm real}$ in terms of Cartesian quasi-isotropic coordinates $x$, $y$, and $z$ (see [Barausse and Buonanno (2010)](https://arxiv.org/abs/0912.3517) Section III).  The spatial coordinates $r$, $\theta$, and $\phi$ referenced throughout are [Boyer-Lindquist coordinates](https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates) (see [Barausse and Buonanno (2010)](https://arxiv.org/abs/0912.3517) Section IV).

Please note that throughout this notebook we adopt the following conventions:

1. $c = 1$ where $c$ is the speed of light in a vacuum,
1. spacial tensor indices are denoted by lowercase Latin letters,
1. repeated indices indicate Einstein summation notation, and
1. we normalize $M=1$ in all terms except for $\eta$ and $\mu$ for agreement with LALSuite.  Nonetheless, $M$ appears in other text cells for comparison with the cited literature.

Running this notebook to completion will generate a file called v4P_Hreal_on_bottom.py.  This file contains the Python function v4P_compute_Hreal(), which takes as input m1, m2 (each in solar masses), the value of the [Euler-Mascheroni constant](https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant), a tortoise coordinate, and values for all twelve dynamic variables (3 components of the separation vector, three components of the momentum vector, and three spin components of each compact object).  Note that the spin components should be dimensionless.

### Citations
Throughout this module, we will refer to
* [Barausse and Buonanno (2010)](https://arxiv.org/abs/0912.3517) as BB2010,
* [Barausse and Buonanno (2011)](https://arxiv.org/abs/1107.2904) as BB2011,
* [Ossokine, Buonanno, Marsat, et al (2020)](https://arxiv.org/abs/2004.09442) as OB2020,
* [Steinhoff, Hinderer, Buonanno, et al (2016)](https://arxiv.org/abs/1608.01907) as SH2016,
* [Bohe, Shao, Taracchini, et al (2017)](https://arxiv.org/pdf/1611.03703.pdf) as BL2017,
* [Pan, Buonanno, Buchman, et. al (2010)](https://arxiv.org/abs/0912.3466v2) as P2010,
* [Taracchini, Buonanno, Pan, et al (2014)](https://arxiv.org/abs/1311.2544) as T2014,
* [Taracchini, Pan, Buonanno, et al (2012)](https://arxiv.org/abs/1202.0790) as T2012, and
* [Damour, Jaranowski, and Schaefer (2000)](https://arxiv.org/abs/gr-qc/0005034) as D2000.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 0](#outputcreation): Creating the output directory for SEOBNR
1. [Step 1](#hreal): The Real Hamiltonian $H_{\rm real}$
1. [Step 2](#heff): The Effective Hamiltonian $H_{\rm eff}$
1. [Step 3](#heff_terms): Terms of $H_{\rm eff}$  
    1. [Step 3.a](#hs): Leading Order Spin Effects $H_{\rm S}$  
    1. [Step 3.b](#hns): The Nonspinning Hamiltonian $H_{\rm NS}$  
    1. [Step 3.c](#hd): The Quadrupole Deformation $H_{\rm D}$
1. [Step 4](#hso): The Spin-Orbit Term $H_{\rm SO}$  
    1. [Step 4.a](#hsoterm1): $H_{\rm SO}$ Term 1  
    1. [Step 4.b](#hsoterm2coeff): $H_{\rm SO}$ Term 2 Coefficient  
    1. [Step 4.c](#hsoterm2): $H_{\rm SO}$ Term 2  
        1. [Step 4.c.i](#hsoterm2a): $H_{\rm SO}$ Term 2a  
        1. [Step 4.c.ii](#hsoterm2b): $H_{\rm SO}$ Term 2b  
        1. [Step 4.c.iii](#hsoterm2c): $H_{\rm SO}$ Term 2c
1. [Step 5](#hss): The Spin-Spin Term $H_{\rm SS}$  
    1. [Step 5.a](#hssterm1): $H_{\rm SS}$ Term 1  
    1. [Step 5.b](#hssterm2coeff): $H_{\rm SS}$ Term 2 coefficient  
    1. [Step 5.c](#hssterm2): $H_{\rm SS}$ Term 2  
    1. [Step 5.d](#hssterm3coeff): $H_{\rm SS}$ Term 3 coefficient  
    1. [Step 5.e](#hssterm3): $H_{\rm SS}$ Term 3
1. [Step 6](#hnsterms): The $H_{\rm NS}$ Terms  
    1. [Step 6.a](#betapsum): $\beta p$ Sum  
    1. [Step 6.b](#alpha): $\alpha$  
    1. [Step 6.c](#hnsradicand): $H_{\rm NS}$ Radicand  
        1. [Step 6.c.i](#gammappsum): $\gamma p$ Sum  
        1. [Step 6.c.ii](#q4): ${\cal Q}_{4}$
1. [Step 7](#hdterms): The $H_{\rm D}$ Terms  
    1. [Step 7.a](#hdcoeff): $H_{\rm D}$ Coefficient  
    1. [Step 7.b](#hdsum): $H_{\rm D}$ Sum  
        1. [Step 7.b.i](#hdsumterm1): $H_{\rm D}$ Sum Term 1  
        1. [Step 7.b.ii](#hdsumterm2): $H_{\rm D}$ Sum Term 2
1. [Step 8](#dotproducts): Common Dot Products
    1. [Step 8.a](#sdotxi): ${\bf S} \cdot \boldsymbol{\xi}$  
    1. [Step 8.b](#sdotv): ${\bf S} \cdot {\bf v}$  
    1. [Step 8.c](#sdotn): ${\bf S} \cdot {\bf n}$  
    1. [Step 8.d](#sdotskerrhat): ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$  
    1. [Step 8.e](#sstardotn): ${\bf S}^{*} \cdot {\bf n}$
1. [Step 9](#hreal_spin_combos): $H_{\rm real}$ Spin Combination ${\bf S}^{*}$  
    1. [Step 9a](#sstar): ${\bf S}^{*}$  
    1. [Step 9b](#deltasigmastar): $\Delta_{\sigma^{*}}$  
    1. [Step 9c](#sigmastarcoeff): $\sigma^{*}$ Coefficient  
        1. [Step 9c i](#sigmastarcoeffterm1): $\sigma^{*}$ Coefficient Term 1  
        1. [Step 9c ii](#sigmastarcoeffterm2): $\sigma^{*}$ Coefficient Term 2   
    1. [Step 9d](#sigmacoeff): $\sigma$ Coefficient  
        1. [Step 9d i](#sigmacoeffterm1): $\sigma$ Coefficient Term 1  
        1. [Step 9d ii](#sigmacoeffterm2): $\sigma$ Coefficient Term 2  
        1. [Step 9d iii](#sigmacoeffterm3): $\sigma$ Coefficient Term 3
1. [Step 10](#metpotderivs): Derivatives of the Metric Potential  
    1. [Step 10.a](#omegar): $\omega_{r}$  
    1. [Step 10.b](#nur): $\nu_{r}$  
    1. [Step 10.c](#mur): $\mu_{r}$  
    1. [Step 10.d](#omegacostheta): $\omega_{\cos\theta}$  
    1. [Step 10.e](#nucostheta): $\nu_{\cos\theta}$  
    1. [Step 10.f](#mucostheta): $\mu_{\cos\theta}$  
    1. [Step 10.g](#lambdatprm): $\Lambda_{t}^{\prime}$  
    1. [Step 10.h](#omegatildeprm): $\tilde{\omega}_{\rm fd}^{\prime}$
1. [Step 11](#metpots): The Deformed and Rescaled Metric Potentials  
    1. [Step 11.a](#omega): $\omega$  
    1. [Step 11.b](#exp2nu): $e^{2 \nu}$  
    1. [Step 11.c](#btilde): $\tilde{B}$  
    1. [Step 11.d](#brtilde): $\tilde{B}_{r}$  
    1. [Step 11.e](#exp2mu): $e^{2 \tilde{\mu}}$  
    1. [Step 11.f](#jtilde): $\tilde{J}$  
    1. [Step 11.g](#q): $Q$  
        1. [Step 11.g.i](#drsipn2): $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$  
        1. [Step 11.g.ii](#qcoeff1): Q Coefficient 1  
        1. [Step 11.g.iii](#qcoeff2): Q Coefficient 2
1. [Step 12](#tort): Tortoise terms  
    1. [Step 12.a](#pphi): $p_{\phi}$  
    1. [Step 12.b](#pdotvr): $\hat{\bf p} \cdot {\bf v} r$  
    1. [Step 12.c](#pdotn): $\hat{\bf p} \cdot {\bf n}$  
    1. [Step 12.d](#pdotxir): $\hat{\bf p} \cdot \boldsymbol{\xi} r$  
    1. [Step 12.e](#hatp): $\hat{\bf p}$  
    1. [Step 12.f](#prt): prT  
    1. [Step 12.g](#csi2): csi2  
    1. [Step 12.h](#csi1): csi1  
    1. [Step 12.i](#csi): csi
1. [Step 13](#metric): Metric Terms  
    1. [Step 13.a](#lambdat): $\Lambda_{t}$  
    1. [Step 13.b](#deltar): $\Delta_{r}$  
    1. [Step 13.c](#deltat): $\Delta_{t}$  
    1. [Step 13.d](#deltatprm): $\Delta_{t}^{\prime}$  
    1. [Step 13.e](#deltau): $\Delta_{u}$  
        1. [Step 13.e.i](#deltaubar): $\bar{\Delta}_{u}$  
        1. [Step 13.e.ii](#deltaucalib): $\Delta_{u}$ Calibration Term  
        1. [Step 13.e.iii](#calib_coeffs): Calibration Coefficients  
        1. [Step 13.e.iv](#k): $K$  
        1. [Step 13.e.v](#chi): $\chi$ 
    1. [Step 13.f](#omegatilde): $\tilde{\omega}_{\rm fd}$  
    1. [Step 13.g](#dinv): $D^{-1}$
1. [Step 14](#coord): Terms Dependent on Coordinates  
    1. [Step 14.a](#usigma): $\Sigma$  
    1. [Step 14.b](#w2): $\varpi^{2}$   
    1. [Step 14.d](#sin2theta): $\sin^{2}\theta$  
    1. [Step 14.e](#costheta): $\cos\theta$
1. [Step 15](#vectors): Important Vectors  
    1. [Step 15.a](#v): ${\bf v}$  
    1. [Step 15.b](#xi): $\boldsymbol{\xi}$  
    1. [Step 15.c](#e3): ${\bf e}_{3}$  
    1. [Step 15.d](#n): ${\bf n}$
    1. [Step 15.e](#sperp): ${\bf S}^{\perp}$
    1. [Step 15.f](#orb_momentum): ${\bf L}$
1. [Step 16](#spin_combos): Spin Combinations $\boldsymbol{\sigma}$, $\boldsymbol{\sigma}^{*}$, and ${\bf S}_{\rm Kerr}$   
    1. [Step 16.a](#a): $a$  
    1. [Step 16.b](#skerrhat): $\hat{\bf S}_{\rm Kerr}$  
    1. [Step 16.c](#skerrmag): $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$  
    1. [Step 16.d](#skerr): ${\bf S}_{\rm Kerr}$  
    1. [Step 16.e](#sigma): $\boldsymbol{\sigma}$  
    1. [Step 16.f](#sigmastar): $\boldsymbol{\sigma}^{*}$
1. [Step 17](#fundquant): Fundamental Quantities  
    1. [Step 17.a](#u): $u$  
    1. [Step 17.b](#r): $r$  
    1. [Step 17.c](#eta): $\eta$  
    1. [Step 17.d](#mu): $\mu$  
    1. [Step 17.e](#m): $M$
1. [Step 18](#validation): Validation
1. [Step 19](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='outputcreation'></a>

# Step 0: Creating the output directory for SEOBNR \[Back to [top](#toc)\]
$$\label{outputcreation}$$

First we create the output directory for SEOBNR (if it does not already exist):


```python
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface

# Create C code output directory:
Ccodesdir = "SEOBNR"
# Then create an output directory in case it does not exist
cmd.mkdir(Ccodesdir)
```

<a id='hreal'></a>

# Step 1: The real Hamiltonian $H_{\textrm{real}}$ \[Back to [top](#toc)\]
$$\label{hreal}$$

The SEOB Hamiltonian $H_{\rm real}$ is given by [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.69):

\begin{equation*}
    H_{\rm real} = M \sqrt{ 1 + 2 \eta \left( \frac{ H_{\rm eff} }{ \mu } - 1 \right) }.
\end{equation*}

Here $H_{\rm eff}$ (defined in [this cell](#heff)) is an *effective* Hamiltonian (see [this cell](#intro)) and $M$ (defined in [this cell](#m)), $\mu$ (defined in [this cell](#mu)), and $\eta$ (defined in [this cell](#eta)) are constants determined by $m_{1}$ and $m_{2}$.


```python
%%writefile $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt
Hreal = sp.sqrt(1 + 2*eta*(Heff - 1))
```

    Overwriting SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='heff'></a>

# Step 2: The Effective Hamiltonian $H_{\rm eff}$ \[Back to [top](#toc)\]
$$\label{heff}$$

The effective Hamiltonian $H_{\rm eff}$ is given by [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.70):

\begin{equation*}
    H_{\rm eff} = H_{\rm S} + \underbrace{ \beta^{i} p_{i} + \alpha \sqrt{ \mu^{2} + \gamma^{ij} p_{i} p_{j} + {\cal Q}_{4} } }_{ H_{\rm NS} } - \underbrace{ \frac{ \mu }{ 2 M r^{3} } \left( \delta^{ij} - 3 n^{i} n^{j} \right) S^{*}_{i} S^{*}_{j} }_{ H_{\rm D} }.
\end{equation*}

Here $H_{\rm S}$ (considered further in [this cell](#hs)) denotes leading order effects of spin-spin and spin-orbit coupling, $H_{\rm NS}$ (considered further in [this cell](#hns)) is the Hamiltonian for a nonspinning test particle, and $H_{\rm D}$ (considered further in [this cell](#hd)) describes quadrupole deformation of the coupling of the particle's spin with itself to leading order.  [T2014](https://arxiv.org/abs/1311.2544) adds to $H_{\rm eff}$ a 3PN spin-spin term given by

\begin{equation*}
    \frac{d_{\rm SS} \eta }{ r^{4} } \left( {\bf S}_{1}^{2} + {\bf S}_{2}^{2} \right)
\end{equation*}

where $d_{\rm SS}$ is an adjustable parameter determined by fitting to numerical relativity results.  Equation (4.13) of [BL2017](https://arxiv.org/pdf/1611.03703.pdf) gives

\begin{equation*}
    d_{\rm SS} = 528.511252 \chi^{3} \eta^{2} - 41.000256 \chi^{3} \eta + 1161.780126 \chi^{2} \eta^{3}  - 326.324859 \chi^{2} \eta^{2} + 37.196389 \chi \eta + 706.958312 \eta^{3} - 36.027203 \eta + 6.068071.
\end{equation*}

We take $u \equiv \frac{1}{r}$ (as described in [this cell](#u)), and define $\eta$ in [this cell](#eta) $\chi$ in [this cell](#chi).  Note that the coefficients for $d_{\rm SS}$ have been rounded to coincide with the LALSuite implementation (see the file LALSimIMRSpinEOBHamiltonian.h).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
dSS = 528.511*chi*chi*chi*eta*eta - 41.0003*chi*chi*chi*eta + 1161.78*chi*chi*eta*eta*eta - 326.325*chi*chi*eta*eta
    + 37.1964*chi*eta + 706.958*eta*eta*eta - 36.0272*eta + 6.06807
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='heff_terms'></a>

# Step 3: Terms of $H_{\rm eff}$ \[Back to [top](#toc)\]
$$\label{heff_terms}$$

In this step, we break down each of the terms $H_{\rm S}$ (defined in [this cell](#hs)), $H_{\rm NS}$ (defined in [this cell](#hns)), and $H_{\rm D}$ (defined in [this cell](#hd)) in $H_{\rm eff}$ (defined in [this cell](#heff)).

<a id='hs'></a>

## Step 3.a: Leading Order Spin Effects $H_{\rm S}$ \[Back to [top](#toc)\]
$$\label{hs}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.17),

\begin{equation*}
    H_{\rm S} = H_{\rm SO} + H_{\rm SS}
\end{equation*}

where $H_{\rm SO}$ (defined in [this cell](#hso)) includes spin-orbit terms and $H_{\rm SS}$ (defined in [this cell](#hss)) includes spin-spin terms.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hs = Hso + Hss
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hns'></a>

## Step 3.b: The Nonspinning Hamiltonian $H_{\rm NS}$ \[Back to [top](#toc)\]
$$\label{hns}$$

We defined $H_{\rm NS}$ in [this cell](#heff) as

\begin{equation*}
    H_{\rm NS} = \underbrace{ \beta^{i} p_{i} }_{ \beta\ p\ \rm sum } + \alpha \sqrt{ \smash[b]{ \underbrace{ \mu^{2} + \gamma^{ij} p_{i} p_{j} + {\cal Q}_{4} }_{ H_{\rm NS}\ \rm radicand } } }.
\end{equation*}

We compute $\beta\ p$ sum in [this cell](#betapsum), $\alpha$ in [this cell](#alpha), and $H_{\rm NS}$ radicand in [this cell](#hnsradicand).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hns = betapsum + alpha*sp.sqrt(Hnsradicand)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hd'></a>

## Step 3.c: The Quadrupole Deformation $H_{\rm D}$ \[Back to [top](#toc)\]
$$\label{hd}$$

We defined $H_{\rm D}$ in [this cell](#heff) as:

\begin{equation*}
    H_{\rm D} = \underbrace{ \frac{ \mu }{ 2 M r^{3} } }_{H_{\rm D}\ {\rm coefficient}} \underbrace{ \left( \delta^{ij} - 3 n^{i} n^{j} \right) S^{*}_{i} S^{*}_{j} }_{H_{\rm D}\ {\rm sum}}
\end{equation*}

We compute $H_{\rm D}$ coefficient in [this cell](#hdcoeff) and $H_{\rm D}$ sum in [this cell](#hdsum).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hd = Hdcoeff*Hdsum
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hso'></a>

# Step 4: The Spin-Orbit Term $H_{\rm SO}$ \[Back to [top](#toc)\]
$$\label{hso}$$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.18) as:

\begin{align*}
    H_{\rm SO} = H_{\rm SO}\ {\rm Term\ 1} + H_{\rm SO}\ {\rm Term\ 2\ coefficient} * H_{\rm SO}\ {\rm Term\ 2}.
\end{align*}

We define and consider $H_{\rm SO}$ Term 1 in [this cell](#hsoterm1), $H_{\rm SO}$ Term 2 coefficient in [this cell](#hsoterm2coeff), and $H_{\rm SO}$ Term 2 in [this cell](#hsoterm2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hso = HsoTerm1 + HsoTerm2coeff*HsoTerm2
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm1'></a>

## Step 4.a: $H_{\rm SO}$ Term 1 \[Back to [top](#toc)\]
$$\label{hsoterm1}$$

Combining our notation $H_{\rm SO}$ (defined in [this cell](#hso)) with [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.18), we have

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 1} = \frac{ e^{2 \nu - \tilde{\mu} } \left( e^{\tilde{\mu} + \nu} - \tilde{B} \right) \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( {\bf S} \cdot \hat{\bf S}_{\rm Kerr} \right) }{ \tilde{B}^{2} \sqrt{Q} \xi^{2} }.
\end{equation*}

We will write

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 1} = \frac{ e^{2 \nu} \left( e^{\tilde{\mu}} e^{\nu} - \tilde{B} \right) \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( {\bf S} \cdot \hat{\bf S}_{\rm Kerr} \right) }{ e^{ \tilde{\mu} } \tilde{B}^{2} \sqrt{Q} \xi^{2} }.
\end{equation*}

We define $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), $\tilde{B}$ in [this cell](#btilde), $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir), ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$ in [this cell](#sdotskerrhat), $Q$ in [this cell](#q), and $\boldsymbol{\xi}^{2}$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm1 = exp2nu*(expmu*expnu - Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*sp.sqrt(Q)*xisq)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm2coeff'></a>

## Step 4.b: $H_{\rm SO}$ Term 2 Coefficient \[Back to [top](#toc)\]
$$\label{hsoterm2coeff}$$

Combining our notation $H_{\rm SO}$ (defined in [this cell](#hso)) with [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.18), we have

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2\ coefficient} = \frac{ e^{\nu - 2 \tilde{\mu}} }{ \tilde{B}^{2} \left( \sqrt{Q} + 1 \right) \sqrt{Q} \xi^{2} }
\end{equation*}

which we write in the form

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2\ coefficient} = \frac{ e^{\nu} }{ e^{2 \tilde{\mu}} \tilde{B}^{2} \left( Q + \sqrt{Q} \right) \xi^{2} }.
\end{equation*}

We define and consider $e^{\nu}$ in [this cell](#exp2nu), $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $\tilde{B}$ in [this cell](#btilde), $Q$ in [this cell](#q), and $\xi^{2}$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm2coeff = expnu/(exp2mu*Btilde*Btilde*(Q + sp.sqrt(Q))*xisq)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm2'></a>

## Step 4.c: $H_{\rm SO}$ Term 2 \[Back to [top](#toc)\]
$$\label{hsoterm2}$$

Combining our notation $H_{\rm SO}$ (defined in [this cell](#hso)) with [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.18), we have

\begin{align*}
    H_{\rm SO}\ {\rm Term\ 2} &= \underbrace{ \left( {\bf S} \cdot \boldsymbol{\xi} \right) \tilde{J} \left[ \mu_r \left( \hat{\bf p} \cdot {\bf v} r \right) \left( \sqrt{Q} + 1 \right) - \mu_{\cos \theta} \left( \hat{\bf p} \cdot {\bf n} \right) \xi^{2} -\sqrt{Q} \left( \nu_r \left( \hat{\bf p} \cdot {\bf v} r \right) + \left( \mu_{\cos \theta} - \nu_{\cos \theta} \right) \left( \hat{\bf p} \cdot {\bf n} \right) \xi^{2} \right) \right] \tilde{B}^{2} }_{H_{\rm SO}\ {\rm Term\ 2a}} \\
        &\ \ \ \ \ + \underbrace{ e^{\tilde{\mu} + \nu} \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( 2 \sqrt{Q} + 1 \right) \left[ \tilde{J} \nu_r \left( {\bf S} \cdot {\bf v} \right) - \nu_{\cos \theta} \left( {\bf S} \cdot {\bf n} \right) \xi^{2} \right] \tilde{B} }_{H_{\rm SO}\ {\rm Term\ 2b}} - \underbrace{ \tilde{J} \tilde{B}_{r} e^{\tilde{\mu} + \nu} \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( \sqrt{Q} + 1 \right) \left( {\bf S} \cdot {\bf v} \right) }_{H_{\rm SO}\ {\rm Term\ 2c}}
\end{align*}

We compute $H_{\rm SO}$ Term 2a in [this cell](#hsoterm2a), $H_{\rm SO}$ Term 2b in [this cell](#hsoterm2b), and $H_{\rm SO}$ Term 2c in [this cell](#hsoterm2c).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm2 = HsoTerm2a + HsoTerm2b - HsoTerm2c
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm2a'></a>

### Step 4.c.i: $H_{\rm SO}$ Term 2a \[Back to [top](#toc)\]
$$\label{hsoterm2a}$$

We defined $H_{\rm S0}$ Term 2a  in [this cell](#hsoterm2) as

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2a} = \left( {\bf S} \cdot \boldsymbol{\xi} \right) \tilde{J} \left[ \mu_r \left( \hat{\bf p} \cdot {\bf v} r \right) \left( \sqrt{Q} + 1 \right) - \mu_{\cos \theta} \left( \hat{\bf p} \cdot {\bf n} \right) \xi^{2} -\sqrt{Q} \left( \nu_r \left( \hat{\bf p} \cdot {\bf v} r \right) + \left( \mu_{\cos \theta} - \nu_{\cos \theta} \right) \left( \hat{\bf p} \cdot {\bf n} \right) \xi^{2} \right) \right] \tilde{B}^{2}.
\end{equation*}

We define ${\bf S} \cdot \boldsymbol{\xi}$ in [this cell](#sdotxi), $\tilde{J}$ in [this cell](#jtilde), $\mu_{r}$ in [this cell](#mur), $\hat{\bf p} \cdot {\bf v} r$ in [this cell](#pdotvr), $Q$ in [this cell](#q), $\mu_{\cos \theta}$ in [this cell](#mucostheta), $\hat{\bf p} \cdot {\bf n}$ in [this cell](#pdotn), $\xi^{2}$ in [this cell](#sin2theta), $\nu_{r}$ in [this cell](#nur), $\nu_{\cos\theta}$ in [this cell](#nucostheta), and $\tilde{B}$ in [this cell](#btilde).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm2a = Sdotxi*Jtilde*(mur*pdotvr*(sp.sqrt(Q) + 1) - mucostheta*pdotn*xisq
                           - sp.sqrt(Q)*(nur*pdotvr + (mucostheta - nucostheta)*pdotn*xisq))*Btilde*Btilde
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm2b'></a>

### Step 4.c.ii: $H_{\rm SO}$ Term 2b \[Back to [top](#toc)\]
$$\label{hsoterm2b}$$

We defined $H_{\rm S0}$ Term 2b  in [this cell](#hsoterm2) as

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2b} = e^{\tilde{\mu} + \nu} \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( 2 \sqrt{Q} + 1 \right) \left[ \tilde{J} \nu_r \left( {\bf S} \cdot {\bf v} \right) - \nu_{\cos \theta} \left( {\bf S} \cdot {\bf n} \right) \xi^{2} \right] \tilde{B}.
\end{equation*}

We define $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), $\hat{\bf p} \cdot \xi r$ in [this cell](#pdotxir), $Q$ in [this cell](#q), $\tilde{J}$ in [this cell](#jtilde), $\nu_{r}$ in [this cell](#nur), ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv), $\nu_{\cos\theta}$ in [this cell](#nucostheta), ${\bf S} \cdot {\bf n}$ in [this cell](#sdotn), $\xi^{2}$ in [this cell](#sin2theta), and $\tilde{B}$ in [this cell](#btilde).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm2b = expmu*expnu*pdotxir*(2*sp.sqrt(Q) + 1)*(Jtilde*nur*Sdotv - nucostheta*Sdotn*xisq)*Btilde
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hsoterm2c'></a>

### Step 4.c.iii: $H_{\rm SO}$ Term 2c \[Back to [top](#toc)\]
$$\label{hsoterm2c}$$

We defined $H_{\rm S0}$ Term 2c in [this cell](#hsoterm2) as

\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2c} = \tilde{J} \tilde{B}_{r} e^{\tilde{\mu} + \nu} \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( \sqrt{Q} + 1 \right) \left( {\bf S} \cdot {\bf v} \right)
\end{equation*}

We define $\tilde{J}$ in [this cell](#jtilde), $\tilde{B}_{r}$ in [this cell](#brtilde), $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), $\hat{\bf p} \cdot \xi r$ in [this cell](#pdotxir), $Q$ in [this cell](#q), and ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HsoTerm2c = Jtilde*Brtilde*expmu*expnu*pdotxir*(sp.sqrt(Q) + 1)*Sdotv
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hss'></a>

# Step 5: The Spin-Spin Term $H_{\rm SS}$ \[Back to [top](#toc)\]
$$\label{hss}$$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

\begin{equation*}
    H_{\rm SS} = H_{\rm SS}\ {\rm Term\ 1} + H_{\rm SS}\ {\rm Term\ 2\ coefficient} * H_{\rm SS}\ {\rm Term\ 2} + H_{\rm SS}\ {\rm Term\ 3\ coefficient} * H_{\rm SS}\ {\rm Term\ 3}.
\end{equation*}

We define $H_{\rm SS}$ Term 1 in [this cell](#hssterm1), $H_{\rm SS}$ Term 2 coefficient in [this cell](#hssterm2coeff), $H_{\rm SS}$ Term 2 in [this cell](#hssterm2), $H_{\rm SS}$ Term 3 coefficient in [this cell](#hssterm3coeff), and $H_{\rm SS}$ Term 3 in [this cell](#hssterm3).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hss = HssTerm1 + HssTerm2coeff*HssTerm2 + HssTerm3coeff*HssTerm3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hssterm1'></a>

## Step 5.a: $H_{\rm SS}$ Term 1 \[Back to [top](#toc)\]
$$\label{hssterm1}$$

Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with our definition of $H_{\rm SS}$ Term 1 in [this cell](#hss), we have

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 1} = \omega \left( {\bf S} \cdot \hat{\bf S}_{\rm Kerr} \right).
\end{equation*}

We define $\omega$ in [this cell](#omega) and ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$ in [this cell](#sdotskerrhat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HssTerm1 = omega*SdotSkerrhat
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hssterm2coeff'></a>

## Step 5.b: $H_{\rm SS}$ Term 2 Coefficient \[Back to [top](#toc)\]
$$\label{hssterm2coeff}$$

Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with ore definition of $H_{\rm SS}$ Term 2 coefficient in [this cell](#hss), we have

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 2\ coefficient} = \frac{ e^{-3 \tilde{\mu} -\nu} \tilde{J} \omega_{r} }{ 2 \tilde{B} \left( \sqrt{Q} + 1 \right) \sqrt{Q} \xi^{2} }
\end{equation*}

which we write as

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 2\ coefficient} = \frac{ \tilde{J} \omega_{r} }{ 2 e^{2 \tilde{\mu}} e^{\tilde{\mu}} e^{\nu} \tilde{B} \left( Q + \sqrt{Q} \right) \xi^{2} }.
\end{equation*}

We define $\tilde{J}$ in [this cell](#jtilde), $\omega_{r}$ in [this cell](#omegar), $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), $\tilde{B}$ in [this cell](#btilde), $Q$ in [this cell](#q),  and $\xi^{2}$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HssTerm2coeff = Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q + sp.sqrt(Q))*xisq)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hssterm2'></a>

## Step 5.c: $H_{\rm SS}$ Term 2 \[Back to [top](#toc)\]
$$\label{hssterm2}$$

Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with our definition of $H_{\rm SS}$ Term 2 in [this cell](#hss), we have

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 2} = -e^{\tilde{\mu} + \nu} \left( {\bf \hat{p}} \cdot {\bf v} r \right) \left( {\bf \hat{p}} \cdot {\bf \xi} r \right) \left( {\bf S} \cdot {\bf \xi} \right)
\tilde{B} + e^{2 \left( \tilde{\mu} + \nu \right)} \left( {\bf \hat{p}} \cdot {\bf \xi} r \right)^2 \left( {\bf S}
\cdot {\bf v} \right) + e^{2 \tilde{\mu}} \left( 1 + \sqrt{Q} \right) \sqrt{Q} \left( {\bf S} \cdot {\bf v} \right)\xi^2 \tilde{B}^{2} + \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left[ \left( {\bf \hat{p}} \cdot {\bf v} r \right)
\left( {\bf S} \cdot {\bf n}\right) - \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf S} \cdot {\bf v} \right)\right] \xi^{2} \tilde{B}^{2}
\end{equation*}

which we write as

\begin{align*}
    H_{\rm SS}\ {\rm Term\ 2} &= e^{\tilde{\mu}} \left( {\bf \hat{p}} \cdot {\bf \xi} r \right) \left[ e^{\tilde{\mu}} e^{2 \nu} \left( {\bf \hat{p}} \cdot {\bf \xi} r \right) \left( {\bf S} \cdot {\bf v} \right) - e^{\nu} \left( {\bf \hat{p}} \cdot {\bf v} r \right) \left( {\bf S} \cdot {\bf \xi} \right)
\tilde{B} \right] \\
&\ \ \ \ \ + \xi^2 \tilde{B}^{2} \left\{ e^{2 \tilde{\mu}} \left( \sqrt{Q} + Q \right) \left( {\bf S} \cdot {\bf v} \right) + \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left[ \left( {\bf \hat{p}} \cdot {\bf v} r \right)
\left( {\bf S} \cdot {\bf n}\right) - \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf S} \cdot {\bf v} \right)\right] \right\}
\end{align*}

We define $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir), $e^{\nu}$ in [this cell](#exp2nu), ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv), $\hat{\bf p} \cdot {\bf v} r$ in [this cell](#pdotvr), ${\bf S} \cdot \boldsymbol{\xi}$ in [this cell](#sdotxi), $\tilde{B}$ in [this cell](#btilde), $Q$ in [this cell](#q), $\tilde{J}$ in [this cell](#jtilde), $\hat{\bf p} \cdot {\bf n}$ in [this cell](#pdotn), ${\bf S} \cdot {\bf n}$ in [this cell](#sdotn), and $\xi^{2}$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HssTerm2 = expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv - expnu*pdotvr*Sdotxi*Btilde)
            + xisq*Btilde*Btilde*(exp2mu*(sp.sqrt(Q) + Q)*Sdotv
            + Jtilde*pdotn*(pdotvr*Sdotn - Jtilde*pdotn*Sdotv))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hssterm3coeff'></a>

## Step 5.d: $H_{\rm SS}$ Term 3 Coefficient \[Back to [top](#toc)\]
$$\label{hssterm3coeff}$$

Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with our definition of $H_{\rm SS}$ Term 3 coefficient in [this cell](#hss), we have

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 3\ coefficient} = \frac{ e^{-3 \tilde{\mu} - \nu} \omega_{\cos\theta} }{ 2 \tilde{B} \left( \sqrt{Q} + 1 \right) \sqrt{Q} }
\end{equation*}

which we write as

\begin{equation*}
    H_{\rm SS}\ {\rm Term\ 3\ coefficient} = \frac{ \omega_{\cos\theta} }{ 2 e^{2 \tilde{\mu}} e^{\tilde{\mu}} e^{\nu} \tilde{B} \left( Q + \sqrt{Q} \right) }.
\end{equation*}

We define $\omega_{\cos\theta}$ in [this cell](#omegacostheta), $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), and $\tilde{B}$ in [this cell](#btilde), $Q$ in [this cell](#q).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HssTerm3coeff = omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q + sp.sqrt(Q)))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hssterm3'></a>

## Step 5.e: $H_{\rm SS}$ Term 3 \[Back to [top](#toc)\]
$$\label{hssterm3}$$

Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with our definition of $H_{\rm SS}$ Term 3 in [this cell](#hss), we have

\begin{align*}
    H_{\rm SS}\ {\rm Term\ 3} &= -e^{2 \left( \tilde{\mu} + \nu \right)} \left( \hat{\bf p} \cdot {\bf \xi} r \right)^{2} \left( {\bf S} \cdot {\bf n} \right) + e^{\tilde{\mu} +\nu} \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf \hat{p}} \cdot {\bf \xi} r \right) \left( {\bf S} \cdot {\bf \xi} \right) \tilde{B} \\
        &\ \ \ \ \ + \left[ \left( {\bf S} \cdot {\bf n} \right) \left( {\bf \hat{p}} \cdot {\bf v} r \right)^{2} - \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf S} \cdot {\bf v} \right) \left( {\bf \hat{p}} \cdot {\bf v} r\right) - e^{2 \tilde{\mu}} \left( 1 + \sqrt{Q} \right) \sqrt{Q} \left( {\bf S} \cdot {\bf n} \right) \xi^{2} \right] \tilde{B}^{2}
\end{align*}

which we write as

\begin{align*}
    H_{\rm SS}\ {\rm Term\ 3} &= e^{\tilde{\mu}} e^{\nu} \left( \hat{\bf p} \cdot {\bf \xi} r \right) \left[ \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf S} \cdot {\bf \xi} \right) \tilde{B} - e^{\tilde{\mu}} e^{\nu} \left( \hat{\bf p} \cdot {\bf \xi} r \right) \left( {\bf S} \cdot {\bf n} \right) \right] \\
        &\ \ \ \ \ + \left\{ \left( {\bf \hat{p}} \cdot {\bf v} r \right) \left[ \left( {\bf S} \cdot {\bf n} \right) \left( {\bf \hat{p}} \cdot {\bf v} r \right) - \tilde{J} \left( {\bf \hat{p}} \cdot {\bf n} \right) \left( {\bf S} \cdot {\bf v} \right) \right] - e^{2 \tilde{\mu}} \left( \sqrt{Q} + Q \right) \left( {\bf S} \cdot {\bf n} \right) \xi^{2} \right\} \tilde{B}^{2}
\end{align*}

We define $e^{\tilde{\mu}}$ in [this cell](#exp2mu), $e^{\nu}$ in [this cell](#exp2nu), $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir), $\tilde{J}$ in [this cell](#jtilde), $\hat{\bf p} \cdot {\bf n}$ in [this cell](#pdotn), ${\bf S} \cdot \boldsymbol{\xi}$ in [this cell](#sdotxi), $\tilde{B}$ in [this cell](#btilde), ${\bf S} \cdot {\bf n}$ in [this cell](#sdotn), $\hat{\bf p} \cdot {\bf v} r$ in [this cell](#pdotvr), ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv), $Q$ in [this cell](#q), and $\xi^{2}$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HssTerm3 = expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde - expmu*expnu*pdotxir*Sdotn)
            + (pdotvr*(Sdotn*pdotvr - Jtilde*pdotn*Sdotv) - exp2mu*(sp.sqrt(Q) + Q)*Sdotn*xisq)*Btilde*Btilde
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hnsterms'></a>

# Step 6: $H_{\rm NS}$ Terms \[Back to [top](#toc)\]
$$\label{hnsterms}$$

We collect here the terms in $H_{\rm NS}$ (defined in [this cell](#hns)).

<a id='betapsum'></a>

## Step 6.a: $\beta p$ sum \[Back to [top](#toc)\]
$$\label{betapsum}$$

We defined the term $\beta p$ sum in [this cell](#hns) as

\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}

Therefore

\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}

We define $\tilde{\omega}_{\rm fd}$ in [this cell](#omegatilde), $\Lambda_{t}$ in [this cell](#lambdat), and $p_{\phi}$ in [this cell](#pphi).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

betapsum = omegatilde*pphi/Lambdat
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='alpha'></a>

## Step 6.b: $\alpha$ \[Back to [top](#toc)\]
$$\label{alpha}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have
\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}

Therefore

\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](#deltat), $\Sigma$ in [this cell](#usigma), and $\Lambda_{t}$ in [this cell](#lambdat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

alpha = sp.sqrt(Deltat*Sigma/Lambdat)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hnsradicand'></a>

## Step 6.c: $H_{\rm NS}$ radicand \[Back to [top](#toc)\]
$$\label{hnsradicand}$$

Recall that we defined $H_{\rm NS}$ radicand in [this cell](#hns) as

\begin{equation*}
    H_{\rm NS}\ {\rm radicand} = \mu^{2} + \underbrace{\gamma^{ij} p_{i} p_{j}}_{\gamma p\ \rm sum} + {\cal Q}_{4}
\end{equation*}

We define $\mu$ in [this cell](#mu), $\gamma p$ sum in [this cell](#gammappsum), and ${\cal Q}_{4}$ in [this cell](#q4).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hnsradicand = 1 + gammappsum + Q4
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='gammappsum'></a>

### Step 6.c.i: $\gamma^{ij} p_{i} p_{j}$ \[Back to [top](#toc)\]
$$\label{gammappsum}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.46), we have

\begin{equation*}
    \gamma^{ij} = g^{ij} - \frac{ g^{ti} g^{tj} }{ g^{tt} }.
\end{equation*}

Combining this result with [BB2010](https://arxiv.org/abs/0912.3517) Equations 5.36, we have

\begin{equation*}
    \gamma^{r\theta} = \gamma^{r\phi} = \gamma^{\theta r} = \gamma^{\theta\phi} = \gamma^{\phi r} = \gamma^{\phi\theta} = 0
\end{equation*}

and

\begin{align*}
    \gamma^{rr} &= g^{rr} = \frac{ \Delta_{r} }{ \Sigma } \\
    \gamma^{\theta\theta} &= g^{\theta\theta} = \frac{ 1 }{ \Sigma } \\
    \gamma^{\phi\phi} &= \frac{ \Sigma }{ \Lambda_{t} \sin^{2} \theta }.
\end{align*}

Therefore

\begin{align*}
    \gamma^{ij} p_{i} p_{j} &= \gamma^{rr} p_{r} p_{r} + \gamma^{\theta\theta} p_{\theta} p_{\theta} + \gamma^{\phi\phi} p_{\phi} p_{\phi} \\
        &= \frac{ \Delta_{r} }{ \Sigma } p_{r}^{2} + \frac{ 1 }{ \Sigma } p_{\theta}^{2} + \frac{ \Sigma }{ \Lambda_{t} \sin^{2} \theta } p_{\phi}^{2}.
\end{align*}

Converting Boyer-Lindquist coordinates to tortoise coordinates (the transformation for which is found in the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2)), we have

\begin{align*}
    p_{r} &= \hat{\bf p} \cdot {\bf n} \\
    p_{\theta} &= \hat{\bf p} \cdot {\bf v} \frac{ r }{ \sin \theta } \\
    p_{\phi} &= \hat{\bf p} \cdot \boldsymbol{\xi} r.
\end{align*}

Therefore

\begin{equation*}
    \gamma^{ij} p_{i} p_{j} = \frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2} + \Sigma^{-1} \left( \hat{\bf p} \cdot {\bf v} \frac{ r }{ \sin \theta } \right)^{2} + \frac{ \Sigma }{ \Lambda_{t} \sin^{2} \theta } \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right)^{2}.
\end{equation*}

We define $\Delta_{r}$ in [this cell](#deltar), $\Sigma$ in [this cell](#sigma), $\hat{\bf p} \cdot {\bf n}$ in [this cell](#pdotn), $\hat{\bf p} \cdot {\bf v} r$ in [this cell](#pdotvr), $\sin^{2} \theta$ in [this cell](#sin2theta), $\Lambda_{t}$ in [this cell](#lambdat), and $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

gammappsum = Deltar/Sigma*pdotn*pdotn + 1/Sigma*pdotvr*pdotvr/sin2theta + Sigma/Lambdat/sin2theta*pdotxir*pdotxir
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='q4'></a>

### Step 6.c.ii: ${\cal Q}_{4}$ \[Back to [top](#toc)\]
$$\label{q4}$$

From [T2012](https://arxiv.org/abs/1202.0790) Equation (15),

\begin{equation*}
    {\cal Q}_{4} \propto \frac{ p_{r^{*}}^{4} }{ r^{2} } \left( r^{2} + \chi_{\rm Kerr}^{2} \right)^{4}.
\end{equation*}

We denote $p_{r^{*}}$ by prT.  Converting from tortoise coordinates to physical coordinates(the transformation for which is found in the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2)), we find

\begin{equation*}
    {\cal Q}_{4} = \frac{ prT^{4} }{ r^{2} } z_{3}
\end{equation*}

where $z_{3}$ is found in [D2000](https://arxiv.org/abs/gr-qc/0005034) Equation (4.34):

\begin{equation*}
    z_{3} = 2 \left( 4 - 3 \nu \right) \nu.
\end{equation*}

In the notation of [BB2010](https://arxiv.org/abs/0912.3517), $\nu = \eta$ (see discussion after [T2012](https://arxiv.org/abs/1202.0790) Equation (2)).  Thus

\begin{equation*}
    {\cal Q}_{4} = 2 prT^{4} u^{2} \left( 4 - 3 \eta \right) \eta.
\end{equation*}

We define prT in [this cell](#prt), $u$ in [this cell](#u), and $\eta$ in [this cell](#eta) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Q4 = 2*prT*prT*prT*prT*u*u*(4 - 3*eta)*eta
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hdterms'></a>

# Step 7: The $H_{\rm D}$ Terms \[Back to [top](#toc)\]
$$\label{hdterms}$$

Recall we defined $H_{\rm D}$ in [this cell](#hd) as

\begin{equation*}
    H_{\rm D} = H_{\rm D}\ {\rm coeffecient} * H_{\rm D}\ {\rm sum}.
\end{equation*}

In this step, we break down each of $H_{\rm D}$ coefficient (defined in [this cell](#hdcoeff)) and $H_{\rm D}$ sum (defined in [this cell](#hdsum)).

<a id='hdcoeff'></a>

## Step 7.a: $H_{\rm D}$ Coefficient \[Back to [top](#toc)\]
$$\label{hdcoeff}$$

From our definition of $H_{\rm D}$ in [this cell](#hd), we have

\begin{equation*}
    H_{\rm D}\ {\rm coefficient} = \frac{ \mu }{ 2 M r^{3} },
\end{equation*}

and recalling the definition of [$\eta$](#eta) we'll write

\begin{equation*}
    H_{\rm D}\ {\rm coefficient} = \frac{ \eta }{ 2 r^{3} }.
\end{equation*}

We define $\eta$ in [this cell](#eta) and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hdcoeff = sp.Rational(1,2)/(r*r*r)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hdsum'></a>

## Step 7.b: $H_{\rm D}$ Sum \[Back to [top](#toc)\]
$$\label{hdsum}$$

From our definition of $H_{\rm D}$ in [this cell](#hd), we have

\begin{align*}
    H_{\rm D}\ {\rm sum} &= \left( \delta^{ij} - 3 n^{i} n^{j} \right) S^{*}_{i} S^{*}_{j} \\
        &= \underbrace{\delta^{ij} S^{*}_{i} S^{*}_{j}}_{\rm Term\ 1} - \underbrace{3 n^{i} n^{j} S^{*}_{i} S^{*}_{j}}_{\rm Term\ 2}.
\end{align*}

We compute $H_{\rm D}$ Term 1 in [this cell](#hdsumterm1) and $H_{\rm D}$ Term 2 in [this cell](#hdsumterm2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Hdsum = HdsumTerm1 - HdsumTerm2
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hdsumterm1'></a>

### Step 7.b.i: $H_{\rm D}$ Sum Term 1 \[Back to [top](#toc)\]
$$\label{hdsumterm1}$$

From our definition of $H_{\rm D}$ sum Term 1 in [this cell](#hdsum), we have

\begin{equation*}
    H_{\rm D}\ {\rm sum\ Term\ 1} = \delta^{ij} S^{*}_{i} S^{*}_{j}
\end{equation*}

where $\delta^{ij}$ is the Kronecker delta:

\begin{equation*}
    \delta_{ij} = \left\{ \begin{array}{cc}
        0, & i \not= j \\
        1, & i = j. \end{array} \right.
\end{equation*}

Thus we have

\begin{equation*}
    H_{\rm D}\ {\rm sum\ Term\ 1} = S^{*}_{1} S^{*}_{1} + S^{*}_{2} S^{*}_{2} + S^{*}_{3} S^{*}_{3}
\end{equation*}

We define ${\bf S}^{*}$ in [this cell](#hreal_spin_combos).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HdsumTerm1 = Sstar1*Sstar1 + Sstar2*Sstar2 + Sstar3*Sstar3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hdsumterm2'></a>

### Step 7.b.ii: $H_{\rm D}$ Sum Term 2 \[Back to [top](#toc)\]
$$\label{hdsumterm2}$$

From our definition of $H_{\rm D}$ sum Term 2 in [this cell](#hdsum), we have

\begin{align*}
    H_{\rm D}\ {\rm sum\ Term\ 2} &= 3 n^{i} n^{j} S^{*}_{i} S^{*}_{j} \\
        &= 3 \left( {\bf S}^{*} \cdot {\bf n} \right)^{2} \\
\end{align*}


We define ${\bf S}^{*} \cdot {\bf n}$ in [this cell](#sstardotn).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

HdsumTerm2 = 3*Sstardotn*Sstardotn
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='dotproducts'></a>

# Step 8: Common Dot Products \[Back to [top](#toc)\]
$$\label{dotproducts}$$

What follows are definitions of many common dot products.

<a id='sdotxi'></a>

## Step 8.a: ${\bf S} \cdot \boldsymbol{\xi}$ \[Back to [top](#toc)\]
$$\label{sdotxi}$$

We have

\begin{equation*}
    {\bf S} \cdot \boldsymbol{\xi} = S^{1} \xi^{1} + S^{2} \xi^{2} + S^{3} \xi^{3}
\end{equation*}

We define $\xi$ in [this cell](#xi).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sdotxi = S1*xi1 + S2*xi2 + S3*xi3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sdotv'></a>

## Step 8.b: ${\bf S} \cdot {\bf v}$ \[Back to [top](#toc)\]
$$\label{sdotv}$$

We have

\begin{equation*}
    {\bf S} \cdot {\bf v} = S^{1} v^{1} + S^{2} v^{2} + S^{3} v^{3}.
\end{equation*}

We define ${\bf v}$ in [this cell](#v).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sdotv = S1*v1 + S2*v2 + S3*v3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sdotn'></a>

## Step 8.c: ${\bf S} \cdot {\bf n}$ \[Back to [top](#toc)\]
$$\label{sdotn}$$

We have

\begin{equation*}
    {\bf S} \cdot {\bf n} = S^{1} n^{1} + S^{2} n^{2} + S^{3} n^{3}.
\end{equation*}

We define ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sdotn = S1*n1 + S2*n2 + S3*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sdotskerrhat'></a>

## Step 8.d: ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$ \[Back to [top](#toc)\]
$$\label{sdotskerrhat}$$

We have

\begin{equation*}
    {\bf S} \cdot \hat{\bf S}_{\rm Kerr} = S^{1} \hat{S}_{\rm Kerr}^{1} + S^{2} \hat{S}_{\rm Kerr}^{2} + S^{3} \hat{S}_{\rm Kerr}^{3}.
\end{equation*}

We define $\hat{\bf S}_{\rm Kerr}$ in [this cell](#skerrhat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

SdotSkerrhat = S1*Skerrhat1 + S2*Skerrhat2 + S3*Skerrhat3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sstardotn'></a>

## Step 8.e: ${\bf S}^{*} \cdot {\bf n}$ \[Back to [top](#toc)\]
$$\label{sstardotn}$$

We have

\begin{equation*}
    {\bf S}^{*} \cdot {\bf n} = {\bf S}^{*}_{1} n_{1} + {\bf S}^{*}_{2} n_{2} + {\bf S}^{*}_{3} n_{3}.
\end{equation*}

We define ${\bf S}^{*}$ in [this cell](#sstar) and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sstardotn = Sstar1*n1 + Sstar2*n2 + Sstar3*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hreal_spin_combos'></a>

# Step 9: $H_{\rm real}$ Spin Combination ${\bf S}^{*}$ \[Back to [top](#toc)\]
$$\label{hreal_spin_combos}$$

We collect here terms defining and containing ${\bf S}^{*}$.

<a id='sstar'></a>

## Step 9.a: ${\bf S}^{*}$ \[Back to [top](#toc)\]
$$\label{sstar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.63):

\begin{equation*}
    {\bf S}^{*} = \boldsymbol{\sigma}^{*} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma^{*}}.
\end{equation*}

We define $\boldsymbol{\sigma}^{*}$ in [this cell](#sigmastar) and $\boldsymbol{\Delta}_{\sigma^{*}}$ in [this cell](#deltasigmastar).

Please note: after normalization, ${\bf S} = {\bf S}^{*}$.  See [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.26).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

S1 = Sstar1
S2 = Sstar2
S3 = Sstar3
Sstar1 = sigmastar1 + Deltasigmastar1
Sstar2 = sigmastar2 + Deltasigmastar2
Sstar3 = sigmastar3 + Deltasigmastar3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltasigmastar'></a>

## Step 9.b: $\boldsymbol{\Delta}_{\sigma^{*}}$ \[Back to [top](#toc)\]
$$\label{deltasigmastar}$$

We can write $\boldsymbol{\Delta}_{\sigma^{*}}$ as

\begin{equation*}
    \boldsymbol{\Delta}_{\sigma^{*}} = \boldsymbol{\sigma}^{*} \left( \boldsymbol{\sigma}^{*}\ {\rm coefficient} \right) + \boldsymbol{\sigma} \left( \boldsymbol{\sigma}\ {\rm coefficient} \right)
\end{equation*}

For further dissection, see $\boldsymbol{\sigma}^{*}$ in [this cell](#sigmastar), $\boldsymbol{\sigma}^{*}$ coefficient in [this cell](#sigmastarcoeff), $\boldsymbol{\sigma}$ in [this cell](#sigma), and $\boldsymbol{\sigma}$ coefficient in [this cell](#sigmacoeff).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff
Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff
Deltasigmastar3 = sigmastar3*sigmastarcoeff + sigma3*sigmacoeff
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmastarcoeff'></a>

## Step 9.c: $\boldsymbol{\sigma}^{*}$ coefficient \[Back to [top](#toc)\]
$$\label{sigmastarcoeff}$$

We will break down $\boldsymbol{\sigma}^{*}\ {\rm coefficient}$ into two terms:

\begin{equation*}
    \boldsymbol{\sigma}^{*}\ {\rm coefficient} = \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 1} + \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 2}
\end{equation*}

We compute $\boldsymbol{\sigma}^{*}$ coefficient Term 1 in [this cell](#sigmastarcoeffterm1) and $\boldsymbol{\sigma}^{*}$ coefficient Term 2 in [this cell](#sigmastarcoeffterm2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmastarcoeff = sigmastarcoeffTerm1 + sigmastarcoeffTerm2
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmastarcoeffterm1'></a>

### Step 9.c.i: $\boldsymbol{\sigma}^{*}$ Coefficient Term 1 \[Back to [top](#toc)\]
$$\label{sigmastarcoeffterm1}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (51) with $b_{0} = 0$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), where what is listed below is the coefficient on $\boldsymbol{\sigma}^{*}$:

\begin{align*}
    \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 1} &= \frac{7}{6} \eta \frac{M}{r} + \frac{1}{3} \eta \left( Q - 1 \right) - \frac{5}{2} \eta \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \\
        &= \frac{ \eta }{ 12 } \left( 14 \frac{ M }{ r } + 4 \left( Q - 1 \right) - 30 \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right)
\end{align*}

We group together and compute $Q-1$ in [this cell](#q) and $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](#drsipn2); we define $r$ in [this cell](#r), $\eta$ in [this cell](#eta), and $M$ in [this cell](#m) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmastarcoeffTerm1 = eta/12*(14/r + 4*Qminus1 - 30*DrSipn2)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmastarcoeffterm2'></a>

### Step 9.c.ii: $\boldsymbol{\sigma}^{*}$ Coefficient Term 2 \[Back to [top](#toc)\]
$$\label{sigmastarcoeffterm2}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (52) with all $b_{i} = 0$, $i \in \left\{0, 1, 2, 3\right\}$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), and just the coefficient on $\boldsymbol{\sigma}^{*}$.  In the LALSuite code this is the variable 'sMultiplier1':

\begin{align*}
    \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 2} &= \frac{1}{36} \left( 353 \eta - 27 \eta^2 \right) \left( \frac{M}{r} \right)^{2} + \frac{5}{3} \left( 3 \eta^2 \right) \frac{ \Delta_{r}^{2} }{ \Sigma^{2} } \left( {\bf n} \cdot \hat{\bf p} \right)^{4} \\
            &\ \ \ \ \ + \frac{1}{72} \left( -23 \eta -3 \eta^{2} \right) \left( Q - 1 \right)^{2} + \frac{1}{36} \left( -103 \eta + 60 \eta^{2} \right) \frac{M}{r} \left( Q - 1 \right) \\
            &\ \ \ \ \ + \frac{1}{12} \left( 16 \eta - 21 \eta^{2} \right) \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \left( Q - 1 \right) + \frac{1}{12} \left( 47 \eta - 54 \eta^{2} \right) \frac{M}{r} \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \\
        &= \frac{ \eta }{ 72 r^{2} } \left[ \left( 706 - 54 \eta \right) M^{2} + 360 \eta r^{2} \frac{ \Delta_{r}^{2} }{ \Sigma^{2} } \left( {\bf n} \cdot \hat{\bf p} \right)^{4} + r^{2} \left( -23 - 3 \eta \right) \left( Q - 1 \right)^{2} + \left( -206 + 120 \eta \right) M r \left( Q - 1 \right) \right. \\
            &\ \ \ \ \ + \left. \left( 96 - 126 \eta \right) r^{2} \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \left( Q - 1 \right) + \left( 282 - 324 \eta \right) M r \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right] \\
        &= \frac{ \eta }{ 72 r^{2} } \left[ 706 + r \left( -206 M \left( Q - 1 \right) + 282 M \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} + r \left( Q -1 \right) \left( 96 \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} - 23 \left( Q - 1 \right) \right) \right) \right. \\
            &\ \ \ \ \ + \left. \eta \left( -54 M^{2} + r \left( 120 M \left( Q -1 \right) - 324 M \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right.\right.\right. \\
            &\ \ \ \ \ + \left.\left.\left. r \left( 360 \frac{ \Delta_{r}^{2} }{ \Sigma^{2} } \left( {\bf n} \cdot \hat{\bf p} \right)^{4} + \left( Q - 1 \right) \left( -126 \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} - 3 \left( Q - 1 \right) \right) \right)\right) \right) \right]
\end{align*}

We define $r$ in [this cell](#r), $\eta$ in [this cell](#eta), and $M$ in [this cell](#m); we group together and define $Q - 1$ in [this cell](#q), and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](#drsipn2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmastarcoeffTerm2 = eta/(72*r*r)*(706 + r*(-206*Qminus1 + 282*DrSipn2 + r*Qminus1*(96*DrSipn2 - 23*Qminus1))
                                    + eta*(-54 + r*(120*Qminus1 - 324*DrSipn2
                                    + r*(360*DrSipn2*DrSipn2 + Qminus1*(-126*DrSipn2 - 3*Qminus1)))))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmacoeff'></a>

## Step 9.d: $\boldsymbol{\sigma}$ coefficient \[Back to [top](#toc)\]
$$\label{sigmacoeff}$$

We will break down $\boldsymbol{\sigma}\ {\rm coefficient}$ into three terms:

\begin{equation*}
    \boldsymbol{\sigma}\ {\rm coefficient} = \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 1} + \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 2} + \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 3}
\end{equation*}

We compute $\boldsymbol{\sigma}$ coefficient Term 1 in [this cell](#sigmacoeffterm1), $\boldsymbol{\sigma}$ coefficient Term 2 in [this cell](#sigmacoeffterm2), and $\boldsymbol{\sigma}$ coefficient Term 3 in [this cell](#sigmacoeffterm3).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmacoeff = sigmacoeffTerm1 + sigmacoeffTerm2 + sigmacoeffTerm3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmacoeffterm1'></a>

### Step 9.d.i: $\boldsymbol{\sigma}$ Coefficient Term 1 \[Back to [top](#toc)\]
$$\label{sigmacoeffterm1}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (51) with $a_{0} = 0$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), where what is listed below is the coefficient on $\boldsymbol{\sigma}$:

\begin{align*}
    \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 1} &= -\frac{2}{3} \eta \frac{ M }{ r } + \frac{1}{4} \eta \left( Q - 1 \right) - 3 \eta \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \\
        &= \frac{ \eta }{ 12 } \left( -8 \frac{ M }{ r } + 3 \left( Q - 1 \right) - 36 \smash[b]{\underbrace{ \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} }_{\rm DrSipn2}} \vphantom{\underbrace{a}_{b}} \right)
\end{align*}

We define $\eta$ in [this cell](#eta), $M$ in [this cell](#m), $Q-1$ in [this cell](#q), and $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](#drsipn2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmacoeffTerm1 = eta/12*(-8/r + 3*Qminus1 - 36*DrSipn2)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmacoeffterm2'></a>

### Step 9.d.ii: $\boldsymbol{\sigma}$ Coefficient Term 2 \[Back to [top](#toc)\]
$$\label{sigmacoeffterm2}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (52) with all $a_{i} = 0$, $i \in \left\{0, 1, 2, 3\right\}$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), and just the coefficient on $\boldsymbol{\sigma}$:

\begin{align*}
    \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 2} &= \frac{1}{9} \left( -56 \eta -21 \eta^{2} \right) \left( \frac{ M }{ r } \right)^{2} + \frac{5}{24} \left( 27 \eta^{2} \right) \frac{ \Delta_r^{2} }{ \Sigma^{2} } \left( {\bf n} \cdot \hat{\bf p} \right)^{4} \\
            &\ \ \ \ \ + \frac{1}{144} \left(-45 \eta \right) \left( Q - 1 \right)^{2} + \frac{1}{36} \left( -109 \eta + 51 \eta^{2} \right) \frac{ M }{ r } \left( Q - 1 \right) \\
        &\ \ \ \ \ + \frac{1}{24} \left( 6 \eta - 39\eta^{2} \right) \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \left( Q - 1 \right) + \frac{1}{24} \left( -16 \eta - 147 \eta^{2} \right) \frac{ M }{ r } \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \\
        &= \frac{ \eta }{ 144 r^{2} } \left[ -896 M^{2} + r \left( -436 M \left( Q - 1 \right) - 96 M \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right.\right. \\
            &\ \ \ \ \ \left.\left. + r \left( -45 \left( Q - 1 \right)^{2} + 36 \left( Q - 1 \right) \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right) \right) + \eta \left( -336 M^{2} + r \left( 204 M \left( Q -1 \right) - 882 M \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right.\right.\right. \\
            &\ \ \ \ \ \left.\left.\left. + r \left( 810 \frac{ \Delta_{r}^{2} }{ \Sigma^{2} } \left( {\bf n} \cdot \hat{\bf p} \right)^{4}  - 234 \left( Q - 1 \right) \frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \right) \right) \right) \right]
\end{align*}

We define $\eta$ in [this cell](#eta), $M$ in [this cell](#M), $Q - 1$ in [this cell](#q), and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](#drsipn2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmacoeffTerm2 = eta/(144*r*r)*(-896 + r*(-436*Qminus1 - 96*DrSipn2 + r*(-45*Qminus1*Qminus1
                                    + 36*Qminus1*DrSipn2)) + eta*(-336 + r*(204*Qminus1 - 882*DrSipn2
                                    + r*(810*DrSipn2*DrSipn2 - 234*Qminus1*DrSipn2))))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmacoeffterm3'></a>

### Step 9.d.iii: $\boldsymbol{\sigma}$ Coefficient Term 3 \[Back to [top](#toc)\]
$$\label{sigmacoeffterm3}$$

From Section II of [T2014)](https://arxiv.org/abs/1311.2544),

\begin{equation*}
    \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 3} = \eta d_{\rm SO} u^{3}.
\end{equation*}

where $d_{\rm SO}$ is a fitting parameter.  Equation (4.13) of [BL2017](https://arxiv.org/pdf/1611.03703.pdf) gives

\begin{equation*}
    d_{\rm SO} = 147.481449 \chi^{3} \eta^{2} - 568.651115 \chi^3 \eta + 66.198703 \chi^{3} - 343.313058 \chi^{2} \eta + 2495.293427 \chi \eta^{2} - 44.532373
\end{equation*}

We define $\eta$ in [this cell](#eta), $u$ in [this cell](#u), and $\chi$ in [this cell](#chi).  Note that the values have been rounded to agree with those in the LALSuite implementation (see the file LALSimIMRSpinEOBHamiltonian.h).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmacoeffTerm3 = eta*dSO*u*u*u
dSO = 147.481*chi*chi*chi*eta*eta - 568.651*chi*chi*chi*eta + 66.1987*chi*chi*chi - 343.313*chi*chi*eta
    + 2495.29*chi*eta*eta - 44.5324
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='metpotderivs'></a>

# Step 10: Derivatives of the Metric Potential \[Back to [top](#toc)\]
$$\label{metpotderivs}$$

We collect here terms dependent on derivatives of the metric potential (see [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.47)).

<a id='omegar'></a>

## Step 10.a: $\omega_{r}$ \[Back to [top](#toc)\]
$$\label{omegar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47b) we have

\begin{equation*}
    \omega_{r} = \frac{ \Lambda_{t} \tilde{\omega}_{\rm fd}^{\prime} - \Lambda_{t}^{\prime} \tilde{\omega}_{\rm fd} }{ \Lambda_{t}^{2} }.
\end{equation*}

We define $\Lambda_{t}$ in [this cell](#lambdat), $\tilde{\omega}_{\rm fd}^{\prime}$ in [this cell](#omegatildeprm), $\Lambda_{t}^{\prime}$ in [this cell](#lambdatprm), and $\tilde{\omega}_{\rm fd}$ in [this cell](#omegatilde).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

omegar = (Lambdat*omegatildeprm - Lambdatprm*omegatilde)/(Lambdat*Lambdat)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='nur'></a>

## Step 10.b: $\nu_{r}$ \[Back to [top](#toc)\]
$$\label{nur}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47c) we have

\begin{equation*}
    \nu_{r} =  \frac{ r }{ \Sigma } + \frac{ \varpi^{2} \left( \varpi^{2} \Delta^{\prime}_{t} - 4 r \Delta_{t} \right) }{ 2 \Lambda_{t} \Delta_{t} }.
\end{equation*}

We define $r$ in [this cell](#r), $\Sigma$ in [this cell](#usigma), $\varpi^{2}$ in [this cell](#w2), $\Delta_{t}^{\prime}$ in [this cell](#deltatprm), $\Delta_{t}$ in [this cell](#deltat), and $\Lambda_{t}$ in [this cell](#lambdat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

nur = r/Sigma + w2*(w2*Deltatprm - 4*r*Deltat)/(2*Lambdat*Deltat)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='mur'></a>

## Step 10.c: $\mu_{r}$ \[Back to [top](#toc)\]
$$\label{mur}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47d) we have

\begin{equation*}
    \mu_{r} =  \frac{ r }{ \Sigma } - \frac{ 1 }{ \sqrt{ \Delta_{r} } }.
\end{equation*}

We define $r$ in [this cell](#r), $\Sigma$ in [this cell](#usigma), and $\Delta_{r}$ in [this cell](#deltar).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

mur = r/Sigma - 1/sp.sqrt(Deltar)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='omegacostheta'></a>

## Step 10.d: $\omega_{\cos\theta}$ \[Back to [top](#toc)\]
$$\label{omegacostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47f), we have

\begin{equation*}
    \omega_{\cos\theta} = -\frac{ 2 a^{2} \cos\theta \Delta_{t} \tilde{\omega}_{\rm fd} }{ \Lambda_{t}^{2} }.
\end{equation*}

We define $a$ in [this cell](#a), $\cos\theta$ in [this cell](#costheta), $\Delta_{t}$ in [this cell](#deltat), $\tilde{\omega}_{\rm fd}$ in [this cell](#omegatilde), and $\Lambda_{t}$ in [this cell](#lambdat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

omegacostheta = -2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='nucostheta'></a>

## Step 10.e: $\nu_{\cos\theta}$ \[Back to [top](#toc)\]
$$\label{nucostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47g) we have

\begin{equation*}
    \nu_{\cos\theta} = \frac{ a^{2} \varpi^{2} \cos\theta \left( \varpi^{2} - \Delta_{t} \right) }{ \Lambda_{t} \Sigma }.
\end{equation*}

We define $a$ in [this cell](#a), $\varpi^{2}$ in [this cell](#w2), $\cos\theta$ in [this cell](#costheta), $\Delta_{t}$ in [this cell](#deltat), $\Lambda_{t}$ in [this cell](#lambdat), and $\Sigma$ in [this cell](#usigma).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

nucostheta = a*a*w2*costheta*(w2 - Deltat)/(Lambdat*Sigma)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='mucostheta'></a>

## Step 10.f: $\mu_{\cos \theta}$ \[Back to [top](#toc)\]
$$\label{mucostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47h) we have

\begin{equation*}
    \mu_{\cos \theta} =  \frac{ a^{2} \cos \theta }{ \Sigma }.
\end{equation*}

We define $a$ in [this cell](#a), $\cos \theta$ in [this cell](#costheta), and $\Sigma$ in [this cell](#usigma) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

mucostheta = a*a*costheta/Sigma
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='lambdatprm'></a>

## Step 10.g: $\Lambda_{t}^{\prime}$ \[Back to [top](#toc)\]
$$\label{lambdatprm}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.47), we know that the prime notation indicates a derivative with respect to $r$.  Using the definition of $\Lambda_{t}$ in [this cell](#lambdat), we have

\begin{equation*}
    \Lambda_{t}^{\prime} = 4 \left( a^{2} + r^{2} \right) r - a^{2} \Delta_{t}^{\prime} \sin^{2} \theta.
\end{equation*}

We define $a$ in [this cell](#a), $r$ in [this cell](#r), $\Delta_{u}$ in [this cell](#deltau), and $\sin^{2}\theta$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Lambdatprm = 4*(a*a + r*r)*r - 2*a*a*Deltatprm*sin2theta
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='omegatildeprm'></a>

## Step 10.h: $\tilde{\omega}_{\rm fd}^{\prime}$ \[Back to [top](#toc)\]
$$\label{omegatildeprm}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that the prime notation indicates a derivative with respect to $r$.  Using the definition of $\tilde{\omega}_{\rm fd}$ in [this cell](#omegatilde), we have

\begin{equation*}
    \tilde{\omega}_{\rm fd}^{\prime} = 2 a M.
\end{equation*}

We define $a$ in [this cell](#a) and $M$ in [this cell](#m).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

omegatildeprm = 2*a
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='metpots'></a>

# Step 11: The Deformed and Rescaled Metric Potentials \[Back to [top](#toc)\]
$$\label{metpots}$$

We collect here terms of the deformed and scaled metric potentials.  See [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.30)--(5.34) and (5.48)--(5.52).

<a id='omega'></a>

## Step 11.a: $\omega$ \[Back to [top](#toc)\]
$$\label{omega}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.31) we have

\begin{equation*}
    \omega = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }.
\end{equation*}

We define $\tilde{\omega}_{\rm fd}$ in [this cell](#omegatilde) and $\Lambda_{t}$ in [this cell](#lambdat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

omega = omegatilde/Lambdat
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='exp2nu'></a>

## Step 11.b: $e^{2\nu}$ and $e^{\nu}$ \[Back to [top](#toc)\]
$$\label{exp2nu}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.32), we have

\begin{equation*}
    e^{2 \nu} = \frac{ \Delta_{t} \Sigma }{ \Lambda_t }.
\end{equation*}

It follows that

\begin{equation*}
    e^{\nu} = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_t } }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](#deltat), $\Sigma$ in [this cell](#usigma), and $\Lambda_{t}$ in [this cell](#lambdat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

expnu = sp.sqrt(exp2nu)
exp2nu = Deltat*Sigma/Lambdat
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='btilde'></a>

## Step 11.c: $\tilde{B}$ \[Back to [top](#toc)\]
$$\label{btilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.48), we have

\begin{equation*}
    \tilde{B} = \sqrt{ \Delta_{t} }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](#deltat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Btilde = sp.sqrt(Deltat)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='brtilde'></a>

## Step 11.d: $\tilde{B}_{r}$ \[Back to [top](#toc)\]
$$\label{brtilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.49), we have

\begin{equation*}
    \tilde{B}_{r} = \frac{ \sqrt{ \Delta_{r} } \Delta_{t}^{\prime} - 2 \Delta_{t} }{ 2 \sqrt{ \Delta_{r} \Delta_{t} } }.
\end{equation*}

We define $\Delta_{r}$ in [this cell](#deltar), $\Delta_{t}^{\prime}$ in [this cell](#deltatprm), and $\Delta_{t}$ in [this cell](#deltat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Brtilde = (sp.sqrt(Deltar)*Deltatprm - 2*Deltat)/(2*sp.sqrt(Deltar*Deltat))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='exp2mu'></a>

## Step 11.e: $e^{2\tilde{\mu}}$ and $e^{\tilde{\mu}}$ \[Back to [top](#toc)\]
$$\label{exp2mu}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.50), we have

\begin{equation*}
    e^{2 \tilde{\mu}} = \Sigma.
\end{equation*}

It follows that

\begin{equation*}
    e^{\tilde{\mu}} = \sqrt{ \Sigma }.
\end{equation*}


We define $\Sigma$ in [this cell](#usigma).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

expmu = sp.sqrt(exp2mu)
exp2mu = Sigma
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='jtilde'></a>

## Step 11.f: $\tilde{J}$ \[Back to [top](#toc)\]
$$\label{jtilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.51) we have

\begin{equation*}
    \tilde{J} = \sqrt{ \Delta_{r} }.
\end{equation*}

We define $\Delta_{r}$ in [this cell](#deltar) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Jtilde = sp.sqrt(Deltar)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='q'></a>

## Step 11.g: $Q$ \[Back to [top](#toc)\]
$$\label{q}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.52),

\begin{equation*}
    Q = 1 + \underbrace{ \frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2} }_{\rm DrSipn2} + \underbrace{ \frac{ \Sigma }{ \Lambda_t \sin^{2} \theta } }_{\rm Q\ coefficient\ 1} \left( \smash[b]{ \underbrace{ \hat{\bf p} \cdot \boldsymbol{\xi} r }_{\rm pdotxir} } \right)^{2} + \underbrace{ \frac{ 1 }{ \Sigma \sin^{2} \theta } }_{\rm Q\ coefficient\ 2} \left( \smash[b]{ \underbrace{ \hat{\bf p} \cdot {\bf v} r }_{\rm pdotvr} } \right)^{2};
\end{equation*}

We group togther and compute $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ in [this cell](#drsipn2), $\frac{ \Sigma }{ \Lambda_t \sin^{2} \theta }$ in [this cell](#qcoeff1), $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir), $\frac{ 1 }{ \Sigma \sin^{2} \theta }$ in [this cell](#qcoeff2), and $\hat{\bf p} \cdot {\bf v} r$ in [this cell](#pdotvr).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Qminus1 = Q - 1
Q = 1 + DrSipn2 + Qcoeff1*pdotxir*pdotxir + Qcoeff2*pdotvr*pdotvr
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='drsipn2'></a>

### Step 11.g.i: $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ \[Back to [top](#toc)\]
$$\label{drsipn2}$$

We define $\Delta_{r}$ in [this cell](#deltar), $\Sigma$ in [this cell](#usigma), and $\hat{\bf p} \cdot {\bf n}$ in [this cell](#pdotn).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

DrSipn2 = Deltar*pdotn*pdotn/Sigma
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='qcoeff1'></a>

### Step 11.g.ii: Q Coefficient 1 \[Back to [top](#toc)\]
$$\label{qcoeff1}$$

We defined $Q$ coefficient 1 in [this cell](#q) as

\begin{equation*}
    Q\ {\rm coefficient\ 1} = \frac{ \Sigma }{ \Lambda_t \sin^{2} \theta }
\end{equation*}

We define $\Sigma$ in [this cell](#usigma), $\Lambda_{t}$ in [this cell](#lambdat), and $\sin^{2} \theta$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Qcoeff1 = Sigma/(Lambdat*sin2theta)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='qcoeff2'></a>

### Step 11.g.iii: Q Coefficient 2 \[Back to [top](#toc)\]
$$\label{qcoeff2}$$

We defined $Q$ coefficient 2 in [this cell](#q) as

\begin{equation*}
    Q\ {\rm coefficient\ 2} = \frac{ 1 }{ \Sigma \sin^{2} \theta }
\end{equation*}

We define $\Sigma$ in [this cell](#usigma) and $\sin^{2} \theta$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Qcoeff2 = 1/(Sigma*sin2theta)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='tort'></a>

# Step 12: Tortoise Terms \[Back to [top](#toc)\]
$$\label{tort}$$

We collect here terms related to the conversion from Boyer-Lindquist coordinates to tortoise coordinates.  Details of the conversation are given in the appendix of [P2010](https://arxiv.org/abs/0912.3466v2).

<a id='pphi'></a>

## Step 12.a: $p_{\phi}$ \[Back to [top](#toc)\]
$$\label{pphi}$$

From the discussion preceding [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.41), the phi component of the tortoise momentum $p_{\phi}$ is given by

\begin{equation*}
    p_{\phi} = \hat{\bf p} \cdot \boldsymbol{\xi} r.
\end{equation*}

We define $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](#pdotxir).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pphi = pdotxir
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotvr'></a>

## Step 12.b: $\hat{\bf p} \cdot {\bf v} r$ \[Back to [top](#toc)\]
$$\label{pdotvr}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot {\bf v} r = \left( \hat{p}_{1} v_{1} + \hat{p}_{2} v_{2} + \hat{p}_{3} v_{3} \right) r
\end{equation*}

We define $\hat{\bf p}$ in [this cell](#hatp), ${\bf v}$ in [this cell](#v), and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotvr = (phat1*v1 + phat2*v2 + phat3*v3)*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotn'></a>

## Step 12.c: $\hat{\bf p} \cdot {\bf n}$ \[Back to [top](#toc)\]
$$\label{pdotn}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot {\bf n} = \hat{p}_{1} n_{1} + \hat{p}_{2} n_{2} + \hat{p}_{3} n_{3}
\end{equation*}

We define $\hat{\bf p}$ in [this cell](#hatp) and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotn = phat1*n1 + phat2*n2 + phat3*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotxir'></a>

## Step 12.d: $\hat{\bf p} \cdot \boldsymbol{\xi} r$ \[Back to [top](#toc)\]
$$\label{pdotxir}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot \boldsymbol{\xi} r = \left( \hat{p}_{1} \xi_{1} + \hat{p}_{2} \xi_{2} + \hat{p}_{3} \xi_{3} \right) r
\end{equation*}

We define $\hat{\bf p}$ in [this cell](#hatp), $\boldsymbol{\xi}$ in [this cell](#xi), and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotxir = (phat1*xi1 + phat2*xi2 + phat3*xi3)*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hatp'></a>

## Step 12.e: $\hat{\bf p}$ \[Back to [top](#toc)\]
$$\label{hatp}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.41), we have $\hat{\bf p} = {\bf p}/m$ where $m$ is the mass of a nonspinning test particle and ${\bf p}$ is *conjugate* momentum.  Following Lines 319--321 of LALSimIMRSpinEOBHamiltonianPrec.c, we convert the Boyer-Lindquist momentum ${\bf p}$ to the tortoise momentum (see the appendix of [P2010](https://arxiv.org/abs/0912.3466v2)) via

\begin{align*}
    \hat{\bf p} = {\bf p} + {\rm prT} \left( 1 - \frac{1}{\rm csi1} \right) {\bf n}
\end{align*}

We define prT in [this cell](#prt), csi1 in [this cell](#csi1), and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

phat1 = p1 + prT*(1 - 1/csi1)*n1
phat2 = p2 + prT*(1 - 1/csi1)*n2
phat3 = p3 + prT*(1 - 1/csi1)*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='prt'></a>

## Step 12.f: prT \[Back to [top](#toc)\]
$$\label{prt}$$

The first component of the momentum vector, after conversion to tortoise coordinates (see the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2)), is

\begin{align*}
    {\rm prT} = {\rm csi2}\left( p_{1} n_{1} + p_{2} n_{2} + p_{3} n_{3} \right)
\end{align*}

We define csi2 in [this cell](#csi2) and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

prT = csi2*(p1*n1 + p2*n2 + p3*n3)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='csi2'></a>

## Step 12.g: csi2 \[Back to [top](#toc)\]
$$\label{csi2}$$

From the transformation to tortoise coordinates in the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2),

\begin{equation*}
    {\rm csi2} = 1 + \left( \frac{1}{2} - \frac{1}{2}{\rm sign}\left( \frac{3}{2} - \tau \right) \right) \left( {\rm csi} - 1 \right)
\end{equation*}

We define csi in [this cell](#csi); $\tau$ is a tortoise coordinate ($\tau \in \left\{ 0, 1 ,2 \right\}$).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

csi2 = 1 + (sp.Rational(1,2) - sp.Rational(1,2)*sp.sign(sp.Rational(3,2) - tortoise))*(csi - 1)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='csi1'></a>

## Step 12.h: csi1 \[Back to [top](#toc)\]
$$\label{csi1}$$

From the transformation to tortoise coordinates in the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2),

\begin{equation*}
    {\rm csi1} = 1 + \left( 1 - \left\lvert 1 - \tau \right\rvert \right) \left( {\rm csi} - 1 \right)
\end{equation*}

We define csi in [this cell](#csi); $\tau$ is a tortoise coordinate ($\tau \in \left\{ 0, 1 ,2 \right\}$).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

csi1 = 1 + (1 - sp.abs(1-tortoise))*(csi - 1)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='csi'></a>

## Step 12.i: csi \[Back to [top](#toc)\]
$$\label{csi}$$

From the transformation to tortoise coordinates in the Appendix of [P2010](https://arxiv.org/abs/0912.3466v2),

\begin{equation*}
    {\rm csi} = \frac{ \sqrt{ \Delta_{t} \Delta_{r} } }{ \varpi^{2} }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](#deltat), $\Delta_{r}$ in [this cell](#deltar), and $\varpi^{2}$ in [this cell](#w2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

csi = sp.sqrt(Deltar*Deltat)/w2
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='metric'></a>

# Step 13: Metric Terms \[Back to [top](#toc)\]
$$\label{metric}$$

We collect here terms used to define the deformed Kerr metric.  See [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.38)--(5.40) and (5.71)--(5.75).

<a id='lambdat'></a>

## Step 13.a: $\Lambda_{t}$ \[Back to [top](#toc)\]
$$\label{lambdat}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.39),

\begin{equation*}
    \Lambda_{t} = \varpi^{4} - a^{2} \Delta_{t} \sin^{2} \theta.
\end{equation*}

We define $\varpi^{2}$ in [this cell](#w2), $a$ in [this cell](#a), $\Delta_{t}$ in [this cell](#deltat), and $\sin^{2}\theta$ in [this cell](#sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Lambdat = w2*w2 - a*a*Deltat*sin2theta
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltar'></a>

## Step 13.b: $\Delta_{r}$ \[Back to [top](#toc)\]
$$\label{deltar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.38),

\begin{equation*}
    \Delta_{r} = \Delta_{t} D^{-1}.
\end{equation*}

We define $\Delta_{t}$ in [this cell](#deltat) and $D^{-1}$ in [this cell](#dinv).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltar = Deltat*Dinv
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltat'></a>

## Step 13.c: $\Delta_{t}$ \[Back to [top](#toc)\]
$$\label{deltat}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.71), we have

\begin{equation*}
    \Delta_{t} = r^{2} \Delta_{u}.
\end{equation*}

We define $\Delta_{u}$ in [this cell](#deltau) and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltat = r*r*Deltau
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltatprm'></a>

## Step 13.d: $\Delta_{t}^{\prime}$ \[Back to [top](#toc)\]
$$\label{deltatprm}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that the prime notation indicates a derivative with respect to $r$.  Using the definition of [$\Delta_{t}$](#deltat), we have

\begin{equation*}
    \Delta_{t}^{\prime} = 2 r \Delta_{u} + r^{2} \Delta_{u}^{\prime}.
\end{equation*}

We define $\Delta_{u}$ in [this cell](#deltau) and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltatprm = 2*r*Deltau + r*r*Deltauprm
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltau'></a>

## Step 13.e: $\Delta_{u}$ \[Back to [top](#toc)\]
$$\label{deltau}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.73), we have

\begin{equation*}
    \Delta_u = \bar{\Delta}_{u} \left[ \smash[b]{\underbrace{ 1 + \eta \Delta_{0} + \eta \log \left( 1 + {\rm logarg} \right) }_{\Delta_{u}\ {\rm calibration\ term}}} \vphantom{\underbrace{1}_{n}} \right]
\end{equation*}

We compute $\bar{\Delta}_{u}$  in [this cell](#deltaubar) and $\Delta_{u}$ calibration term and logarg in [this cell](#deltaucalib).  From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that primes denote derivatives with respect to $r$.  We have

\begin{equation*}
    \Delta_u = \bar{\Delta}^{\prime}_{u} \left( \Delta_{u}\ {\rm calibration\ term} \right) + \bar{\Delta}_{u} \left( \Delta_{u}\ {\rm calibration\ term} \right)^{\prime}
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltauprm = Deltaubarprm*Deltaucalib + Deltaubar*Deltaucalibprm
Deltau = Deltaubar*Deltaucalib
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltaubar'></a>

### Step 13.e.i: $\bar{\Delta}_{u}$ \[Back to [top](#toc)\]
$$\label{deltaubar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.75), we have

\begin{equation*}
    \bar{\Delta}_u = \frac{ a^{2} u^{2} }{ M^{2} } + \frac{ 1 }{ \eta K - 1 } \left( 2 u + \frac{ 1 }{ \eta K - 1 } \right).
\end{equation*}

We define $a$ in [this cell](#a), $u$ in [this cell](#u), $M$ in [this cell](#m), $\eta$ in [this cell](#eta), and $K$ in [this cell](#k).  From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that primes denote derivatives with respect to $r$.  We have

\begin{equation*}
    \bar{\Delta}^{\prime}_u = \frac{ -2 a^{2} u^{3} }{ M^{2} } - \frac{ 2 u^{2} }{ \eta K - 1 }.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltaubarprm = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
Deltaubar = a*a*u*u + (2*u + 1/etaKminus1)/etaKminus1
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='deltaucalib'></a>

### Step 13.e.ii: $\Delta_{u}$ Calibration Term \[Back to [top](#toc)\]
$$\label{deltaucalib}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.73), we have

\begin{align*}
    \Delta_u\ {\rm calibration\ term} &= 1 + \eta \Delta_{0} + \eta \log \left( 1 + \Delta_{1} u + \Delta_{2} u^{2} + \Delta_{3} u^{3} + \Delta_{4} u^{4} \right) \\
        &= 1 + \eta \left[ \Delta_{0} + \log \left( 1 + \Delta_{1} u + \Delta_{2} u^{2} + \Delta_{3} u^{3} + \Delta_{4} u^{4} \right) \right].
\end{align*}

In [T2014](https://arxiv.org/pdf/1311.2544.pdf) Equation (2) an additional term is and is defined in Equation (A2) of [this paper](https://arxiv.org/abs/1608.01907v2).  We then have

\begin{equation*}
    \Delta_u\ {\rm calibration\ term} = 1 + \eta \left[ \Delta_{0} + \log \left( 1 + \Delta_{1} u + \Delta_{2} u^{2} + \Delta_{3} u^{3} + \Delta_{4} u^{4} + \Delta_{5} u^{5} \right) \right].
\end{equation*}

<font color='red'>In the LALSuite code itself (see LALSimIMRSpinEOBHamiltonianPrec.c line 274 on Git commit a70b43d), there's one more term ($\Delta_{5\ell}$), for which documentation is elusive.</font>  That brings us to

\begin{equation*}
    \Delta_u\ {\rm calibration\ term} = 1 + \eta \left[ \Delta_{0} + \log \left( 1 + \underbrace{ \Delta_{1} u + \Delta_{2} u^{2} + \Delta_{3} u^{3} + \Delta_{4} u^{4} + \Delta_{5} u^{5} + \Delta_{5\ell} u^{5} \ln\left(u\right) }_{ \rm logarg } \right) \right].
\end{equation*}

Note our notation for logarg.  We define $u$ in [this cell](#u), $\eta$ in [this cell](#eta), and the calibration coefficients $\Delta_{i}$, $i \in \left\{0, 1, 2, 3, 4\right\}$, in [this cell](#calib_coeffs).

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that primes denote derivatives with respect to $r$.  We have
\begin{equation*}
    \left( \Delta_u\ {\rm calibration\ term} \right)^{\prime} = \frac{ -\eta u^{2} \left( \Delta_{1} + 2 \Delta_{2} u + 3 \Delta_{3} u^{2} + 4 \Delta_{4} u^{3} + 5 \Delta_{5} u^{4} + 5 \Delta_{5\ell} u^{4} \ln\left( u \right) + \Delta_{5\ell} u^{5} u^{-1} \right) }{ 1 + {\rm logarg} }.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Deltaucalibprm = -eta*u*u*(Delta1 + u*(2*Delta2 + u*(3*Delta3
                            + u*(4*Delta4 + u*(5*(Delta5 + Delta5l*sp.log(u)))))))/(1 + logarg)
Deltaucalib = 1 + eta*(Delta0 + sp.log(1 + logarg))
logarg = u*(Delta1 + u*(Delta2 + u*(Delta3 + u*(Delta4 + u*(Delta5 + Delta5l*sp.log(u))))))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='calib_coeffs'></a>

### Step 13.e.iii: Calibration Coefficients $\Delta_{i}$, $i \in \left\{0, 1, 2, 3, 4\right\}$ \[Back to [top](#toc)\]
$$\label{calib_coeffs}$$

The first term in the brackets of [SH2016](https://arxiv.org/abs/1608.01907) Equation (A2c) is

\begin{equation*}
        \Delta_{5\ell} = \frac{64}{5} \left( \eta K - 1 \right)^{2}.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta5l = etaKminus1*etaKminus1*sp.Rational(64,5)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


We combine [SH2016](https://arxiv.org/abs/1608.01907) Equation (A2c) with [BL2017](https://arxiv.org/abs/1611.03703) Equation (2.3) (the last term in our expression) to find

\begin{align*}
    \Delta_{5} &= \left( \eta K - 1 \right)^{2} \left( -\frac{4237}{60} + \frac{128}{5}\gamma + \frac{2275}{512} \pi^{2} - \frac{1}{3} a^{2} \left\{ \Delta_{1}^{3} - 3 \Delta_{1} \Delta_{2} + 3 \Delta_{3} \right\} \right. \\
        &\ \ \ \ \ - \frac{ \Delta_{1}^{5} - 5 \Delta_{1}^{3} \Delta_{2} + 5 \Delta_{1} \Delta_{2}^{2} + 5 \Delta_{1}^{2} \Delta_{3} - 5 \Delta_{2} \Delta_{3} - 5 \Delta_{1} \Delta_{4} }{ 5 \left( \eta K - 1 \right)^{2} } \\
         &\left.\ \ \ \ \  + \frac{ \Delta_{1}^{4} - 4 \Delta_{1}^{2}  \Delta_{2} + 2 \Delta_{2}^{2} + 4 \Delta_{1} \Delta_{3} - 4 \Delta_{4} }{ 2\left( \eta K - 1 \right) } + \frac{256}{5} \log(2) + \left\{ \frac{41\pi^2}{32} - \frac{221}{6} \right\} \eta \right)
\end{align*}

Note that we have exlcuded the first term in the brackets in (A2c); this is the term $\Delta_{5\ell}$ which we defined in [this cell](#calib_coeffs).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta5 = etaKminus1*etaKminus1*(sp.Rational(-4237,60) + sp.Rational(128,5)*EMgamma
                        + sp.Rational(2275,512)*sp.pi*sp.pi - sp.Rational(1,3)*a*a*(Delta1*Delta1*Delta1
                        - 3*Delta1*Delta2 + 3*Delta3) - (Delta1*Delta1*Delta1*Delta1*Delta1
                        - 5*Delta1*Delta1*Delta1*Delta2 + 5*Delta1*Delta2*Delta2 + 5*Delta1*Delta1*Delta3
                        - 5*Delta2*Delta3 - 5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)
                        + (Delta1*Delta1*Delta1*Delta1 - 4*Delta1*Delta1*Delta2 + 2*Delta2*Delta2
                        + 4*Delta1*Delta3 - 4*Delta4)/(2*etaKminus1) + sp.Rational(256,5)*sp.log(2)
                        + (sp.Rational(41,32)*sp.pi*sp.pi - sp.Rational(221,6))*eta)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.81), we have

\begin{align*}
    \Delta_{4} &= \frac{1}{12} \left\{ 6 \frac{ a^{2} }{ M^{2} } \left( \Delta_{1}^{2} - 2 \Delta_{2} \right) \left( \eta K - 1 \right)^{2} + 3 \Delta_{1}^{4} - 8 \left( \eta K - 1 \right) \Delta_{1}^{3} - 12 \Delta_{2} \Delta_{1}^{2} + 12 \left[ 2 \left( \eta K - 1 \right) \Delta_{2} + \Delta_{3} \right] \Delta_{1} \right.\\
        &\left.\ \ \ \ \ + 12 \left( \frac{94}{3} - \frac{41}{32} \pi^{2} \right) \left( \eta K - 1 \right)^{2} + 6 \left[ \Delta_{2}^{2} - 4 \Delta_{3} \left( \eta K - 1 \right) \right] \right\} \\
\end{align*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta4 = sp.Rational(1,12)*(6*a*a*(Delta1*Delta1 - 2*Delta2)*etaKminus1*etaKminus1 + 3*Delta1*Delta1*Delta1*Delta1
                        - 8*etaKminus1*Delta1*Delta1*Delta1 -12*Delta2*Delta1*Delta1 + 12*(2*etaKminus1*Delta2
                        + Delta3)*Delta1 + 12*(sp.Rational(94,3)
                        - sp.Rational(41,32)*sp.pi*sp.pi)*etaKminus1*etaKminus1 + 6*(Delta2*Delta2
                        - 4*Delta3*etaKminus1))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.80), we have

\begin{align*}
    \Delta_{3} &= \frac{1}{3} \left[ -\Delta_{1}^{3} + 3 \left( \eta K - 1 \right) \Delta_{1}^{2} + 3 \Delta_{2} \Delta_{1} - 6 \left( \eta K - 1 \right) \left( -\eta K + \Delta_{2} + 1 \right) - 3 \frac{ a^{2} }{ M^{2} } \left( \eta K - 1 \right)^{2} \Delta_{1} \right] \\
        &= -\frac{1}{3}\Delta_{1}^{3} + \left( \eta K - 1 \right) \Delta_{1}^{2} + \Delta_{2} \Delta_{1} - 2 \left( \eta K - 1 \right) \left( \Delta_{2}- \left( \eta K - 1 \right) \right) - \frac{ a^{2} }{ M^{2} } \left( \eta K - 1 \right)^{2} \Delta_{1}
\end{align*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta3 = -sp.Rational(1,3)*Delta1*Delta1*Delta1 + etaKminus1*Delta1*Delta1 + Delta2*Delta1
                        -2*etaKminus1*(Delta2 - etaKminus1) - a*a*etaKminus1*etaKminus1*Delta1
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.79, we have

\begin{equation*}
    \Delta_{2} = \frac{1}{2} \Delta_{1} \left( -4 \eta K + \Delta_{1} + 4 \right) - \frac{ a^{2} }{ M^{2} } \left( \eta K - 1 \right)^{2} \Delta_{0}\\
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta2 = sp.Rational(1,2)*Delta1*(Delta1 - 4*etaKminus1) - a*a*etaKminus1*etaKminus1*Delta0
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.78), we have

\begin{equation*}
    \Delta_{1} = -2 \left( \eta K - 1 \right) \left( K + \Delta_{0} \right)
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta1 = -2*etaKminus1*(K + Delta0)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.77), we have

\begin{equation*}
    \Delta_{0} = K \left( \eta K - 2 \right)
\end{equation*}

We define $K$ and $\eta K-1$ in [this cell](#k), $\eta$ in [this cell](#eta), $a$ in [this cell](#a), and $M$ in [this cell](#m).  Note that the constant $\gamma$ is the Euler-Mascheroni, and the value is taken from the [LALSuite documentation](https://lscsoft.docs.ligo.org/lalsuite/lal/group___l_a_l_constants__h.html).  In the Python code, we denote $\gamma$ by EMgamma.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Delta0 = K*(eta*K - 2)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='k'></a>

### Step 13.e.iv: $K$ \[Back to [top](#toc)\]
$$\label{k}$$

The calibration constant $K$ is defined in [BL2017](https://arxiv.org/abs/1611.03703) Section IV, Equations (4.8) and (4.12) (equations copied from BL2017 source code):

\begin{equation*}
  \left.K\right|_{\chi=0} = 267.788247\,  \nu^3 -126.686734\,  \nu^2 + 10.257281\,\nu  + 1.733598,
\end{equation*}

\begin{align*}
  K =& - 59.165806\,\chi^3\nu^3 - 0.426958\,\chi^3\nu + 1.436589\,\chi^3 + 31.17459\,\chi^2\nu^3 + 6.164663\,\chi^2\nu^2 - 1.380863\,\chi^2 \\
    & - 27.520106\,\chi \nu^3 + 17.373601\,\chi\nu^2 + 2.268313\,\chi\nu - 1.62045\,\chi +\left.K\right|_{\chi=0}
\end{align*}

Here $\left.K\right|_{\chi=0}$ denotes the nonspining fits and $\nu = \eta$ (see discussion after [BL2017](https://arxiv.org/pdf/1611.03703.pdf) Equation (2.1)).  Furthermore, most coefficients are rounded in the LALSuite code; below we round to match the values therein (see the file LALSimIMRSpinEOBHamiltonian.h).  The term $\eta K - 1$ is sufficiently common that we also define it:

\begin{equation*}
    {\rm etaKminus1} = \eta K - 1.
\end{equation*}

We define $\eta$ in [this cell](#eta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

etaKminus1 = eta*K - 1
K = - 59.1658*chi*chi*chi*eta*eta*eta - 0.426958*chi*chi*chi*eta + 1.43659*chi*chi*chi
    + 31.1746*chi*chi*eta*eta*eta + 6.16466*chi*chi*eta*eta - 1.38086*chi*chi - 27.5201*chi*eta*eta*eta
    + 17.3736*chi*eta*eta + 2.26831*chi*eta - 1.62045*chi + Kchi0
Kchi0 = 267.788*eta*eta*eta -126.687*eta*eta + 10.2573*eta  + 1.7336
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='chi'></a>

### Step 13.e.v: $\chi$ \[Back to [top](#toc)\]
$$\label{chi}$$

The augmented spin $\chi$ is defined in [OB2020](https://arxiv.org/abs/2004.09442) Equation (3.7) (where it is denoted $\tilde{\chi}$; note that $\nu = \eta$ from the first paragraph of Section II):

\begin{equation*}
    \chi = \frac{{\bf S}_{\rm Kerr} \cdot \hat{\bf L}}{1 - 2\eta} + \alpha\frac{({\bf S}_{1}^{\perp} + {\bf S}_{2}^{\perp}) \cdot {\bf S}_{\rm Kerr}}{ \left\lvert{\bf S}_{\rm Kerr} \right\rvert (1 - 2\eta)}.
\end{equation*}

From the discussion after this equation, we take $\alpha = \frac{1}{2}$.  We define ${\bf S}^{\perp} = {\bf S}_{1}^{\perp} + {\bf S}_{2}^{\perp}$ in [this cell](#sperp), ${\bf L}$ in [this cell](#orb_momentum), $\left\lvert{\bf S}_{\rm Kerr} \right\rvert$ in [this cell](#skerrmag), and $\eta$ in [this cell](#eta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

chi = (Skerr1*Lhat1 + Skerr2*Lhat2 + Skerr3*Lhat3)/(1 - 2*eta)
        + sp.Rational(1,2)*(Sperp1*Skerr1 + Sperp2*Skerr2 + Sperp3*Skerr3)/(Skerrmag*(1. - 2.*eta))
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='omegatilde'></a>

## Step 13.f: $\tilde{\omega}_{\rm fd}$ \[Back to [top](#toc)\]
$$\label{omegatilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.40), we have

\begin{equation*}
    \tilde{\omega}_{\rm fd} = 2 a M r + \omega_{1}^{\rm fd} \eta \frac{ a M^{3} }{ r } + \omega_{2}^{\rm fd} \eta \frac{ M a^{3} }{ r }.
\end{equation*}

From discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (6.7), we set $\omega_{1}^{\rm fd} = \omega_{2}^{\rm fd} = 0$.  Thus

\begin{equation*}
    \tilde{\omega}_{\rm fd} = 2 a M r.
\end{equation*}

We define $a$ in [this cell](#a), $M$ in [this cell](#m), and $r$ in [this cell](#r) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

omegatilde = 2*a*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='dinv'></a>

## Step 13.g: $D^{-1}$ \[Back to [top](#toc)\]
$$\label{dinv}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.83),

\begin{equation*}
    D^{-1} = 1 + \log \left[ 1 + 6 \eta u^{2} + 2 \left( 26 - 3 \eta \right) \eta u^{3} \right].
\end{equation*}

We define $\eta$ in [this cell](#eta) and $u$ in [this cell](#u).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Dinv = 1 + sp.log(1 + 6*eta*u*u + 2*(26 - 3*eta)*eta*u*u*u)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='coord'></a>

# Step 14: Terms Dependent on Coordinates \[Back to [top](#toc)\]
$$\label{coord}$$

We collect here terms directly dependent on the coordinates.  See [BB2010](https://arxiv.org/abs/0912.3517) Equations (4.5) and (4.6).

<a id='usigma'></a>

## Step 14.a: $\Sigma$ \[Back to [top](#toc)\]
$$\label{usigma}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.5), we have

\begin{equation*}
    \Sigma = r^{2} + a^{2} \cos^{2} \theta.
\end{equation*}

We define $r$ in [this cell](#r), $a$ in [this cell](#a), and $\cos \theta$ in [this cell](#costheta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sigma = r*r + a*a*costheta*costheta
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='w2'></a>

## Step 14.b: $\varpi^{2}$ \[Back to [top](#toc)\]
$$\label{w2}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.7),

\begin{equation*}
    \varpi^{2} = a^{2} + r^{2}.
\end{equation*}

We define $a$ in [this cell](#a) and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

w2 = a*a + r*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sin2theta'></a>

## Step 14.d: $\sin^{2} \theta$ \[Back to [top](#toc)\]
$$\label{sin2theta}$$

Using a common trigonometric idenitity,

\begin{equation*}
    \sin^{2} \theta = 1 - \cos^{2} \theta.
\end{equation*}

We define $\cos \theta$ in [this cell](#costheta).  Note that by construction (from discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.52))

\begin{equation*}
    \xi^{2} = \sin^{2} \theta.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

xisq = sin2theta
sin2theta = 1 - costheta*costheta
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='costheta'></a>

## Step 14.e: $\cos \theta$ \[Back to [top](#toc)\]
$$\label{costheta}$$

From the discussion in [BB2010](https://arxiv.org/abs/0912.3517) after Equation (5.52) (noting that ${\bf e}_{3} = \hat{\bf S}_{\rm Kerr}$),

\begin{equation*}
    \cos \theta = {\bf e}_{3} \cdot {\bf n} = {\bf e}_{3}^{1} n^{1} + {\bf e}_{3}^{2} n^{2} + {\bf e}_{3}^{3} n^{3}.
\end{equation*}

We define ${\bf e}_{3}$ in [this cell](#e3) and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

costheta = e31*n1 + e32*n2 + e33*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='vectors'></a>

# Step 15: Important Vectors \[Back to [top](#toc)\]
$$\label{vectors}$$

We collect the vectors common for computing $H_{\rm real}$ (defined in [this cell](#hreal)) below.

<a id='v'></a>

## Step 15.a: ${\bf v}$ \[Back to [top](#toc)\]
$$\label{v}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.39), we have

\begin{equation*}
    {\bf v} = {\bf n} \times \boldsymbol{\xi}.
\end{equation*}

We define ${\bf n}$ in [this cell](#n) and $\boldsymbol{\xi}$ in [this cell](#xi).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

v1 = n2*xi3 - n3*xi2
v2 = n3*xi1 - n1*xi3
v3 = n1*xi2 - n2*xi1
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='xi'></a>

## Step 15.b: $\boldsymbol{\xi}$ \[Back to [top](#toc)\]
$$\label{xi}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.38), we have

\begin{equation*}
    \boldsymbol{\xi} = {\bf e}_{3} \times {\bf n}.
\end{equation*}

We define ${\bf e}_{3}$ in [this cell](#e3) and ${\bf n}$ in [this cell](#n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

xi1 = e32*n3 - e33*n2
xi2 = e31*n3 + e33*n1
xi3 = e31*n2 - e32*n1
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='e3'></a>

## Step 15.c: ${\bf e}_{3}$ \[Back to [top](#toc)\]
$$\label{e3}$$

From the discussion in [BB2010](https://arxiv.org/abs/0912.3517) after Equation (5.52),

\begin{equation*}
    {\bf e}_{3} = \hat{\bf S}_{\rm Kerr}.
\end{equation*}

We define $\hat{\bf S}_{\rm Kerr}$ in [this cell](#skerrhat).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

e31 = Skerrhat1
e32 = Skerrhat2
e33 = Skerrhat3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='n'></a>

## Step 15.d: ${\bf n}$ \[Back to [top](#toc)\]
$$\label{n}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.37), we have

\begin{equation*}
    {\bf n} = \frac{\bf x }{ r }
\end{equation*}

where ${\bf x} = (x, y, z)$.  We define $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

n1 = x/r
n2 = y/r
n3 = z/r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sperp'></a>

## Step 15.e: ${\bf S}_{\rm perp}$ \[Back to [top](#toc)\]
$$\label{sperp}$$

In the definition of [$\chi$](#chi) we denoted ${\bf S}^{\perp} = {\bf S}_{1}^{\perp} + {\bf S}_{2}^{\perp}$.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Sperp1 = S1perp1 + S2perp1
Sperp2 = S1perp2 + S2perp2
Sperp3 = S1perp3 + S2perp3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


From the discussion after [OB2020](https://arxiv.org/abs/2004.09442) Equation (3.7), we find

\begin{align*}
    {\bf S}_{1}^{\perp} &= {\bf S}_{1} - \left( {\bf S}_{1} \cdot \hat{\bf L} \right) \hat{\bf L} \\
    {\bf S}_{2}^{\perp} &= {\bf S}_{2} - \left( {\bf S}_{2} \cdot \hat{\bf L} \right) \hat{\bf L}
\end{align*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

S2perp1 = S2x - S2dotLhat*Lhat1
S2perp2 = S2y - S2dotLhat*Lhat2
S2perp3 = S2z - S2dotLhat*Lhat3

S1perp1 = S1x - S1dotLhat*Lhat1
S1perp2 = S1y - S1dotLhat*Lhat2
S1perp3 = S1z - S1dotLhat*Lhat3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


For convenience, we here compute ${\bf S}_{1} \cdot \hat{\bf L}$ and ${\bf S}_{2} \cdot \hat{\bf L}$.  We define $\hat{\bf L}$ in [this cell](#orb_momentum).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

S1dotLhat = S1x*Lhat1 + S1y*Lhat2 + S1z*Lhat3
S2dotLhat = S2x*Lhat1 + S2y*Lhat2 + S2z*Lhat3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='orb_momentum'></a>

## Step 15.f: ${\bf L}$ \[Back to [top](#toc)\]
$$\label{orb_momentum}$$

From the discussion after [P2010](https://arxiv.org/abs/0912.3466v2), the orbital angular momentum ${\bf L}$ of the system is given by

\begin{equation*}
    {\bf L} = {\bf x }\times{\bf p }
\end{equation*}

where ${\bf x} = (x, y, z)$ is the position vector and ${\bf p} = (p_{1}, p_{2}, p_{3})$ is the momentum vector.  We denote by $\hat{\bf L}$ the normed orbital angular momentum.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Lhat1 = L1/Lnorm
Lhat2 = L2/Lnorm
Lhat3 = L3/Lnorm

Lnorm = sp.sqrt(L1*L1 + L2*L2 + L3*L3)

L1 = y*p3 - z*p2
L2 = z*p1 - x*p3
L3 = x*p2 - y*p1
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='spin_combos'></a>

# Step 16: Spin Combinations $\boldsymbol{\sigma}$, $\boldsymbol{\sigma}^{*}$, and ${\bf S}_{\rm Kerr}$ \[Back to [top](#toc)\]
$$\label{spin_combos}$$

We collect here various combinations of the spins.

<a id='a'></a>

## Step 16.a: $a$ \[Back to [top](#toc)\]
$$\label{a}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.9), we have

\begin{equation*}
    a = \frac{ \left\lvert {\bf S}_{\rm Kerr} \right\rvert }{ M }.
\end{equation*}

We define $\left\lvert{\bf S}_{\rm Kerr}\right\rvert$ in [this cell](#skerrmag) and $M$ in [this cell](#m).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

a = Skerrmag
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='skerrhat'></a>

## Step 16.b: $\hat{\bf S}_{\rm Kerr}$ \[Back to [top](#toc)\]
$$\label{skerrhat}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.24), we have

\begin{equation*}
    \hat{\bf S}_{\rm Kerr} = \frac{ {\bf S}_{\rm Kerr} }{ \left\lvert {\bf S}_{\rm Kerr} \right\rvert }.
\end{equation*}

We define ${\bf S}_{\rm Kerr}$ in [this cell](#skerr).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Skerrhat1 = Skerr1/Skerrmag
Skerrhat2 = Skerr2/Skerrmag
Skerrhat3 = Skerr3/Skerrmag
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='skerrmag'></a>

## Step 16.c: $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$ \[Back to [top](#toc)\]
$$\label{skerrmag}$$

We have

\begin{equation*}
    \left\lvert {\bf S}_{\rm Kerr} \right\rvert = \sqrt{ {\bf S}_{\rm Kerr}^{1} {\bf S}_{\rm Kerr}^{1} + {\bf S}_{\rm Kerr}^{2} {\bf S}_{\rm Kerr}^{2} + {\bf S}_{\rm Kerr}^{3} {\bf S}_{\rm Kerr}^{3} }.
\end{equation*}

We define ${\bf S}_{\rm Kerr}$ in [this cell](#skerr).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Skerrmag = sp.sqrt(Skerr1*Skerr1 + Skerr2*Skerr2 + Skerr3*Skerr3)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='skerr'></a>

## Step 16.d: ${\bf S}_{\rm Kerr}$ \[Back to [top](#toc)\]
$$\label{skerr}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.64):

\begin{equation*}
    {\bf S}_{\rm Kerr} = \boldsymbol{\sigma} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma}.
\end{equation*}

In [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.67), $\boldsymbol{\Delta}_{\sigma} = 0$.  Thus

\begin{equation*}
    {\bf S}_{\rm Kerr} = \boldsymbol{\sigma}.
\end{equation*}

We define $\boldsymbol{\sigma}$ in [this cell](#sigma).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Skerr1 = sigma1
Skerr2 = sigma2
Skerr3 = sigma3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigma'></a>

## Step 16.e: $\boldsymbol{\sigma}$ \[Back to [top](#toc)\]
$$\label{sigma}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.2):

\begin{equation*}
    \boldsymbol{\sigma} = {\bf S}_{1} + {\bf S}_{2}.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigma1 = S1x + S2x
sigma2 = S1y + S2y
sigma3 = S1z + S2z
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='sigmastar'></a>

## Step 16.f: $\boldsymbol{\sigma}^{*}$ \[Back to [top](#toc)\]
$$\label{sigmastar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.3):

\begin{equation*}
    \boldsymbol{\sigma}^{*} = \frac{ m_{2} }{ m_{1} } {\bf S}_{1} + \frac{ m_{1} }{ m_{2} }{\bf S}_{2}.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

sigmastar1 = m2/m1*S1x + m1/m2*S2x
sigmastar2 = m2/m1*S1y + m1/m2*S2y
sigmastar3 = m2/m1*S1z + m1/m2*S2z
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='fundquant'></a>

# Step 17: Fundamental Quantities \[Back to [top](#toc)\]
$$\label{fundquant}$$

 We collect here fundamental quantities from which we build $H_{\rm real}$ (defined in [this cell](#Hreal)).

<a id='u'></a>

## Step 17.a: $u$ \[Back to [top](#toc)\]
$$\label{u}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.40),

\begin{equation*}
    u = \frac{ M }{ r }.
\end{equation*}

We define $M$ in [this cell](#m) and $r$ in [this cell](#r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

u = 1/r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='r'></a>

## Step 17.b: $r$ \[Back to [top](#toc)\]
$$\label{r}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.52),

\begin{equation*}
    r = \sqrt{ x^{2} + y^{2} + z^{2} }.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

r = sp.sqrt(x*x + y*y + z*z)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='eta'></a>

## Step 17.c: $\eta$ \[Back to [top](#toc)\]
$$\label{eta}$$

From the discussion preceding [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.1),

\begin{equation*}
    \eta = \frac{ \mu }{ M }.
\end{equation*}

We define $\mu$ in [this cell](#mu).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

eta = mu/M
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='mu'></a>

## Step 17.d: $\mu$ \[Back to [top](#toc)\]
$$\label{mu}$$

From the discussion preceding [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.1),

\begin{equation*}
    \mu = \frac{ m_{1} m_{2} }{ M }.
\end{equation*}

We define $M$ in [this cell](#m).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

mu = m1*m2/M
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='m'></a>

## Step 17.e: $M$ \[Back to [top](#toc)\]
$$\label{m}$$

From the discussion preceding [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.1),

\begin{equation*}
    M = m_{1} + m_{2}.
\end{equation*}


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

M = m1 + m2
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='validation'></a>

# Step 18: Validation \[Back to [top](#toc)\]
$$\label{validation}$$

The following code cell reverses the order of the expressions output to SEOBNR/Hamiltonian_on_top.txt and creates a Python function to validate the value of $H_{\rm real}$ against the SEOBNRv3 Hamiltonian value computed in LALSuite git commit bba40f21e9 for command-line input parameters

-M 23 -m 10 -f 20 -X 0.01 -Y 0.02 -Z -0.03 -x 0.04 -y -0.05 -z 0.06.


```python
import numpy as np
import difflib, sys, os

# The subterms in the Hamiltonian expression are sometimes written on more than
# one line for readability in this Jupyter notebook.  We first create a file of
# one-line expressions, Hamiltonian-Hreal_one_line_expressions.txt.
with open(os.path.join(Ccodesdir,"v4P_Hamiltonian-Hreal_one_line_expressions.txt"), "w") as output:
    count = 0
    # Read output of this notebook
    for line in list(open("SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt")):
        # Read the first line
        if count == 0:
            prevline=line
        #Check if prevline is a complete expression
        elif "=" in prevline and "=" in line:
            output.write("%s\n" % prevline.strip('\n'))
            prevline=line
        # Check if line needs to be adjoined to prevline
        elif "=" in prevline and not "=" in line:
            prevline = prevline.strip('\n')
            prevline = (prevline+line).replace(" ","")
        # Be sure to print the last line.
        if count == len(list(open("SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt")))-1:
            if not "=" in line:
                print("ERROR. Algorithm not robust if there is no equals sign on the final line. Sorry.")
                sys.exit(1)
            else:
                output.write("%s" % line)
        count = count + 1

# Now reverse the expressions and write them in a function
# This formulation is used to check that we get a reasonable H_real value
with open(os.path.join(Ccodesdir,"v4P_Hreal_on_bottom.py"), "w") as output:
    output.write("import numpy as np\ndef compute_v4P_Hreal(m1=23., m2=10., EMgamma=0.577215664901532860606512090082402431, tortoise=1, x=2.1242072413581923e+01, y=0., z=0., p1=0., p2=2.1696072000958128e-01, p3=1.0000000000000000e-03, S1x=4.8576675849403119e-03, S1y=9.7153351698806237e-03, S1z=-1.4573002754820936e-02, S2x=3.6730945821854912e-03, S2y=-4.5913682277318639e-03, S2z=5.5096418732782371e-03):\n")
    for line in reversed(list(open("SEOBNR/v4P_Hamiltonian-Hreal_one_line_expressions.txt"))):
        output.write("    %s\n" % line.rstrip().replace("sp.sqrt", "np.sqrt").replace("sp.Rational",
                                "np.divide").replace("sp.abs", "np.abs").replace("sp.log",
                                "np.log").replace("sp.sign", "np.sign").replace("sp.pi",
                                "np.pi"))
    output.write("    return Hreal")

# Now reverse the expressions in a standalone text file
# This formulation is used as a harsher validation check that all expressions agree with a trusted list
with open(os.path.join(Ccodesdir,"v4P_Hamiltonian_expressions.txt-VALIDATION"), "w") as output:
    for line in reversed(list(open("SEOBNR/v4P_Hamiltonian-Hreal_one_line_expressions.txt"))):
        output.write("%s\n" % line.rstrip().replace("sp.sqrt", "np.sqrt").replace("sp.Rational",
                                "np.divide").replace("sp.abs", "np.abs").replace("sp.log",
                                "np.log").replace("sp.sign", "np.sign").replace("sp.pi",
                                "np.pi"))

print("Printing difference between notebook output and a trusted list of expressions...")
# Open the files to compare
file = "v4P_Hamiltonian_expressions.txt"
outfile = "v4P_Hamiltonian_expressions.txt-VALIDATION"

print("Checking file " + outfile)
with open(os.path.join(Ccodesdir,file), "r") as file1, open(os.path.join(Ccodesdir,outfile), "r") as file2:
    # Read the lines of each file
    file1_lines=[]
    file2_lines=[]
    for line in file1.readlines():
        file1_lines.append(line.replace(" ", ""))
    for line in file2.readlines():
        file2_lines.append(line.replace(" ", ""))
    num_diffs = 0
    for line in difflib.unified_diff(file1_lines, file2_lines, fromfile=os.path.join(Ccodesdir,file), tofile=os.path.join(Ccodesdir,outfile)):
        sys.stdout.writelines(line)
        num_diffs = num_diffs + 1
    if num_diffs == 0:
        print("No difference. TEST PASSED!")
    else:
        print("ERROR: Disagreement found with the trusted file. See differences above.")
        sys.exit(1)

# Import the new Hamiltonian function and the trusted Hamiltonian function
import SEOBNR.SEOBNR_v4P_Hamiltonian as Hreal_trusted
import SEOBNR.v4P_Hreal_on_bottom as Hreal_new

# Compute the trusted and new Hamiltonian values; compare; exit if they disagree!
Hreal = Hreal_trusted.compute_v4P_Hreal()
Hreal_temp = Hreal_new.compute_v4P_Hreal()

if(np.abs(Hreal-Hreal_temp)>1e-14):
    print("ERROR. You have broken the Hamiltonian computation!")
    print("Hreal_trusted was ",Hreal)
    print("...and Hreal is now ", Hreal_temp)
    sys.exit(1)
```

    Printing difference between notebook output and a trusted list of expressions...
    Checking file v4P_Hamiltonian_expressions.txt-VALIDATION
    No difference. TEST PASSED!


<a id='latex_pdf_output'></a>

# Step 19: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-SEOBNR_Documentation.pdf](Tutorial-SEOBNR_Documentation.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-SEOBNR_Documentation")
```

    [NbConvertApp] WARNING | pattern 'Tutorial-SEOBNR_Documentation.ipynb' matched no files
    Created Tutorial-SEOBNR_Documentation.tex, and compiled LaTeX file to PDF
        file Tutorial-SEOBNR_Documentation.pdf



```python

```
