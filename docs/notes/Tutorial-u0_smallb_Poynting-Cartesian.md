<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Computing the 4-Velocity Time-Component $u^0$, the Magnetic Field Measured by a Comoving Observer $b^{\mu}$, and the Poynting Vector $S^i$

## Authors: Zach Etienne & Patrick Nelson

## This notebook computes key parameters in relativistic MHD, utilizing the [BSSN.ADMBSSN_tofrom_4metric](../edit/BSSN/ADMBSSN_tofrom_4metric.py) NRPy+ module and the Valencia 3-velocity.  The computations are then validated against the [`u0_smallb_Poynting__Cartesian.py`](../edit/u0_smallb_Poynting__Cartesian/u0_smallb_Poynting__Cartesian.py) NRPy+ module.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated against a trusted code (the hand-written smallbPoynET in WVUThorns_diagnostics, which itself is based on expressions in IllinoisGRMHD... which was validated against the original GRMHD code of the Illinois NR group)

### NRPy+ Source Code for this module: [u0_smallb_Poynting__Cartesian.py](../edit/u0_smallb_Poynting__Cartesian/u0_smallb_Poynting__Cartesian.py)

[comment]: <> (Introduction: TODO)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#u0bu): Computing $u^0$ and $b^{\mu}$
    1. [Step 1.a](#4metric): Compute the 4-metric $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ from the ADM 3+1 variables, using the [`BSSN.ADMBSSN_tofrom_4metric`](../edit/BSSN/ADMBSSN_tofrom_4metric.py) ([**tutorial**](Tutorial-ADMBSSN_tofrom_4metric.ipynb)) NRPy+ module
    1. [Step 1.b](#u0): Compute $u^0$ from the Valencia 3-velocity
    1. [Step 1.c](#uj): Compute $u_j$ from $u^0$, the Valencia 3-velocity, and $g_{\mu\nu}$
    1. [Step 1.d](#gamma): Compute $\gamma=$ `gammaDET` from the ADM 3+1 variables
    1. [Step 1.e](#beta): Compute $b^\mu$
1. [Step 2](#poynting_flux): Defining the Poynting Flux Vector $S^{i}$
    1. [Step 2.a](#g): Computing $g^{i\nu}$
    1. [Step 2.b](#s): Computing $S^{i}$
1. [Step 3](#code_validation): Code Validation against `u0_smallb_Poynting__Cartesian` NRPy+ module
1. [Step 4](#appendix): Appendix: Proving Eqs. 53 and 56 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='u0bu'></a>

# Step 1: Computing $u^0$ and $b^{\mu}$ \[Back to [top](#toc)\]
$$\label{u0bu}$$

First some definitions. The spatial components of $b^{\mu}$ are simply the magnetic field as measured by an observer comoving with the plasma $B^{\mu}_{\rm (u)}$, divided by $\sqrt{4\pi}$. In addition, in the ideal MHD limit, $B^{\mu}_{\rm (u)}$ is orthogonal to the plasma 4-velocity $u^\mu$, which sets the $\mu=0$ component. 

Note also that $B^{\mu}_{\rm (u)}$ is related to the magnetic field as measured by a *normal* observer $B^i$ via a simple projection (Eq 21 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)), which results in the expressions (Eqs 23 and 24 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)):

\begin{align}
\sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
\sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}\\
\end{align}

$B^i$ is related to the actual magnetic field evaluated in IllinoisGRMHD, $\tilde{B}^i$ via

$$B^i = \frac{\tilde{B}^i}{\gamma},$$

where $\gamma$ is the determinant of the spatial 3-metric.

The above expressions will require that we compute
1. the 4-metric $g_{\mu\nu}$ from the ADM 3+1 variables
1. $u^0$ from the Valencia 3-velocity
1. $u_j$ from $u^0$, the Valencia 3-velocity, and $g_{\mu\nu}$
1. $\gamma$ from the ADM 3+1 variables

<a id='4metric'></a>

## Step 1.a: Compute the 4-metric $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ from the ADM 3+1 variables, using the [`BSSN.ADMBSSN_tofrom_4metric`](../edit/BSSN/ADMBSSN_tofrom_4metric.py) ([**tutorial**](Tutorial-ADMBSSN_tofrom_4metric.ipynb)) NRPy+ module \[Back to [top](#toc)\]
$$\label{4metric}$$

We are given $\gamma_{ij}$, $\alpha$, and $\beta^i$ from ADMBase, so let's first compute 

$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$


```python
# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp                         # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par             # NRPy+: Parameter interface
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from outputC import outputC                # NRPy+: Basic C code output functionality
import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions

# Set spatial dimension = 3
DIM=3

thismodule = "smallbPoynET"

# Step 1.a: Compute the 4-metric $g_{\mu\nu}$ and its inverse
#           $g^{\mu\nu}$ from the ADM 3+1 variables, using the
#           BSSN.ADMBSSN_tofrom_4metric NRPy+ module
import BSSN.ADMBSSN_tofrom_4metric as AB4m
gammaDD,betaU,alpha = AB4m.setup_ADM_quantities("ADM")
AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
g4DD = AB4m.g4DD
AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
g4UU = AB4m.g4UU
```

<a id='u0'></a>

## Step 1.b: Compute $u^0$ from the Valencia 3-velocity \[Back to [top](#toc)\]
$$\label{u0}$$

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

In order for the bottom expression to hold true, the left-hand side must be between 0 and 1. Again, this is not guaranteed due to the appearance of numerical errors. In fact, a robust algorithm will not allow $\Gamma^2$ to become too large (which might contribute greatly to the stress energy of a given gridpoint), so let's define $\Gamma_{\rm max}$, the largest allowed Lorentz factor. 

Then our algorithm for computing $u^0$ is as follows:

If
$$R=\gamma_{ij}v^i_{(n)}v^j_{(n)}>1 - \frac{1}{\Gamma_{\rm max}^2},$$ 
then adjust the 3-velocity $v^i$ as follows:

$$v^i_{(n)} = \sqrt{\frac{1 - \frac{1}{\Gamma_{\rm max}^2}}{R}}v^i_{(n)}.$$

After this rescaling, we are then guaranteed that if $R$ is recomputed, it will be set to its ceiling value $R=R_{\rm max} = 1 - \frac{1}{\Gamma_{\rm max}^2}$.

Then, regardless of whether the ceiling on $R$ was applied, $u^0$ can be safely computed via

$$
u^0 = \frac{1}{\alpha \sqrt{1-R}}.
$$


```python
ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU",DIM=3)

# Step 1: Compute R = 1 - 1/max(Gamma)
R = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

GAMMA_SPEED_LIMIT = par.Cparameters("REAL",thismodule,"GAMMA_SPEED_LIMIT",10.0) # Default value based on
                                                                                # IllinoisGRMHD.
                                                                                # GiRaFFE default = 2000.0
Rmax = 1 - 1/(GAMMA_SPEED_LIMIT*GAMMA_SPEED_LIMIT)

rescaledValenciavU = ixp.zerorank1()
for i in range(DIM):
    rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rmax/R)

rescaledu0 = 1/(alpha*sp.sqrt(1-Rmax))
regularu0 =  1/(alpha*sp.sqrt(1-R))

computeu0_Cfunction  = """
/* Function for computing u^0 from Valencia 3-velocity. */
/* Inputs: ValenciavU[], alpha, gammaDD[][], GAMMA_SPEED_LIMIT (C parameter) */
/* Output: u0=u^0 and velocity-limited ValenciavU[] */\n\n"""

computeu0_Cfunction += outputC([R,Rmax],["const double R","const double Rmax"],"returnstring",
                               params="includebraces=False,CSE_varprefix=tmpR,outCverbose=False")

computeu0_Cfunction += "if(R <= Rmax) "
computeu0_Cfunction += outputC(regularu0,"u0","returnstring",
                               params="includebraces=True,CSE_varprefix=tmpnorescale,outCverbose=False")
computeu0_Cfunction += " else "
computeu0_Cfunction += outputC([rescaledValenciavU[0],rescaledValenciavU[1],rescaledValenciavU[2],rescaledu0],
                               ["ValenciavU0","ValenciavU1","ValenciavU2","u0"],"returnstring",
                               params="includebraces=True,CSE_varprefix=tmprescale,outCverbose=False")

print(computeu0_Cfunction)
```

    
    /* Function for computing u^0 from Valencia 3-velocity. */
    /* Inputs: ValenciavU[], alpha, gammaDD[][], GAMMA_SPEED_LIMIT (C parameter) */
    /* Output: u0=u^0 and velocity-limited ValenciavU[] */
    
    const double R = ((ValenciavU0)*(ValenciavU0))*gammaDD00 + 2*ValenciavU0*ValenciavU1*gammaDD01 + 2*ValenciavU0*ValenciavU2*gammaDD02 + ((ValenciavU1)*(ValenciavU1))*gammaDD11 + 2*ValenciavU1*ValenciavU2*gammaDD12 + ((ValenciavU2)*(ValenciavU2))*gammaDD22;
    const double Rmax = 1 - 1/((GAMMA_SPEED_LIMIT)*(GAMMA_SPEED_LIMIT));
    if(R <= Rmax) {
      u0 = 1/(alpha*sqrt(-((ValenciavU0)*(ValenciavU0))*gammaDD00 - 2*ValenciavU0*ValenciavU1*gammaDD01 - 2*ValenciavU0*ValenciavU2*gammaDD02 - ((ValenciavU1)*(ValenciavU1))*gammaDD11 - 2*ValenciavU1*ValenciavU2*gammaDD12 - ((ValenciavU2)*(ValenciavU2))*gammaDD22 + 1));
    }
     else {
      const double tmprescale_1 = sqrt((1 - 1/((GAMMA_SPEED_LIMIT)*(GAMMA_SPEED_LIMIT)))/(((ValenciavU0)*(ValenciavU0))*gammaDD00 + 2*ValenciavU0*ValenciavU1*gammaDD01 + 2*ValenciavU0*ValenciavU2*gammaDD02 + ((ValenciavU1)*(ValenciavU1))*gammaDD11 + 2*ValenciavU1*ValenciavU2*gammaDD12 + ((ValenciavU2)*(ValenciavU2))*gammaDD22));
      ValenciavU0 = ValenciavU0*tmprescale_1;
      ValenciavU1 = ValenciavU1*tmprescale_1;
      ValenciavU2 = ValenciavU2*tmprescale_1;
      u0 = fabs(GAMMA_SPEED_LIMIT)/alpha;
    }
    


<a id='uj'></a>

## Step 1.c: Compute $u_j$ from $u^0$, the Valencia 3-velocity, and $g_{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{uj}$$

The basic equation is

\begin{align}
u_j &= g_{\mu j} u^{\mu} \\
&= g_{0j} u^0 + g_{ij} u^i \\
&= \beta_j u^0 + \gamma_{ij} u^i \\
&= \beta_j u^0 + \gamma_{ij} u^0 \left(\alpha v^i_{(n)} - \beta^i\right) \\
&= u^0 \left(\beta_j + \gamma_{ij} \left(\alpha v^i_{(n)} - \beta^i\right) \right)\\
&= \alpha u^0 \gamma_{ij} v^i_{(n)} \\
\end{align}


```python
u0 = par.Cparameters("REAL",thismodule,"u0",1e300) # Will be overwritten in C code. Set to crazy value to ensure this.

uD = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        uD[j] += alpha*u0*gammaDD[i][j]*ValenciavU[i]
```

<a id='beta'></a>

## Step 1.d: Compute $b^\mu$ \[Back to [top](#toc)\]
$$\label{beta}$$

We compute $b^\mu$ from the above expressions:

\begin{align}
\sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
\sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}\\
\end{align}

$B^i$ is exactly equal to the $B^i$ evaluated in IllinoisGRMHD/GiRaFFE.

Pulling this together, we currently have available as input:
+ $\tilde{B}^i$
+ $u_j$
+ $u^0$,

with the goal of outputting now $b^\mu$ and $b^2$:


```python
M_PI = par.Cparameters("#define",thismodule,"M_PI","")
BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU",DIM=3)

# uBcontraction = u_i B^i
uBcontraction = sp.sympify(0)
for i in range(DIM):
    uBcontraction += uD[i]*BU[i]

# uU = 3-vector representing u^i = u^0 \left(\alpha v^i_{(n)} - \beta^i\right)
uU = ixp.zerorank1()
for i in range(DIM):
    uU[i] = u0*(alpha*ValenciavU[i] - betaU[i])

smallb4U = ixp.zerorank1(DIM=4)
smallb4U[0] = uBcontraction/(alpha*sp.sqrt(4*M_PI))
for i in range(DIM):
    smallb4U[1+i] = (BU[i] + uBcontraction*uU[i])/(alpha*u0*sp.sqrt(4*M_PI))
```

<a id='poynting_flux'></a>

# Step 2: Defining the Poynting Flux Vector $S^{i}$ \[Back to [top](#toc)\]
$$\label{poynting_flux}$$

The Poynting flux is defined in Eq. 11 of [Kelly *et al.*](https://arxiv.org/pdf/1710.02132.pdf) (note that we choose the minus sign convention so that the Poynting luminosity across a spherical shell is $L_{\rm EM} = \int (-\alpha T^i_{\rm EM\ 0}) \sqrt{\gamma} d\Omega = \int S^r \sqrt{\gamma} d\Omega$, as in [Farris *et al.*](https://arxiv.org/pdf/1207.3354.pdf):

$$
S^i = -\alpha T^i_{\rm EM\ 0} = -\alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)
$$



<a id='s'></a>

## Step 2.a: Computing $S^{i}$ \[Back to [top](#toc)\]
$$\label{s}$$

Given $g^{\mu\nu}$ computed above, we focus first on the $g^i{}_{0}$ term by computing 
$$
g^\mu{}_\delta = g^{\mu\nu} g_{\nu \delta},
$$
and then the rest of the Poynting flux vector can be immediately computed from quantities defined above:
$$
S^i = -\alpha T^i_{\rm EM\ 0} = -\alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)
$$


```python
# Step 2.a.i: compute g^\mu_\delta:
g4UD = ixp.zerorank2(DIM=4)
for mu in range(4):
    for delta in range(4):
        for nu in range(4):
            g4UD[mu][delta] += g4UU[mu][nu]*g4DD[nu][delta]

# Step 2.a.ii: compute b_{\mu}
smallb4D = ixp.zerorank1(DIM=4)
for mu in range(4):
    for nu in range(4):
        smallb4D[mu] += g4DD[mu][nu]*smallb4U[nu]

# Step 2.a.iii: compute u_0 = g_{mu 0} u^{mu} = g4DD[0][0]*u0 + g4DD[i][0]*uU[i]
u_0 = g4DD[0][0]*u0
for i in range(DIM):
    u_0 += g4DD[i+1][0]*uU[i]

# Step 2.a.iv: compute b^2, setting b^2 = smallb2etk, as gridfunctions with base names ending in a digit
#          are forbidden in NRPy+.
smallb2etk = sp.sympify(0)
for mu in range(4):
    smallb2etk += smallb4U[mu]*smallb4D[mu]

# Step 2.a.v: compute S^i
PoynSU = ixp.zerorank1()
for i in range(DIM):
    PoynSU[i] = -alpha * (smallb2etk*uU[i]*u_0 + sp.Rational(1,2)*smallb2etk*g4UD[i+1][0] - smallb4U[i+1]*smallb4D[0])
```

<a id='code_validation'></a>

# Step 3: Code Validation against `u0_smallb_Poynting__Cartesian` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for u0, smallbU, smallb2etk, and PoynSU between

1. this tutorial and 
2. the NRPy+ [u0_smallb_Poynting__Cartesian module](../edit/u0_smallb_Poynting__Cartesian/u0_smallb_Poynting__Cartesian.py).


```python
import sys
import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0etc
u0etc.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)

if u0etc.computeu0_Cfunction != computeu0_Cfunction:
    print("FAILURE: u0 C code has changed!")
    sys.exit(1)
else:
    print("PASSED: u0 C code matches!")

for i in range(4):
    print("u0etc.smallb4U["+str(i)+"] - smallb4U["+str(i)+"] = "
              + str(u0etc.smallb4U[i]-smallb4U[i]))

print("u0etc.smallb2etk - smallb2etk = " + str(u0etc.smallb2etk-smallb2etk))

for i in range(DIM):
    print("u0etc.PoynSU["+str(i)+"] - PoynSU["+str(i)+"] = "
              + str(u0etc.PoynSU[i]-PoynSU[i]))
```

    PASSED: u0 C code matches!
    u0etc.smallb4U[0] - smallb4U[0] = 0
    u0etc.smallb4U[1] - smallb4U[1] = 0
    u0etc.smallb4U[2] - smallb4U[2] = 0
    u0etc.smallb4U[3] - smallb4U[3] = 0
    u0etc.smallb2etk - smallb2etk = 0
    u0etc.PoynSU[0] - PoynSU[0] = 0
    u0etc.PoynSU[1] - PoynSU[1] = 0
    u0etc.PoynSU[2] - PoynSU[2] = 0


<a id='appendix'></a>

# Step 4: Appendix: Proving Eqs. 53 and 56 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf)
$$\label{appendix}$$

$u^\mu u_\mu = -1$ implies

\begin{align}
g^{\mu\nu} u_\mu u_\nu &= g^{00} \left(u_0\right)^2 + 2 g^{0i} u_0 u_i + g^{ij} u_i u_j = -1 \\
\implies &g^{00} \left(u_0\right)^2 + 2 g^{0i} u_0 u_i + g^{ij} u_i u_j + 1 = 0\\
& a x^2 + b x + c = 0
\end{align}

Thus we have a quadratic equation for $u_0$, with solution given by

\begin{align}
u_0 &= \frac{-b \pm \sqrt{b^2 - 4 a c}}{2 a} \\
&= \frac{-2 g^{0i}u_i \pm \sqrt{\left(2 g^{0i} u_i\right)^2 - 4 g^{00} (g^{ij} u_i u_j + 1)}}{2 g^{00}}\\
&= \frac{-g^{0i}u_i \pm \sqrt{\left(g^{0i} u_i\right)^2 - g^{00} (g^{ij} u_i u_j + 1)}}{g^{00}}\\
\end{align}

Notice that (Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf))
$$
g^{\mu\nu} = \begin{pmatrix} 
-\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
\frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
\end{pmatrix},
$$
so we have

\begin{align}
u_0 &= \frac{-\beta^i u_i/\alpha^2 \pm \sqrt{\left(\beta^i u_i/\alpha^2\right)^2 + 1/\alpha^2 (g^{ij} u_i u_j + 1)}}{1/\alpha^2}\\
&= -\beta^i u_i \pm \sqrt{\left(\beta^i u_i\right)^2 + \alpha^2 (g^{ij} u_i u_j + 1)}\\
&= -\beta^i u_i \pm \sqrt{\left(\beta^i u_i\right)^2 + \alpha^2 \left(\left[\gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}\right] u_i u_j + 1\right)}\\
&= -\beta^i u_i \pm \sqrt{\left(\beta^i u_i\right)^2 + \alpha^2 \left(\gamma^{ij}u_i u_j + 1\right) - \beta^i\beta^j u_i u_j}\\
&= -\beta^i u_i \pm \sqrt{\alpha^2 \left(\gamma^{ij}u_i u_j + 1\right)}\\
\end{align}

Now, since 

$$
u^0 = g^{\alpha 0} u_\alpha = -\frac{1}{\alpha^2} u_0 + \frac{\beta^i u_i}{\alpha^2},
$$

we get

\begin{align}
u^0 &= \frac{1}{\alpha^2} \left(u_0 + \beta^i u_i\right) \\
&= \pm \frac{1}{\alpha^2} \sqrt{\alpha^2 \left(\gamma^{ij}u_i u_j + 1\right)}\\
&= \pm \frac{1}{\alpha} \sqrt{\gamma^{ij}u_i u_j + 1}\\
\end{align}

By convention, the relativistic Gamma factor is positive and given by $\alpha u^0$, so we choose the positive root. Thus we have derived Eq. 53 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf):

$$
u^0 = \frac{1}{\alpha} \sqrt{\gamma^{ij}u_i u_j + 1}.
$$

Next we evaluate 

\begin{align}
u^i &= u_\mu g^{\mu i} \\
&= u_0 g^{0 i} + u_j g^{i j}\\
&= u_0 \frac{\beta^i}{\alpha^2} + u_j \left(\gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}\right)\\
&= \gamma^{ij} u_j + u_0 \frac{\beta^i}{\alpha^2} - u_j \frac{\beta^i\beta^j}{\alpha^2}\\
&= \gamma^{ij} u_j + \frac{\beta^i}{\alpha^2} \left(u_0  - u_j \beta^j\right)\\
&= \gamma^{ij} u_j - \beta^i u^0,\\
\implies v^i &= \frac{\gamma^{ij} u_j}{u^0} - \beta^i
\end{align}

which is equivalent to Eq. 56 in [Duez *et al* (2005)](https://arxiv.org/pdf/astro-ph/0503420.pdf). Notice in the last step, we used the above definition of $u^0$.

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-u0_smallb_Poynting-Cartesian.pdf](Tutorial-u0_smallb_Poynting-Cartesian.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-u0_smallb_Poynting-Cartesian")
```

    [pandoc warning] Duplicate link reference `[comment]' "source" (line 22, column 1)
    Created Tutorial-u0_smallb_Poynting-Cartesian.tex, and compiled LaTeX file
        to PDF file Tutorial-u0_smallb_Poynting-Cartesian.pdf

