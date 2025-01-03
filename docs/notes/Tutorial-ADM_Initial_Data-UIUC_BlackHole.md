<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [UIUC Black Hole](https://arxiv.org/abs/1001.4077) Initial data

## Authors: Terrence Pierre Jacques, Zach Etienne, & Ian Ruchlin

### Formatting improvements courtesy Brandon Clark

## This notebook introduces a module for setting up UIUC Black Hole initial data for studying highly spinning black holes ([Liu, Etienne, & Shapiro, PRD 80 121503, 2009](https://arxiv.org/abs/1001.4077)). Using the [Exact ADM Spherical-or-Cartesian-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_Exact_ADM_Spherical_or_Cartesian_to_BSSNCurvilinear.ipynb), it supports transitions from spherical to any defined coordinate system.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated to exhibit convergence to zero of the Hamiltonian and momentum constraint violation at the expected order to the exact solution (see plots at bottom of [the exact initial data validation start-to-finish tutorial notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb); momentum constraint violation in non-$\phi$ directions is zero), and all quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).

### NRPy+ Source Code for this module: [BSSN/UIUCBlackHole.py](../edit/BSSN/UIUCBlackHole.py)

## Introduction:
UIUC black holes have the advantage of finite coordinate radius in the maximal spin limit. It is therefore excellent for studying very highly spinning black holes. This module sets the UIUC black hole at the origin. 

<a id='toc'></a>

# Table of Contents:  
$$\label{toc}$$

1. [Step 1](#initializenrpy): Set up the needed NRPy+ infrastructure and declare core gridfunctions
1. [Step 2](#bl_radius): The Boyer-Lindquist Radius
    1. [Step 2.a](#define_inner_outer_radii): Define the inner and outer radii
    1. [Step 2.b](#define_bl_radius): Define the Boyer-Lindquist radius
1. [Step 3](#line_element): Define the line element, and extract components of $\gamma_{ij}$
1. [Step 4](#extrinsic_curvature): Define and construct nonzero components of the extrinsic curvature $K_{ij}$
1. [Step 5](#lapse_shift): Construct Lapse function $\alpha$ and components of shift vector $\beta$
1. [Step 6](#code_validation): Code Validation against `BSSN.UIUCBlackHole` NRPy+ module
1. [Step 7](#latex_pdf_output) Output this notebook to $\LaTeX$-formatted PDF file



<a id='initializenrpy'></a>

# Step 1: Set up the needed NRPy+ infrastructure and declare core gridfunctions \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

First, we will import the core modules of Python/NRPy+ and specify the main gridfunctions that we will need.
Second, we set some basic NRPy+ parameters. E.g., set the spatial dimension parameter to 3.

**Inputs for initial data**:

* The black hole mass $M$.
* The dimensionless spin parameter $\chi = a/M$

**Additional variables needed for spacetime evolution**:

* Desired coordinate system Boyer-Lindquist coordinates $(r_{BL}, \theta, \phi)$
<br>
* Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=1$ and $\beta^i=B^i=0$. $\alpha = \psi^{-2}$ will yield much better behavior, but the conformal factor $\psi$ depends on the desired *destination* coordinate system (which may not be spherical coordinates).


```python
# Step P0: Load needed modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
r,th,ph = sp.symbols('r th ph', real=True)

thismodule = "UIUCBlackHole"

# Step 0: Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Step 1: Set psi, the conformal factor:

# The UIUC initial data represent a Kerr black hole with mass M
#  and dimensionless spin chi in UIUC quasi-isotropic coordinates,
#   see https://arxiv.org/abs/1001.4077
# Input parameters:
M,chi = par.Cparameters("REAL", thismodule, ["M","chi"],[1.0,0.99])

# Spin per unit mass
a = M*chi
```

<a id='bl_radius'></a>

# Step 2: The Boyer-Lindquist Radius \[Back to [top](#toc)\]
$$\label{bl_radius}$$


<a id='define_inner_outer_radii'></a>

## Step 2.a: Defining the Inner and Outer Radii \[Back to [top](#toc)\]
$$\label{define_inner_outer_radii}$$



Boyer-Lindquist radii of the outer (+) and inner (−) horizons of the BH, defined under equation 1 in [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077) as 
$$ r_{\pm} = M \pm \sqrt{M^2 - a^2}$$


```python
# Defined under equation 1 in Liu, Etienne, & Shapiro (2009)
# https://arxiv.org/pdf/1001.4077.pdf

# Boyer - Lindquist outer horizon
rp = M + sp.sqrt(M**2 - a**2)
# Boyer - Lindquist inner horizon
rm = M - sp.sqrt(M**2 - a**2)
```

<a id='define_bl_radius'></a>

## Step 2.b: Define the Boyer-Lindquist Radius  \[Back to [top](#toc)\]
$$\label{define_bl_radius}$$

Define $r_{BL}$, equation 11 of [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077), using the radial coordinate $r$:

$$  r_{BL} = r \left( 1 + \frac{r_+}{4r}\right)^2.  $$


```python
# Boyer - Lindquist radius in terms of UIUC radius
# Eq. 11
# r_{BL} = r * ( 1 + r_+ / 4r )^2
rBL = r*(1 + rp / (4*r))**2
```

Quantities used to calculate the spatial metric $\gamma_{ij}$, found under equation 2 of [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077):
$$  \Sigma = r_{BL}^2 + a^2 \cos^2 \theta, $$

$$ \Delta = r_{BL}^2 - 2Mr_{BL} + a^2,  $$ 

$$  A = \left(r_{BL}^2 + a^2\right)^2 - \Delta a^2 \sin^2 \theta.  $$ 


```python
# Expressions found below Eq. 2
# Sigma = r_{BL}^2 + a^2 cos^2 theta
SIG = rBL**2 + a**2*sp.cos(th)**2

# Delta = r_{BL}^2 - 2Mr_{BL} + a^2
DEL = rBL**2 - 2*M*rBL + a**2

# A = (r_{BL}^2 + a^2)^2 - Delta a^2 sin^2 theta
AA = (rBL**2 + a**2)**2 - DEL*a**2*sp.sin(th)**2
```

<a id='line_element'></a>

# Step 3: Define the Line element and extract components of $\gamma_{ij}$ \[Back to [top](#toc)\]
$$\label{line_element}$$

The line element, defined in equation 13 of [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077):

$$ ds^2 = \frac{\Sigma\left(r + \frac{r_+}{4}\right)^2 } {r^3 \left(r_{BL} - r_- \right)} dr^2 + \Sigma d\theta^2  +  \frac{ A \sin^2 \theta  }{\Sigma} d\phi^2  $$


```python
# *** The ADM 3-metric in spherical basis ***
gammaDD = ixp.zerorank2()
# Declare the nonzero components of the 3-metric (Eq. 13):

# ds^2 = Sigma (r + r_+/4)^2 / ( r^3 (r_{BL} - r_- ) * dr^2 +
# Sigma d theta^2  +  (A sin^2 theta) / Sigma  *  d\phi^2

gammaDD[0][0] = ((SIG*(r + rp/4)**2)/(r**3*(rBL - rm)))
gammaDD[1][1] = SIG
gammaDD[2][2] = AA/SIG*sp.sin(th)**2
```

<a id='extrinsic_curvature'></a>

# Step 4: Define and construct nonzero components of extrinsic curvature $K_{ij}$ \[Back to [top](#toc)\]
$$\label{extrinsic_curvature}$$



Nonzero components of the extrinsic curvature, equation 14 of [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077):

$$ K_{r\phi} = K_{\phi r} = \frac{Ma\sin^2\theta}{\Sigma\sqrt{A\Sigma}} \ 
    \left[3r^4_{BL} + 2a^2 r^2_{BL} - a^4 - a^2 \left(r^2_{BL} - a^2\right) \sin^2 \theta\right] \
    \left(1 + \frac{r_+}{4r}\right) \frac{1}{\sqrt{r(r_{BL} - r_-)}}  $$


```python
# *** The physical trace-free extrinsic curvature in spherical basis ***
# Nonzero components of the extrinsic curvature K, given by
# Eq. 14 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:
KDD     = ixp.zerorank2() # K_{ij} = 0 for these initial data


# K_{r phi} = K_{phi r} = (Ma sin^2 theta) / (Sigma sqrt{A Sigma}) *
#     [3r^4_{BL} + 2a^2 r^2_{BL} - a^4 - a^2 (r^2_{BL} - a^2) sin^2 theta] *
#     (1 + r_+ / 4r) (1 / sqrt{r(r_{BL} - r_-)})

KDD[0][2] = KDD[2][0] = (M*a*sp.sin(th)**2)/(SIG*sp.sqrt(AA*SIG))*\
                (3*rBL**4 + 2*a**2*rBL**2 - a**4- a**2*(rBL**2 - a**2)*\
                 sp.sin(th)**2)*(1 + rp/(4*r))*1/sp.sqrt(r*(rBL - rm))
```

Nonzero components of the extrinsic curvature, equation 15 of [Liu, Etienne, & Shapiro (2009)](https://arxiv.org/abs/1001.4077):

$$ K_{\theta\phi} = K_{\phi\theta} = -\frac{2a^3 Mr_{BL}\cos\theta \sin^3\theta} {\Sigma \sqrt{A\Sigma} } \left(r - \frac{r_+}{4}\right) \sqrt {\frac{r_{BL} - r_-}{r} }  $$


```python
# Components of the extrinsic curvature K, given by
# Eq. 15 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:

# K_{theta phi} = K_{phi theta} = -(2a^3 Mr_{BL} cos theta sin^3 theta) /
#         (Sigma sqrt{A Sigma}) x (r - r_+ / 4) sqrt{(r_{BL} - r_-) / r }

KDD[1][2] = KDD[2][1] = -((2*a**3*M*rBL*sp.cos(th)*sp.sin(th)**3)/ \
                (SIG*sp.sqrt(AA*SIG)))*(r - rp/4)*sp.sqrt((rBL - rm)/r)
```

<a id='lapse_shift'></a>

# Step 5: Construct Lapse function $\alpha$ and components of shift vector $\beta$ \[Back to [top](#toc)\]
$$\label{lapse_shift}$$

$$\alpha=W=e^{-2\phi}$$, and

$$\beta^i=B^i=0$$

where $\phi$ is the BSSN conformal factor.


```python
betaU = ixp.zerorank1() # We generally choose \beta^i = 0 for these initial data
BU    = ixp.zerorank1() # We generally choose B^i = 0 for these initial data

# Finally set alpha. We generally choose alpha = 1/psi**2 (psi = BSSN conformal factor)
#                    for these initial data
import BSSN.BSSN_quantities as Bq  # Sets default for EvolvedConformalFactor_cf
import BSSN.BSSN_in_terms_of_ADM as BitoA
import reference_metric as rfm
rfm.reference_metric()
try:
    cf_type = par.parval_from_str("EvolvedConformalFactor_cf")
except:
    print("UIUCBlackHole Error: Must set BSSN_quantities::EvolvedConformalFactor_cf;")
    print("                     the lapse is set in terms of the BSSN conformal factor")
    sys.exit(1)

BitoA.cf_from_gammaDD(gammaDD)
cf = BitoA.cf

# Let's choose alpha = 1/psi**2 (psi = BSSN conformal factor) for these initial data,
# where psi = exp(phi); chi = 1/psi**4; W = 1/psi**2
if cf_type == "phi":
    alpha = sp.exp(-2*cf)
elif cf_type == "chi":
    alpha = sp.sqrt(cf)
elif cf_type == "W":
    alpha = cf
else:
    print("Error EvolvedConformalFactor_cf type = \""+cf_type+"\" unknown.")
    sys.exit(1)

# Validated against original SENR: KDD[0][2], KDD[1][2], gammaDD[2][2], gammaDD[0][0], gammaDD[1][1]
# print(sp.mathematica_code(gammaDD[1][1]))
```

<a id='code_validation'></a>

# Step 6: Code Validation against `BSSN.UIUCBlackHole` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$


Here, as a code validation check, we verify agreement in the SymPy expressions for UIUC black hole initial data between

1. this tutorial and 
2. the NRPy+ [BSSN.UIUCBlackHole](../edit/BSSN/UIUCBlackHole.py) module.


```python
# Step 3: Code Validation against BSSN.UIUCBlackHole NRPy+ module

import BSSN.UIUCBlackHole as uibh
uibh.UIUCBlackHole()

def compare(q1, q2, q1name, q2name):
    if sp.simplify(q1 - q2) != 0:
        print("Error: "+q1name+" - "+q2name+" = "+str(sp.simplify(q1 - q2)))
        sys.exit(1) 

print("Consistency check between UIUCBlackHole tutorial and NRPy+ BSSN.UIUCBlackHole module.")
compare(alpha, uibh.alpha, "alpha", "uibh.alpha")
for i in range(DIM):
    compare(betaU[i], uibh.betaU[i], "betaU"+str(i), "uibh.betaU"+str(i))
    compare(BU[i], uibh.BU[i], "BU"+str(i), "uibh.BU"+str(i))
    for j in range(DIM):
        compare(gammaDD[i][j], uibh.gammaDD[i][j], "gammaDD"+str(i)+str(j), "uibh.gammaDD"+str(i)+str(j))
        compare(KDD[i][j], uibh.KDD[i][j], "KDD"+str(i)+str(j), "uibh.KDD"+str(i)+str(j))
print("ALL TESTS PASS")
```

    Consistency check between UIUCBlackHole tutorial and NRPy+ BSSN.UIUCBlackHole module.
    ALL TESTS PASS


<a id='latex_pdf_output'></a>

# Step 7: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_Initial_Data-UIUC_BlackHole.pdf](Tutorial-ADM_Initial_Data-UIUC_BlackHole.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data-UIUC_BlackHole")
```

    Created Tutorial-ADM_Initial_Data-UIUC_BlackHole.tex, and compiled LaTeX
        file to PDF file Tutorial-ADM_Initial_Data-UIUC_BlackHole.pdf



```python

```
