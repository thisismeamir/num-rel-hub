<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [Shifted Kerr-Schild Solution](https://arxiv.org/pdf/1704.00599.pdf) Initial Data

## Authors: George Vopal & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook features a module that sets up Shifted Kerr-Schild initial data as per [Etienne et al., 2017 GiRaFFE](https://arxiv.org/pdf/1704.00599.pdf). The module confirms the expected order convergence to zero of the Hamiltonian and momentum constraint violations while utilizing Shifted Kerr-Schild coordinates to promote stability within the black hole horizon during the evolution of hydrodynamic, MHD, and FFE fields.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated to exhibit convergence to zero of the Hamiltonian and momentum constraint violations at the expected order to the exact solution (see plots at bottom of [the exact initial data validation start-to-finish tutorial notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb)).

### NRPy+ Source Code for this module: [BSSN/ShiftedKerrSchild.py](../edit/BSSN/ShiftedKerrSchild.py)

## Introduction:
Shifted Kerr-Schild coordinates are similar to the trumpet spacetime, in that $r=0$ maps to some finite radius surface in Kerr-Schild coordinates. The radial shift $r_0$ both reduces the black hole's coordinate size and causes the very strongly-curved spacetime fields at $r<r_{0}$ to vanish deep inside the horizon, which aids in numerical stability, e.g., when evolving hydrodynamic, MHD, and FFE fields inside the horizon.

<a id='toc'></a>

# Table of Contents:  
$$\label{toc}$$

1. [Step 1](#initialize_nrpy): Set up the needed NRPy+ infrastructure and declare core gridfunctions
1. [Step 2](#kerr_schild_lapse): The Kerr-Schild Lapse, Shift, and 3-Metric
    1. [Step 2.a](#define_rho): Define $\rho^{2}$, $\alpha$, $\beta^{r}$, $\beta^{\theta}$, $\beta^{\phi}$, $\gamma_{r\theta}$, $\gamma_{\theta\phi}$
    1. [Step 2.b](#nonzero_gamma): Define and construct nonzero components of $\gamma_{ij}$
1. [Step 3](#extrinsic_curvature): The extrinsic curvature $K_{ij}$
    1. [Step 3.a](#abc): Define useful quantities $A$, $B$, $C$
    1. [Step 3.b](#nonzero_k): Define and construct nonzero components of $K_{ij}$
1. [Step 4](#code_validation): Code Validation against `BSSN.ShiftedKerrSchild` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initialize_nrpy'></a>

# Step 1: Set up the needed NRPy+ infrastructure and declare core gridfunctions \[Back to [top](#toc)\]
$$\label{initialize_nrpy}$$

First, we will import the core modules from Python/NRPy+ and specify the main gridfunctions we will need.

**Input for initial data**:

* The black hole mass $M$.
* The black hole spin parameter $a$
* The radial offset $r_0$



```python
# Step P0: Load needed modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
r,th,ph = sp.symbols('r th ph', real=True)

thismodule = "ShiftedKerrSchild"

DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Input parameters:
M, a, r0 = par.Cparameters("REAL", thismodule,
                           ["M", "a", "r0"],
                           [1.0, 0.9,  1.0])

# Auxiliary variables:
rho2 = sp.symbols('rho2', real=True)
```

<a id='kerr_schild_lapse'></a>

# Step 2: The Kerr-Schild Lapse, Shift, and 3-Metric \[Back to [top](#toc)\]
$$\label{kerr_schild_lapse}$$

<a id='define_rho'></a>

## Step 2.a: Define $\rho^{2}$, $\alpha$, $\beta^{r_{\rm KS}}$, $\beta^{\theta}$, $\beta^{\phi}$, $\gamma_{r_{\rm KS}\theta}$, $\gamma_{\theta\phi}$ \[Back to [top](#toc)\]
$$\label{define_rho}$$

The relationship between the Kerr-Schild radius $r_{\rm KS}$ and the radial coordinate used on our numerical grid $r$, is given by

$$
r_{\rm KS} = r + r_0,
$$
where $r_0\ge 0$ is the radial shift.

Notice that the radial shift has no impact on Jacobians since $\frac{\partial{r_{\rm KS}}}{\partial{r}}=1$. $r_0$ must be set to a value less than the horizon radius $R$, but not so close to $R$ that finite-difference stencils from outside the horizon cross $r=0$. Thus $r_0$ must be set with consideration of the numerical grid structure in mind, as nonzero values of $r_0$ will shrink the coordinate size of the black hole by exactly $r_0$.

All of these equations are as defined in the appendix of the original GiRaFFE paper ([Etienne et al., 2017 GiRaFFE](https://arxiv.org/pdf/1704.00599.pdf)).
<br>
First, we define $\rho^{2}$ as

<br>

$$ \rho^2 = r_{\rm KS} + a^{2}\cos^{2}(\theta) $$

<br>

And we then define the Kerr-Schild lapse $\alpha$ from equation (A.1)

<br>

$$ \alpha = \frac{1}{\sqrt{1 + \frac{2Mr_{\rm KS}}{\rho^2}}} $$

<br>

And the shift $\beta$ from equations (A.2) & (A.3)

<br>

$$ \beta^{r_{\rm KS}} = \alpha^2\frac{2Mr_{\rm KS}}{\rho^2} $$

<br>

$$ \beta^{\theta} = \beta^{\phi} = \gamma_{r_{\rm KS}\theta} = \gamma_{\theta\phi} = 0 $$


```python
# Step 1: Define rho^2, alpha, beta^(r_{KS}), beta^(theta), beta^(phi), gamma_{r_{KS}theta}, gamma_{theta\phi}

# r_{KS} = r + r0
rKS = r+r0

# rho^2 = rKS^2 + a^2*cos^2(theta)
rho2 = rKS*rKS + a*a*sp.cos(th)**2

# alpha = 1/sqrt{1 + M*rKS/rho^2}
alpha = 1/(sp.sqrt(1 + 2*M*rKS/rho2))

# Initialize the shift vector, \beta^i, to zero.
betaU = ixp.zerorank1()
# beta^r = alpha^2*2Mr/rho^2
betaU[0] = alpha*alpha*2*M*rKS/rho2

# Time derivative of shift vector beta^i, B^i, is zero.
BU = ixp.zerorank1()
```

<a id='nonzero_gamma'></a>

## Step 2.b: Define and construct nonzero components $\gamma_{r_{\rm KS}r_{\rm KS}}$, $\gamma_{r_{\rm KS}\phi}$, $\gamma_{\theta\theta}$, $\gamma_{\phi\phi}$ \[Back to [top](#toc)\]
$$\label{nonzero_gamma}$$

From equations (A.4)-(A.7) of [Etienne et al., 2017](https://arxiv.org/pdf/1704.00599.pdf) we define the nonzero components of the 3-metric:

<br>

$$ \gamma_{r_{\rm KS}r_{\rm KS}} = 1 + \frac{2Mr_{\rm KS}}{\rho^2} $$

<br>

$$ \gamma_{r_{\rm KS}\phi} = -a\gamma_{r_{\rm KS}r_{\rm KS}}\sin^2(\theta) $$

<br>

$$ \gamma_{\theta\theta} = \rho^2 $$

<br>

$$ \gamma_{\phi\phi} = \left(r_{\rm KS}^2 + a^2 + \frac{2Mr_{\rm KS}}{\rho^2}a^{2}\sin^{2}(\theta)\right)\sin^{2}(\theta) $$


```python
# Step 2: Define and construct nonzero components gamma_{r_{KS}r_{KS}}$, gamma_{r_{KS}phi},
#         gamma_{thetatheta}, gamma_{phiphi}

# Initialize \gamma_{ij} to zero.
gammaDD = ixp.zerorank2()

# gammaDD{rKS rKS} = 1 +2M*rKS/rho^2
gammaDD[0][0] = 1 + 2*M*rKS/rho2

# gammaDD{rKS phi} = -a*gammaDD{r r}*sin^2(theta)
gammaDD[0][2] = gammaDD[2][0] = -a*gammaDD[0][0]*sp.sin(th)**2

# gammaDD{theta theta} = rho^2
gammaDD[1][1] = rho2

# gammaDD{phi phi} = (rKS^2 + a^2 + 2Mr/rho^2*a^2*sin^2(theta))*sin^2(theta)
gammaDD[2][2] = (rKS*rKS + a*a + 2*M*rKS*a*a*sp.sin(th)**2/rho2)*sp.sin(th)**2
```

<a id='extrinsic_curvature'></a>

# Step 3: The extrinsic curvature $K_{ij}$ \[Back to [top](#toc)\]
$$\label{extrinsic_curvature}$$

<a id='abc'></a>

## Step 3.a: Define useful quantities $A$, $B$, $C$ \[Back to [top](#toc)\]
$$\label{abc}$$

From equations (A.8)-(A.10) of [Etienne et al., 2017](https://arxiv.org/pdf/1704.00599.pdf) we define the following expressions which will help simplify the nonzero extrinsic curvature components:

<br>

$$ A = \left(a^{2}\cos(2\theta) + a^{2} + 2r_{\rm KS}^{2}\right) $$

<br>

$$ B = A + 4Mr_{\rm KS} $$

<br>

$$ D = \sqrt{\frac{2Mr_{\rm KS}}{a^{2}\cos^{2}(\theta) + r_{\rm KS}^2} + 1} $$



```python
# Step 3: Define useful quantities A, B, C

# A = (a^2*cos^2(2theta) + a^2 + 2r^2)
A = (a*a*sp.cos(2*th) + a*a + 2*rKS*rKS)

# B = A + 4M*rKS
B = A + 4*M*rKS

# D = \sqrt(2M*rKS/(a^2cos^2(theta) + rKS^2) + 1)
D = sp.sqrt(2*M*rKS/(a*a*sp.cos(th)**2 + rKS*rKS) + 1)
```

<a id='nonzero_k'></a>

## Step 3.b: Define and construct nonzero components of $K_{ij}$ \[Back to [top](#toc)\]
$$\label{nonzero_k}$$

We will now express the extrinsic curvature $K_{ij}$ in spherical polar coordinates.

From equations (A.11) - (A.13) of [Etienne et al., 2017](https://arxiv.org/pdf/1704.00599.pdf) we define the following:

$$ K_{r_{\rm KS}r_{\rm KS}} = \frac{D(A + 2Mr_{\rm KS})}{A^{2}B}\left[4M\left(a^{2}\cos(2\theta) + a^{2} - 2r_{\rm KS}^{2}\right)\right] $$

<br>

$$ K_{r_{\rm KS}\theta} = \frac{D}{AB}\left[8a^{2}Mr_{\rm KS}\sin(\theta)\cos(\theta)\right] $$

<br>

$$ K_{r_{\rm KS}\phi} = \frac{D}{A^2}\left[-2aM\sin^{2}(\theta)\left(a^{2}\cos(2\theta) + a^{2} - 2r_{\rm KS}^{2}\right)\right] $$


```python
# Step 4: Define the extrinsic curvature in spherical polar coordinates
# Establish the 3x3 zero-matrix
KDD = ixp.zerorank2()

# *** Fill in the nonzero components ***
# *** This will create an upper-triangular matrix ***
# K_{r r} = D(A+2Mr)/(A^2*B)[4M(a^2*cos(2theta) + a^2 - 2r^2)]
KDD[0][0] = D*(A+2*M*rKS)/(A*A*B)*(4*M*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))

# K_{r theta} = D/(AB)[8a^2*Mr*sin(theta)cos(theta)]
KDD[0][1] = KDD[1][0] = D/(A*B)*(8*a*a*M*rKS*sp.sin(th)*sp.cos(th))

# K_{r phi} = D/A^2[-2aMsin^2(theta)(a^2cos(2theta)+a^2-2r^2)]
KDD[0][2] = KDD[2][0] =  D/(A*A)*(-2*a*M*sp.sin(th)**2*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))
```

And from equations (A.14) - (A.17) of [Etienne et al., 2017](https://arxiv.org/pdf/1704.00599.pdf) we define the following expressions to complete the upper-triangular matrix $K_{ij}$:

$$ K_{\theta\theta} = \frac{D}{B}\left[4Mr_{\rm KS}^{2}\right] $$

<br>

$$ K_{\theta\phi} = \frac{D}{AB}\left[-8a^{3}Mr_{\rm KS}\sin^{3}(\theta)\cos(\theta)\right] $$

<br>

$$ K_{\phi\phi} = \frac{D}{A^{2}B}\left[2Mr_{\rm KS}\sin^{2}(\theta)\left(a^{4}(r_{\rm KS}-M)\cos(4\theta) + a^{4}(M + 3r_{\rm KS}) + 4a^{2}r_{\rm KS}^{2}(2r_{\rm KS} - M) + 4a^{2}r_{\rm KS}\cos(2\theta)\left(a^{2} + r_{\rm KS}(M + 2r_{\rm KS})\right) + 8r_{\rm KS}^{5}\right)\right] $$


```python
# K_{theta theta} = D/B[4Mr^2]
KDD[1][1] = D/B*(4*M*rKS*rKS)

# K_{theta phi} = D/(AB)*(-8*a^3*Mr*sin^3(theta)cos(theta))
KDD[1][2] = KDD[2][1] = D/(A*B)*(-8*a**3*M*rKS*sp.sin(th)**3*sp.cos(th))

# K_{phi phi} = D/(A^2*B)[2Mr*sin^2(theta)(a^4(M+3r)
#   +4a^2r^2(2r-M)+4a^2r*cos(2theta)(a^2+r(M+2r))+8r^5)]
KDD[2][2] = D/(A*A*B)*(2*M*rKS*sp.sin(th)**2*(a**4*(rKS-M)*sp.cos(4*th)\
                        + a**4*(M+3*rKS)+4*a*a*rKS*rKS*(2*rKS-M)\
                        + 4*a*a*rKS*sp.cos(2*th)*(a*a + rKS*(M + 2*rKS)) + 8*rKS**5))
```

<a id='code_validation'></a>

# Step 4: Code Validation against `BSSN.ShiftedKerrSchild` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for Shifted Kerr-Schild initial data between

1. this tutorial and 
2. the NRPy+ [BSSN.ShiftedKerrSchild](../edit/BSSN/ShiftedKerrSchild.py) module.


```python
# Step 3: Code Validation against BSSN.ShiftedKerrSchild NRPy+ module

import BSSN.ShiftedKerrSchild as sks
sks.ShiftedKerrSchild()

def compare(q1, q2, q1name, q2name):
    if sp.simplify(q1 - q2) != 0:
        print("Error: "+q1name+" - "+q2name+" = "+str(sp.simplify(q1 - q2)))
        sys.exit(1) 

print("Consistency check between ShiftedKerrSchild tutorial and NRPy+ BSSN.ShifedKerrSchild module.")
compare(alpha, sks.alpha, "alpha", "sks.alpha")
for i in range(DIM):
    compare(betaU[i], sks.betaU[i], "betaU"+str(i), "sks.betaU"+str(i))
    compare(BU[i], sks.BU[i], "BU"+str(i), "sks.BU"+str(i))
    for j in range(DIM):
        compare(gammaDD[i][j], sks.gammaDD[i][j], "gammaDD"+str(i)+str(j), "sks.gammaDD"+str(i)+str(j))
        compare(KDD[i][j], sks.KDD[i][j], "KDD"+str(i)+str(j), "sks.KDD"+str(i)+str(j))
print("ALL TESTS PASS")
```

    Consistency check between ShiftedKerrSchild tutorial and NRPy+ BSSN.ShifedKerrSchild module.
    ALL TESTS PASS


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_Initial_Data-ShiftedKerrSchild.pdf](Tutorial-ADM_Initial_Data-ShiftedKerrSchild.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data-ShiftedKerrSchild")
```

    Created Tutorial-ADM_Initial_Data-ShiftedKerrSchild.tex, and compiled LaTeX
        file to PDF file Tutorial-ADM_Initial_Data-ShiftedKerrSchild.pdf

