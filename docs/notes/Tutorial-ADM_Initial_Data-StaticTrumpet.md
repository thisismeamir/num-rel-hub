<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [Static Trumpet Black Hole](https://arxiv.org/abs/1403.5484) Initial data

## Authors: Terrence Pierre Jacques, Zach Etienne & Ian Ruchlin
### Formatting improvements courtesy Brandon Clark

## This notebook introduces a module establishing Static Trumpet Black Hole initial data, based on research by [Dennison and Baumgarte, 2014 Class. Quantum Grav. 31 117001](https://arxiv.org/abs/1403.5484). It allows conversion from Schwarzschild spacetime to any coordinate system detailed in [reference_metric.py](../edit/reference_metric.py), using the [Exact ADM Spherical-or-Cartesian-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_Exact_ADM_Spherical_or_Cartesian_to_BSSNCurvilinear.ipynb). Emphasizing isotropic coordinates, it validates against the original SENR code.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated to exhibit convergence to zero of the Hamiltonian constraint violation at the expected order to the exact solution (see plot at bottom of [the exact initial data validation start-to-finish tutorial notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb); momentum constraint is zero), and all quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).

### NRPy+ Source Code for this module: [BSSN/StaticTrumpet.py](../edit/BSSN/StaticTrumpet.py)

## Introduction:
These initial data are derived from a family of analytical coordinate systems representing the Schwarzschild spacetime. The coordinates extend smoothly through the black hole event horizon, the spatial coordinates are isotropic (so that the spatial metric can be written as a conformal factor to some power times a flat spatial metric), and, for almost all members of the family, the spatial slices take a so-called $\textit{trumpet geometry}$. Moreover, all expressions are surprisingly simple. This module sets the static trumpet black hole at the origin.

<a id='toc'></a>

# Table of Contents:  
$$\label{toc}$$

1. [Step 1](#initialize_nrpy): Set up the needed NRPy+ infrastructure and declare core gridfunctions
1. [Step 2](#conformal_factor_psi): The conformal factor $\psi$
    1. [Step 2.a](#define_psi): Define the conformal factor $\psi$
    1. [Step 2.b](#nonzero_gamma): Define and construct nonzero components of $\gamma_{ij}$
1. [Step 3](#extrinsic_curvature): Define and construct nonzero components of the extrinsic curvature $K_{ij}$, at the radius $R_0 = M$
1. [Step 4](#lapse_shift): Construct Lapse function $\alpha$ and components of shift vector $\beta$
1. [Step 5](#code_validation): Code Validation against `BSSN.StaticTrumpet` NRPy+ module
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initialize_nrpy'></a>

# Step 1: Set up the needed NRPy+ infrastructure and declare core gridfunctions \[Back to [top](#toc)\]
$$\label{initialize_nrpy}$$

First, we will import the core modules of Python/NRPy+ and specify the main gridfunctions that we will need. 

**Input for initial data**:

* The black hole mass $M$.



```python
# Step P0: Load needed modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
r,th,ph = sp.symbols('r th ph', real=True)

thismodule = "StaticTrumpet"

# Step 0: Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Step 1: Set psi, the conformal factor:
# Input parameters:
M = par.Cparameters("REAL", thismodule, ["M"], 1.0)
```

<a id='conformal_factor_psi'></a>

# Step 2: The conformal factor $\psi$ \[Back to [top](#toc)\]
$$\label{conformal_factor_psi}$$

<a id='define_psi'></a>

## Step 2.a: Define the conformal factor $\psi$ \[Back to [top](#toc)\]
$$\label{define_psi}$$


The conformal factor, defined in equation 13 of [Dennison and Baumgarte (2014)](https://arxiv.org/abs/1403.5484), setting $R_0 = M$,
$$ \psi = \sqrt{1 + \frac{M}{r}}. $$


```python
# Step 1: Set psi, the StaticTrumpet conformal factor
# Dennison and Baumgarte (2014) Eq. 13
# https://arxiv.org/pdf/1403.5484.pdf

# psi = sqrt{1 + M/r }
psi0 = sp.sqrt(1 + M/r)
```

<a id='nonzero_gamma'></a>

## Step 2.b: Define and construct nonzero components of $\gamma_{ij}$ \[Back to [top](#toc)\]
$$\label{nonzero_gamma}$$

The spatial metric, defined in equation 15 of [Dennison and Baumgarte (2014)](https://arxiv.org/abs/1403.5484),
$$ \gamma_{ij} = \psi^4 \eta_{ij}, $$

where $\eta_{ij}$ is the flat metric in spherical polar coordinates.


```python
# *** The physical spatial metric in spherical basis ***
# Set the upper-triangle of the matrix...
# Eq. 15
# eta_00 = 1, eta_11 = r^2, eta_22 = r^2 * sin^2 (theta)
gammaDD = ixp.zerorank2()
gammaDD[0][0] = psi0**4
gammaDD[1][1] = psi0**4 * r**2
gammaDD[2][2] = psi0**4 * r**2*sp.sin(th)**2
```

<a id='extrinsic_curvature'></a>

# Step 3: Define and construct nonzero components of the extrinsic curvature $K_{ij}$, at the radius $R_0 = M$ \[Back to [top](#toc)\]
$$\label{extrinsic_curvature}$$

Components of the extrinsic curvature in spherical basis, defined in equations 19 and 20 of [Dennison and Baumgarte (2014)](https://arxiv.org/abs/1403.5484),

$$ K_{rr} = - \frac{r \left( M-R_0 \right) + MR_0}{r^2 f_1}, $$

<br>

$$ K_{\theta\theta} = \frac{K_{\phi\phi}}{\sin^2 \theta} = f_1, $$

<br>

where $f_1 = \sqrt{2r \left( M-R_0 \right) + R_0 \left( 2M-R_0 \right).}$ Setting $R_0 = M$, these equations reduce to

<br>

$$ K_{rr} = -\frac{M}{r^2}, $$

<br>

$$ K_{\theta\theta} = \frac{K_{\phi\phi}}{\sin^2 \theta} = M. $$



```python
# *** The physical trace-free extrinsic curvature in spherical basis ***
# Set the upper-triangle of the matrix...

# Eq.19 and 20
KDD = ixp.zerorank2()

# K_{rr} = - M / r^2
KDD[0][0] = -M / r**2

# K_{theta theta} = K_{phi phi} / sin^2 theta = M
KDD[1][1] = M

KDD[2][2] = M * sp.sin(th)**2
```

<a id='lapse_shift'></a>

# Step 4: Construct Lapse function $\alpha$ and components of shift vector $\beta$\[Back to [top](#toc)\]
$$\label{lapse_shift}$$

Laspe function and shift vector components, equation 15 of [Dennison and Baumgarte (2014)](https://arxiv.org/abs/1403.5484), setting $R_0 = M$,
$$ \alpha = \frac{r}{r+M}, $$
<br>
$$ \beta^r = \frac{Mr}{\left (r+M \right)^2}, $$
<br>
$$\beta^\theta = \beta^\phi = 0. $$


```python
# Eq. 15
# alpha = r / (r+M)
alpha = r / (r + M)

betaU = ixp.zerorank1()
# beta^r = Mr / (r + M)^2
betaU[0] = M*r / (r + M)**2

BU = ixp.zerorank1()
```

<a id='code_validation'></a>

# Step 5: Code Validation against `BSSN.StaticTrumpet` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for static trumpet black hole initial data between

1. this tutorial and 
2. the NRPy+ [BSSN.StaticTrumpet](../edit/BSSN/StaticTrumpet.py) module.


```python
# Step 3: Code Validation against BSSN.StaticTrumpet NRPy+ module

import BSSN.StaticTrumpet as st
st.StaticTrumpet()

def compare(q1, q2, q1name, q2name):
    if sp.simplify(q1 - q2) != 0:
        print("Error: "+q1name+" - "+q2name+" = "+str(sp.simplify(q1 - q2)))
        sys.exit(1) 

print("Consistency check between StaticTrumpet tutorial and NRPy+ BSSN.StaticTrumpet module.")
compare(alpha, st.alpha, "alpha", "st.alpha")
for i in range(DIM):
    compare(betaU[i], st.betaU[i], "betaU"+str(i), "st.betaU"+str(i))
    compare(BU[i], st.BU[i], "BU"+str(i), "st.BU"+str(i))
    for j in range(DIM):
        compare(gammaDD[i][j], st.gammaDD[i][j], "gammaDD"+str(i)+str(j), "st.gammaDD"+str(i)+str(j))
        compare(KDD[i][j], st.KDD[i][j], "KDD"+str(i)+str(j), "st.KDD"+str(i)+str(j))
print("ALL TESTS PASS")
```

    Consistency check between StaticTrumpet tutorial and NRPy+ BSSN.StaticTrumpet module.
    ALL TESTS PASS


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_Initial_Data-StaticTrumpet.pdf](Tutorial-ADM_Initial_Data-StaticTrumpet.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data-StaticTrumpet")
```

    Created Tutorial-ADM_Initial_Data-StaticTrumpet.tex, and compiled LaTeX
        file to PDF file Tutorial-ADM_Initial_Data-StaticTrumpet.pdf

