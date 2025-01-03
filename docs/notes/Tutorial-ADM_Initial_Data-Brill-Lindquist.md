<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [Brill-Lindquist](https://journals.aps.org/pr/abstract/10.1103/PhysRev.131.471) Initial data
## Author: Zach Etienne
###  Formatting improvements courtesy Brandon Clark

## This notebook presents Brill-Lindquist initial data for merging black holes in Cartesian coordinates, using Python/NRPy+ modules. It establishes the conformal factor, $\psi$, defines ADM variables in Cartesian coordinates, and verifies the setup against NRPy+'s BrillLindquist module. 

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated to exhibit convergence to zero of the Hamiltonian constraint violation at the expected order to the exact solution (see plot at bottom of [the exact initial data validation start-to-finish tutorial notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb); momentum constraint is zero), and all quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).

### NRPy+ Source Code for this module: [BrillLindquist.py](../edit/BSSN/BrillLindquist.py)

## Introduction:
This module sets up initial data for a merging black hole system in Cartesian coordinates.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#initialdata): Setting up Brill-Lindquist initial data
1. [Step 3](#code_validation): Code Validation against BSSN/BrillLindquist NRPy+ module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from Python/NRPy+:


```python
# Step 1: Initialize core Python/NRPy+ modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
```

<a id='initialdata'></a>

# Step 2: Setting up Brill-Lindquist initial data \[Back to [top](#toc)\]
$$\label{initialdata}$$ 

Here we set up Brill-Lindquist initial data ([Brill & Lindquist, Phys. Rev. 131, 471, 1963](https://journals.aps.org/pr/abstract/10.1103/PhysRev.131.471); see also Eq. 1 of [Brandt & Brügmann, arXiv:gr-qc/9711015v1](https://arxiv.org/pdf/gr-qc/9711015v1.pdf)):

The conformal factor $\psi$ for Brill-Lindquist initial data is given by
$$\psi = e^{\phi} = 1 + \sum_{i=1}^N \frac{m_{(i)}}{2 \left|\vec{r}_{(i)} - \vec{r}\right|};\quad K_{ij}=0,$$

where $\psi$ is written in terms of the 3-metric $\gamma_{ij}$ as

$$
\gamma_{ij} = \psi^4 \delta_{ij}.
$$

The extrinsic curvature is zero:
$$
K_{ij} = 0
$$

These data consist of $N$ nonspinning black holes initially at rest. This module restricts to the case of two such black holes, positioned anywhere. Here, we implement $N=2$.

**Inputs for $\psi$**:
* The position and (bare) mass of black hole 1: $\left(x_{(1)},y_{(1)},z_{(1)}\right)$ and $m_{(1)}$, respectively
* The position and (bare) mass of black hole 2: $\left(x_{(2)},y_{(2)},z_{(2)}\right)$ and $m_{(2)}$, respectively

**Additional variables needed for spacetime evolution**:
* Desired coordinate system
* Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=\psi^{-2}$ and $\beta^i=B^i=0$.


```python
# Step 2: Setting up Brill-Lindquist initial data

thismodule = "Brill-Lindquist"
BH1_posn_x,BH1_posn_y,BH1_posn_z = par.Cparameters("REAL", thismodule,
                                                   ["BH1_posn_x","BH1_posn_y","BH1_posn_z"],
                                                   [         0.0,         0.0,        +0.5])
BH1_mass = par.Cparameters("REAL", thismodule, ["BH1_mass"],1.0)
BH2_posn_x,BH2_posn_y,BH2_posn_z = par.Cparameters("REAL", thismodule,
                                                   ["BH2_posn_x","BH2_posn_y","BH2_posn_z"],
                                                   [         0.0,         0.0,        -0.5])
BH2_mass = par.Cparameters("REAL", thismodule, ["BH2_mass"],1.0)

# Step 2.a: Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

Cartxyz = ixp.declarerank1("Cartxyz")

# Step 2.b: Set psi, the conformal factor:
psi = sp.sympify(1)
psi += BH1_mass / ( 2 * sp.sqrt((Cartxyz[0]-BH1_posn_x)**2 + (Cartxyz[1]-BH1_posn_y)**2 + (Cartxyz[2]-BH1_posn_z)**2) )
psi += BH2_mass / ( 2 * sp.sqrt((Cartxyz[0]-BH2_posn_x)**2 + (Cartxyz[1]-BH2_posn_y)**2 + (Cartxyz[2]-BH2_posn_z)**2) )

# Step 2.c: Set all needed ADM variables in Cartesian coordinates
gammaDD = ixp.zerorank2()
KDD     = ixp.zerorank2() # K_{ij} = 0 for these initial data
for i in range(DIM):
    gammaDD[i][i] = psi**4

alpha = 1/psi**2
betaU = ixp.zerorank1() # We generally choose \beta^i = 0 for these initial data
BU    = ixp.zerorank1() # We generally choose B^i = 0 for these initial data
```

<a id='code_validation'></a>

# Step 3: Code Validation against BSSN/BrillLindquist NRPy+ module  \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for Brill-Lindquist initial data between

1. this tutorial and 
2. the NRPy+ [BSSN/BrillLindquist.py](../edit/BSSN/BrillLindquist.py) module.


```python
# Step 3: Code Validation against BSSN.BrillLindquist NRPy+ module

import BSSN.BrillLindquist as bl
bl.BrillLindquist()

def compare(q1, q2, q1name, q2name):
    if sp.simplify(q1 - q2) != 0:
        print("Error: "+q1name+" - "+q2name+" = "+str(sp.simplify(q1 - q2)))
        sys.exit(1) 

print("Consistency check between Brill-Lindquist tutorial and NRPy+ BSSN.BrillLindquist module.")
compare(alpha, bl.alpha, "alpha", "bl.alpha")
for i in range(DIM):
    compare(betaU[i], bl.betaU[i], "betaU"+str(i), "bl.betaU"+str(i))
    compare(BU[i], bl.BU[i], "BU"+str(i), "bl.BU"+str(i))
    for j in range(DIM):
        compare(gammaDD[i][j], bl.gammaDD[i][j], "gammaDD"+str(i)+str(j), "bl.gammaDD"+str(i)+str(j))
        compare(KDD[i][j], bl.KDD[i][j], "KDD"+str(i)+str(j), "bl.KDD"+str(i)+str(j))
print("ALL TESTS PASS")
```

    Consistency check between Brill-Lindquist tutorial and NRPy+ BSSN.BrillLindquist module.
    ALL TESTS PASS


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_Initial_Data-Brill-Lindquist](Tutorial-ADM_Initial_Data-Brill-Lindquist.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)



```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data-Brill-Lindquist")
```

    Created Tutorial-ADM_Initial_Data-Brill-Lindquist.tex, and compiled LaTeX
        file to PDF file Tutorial-ADM_Initial_Data-Brill-Lindquist.pdf

