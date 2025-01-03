<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Tetrads for Evaluating the Outgoing Gravitational Wave Weyl scalar $\psi_4$

## Authors: Patrick Nelson & Zach Etienne

## This tutorial demonstrates the construction of quasi-Kinnersley tetrads for the Weyl scalar $\Psi_4$. It details the process of initializing core NRPy+ modules, setting the numerical grid's coordinate system, defining vectors, and orthogonalizing them using the Levi-Civita symbol.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** See [$\psi_4$ notebook](Tutorial-Psi4.ipynb), whose Python module depends on the one presented here, for all validation notes.

### NRPy+ Source Code for this module: [BSSN/Psi4_tetrads.py](../edit/BSSN/Psi4_tetrads.py)

## Introduction: 
This module constructs tetrad vectors $l^\mu$, $m^\mu$, and $n^\mu$ for the $\psi_4$ Weyl scalar, a quantity that is immensely useful when extracting gravitational wave content from a numerical relativity simulation. $\psi_4$ is related to the gravitational wave strain via

$$
\psi_4 = \ddot{h}_+ - i \ddot{h}_\times.
$$

We construct $\psi_4$ from the standard ADM spatial metric $\gamma_{ij}$ and extrinsic curvature $K_{ij}$, and their derivatives. The full expression is given by Eq. 5.1 in [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf):

\begin{align}
\psi_4 &= \left[ {R}_{ijkl}+2K_{i[k}K_{l]j}\right]
{n}^i\bar{m}^j{n}^k\bar{m}^l  \\
& -8\left[ K_{j[k,l]}+{\Gamma }_{j[k}^pK_{l]p}\right]
{n}^{[0}\bar{m}^{j]}{n}^k\bar{m}^l \\
& +4\left[ {R}_{jl}-K_{jp}K_l^p+KK_{jl}\right]
{n}^{[0}\bar{m}^{j]}{n}^{[0}\bar{m}^{l]},
\end{align}

Note that $\psi_4$ is complex, with the imaginary components originating from the tetrad vector $m^\mu$. This module does not specify a tetrad; instead, it only constructs the above expression leaving $m^\mu$ and $n^\mu$ unspecified. This module defines these tetrad quantities, implementing the quasi-Kinnersley tetrad of [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf), also referred to as "***the BCL paper***".

### A Note on Notation:

As is standard in NRPy+, 

* Greek indices range from 0 to 3, inclusive, with the zeroth component denoting the temporal (time) component.
* Latin indices range from 0 to 2, inclusive, with the zeroth component denoting the first spatial component.

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This tutorial notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize needed NRPy+ modules
1. [Step 2](#quasikinnersley): The quasi-Kinnersley tetrad
1. [Step 3](#code_validation): Code Validation against `BSSN.Psi4_tetrads` NRPy+ module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1.a: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python modules for multiplatform OS-level functions

# Step 1.b: Set the coordinate system for the numerical grid
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")

# Step 1.c: Given the chosen coordinate system, set up
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
rfm.reference_metric()

# Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
#           a 3+1-dimensional decomposition of the general
#           relativistic field equations)
DIM = 3

# Step 1.e: Import all ADM quantities as written in terms of BSSN quantities
import BSSN.ADM_in_terms_of_BSSN as AB
AB.ADM_in_terms_of_BSSN()

# Step 1.f: Initialize TetradChoice parameter
thismodule = __name__
# Current option: QuasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)
par.initialize_param(par.glb_param("char", thismodule, "TetradChoice", "QuasiKinnersley"))
```

<a id='quasikinnersley'></a>

# Step 2: The quasi-Kinnersley tetrad of [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf) \[Back to [top](#toc)\]
$$\label{quasikinnersley}$$

To define the Weyl scalars, first, a tetrad must be chosen. Below, for compatibility with the [WeylScal4 diagnostic module](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src/master/WeylScal4/), we implement the quasi-Kinnersley tetrad of [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf).

We begin with the vectors given in eqs. 5.6 and 5.7 of the BCL paper, which are orthogonal to each other in flat spacetime; one is in the $\phi$ direction, one is in $r$, and the third is the cross product of the first two:
\begin{align}
    v_1^a &= [-y,x,0] \\
    v_2^a &= [x,y,z] \\
    v_3^a &= {\rm det}(\gamma)^{1/2} \gamma^{ad} \epsilon_{dbc} v_1^b v_2^c,
\end{align}

Notice that $v_1^a$ and $v_2^a$ assume the Cartesian basis, but $\gamma^{ad}$ will be in the $xx^i$ basis given by the chosen `reference_metric::CoordSystem`. Thus to construct $v_3^a$, we must first perform a change of basis on $v_1^a$ and $v_2^a$:

$$
v_{1,{\rm xx}}^a = \frac{\partial xx^a}{\partial x_{\rm Cart}^b} v_{1,{\rm Cart}}^b.
$$
This equation is problematic because we generally do not have a closed-form expression for components of the $xx^a$ vector as functions of the Cartesian coordinate vector components $x_{\rm Cart}^a$. However, we do have closed-form expressions for components of $x_{\rm Cart}^a$ as functions of $xx^a$. Thus we can construct the needed Jacobian matrix $\frac{\partial xx^a}{\partial x_{\rm Cart}^b}$ by evaluating the derivative $\frac{\partial x_{\rm Cart}^b}{\partial xx^a}$ and performing a simple matrix inversion:
$$
\frac{\partial xx^a}{\partial x_{\rm Cart}^b} = \left(\frac{\partial x_{\rm Cart}^b}{\partial xx^a} \right)^{-1}.
$$


```python
# Step 2.a: Declare the Cartesian x,y,z in terms of
#           xx0,xx1,xx2.
x = rfm.xx_to_Cart[0]
y = rfm.xx_to_Cart[1]
z = rfm.xx_to_Cart[2]

# Step 2.b: Declare detgamma and gammaUU from
#           BSSN.ADM_in_terms_of_BSSN;
#           simplify detgamma & gammaUU expressions,
#           which expedites Psi4 codegen.
detgamma = sp.simplify(AB.detgamma)
gammaUU = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        gammaUU[i][j] = sp.simplify(AB.gammaUU[i][j])

# Step 2.c: Define v1U and v2U
v1UCart = [-y, x, sp.sympify(0)]
v2UCart = [x, y, z]

# Step 2.d: Construct the Jacobian d x_Cart^i / d xx^j
Jac_dUCart_dDrfmUD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        Jac_dUCart_dDrfmUD[i][j] = sp.simplify(sp.diff(rfm.xx_to_Cart[i], rfm.xx[j]))

# Step 2.e: Invert above Jacobian to get needed d xx^j / d x_Cart^i
Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)

# Step 2.e.i: Simplify expressions for d xx^j / d x_Cart^i:
for i in range(DIM):
    for j in range(DIM):
        Jac_dUrfm_dDCartUD[i][j] = sp.simplify(Jac_dUrfm_dDCartUD[i][j])

# Step 2.f: Transform v1U and v2U from the Cartesian to the xx^i basis
v1U = ixp.zerorank1()
v2U = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        v1U[i] += Jac_dUrfm_dDCartUD[i][j] * v1UCart[j]
        v2U[i] += Jac_dUrfm_dDCartUD[i][j] * v2UCart[j]
```

... next we construct the third tetrad vector $v_3^a={\rm det}(\gamma)^{1/2} \gamma^{ad} \epsilon_{dbc} v_1^b v_2^c$, using the Levi-Civita symbol $\epsilon_{dbc}$ as defined in [indexedexp.py](../edit/indexedexp.py):


```python
# Step 2.g: Define v3U
v3U = ixp.zerorank1()
LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()
for a in range(DIM):
    for b in range(DIM):
        for c in range(DIM):
            for d in range(DIM):
                v3U[a] += sp.sqrt(detgamma)*gammaUU[a][d]*LeviCivitaSymbolDDD[d][b][c]*v1U[b]*v2U[c]

# Step 2.g.i: Simplify expressions for v1U,v2U,v3U. This greatly expedites the C code generation (~10x faster)
#             Drat. Simplification with certain versions of SymPy & coord systems results in a hang. Let's just
#             evaluate the expressions so the most trivial optimizations can be performed.
for a in range(DIM):
    v1U[a] = v1U[a].doit() #sp.simplify(v1U[a])
    v2U[a] = v2U[a].doit() #sp.simplify(v2U[a])
    v3U[a] = v3U[a].doit() #sp.simplify(v3U[a])
```

As our next step, we carry out the Gram-Schmidt orthonormalization process. The vectors $v_i^a$ are placeholders in the code; the final product of the orthonormalization is the vectors $e_i^a$. So,
\begin{align}
e_1^a &= \frac{v_1^a}{\sqrt{\omega_{11}}} \\
e_2^a &= \frac{v_2^a - \omega_{12} e_1^a}{\sqrt{\omega_{22}}} \\
e_3^a &= \frac{v_3^a - \omega_{13} e_1^a - \omega_{23} e_2^a}{\sqrt{\omega_{33}}}, \text{ where}\\
\omega_{ij} &= v_i^a v_j^b \gamma_{ab}
\end{align}

Note that the above expressions must be evaluated with the numerators first so that the denominators generate the proper normalization.


```python
# Step 2.h: Define omega_{ij}
omegaDD = ixp.zerorank2()
gammaDD = AB.gammaDD
def v_vectorDU(v1U,v2U,v3U,  i,a):
    if i==0:
        return v1U[a]
    if i==1:
        return v2U[a]
    if i==2:
        return v3U[a]
    print("ERROR: unknown vector!")
    sys.exit(1)

def update_omega(omegaDD, i,j, v1U,v2U,v3U,gammaDD):
    omegaDD[i][j] = sp.sympify(0)
    for a in range(DIM):
        for b in range(DIM):
            omegaDD[i][j] += v_vectorDU(v1U,v2U,v3U, i,a)*v_vectorDU(v1U,v2U,v3U, j,b)*gammaDD[a][b]

# Step 2.i: Define e^a_i. Note that:
#           omegaDD[0][0] = \omega_{11} above;
#           omegaDD[1][1] = \omega_{22} above, etc.
e1U = ixp.zerorank1()
e2U = ixp.zerorank1()
e3U = ixp.zerorank1()
# First e_1^a: Orthogonalize & normalize:
update_omega(omegaDD, 0,0, v1U,v2U,v3U,gammaDD)
for a in range(DIM):
    e1U[a] = v1U[a]/sp.sqrt(omegaDD[0][0])

# Next e_2^a: First orthogonalize:
update_omega(omegaDD, 0,1, e1U,v2U,v3U,gammaDD)
for a in range(DIM):
    e2U[a] = (v2U[a] - omegaDD[0][1]*e1U[a])
# Then normalize:
update_omega(omegaDD, 1,1, e1U,e2U,v3U,gammaDD)
for a in range(DIM):
    e2U[a] /= sp.sqrt(omegaDD[1][1])

# Next e_3^a: First orthogonalize:
update_omega(omegaDD, 0,2, e1U,e2U,v3U,gammaDD)
update_omega(omegaDD, 1,2, e1U,e2U,v3U,gammaDD)
for a in range(DIM):
    e3U[a] = (v3U[a] - omegaDD[0][2]*e1U[a] - omegaDD[1][2]*e2U[a])
# Then normalize:
update_omega(omegaDD, 2,2, e1U,e2U,e3U,gammaDD)
for a in range(DIM):
    e3U[a] /= sp.sqrt(omegaDD[2][2])
```

Once we have orthogonal, normalized vectors, we can construct the tetrad itself, again drawing on eqs. 5.6. We can draw on SymPy's built-in tools for complex numbers to build the complex vector $m^a$:
\begin{align}
    l^\mu &= \frac{1}{\sqrt{2}} \left(u^\mu + r^\mu\right) \\
    n^\mu &= \frac{1}{\sqrt{2}} \left(u^\mu - r^\mu\right) \\
    \Re(m^\mu) &= \frac{1}{\sqrt{2}} \theta^\mu \\
    \Im(m^\mu) &= \frac{1}{\sqrt{2}} \phi^\mu,
\end{align}
where $r^\mu=\{0,e_2^i\}$, $\theta^\mu=\{0,e_3^i\}$, $\phi^\mu=\{0,e_1^i\}$, and $u^\mu$ is the time-like unit normal to the hypersurface.


```python
# Step 2.j: Construct l^mu, n^mu, and m^mu, based on r^mu, theta^mu, phi^mu, and u^mu:
r4U     = ixp.zerorank1(DIM=4)
u4U     = ixp.zerorank1(DIM=4)
theta4U = ixp.zerorank1(DIM=4)
phi4U   = ixp.zerorank1(DIM=4)

for a in range(DIM):
    r4U[    a+1] = e2U[a]
    theta4U[a+1] = e3U[a]
    phi4U[  a+1] = e1U[a]

# FIXME? assumes alpha=1, beta^i = 0
u4U[0] = 1

l4U = ixp.zerorank1(DIM=4)
n4U = ixp.zerorank1(DIM=4)
mre4U  = ixp.zerorank1(DIM=4)
mim4U  = ixp.zerorank1(DIM=4)

# M_SQRT1_2 = 1 / sqrt(2) (defined in math.h on Linux)
M_SQRT1_2 = par.Cparameters("#define",thismodule,"M_SQRT1_2","")
isqrt2 = M_SQRT1_2 #1/sp.sqrt(2) <- SymPy drops precision to 15 sig. digits in unit tests
for mu in range(4):
    l4U[mu]   = isqrt2*(u4U[mu] + r4U[mu])
    n4U[mu]   = isqrt2*(u4U[mu] - r4U[mu])
    mre4U[mu] = isqrt2*theta4U[mu]
    mim4U[mu] = isqrt2*  phi4U[mu]

# ltetU,ntetU,remtetU,immtetU,e1U,e2U,e3U
for mu in range(4):
    l4U[mu]   = isqrt2*(u4U[mu] + r4U[mu])
    n4U[mu]   = isqrt2*(u4U[mu] - r4U[mu])
    mre4U[mu] = isqrt2*theta4U[mu]
    mim4U[mu] = isqrt2*  phi4U[mu]
```

<a id='code_validation'></a>

# Step 3: Code validation against `BSSN.Psi4_tetrads` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.Psi4_tetrads](../edit/BSSN/Psi4_tetrads.py) module.

By default, we compare all quantities in Spherical coordinates, though other coordinate systems may be chosen.


```python
def comp_func(expr1,expr2,basename,prefixname2="BP4T."):
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

import BSSN.Psi4_tetrads as BP4T
BP4T.Psi4_tetrads()

for mu in range(4):
    namecheck_list.extend([gfnm("l4U",mu),gfnm("n4U",mu),gfnm("mre4U",mu),gfnm("mim4U",mu)])
    exprcheck_list.extend([BP4T.l4U[mu],BP4T.n4U[mu],BP4T.mre4U[mu],BP4T.mim4U[mu]])
    expr_list.extend([l4U[mu],n4U[mu],mre4U[mu],mim4U[mu]])

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

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Psi4_tetrads.pdf](Tutorial-Psi4_tetrads.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Psi4_tetrads")
```

    Created Tutorial-Psi4_tetrads.tex, and compiled LaTeX file to PDF file
        Tutorial-Psi4_tetrads.pdf

