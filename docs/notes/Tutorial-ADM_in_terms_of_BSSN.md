<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# ADM Quantities in terms of BSSN Quantities

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook demonstrates the conversion between the 4-metric $g_{\mu\nu}$ and the ADM or BSSN variables. Initially, ADM quantities are defined directly or in terms of BSSN quantities, followed by writing the 4-metric and its inverse using these quantities. The notebook then details the conversion of the ADM/BSSN metric quantities back into the 4-metric and validates the result against the NRPy+ BSSN.ADMBSSN_tofrom_4metric module.

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Source Code for this module: [ADM_in_terms_of_BSSN.py](../edit/BSSN/ADM_in_terms_of_BSSN.py)

## Introduction:
This tutorial notebook constructs all quantities in the [ADM formalism](https://en.wikipedia.org/wiki/ADM_formalism) (see also Chapter 2 in Baumgarte & Shapiro's book *Numerical Relativity*) in terms of quantities in our adopted (covariant, tensor-rescaled) BSSN formalism. That is to say, we will write the ADM quantities $\left\{\gamma_{ij},K_{ij},\alpha,\beta^i\right\}$ and their derivatives in terms of the BSSN quantities $\left\{\bar{\gamma}_{ij},\text{cf},\bar{A}_{ij},\text{tr}K,\alpha,\beta^i\right\}$ and their derivatives.

### A Note on Notation:

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates the temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

As a corollary, any expressions in NRPy+ involving mixed Greek and Latin indices will need to offset one set of indices by one; a Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial notebook).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#threemetric): The ADM three-metric $\gamma_{ij}$ and its derivatives in terms of rescaled BSSN quantities
    1. [Step 2.a](#derivatives_e4phi): Derivatives of $e^{4\phi}$
    1. [Step 2.b](#derivatives_adm_3metric): Derivatives of the ADM three-metric: $\gamma_{ij,k}$ and $\gamma_{ij,kl}$
    1. [Step 2.c](#christoffel): Christoffel symbols  $\Gamma^i_{jk}$  associated with the ADM 3-metric $\gamma_{ij}$
1. [Step 3](#extrinsiccurvature): The ADM extrinsic curvature $K_{ij}$ and its derivatives in terms of rescaled BSSN quantities
1. [Step 4](#code_validation): Code Validation against `BSSN.ADM_in_terms_of_BSSN` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from Python/NRPy+:


```python
# Step 1.a: Import all needed modules from NRPy+
import NRPy_param_funcs as par    # NRPy+: parameter interface
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python module for multiplatform OS-level functions

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

# Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
import BSSN.BSSN_quantities as Bq
Bq.BSSN_basic_tensors()
gammabarDD = Bq.gammabarDD
cf         = Bq.cf
AbarDD     = Bq.AbarDD
trK        = Bq.trK

Bq.gammabar__inverse_and_derivs()
gammabarDD_dD  = Bq.gammabarDD_dD
gammabarDD_dDD = Bq.gammabarDD_dDD

Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()
AbarDD_dD = Bq.AbarDD_dD
```

<a id='threemetric'></a>

# Step 2: The ADM three-metric $\gamma_{ij}$ and its derivatives in terms of rescaled BSSN quantities. \[Back to [top](#toc)\]
$$\label{threemetric}$$

The ADM three-metric is written in terms of the covariant BSSN three-metric tensor as (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
$$
\gamma_{ij} = \left(\frac{\gamma}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{i j},
$$
where $\gamma=\det{\gamma_{ij}}$ and $\bar{\gamma}=\det{\bar{\gamma}_{ij}}$. 

The "standard" BSSN conformal factor $\phi$ is given by (Eq. 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

\begin{align}
\phi &= \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right) \\
\implies e^{\phi} &= \left(\frac{\gamma}{\bar{\gamma}}\right)^{1/12} \\
\implies e^{4 \phi} &= \left(\frac{\gamma}{\bar{\gamma}}\right)^{1/3}
\end{align}

Thus the ADM three-metric may be written in terms of the BSSN three-metric and conformal factor $\phi$ as

$$
\gamma_{ij} = e^{4 \phi} \bar{\gamma}_{i j}.
$$

NRPy+'s implementation of BSSN allows for $\phi$ and two other alternative conformal factors to be defined:

\begin{align}
\chi &= e^{-4\phi} \\
W &= e^{-2\phi},
\end{align}

Thus if `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"chi"`, then

\begin{align}
\gamma_{ij} &= \frac{1}{\chi} \bar{\gamma}_{i j} \\
&= \frac{1}{\text{cf}} \bar{\gamma}_{i j},
\end{align}

and if `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"W"`, then
\begin{align}
\gamma_{ij} &= \frac{1}{W^2} \bar{\gamma}_{i j} \\
&= \frac{1}{\text{cf}^2} \bar{\gamma}_{i j}.
\end{align}


```python
# Step 2: The ADM three-metric gammaDD and its
#         derivatives in terms of BSSN quantities.
gammaDD = ixp.zerorank2()

exp4phi = sp.sympify(0)
if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
    exp4phi    = sp.exp(4*cf)
elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
    exp4phi    = (1 / cf)
elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
    exp4phi    = (1 / cf**2)
else:
    print("Error EvolvedConformalFactor_cf type = \""+par.parval_from_str("EvolvedConformalFactor_cf")+"\" unknown.")
    sys.exit(1)

for i in range(DIM):
    for j in range(DIM):
        gammaDD[i][j] = exp4phi*gammabarDD[i][j]
```

<a id='derivatives_e4phi'></a>

## Step 2.a: Derivatives of $e^{4\phi}$ \[Back to [top](#toc)\]
$$\label{derivatives_e4phi}$$

To compute derivatives of $\gamma_{ij}$ in terms of BSSN variables and their derivatives, we will first need derivatives of $e^{4\phi}$ in terms of the conformal BSSN variable `cf`.

\begin{align}
\frac{\partial}{\partial x^i} e^{4\phi} &= 4 e^{4\phi} \phi_{,i} \\
\implies \frac{\partial}{\partial x^j} \frac{\partial}{\partial x^i} e^{4\phi} &= \frac{\partial}{\partial x^j} \left(4 e^{4\phi} \phi_{,i}\right) \\
&= 16 e^{4\phi} \phi_{,i} \phi_{,j} + 4 e^{4\phi}  \phi_{,ij}
\end{align}

Thus computing first and second derivatives of $e^{4\phi}$ in terms of the BSSN quantity  `cf` requires only that we evaluate $\phi_{,i}$ and $\phi_{,ij}$ in terms of $e^{4\phi}$ (computed above in terms of `cf`) and derivatives of `cf`:

If `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"phi"`, then
\begin{align}
\phi_{,i} &= \text{cf}_{,i} \\
\phi_{,ij} &= \text{cf}_{,ij}
\end{align}

If `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"chi"`, then
\begin{align}
\text{cf} = e^{-4\phi} \implies \text{cf}_{,i} &= -4 e^{-4\phi} \phi_{,i} \\
\implies \phi_{,i} &= -\frac{e^{4\phi}}{4} \text{cf}_{,i} \\
\implies \phi_{,ij} &= -e^{4\phi} \phi_{,j} \text{cf}_{,i} -\frac{e^{4\phi}}{4} \text{cf}_{,ij}\\
&= -e^{4\phi} \left(-\frac{e^{4\phi}}{4} \text{cf}_{,j}\right) \text{cf}_{,i} -\frac{e^{4\phi}}{4} \text{cf}_{,ij} \\
&= \frac{1}{4} \left[\left(e^{4\phi}\right)^2 \text{cf}_{,i} \text{cf}_{,j} -e^{4\phi} \text{cf}_{,ij}\right] \\
\end{align}

If `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"W"`, then
\begin{align}
\text{cf} = e^{-2\phi} \implies \text{cf}_{,i} &= -2 e^{-2\phi} \phi_{,i} \\
\implies \phi_{,i} &= -\frac{e^{2\phi}}{2} \text{cf}_{,i} \\
\implies \phi_{,ij} &= -e^{2\phi} \phi_{,j} \text{cf}_{,i} -\frac{e^{2\phi}}{2} \text{cf}_{,ij}\\
&= -e^{2\phi} \left(-\frac{e^{2\phi}}{2} \text{cf}_{,j}\right) \text{cf}_{,i} -\frac{e^{2\phi}}{2} \text{cf}_{,ij} \\
&= \frac{1}{2} \left[e^{4\phi} \text{cf}_{,i} \text{cf}_{,j} -e^{2\phi} \text{cf}_{,ij}\right] \\
\end{align}


```python
# Step 2.a: Derivatives of $e^{4\phi}$
phidD = ixp.zerorank1()
phidDD = ixp.zerorank2()
cf_dD  = ixp.declarerank1("cf_dD")
cf_dDD = ixp.declarerank2("cf_dDD","sym01")
if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
    for i in range(DIM):
        phidD[i]  = cf_dD[i]
        for j in range(DIM):
            phidDD[i][j] = cf_dDD[i][j]
elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
    for i in range(DIM):
        phidD[i]  = -sp.Rational(1,4)*exp4phi*cf_dD[i]
        for j in range(DIM):
            phidDD[i][j] = sp.Rational(1,4)*( exp4phi**2*cf_dD[i]*cf_dD[j] - exp4phi*cf_dDD[i][j] )
elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
    exp2phi = (1 / cf)
    for i in range(DIM):
        phidD[i]  = -sp.Rational(1,2)*exp2phi*cf_dD[i]
        for j in range(DIM):
            phidDD[i][j] = sp.Rational(1,2)*( exp4phi*cf_dD[i]*cf_dD[j] - exp2phi*cf_dDD[i][j] )
else:
    print("Error EvolvedConformalFactor_cf type = \""+par.parval_from_str("EvolvedConformalFactor_cf")+"\" unknown.")
    sys.exit(1)

exp4phidD  = ixp.zerorank1()
exp4phidDD = ixp.zerorank2()
for i in range(DIM):
    exp4phidD[i] = 4*exp4phi*phidD[i]
    for j in range(DIM):
        exp4phidDD[i][j] = 16*exp4phi*phidD[i]*phidD[j] + 4*exp4phi*phidDD[i][j]
```

<a id='derivatives_adm_3metric'></a>

## Step 2.b: Derivatives of the ADM three-metric: $\gamma_{ij,k}$ and $\gamma_{ij,kl}$ \[Back to [top](#toc)\]
$$\label{derivatives_adm_3metric}$$

Recall the relation between the ADM three-metric $\gamma_{ij}$, the BSSN conformal three-metric $\bar{\gamma}_{i j}$, and the BSSN conformal factor $\phi$:

$$
\gamma_{ij} = e^{4 \phi} \bar{\gamma}_{i j}.
$$

Now that we have constructed derivatives of $e^{4 \phi}$ in terms of the chosen BSSN conformal factor `cf`, and the [BSSN.BSSN_quantities module](../edit/BSSN/BSSN_quantities.py) ([**tutorial**](Tutorial-BSSN_quantities.ipynb)) defines derivatives of $\bar{\gamma}_{ij}$ in terms of rescaled BSSN variables, derivatives of $\gamma_{ij}$ can be immediately constructed using the product rule:

\begin{align}
\gamma_{ij,k} &= \left(e^{4 \phi}\right)_{,k} \bar{\gamma}_{i j} + e^{4 \phi} \bar{\gamma}_{ij,k} \\
\gamma_{ij,kl} &= \left(e^{4 \phi}\right)_{,kl} \bar{\gamma}_{i j} + \left(e^{4 \phi}\right)_{,k} \bar{\gamma}_{i j,l} + \left(e^{4 \phi}\right)_{,l} \bar{\gamma}_{ij,k} + e^{4 \phi} \bar{\gamma}_{ij,kl}
\end{align}


```python
# Step 2.b: Derivatives of gammaDD, the ADM three-metric
gammaDDdD  = ixp.zerorank3()
gammaDDdDD = ixp.zerorank4()

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            gammaDDdD[i][j][k] = exp4phidD[k]*gammabarDD[i][j] + exp4phi*gammabarDD_dD[i][j][k]
            for l in range(DIM):
                gammaDDdDD[i][j][k][l] = exp4phidDD[k][l]*gammabarDD[i][j] + \
                                         exp4phidD[k]*gammabarDD_dD[i][j][l] + \
                                         exp4phidD[l]*gammabarDD_dD[i][j][k] + \
                                         exp4phi*gammabarDD_dDD[i][j][k][l]
```

<a id='christoffel'></a>

## Step 2.c: Christoffel symbols  $\Gamma^i_{jk}$  associated with the ADM 3-metric $\gamma_{ij}$ \[Back to [top](#toc)\]
$$\label{christoffel}$$

The 3-metric analog to the definition of Christoffel symbol (Eq. 1.18) in Baumgarte & Shapiro's *Numerical Relativity* is given by
$$
\Gamma^i_{jk} = \frac{1}{2} \gamma^{il} \left(\gamma_{lj,k} + \gamma_{lk,j} - \gamma_{jk,l} \right),
$$
which we implement here:


```python
# Step 2.c: 3-Christoffel symbols associated with ADM 3-metric gammaDD
# Step 2.c.i: First compute the inverse 3-metric gammaUU:
gammaUU, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

GammaUDD = ixp.zerorank3()

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                GammaUDD[i][j][k] += sp.Rational(1,2)*gammaUU[i][l]* \
                                (gammaDDdD[l][j][k] + gammaDDdD[l][k][j] - gammaDDdD[j][k][l])
```

<a id='extrinsiccurvature'></a>

# Step 3: The ADM extrinsic curvature $K_{ij}$ and its derivatives in terms of rescaled BSSN quantities. \[Back to [top](#toc)\]
$$\label{extrinsiccurvature}$$

The ADM extrinsic curvature may be written in terms of the BSSN trace-free extrinsic curvature tensor $\bar{A}_{ij}$ and the trace of the ADM extrinsic curvature $K$:

\begin{align}
K_{ij} &= \left(\frac{\gamma}{\bar{\gamma}}\right)^{1/3} \bar{A}_{ij} + \frac{1}{3} \gamma_{ij} K \\
&= e^{4\phi} \bar{A}_{ij} + \frac{1}{3} \gamma_{ij} K \\
\end{align}

We only compute first spatial derivatives of $K_{ij}$, as higher-derivatives are generally not needed:
$$
K_{ij,k} = \left(e^{4\phi}\right)_{,k} \bar{A}_{ij} + e^{4\phi} \bar{A}_{ij,k} + \frac{1}{3} \left(\gamma_{ij,k} K + \gamma_{ij} K_{,k}\right)
$$
which is expressed in terms of quantities already defined.


```python
# Step 3: Define ADM extrinsic curvature KDD and
#         its first spatial derivatives KDDdD
#         in terms of BSSN quantities
KDD = ixp.zerorank2()

for i in range(DIM):
    for j in range(DIM):
        KDD[i][j] = exp4phi*AbarDD[i][j] + sp.Rational(1,3)*gammaDD[i][j]*trK

KDDdD = ixp.zerorank3()
trK_dD = ixp.declarerank1("trK_dD")
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            KDDdD[i][j][k] = exp4phidD[k]*AbarDD[i][j] + exp4phi*AbarDD_dD[i][j][k] + \
                             sp.Rational(1,3)*(gammaDDdD[i][j][k]*trK + gammaDD[i][j]*trK_dD[k])
```

<a id='code_validation'></a>

# Step 4: Code Validation against `BSSN.ADM_in_terms_of_BSSN` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions between

1. this tutorial and 
2. the NRPy+ [BSSN.ADM_in_terms_of_BSSN](../edit/BSSN/ADM_in_terms_of_BSSN.py) module.



```python
def comp_func(expr1,expr2,basename,prefixname2="Bq."):
    if str(expr1-expr2)!="0":
        print(basename+" - "+prefixname2+basename+" = "+ str(expr1-expr2))
        return 1
    return 0

def gfnm(basename,idx1,idx2=None,idx3=None,idx4=None):
    if idx2 is None:
        return basename+"["+str(idx1)+"]"
    if idx3 is None:
        return basename+"["+str(idx1)+"]["+str(idx2)+"]"
    if idx4 is None:
        return basename+"["+str(idx1)+"]["+str(idx2)+"]["+str(idx3)+"]"
    return basename+"["+str(idx1)+"]["+str(idx2)+"]["+str(idx3)+"]["+str(idx4)+"]"

expr_list = []
exprcheck_list = []
namecheck_list = []

import BSSN.ADM_in_terms_of_BSSN as AB
AB.ADM_in_terms_of_BSSN()

namecheck_list.extend(["detgamma"])
exprcheck_list.extend([AB.detgamma])
expr_list.extend([detgamma])
for i in range(DIM):
    for j in range(DIM):
        namecheck_list.extend([gfnm("gammaDD",i,j),gfnm("gammaUU",i,j),gfnm("KDD",i,j)])
        exprcheck_list.extend([AB.gammaDD[i][j],AB.gammaUU[i][j],AB.KDD[i][j]])
        expr_list.extend([gammaDD[i][j],gammaUU[i][j],KDD[i][j]])
        for k in range(DIM):
            namecheck_list.extend([gfnm("gammaDDdD",i,j,k),gfnm("GammaUDD",i,j,k),gfnm("KDDdD",i,j,k)])
            exprcheck_list.extend([AB.gammaDDdD[i][j][k],AB.GammaUDD[i][j][k],AB.KDDdD[i][j][k]])
            expr_list.extend([gammaDDdD[i][j][k],GammaUDD[i][j][k],KDDdD[i][j][k]])
            for l in range(DIM):
                namecheck_list.extend([gfnm("gammaDDdDD",i,j,k,l)])
                exprcheck_list.extend([AB.gammaDDdDD[i][j][k][l]])
                expr_list.extend([gammaDDdDD[i][j][k][l]])

num_failures = 0
for i in range(len(expr_list)):
    num_failures += comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])


if num_failures == 0:
    print("ALL TESTS PASSED!")
else:
    print("ERROR. " + str(num_failures) + " TESTS FAILED")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_in_terms_of_BSSN.pdf](Tutorial-ADM_in_terms_of_BSSN.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_in_terms_of_BSSN")
```

    Created Tutorial-ADM_in_terms_of_BSSN.tex, and compiled LaTeX file to PDF
        file Tutorial-ADM_in_terms_of_BSSN.pdf

