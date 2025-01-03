<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# BSSN Quantities in terms of ADM Quantities
## Author: Zach Etienne

## This notebook showcases the conversion of ADM to BSSN variables using NRPy+, preparing them for solving Einstein's equations with the BSSN formulation. It covers the rescaling of 3-metric, extrinsic curvature, and gauge quantities, with an outline of the rescaling processes from covariant BSSN formulation and the underlying references.

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Source Code for this module: [BSSN_in_terms_of_ADM.py](../edit/BSSN/BSSN_in_terms_of_ADM.py)

## Introduction:
This module documents the conversion of ADM variables:

$$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i\right\}$$

into BSSN variables

$$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$ 

in the desired curvilinear basis (given by `reference_metric::CoordSystem`). Then it rescales the resulting BSSNCurvilinear variables (as defined in [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb)) into the form needed for solving Einstein's equations with the BSSN formulation:

$$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}.$$

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules; set desired output BSSN Curvilinear coordinate system set to Spherical
1. [Step 2](#adm2bssn): Perform the ADM-to-BSSN conversion for 3-metric, extrinsic curvature, and gauge quantities
    1. [Step 2.a](#adm2bssn_gamma): Convert ADM $\gamma_{ij}$ to BSSN $\bar{\gamma}_{ij}$; rescale to get $h_{ij}$
    1. [Step 2.b](#admexcurv_convert): Convert the ADM extrinsic curvature $K_{ij}$ to BSSN $\bar{A}_{ij}$ and $K$; rescale to get $a_{ij}$, $K$.
    1. [Step 2.c](#lambda): Define $\bar{\Lambda}^i$
    1. [Step 2.d](#conformal): Define the conformal factor variable `cf`
1. [Step 3](#code_validation): Code Validation against `BSSN.BSSN_in_terms_of_ADM` NRPy+ module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$



```python
# Step 1: Import needed core NRPy+ modules
import sympy as sp                 # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par     # NRPy+: Parameter interface
import indexedexp as ixp           # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm     # NRPy+: Reference metric support
import sys                         # Standard Python modules for multiplatform OS-level functions
import BSSN.BSSN_quantities as Bq  # NRPy+: This module depends on the parameter EvolvedConformalFactor_cf,
                                   #        which is defined in BSSN.BSSN_quantities

# Step 1.a: Set DIM=3, as we're using a 3+1 decomposition of Einstein's equations
DIM=3
```

<a id='adm2bssn'></a>

# Step 2: Perform the ADM-to-BSSN conversion for 3-metric, extrinsic curvature, and gauge quantities \[Back to [top](#toc)\]
$$\label{adm2bssn}$$

Here we convert ADM quantities to their BSSN Curvilinear counterparts.

<a id='adm2bssn_gamma'></a>

## Step 2.a: Convert ADM $\gamma_{ij}$ to BSSN $\bar{\gamma}_{ij}$; rescale to get $h_{ij}$ \[Back to [top](#toc)\]
$$\label{adm2bssn_gamma}$$

We have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
$$
\bar{\gamma}_{i j} = \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \gamma_{ij},
$$
where we always make the choice $\bar{\gamma} = \hat{\gamma}$.

After constructing $\bar{\gamma}_{ij}$, we rescale to get $h_{ij}$ according to the prescription described in the [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

$$
h_{ij} = (\bar{\gamma}_{ij} - \hat{\gamma}_{ij})/\text{ReDD[i][j]}.
$$


```python
# Step 2: All ADM quantities were input into this function in the Spherical or Cartesian
#         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
#         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
#         Here we convert ADM quantities to their BSSN Curvilinear counterparts:

# Step 2.a: Convert ADM $\gamma_{ij}$ to BSSN $\bar{gamma}_{ij}$:
#           We have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
def gammabarDD_hDD(gammaDD):
    global gammabarDD,hDD
    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.hDD_given_ADM(): Must call reference_metric() first!")
        sys.exit(1)
    # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
    _gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    gammabarDD = ixp.zerorank2()
    hDD        = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*gammaDD[i][j]
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
```

<a id='admexcurv_convert'></a>

## Step 2.b: Convert the ADM extrinsic curvature $K_{ij}$ to BSSN quantities $\bar{A}_{ij}$ and $K={\rm tr}(K_{ij})$; rescale $\bar{A}_{ij}$ to get $a_{ij}$ \[Back to [top](#toc)\]
$$\label{admexcurv_convert}$$

Convert the ADM extrinsic curvature $K_{ij}$ to the trace-free extrinsic curvature $\bar{A}_{ij}$, plus the trace of the extrinsic curvature $K$, where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
\begin{align}
K &= \gamma^{ij} K_{ij} \\
\bar{A}_{ij} &= \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \left(K_{ij} - \frac{1}{3} \gamma_{ij} K \right)
\end{align}

After constructing $\bar{A}_{ij}$, we rescale to get $a_{ij}$ according to the prescription described in the [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

$$
a_{ij} = \bar{A}_{ij}/\text{ReDD[i][j]}.
$$


```python
# Step 2.b: Convert the extrinsic curvature K_{ij} to the trace-free extrinsic
#           curvature \bar{A}_{ij}, plus the trace of the extrinsic curvature K,
#           where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
def trK_AbarDD_aDD(gammaDD,KDD):
    global trK,AbarDD,aDD
    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.trK_AbarDD(): Must call reference_metric() first!")
        sys.exit(1)
    # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    # K = gamma^{ij} K_{ij}, and
    # \bar{A}_{ij} &= (\frac{\bar{gamma}}{gamma})^{1/3}*(K_{ij} - \frac{1}{3}*gamma_{ij}*K)
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j]*KDD[i][j]

    AbarDD = ixp.zerorank2()
    aDD    = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AbarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*(KDD[i][j] - sp.Rational(1,3)*gammaDD[i][j]*trK)
            aDD[i][j]    = AbarDD[i][j] / rfm.ReDD[i][j]
```

<a id='lambda'></a>

## Step 2.c: Assuming the ADM 3-metric $\gamma_{ij}$ is given as an explicit function of `(xx0,xx1,xx2)`, convert to BSSN $\bar{\Lambda}^i$; rescale to compute $\lambda^i$ \[Back to [top](#toc)\]
$$\label{lambda}$$

To define $\bar{\Lambda}^i$ we implement Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf):
$$
\bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).
$$

The [reference_metric.py](../edit/reference_metric.py) module provides us with exact, closed-form expressions for $\hat{\Gamma}^i_{jk}$, so here we need only compute exact expressions for $\bar{\Gamma}^i_{jk}$, based on $\gamma_{ij}$ given as an explicit function of `(xx0,xx1,xx2)`. This is particularly useful when setting up initial data.

After constructing $\bar{\Lambda}^i$, we rescale to get $\lambda^i$ according to the prescription described in the [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

$$
\lambda^i = \bar{\Lambda}^i/\text{ReU[i]}.
$$


```python
# Step 2.c: Define \bar{Lambda}^i (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
def LambdabarU_lambdaU__exact_gammaDD(gammaDD):
    global LambdabarU,lambdaU

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    gammabarUU, _gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

    # First compute Christoffel symbols \bar{Gamma}^i_{jk}, with respect to barred metric:
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1,2)*gammabarUU[i][l]*( sp.diff(gammabarDD[l][j],rfm.xx[k]) +
                                                                                sp.diff(gammabarDD[l][k],rfm.xx[j]) -
                                                                                sp.diff(gammabarDD[j][k],rfm.xx[l]) )
    # Next evaluate \bar{Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])
    for i in range(DIM):
        # We evaluate LambdabarU[i] here to ensure proper cancellations. If these cancellations
        #   are not applied, certain expressions (e.g., lambdaU[0] in StaticTrumpet) will
        #   cause SymPy's (v1.5+) CSE algorithm to hang
        LambdabarU[i] = LambdabarU[i].doit()
    lambdaU    = ixp.zerorank1()
    for i in range(DIM):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]
```

<a id='conformal'></a>

## Step 2.d: Define the conformal factor variable `cf` \[Back to [top](#toc)\]
$$\label{conformal}$$

We define the conformal factor variable `cf` based on the setting of the `"BSSN_quantities::EvolvedConformalFactor_cf"` parameter.

For example if `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"phi"`, we can use Eq. 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf), which in arbitrary coordinates is written:

$$
\phi = \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right).
$$

Alternatively if `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"chi"`, then
$$
\chi = e^{-4 \phi} = \exp\left(-4 \frac{1}{12} \left(\frac{\gamma}{\bar{\gamma}}\right)\right) 
= \exp\left(-\frac{1}{3} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/3}.
$$

Finally if `"BSSN_quantities::EvolvedConformalFactor_cf"` is set to `"W"`, then
$$
W = e^{-2 \phi} = \exp\left(-2 \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = 
\exp\left(-\frac{1}{6} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = 
\left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/6}.
$$


```python
# Step 2.d: Set the conformal factor variable cf, which is set
#           by the "BSSN_quantities::EvolvedConformalFactor_cf" parameter. For example if
#           "EvolvedConformalFactor_cf" is set to "phi", we can use Eq. 3 of
#           [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf),
#           which in arbitrary coordinates is written:
def cf_from_gammaDD(gammaDD):
    global cf

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    _gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)
    _gammaUU, gammaDET       = ixp.symm_matrix_inverter3x3(gammaDD)

    cf = sp.sympify(0)

    if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
        # phi = \frac{1}{12} log(\frac{gamma}{\bar{gamma}}).
        cf = sp.Rational(1,12)*sp.log(gammaDET/gammabarDET)
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
        # chi = exp(-4*phi) = exp(-4*\frac{1}{12}*(\frac{gamma}{\bar{gamma}}))
        #      = exp(-\frac{1}{3}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{\bar{gamma}})^{-1/3}.
        #
        cf = (gammaDET/gammabarDET)**(-sp.Rational(1,3))
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
        # W = exp(-2*phi) = exp(-2*\frac{1}{12}*log(\frac{gamma}{\bar{gamma}}))
        #   = exp(-\frac{1}{6}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{bar{gamma}})^{-1/6}.
        cf = (gammaDET/gammabarDET)**(-sp.Rational(1,6))
    else:
        print("Error EvolvedConformalFactor_cf type = \""+par.parval_from_str("EvolvedConformalFactor_cf")+"\" unknown.")
        sys.exit(1)
```

<a id='betvet'></a>

## Step 2.e: Rescale $\beta^i$ and $B^i$ to compute $\mathcal{V}^i={\rm vet}^i$ and $\mathcal{B}^i={\rm bet}^i$, respectively \[Back to [top](#toc)\]
$$\label{betvet}$$

We rescale $\beta^i$ and $B^i$ according to the prescription described in the [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
\begin{align}
\mathcal{V}^i &= \beta^i/\text{ReU[i]}\\
\mathcal{B}^i &= B^i/\text{ReU[i]}.
\end{align}


```python
# Step 2.e: Rescale beta^i and B^i according to the prescription described in
#         the [BSSN in curvilinear coordinates tutorial notebook](Tutorial-BSSNCurvilinear.ipynb)
#         (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
#
# \mathcal{V}^i &= beta^i/(ReU[i])
# \mathcal{B}^i &= B^i/(ReU[i])
def betU_vetU(betaU,BU):
    global vetU,betU

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.bet_vet(): Must call reference_metric() first!")
        sys.exit(1)
    vetU    = ixp.zerorank1()
    betU    = ixp.zerorank1()
    for i in range(DIM):
        vetU[i]    =      betaU[i] / rfm.ReU[i]
        betU[i]    =         BU[i] / rfm.ReU[i]
```

<a id='code_validation'></a>

# Step 3: Code Validation against `BSSN.BSSN_in_terms_of_ADM` module \[Back to [top](#toc)\] 
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for [UIUC initial data](Tutorial-ADM_Initial_Data-UIUC_BlackHole.ipynb) between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_in_terms_of_ADM](../edit/BSSN/BSSN_in_terms_of_ADM.py) module.

As no basis transformation is performed, we analyze these expressions in their native, Spherical coordinates.


```python
# Step 3.a: Set the desired *output* coordinate system to Spherical:
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

# Step 3.b: Set up initial data; assume UIUC spinning black hole initial data
import BSSN.UIUCBlackHole as uibh
uibh.UIUCBlackHole()

# Step 3.c: Call above functions to convert ADM to BSSN curvilinear
gammabarDD_hDD(                   uibh.gammaDD)
trK_AbarDD_aDD(                   uibh.gammaDD,uibh.KDD)
LambdabarU_lambdaU__exact_gammaDD(uibh.gammaDD)
cf_from_gammaDD(                  uibh.gammaDD)
betU_vetU(                        uibh.betaU,uibh.BU)

# Step 3.d: Now load the BSSN_in_terms_of_ADM module and perform the same conversion
import BSSN.BSSN_in_terms_of_ADM as BitoA
BitoA.gammabarDD_hDD(                   uibh.gammaDD)
BitoA.trK_AbarDD_aDD(                   uibh.gammaDD,uibh.KDD)
BitoA.LambdabarU_lambdaU__exact_gammaDD(uibh.gammaDD)
BitoA.cf_from_gammaDD(                  uibh.gammaDD)
BitoA.betU_vetU(                        uibh.betaU,uibh.BU)
```


```python
# Step 3.e: Perform the consistency check
def compare(q1, q2, q1name, q2name):
    if sp.simplify(q1 - q2) != 0:
        print("Error: "+q1name+" - "+q2name+" = "+str(sp.simplify(q1 - q2)))
        sys.exit(1) 

print("Consistency check between this tutorial notebook and BSSN.BSSN_in_terms_of_ADM NRPy+ module:")

compare(cf, BitoA.cf, "cf", "BitoA.cf")
compare(trK, BitoA.trK, "trK", "BitoA.trK")
# alpha is the only variable that remains unchanged:
# print("alpha - BitoA.alpha = " + str(alpha - BitoA.alpha))
for i in range(DIM):
    compare(vetU[i], BitoA.vetU[i], "vetU"+str(i), "BitoA.vetU"+str(i))
    compare(betU[i], BitoA.betU[i], "betU"+str(i), "BitoA.betU"+str(i))
    compare(lambdaU[i], BitoA.lambdaU[i], "lambdaU"+str(i), "BitoA.lambdaU"+str(i))
    for j in range(DIM):
        compare(hDD[i][j], BitoA.hDD[i][j], "hDD"+str(i)+str(j), "BitoA.hDD"+str(i)+str(j))
        compare(aDD[i][j], BitoA.aDD[i][j], "aDD"+str(i)+str(j), "BitoA.aDD"+str(i)+str(j))

print("ALL TESTS PASS")
```

    Consistency check between this tutorial notebook and BSSN.BSSN_in_terms_of_ADM NRPy+ module:
    ALL TESTS PASS


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename [Tutorial-BSSN_in_terms_of_ADM.pdf](Tutorial-BSSN_in_terms_of_ADM.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_in_terms_of_ADM")
```

    Created Tutorial-BSSN_in_terms_of_ADM.tex, and compiled LaTeX file to PDF
        file Tutorial-BSSN_in_terms_of_ADM.pdf

