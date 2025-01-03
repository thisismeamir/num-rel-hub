<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Converting between the 4-metric $g_{\mu\nu}$ and ADM variables $\left\{\gamma_{ij}, \alpha, \beta^i\right\}$ or BSSN variables $\left\{h_{ij}, {\rm cf}, \alpha, {\rm vet}^i\right\}$
## Author: Zach Etienne

## This notebook presents the NRPy+ Python module [BSSN/ADMBSSN_tofrom_4metric.py](../edit/BSSN/ADMBSSN_tofrom_4metric.py), which facilitates conversion between the 4-metric and ADM or BSSN variables. The notebook also outlines the steps to write ADM/BSSN metric quantities in terms of 4-metric $g_{\mu\nu}$, while excluding extrinsic curvature $K_{ij}$, or the BSSN $\bar{A}_{ij}$, $K$.

**Notebook Status:** <font color='orange'><b> Self-validated, some additional tests performed </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). In addition, the construction of $g_{\mu\nu}$ and $g^{\mu\nu}$ from BSSN variables has passed the test $g^{\mu\nu}g_{\mu\nu}=4$ [below](#validationcontraction). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Source Code for this module: [BSSN/ADMBSSN_tofrom_4metric.py](../edit/BSSN/ADMBSSN_tofrom_4metric.py)

[comment]: <> (Introduction: TODO)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#setup_ADM_quantities): `setup_ADM_quantities(inputvars)`: If `inputvars="ADM"` declare ADM quantities $\left\{\gamma_{ij},\beta^i,\alpha\right\}$; if `inputvars="ADM"` define ADM quantities in terms of BSSN quantities
1. [Step 2](#admbssn_to_fourmetric): Write 4-metric $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ in terms of ADM or BSSN quantities
    1. [Step 2.a](#admbssn_to_fourmetric_lower): 4-metric $g_{\mu\nu}$ in terms of ADM or BSSN quantities
    1. [Step 2.b](#admbssn_to_fourmetric_inv): 4-metric inverse $g^{\mu\nu}$ in terms of ADM or BSSN quantities
    1. [Step 2.c](#validationcontraction): Validation check: Confirm $g_{\mu\nu}g^{\mu\nu}=4$
1. [Step 3](#fourmetric_to_admbssn): Write ADM/BSSN metric quantities in terms of 4-metric $g_{\mu\nu}$ (Excludes extrinsic curvature $K_{ij}$ or the BSSN $\bar{A}_{ij}$, $K$)
    1. [Step 3.a](#adm_ito_fourmetric_validate): ADM in terms of 4-metric validation: Confirm $\gamma_{ij}\gamma^{ij}=3$
    1. [Step 3.b](#bssn_ito_fourmetric_validate): BSSN in terms of 4-metric validation: Confirm $\bar{\gamma}_{ij}\bar{\gamma}^{ij}=3$
1. [Step 4](#code_validation): Code Validation against `BSSN.ADMBSSN_tofrom_4metric` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='setup_ADM_quantities'></a>

# Step 1: `setup_ADM_quantities(inputvars)`: If `inputvars="ADM"` declare ADM quantities $\left\{\gamma_{ij},\beta^i,\alpha\right\}$; if `inputvars="ADM"` define ADM quantities in terms of BSSN quantities \[Back to [top](#toc)\]
$$\label{setup_ADM_quantities}$$


```python
import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import sys

def setup_ADM_quantities(inputvars):
    if inputvars == "ADM":
        gammaDD = ixp.declarerank2("gammaDD", "sym01")
        betaU = ixp.declarerank1("betaU")
        alpha = sp.symbols("alpha", real=True)
    elif inputvars == "BSSN":
        import BSSN.ADM_in_terms_of_BSSN as AitoB

        # Construct gamma_{ij} in terms of cf & gammabar_{ij}
        AitoB.ADM_in_terms_of_BSSN()
        gammaDD = AitoB.gammaDD
        # Next construct beta^i in terms of vet^i and reference metric quantities
        import BSSN.BSSN_quantities as Bq

        Bq.BSSN_basic_tensors()
        betaU = Bq.betaU
        alpha = sp.symbols("alpha", real=True)
    else:
        print("inputvars = " + str(inputvars) + " not supported. Please choose ADM or BSSN.")
        sys.exit(1)
    return gammaDD,betaU,alpha
```

<a id='admbssn_to_fourmetric'></a>

# Step 2: Write 4-metric $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ in terms of ADM or BSSN variables \[Back to [top](#toc)\]
$$\label{admbssn_to_fourmetric}$$

<a id='admbssn_to_fourmetric_lower'></a>

## Step 2.a: 4-metric $g_{\mu\nu}$ in terms of ADM or BSSN variables \[Back to [top](#toc)\]
$$\label{admbssn_to_fourmetric_lower}$$

Given ADM variables $\left\{\gamma_{ij},\beta^i,\alpha \right\}$, which themselves may be written in terms of the rescaled BSSN curvilinear variables $\left\{h_{ij},{\rm cf},\mathcal{V}^i,\alpha \right\}$ for our chosen reference metric via simple function calls to `ADM_in_terms_of_BSSN()` and `BSSN_quantities.BSSN_basic_tensors()`, we are to construct the 4-metric $g_{\mu\nu}$. 

We accomplish this via Eq. 2.122 (which can be trivially derived from the ADM 3+1 line element) of Baumgarte & Shapiro's *Numerical Relativity* (henceforth B&S):
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix},
$$
where the shift vector $\beta^i$ is lowered via (Eq. 2.121):

$$\beta_k = \gamma_{ik} \beta^i.$$


```python
def g4DD_ito_BSSN_or_ADM(inputvars):
    # Step 0: Declare g4DD as globals, to make interfacing with other modules/functions easier
    global g4DD

    # Step 1: Check that inputvars is set to a supported value
    gammaDD,betaU,alpha = setup_ADM_quantities(inputvars)

    # Step 2: Compute g4DD = g_{mu nu}:
    # To get \gamma_{\mu \nu} = gamma4DD[mu][nu], we'll need to construct the 4-metric, using Eq. 2.122 in B&S:
    g4DD = ixp.zerorank2(DIM=4)

    # Step 2.a: Compute beta_i via Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j] * betaU[j]

    # Step 2.b: Compute beta_i beta^i, the beta contraction.
    beta2 = sp.sympify(0)
    for i in range(3):
        beta2 += betaU[i] * betaD[i]

    # Step 2.c: Construct g4DD via Eq. 2.122 in B&S
    g4DD[0][0] = -alpha ** 2 + beta2
    for mu in range(1, 4):
        g4DD[mu][0] = g4DD[0][mu] = betaD[mu - 1]
    for mu in range(1, 4):
        for nu in range(1, 4):
            g4DD[mu][nu] = gammaDD[mu - 1][nu - 1]
```

<a id='admbssn_to_fourmetric_inv'></a>

## Step 2.b: Inverse 4-metric $g^{\mu\nu}$ in terms of ADM or BSSN variables \[Back to [top](#toc)\]
$$\label{admbssn_to_fourmetric_inv}$$ 

B&S also provide a convenient form for the inverse 4-metric (Eq. 2.119; also Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf)):
$$
g^{\mu\nu} = \gamma^{\mu\nu} - n^\mu n^\nu = 
\begin{pmatrix} 
-\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
\frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
\end{pmatrix},
$$
where the unit normal vector to the hypersurface is given by $n^{\mu} = \left(\alpha^{-1},-\beta^i/\alpha\right)$.


```python
def g4UU_ito_BSSN_or_ADM(inputvars):
    # Step 0: Declare g4UU as globals, to make interfacing with other modules/functions easier
    global g4UU

    # Step 1: Check that inputvars is set to a supported value
    gammaDD,betaU,alpha = setup_ADM_quantities(inputvars)

    # Step 2: Compute g4UU = g_{mu nu}:
    # To get \gamma^{\mu \nu} = gamma4UU[mu][nu], we'll need to use Eq. 2.119 in B&S.
    g4UU = ixp.zerorank2(DIM=4)

    # Step 3: Construct g4UU = g^{mu nu}
    # Step 3.a: Compute gammaUU based on provided gammaDD:
    gammaUU, _gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)

    # Then evaluate g4UU:
    g4UU = ixp.zerorank2(DIM=4)

    g4UU[0][0] = -1 / alpha**2
    for mu in range(1,4):
        g4UU[0][mu] = g4UU[mu][0] = betaU[mu-1]/alpha**2
    for mu in range(1,4):
        for nu in range(1,4):
            g4UU[mu][nu] = gammaUU[mu-1][nu-1] - betaU[mu-1]*betaU[nu-1]/alpha**2
```

<a id='validationcontraction'></a>

## Step 2.c: Validation check: Confirm $g_{\mu\nu}g^{\mu\nu}=4$ \[Back to [top](#toc)\]
$$\label{validationcontraction}$$ 

Next we compute $g^{\mu\nu} g_{\mu\nu}$ as a validation check. It should equal 4:


```python
g4DD_ito_BSSN_or_ADM("BSSN")
g4UU_ito_BSSN_or_ADM("BSSN")
contraction = 0
for mu in range(4):
    for nu in range(4):
        contraction += g4DD[mu][nu]*g4UU[mu][nu]
if sp.simplify(contraction) == sp.sympify(4):
    print("TEST PASSED!")
else:
    print("TEST FAILED: "+str(contraction)+" does not apparently equal 4.")
    sys.exit(1)
```

    TEST PASSED!


<a id='fourmetric_to_admbssn'></a>

# Step 3: Write ADM/BSSN metric quantities in terms of 4-metric $g_{\mu\nu}$ (Excludes extrinsic curvature $K_{ij}$, the BSSN $a_{ij}$, $K$, and $\lambda^i$) \[Back to [top](#toc)\]
$$\label{fourmetric_to_admbssn}$$ 

Given $g_{\mu\nu}$, we now compute ADM/BSSN metric quantities, excluding extrinsic curvature. 

Let's start by computing the ADM quantities in terms of the 4-metric $g_{\mu\nu}$.

Recall that
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$

From this equation, we immediately obtain $\gamma_{ij}$. However we need $\beta^i$ and $\alpha$. After computing the inverse of $\gamma_{ij}$, $\gamma^{ij}$, we raise $\beta_j$ via $\beta^i=\gamma^{ij} \beta_j$ and then compute $\alpha$ via $\alpha = \sqrt{\beta^k \beta_k - g_{00}}$. To convert to BSSN variables $\left\{h_{ij},{\rm cf},\mathcal{V}^i,\alpha \right\}$, we need only convert from ADM via function calls to [`BSSN.BSSN_in_terms_of_ADM`](../edit/BSSN/BSSN_in_terms_of_ADM.py) ([**tutorial**](Tutorial-BSSN_in_terms_of_ADM.ipynb)).


```python
def BSSN_or_ADM_ito_g4DD(inputvars):
    # Step 0: Declare output variables as globals, to make interfacing with other modules/functions easier
    if inputvars == "ADM":
        global gammaDD,betaU,alpha
    elif inputvars == "BSSN":
        global hDD,cf,vetU,alpha
    else:
        print("inputvars = " + str(inputvars) + " not supported. Please choose ADM or BSSN.")
        sys.exit(1)

    # Step 1: declare g4DD as symmetric rank-4 tensor:
    g4DD    = ixp.declarerank2("g4DD","sym01",DIM=4)

    # Step 2: Compute gammaDD & betaD
    betaD   = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    for i in range(3):
        betaD[i] = g4DD[0][i]
        for j in range(3):
            gammaDD[i][j] = g4DD[i+1][j+1]

    # Step 3: Compute betaU
    # Step 3.a: Compute gammaUU based on provided gammaDD
    gammaUU, _gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)

    # Step 3.b: Use gammaUU to raise betaU
    betaU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaU[i] += gammaUU[i][j]*betaD[j]

    # Step 4:   Compute alpha = sqrt(beta^2 - g_{00}):
    # Step 4.a: Compute beta^2 = beta^k beta_k:
    beta_squared = sp.sympify(0)
    for k in range(3):
        beta_squared += betaU[k]*betaD[k]

    # Step 4.b: alpha = sqrt(beta^2 - g_{00}):
    alpha = sp.sqrt(sp.simplify(beta_squared) - g4DD[0][0])

    # Step 5: If inputvars == "ADM", we are finished. Return.
    if inputvars == "ADM":
        return

    # Step 6: If inputvars == "BSSN", convert ADM to BSSN & return hDD, cf,
    import BSSN.BSSN_in_terms_of_ADM as BitoA
    dummyBU = ixp.zerorank1()
    BitoA.gammabarDD_hDD( gammaDD)
    BitoA.cf_from_gammaDD(gammaDD)
    BitoA.betU_vetU(      betaU,dummyBU)
    hDD  = BitoA.hDD
    cf   = BitoA.cf
    vetU = BitoA.vetU
```

<a id='adm_ito_fourmetric_validate'></a>

## Step 3.a: ADM in terms of 4-metric validation: Confirm $\gamma_{ij}\gamma^{ij}=3$  \[Back to [top](#toc)\]
$$\label{adm_ito_fourmetric_validate}$$

Next we compute $\gamma^{ij} \gamma_{ij}$ as a validation check. It should equal 3:


```python
BSSN_or_ADM_ito_g4DD("ADM")
gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)

contraction = sp.sympify(0)
for i in range(3):
    for j in range(3):
        contraction += gammaDD[i][j]*gammaUU[i][j]
if sp.simplify(contraction) == sp.sympify(3):
    print("TEST PASSED!")
else:
    print("TEST FAILED: "+str(contraction)+" does not apparently equal 3.")
    sys.exit(1)
```

    TEST PASSED!


<a id='bssn_ito_fourmetric_validate'></a>

## Step 3.b: BSSN in terms of 4-metric validation: Confirm $\bar{\gamma}_{ij}\bar{\gamma}^{ij}=3$  \[Back to [top](#toc)\]
$$\label{bssn_ito_fourmetric_validate}$$

Next we compute $\bar{\gamma}_{ij}\bar{\gamma}^{ij}$ as a validation check. It should equal 3:


```python
import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","SinhCylindrical")
rfm.reference_metric()

BSSN_or_ADM_ito_g4DD("BSSN")
gammabarDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
        gammabarDD[i][j] = hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]

gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

contraction = sp.sympify(0)
for i in range(3):
    for j in range(3):
        contraction += gammabarDD[i][j]*gammabarUU[i][j]
if sp.simplify(contraction) == sp.sympify(3):
    print("TEST PASSED!")
else:
    print("TEST FAILED: "+str(contraction)+" does not apparently equal 3.")
    sys.exit(1)
```

    TEST PASSED!


<a id='code_validation'></a>

## Step 4: Code Validation against `BSSN.ADMBSSN_tofrom_4metric` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for BrillLindquist initial data between
1. this tutorial and 
2. the NRPy+ [BSSN.ADMBSSN_tofrom_4metric](../edit/BSSN/ADMBSSN_tofrom_4metric.py) module.

By default, we analyze these expressions in SinhCylindrical coordinates, though other coordinate systems may be chosen.


```python
par.set_parval_from_str("reference_metric::CoordSystem","SinhCylindrical")
rfm.reference_metric()

import BSSN.ADMBSSN_tofrom_4metric as AB4m
for inputvars in ["BSSN","ADM"]:
    g4DD_ito_BSSN_or_ADM(inputvars)
    AB4m.g4DD_ito_BSSN_or_ADM(inputvars)
    for i in range(4):
        for j in range(4):
            print(inputvars+" input: g4DD["+str(i)+"]["+str(j)+"] - g4DD_mod["+str(i)+"]["
                  +str(j)+"] = "+str(g4DD[i][j]-AB4m.g4DD[i][j]))

    g4UU_ito_BSSN_or_ADM(inputvars)
    AB4m.g4UU_ito_BSSN_or_ADM(inputvars)
    for i in range(4):
        for j in range(4):
            print(inputvars+" input: g4UU["+str(i)+"]["+str(j)+"] - g4UU_mod["+str(i)+"]["
                  +str(j)+"] = "+str(g4UU[i][j]-AB4m.g4UU[i][j]))

BSSN_or_ADM_ito_g4DD("BSSN")
AB4m.BSSN_or_ADM_ito_g4DD("BSSN")
print("BSSN QUANTITIES (ito 4-metric g4DD)")
print("cf - mod_cf = " + str(cf - AB4m.cf))
print("alpha - mod_alpha = " + str(alpha - AB4m.alpha))
for i in range(3):
    print("vetU["+str(i)+"] - mod_vetU["+str(i)+"] = " + str(vetU[i] - AB4m.vetU[i]))
    for j in range(3):
        print("hDD["+str(i)+"]["+str(j)+"] - mod_hDD["+str(i)+"]["+str(j)+"] = "
              + str(hDD[i][j] - AB4m.hDD[i][j]))

BSSN_or_ADM_ito_g4DD("ADM")
AB4m.BSSN_or_ADM_ito_g4DD("ADM")
print("ADM QUANTITIES (ito 4-metric g4DD)")
print("alpha - mod_alpha = " + str(alpha - AB4m.alpha))
for i in range(3):
    print("betaU["+str(i)+"] - mod_betaU["+str(i)+"] = " + str(betaU[i] - AB4m.betaU[i]))
    for j in range(3):
        print("gammaDD["+str(i)+"]["+str(j)+"] - mod_gammaDD["+str(i)+"]["+str(j)+"] = "
              + str(gammaDD[i][j] - AB4m.gammaDD[i][j]))
```

    BSSN input: g4DD[0][0] - g4DD_mod[0][0] = 0
    BSSN input: g4DD[0][1] - g4DD_mod[0][1] = 0
    BSSN input: g4DD[0][2] - g4DD_mod[0][2] = 0
    BSSN input: g4DD[0][3] - g4DD_mod[0][3] = 0
    BSSN input: g4DD[1][0] - g4DD_mod[1][0] = 0
    BSSN input: g4DD[1][1] - g4DD_mod[1][1] = 0
    BSSN input: g4DD[1][2] - g4DD_mod[1][2] = 0
    BSSN input: g4DD[1][3] - g4DD_mod[1][3] = 0
    BSSN input: g4DD[2][0] - g4DD_mod[2][0] = 0
    BSSN input: g4DD[2][1] - g4DD_mod[2][1] = 0
    BSSN input: g4DD[2][2] - g4DD_mod[2][2] = 0
    BSSN input: g4DD[2][3] - g4DD_mod[2][3] = 0
    BSSN input: g4DD[3][0] - g4DD_mod[3][0] = 0
    BSSN input: g4DD[3][1] - g4DD_mod[3][1] = 0
    BSSN input: g4DD[3][2] - g4DD_mod[3][2] = 0
    BSSN input: g4DD[3][3] - g4DD_mod[3][3] = 0
    BSSN input: g4UU[0][0] - g4UU_mod[0][0] = 0
    BSSN input: g4UU[0][1] - g4UU_mod[0][1] = 0
    BSSN input: g4UU[0][2] - g4UU_mod[0][2] = 0
    BSSN input: g4UU[0][3] - g4UU_mod[0][3] = 0
    BSSN input: g4UU[1][0] - g4UU_mod[1][0] = 0
    BSSN input: g4UU[1][1] - g4UU_mod[1][1] = 0
    BSSN input: g4UU[1][2] - g4UU_mod[1][2] = 0
    BSSN input: g4UU[1][3] - g4UU_mod[1][3] = 0
    BSSN input: g4UU[2][0] - g4UU_mod[2][0] = 0
    BSSN input: g4UU[2][1] - g4UU_mod[2][1] = 0
    BSSN input: g4UU[2][2] - g4UU_mod[2][2] = 0
    BSSN input: g4UU[2][3] - g4UU_mod[2][3] = 0
    BSSN input: g4UU[3][0] - g4UU_mod[3][0] = 0
    BSSN input: g4UU[3][1] - g4UU_mod[3][1] = 0
    BSSN input: g4UU[3][2] - g4UU_mod[3][2] = 0
    BSSN input: g4UU[3][3] - g4UU_mod[3][3] = 0
    ADM input: g4DD[0][0] - g4DD_mod[0][0] = 0
    ADM input: g4DD[0][1] - g4DD_mod[0][1] = 0
    ADM input: g4DD[0][2] - g4DD_mod[0][2] = 0
    ADM input: g4DD[0][3] - g4DD_mod[0][3] = 0
    ADM input: g4DD[1][0] - g4DD_mod[1][0] = 0
    ADM input: g4DD[1][1] - g4DD_mod[1][1] = 0
    ADM input: g4DD[1][2] - g4DD_mod[1][2] = 0
    ADM input: g4DD[1][3] - g4DD_mod[1][3] = 0
    ADM input: g4DD[2][0] - g4DD_mod[2][0] = 0
    ADM input: g4DD[2][1] - g4DD_mod[2][1] = 0
    ADM input: g4DD[2][2] - g4DD_mod[2][2] = 0
    ADM input: g4DD[2][3] - g4DD_mod[2][3] = 0
    ADM input: g4DD[3][0] - g4DD_mod[3][0] = 0
    ADM input: g4DD[3][1] - g4DD_mod[3][1] = 0
    ADM input: g4DD[3][2] - g4DD_mod[3][2] = 0
    ADM input: g4DD[3][3] - g4DD_mod[3][3] = 0
    ADM input: g4UU[0][0] - g4UU_mod[0][0] = 0
    ADM input: g4UU[0][1] - g4UU_mod[0][1] = 0
    ADM input: g4UU[0][2] - g4UU_mod[0][2] = 0
    ADM input: g4UU[0][3] - g4UU_mod[0][3] = 0
    ADM input: g4UU[1][0] - g4UU_mod[1][0] = 0
    ADM input: g4UU[1][1] - g4UU_mod[1][1] = 0
    ADM input: g4UU[1][2] - g4UU_mod[1][2] = 0
    ADM input: g4UU[1][3] - g4UU_mod[1][3] = 0
    ADM input: g4UU[2][0] - g4UU_mod[2][0] = 0
    ADM input: g4UU[2][1] - g4UU_mod[2][1] = 0
    ADM input: g4UU[2][2] - g4UU_mod[2][2] = 0
    ADM input: g4UU[2][3] - g4UU_mod[2][3] = 0
    ADM input: g4UU[3][0] - g4UU_mod[3][0] = 0
    ADM input: g4UU[3][1] - g4UU_mod[3][1] = 0
    ADM input: g4UU[3][2] - g4UU_mod[3][2] = 0
    ADM input: g4UU[3][3] - g4UU_mod[3][3] = 0
    BSSN QUANTITIES (ito 4-metric g4DD)
    cf - mod_cf = 0
    alpha - mod_alpha = 0
    vetU[0] - mod_vetU[0] = 0
    hDD[0][0] - mod_hDD[0][0] = 0
    hDD[0][1] - mod_hDD[0][1] = 0
    hDD[0][2] - mod_hDD[0][2] = 0
    vetU[1] - mod_vetU[1] = 0
    hDD[1][0] - mod_hDD[1][0] = 0
    hDD[1][1] - mod_hDD[1][1] = 0
    hDD[1][2] - mod_hDD[1][2] = 0
    vetU[2] - mod_vetU[2] = 0
    hDD[2][0] - mod_hDD[2][0] = 0
    hDD[2][1] - mod_hDD[2][1] = 0
    hDD[2][2] - mod_hDD[2][2] = 0
    ADM QUANTITIES (ito 4-metric g4DD)
    alpha - mod_alpha = 0
    betaU[0] - mod_betaU[0] = 0
    gammaDD[0][0] - mod_gammaDD[0][0] = 0
    gammaDD[0][1] - mod_gammaDD[0][1] = 0
    gammaDD[0][2] - mod_gammaDD[0][2] = 0
    betaU[1] - mod_betaU[1] = 0
    gammaDD[1][0] - mod_gammaDD[1][0] = 0
    gammaDD[1][1] - mod_gammaDD[1][1] = 0
    gammaDD[1][2] - mod_gammaDD[1][2] = 0
    betaU[2] - mod_betaU[2] = 0
    gammaDD[2][0] - mod_gammaDD[2][0] = 0
    gammaDD[2][1] - mod_gammaDD[2][1] = 0
    gammaDD[2][2] - mod_gammaDD[2][2] = 0


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename [Tutorial-ADMBSSN_tofrom_4metric.pdf](Tutorial-ADMBSSN_tofrom_4metric.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADMBSSN_tofrom_4metric")
```

    Created Tutorial-ADMBSSN_tofrom_4metric.tex, and compiled LaTeX file to PDF
        file Tutorial-ADMBSSN_tofrom_4metric.pdf

