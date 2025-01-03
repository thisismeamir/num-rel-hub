<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [BSSN](http://www2.yukawa.kyoto-u.ac.jp/~yuichiro.sekiguchi/3+1.pdf) Stress-Energy Source Terms

## Author: Zach Etienne

## This module constructs the BSSN stress-energy source terms, in terms of $T^{\mu\nu}$, as prescribed in the reference metric approach of [Baumgarte, Montero, Cordero-Carrión, and Müller (2012)](https://arxiv.org/abs/1211.6632), which builds upon the covariant Lagrangian BSSN formalism of [Brown (2009)](https://arxiv.org/abs/0902.3652).

**Notebook Status:** <font color='orange'><b> Self-validated </b></font>

**Validation Notes:** None yet.

### NRPy+ Source Code for this module: [BSSN/BSSN_stress_energy_source_terms.py](../edit/BSSN/BSSN_stress_energy_source_terms.py)


## Introduction

In [the NRPy+ tutorial on the BSSN formulation](Tutorial-BSSN_formulation.ipynb) we outlined the BSSN formulation of Einstein's equations *in the absence of stress-energy* (i.e., in a vacuum where Einstein's equations reduce to $G^{\mu\nu}=0$). When $T^{\mu\nu}$ is nonzero, stress-energy source terms must appear on the right-hand sides of the BSSN equations in order to ensure Einstein's equations of general relativity are satisfied.

Analyzing Eqs. 9 of [Baumgarte, Montero, Cordero-Carrión, and Müller](https://arxiv.org/pdf/1211.6632.pdf), we see that adding stress-energy source terms $T_{\mu\nu}$ to Einstein's equations in vacuum simply adjusts the right-hand sides of the $\partial_t \bar{A}_{ij}$, $\partial_t K$, and $\partial_t \bar{\Lambda}^i$ equations as follows:


\begin{eqnarray}
\ \partial_t \bar{A}_{ij} &=& \left[\partial_t \bar{A}_{ij}\right]_{\rm vacuum}\ {\color{blue}{-\ 8\pi \alpha e^{-4\phi} \left(S_{ij}\right)^{\rm TF}}} \\
\partial_t K &=& \left[\partial_t K\right]_{\rm vacuum}\ {\color{blue}{+\ 4\pi \alpha (\rho + S)}} \\
\partial_t \bar{\Lambda}^i &=& \left[\partial_t \bar{\Lambda}^{i}\right]_{\rm vacuum}\ {\color{blue}{-\ 16\pi \alpha \bar{\gamma}^{ij} S_j}},
\end{eqnarray}

where $\rho$, $S$, $S_i$, and $S_{ij}$ are related to the stress-energy tensor $T^{\mu\nu}$ as follows (Eq. 10 of [Baumgarte, Montero, Cordero-Carrión, and Müller](https://arxiv.org/pdf/1211.6632.pdf)):

\begin{eqnarray}
\ S_{ij} &=& \gamma_{i \mu} \gamma_{j \nu} T^{\mu \nu} \\
S_{i} &=& -\gamma_{i\mu} n_\nu T^{\mu\nu} \\
S &=& \gamma^{ij} S_{ij} \\
\rho &=& n_\mu n_\nu T^{\mu\nu},
\end{eqnarray}

the unit normal one-form on each spatial slice $n_{\mu}$ is given by Eq. 10 of [Baumgarte, Montero, Cordero-Carrión, and Müller](https://arxiv.org/pdf/1211.6632.pdf)):

$$
n_\mu = (-\alpha,0,0,0),
$$

and Baumgarte & Shapiro Eq. 2.27 gives $\gamma_{\mu\nu}$:

$$\gamma_{\mu\nu} = g_{\mu\nu} + n_\mu n_\nu.$$

Further, analyzing Eqs. 13 & 14 of [Baumgarte, Montero, Cordero-Carrión, and Müller](https://arxiv.org/pdf/1211.6632.pdf) we find that adding stress-energy source terms $T_{\mu\nu}$ to Einstein's equations in vacuum adjusts the BSSN constraint equations as follows:
\begin{eqnarray}
\ \mathcal{H} &=& \left[\mathcal{H}\right]_{\rm vacuum}\ {\color{blue}{-\ 16\pi \rho}} \\
\mathcal{M}^i &=& \left[\mathcal{M}^i\right]_{\rm vacuum}\ {\color{blue}{-\ 8\pi S^i}},
\end{eqnarray}

This module will construct expressions for $S_{ij}$, $S_i$, $S$, and $\rho$ in terms of $T^{\mu\nu}$, and also add the necessary terms to the BSSN RHSs and constraints.

### A Note on Notation

As is standard in NRPy+, 

* Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
* Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.

As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
1. [Step 2](#bssn_sourceterms): BSSN source terms, in terms of $T^{\mu\nu}$
    1. [Step 2.a](#gamma4dd): Define `gamma4DD[mu][nu]` = $g_{\mu \nu} + n_{\mu} n_{\nu}$
    1. [Step 2.b](#t4uu): Declare `T4UU[mu][nu]`=$T^{\mu\nu}$
    1. [Step 2.c](#define_bssn_sourceterms): Define BSSN source terms    
1. [Step 3](#add_bssn_sourceterms_to_rhss): Add BSSN source terms to BSSN RHSs
1. [Step 4](#add_bssn_sourceterms_to_constraints): Add BSSN source terms to BSSN Constraints
1. [Step 5](#code_validation): Code Validation against `BSSN.BSSN_stress_energy_source_terms` NRPy+ module
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules & set up reference metric \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp                         # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par             # NRPy+: Parameter interface
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm             # NRPy+: Reference metric support
import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions
import BSSN.ADM_in_terms_of_BSSN as AitoB  # NRPy+: ADM quantities in terms of BSSN quantities

# Step 1.a: Set up reference metric. We'll choose SinhSpherical here, but
#           could choose any CoordSystem defined in reference_metric.py:
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

thismodule = "BSSN.BSSN_stress_energy_source_terms"
```

<a id='bssn_sourceterms'></a>

# Step 2: BSSN source terms, in terms of $T^{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{bssn_sourceterms}$$

<a id='gamma4dd'></a>

## Step 2.a: Define `gamma4DD[mu][nu]` = $g_{\mu \nu} + n_{\mu} n_{\nu}$ \[Back to [top](#toc)\]
$$\label{gamma4dd}$$


```python
# Step 2.a: Define gamma4DD[mu][nu] = g_{mu nu} + n_{mu} n_{nu}
alpha = sp.symbols("alpha",real=True)
zero  = sp.sympify(0)
n4D      = [-alpha, zero, zero ,zero]
AB4m.g4DD_ito_BSSN_or_ADM("BSSN")

gamma4DD = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        gamma4DD[mu][nu] = AB4m.g4DD[mu][nu] + n4D[mu]*n4D[nu]
```

<a id='t4uu'></a>

## Step 2.b: Declare `T4UU[mu][nu]`=$T^{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{t4uu}$$


```python
# Step 2.b: Declare T4UU
T4UU = ixp.declarerank2("T4UU","sym01",DIM=4)
```

<a id='define_bssn_sourceterms'></a>

## Step 2.c: Define BSSN source terms \[Back to [top](#toc)\]
$$\label{define_bssn_sourceterms}$$

Recall from above, we have:
\begin{eqnarray}
\ S_{ij} &=& \gamma_{i \mu} \gamma_{j \nu} T^{\mu \nu} \\
S_{i} &=& -\gamma_{i\mu} n_\nu T^{\mu\nu} \\
S &=& \gamma^{ij} S_{ij} \\
\rho &=& n_\mu n_\nu T^{\mu\nu}.
\end{eqnarray}


```python
# Step 2.c: Define BSSN source terms
# Step 2.c.i: S_{ij} = gamma_{i mu} gamma_{j nu} T^{mu nu}
SDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        for mu in range(4):
            for nu in range(4):
                SDD[i][j] += gamma4DD[i+1][mu] * gamma4DD[j+1][nu] * T4UU[mu][nu]
# Step 2.c.ii: S_{i} = -gamma_{i mu} n_{nu} T^{mu nu}
SD = ixp.zerorank1()
for i in range(3):
    for mu in range(4):
        for nu in range(4):
            SD[i] += - gamma4DD[i+1][mu] * n4D[nu] * T4UU[mu][nu]
# Step 2.c.iii: S = gamma^{ij} S_{ij}
AitoB.ADM_in_terms_of_BSSN()
S = zero
for i in range(3):
    for j in range(3):
        S += AitoB.gammaUU[i][j]*SDD[i][j]
# Step 2.c.iv: rho = n_{mu} n_{nu} T^{mu nu}
rho = zero
for mu in range(4):
    for nu in range(4):
        rho += n4D[mu]*n4D[nu]*T4UU[mu][nu]
```

<a id='add_bssn_sourceterms_to_rhss'></a>

# Step 3: Add BSSN source terms to BSSN RHSs \[Back to [top](#toc)\]
$$\label{add_bssn_sourceterms_to_rhss}$$

Recall from above we need to make the following modifications:
\begin{eqnarray}
\ \partial_t \bar{A}_{ij} &=& \left[\partial_t \bar{A}_{ij}\right]_{\rm vacuum}\ {\color{blue}{-\ 8\pi \alpha e^{-4\phi} \left(S_{ij}\right)^{\rm TF}}} \\
\partial_t K &=& \left[\partial_t K\right]_{\rm vacuum}\ {\color{blue}{+\ 4\pi \alpha (\rho + S)}} \\
\partial_t \bar{\Lambda}^i &=& \left[\partial_t \bar{\Lambda}^{i}\right]_{\rm vacuum}\ {\color{blue}{-\ 16\pi \alpha \bar{\gamma}^{ij} S_j}},
\end{eqnarray}

where $$\left(S_{ij}\right)^{\rm TF} = S_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{km} S_{km}.$$

*Exercise to student:* Prove that replacing the $\bar{\gamma}_{ij}$ and $\bar{\gamma}^{km}$ with $\gamma_{ij}$ and $\gamma^{km}$, respectively, results in exactly the same expression for $\left(S_{ij}\right)^{\rm TF}$.


```python
# Step 3: Add BSSN stress-energy source terms to BSSN RHSs
import BSSN.BSSN_quantities as Bq
# Can't #declare M_PI here, as it is not SIMD-compatible.
PI = par.Cparameters("REAL", thismodule, ["PI"], "3.14159265358979323846264338327950288")
alpha = sp.symbols("alpha",real=True)
zero  = sp.sympify(0)

# Step 3.a: Initialize RHS source terms to zero.
sourceterm_trK_rhs     = zero
sourceterm_a_rhsDD     = ixp.zerorank2()
sourceterm_lambda_rhsU = ixp.zerorank1()

# Step 3.b: trK_rhs
sourceterm_trK_rhs = 4*PI*alpha*(rho + S)

# Step 3.c: Abar_rhsDD:
# Step 3.c.i: Compute trace-free part of S_{ij}:
Bq.BSSN_basic_tensors() # Sets gammabarDD
gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(Bq.gammabarDD) # Set gammabarUU
tracefree_SDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        tracefree_SDD[i][j] = SDD[i][j]
for i in range(3):
    for j in range(3):
        for k in range(3):
            for m in range(3):
                tracefree_SDD[i][j] += -sp.Rational(1,3)*Bq.gammabarDD[i][j]*gammabarUU[k][m]*SDD[k][m]

# Step 3.c.ii: Define exp_m4phi = e^{-4 phi}
Bq.phi_and_derivs()
# Step 3.c.iii: Evaluate RHS
for i in range(3):
    for j in range(3):
        Abar_rhsDDij             = -8*PI*alpha*Bq.exp_m4phi*tracefree_SDD[i][j]
        sourceterm_a_rhsDD[i][j] =  Abar_rhsDDij / rfm.ReDD[i][j]

# Step 3.d: Stress-energy part of Lambdabar_rhsU = stressenergy_Lambdabar_rhsU
sourceterm_Lambdabar_rhsU = ixp.zerorank1()
for i in range(3):
    for j in range(3):
        sourceterm_Lambdabar_rhsU[i] += -16*PI*alpha*gammabarUU[i][j]*SD[j]
for i in range(3):
    sourceterm_lambda_rhsU[i] = sourceterm_Lambdabar_rhsU[i] / rfm.ReU[i]
```

<a id='add_bssn_sourceterms_to_constraints'></a>

# Step 4: Add BSSN source terms to BSSN Constraints \[Back to [top](#toc)\]
$$\label{add_bssn_sourceterms_to_constraints}$$

Recall from above we need to make the following modifications:
\begin{eqnarray}
\ \mathcal{H} &=& \left[\mathcal{H}\right]_{\rm vacuum}\ {\color{blue}{-\ 16\pi \rho}} \\
\mathcal{M}^i &=& \left[\mathcal{M}^i\right]_{\rm vacuum}\ {\color{blue}{-\ 8\pi S^i}},
\end{eqnarray}

where 
$$
S^i = \gamma^{ij} S_j,
$$
and $\gamma^{ij}$ is the inverse ADM 3-metric.


```python
# Step 4: Add BSSN stress-energy source terms to BSSN constraints
# Step 4.a: Initialize constraint source terms to zero.
sourceterm_H  = sp.sympify(0)
sourceterm_MU = ixp.zerorank1()

# Step 4.b: Add source term to the Hamiltonian constraint H
sourceterm_H = -16*PI*rho

# Step 4.c: Add source term to the momentum constraint M^i
# Step 4.c.i: Compute gammaUU in terms of BSSN quantities
import BSSN.ADM_in_terms_of_BSSN as AitoB
AitoB.ADM_in_terms_of_BSSN() # Provides gammaUU
# Step 4.c.ii: Raise S_i
SU = ixp.zerorank1()
for i in range(3):
    for j in range(3):
        SU[i] += AitoB.gammaUU[i][j]*SD[j]
# Step 4.c.iii: Add source term to momentum constraint & rescale:
for i in range(3):
    sourceterm_MU[i] = -8 * PI * SU[i] / rfm.ReU[i]
```

<a id='code_validation'></a>

# Step 5: Code Validation against `BSSN.BSSN_stress_energy_source_terms` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_stressenergy_source_terms](../edit/BSSN/BSSN_stressenergy_source_terms.py) module.

By default, we analyze these expressions in SinhSpherical coordinates, though other coordinate systems may be chosen.


```python
# Step 5: Code Validation against BSSN.BSSN_stress_energy_source_terms NRPy+ module

# We already have SymPy expressions for BSSN source terms
#         in terms of other SymPy variables.
#
#         Here, we will use the above-defined BSSN stress-energy source term expressions
#         to validate against the same expressions in the
#         BSSN/BSSN_stress_energy_source_terms.py file, to ensure consistency between
#         this tutorial and the module itself.
import BSSN.BSSN_stress_energy_source_terms as Bsest

print("Consistency check between BSSN_stress_energy_source_terms tutorial and NRPy+ module: ALL SHOULD BE ZERO.")

print("STRESS-ENERGY SOURCE TERMS:")
Bsest.stress_energy_source_terms_ito_T4UU_and_ADM_or_BSSN_metricvars("BSSN")
print("rho - Bsest.rho = " + str(rho - Bsest.rho))
print("S   - Bsest.S   = " + str(S   - Bsest.S))
for i in range(3):
    print("SD["+str(i)+"] - Bsest.SD["+str(i)+"] = " + str(SD[i] - Bsest.SD[i]))
for i in range(3):
    for j in range(3):
        print("SDD["+str(i)+"]["+str(j)+"] - Bsest.SDD["+str(i)+"]["+str(j)+"] = " + str(SDD[i][j] - Bsest.SDD[i][j]))

print("\nBSSN RHSs SOURCE TERMS:")
Bsest.BSSN_source_terms_for_BSSN_RHSs()
print("sourceterm_trK_rhs - Bsest.sourceterm_trK_rhs = " + str(sourceterm_trK_rhs - Bsest.sourceterm_trK_rhs))
for i in range(3):
    for j in range(3):
        print("sourceterm_a_rhsDD["+str(i)+"]["+str(j)+"] - Bsest.sourceterm_a_rhsDD["+str(i)+"]["+str(j)+"] = " +
              str(sourceterm_a_rhsDD[i][j] - Bsest.sourceterm_a_rhsDD[i][j]))
for i in range(3):
    print("sourceterm_lambda_rhsU["+str(i)+"] - Bsest.sourceterm_lambda_rhsU["+str(i)+"] = " +
          str(sourceterm_lambda_rhsU[i] - Bsest.sourceterm_lambda_rhsU[i]))


print("\nBSSN CONSTRAINTS SOURCE TERMS:")
Bsest.BSSN_source_terms_for_BSSN_constraints()
print("sourceterm_H - Bsest.sourceterm_H = " + str(sourceterm_H - Bsest.sourceterm_H))
for i in range(3):
    print("sourceterm_MU["+str(i)+"] - Bsest.sourceterm_MU["+str(i)+"] = " +
          str(sourceterm_MU[i] - Bsest.sourceterm_MU[i]))
```

    Consistency check between BSSN_stress_energy_source_terms tutorial and NRPy+ module: ALL SHOULD BE ZERO.
    STRESS-ENERGY SOURCE TERMS:
    rho - Bsest.rho = 0
    S   - Bsest.S   = 0
    SD[0] - Bsest.SD[0] = 0
    SD[1] - Bsest.SD[1] = 0
    SD[2] - Bsest.SD[2] = 0
    SDD[0][0] - Bsest.SDD[0][0] = 0
    SDD[0][1] - Bsest.SDD[0][1] = 0
    SDD[0][2] - Bsest.SDD[0][2] = 0
    SDD[1][0] - Bsest.SDD[1][0] = 0
    SDD[1][1] - Bsest.SDD[1][1] = 0
    SDD[1][2] - Bsest.SDD[1][2] = 0
    SDD[2][0] - Bsest.SDD[2][0] = 0
    SDD[2][1] - Bsest.SDD[2][1] = 0
    SDD[2][2] - Bsest.SDD[2][2] = 0
    
    BSSN RHSs SOURCE TERMS:
    sourceterm_trK_rhs - Bsest.sourceterm_trK_rhs = 0
    sourceterm_a_rhsDD[0][0] - Bsest.sourceterm_a_rhsDD[0][0] = 0
    sourceterm_a_rhsDD[0][1] - Bsest.sourceterm_a_rhsDD[0][1] = 0
    sourceterm_a_rhsDD[0][2] - Bsest.sourceterm_a_rhsDD[0][2] = 0
    sourceterm_a_rhsDD[1][0] - Bsest.sourceterm_a_rhsDD[1][0] = 0
    sourceterm_a_rhsDD[1][1] - Bsest.sourceterm_a_rhsDD[1][1] = 0
    sourceterm_a_rhsDD[1][2] - Bsest.sourceterm_a_rhsDD[1][2] = 0
    sourceterm_a_rhsDD[2][0] - Bsest.sourceterm_a_rhsDD[2][0] = 0
    sourceterm_a_rhsDD[2][1] - Bsest.sourceterm_a_rhsDD[2][1] = 0
    sourceterm_a_rhsDD[2][2] - Bsest.sourceterm_a_rhsDD[2][2] = 0
    sourceterm_lambda_rhsU[0] - Bsest.sourceterm_lambda_rhsU[0] = 0
    sourceterm_lambda_rhsU[1] - Bsest.sourceterm_lambda_rhsU[1] = 0
    sourceterm_lambda_rhsU[2] - Bsest.sourceterm_lambda_rhsU[2] = 0
    
    BSSN CONSTRAINTS SOURCE TERMS:
    sourceterm_H - Bsest.sourceterm_H = 0
    sourceterm_MU[0] - Bsest.sourceterm_MU[0] = 0
    sourceterm_MU[1] - Bsest.sourceterm_MU[1] = 0
    sourceterm_MU[2] - Bsest.sourceterm_MU[2] = 0


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_stress_energy_source_terms.pdf](Tutorial-BSSN_stress_energy_source_terms.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_stress_energy_source_terms")
```

    Created Tutorial-BSSN_stress_energy_source_terms.tex, and compiled LaTeX
        file to PDF file Tutorial-BSSN_stress_energy_source_terms.pdf

