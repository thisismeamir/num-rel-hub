<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Using NRPy+ to Construct SymPy expressions for Weyl scalars and invariants in Cartesian coordinates

## Author: Patrick Nelson & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This module adopts the prescription of [Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)](https://arxiv.org/abs/gr-qc/0104063) (henceforth, the "BCL paper") to construct the Weyl scalars $\psi_0$, $\psi_1$, $\psi_2$, $\psi_3$, and $\psi_4$ from an approximate Kinnersley tetrad. 

### It also constructs the corresponding Weyl invariants adopting the same approach as Einstein Toolkit's [Kranc-generated](http://kranccode.org/) [WeylScal4 diagnostic module](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src). We will also follow that thorn's approach for other parts of this code.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** The numerical implementation of expressions constructed in this module have been validated against a trusted code (the WeylScal4 Einstein Toolkit thorn).

### NRPy+ Source Code for this module: 
* [WeylScal4NRPy/WeylScalarInvariants_Cartesian.py](../edit/WeylScal4NRPy/WeylScalarInvariants_Cartesian.py)
* [WeylScal4NRPy/WeylScalars_Cartesian.py](../edit/WeylScal4NRPy/WeylScalars_Cartesian.py)

## Introduction: 
As this module is meant for Cartesian coordinates, all quantities are already rescaled. Further, we assume physical (as opposed to conformal) quantities including the 3-metric $\gamma_{ij}$ and 3-extrinsic curvature $K_{ij}$ are provided as input gridfunctions.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Set core NRPy+ parameters for numerical grids and reference metric
1. [Step 2](#register_gridfunctions): Defining the Levi-Civita Symbol
1. [Step 3](#qktetrad): Defining the Approximate Quasi-Kinnersley Tetrad
1. [Step 4](#tensor): Building the Riemann and extrinsic curvature tensors
1. [Step 5](#psi4): Putting it all together and calculating $\psi_4$
    1. [Step 5.a](#code_validation1): Code Validation against `WeylScal4NRPy.WeylScalars_Cartesian` NRPy+ Module
1. [Step 6](#invariant_scalars): The Invariant Scalars
    1. [Step 6.a](#code_validation2): Code Validation against `WeylScal4NRPy.WeylScalarInvariants_Cartesian` NRPy+ Module
1. [Step 7](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Set core NRPy+ parameters for numerical grids and reference metric \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri                # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import sys                        # Standard Python modules for multiplatform OS-level functions
```


```python
# Step 2: Initialize WeylScalar parameters
thismodule = __name__
# Currently only one option: Approx_QuasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)
par.initialize_param(par.glb_param("char", thismodule, "TetradChoice", "Approx_QuasiKinnersley"))
# The default value will output all psis
par.initialize_param(par.glb_param("bool", thismodule, "output_all_psis", False))
```

<a id='register_gridfunctions'></a>

# Step 2: Declare input and output variable gridfunctions \[Back to [top](#toc)\]
$$\label{register_gridfunctions}$$

We declare gridfunctions quantities defined here depend on, including:

* the physical metric $\gamma_{ij}$, 
* the extrinsic curvature $K_{ij}$, 
* the Cartesian coordinates $(x,y,z)$,

Also, the output gridfunctions are as follows:

* the real and imaginary components of $\psi_4$, and
* the Weyl curvature invariants


```python
# Step 3.a: Set spatial dimension (must be 3 for BSSN)
par.set_parval_from_str("grid::DIM",3)

# Step 3.b: declare the additional gridfunctions (i.e., functions whose values are declared
#          at every grid point, either inside or outside of our SymPy expressions) needed
#         for this thorn:
#           * the physical metric $\gamma_{ij}$,
#           * the extrinsic curvature $K_{ij}$,
#           * the real and imaginary components of $\psi_4$, and
#           * the Weyl curvature invariants:
gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01") # The AUX or EVOL designation is *not*
                                                                                # used in diagnostic modules.
kDD = ixp.register_gridfunctions_for_single_rank2("AUX","kDD", "sym01")
x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])
psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i = gri.register_gridfunctions("AUX",["psi4r","psi4i",
                                                                                                "psi3r","psi3i",
                                                                                                "psi2r","psi2i",
                                                                                                "psi1r","psi1i",
                                                                                                "psi0r","psi0i"])
curvIr,curvIi,curvJr,curvJi,J1curv,J2curv,J3curv,J4curv = gri.register_gridfunctions("AUX",["curvIr","curvIi",
                                                                                            "curvJr","curvJi",
                                                                                            "J1curv","J2curv",
                                                                                            "J3curv","J4curv"])
```

<a id='qktetrad'></a>

# Step 3: Defining the Approximate Quasi-Kinnersley Tetrad \[Back to [top](#toc)\]
$$\label{qktetrad}$$

To define the Weyl scalars, first, a tetrad must be chosen. Below, for compatibility with the [WeylScal4 diagnostic module](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src), we implement the approximate  quasi-Kinnersley tetrad of the [BCL paper](https://arxiv.org/abs/gr-qc/0104063).

We begin with the vectors given in eqs. 5.6 and 5.7 of the [BCL paper](https://arxiv.org/abs/gr-qc/0104063),
\begin{align}
    v_1^a &= [-y,x,0] \\
    v_2^a &= [x,y,z] \\
    v_3^a &= {\rm det}(g)^{1/2} g^{ad} \epsilon_{dbc} v_1^b v_2^c,
\end{align}
and carry out the Gram-Schmidt orthonormalization process. Note that these vectors are initially orthogonal to each other; one is in the $\phi$ direction, one is in $r$, and the third is the cross product of the first two. The vectors $w_i^a$ are placeholders in the code; the final product of the orthonormalization is the vectors $e_i^a$. So,
\begin{align}
e_1^a &= \frac{v_1^a}{\sqrt{\omega_{11}}} \\
e_2^a &= \frac{v_2^a - \omega_{12} e_1^a}{\sqrt{\omega_{22}}} \\
e_3^a &= \frac{v_3^a - \omega_{13} e_1^a - \omega_{23} e_2^a}{\sqrt{\omega_{33}}}, \\
\end{align}
where $\omega_{ij} = v_i^a v_j^b \gamma_{ab}$ needs to be updated between steps (to save resources, we can get away with only calculating components as needed) and uses $e_i^a$ instead of $v_i^a$ if it has been calculated. Recall that $\gamma_{ab}$ was declared as a gridfunction above.

Once we have orthogonal, normalized vectors, we can construct the tetrad itself, again drawing on eqs. 5.6. We could draw on SymPy's built-in tools for complex numbers to build the complex vectors $m^a$ and $(m^*)^a$; however, the final expressions for the Weyl scalars are complex enough that `sp.re()` and `sp.im()` are prohibitively time-consuming. To get around this, we will define the real and imaginary components of $\overset{*}{m}{}^a$, and do the complex algebra by hand. Thus,
\begin{align}
    l^a &= \frac{1}{\sqrt{2}} e_2^a \\
    n^a &= -\frac{1}{\sqrt{2}} e_2^a \\
    m^a &= \frac{1}{\sqrt{2}} (e_3^a + i e_1^a) \\
    \overset{*}{m}{}^a &= \frac{1}{\sqrt{2}} (e_3^a - i e_1^a).
\end{align}

In coding this procedure, we will follow the code from [WeylScal4](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src) very closely. We will also assume that $l^0 = n^0 = \frac{1}{\sqrt{2}}$ and that $m^0 = \overset{*}{m}{}^0 = 0$ (again, following the example of the Kranc-generated WeylScal4). This last assumption in particular will significantly reduce the terms needed to find $\psi_4$.


```python
# Step 4: Set which tetrad is used; at the moment, only one supported option

# The tetrad depends in general on the inverse 3-metric gammaUU[i][j]=\gamma^{ij}
#          and the determinant of the 3-metric (detgamma), which are defined in
#          the following line of code from gammaDD[i][j]=\gamma_{ij}.
tmpgammaUU, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)
detgamma = sp.simplify(detgamma)
gammaUU = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        gammaUU[i][j] = sp.simplify(tmpgammaUU[i][j])

if par.parval_from_str("TetradChoice") == "Approx_QuasiKinnersley":
    # Eqs 5.6 in https://arxiv.org/pdf/gr-qc/0104063.pdf
    xmoved = x# - xorig
    ymoved = y# - yorig
    zmoved = z# - zorig

    # Step 5.a: Choose 3 orthogonal vectors. Here, we choose one in the azimuthal
    #          direction, one in the radial direction, and the cross product of the two.
    # Eqs 5.7
    v1U = ixp.zerorank1()
    v2U = ixp.zerorank1()
    v3U = ixp.zerorank1()
    v1U[0] = -ymoved
    v1U[1] = xmoved# + offset
    v1U[2] = sp.sympify(0)
    v2U[0] = xmoved# + offset
    v2U[1] = ymoved
    v2U[2] = zmoved
    LeviCivitaSymbol_rank3 = ixp.LeviCivitaSymbol_dim3_rank3()
    for a in range(3):
        for b in range(3):
            for c in range(3):
                for d in range(3):
                    v3U[a] += sp.sqrt(detgamma) * gammaUU[a][d] * LeviCivitaSymbol_rank3[d][b][c] * v1U[b] *v2U[c]

    for a in range(3):
        v3U[a] = sp.simplify(v3U[a])

    # Step 5.b: Gram-Schmidt orthonormalization of the vectors.
    # The w_i^a vectors here are used to temporarily hold values on the way to the final vectors e_i^a

    # e_1^a &= \frac{v_1^a}{\omega_{11}}
    # e_2^a &= \frac{v_2^a - \omega_{12} e_1^a}{\omega_{22}}
    # e_3^a &= \frac{v_3^a - \omega_{13} e_1^a - \omega_{23} e_2^a}{\omega_{33}},

    # Normalize the first vector
    w1U = ixp.zerorank1()
    for a in range(3):
        w1U[a] = v1U[a]
    omega11 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega11 += w1U[a] * w1U[b] * gammaDD[a][b]
    e1U = ixp.zerorank1()
    for a in range(3):
        e1U[a] = w1U[a] / sp.sqrt(omega11)

    # Subtract off the portion of the first vector along the second, then normalize
    omega12 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega12 += e1U[a] * v2U[b] * gammaDD[a][b]
    w2U = ixp.zerorank1()
    for a in range(3):
        w2U[a] = v2U[a] - omega12*e1U[a]
    omega22 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega22 += w2U[a] * w2U[b] *gammaDD[a][b]
    e2U = ixp.zerorank1()
    for a in range(3):
        e2U[a] = w2U[a] / sp.sqrt(omega22)

    # Subtract off the portion of the first and second vectors along the third, then normalize
    omega13 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega13 += e1U[a] * v3U[b] * gammaDD[a][b]
    omega23 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega23 += e2U[a] * v3U[b] * gammaDD[a][b]
    w3U = ixp.zerorank1()
    for a in range(3):
        w3U[a] = v3U[a] - omega13*e1U[a] - omega23*e2U[a]
    omega33 = sp.sympify(0)
    for a in range(3):
        for b in range(3):
            omega33 += w3U[a] * w3U[b] * gammaDD[a][b]
    e3U = ixp.zerorank1()
    for a in range(3):
        e3U[a] = w3U[a] / sp.sqrt(omega33)

    # Step 5.c: Construct the tetrad itself.
    # Eqs. 5.6:
    # l^a &= \frac{1}{\sqrt{2}} e_2^a
    # n^a &= -\frac{1}{\sqrt{2}} e_2^a
    # m^a &= \frac{1}{\sqrt{2}} (e_3^a + i e_1^a)
    # \overset{*}{m}{}^a &= \frac{1}{\sqrt{2}} (e_3^a - i e_1^a)
    isqrt2 = 1/sp.sqrt(2)
    ltetU = ixp.zerorank1()
    ntetU = ixp.zerorank1()
    #mtetU = ixp.zerorank1()
    #mtetccU = ixp.zerorank1()
    remtetU = ixp.zerorank1() # SymPy did not like trying to take the real/imaginary parts of such a
    immtetU = ixp.zerorank1() # complicated expression, so we do it ourselves.
    for i in range(3):
        ltetU[i] = isqrt2 * e2U[i]
        ntetU[i] = -isqrt2 * e2U[i]
        remtetU[i] = isqrt2 * e3U[i]
        immtetU[i] = isqrt2 * e1U[i]
    nn = isqrt2

else:
    print("Error: TetradChoice == "+par.parval_from_str("TetradChoice")+" unsupported!")
    sys.exit(1)
```

<a id='tensor'></a>

# Step 4: Building the Riemann and extrinsic curvature tensors \[Back to [top](#toc)\]
$$\label{tensor}$$

Now that we have the tetrad in place, we can contract it with the Weyl tensor to obtain the Weyl scalars. Naturally, we must first construct the Weyl tensor to do so; we will not do this directly, instead following the example of the [BCL paper](https://arxiv.org/abs/gr-qc/0104063) and the original WeylScal4. We will first build the Christoffel symbols,
$$\Gamma^i_{kl} = \frac{1}{2} \gamma^{im} (\gamma_{mk,l} + \gamma_{ml,k} - \gamma_{kl,m}).
$$


```python
#Step 5: Declare and construct the second derivative of the metric.
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")

# Define the Christoffel symbols
GammaUDD = ixp.zerorank3(3)
for i in range(3):
    for k in range(3):
        for l in range(3):
            for m in range(3):
                GammaUDD[i][k][l] += (sp.Rational(1,2))*gammaUU[i][m]*\
                                     (gammaDD_dD[m][k][l] + gammaDD_dD[m][l][k] - gammaDD_dD[k][l][m])
```

We will also need the Riemann curvature tensor,
$$
R_{abcd} = \frac{1}{2} (\gamma_{ad,cb}+\gamma_{bc,da}-\gamma_{ac,bd}-\gamma_{bd,ac}) + \gamma_{je} \Gamma^{j}_{bc}\Gamma^{e}_{ad} - \gamma_{je} \Gamma^{j}_{bd} \Gamma^{e}_{ac},
$$
since several terms in our expression for $\psi_4$ are contractions of this tensor.
To do this, we need second derivatives of the metric tensor, $\gamma_{ab,cd}$, using the finite differencing functionality in NRPy+. 


```python
# Step 6.a: Declare and construct the Riemann curvature tensor:
gammaDD_dDD = ixp.declarerank4("gammaDD_dDD","sym01_sym23")
RiemannDDDD = ixp.zerorank4()
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                RiemannDDDD[a][b][c][d] += (gammaDD_dDD[a][d][c][b] +
                                            gammaDD_dDD[b][c][d][a] -
                                            gammaDD_dDD[a][c][b][d] -
                                            gammaDD_dDD[b][d][a][c]) * sp.Rational(1,2)
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                for e in range(3):
                    for j in range(3):
                        RiemannDDDD[a][b][c][d] +=  gammaDD[j][e] * GammaUDD[j][b][c] * GammaUDD[e][a][d] - \
                                                    gammaDD[j][e] * GammaUDD[j][b][d] * GammaUDD[e][a][c]
```

We will need the trace of the extrinsic curvature tensor $K_{ij}$, which can be computed as usual: $K = \gamma^{ij} K_{ij}$.


```python
# Step 6.b: We also need the extrinsic curvature tensor $K_{ij}$. This can be built from quantities from BSSN_RHSs
#           For now, we assume this is a gridfunction (We assume the ADM formalism for now).
#extrinsicKDD = ixp.zerorank2()
#for i in range(3):
#    for j in range(3):
#        extrinsicKDD[i][j] = (bssn.AbarDD[i][j] + sp.Rational(1,3)*gammaDD[i][j]*bssn.trK)/bssn.exp_m4phi
# We will, however, need to calculate the trace of K seperately:
trK = sp.sympify(0)
for i in range(3):
    for j in range(3):
        trK += gammaUU[i][j] * kDD[i][j]
```

<a id='psi4'></a>

# Step 5: Putting it all together and calculating $\psi_4$ \[Back to [top](#toc)\]
$$\label{psi4}$$

We do not need to explicitly build the Weyl tensor itself, because the [BCL paper](https://arxiv.org/abs/gr-qc/0104063) shows that, for the Weyl tensor $C_{ijkl}$,
\begin{align}
\psi_4 =&\  C_{ijkl} n^i \overset{*}{m}{}^j n^k \overset{*}{m}{}^l\\
=&\ (R_{ijkl} + 2K_{i[k}K_{l]j}) n^i \overset{*}{m}{}^j n^k \overset{*}{m}{}^l \\
&- 8 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}) n^{[0} \overset{*}{m}{}^{j]} n^k \overset{*}{m}{}^l \\
&+ 4 (R_{jl} - K_{jp} K^p_l + KK_{jl}) n^{[0} \overset{*}{m}{}^{j]} n^{[0} \overset{*}{m}{}^{l]}.
\end{align}

Note here the brackets around pairs of indices. This indicates the antisymmetric part of a tensor; that is, for some arbitrary tensor $A_{ij}$, $A_{[ij]} = \frac{1}{2}(A_{ij}-A_{ji})$. This applies identically for indices belonging to separate tensors as well as superscripts in place of subscripts.

The other Weyl scalars from the [BCL paper](https://arxiv.org/abs/gr-qc/0104063), appendix A, that we may want to consider are 
\begin{align}
\psi_3 =&\ (R_{ijkl} + 2K_{i[k}K_{l]j}) l^i n^j \overset{*}{m}{}^k n^l \\
&- 4 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}) (l^{[0} n^{j]} \overset{*}{m}{}^k n^l + l^k n^j\overset{*}{m}{}^{[0} n^{l]}) \\
&+ 4 (R_{jl} - K_{jp} K^p_l + KK_{jl}) l^{[0} n^{j]} \overset{*}{m}{}^{[0} n^{l]} \\
\psi_2 =&\ (R_{ijkl} + 2K_{i[k}K_{l]j}) l^i m^j \overset{*}{m}{}^k n^l \\
&- 4 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}) (l^{[0} m^{j]} \overset{*}{m}{}^k n^l + l^k m^l \overset{*}{m}{}^{[0} n^{j]}) \\
&+ 4 (R_{jl} - K_{jp} K^p_l + KK_{jl}) l^{[0} m^{j]} \overset{*}{m}{}^{[0} n^{l]} \\
\psi_1 =&\ (R_{ijkl} + 2K_{i[k}K_{l]j}) n^i l^j m^k l^l \\
&- 4 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}) (n^{[0} l^{j]} m^k l^l + n^k l^l m^{[0} l^{j]}) \\
&+ 4 (R_{jl} - K_{jp} K^p_l + KK_{jl}) n^{[0} l^{j]} m^{[0} l^{l]} \\
\psi_0 =&\ (R_{ijkl} + 2K_{i[k}K_{l]j}) l^i m^j l^k m^l \\
&- 8 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}) l^{[0} m^{j]} l^k m^l \\
&+ 4 (R_{jl} - K_{jp} K^p_l + KK_{jl}) l^{[0} m^{j]} l^{[0} m^{l]}. \\
\end{align}

To make it easier to track the construction of this expression, we will break it down into three parts, by first defining each of the parenthetical terms above separately. This is effectively identical to the procedure used in the Mathematica notebook that generates the original [WeylScal4](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src). That is, let 
\begin{align}
\text{GaussDDDD[i][j][k][l]} =& R_{ijkl} + 2K_{i[k}K_{l]j},
\end{align}


```python
# Step 7: Build the formula for \psi_4.
# Gauss equation: involving the Riemann tensor and extrinsic curvature.
GaussDDDD = ixp.zerorank4()
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                GaussDDDD[i][j][k][l] += RiemannDDDD[i][j][k][l] + kDD[i][k]*kDD[l][j] - kDD[i][l]*kDD[k][j]
```

\begin{align}
\text{CodazziDDD[j][k][l]} =& -2 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p}),
\end{align}                                                       


```python
# Codazzi equation: involving partial derivatives of the extrinsic curvature.
# We will first need to declare derivatives of kDD
kDD_dD = ixp.declarerank3("kDD_dD","sym01")
CodazziDDD = ixp.zerorank3()
for j in range(3):
    for k in range(3):
        for l in range(3):
            CodazziDDD[j][k][l] += kDD_dD[j][l][k] - kDD_dD[j][k][l]

for j in range(3):
    for k in range(3):
        for l in range(3):
            for p in range(3):
                CodazziDDD[j][k][l] += GammaUDD[p][j][l]*kDD[k][p] - GammaUDD[p][j][k]*kDD[l][p]
```

and
\begin{align}
\text{RojoDD[j][l]} = & R_{jl} - K_{jp} K^p_l + KK_{jl} \\
= & \gamma^{pd} R_{jpld} - K_{jp} K^p_l + KK_{jl},
\end{align}


```python
# Another piece. While not associated with any particular equation,
# this is still useful for organizational purposes.
RojoDD = ixp.zerorank2()
for j in range(3):
    for l in range(3):
        RojoDD[j][l] += trK*kDD[j][l]

for j in range(3):
    for l in range(3):
        for p in range(3):
            for d in range(3):
                RojoDD[j][l] += gammaUU[p][d]*RiemannDDDD[j][p][l][d] - kDD[j][p]*gammaUU[p][d]*kDD[d][l]
```

where these quantities are so named because of their relation to the Gauss-Codazzi equations. Then, we simply contract these with the tetrad we chose earlier to arrive at an expression for $\psi_4$. The barred Christoffel symbols and barred Ricci tensor have already been calculated by [BSSN/BSSN_RHSs.py](../edit/BSSN/BSSN_RHSs.py), so we use those values. So, our expression for $\psi_4$ has become 
\begin{align}
\psi_4 =&\ (\text{GaussDDDD[i][j][k][l]}) n^i \overset{*}{m}{}^j n^k \overset{*}{m}{}^l \\
&+2 (\text{CodazziDDD[j][k][l]}) n^{0} \overset{*}{m}{}^{j} n^k \overset{*}{m}{}^l \\
&+ (\text{RojoDD[j][l]}) n^{0} \overset{*}{m}{}^{j} n^{0} \overset{*}{m}{}^{l}.
\end{align}

Likewise, we can rewrite the other Weyl scalars:
\begin{align}
\psi_3 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i n^j \overset{*}{m}{}^k n^l \\
&+ (\text{CodazziDDD[j][k][l]}) (l^{0} n^{j} \overset{*}{m}{}^k n^l - l^{j} n^{0} \overset{*}{m}{}^k n^l - l^k n^j\overset{*}{m}{}^l n^0) \\
&- (\text{RojoDD[j][l]}) l^{0} n^{j} \overset{*}{m}{}^l n^0 - l^{j} n^{0} \overset{*}{m}{}^l n^0 \\
\psi_2 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j \overset{*}{m}{}^k n^l \\
&+ (\text{CodazziDDD[j][k][l]}) (l^{0} m^{j} \overset{*}{m}{}^k n^l - l^{j} m^{0} \overset{*}{m}{}^k n^l - l^k m^l \overset{*}{m}{}^l n^0) \\
&- (\text{RojoDD[j][l]}) l^0 m^j \overset{*}{m}{}^l n^0 \\
\psi_1 =&\ (\text{GaussDDDD[i][j][k][l]}) n^i l^j m^k l^l \\
&+ (\text{CodazziDDD[j][k][l]}) (n^{0} l^{j} m^k l^l - n^{j} l^{0} m^k l^l - n^k l^l m^j l^0) \\
&- (\text{RojoDD[j][l]}) (n^{0} l^{j} m^l l^0 - n^{j} l^{0} m^l l^0) \\
\psi_0 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j l^k m^l \\
&+2 (\text{CodazziDDD[j][k][l]}) (l^0 m^j l^k m^l + l^k m^l l^0 m^j) \\
&+ (\text{RojoDD[j][l]}) l^0 m^j l^0 m^j. \\
\end{align}

We will start by setting the scalars to SymPy's $0$ (this is done so that Python knows that the scalars are symbolic, not numeric, avoiding some potential bugs later on) and then performing the needed contractions of `RojoDD[j][l]`. Recall that the tetrad vectors were defined above and that we just built `RojoDD[j][l]` from the Ricci tensor and extrinsic curvature.

The relevant terms here are:
\begin{align}
\psi_4:&\ \ \ (\text{RojoDD[j][l]}) n^{0} \overset{*}{m}{}^{j} n^{0} \overset{*}{m}{}^{l} \\
\psi_3:&\ -(\text{RojoDD[j][l]}) (l^{0} n^{j} \overset{*}{m}{}^l n^0 - l^{j} n^{0} \overset{*}{m}{}^l n^0) \\
\psi_2:&\ - (\text{RojoDD[j][l]}) l^0 m^j \overset{*}{m}{}^l n^0 \\
\psi_1:&\ - (\text{RojoDD[j][l]}) (n^{0} l^{j} m^l l^0 - n^{j} l^{0} m^l l^0) \\
\psi_0:&\ \ \ (\text{RojoDD[j][l]}) l^0 m^j l^0 m^j. \\
\end{align}


```python
# Now we can calculate $\psi_4$ itself!
psi4r = sp.sympify(0)
psi4i = sp.sympify(0)
psi3r = sp.sympify(0)
psi3i = sp.sympify(0)
psi2r = sp.sympify(0)
psi2i = sp.sympify(0)
psi1r = sp.sympify(0)
psi1i = sp.sympify(0)
psi0r = sp.sympify(0)
psi0i = sp.sympify(0)
for l in range(3):
    for j in range(3):
        psi4r += RojoDD[j][l] * nn * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
        psi4i += RojoDD[j][l] * nn * nn * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
        psi3r +=-RojoDD[j][l] * nn * nn * (ntetU[j]-ltetU[j]) * remtetU[l]
        psi3i += RojoDD[j][l] * nn * nn * (ntetU[j]-ltetU[j]) * immtetU[l]
        psi2r +=-RojoDD[j][l] * nn * nn * (remtetU[l]*remtetU[j]+immtetU[j]*immtetU[l])
        psi2i +=-RojoDD[j][l] * nn * nn * (immtetU[l]*remtetU[j]-remtetU[j]*immtetU[l])
        psi1r += RojoDD[j][l] * nn * nn * (ntetU[j]*remtetU[l]-ltetU[j]*remtetU[l])
        psi1i += RojoDD[j][l] * nn * nn * (ntetU[j]*immtetU[l]-ltetU[j]*immtetU[l])
        psi0r += RojoDD[j][l] * nn * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
        psi0i += RojoDD[j][l] * nn * nn * (remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])
```

Now, we will add the contractions of `CodazziDDD[j][k][l]` to the Weyl Scalars. Again, we use the null tetrad we constructed and the tensor `CodazziDDD[j][k][l]` we constructed from the extrinsic curvature and Christoffel symbols.

The relevant terms here are:
\begin{align}
\psi_4:&\ 2 (\text{CodazziDDD[j][k][l]}) n^{0} \overset{*}{m}{}^{j} n^k \overset{*}{m}{}^l \\
\psi_3:&\ \ \ (\text{CodazziDDD[j][k][l]}) (l^{0} n^{j} \overset{*}{m}{}^k n^l - l^{j} n^{0} \overset{*}{m}{}^k n^l - l^k n^j\overset{*}{m}{}^l n^0) \\
\psi_2:&\ \ \ (\text{CodazziDDD[j][k][l]}) (l^{0} m^{j} \overset{*}{m}{}^k n^l - l^{j} m^{0} \overset{*}{m}{}^k n^l - l^k m^l \overset{*}{m}{}^l n^0) \\
\psi_1:&\ \ \ (\text{CodazziDDD[j][k][l]}) (n^{0} l^{j} m^k l^l - n^{j} l^{0} m^k l^l - n^k l^l m^j l^0) \\
\psi_0:&\ 2 (\text{CodazziDDD[j][k][l]}) (l^0 m^j l^k m^l + l^k m^l l^0 m^j). \\
\end{align}


```python
for l in range(3):
    for j in range(3):
        for k in range(3):
            psi4r += 2 * CodazziDDD[j][k][l] * ntetU[k] * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
            psi4i += 2 * CodazziDDD[j][k][l] * ntetU[k] * nn * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
            psi3r += 1 * CodazziDDD[j][k][l] * nn * ((ntetU[j]-ltetU[j])*remtetU[k]*ntetU[l]-remtetU[j]*ltetU[k]*ntetU[l])
            psi3i +=-1 * CodazziDDD[j][k][l] * nn * ((ntetU[j]-ltetU[j])*immtetU[k]*ntetU[l]-immtetU[j]*ltetU[k]*ntetU[l])
            psi2r += 1 * CodazziDDD[j][k][l] * nn * (ntetU[l]*(remtetU[j]*remtetU[k]+immtetU[j]*immtetU[k])-ltetU[k]*(remtetU[j]*remtetU[l]+immtetU[j]*immtetU[l]))
            psi2i += 1 * CodazziDDD[j][k][l] * nn * (ntetU[l]*(immtetU[j]*remtetU[k]-remtetU[j]*immtetU[k])-ltetU[k]*(remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l]))
            psi1r += 1 * CodazziDDD[j][k][l] * nn * (ltetU[j]*remtetU[k]*ltetU[l]-remtetU[j]*ntetU[k]*ltetU[l]-ntetU[j]*remtetU[k]*ltetU[l])
            psi1i += 1 * CodazziDDD[j][k][l] * nn * (ltetU[j]*immtetU[k]*ltetU[l]-immtetU[j]*ntetU[k]*ltetU[l]-ntetU[j]*immtetU[k]*ltetU[l])
            psi0r += 2 * CodazziDDD[j][k][l] * nn * ltetU[k]*(remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
            psi0i += 2 * CodazziDDD[j][k][l] * nn * ltetU[k]*(remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])
```

Finally, we will add the contractions of `GaussDDDD[i][j][k][l]` (from the Riemann tensor and extrinsic curvature, above) with the null tetrad.

The relevant terms here are:
\begin{align}
\psi_4:&\ (\text{GaussDDDD[i][j][k][l]}) n^i \overset{*}{m}{}^j n^k \overset{*}{m}{}^l \\
\psi_3:&\ (\text{GaussDDDD[i][j][k][l]}) l^i n^j \overset{*}{m}{}^k n^l \\
\psi_2:&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j \overset{*}{m}{}^k n^l \\
\psi_1:&\ (\text{GaussDDDD[i][j][k][l]}) n^i l^j m^k l^l \\
\psi_0:&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j l^k m^l. \\
\end{align}


```python
for l in range(3):
    for j in range(3):
        for k in range(3):
            for i in range(3):
                psi4r += GaussDDDD[i][j][k][l] * ntetU[i] * ntetU[k] * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                psi4i += GaussDDDD[i][j][k][l] * ntetU[i] * ntetU[k] * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
                psi3r += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[j] * remtetU[k] * ntetU[l]
                psi3i +=-GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[j] * immtetU[k] * ntetU[l]
                psi2r += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[l] * (remtetU[j]*remtetU[k]+immtetU[j]*immtetU[k])
                psi2i += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[l] * (immtetU[j]*remtetU[k]-remtetU[j]*immtetU[k])
                psi1r += GaussDDDD[i][j][k][l] * ntetU[i] * ltetU[j] * remtetU[k] * ltetU[l]
                psi1i += GaussDDDD[i][j][k][l] * ntetU[i] * ltetU[j] * immtetU[k] * ltetU[l]
                psi0r += GaussDDDD[i][j][k][l] * ltetU[i] * ltetU[k] * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                psi0i += GaussDDDD[i][j][k][l] * ltetU[i] * ltetU[k] * (remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])
```

<a id='code_validation1'></a>

## Step 5.a: Code Validation against `WeylScal4NRPy.WeylScalars_Cartesian` NRPy+ Module  \[Back to [top](#toc)\]
$$\label{code_validation1}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for Weyl invariants between

1. this tutorial and 
2. the NRPy+ [WeylScal4NRPy.WeylScalars_Cartesian](../edit/WeylScal4NRPy/WeylScalars_Cartesian.py) module.


```python
#psi4rb,psi4ib,psi3rb,psi3ib,psi2rb,psi2ib,psi1rb,psi1ib,psi0rb,psi0ib = psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i
gri.glb_gridfcs_list = []
import WeylScal4NRPy.WeylScalars_Cartesian as weyl
par.set_parval_from_str("WeylScal4NRPy.WeylScalars_Cartesian::output_scalars","all_psis")
weyl.WeylScalars_Cartesian()

print("Consistency check between WeylScalars_Cartesian tutorial and NRPy+ module: ALL SHOULD BE ZERO.")

print("psi4r - weyl.psi4r = " + str(psi4r - weyl.psi4r))
print("psi4i - weyl.psi4i = " + str(psi4i - weyl.psi4i))
print("psi3r - weyl.psi3r = " + str(psi3r - weyl.psi3r))
print("psi3i - weyl.psi3i = " + str(psi3i - weyl.psi3i))
print("psi2r - weyl.psi2r = " + str(psi2r - weyl.psi2r))
print("psi2i - weyl.psi2i = " + str(psi2i - weyl.psi2i))
print("psi1r - weyl.psi1r = " + str(psi1r - weyl.psi1r))
print("psi1i - weyl.psi1i = " + str(psi1i - weyl.psi1i))
print("psi0r - weyl.psi0r = " + str(psi0r - weyl.psi0r))
print("psi0i - weyl.psi0i = " + str(psi0i - weyl.psi0i))
```

    Consistency check between WeylScalars_Cartesian tutorial and NRPy+ module: ALL SHOULD BE ZERO.
    psi4r - weyl.psi4r = 0
    psi4i - weyl.psi4i = 0
    psi3r - weyl.psi3r = 0
    psi3i - weyl.psi3i = 0
    psi2r - weyl.psi2r = 0
    psi2i - weyl.psi2i = 0
    psi1r - weyl.psi1r = 0
    psi1i - weyl.psi1i = 0
    psi0r - weyl.psi0r = 0
    psi0i - weyl.psi0i = 0


<a id='invariant_scalars'></a>

# Step 6: The Invariant Scalars \[Back to [top](#toc)\]
$$\label{invariant_scalars}$$

We may also wish to compute the invariant scalars, whose value does not depend on the choice of the null tetrad. While they are defined using the Weyl tensor, they can also be expressed in terms of the Weyl scalars. We will use those expressions for simplicity.

Following after the method used in the Kranc code, we will read in the already-computed values of the Weyl scalars to find the invariants instead of trying to make NRPy output a very large expression in terms of the metric and extrinsic curvature.

We will start with the invariants $I$ and $J$, as defined in equations (2.3a) and (2.3b) of [arXiv:gr-qc/0407013](https://arxiv.org/abs/gr-qc/0407013). They are as follows.
\begin{align}
I &= 3 \psi_2^2 - 4 \psi_1 \psi_3 + \psi_4 \psi_0 \\
J &= 
\begin{vmatrix}
\psi_4 & \psi_3 & \psi_2 \\
\psi_3 & \psi_2 & \psi_1 \\
\psi_2 & \psi_1 & \psi_0 \\
\end{vmatrix}
\end{align}
Here, since we can work in terms of the Weyl scalars themselves, we will use SymPy's built-in tools for handling complex numbers, which will not get overwhelmed as they did when computing the Weyl scalars.


```python
gri.glb_gridfcs_list = []
psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i = gri.register_gridfunctions("AUX",["psi4r","psi4i",
                                                                                                "psi3r","psi3i",
                                                                                                "psi2r","psi2i",
                                                                                                "psi1r","psi1i",
                                                                                                "psi0r","psi0i"])

psi4 = psi4r + sp.I * psi4i
psi3 = psi3r + sp.I * psi3i
psi2 = psi2r + sp.I * psi2i
psi1 = psi1r + sp.I * psi1i
psi0 = psi0r + sp.I * psi0i

curvIr = sp.re(3*psi2*psi2 - 4*psi1*psi3 + psi4*psi0)
curvIi = sp.im(3*psi2*psi2 - 4*psi1*psi3 + psi4*psi0)
curvJr = sp.re(psi4 * (psi2*psi0 - psi1*psi1) - \
               psi3 * (psi3*psi0 - psi1*psi2) +\
               psi2 * (psi3*psi1 - psi2*psi2) )
curvJi = sp.im(psi4 * (psi2*psi0 - psi1*psi1) - \
               psi3 * (psi3*psi0 - psi1*psi2) +\
               psi2 * (psi3*psi1 - psi2*psi2) )
```

We will now code the invariants $J_1$, $J_2$, $J_3$, and $J_4$, as found in equations B5-B8 of [arXiv:0704.1756](https://arxiv.org/abs/0704.1756). As with the other invariants, we will simply read in the values of the gridfunctions that they already calculated (that is, the Weyl scalars). These equations are based directly on those used in the Mathematica notebook that generates `WeylScal4` (available at [this](https://bitbucket.org/einsteintoolkit/einsteinanalysis/src) repository), modified so that Python can interpret them. Those equations were generated in turn using `xTensor` from equations B5-B8.


```python
J1curv =-16*(3*psi2i**2-3*psi2r**2-4*psi1i*psi3i+4*psi1r*psi3r+psi0i*psi4i-psi0r*psi4r)

J2curv = 96*(-3*psi2i**2*psi2r+psi2r**3+2*psi1r*psi2i*psi3i+2*psi1i*psi2r*psi3i-psi0r*psi3i**2+
            2*psi1i*psi2i*psi3r-2*psi1r*psi2r*psi3r-2*psi0i*psi3i*psi3r+psi0r*psi3r**2-
            2*psi1i*psi1r*psi4i+psi0r*psi2i*psi4i+psi0i*psi2r*psi4i-psi1i**2*psi4r+psi1r**2*psi4r+
            psi0i*psi2i*psi4r-psi0r*psi2r*psi4r)

J3curv = 64*(9*psi2i**4-54*psi2i**2*psi2r**2+9*psi2r**4-24*psi1i*psi2i**2*psi3i+48*psi1r*psi2i*psi2r*psi3i+
            24*psi1i*psi2r**2*psi3i+16*psi1i**2*psi3i**2-16*psi1r**2*psi3i**2+
            24*psi1r*psi2i**2*psi3r+48*psi1i*psi2i*psi2r*psi3r-24*psi1r*psi2r**2*psi3r-64*psi1i*psi1r*psi3i*psi3r-
            16*psi1i**2*psi3r**2+16*psi1r**2*psi3r**2+6*psi0i*psi2i**2*psi4i-12*psi0r*psi2i*psi2r*psi4i-
            6*psi0i*psi2r**2*psi4i-8*psi0i*psi1i*psi3i*psi4i+8*psi0r*psi1r*psi3i*psi4i+8*psi0r*psi1i*psi3r*psi4i+
            8*psi0i*psi1r*psi3r*psi4i+psi0i**2*psi4i**2-psi0r**2*psi4i**2-6*psi0r*psi2i**2*psi4r-
            12*psi0i*psi2i*psi2r*psi4r+6*psi0r*psi2r**2*psi4r+8*psi0r*psi1i*psi3i*psi4r+8*psi0i*psi1r*psi3i*psi4r+
            8*psi0i*psi1i*psi3r*psi4r-8*psi0r*psi1r*psi3r*psi4r-4*psi0i*psi0r*psi4i*psi4r-psi0i**2*psi4r**2+
            psi0r**2*psi4r**2)

J4curv = -640*(-15*psi2i**4*psi2r+30*psi2i**2*psi2r**3-3*psi2r**5+10*psi1r*psi2i**3*psi3i+
              30*psi1i*psi2i**2*psi2r*psi3i-30*psi1r*psi2i*psi2r**2*psi3i-10*psi1i*psi2r**3*psi3i-
              16*psi1i*psi1r*psi2i*psi3i**2-3*psi0r*psi2i**2*psi3i**2-8*psi1i**2*psi2r*psi3i**2+
              8*psi1r**2*psi2r*psi3i**2-6*psi0i*psi2i*psi2r*psi3i**2+3*psi0r*psi2r**2*psi3i**2+
              4*psi0r*psi1i*psi3i**3+4*psi0i*psi1r*psi3i**3+10*psi1i*psi2i**3*psi3r-
              30*psi1r*psi2i**2*psi2r*psi3r-30*psi1i*psi2i*psi2r**2*psi3r+10*psi1r*psi2r**3*psi3r-
              16*psi1i**2*psi2i*psi3i*psi3r+16*psi1r**2*psi2i*psi3i*psi3r-6*psi0i*psi2i**2*psi3i*psi3r+
              32*psi1i*psi1r*psi2r*psi3i*psi3r+12*psi0r*psi2i*psi2r*psi3i*psi3r+6*psi0i*psi2r**2*psi3i*psi3r+
              12*psi0i*psi1i*psi3i**2*psi3r-12*psi0r*psi1r*psi3i**2*psi3r+16*psi1i*psi1r*psi2i*psi3r**2+
              3*psi0r*psi2i**2*psi3r**2+8*psi1i**2*psi2r*psi3r**2-8*psi1r**2*psi2r*psi3r**2+
              6*psi0i*psi2i*psi2r*psi3r**2-3*psi0r*psi2r**2*psi3r**2-12*psi0r*psi1i*psi3i*psi3r**2-
              12*psi0i*psi1r*psi3i*psi3r**2-4*psi0i*psi1i*psi3r**3+4*psi0r*psi1r*psi3r**3-
              6*psi1i*psi1r*psi2i**2*psi4i+2*psi0r*psi2i**3*psi4i-6*psi1i**2*psi2i*psi2r*psi4i+
              6*psi1r**2*psi2i*psi2r*psi4i+6*psi0i*psi2i**2*psi2r*psi4i+6*psi1i*psi1r*psi2r**2*psi4i-
              6*psi0r*psi2i*psi2r**2*psi4i-2*psi0i*psi2r**3*psi4i+12*psi1i**2*psi1r*psi3i*psi4i-
              4*psi1r**3*psi3i*psi4i-2*psi0r*psi1i*psi2i*psi3i*psi4i-2*psi0i*psi1r*psi2i*psi3i*psi4i-
              2*psi0i*psi1i*psi2r*psi3i*psi4i+2*psi0r*psi1r*psi2r*psi3i*psi4i-2*psi0i*psi0r*psi3i**2*psi4i+
              4*psi1i**3*psi3r*psi4i-12*psi1i*psi1r**2*psi3r*psi4i-2*psi0i*psi1i*psi2i*psi3r*psi4i+
              2*psi0r*psi1r*psi2i*psi3r*psi4i+2*psi0r*psi1i*psi2r*psi3r*psi4i+2*psi0i*psi1r*psi2r*psi3r*psi4i-
              2*psi0i**2*psi3i*psi3r*psi4i+2*psi0r**2*psi3i*psi3r*psi4i+2*psi0i*psi0r*psi3r**2*psi4i-
              psi0r*psi1i**2*psi4i**2-2*psi0i*psi1i*psi1r*psi4i**2+psi0r*psi1r**2*psi4i**2+
              2*psi0i*psi0r*psi2i*psi4i**2+psi0i**2*psi2r*psi4i**2-psi0r**2*psi2r*psi4i**2-3*psi1i**2*psi2i**2*psi4r+
              3*psi1r**2*psi2i**2*psi4r+2*psi0i*psi2i**3*psi4r+12*psi1i*psi1r*psi2i*psi2r*psi4r-
              6*psi0r*psi2i**2*psi2r*psi4r+3*psi1i**2*psi2r**2*psi4r-3*psi1r**2*psi2r**2*psi4r-
              6*psi0i*psi2i*psi2r**2*psi4r+2*psi0r*psi2r**3*psi4r+4*psi1i**3*psi3i*psi4r-12*psi1i*psi1r**2*psi3i*psi4r-
              2*psi0i*psi1i*psi2i*psi3i*psi4r+2*psi0r*psi1r*psi2i*psi3i*psi4r+2*psi0r*psi1i*psi2r*psi3i*psi4r+
              2*psi0i*psi1r*psi2r*psi3i*psi4r-psi0i**2*psi3i**2*psi4r+psi0r**2*psi3i**2*psi4r-
              12*psi1i**2*psi1r*psi3r*psi4r+4*psi1r**3*psi3r*psi4r+2*psi0r*psi1i*psi2i*psi3r*psi4r+
              2*psi0i*psi1r*psi2i*psi3r*psi4r+2*psi0i*psi1i*psi2r*psi3r*psi4r-2*psi0r*psi1r*psi2r*psi3r*psi4r+
              4*psi0i*psi0r*psi3i*psi3r*psi4r+psi0i**2*psi3r**2*psi4r-psi0r**2*psi3r**2*psi4r-
              2*psi0i*psi1i**2*psi4i*psi4r+4*psi0r*psi1i*psi1r*psi4i*psi4r+2*psi0i*psi1r**2*psi4i*psi4r+
              2*psi0i**2*psi2i*psi4i*psi4r-2*psi0r**2*psi2i*psi4i*psi4r-4*psi0i*psi0r*psi2r*psi4i*psi4r+
              psi0r*psi1i**2*psi4r**2+2*psi0i*psi1i*psi1r*psi4r**2-psi0r*psi1r**2*psi4r**2-
              2*psi0i*psi0r*psi2i*psi4r**2-psi0i**2*psi2r*psi4r**2+psi0r**2*psi2r*psi4r**2)

#cse_output = sp.cse(psi0i,sp.numbered_symbols("tmp"))
#for commonsubexpression in cse_output:
#    print("hello?",commonsubexpression)
#for commonsubexpression in cse_output[0]:
#    print((str(commonsubexpression[0])+" = "+str(commonsubexpression[1])+";").replace("**","^").replace("_d","d"))
#for i,result in enumerate(cse_output[1]):
#   print(("psi0iPy = "+str(result)+";").replace("**","^").replace("_d","d"))
# These replace commands are used to allow us to validate against Einstein Toolkit's WeylScal4 thorn in Mathematica.
# Specifically, the first changes exponentiation to Mathematica's format, and the second strips the underscores
# that have a very specific meaning in Mathematica and thus cannot be used in variable names.
```

<a id='code_validation2'></a>

## Step 6.a: Code Validation against `WeylScal4NRPy.WeylScalarInvariants_Cartesian` NRPy+ Module  \[Back to [top](#toc)\]
$$\label{code_validation2}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for Weyl invariants between

1. this tutorial and 
2. the NRPy+ [WeylScal4NRPy.WeylScalarInvariants_Cartesian](../edit/WeylScal4NRPy/WeylScalarInvariants_Cartesian.py) module.


```python
# Reset the list of gridfunctions, as registering a gridfunction
#   twice will spawn an error.
gri.glb_gridfcs_list = []

import WeylScal4NRPy.WeylScalarInvariants_Cartesian as invar
invar.WeylScalarInvariants_Cartesian()

num_failures = 0
if curvIr - invar.curvIr != 0: num_failures += 1
if curvIi - invar.curvIi != 0: num_failures += 1
if curvJr - invar.curvJr != 0: num_failures += 1
if curvJi - invar.curvJi != 0: num_failures += 1
if J1curv - invar.J1curv != 0: num_failures += 1
if J2curv - invar.J2curv != 0: num_failures += 1
if J3curv - invar.J3curv != 0: num_failures += 1
if J4curv - invar.J4curv != 0: num_failures += 1

import sys
if num_failures == 0:
    print("ScalarInvariants_Cartesian tutorial vs NRPy+ module: TESTS PASSED.")
else:
    print("ScalarInvariants_Cartesian tutorial vs NRPy+ module: TESTS FAILED, WITH "+str(num_failures)+" FAILURES")
    sys.exit(1)
```

    Consistency check between ScalarInvariants_Cartesian tutorial and NRPy+ module for invariant scalars: PASSED.


<a id='latex_pdf_output'></a>

# Step 7: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-WeylScalarsInvariants-Cartesian.pdf](Tutorial-WeylScalarsInvariants-Cartesian.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-WeylScalarsInvariants-Cartesian")
```

    Created Tutorial-WeylScalarsInvariants-Cartesian.tex, and compiled LaTeX
        file to PDF file Tutorial-WeylScalarsInvariants-Cartesian.pdf

