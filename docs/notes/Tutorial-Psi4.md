<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# The Outgoing Gravitational Wave Weyl scalar $\psi_4$

## Author: Zach Etienne

## This notebook presents the construction of $\psi_4$, a complex scalar for gravitational wave analysis. Using the ADM spatial metric, extrinsic curvature, and arbitrary tetrad vectors, a detailed process is outlined to form $\psi_4$ following [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf).

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated as follows:

* Agreement (to roundoff error) with the WeylScal4 ETK thorn in Cartesian coordinates (as it agrees to roundoff error with Patrick Nelson's [Cartesian Weyl Scalars & Invariants NRPy+ tutorial notebook](Tutorial-WeylScalarsInvariants-Cartesian.ipynb), which itself was validated against WeylScal4). 
* In SinhSpherical coordinates this module yields amplitude falloff and phase agreement with black hole perturbation theory for the $l=2,m=0$ (dominant, spin-weight -2 spherical harmonic) mode of a ringing Brill-Lindquist black hole remnant to more than 7 decades in amplitude, surpassing the agreement seen in Fig. 6 of [Ruchlin, Etienne, & Baumgarte](https://arxiv.org/pdf/1712.07658.pdf). (as shown in [corresponding start-to-finish notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide-Psi4.ipynb); must choose EvolOption = "high resolution").
* Above head-on collision calculation was performed with the [Einstein Toolkit](https://einsteintoolkit.org/), and results were found to match in the case of all spherical harmonic modes.
* An equal-mass binary black hole calculation (QC-0) was performed using this module for $\psi_4$ (SinhSpherical coordinates in a region where gravitational waves extracted), and excellent agreement was observed for all (spin-weight -2 spherical harmonic) modes up to and including $l=4$, when compared with the same calculation performed with the [Einstein Toolkit](https://einsteintoolkit.org/).

### NRPy+ Source Code for this module: [BSSN/Psi4.py](../edit/BSSN/Psi4.py)

## Introduction:
This module constructs $\psi_4$, a quantity that is immensely useful when extracting gravitational wave content from a numerical relativity simulation. $\psi_4$ is related to the gravitational wave strain via

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

Note that $\psi_4$ is complex, with the imaginary components originating from the tetrad vector $m^\mu$. This module does not specify a tetrad; instead, it only constructs the above expression leaving $m^\mu$ and $n^\mu$ unspecified. The [next module on tetrads defines these tetrad quantities](Tutorial-Psi4_tetrads.ipynb) (currently only a quasi-Kinnersley tetrad is supported).

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
1. [Step 2](#riemann): Constructing the 3-Riemann tensor $R_{ik\ell m}$
1. [Step 3](#rank4termone): Constructing the rank-4 tensor in Term 1 of $\psi_4$: $R_{ijkl} + 2 K_{i[k} K_{l]j}$
1. [Step 4](#rank3termtwo): Constructing the rank-3 tensor in Term 2 of $\psi_4$: $-8 \left(K_{j[k,l]} + \Gamma^{p}_{j[k} K_{l]p} \right)$
1. [Step 5](#rank2termthree): Constructing the rank-2 tensor in term 3 of $\psi_4$: $+4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right)$
1. [Step 6](#psifour): Constructing $\psi_4$ through contractions of the above terms with arbitrary tetrad vectors $n^\mu$ and $m^\mu$
1. [Step 7](#code_validation): Code Validation against `BSSN.Psi4` NRPy+ module
1. [Step 8](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from NRPy+:


```python
# Step 1.a: import all needed modules from NRPy+:
import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import reference_metric as rfm

# Step 1.b: Set the coordinate system for the numerical grid
#           Note that this parameter is assumed to be set
#           prior to calling the Python Psi4.py module,
#           so this Step will not appear there.
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

# Step 1.f: Initialize tetrad vectors.
#           mre4U = $\text{Re}(m^\mu)$
#           mim4U = $\text{Im}(m^\mu)$, and
#           n4U   = $n^\mu$
#           Note that in the separate Python Psi4.py
#           module, these will be set to the tetrad
#           chosen within the Psi4_tetrads.py module.

#           We choose the most general form for the
#           tetrad vectors here instead, to ensure complete
#           code validation.
mre4U = ixp.declarerank1("mre4U",DIM=4)
mim4U = ixp.declarerank1("mim4U",DIM=4)
n4U   = ixp.declarerank1("n4U"  ,DIM=4)
```

<a id='riemann'></a>

# Step 2: Constructing the 3-Riemann tensor $R_{ik\ell m}$ \[Back to [top](#toc)\]
$$\label{riemann}$$

Analogously to Christoffel symbols, the Riemann tensor is a measure of the curvature of an $N$-dimensional manifold. Thus the 3-Riemann tensor is not simply a projection of the 4-Riemann tensor (see e.g., Eq. 2.7 of [Campanelli *et al* (1998)](https://arxiv.org/pdf/gr-qc/9803058.pdf) for the relation between 4-Riemann and 3-Riemann), as $N$-dimensional Riemann tensors are meant to define a notion of curvature given only the associated $N$-dimensional metric. 

So, given the ADM 3-metric, the Riemann tensor in arbitrary dimension is given by the 3-dimensional version of Eq. 1.19 in Baumgarte & Shapiro's *Numerical Relativity*. I.e.,

$$
R^i_{jkl} = \partial_k \Gamma^{i}_{jl} - \partial_l \Gamma^{i}_{jk} + \Gamma^i_{mk} \Gamma^m_{jl} - \Gamma^{i}_{ml} \Gamma^{m}_{jk},
$$
where $\Gamma^i_{jk}$ is the Christoffel symbol associated with the 3-metric $\gamma_{ij}$:

$$
\Gamma^l_{ij} = \frac{1}{2} \gamma^{lk} \left(\gamma_{ki,j} + \gamma_{kj,i} - \gamma_{ij,k} \right) 
$$

Notice that this equation for the Riemann tensor is equivalent to the equation given in the Wikipedia article on [Formulas in Riemannian geometry](https://en.wikipedia.org/w/index.php?title=List_of_formulas_in_Riemannian_geometry&oldid=882667524):

$$
R^\ell{}_{ijk}=
\partial_j \Gamma^\ell{}_{ik}-\partial_k\Gamma^\ell{}_{ij}
+\Gamma^\ell{}_{js}\Gamma_{ik}^s-\Gamma^\ell{}_{ks}\Gamma^s{}_{ij},
$$
with the replacements $i\to \ell$, $j\to i$, $k\to j$, $l\to k$, and $s\to m$. Wikipedia also provides a simpler form in terms of second-derivatives of three-metric itself (using the definition of the Christoffel symbol), so that we need not define derivatives of the Christoffel symbol:

$$
R_{ik\ell m}=\frac{1}{2}\left(
\gamma_{im,k\ell} 
+ \gamma_{k\ell,im}
- \gamma_{i\ell,km}
- \gamma_{km,i\ell} \right)
+\gamma_{np} \left(
\Gamma^n{}_{k\ell} \Gamma^p{}_{im} - 
\Gamma^n{}_{km} \Gamma^p{}_{i\ell} \right).
$$

First, we construct the term on the left:


```python
# Step 2: Construct the (rank-4) Riemann curvature tensor associated with the ADM 3-metric:
RDDDD = ixp.zerorank4()
gammaDDdDD = AB.gammaDDdDD

for i in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            for m in range(DIM):
                RDDDD[i][k][l][m] = sp.Rational(1,2) * \
                    (gammaDDdDD[i][m][k][l] + gammaDDdDD[k][l][i][m] - gammaDDdDD[i][l][k][m] - gammaDDdDD[k][m][i][l])
```

... then we add the term on the right:


```python
# ... then we add the term on the right:
gammaDD = AB.gammaDD
GammaUDD = AB.GammaUDD

for i in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            for m in range(DIM):
                for n in range(DIM):
                    for p in range(DIM):
                        RDDDD[i][k][l][m] += gammaDD[n][p] * \
                            (GammaUDD[n][k][l]*GammaUDD[p][i][m] - GammaUDD[n][k][m]*GammaUDD[p][i][l])
```

<a id='rank4termone'></a>

# Step 3: Constructing the rank-4 tensor in Term 1 of $\psi_4$: $R_{ijkl} + 2 K_{i[k} K_{l]j}$ \[Back to [top](#toc)\]
$$\label{rank4termone}$$

Following Eq. 5.1 in [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf), the rank-4 tensor in the first term of $\psi_4$ is given by

$$
R_{ijkl} + 2 K_{i[k} K_{l]j} = R_{ijkl} + K_{ik} K_{lj} - K_{il} K_{kj}
$$


```python
# Step 3: Construct the (rank-4) tensor in term 1 of psi_4 (referring to Eq 5.1 in
#   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
rank4term1DDDD = ixp.zerorank4()
KDD = AB.KDD

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                rank4term1DDDD[i][j][k][l] = RDDDD[i][j][k][l] + KDD[i][k]*KDD[l][j] - KDD[i][l]*KDD[k][j]
```

<a id='rank3termtwo'></a>

# Step 4: Constructing the rank-3 tensor in Term 2 of $\psi_4$: $-8 \left(K_{j[k,l]} + \Gamma^{p}_{j[k} K_{l]p} \right)$ \[Back to [top](#toc)\]
$$\label{rank3termtwo}$$

Following Eq. 5.1 in [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf), the rank-3 tensor in the second term of $\psi_4$ is given by

$$
-8 \left(K_{j[k,l]} + \Gamma^{p}_{j[k} K_{l]p} \right)
$$
First let's construct the first term in this sum: $K_{j[k,l]} = \frac{1}{2} (K_{jk,l} - K_{jl,k})$:


```python
# Step 4: Construct the (rank-3) tensor in term 2 of psi_4 (referring to Eq 5.1 in
#   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
rank3term2DDD = ixp.zerorank3()
KDDdD = AB.KDDdD

for j in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            rank3term2DDD[j][k][l] = sp.Rational(1,2)*(KDDdD[j][k][l] - KDDdD[j][l][k])
```

... then we construct the second term in this sum: $\Gamma^{p}_{j[k} K_{l]p} = \frac{1}{2} (\Gamma^{p}_{jk} K_{lp}-\Gamma^{p}_{jl} K_{kp})$:


```python
# ... then we construct the second term in this sum:
#  \Gamma^{p}_{j[k} K_{l]p} = \frac{1}{2} (\Gamma^{p}_{jk} K_{lp}-\Gamma^{p}_{jl} K_{kp}):
for j in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            for p in range(DIM):
                rank3term2DDD[j][k][l] += sp.Rational(1,2)*(GammaUDD[p][j][k]*KDD[l][p] - GammaUDD[p][j][l]*KDD[k][p])
```

Finally, we multiply the term by $-8$:


```python
# Finally, we multiply the term by $-8$:
for j in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            rank3term2DDD[j][k][l] *= sp.sympify(-8)
```

<a id='rank2termthree'></a>

# Step 5: Constructing the rank-2 tensor in term 3 of $\psi_4$: $+4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right)$ \[Back to [top](#toc)\]
$$\label{rank2termthree}$$

Following Eq. 5.1 in [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf), the rank-2 tensor in the third term of $\psi_4$ is given by

$$
+4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right),
$$
where
\begin{align}
R_{jl} &= R^i_{jil} \\
&= \gamma^{im} R_{ijml} \\
K &= K^i_i \\
&= \gamma^{im} K_{im}
\end{align}

Let's build the components of this term: $R_{jl}$, $K^p_l$, and $K$, as defined above:


```python
# Step 5: Construct the (rank-2) tensor in term 3 of psi_4 (referring to Eq 5.1 in
#   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf

# Step 5.1: Construct 3-Ricci tensor R_{ij} = gamma^{im} R_{ijml}
RDD = ixp.zerorank2()
gammaUU = AB.gammaUU
for j in range(DIM):
    for l in range(DIM):
        for i in range(DIM):
            for m in range(DIM):
                RDD[j][l] += gammaUU[i][m]*RDDDD[i][j][m][l]

# Step 5.2: Construct K^p_l = gamma^{pi} K_{il}
KUD = ixp.zerorank2()
for p in range(DIM):
    for l in range(DIM):
        for i in range(DIM):
            KUD[p][l] += gammaUU[p][i]*KDD[i][l]

# Step 5.3: Construct trK = gamma^{ij} K_{ij}
trK = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        trK += gammaUU[i][j]*KDD[i][j]
```

Next we put these terms together to construct the entire term:
$$
+4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right),
$$


```python
# Next we put these terms together to construct the entire term in parentheses:
# +4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right),
rank2term3DD = ixp.zerorank2()
for j in range(DIM):
    for l in range(DIM):
        rank2term3DD[j][l] = RDD[j][l] + trK*KDD[j][l]
        for p in range(DIM):
            rank2term3DD[j][l] += - KDD[j][p]*KUD[p][l]
# Finally we multiply by +4:
for j in range(DIM):
    for l in range(DIM):
        rank2term3DD[j][l] *= sp.sympify(4)
```

<a id='psifour'></a>

# Step 6: Constructing $\psi_4$ through contractions of the above terms with arbitrary tetrad vectors $m^\mu$ and $n^\mu$ \[Back to [top](#toc)\]
$$\label{psifour}$$

Eq. 5.1 in [Baker, Campanelli, Lousto (2001)](https://arxiv.org/pdf/gr-qc/0104063.pdf) writes $\psi_4$ (which is complex) as the contraction of each of the above terms with products of tetrad vectors:

\begin{align}
\psi_4 &= \left[ {R}_{ijkl}+2K_{i[k}K_{l]j}\right]
{n}^i\bar{m}^j{n}^k\bar{m}^l  \\
& -8\left[ K_{j[k,l]}+{\Gamma }_{j[k}^pK_{l]p}\right]
{n}^{[0}\bar{m}^{j]}{n}^k\bar{m}^l \\
& +4\left[ {R}_{jl}-K_{jp}K_l^p+KK_{jl}\right]
{n}^{[0}\bar{m}^{j]}{n}^{[0}\bar{m}^{l]},
\end{align}
where $\bar{m}^\mu$ is the complex conjugate of $m^\mu$, and $n^\mu$ is real. The third term is given by
\begin{align}
{n}^{[0}\bar{m}^{j]}{n}^{[0}\bar{m}^{l]}
&= \frac{1}{2}({n}^{0}\bar{m}^{j} - {n}^{j}\bar{m}^{0} )\frac{1}{2}({n}^{0}\bar{m}^{l} - {n}^{l}\bar{m}^{0} )\\
&= \frac{1}{4}({n}^{0}\bar{m}^{j} - {n}^{j}\bar{m}^{0} )({n}^{0}\bar{m}^{l} - {n}^{l}\bar{m}^{0} )\\
&= \frac{1}{4}({n}^{0}\bar{m}^{j}{n}^{0}\bar{m}^{l} - {n}^{j}\bar{m}^{0}{n}^{0}\bar{m}^{l} - {n}^{0}\bar{m}^{j}{n}^{l}\bar{m}^{0} +  {n}^{j}\bar{m}^{0}{n}^{l}\bar{m}^{0})
\end{align}

Only $m^\mu$ is complex, so we can separate the real and imaginary parts of $\psi_4$ by hand, defining $M^\mu$ to now be the real part of $\bar{m}^\mu$ and $\mathcal{M}^\mu$ to be the imaginary part. All of the above products are of the form ${n}^\mu\bar{m}^\nu{n}^\eta\bar{m}^\delta$, so let's evaluate the real and imaginary parts of this product once, for all such terms:

\begin{align}
{n}^\mu\bar{m}^\nu{n}^\eta\bar{m}^\delta
&= {n}^\mu(M^\nu - i \mathcal{M}^\nu){n}^\eta(M^\delta - i \mathcal{M}^\delta) \\
&= \left({n}^\mu M^\nu {n}^\eta M^\delta -
{n}^\mu \mathcal{M}^\nu {n}^\eta \mathcal{M}^\delta \right)+
i \left(
-{n}^\mu M^\nu {n}^\eta \mathcal{M}^\delta
-{n}^\mu \mathcal{M}^\nu {n}^\eta M^\delta
\right)
\end{align}


```python
# Step 6: Construct real & imaginary parts of psi_4
#         by contracting constituent rank 2, 3, and 4
#         tensors with input tetrads mre4U, mim4U, & n4U.

def tetrad_product__Real_psi4(n,Mre,Mim,  mu,nu,eta,delta):
    return +n[mu]*Mre[nu]*n[eta]*Mre[delta] - n[mu]*Mim[nu]*n[eta]*Mim[delta]

def tetrad_product__Imag_psi4(n,Mre,Mim,  mu,nu,eta,delta):
    return -n[mu]*Mre[nu]*n[eta]*Mim[delta] - n[mu]*Mim[nu]*n[eta]*Mre[delta]


# We split psi_4 into three pieces, to expedite & possibly parallelize C code generation.
psi4_re_pt = [sp.sympify(0),sp.sympify(0),sp.sympify(0)]
psi4_im_pt = [sp.sympify(0),sp.sympify(0),sp.sympify(0)]

# First term:
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                psi4_re_pt[0] += rank4term1DDDD[i][j][k][l] * \
                                 tetrad_product__Real_psi4(n4U,mre4U,mim4U, i+1,j+1,k+1,l+1)
                psi4_im_pt[0] += rank4term1DDDD[i][j][k][l] * \
                                 tetrad_product__Imag_psi4(n4U,mre4U,mim4U, i+1,j+1,k+1,l+1)

# Second term:
for j in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            psi4_re_pt[1] += rank3term2DDD[j][k][l] * \
                       sp.Rational(1,2)*(+tetrad_product__Real_psi4(n4U,mre4U,mim4U, 0,j+1,k+1,l+1)
                                         -tetrad_product__Real_psi4(n4U,mre4U,mim4U, j+1,0,k+1,l+1) )
            psi4_im_pt[1] += rank3term2DDD[j][k][l] * \
                       sp.Rational(1,2)*(+tetrad_product__Imag_psi4(n4U,mre4U,mim4U, 0,j+1,k+1,l+1)
                                         -tetrad_product__Imag_psi4(n4U,mre4U,mim4U, j+1,0,k+1,l+1) )
# Third term:
for j in range(DIM):
    for l in range(DIM):
        psi4_re_pt[2] += rank2term3DD[j][l] * \
                       (sp.Rational(1,4)*(+tetrad_product__Real_psi4(n4U,mre4U,mim4U, 0,j+1,0,l+1)
                                          -tetrad_product__Real_psi4(n4U,mre4U,mim4U, j+1,0,0,l+1)
                                          -tetrad_product__Real_psi4(n4U,mre4U,mim4U, 0,j+1,l+1,0)
                                          +tetrad_product__Real_psi4(n4U,mre4U,mim4U, j+1,0,l+1,0)))
        psi4_im_pt[2] += rank2term3DD[j][l] * \
                       (sp.Rational(1,4)*(+tetrad_product__Imag_psi4(n4U,mre4U,mim4U, 0,j+1,0,l+1)
                                          -tetrad_product__Imag_psi4(n4U,mre4U,mim4U, j+1,0,0,l+1)
                                          -tetrad_product__Imag_psi4(n4U,mre4U,mim4U, 0,j+1,l+1,0)
                                          +tetrad_product__Imag_psi4(n4U,mre4U,mim4U, j+1,0,l+1,0)))
```

<a id='code_validation'></a>

# Step 7: Code validation against `BSSN.Psi4` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As a code validation check, we verify agreement in the SymPy expressions for the RHSs of the BSSN equations between
1. this tutorial and 
2. the NRPy+ [BSSN.Psi4](../edit/BSSN/Psi4.py) module.

By default, we compare all quantities in Spherical coordinates, though other coordinate systems may be chosen.


```python
# Call the BSSN_RHSs() function from within the
#          BSSN/BSSN_RHSs.py module,
#          which should do exactly the same as in Steps 1-16 above.
import BSSN.Psi4 as BP4
BP4.Psi4(specify_tetrad=False)

print("Consistency check between this tutorial and BSSN.Psi4 NRPy+ module: ALL SHOULD BE ZERO.")

for part in range(3):
    print("psi4_im_pt["+str(part)+"] - BP4.psi4_im_pt["+str(part)+"] = " + str(psi4_im_pt[part] - BP4.psi4_im_pt[part]))
    print("psi4_re_pt["+str(part)+"] - BP4.psi4_re_pt["+str(part)+"] = " + str(psi4_re_pt[part] - BP4.psi4_re_pt[part]))
```

    Consistency check between this tutorial and BSSN.Psi4 NRPy+ module: ALL SHOULD BE ZERO.
    psi4_im_pt[0] - BP4.psi4_im_pt[0] = 0
    psi4_re_pt[0] - BP4.psi4_re_pt[0] = 0
    psi4_im_pt[1] - BP4.psi4_im_pt[1] = 0
    psi4_re_pt[1] - BP4.psi4_re_pt[1] = 0
    psi4_im_pt[2] - BP4.psi4_im_pt[2] = 0
    psi4_re_pt[2] - BP4.psi4_re_pt[2] = 0


<a id='latex_pdf_output'></a>

# Step 8: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Psi4.pdf](Tutorial-Psi4.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Psi4")
```

    Created Tutorial-Psi4.tex, and compiled LaTeX file to PDF file Tutorial-
        Psi4.pdf

