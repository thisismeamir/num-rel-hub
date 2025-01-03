<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Transforming BSSN Variables between Two Bases
## Author: Zach Etienne


## This notebook demonstrates the use of NRPy+ for the transformation of BSSN variables from one basis to another. The tutorial provides the procedure of Jacobian transformations on vectors and tensors to migrate rescaled BSSN variables from a source grid to a destination grid. This process is outlined, including code validation against the `BSSN.BSSN_basis_transforms` NRPy+ module.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). In addition, if the same basis & grid are used for both source and destination, source and destination tensors are identical. **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Source Code for this module: [BSSN/BSSN_basis_transforms.py](../edit/BSSN/BSSN_basis_transforms.py)



## Introduction:

Given the rescaled BSSN variables:

$$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\},$$ 

we perform needed Jacobian transformations to all vectors and tensors to migrate to another basis. This is a four-step process:

1. Un-rescale all BSSN variables on source grid
1. Transform source grid basis to Cartesian, using center of source grid as origin
1. Basis transform from Cartesian to destination basis at point on destination grid (`xx0,xx1,xx2`)${}_{\rm dst}$
1. Compute rescaled BSSN quantities in destination basis

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules and BSSN variables; Unrescale all BSSN variables on source grid
1. [Step 2](#srctocart): Transform source grid basis to Cartesian, using center of source grid as origin
1. [Step 3](#carttodst): Basis transform from Cartesian to destination grid basis at point on destination grid (`xx0,xx1,xx2`)${}_{\rm dst}$
1. [Step 4](#rescaleindstbasis): Compute rescaled BSSN quantities in destination basis
1. [Step 5](#code_validation): Code Validation Tests
    1. [Step 5.a](#nrpy_module_validate): Confirm identical output to `BSSN.BSSN_basis_transforms` NRPy+ module
    1. [Step 5.b](#confirm_same_basis_identity): Confirm that if same basis chosen for input and output, at the same points `(xx0,xx1,xx2)`, the output is identical to the input
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules and BSSN variables; Unrescale all BSSN variables \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Here we declare BSSN quantities on the source grid $$\left\{h_{i j},a_{i j},\lambda^{i}, \mathcal{V}^i, \mathcal{B}^i\right\}$$ 

Unrescaling is documented in the [BSSN quantities tutorial notebook](Tutorial-BSSN_quantities.ipynb).


```python
# Step P1: Import needed NRPy+ core modules:
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import NRPy_param_funcs as par    # NRPy+: Parameter interface

# Step P2: Declare inputs
# Input BSSN variables on source ("src") grid:
src_hDD = ixp.declarerank2("src_hDD","sym01")
src_aDD = ixp.declarerank2("src_aDD","sym01")
src_lambdaU = ixp.declarerank1("src_lambdaU")
src_vetU = ixp.declarerank1("src_vetU")
src_betU = ixp.declarerank1("src_betU")

# Source ("src") grid basis and coordinate point (xx0,xx1,xx2)_{src} = src_xx[i]
src_basis = "SinhCylindrical"
src_xx = ixp.declarerank1("src_xx")
# Destination ("dst") grid basis and coordinate point (xx0,xx1,xx2)_{dst} = dst_xx[i]
dst_basis = "SinhSpherical"
dst_xx = ixp.declarerank1("dst_xx")

# Step 1: Unrescale all BSSN variables

par.set_parval_from_str("reference_metric::CoordSystem",src_basis)
rfm.reference_metric()

# STOLEN FROM BSSN/BSSN_quantities.py:
# Step 1.a: gammabarDD and AbarDD:
src_gammabarDD = ixp.zerorank2()
src_AbarDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
        src_gammabarDD[i][j] = src_hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
        # Abar_{ij}      = a_{ij}*ReDD[i][j]
        src_AbarDD[i][j]     = src_aDD[i][j] * rfm.ReDD[i][j]

# Step 1.b: LambdabarU, betaU, and BU:
src_LambdabarU = ixp.zerorank1()
src_betaU = ixp.zerorank1()
src_BU = ixp.zerorank1()
for i in range(3):
    src_LambdabarU[i] = src_lambdaU[i] * rfm.ReU[i]
    src_betaU[i]      =    src_vetU[i] * rfm.ReU[i]
    src_BU[i]         =    src_betU[i] * rfm.ReU[i]
```

<a id='srctocart'></a>

# Step 2: Transform source grid basis to Cartesian, using center of source grid as origin \[Back to [top](#toc)\]
$$\label{srctocart}$$

Within [`reference_metric.py`](../edit/reference_metric.py), the `compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()` function defines Jacobians relative to the center of the source (reference metric) grid, at a point $x^j_{\rm src}=$(`xx0,xx1,xx2`)${}_{\rm src}$ on the source grid:
$$
{\rm Jac\_dUCart\_dDsrcUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm src}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUsrc\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm src}}{\partial x^j_{\rm Cart}},
$$

using NRPy+'s `generic_matrix_inverter3x3()` function. 

In terms of these, the transformation of BSSN tensors from `"reference_metric::CoordSystem"` coordinates to Cartesian may be written:

\begin{align}
\bar{\Lambda}^i_{\rm Cart} &= \frac{\partial x^i_{\rm Cart}}{\partial x^\ell_{\rm src}} \bar{\Lambda}^\ell_{\rm src}\\
\beta^i_{\rm Cart} &= \frac{\partial x^i_{\rm Cart}}{\partial x^\ell_{\rm src}} \beta^\ell_{\rm src}\\
B^i_{\rm Cart} &= \frac{\partial x^i_{\rm Cart}}{\partial x^\ell_{\rm src}} B^\ell_{\rm src}\\
\bar{\gamma}^{\rm Cart}_{ij} &= 
\frac{\partial x^\ell_{\rm src}}{\partial x^i_{\rm Cart}}
\frac{\partial x^m_{\rm src}}{\partial x^j_{\rm Cart}} \bar{\gamma}^{\rm src}_{\ell m}\\
\end{align}

The transformation for vectors is provided by the [`reference_metric.py`](../edit/reference_metric.py) function `basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_vectorU)`, and the transformation for rank-2 covariant tensors is provided by `basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_tensorDD)`, also found within `reference_metric.py`.

After performing the basis transformation to Cartesian, we relabel `(xx0,xx1,xx2)` by `(src_xx0,src_xx1,src_xx2)` to avoid ambiguity.


```python
# Step 2: Transform source grid basis to Cartesian, using center of source grid as origin

# Step 2.a: Construct Jacobian & Inverse Jacobians:
Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

# Step 2.b: Convert basis of all BSSN *vectors* to Cartesian
CartLambdabarU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_LambdabarU)
CartbetaU      = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_betaU)
CartBU         = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_BU)

# Step 2.c: Convert basis of all BSSN *tensors* to Cartesian
CartgammabarDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_gammabarDD)
CartAbarDD     = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_AbarDD)

# Step 2.d: All BSSN tensor/vector quantities are written in terms of
#           rescaled quantities and (xx0,xx1,xx2) on the SOURCE grid.
#           To avoid confusion with (xx0,xx1,xx2) on the DESTINATION grid,
#           we replace (xx0,xx1,xx2) with (src_xx0,src_xx1,src_xx2) here:
for i in range(3):
    for k in range(3):
        CartLambdabarU[i] = CartLambdabarU[i].subs(rfm.xx[k],src_xx[k])
        CartbetaU[i]      =      CartbetaU[i].subs(rfm.xx[k],src_xx[k])
        CartBU[i]         =         CartBU[i].subs(rfm.xx[k],src_xx[k])
for i in range(3):
    for j in range(3):
        for k in range(3):
            CartgammabarDD[i][j] = CartgammabarDD[i][j].subs(rfm.xx[k],src_xx[k])
            CartAbarDD[i][j]     =     CartAbarDD[i][j].subs(rfm.xx[k],src_xx[k])
```

<a id='carttodst'></a>

# Step 3:  Basis transform from Cartesian to destination grid basis at point on destination grid (`xx0,xx1,xx2`)${}_{\rm dst}$ \[Back to [top](#toc)\]
$$\label{carttodst}$$

We define Jacobians relative to the center of the destination grid, at a point $x^j_{\rm dst}=$(`xx0,xx1,xx2`)${}_{\rm dst}$ on the destination grid:
$$
{\rm Jac\_dUCart\_dDdstUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm dst}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUdst\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm dst}}{\partial x^j_{\rm Cart}},
$$

using NRPy+'s `generic_matrix_inverter3x3()` function. In terms of these, the transformation of BSSN tensors from Cartesian to the destination grid's `"reference_metric::CoordSystem"` coordinates may be written:

\begin{align}
\bar{\Lambda}^i_{\rm dst} &= \frac{\partial x^i_{\rm dst}}{\partial x^\ell_{\rm Cart}} \bar{\Lambda}^\ell_{\rm Cart}\\
\beta^i_{\rm dst} &= \frac{\partial x^i_{\rm dst}}{\partial x^\ell_{\rm Cart}} \beta^\ell_{\rm Cart}\\
B^i_{\rm dst} &= \frac{\partial x^i_{\rm dst}}{\partial x^\ell_{\rm Cart}} B^\ell_{\rm Cart}\\
\bar{\gamma}^{\rm dst}_{ij} &= 
\frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm dst}}
\frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm dst}} \bar{\gamma}^{\rm Cart}_{\ell m}\\
\end{align}


```python
# Step 3: Transform BSSN tensors in Cartesian basis to destination grid basis, using center of dest. grid as origin

# Step 3.a: Set up destination grid coordinate system
par.set_parval_from_str("reference_metric::CoordSystem",dst_basis)
rfm.reference_metric()

# Step 3.b: Next construct Jacobian and inverse Jacobian matrices:
Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

# Step 3.c: Convert basis of all BSSN *vectors* from Cartesian to destination basis
dst_LambdabarU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartLambdabarU)
dst_betaU      = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartbetaU)
dst_BU         = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartBU)

# Step 3.d: Convert basis of all BSSN *tensors* from Cartesian to destination basis
dst_gammabarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, CartgammabarDD)
dst_AbarDD     = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, CartAbarDD)
```

<a id='rescaleindstbasis'></a>

# Step 4:  Compute rescaled BSSN quantities in destination basis \[Back to [top](#toc)\]
$$\label{rescaleindstbasis}$$

Rescaling is documented in the [BSSN quantities tutorial notebook](Tutorial-BSSN_quantities.ipynb).


```python
# Step 4: Rescale all BSSN quantities
# BASED ON BSSN/BSSN_quantities.py:
# Step 4.a: hDD and aDD:
dst_hDD = ixp.zerorank2()
dst_aDD = ixp.zerorank2()
for i in range(3):
    for j in range(3):
        # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
        # ==>     h_{ij} = (gammabar_{ij} - gammahat_{ij}) / ReDD[i][j]
        dst_hDD[i][j] = (dst_gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
        # Abar_{ij}      = a_{ij}*ReDD[i][j]
        # ==>     a_{ij} = Abar_{ij}/ReDD[i][j]
        dst_aDD[i][j]     = dst_AbarDD[i][j] / rfm.ReDD[i][j]

# Step 4.b: lambdaU, vetU, and betU:
dst_lambdaU = ixp.zerorank1()
dst_vetU    = ixp.zerorank1()
dst_betU    = ixp.zerorank1()
for i in range(3):
    # Lambdabar^i = \lambda^i * ReU[i]
    # ==>  \lambda^i = Lambdabar^i / ReU[i]
    dst_lambdaU[i] = dst_LambdabarU[i] / rfm.ReU[i]
    dst_vetU[i]    =      dst_betaU[i] / rfm.ReU[i]
    dst_betU[i]    =         dst_BU[i] / rfm.ReU[i]

# Step 4.c: All BSSN tensor/vector quantities are written in terms of
#           rescaled quantities and (xx0,xx1,xx2) on the DESTINATION grid.
#           To avoid confusion with (xx0,xx1,xx2) on the SOURCE grid,
#           we replace (xx0,xx1,xx2) with (dst_xx0,dst_xx1,dst_xx2) here:
for i in range(3):
    for k in range(3):
        dst_lambdaU[i] = dst_lambdaU[i].subs(rfm.xx[k],dst_xx[k])
        dst_vetU[i]    =    dst_vetU[i].subs(rfm.xx[k],dst_xx[k])
        dst_betU[i]    =    dst_betU[i].subs(rfm.xx[k],dst_xx[k])
for i in range(3):
    for j in range(3):
        for k in range(3):
            dst_hDD[i][j] = dst_hDD[i][j].subs(rfm.xx[k],dst_xx[k])
            dst_aDD[i][j] = dst_aDD[i][j].subs(rfm.xx[k],dst_xx[k])
```

<a id='code_validation'></a>

# Step 5: Code Validation Tests \[Back to [top](#toc)\] 
$$\label{code_validation}$$

<a id='nrpy_module_validate'></a>

## Step 5.a: Confirm identical output to `BSSN.BSSN_basis_transforms` NRPy+ module \[Back to [top](#toc)\] 
$$\label{nrpy_module_validate}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for BrillLindquist initial data between
1. this tutorial and 
2. the NRPy+ [BSSN.BSSN_basis_transforms](../edit/BSSN/BSSN_basis_transforms.py) module.

By default, we analyze these expressions in Spherical coordinates, though other coordinate systems may be chosen.


```python
import BSSN.BSSN_basis_transforms as Bbt

# Set up expressions from separate BSSN.BSSN_basis_transforms Python module
Bbt.BSSN_basis_transform(src_basis,src_xx, dst_basis,dst_xx,
                         src_hDD,src_aDD,src_lambdaU,src_vetU,src_betU)

# Define functions for comparisons between this Jupyter notebook & associated Python module
def comp_func(expr1,expr2,basename,prefixname2="Bq."):
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

# Set up expression lists for comparisons between this Jupyter notebook & Python module
for i in range(3):
    namecheck_list.extend([gfnm("dst_lambdaU",i),gfnm("dst_vetU",i),gfnm("dst_betU",i)])
    exprcheck_list.extend([Bbt.dst_lambdaU[i],Bbt.dst_vetU[i],Bbt.dst_betU[i]])
    expr_list.extend([dst_lambdaU[i],dst_vetU[i],dst_betU[i]])
    for j in range(3):
        namecheck_list.extend([gfnm("dst_hDD",i,j),gfnm("dst_aDD",i,j)])
        exprcheck_list.extend([Bbt.dst_hDD[i][j],Bbt.dst_aDD[i][j]])
        expr_list.extend([dst_hDD[i][j],dst_aDD[i][j]])

# Compare all SymPy expressions
num_failures=0
for i in range(len(expr_list)):
    num_failures += comp_func(expr_list[i],exprcheck_list[i],namecheck_list[i])

if num_failures == 0:
    print("ALL TESTS PASSED!")
else:
    print(str(num_failures) + " FAILURES.")
    sys.exit(1)
```

    ALL TESTS PASSED!


<a id='confirm_same_basis_identity'></a>

## Step 5.b: Confirm that if same basis chosen for input and output, at the same points `(xx0,xx1,xx2)`, the output is identical to the input \[Back to [top](#toc)\] 
$$\label{confirm_same_basis_identity}$$

Next we verify that if the same basis is chosen for input and output, at the same points `(xx0,xx1,xx2)`, the results are identical.


```python
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
Bbt.BSSN_basis_transform("Spherical",src_xx, "Spherical",src_xx,
                         src_hDD,src_aDD,src_lambdaU,src_vetU,src_betU)

all_passed = True
for i in range(3):
    if sp.simplify(Bbt.dst_lambdaU[i])-src_lambdaU[i] != 0:
        print("Error in lambdaU["+str(i)+"]: "+str(sp.simplify(Bbt.dst_lambdaU[i])-src_lambdaU[i])+" != 0")
        all_passed = False
    for j in range(3):
        if sp.simplify(Bbt.dst_hDD[i][j])-src_hDD[i][j] != 0:
            print("Error in hDD["+str(i)+"]["+str(j)+"]: "+sp.simplify(Bbt.dst_hDD[i][j])-src_hDD[i][j]+" != 0")
            all_passed = False
if all_passed:
    print("ALL TESTS PASSED!")
```

    ALL TESTS PASSED!


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename [Tutorial-BSSN-basis_transforms.pdf](Tutorial-BSSN-basis_transforms.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN-basis_transforms")
```

    Created Tutorial-BSSN-basis_transforms.tex, and compiled LaTeX file to PDF
        file Tutorial-BSSN-basis_transforms.pdf

