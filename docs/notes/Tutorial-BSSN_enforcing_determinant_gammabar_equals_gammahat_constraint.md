<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Enforce conformal 3-metric $\det{\bar{\gamma}_{ij}}=\det{\hat{\gamma}_{ij}}$ constraint (Eq. 53 of [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658))

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook uses a method for enhancing numerical relativity simulation stability, based on [Brown](https://arxiv.org/abs/0902.3652)'s covariant BSSN formalism. It resolves numerical errors causing temporal divergence in the conformal 3-metric $\bar{\gamma}{ij}$, which disrupts the underlying PDEs' hyperbolicity. A correction mechanism aligns the determinant of $\bar{\gamma}{ij}$ with the reference metric $\hat{\gamma}_{ij}$ at each Runge-Kutta timestep, safeguarding PDE hyperbolicity. 

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). In addition, its output has been 

### NRPy+ Source Code for this module: [BSSN/Enforce_Detgammahat_Constraint.py](../edit/BSSN/Enforce_Detgammahat_Constraint.py)

## Introduction:
[Brown](https://arxiv.org/abs/0902.3652)'s covariant Lagrangian formulation of BSSN, which we adopt, requires that $\partial_t \bar{\gamma} = 0$, where $\bar{\gamma}=\det \bar{\gamma}_{ij}$. Further, all initial data we choose satisfies $\bar{\gamma}=\hat{\gamma}$. 

However, numerical errors will cause $\bar{\gamma}$ to deviate from a constant in time. This actually disrupts the hyperbolicity of the PDEs, so to cure this, we adjust $\bar{\gamma}_{ij}$ at the end of each Runge-Kutta timestep, so that its determinant satisfies $\bar{\gamma}=\hat{\gamma}$ at all times. We adopt the following, rather standard prescription (Eq. 53 of [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658)):

$$
\bar{\gamma}_{ij} \to \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{ij}.
$$

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 1](#initializenrpy): Initialize needed NRPy+ modules
1. [Step 2](#enforcegammaconstraint): Enforce the $\det{\bar{\gamma}_{ij}}=\det{\hat{\gamma}_{ij}}$ constraint
1. [Step 3](#code_validation): Code Validation against `BSSN.Enforce_Detgammahat_Constraint` NRPy+ module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step P1: import all needed modules from NRPy+:
from outputC import nrpyAbs,lhrh,outCfunction # NRPy+: Core C code output module
import finite_difference as fin   # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par    # NRPy+: parameter interface
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import cmdline_helper as cmd      # NRPy+: Multi-platform Python command-line interface
import sympy as sp                # SymPy, Python's core symbolic algebra package
import BSSN.BSSN_quantities as Bq # NRPy+: BSSN quantities
import os,shutil,sys              # Standard Python modules for multiplatform OS-level functions

# Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Then we set the coordinate system for the numerical grid
par.set_parval_from_str("reference_metric::CoordSystem","SinhSpherical")
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.
```

<a id='enforcegammaconstraint'></a>

# Step 2: Enforce the $\det{\bar{\gamma}_{ij}}=\det{\hat{\gamma}_{ij}}$ constraint \[Back to [top](#toc)\]
$$\label{enforcegammaconstraint}$$

Recall that we wish to make the replacement:
$$
\bar{\gamma}_{ij} \to \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{ij}.
$$
Notice the expression on the right is guaranteed to have determinant equal to $\hat{\gamma}$.

$\bar{\gamma}_{ij}$ is not a gridfunction, so we must rewrite the above in terms of $h_{ij}$:
\begin{align}
\left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{ij} &= \bar{\gamma}'_{ij}  \\
&= \hat{\gamma}_{ij} + \varepsilon'_{ij} \\
&= \hat{\gamma}_{ij} + \text{Re[i][j]} h'_{ij} \\
\implies h'_{ij} &= \left[\left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \bar{\gamma}_{ij} - \hat{\gamma}_{ij}\right] / \text{Re[i][j]} \\
&= \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \frac{\bar{\gamma}_{ij}}{\text{Re[i][j]}} - \delta_{ij}\\
&= \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \frac{\hat{\gamma}_{ij} + \text{Re[i][j]} h_{ij}}{\text{Re[i][j]}} - \delta_{ij}\\
&= \left(\frac{\hat{\gamma}}{\bar{\gamma}}\right)^{1/3} \left(\delta_{ij} + h_{ij}\right) - \delta_{ij}
\end{align}

Upon inspection, when expressing $\hat{\gamma}$ SymPy generates expressions like `(xx0)^{4/3} = pow(xx0, 4./3.)`, which can yield $\text{NaN}$s when `xx0 < 0` (i.e., in the `xx0` ghost zones). To prevent this, we know that $\hat{\gamma}\ge 0$ for all reasonable coordinate systems, so we make the replacement $\hat{\gamma}\to |\hat{\gamma}|$ below:


```python
# We will need the h_{ij} quantities defined within BSSN_RHSs
#    below when we enforce the gammahat=gammabar constraint
# Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
#         which we declare here in case they haven't yet been declared elsewhere.

Bq.declare_BSSN_gridfunctions_if_not_declared_already()
hDD = Bq.hDD
Bq.BSSN_basic_tensors()
gammabarDD = Bq.gammabarDD

# First define the Kronecker delta:
KroneckerDeltaDD = ixp.zerorank2()
for i in range(DIM):
    KroneckerDeltaDD[i][i] = sp.sympify(1)

# The detgammabar in BSSN_RHSs is set to detgammahat when BSSN_RHSs::detgbarOverdetghat_equals_one=True (default),
#    so we manually compute it here:
dummygammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(gammabarDD)

# Next apply the constraint enforcement equation above.
hprimeDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        # Using nrpyAbs here, as it directly translates to fabs() without additional SymPy processing.
        #    This acts to simplify the final expression somewhat.
        hprimeDD[i][j] = \
        (nrpyAbs(rfm.detgammahat)/detgammabar)**(sp.Rational(1,3)) * (KroneckerDeltaDD[i][j] + hDD[i][j]) \
        - KroneckerDeltaDD[i][j]
```

<a id='code_validation'></a>

# Step 3: Code Validation against `BSSN.Enforce_Detgammahat_Constraint` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the C code output between

1. this tutorial and 
2. the NRPy+ [BSSN.Enforce_Detgammahat_Constraint](../edit/BSSN/Enforce_Detgammahat_Constraint.py) module.


```python
##########
# Step 1: Generate enforce_detgammahat_constraint() using functions in this tutorial notebook:

Ccodesdir = os.path.join("enforce_detgammahat_constraint")
# First remove C code output directory if it exists
# Courtesy https://stackoverflow.com/questions/303200/how-do-i-remove-delete-a-folder-that-is-not-empty
shutil.rmtree(Ccodesdir, ignore_errors=True)
# Then create a fresh directory
cmd.mkdir(Ccodesdir)

enforce_detg_constraint_vars = [lhrh(lhs=gri.gfaccess("in_gfs","hDD00"),rhs=hprimeDD[0][0]),
                                lhrh(lhs=gri.gfaccess("in_gfs","hDD01"),rhs=hprimeDD[0][1]),
                                lhrh(lhs=gri.gfaccess("in_gfs","hDD02"),rhs=hprimeDD[0][2]),
                                lhrh(lhs=gri.gfaccess("in_gfs","hDD11"),rhs=hprimeDD[1][1]),
                                lhrh(lhs=gri.gfaccess("in_gfs","hDD12"),rhs=hprimeDD[1][2]),
                                lhrh(lhs=gri.gfaccess("in_gfs","hDD22"),rhs=hprimeDD[2][2]) ]

enforce_gammadet_string = fin.FD_outputC("returnstring",enforce_detg_constraint_vars,
                                         params="outCverbose=False,preindent=1,includebraces=False")

desc = "Enforce det(gammabar) = det(gammahat) constraint."
name = "enforce_detgammahat_constraint"
outCfunction(
    outfile=os.path.join(Ccodesdir, name + ".h-validation"), desc=desc, name=name,
    params="const rfm_struct *restrict rfmstruct, const paramstruct *restrict params, REAL *restrict in_gfs",
    body=enforce_gammadet_string,
    loopopts="AllPoints,enable_rfm_precompute")

##########
# Step 2: Generate enforce_detgammahat_constraint() using functions in BSSN.Enforce_Detgammahat_Constraint

gri.glb_gridfcs_list = []

import BSSN.Enforce_Detgammahat_Constraint as EGC
EGC.output_Enforce_Detgammahat_Constraint_Ccode(outdir=Ccodesdir,
                                                exprs=EGC.Enforce_Detgammahat_Constraint_symb_expressions())

import filecmp
for file in [os.path.join(Ccodesdir,"enforce_detgammahat_constraint.h")]:
    if filecmp.cmp(file,file+"-validation") == False:
        print("VALIDATION TEST FAILED on file: "+file+".")
        sys.exit(1)
    else:
        print("Validation test PASSED on file: "+file)
##########
```

    Output C function enforce_detgammahat_constraint() to file enforce_detgammahat_constraint/enforce_detgammahat_constraint.h-validation
    Output C function enforce_detgammahat_constraint() to file enforce_detgammahat_constraint/enforce_detgammahat_constraint.h
    Validation test PASSED on file: enforce_detgammahat_constraint/enforce_detgammahat_constraint.h


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.pdf](Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint")
```

    Created Tutorial-
        BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.tex, and
        compiled LaTeX file to PDF file Tutorial-
        BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.pdf

