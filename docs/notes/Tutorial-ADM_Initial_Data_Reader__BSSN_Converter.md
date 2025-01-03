<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# BSSN in Curvilinear Coordinates Initial Data Reader
## Author: Zach Etienne

## This notebook introduces a numerical relativity module specializing in initial data reading for curvilinear coordinates and BSSN conversion. The module enables initial data reading, sets quantities for BSSN-based temporal evolution, and converts ADM initial data from spherical or Cartesian to a Cartesian basis, before converting to BSSN variables in the desired reference-metric basis. It also facilitates BSSN vector/tensor transformations from Cartesian to reference-metric basis, including necessary rescaling.

**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Source Code for this module: [BSSN/Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.py](../edit/BSSN/Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.py)

## Introduction:

These C functions take as input the ADM variables

$$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i, B^i\right\}$$

and optionally stress-energy tensor

$$T^{\mu\nu}$$

in the spherical or Cartesian basis. It is up to the user to either provide the needed function for reading in these data or use one provided within NRPy+. The function takes the form

```C
void ID_function(const paramstruct *params, const REAL xCart[3],
                 const ID_persist_struct *restrict ID_persist,
                 initial_data_struct *restrict initial_data)

```

where

* `params` provides NRPy+ input parameters
* `xCart[3]` is the input 3-element array storing the Cartesian coordinate `x,y,z` at which `ID_function()` will set the initial data
* `ID_persist` is an input `struct` storing any data needed by the initial data solver/interpolator. For example if `ID_function()` performs interpolations on pseudospectral grids, `ID_persist` might store the pseudospectral coefficients.
* `initial_data` is the struct that `ID_function()` outputs, storing the ADM variables and optionally the stress-energy tensor in either the spherical or Cartesian basis.

`ID_function()` is called from the main driver routine for reading and setting initial data:

* in the spherical basis, the main driver routine is `initial_data_reader__convert_to_BSSN_from_ADM_spherical()`, and
* in the Cartesian basis, the main driver routine is `initial_data_reader__convert_to_BSSN_from_ADM_Cartesian()`.

These driver routines loop over all gridpoints, and at each gridpoint:

1. `ID_function()` is first called, setting $\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i, B^i\right\}$, and optionally $T^{\mu\nu}$ in the spherical or Cartesian basis.
1. `ADM_SphorCart_to_Cart()` then converts ADM variables from the spherical or Cartesian basis to the Cartesian basis.
1. Next, `ADM_Cart_to_BSSN_Cart()` converts ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.
1. Finally, `initial_data_BSSN_basis_transform_Cartesian_to_rfm()` converts BSSN vectors/tensors from Cartesian to reference-metric basis

The above steps will set all the BSSN quantities needed for full evolutions, at all gridpoints, *except* the *derived* BSSN quantity $\lambda^i$. $\lambda^i$ is set in terms of finite-difference derivatives of other BSSN quantities.

As such, after the other BSSN quantities are set at all gridpoints, the function `initial_data_lambdaU_grid_interior()` is called to compute $\lambda^i$ at all points in the grid *interior*.

**It is up to the user** to then call the inner boundary and extrapolation-outer boundary condition routine from the [Curvilinear boundary conditions driver](Tutorial-Start_to_Finish-Curvilinear_BCs_new_way.ipynb), so that all BSSN quantities are set at all points on the numerical grid, *including* the grid boundaries.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#id_reader_driver): Overview of `initial_data_reader__convert_ADM_..._to_BSSN()`: Driver routine that computes/reads ADM variables at all points on all grids, and converts them to BSSN quantities in chosen curvilinear reference metric
1. [Step 2.a](#example_id_func): Example `ID_function()` generator, for "exact" initial data (i.e., initial data in which all quantities are specified with closed-form expressions)
1. [Step 2.b](#adm_sphorcart_to_cart): `ADM_SphorCart_to_Cart()`: Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis
1. [Step 2.c](#adm_cart_to_bssn_cart) `ADM_Cart_to_BSSN_Cart()`: Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis
1. [Step 2.d](#bssn_cart_to_rescaled_bssn_rfm): `BSSN_Cart_to_rescaled_BSSN_rfm()`: Convert Cartesian-basis BSSN vectors/tensors *except* $\lambda^i$, to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities.
1. [Step 2.e](#lambda): `initial_data_lambdaU_grid_interior()`: Compute $\lambda^i$ from finite-difference derivatives of rescaled metric quantities
1. [Step 2.f](#register_driver_function): Putting it all together: Register `initial_data_reader__convert_ADM_..._to_BSSN()` C function 

1. [Step 3](#nbd): `register_NRPy_basic_defines()`: Register `ID_data_struct` within `NRPy_basic_defines.h`
1. [Step 4](#code_validation):Code Validation against  `BSSN.ADM_Initial_Data_Reader__BSSN_Converter` NRPy+ module
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1: Initialize core Python/NRPy+ modules
from outputC import outputC,lhrh,add_to_Cfunction_dict, Cfunction # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import finite_difference as fin   # NRPy+: Finite difference C code generation module
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
import cmdline_helper as cmd      # NRPy+: Multi-platform Python command-line interface
from pickling import pickle_NRPy_env # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import os, sys                    # Standard Python modules for multiplatform OS-level functions
```


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    Cell In[1], line 2
          1 # Step 1: Initialize core Python/NRPy+ modules
    ----> 2 from outputC import outputC,lhrh,add_to_Cfunction_dict, Cfunction # NRPy+: Core C code output module
          3 from outputC import outC_NRPy_basic_defines_h_dict
          4 import NRPy_param_funcs as par    # NRPy+: Parameter interface


    File ~/Documents/VS_Code/Etienne/nrpytutorial/outputC.py:18
         12 __all__ = ['lhrh', 'outCparams', 'nrpyAbs', 'superfast_uniq', 'check_if_string__error_if_not',
         13            'outputC', 'parse_outCparams_string',
         14            'outC_NRPy_basic_defines_h_dict',
         15            'outC_function_prototype_dict', 'outC_function_dict', 'Cfunction', 'add_to_Cfunction_dict', 'outCfunction']
         17 import loop as lp                             # NRPy+: C code loop interface
    ---> 18 import NRPy_param_funcs as par                # NRPy+: parameter interface
         19 from SIMD import expr_convert_to_SIMD_intrins # NRPy+: SymPy expression => SIMD intrinsics interface
         20 from cse_helpers import cse_preprocess,cse_postprocess  # NRPy+: CSE preprocessing and postprocessing


    File ~/Documents/VS_Code/Etienne/nrpytutorial/NRPy_param_funcs.py:10
          1 # As documented in the NRPy+ tutorial module
          2 #   Tutorial-Coutput__Parameter_Interface.ipynb
          3 #   this core NRPy+ module is used for
       (...)
          7 # Author: Zachariah B. Etienne
          8 #         zachetie **at** gmail **dot* com
    ---> 10 import sympy as sp                   # Import SymPy
         11 import sys                           # Standard Python: OS-independent system functions
         12 from collections import namedtuple   # Standard Python: Enable namedtuple data type


    ModuleNotFoundError: No module named 'sympy'


<a id='id_reader_driver'></a>

# Step 2: Overview of `initial_data_reader__convert_ADM_..._to_BSSN()`: Driver routine that computes/reads ADM variables at all points on all grids, and converts them to BSSN quantities in chosen curvilinear reference metric \[Back to [top](#toc)\]
$$\label{id_reader_driver}$$

Initial data for constructing spacetimes in NRPy+ are provided in either the spherical or Cartesian basis, and the corresponding reader routines are `initial_data_reader__convert_ADM_spherical_to_BSSN()` or `initial_data_reader__convert_ADM_Cartesian_to_BSSN()`, respectively. These routines do the following:

1. Call `ID_function()`: Read or compute ADM initial data at Cartesian point `xCart[3]`=$(x,y,z)$, in the spherical basis inside `initial_data_reader__convert_to_BSSN_from_ADM_Spherical()` or  Cartesian basis inside `initial_data_reader__convert_to_BSSN_from_ADM_Cartesian()`.
1. Call `ADM_SphorCart_to_Cart()`: Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis.
1. Call `ADM_Cart_to_BSSN_Cart()`: Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.
1. Call `BSSN_Cart_to_rescaled_BSSN_rfm()`: Convert all BSSN vectors/tensors *except* $\lambda^i$ in the Cartesian basis, to the basis specified by `reference_metric::CoordSystem`. Rescale BSSN quantities.
1. Call `initial_data_lambdaU_grid_interior()` to compute $\lambda^i$ in the `reference_metric::CoordSystem` basis, *in the grid interior only*.

**Important Note**: After `initial_data_reader__convert_to_BSSN_rfm_from_ADM_sph_or_Cart()` is called, inner/outer boundary conditions must be applied to $\lambda^i$ to ensure it is specified on the grid boundaries.

<a id='example_id_func'></a>

## Step 2.a: Example `ID_function()` generator, for "exact" initial data (i.e., initial data in which all quantities are specified with closed-form expressions) \[Back to [top](#toc)\]
$$\label{example_id_func}$$

This function converts closed-form SymPy expressions for ADM quantities $\alpha$, $\beta^i$, $B^i$, $\gamma_{ij}$ and $K_{ij}$ in the spherical or Cartesian basis at a given point $(x,y,z)$ to optimized C code. It is meant to be passed as the `ID_function()` argument into `initial_data_reader__convert_ADM_..._to_BSSN()`. By not calling it `ID_function()`, we can easily have multiple initial data functions within the same C executable.

Regarding the input quantities, a number of initial data sets are provided within NRPy+'s `BSSN` module, including for example:

* [`BSSN.UIUCBlackHole`](BSSN/UIUCBlackHole.py) for spinning single black hole initial data, in which the coordinate size of the BH does not shrink to zero in the limit of maximum spin
* [`BSSN.ShiftedKerrSchild`](BSSN/ShiftedKerrSchild.py) for spinning single black hole initial data
* [`BSSN.StaticTrumpet`](BSSN/StaticTrumpet.py) for static trumpet nonspinning black hole initial data
* [`BSSN.BrillLindquist`](BSSN/BrillLindquist.py) for initial data of two nonspinning black holes at rest

When any of these initial data set functions is called, it will export the ADM quantities as global variables. For example

```python
import BSSN.UIUCBlackHole as UIB
UIB.UIUCBlackHole()
pickled_NRPy_env = \
    add_to_Cfunction_dict_exact_ADM_ID_function("UIUCBlackHole", "Spherical", UIB.alphaSph, UIB.betaSphU,
                                                UIB.BSphU, UIB.betaSphU, UIB.gammaSphDD, UIB.KSphDD)
```

where the returned pickled NRPy+ environment `pickled_NRPy_env` can be used to run this function in parallel with other functions (parallelized codegen).


```python
def add_to_Cfunction_dict_exact_ADM_ID_function(IDtype, IDCoordSystem, alpha, betaU, BU, gammaDD, KDD):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = IDtype + " initial data"
    c_type = "void"
    name = IDtype
    params = "const paramstruct *params, const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"
    desired_rfm_coord = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", IDCoordSystem)
    rfm.reference_metric()
    body = ""
    if IDCoordSystem == "Spherical":
        body += r"""
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
  const REAL r  = xx0; // Some ID only specify r,th,ph.
  const REAL th = xx1;
  const REAL ph = xx2;
"""
    elif IDCoordSystem == "Cartesian":
        body += r"""  const REAL Cartxyz0=xCart[0], Cartxyz1=xCart[1], Cartxyz2=xCart[2];
"""
    else:
        print("add_to_Cfunction_dict_exact_ADM_ID_function() Error: IDCoordSystem == " + IDCoordSystem + " unsupported")
        sys.exit(1)
    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["initial_data->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += ["initial_data->betaSphorCartU" + str(i)]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += ["initial_data->BSphorCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += ["initial_data->gammaSphorCartDD" + str(i) + str(j)]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += ["initial_data->KSphorCartDD" + str(i) + str(j)]
    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore CoordSystem:
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_coord)
    rfm.reference_metric()
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc, c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True)
    return pickle_NRPy_env()
```

<a id='adm_sphorcart_to_cart'></a>

## Step 2.b: `ADM_SphorCart_to_Cart()`: Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis \[Back to [top](#toc)\]
$$\label{adm_sphorcart_to_cart}$$


```python
def Cfunction_ADM_SphorCart_to_Cart(input_Coord="Spherical", include_T4UU=False):
    includes = []

    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    desc = "Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis"
    c_type = "static void"
    name = "ADM_SphorCart_to_Cart"
    params = """paramstruct *restrict params, const REAL xCart[3], const initial_data_struct *restrict initial_data,
                                  ADM_Cart_basis_struct *restrict ADM_Cart_basis"""

    body = r"""
  // Unpack initial_data for ADM vectors/tensors
"""
    for i in ["betaSphorCartU", "BSphorCartU"]:
        for j in range(3):
            varname = i + str(j)
            body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
        body += "\n"
    for i in ["gammaSphorCartDD", "KSphorCartDD"]:
        for j in range(3):
            for k in range(j, 3):
                varname = i + str(j) + str(k)
                body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
        body += "\n"
    # Read stress-energy tensor in spherical or Cartesian basis if desired.
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                varname = "T4SphorCartUU" + str(mu) + str(nu)
                body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
        body += "\n"

    # Set reference_metric to the input_Coord
    par.set_parval_from_str("reference_metric::CoordSystem", input_Coord)
    rfm.reference_metric()

    body += r"""
  // Perform the basis transform on ADM vectors/tensors from """ + input_Coord + """ to Cartesian:

  // Set destination xx[3] based on desired xCart[3]
  REAL xx0,xx1,xx2;
  """ + outputC(rfm.Cart_to_xx[0:3], ["xx0", "xx1", "xx2"],
             filename='returnstring', params="includebraces=True").\
                                      replace("Cartx","xCart[0]").\
                                      replace("Carty","xCart[1]").\
                                      replace("Cartz","xCart[2]")

    # Define the input variables:
    gammaSphorCartDD = ixp.declarerank2("gammaSphorCartDD", "sym01")
    KSphorCartDD     = ixp.declarerank2("KSphorCartDD", "sym01")
    betaSphorCartU = ixp.declarerank1("betaSphorCartU")
    BSphorCartU    = ixp.declarerank1("BSphorCartU")
    T4SphorCartUU = ixp.declarerank2("T4SphorCartUU", "sym01", DIM=4)

    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammaCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, gammaSphorCartDD)
    KCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, KSphorCartDD)
    betaCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, betaSphorCartU)
    BCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, BSphorCartU)
    T4CartUU = ixp.zerorank2(DIM=4)
    if include_T4UU:
        T4CartUU = rfm.basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD,
                                                                                       T4SphorCartUU)

    alpha = sp.symbols("initial_data->alpha", real=True)
    list_of_output_exprs    = [alpha]
    list_of_output_varnames = ["ADM_Cart_basis->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaCartU[i]]
        list_of_output_varnames += ["ADM_Cart_basis->betaU" + str(i)]
        list_of_output_exprs += [BCartU[i]]
        list_of_output_varnames += ["ADM_Cart_basis->BU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaCartDD[i][j]]
            list_of_output_varnames += ["ADM_Cart_basis->gammaDD" + str(i) + str(j)]
            list_of_output_exprs += [KCartDD[i][j]]
            list_of_output_varnames += ["ADM_Cart_basis->KDD" + str(i) + str(j)]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += ["ADM_Cart_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore reference metric globals to coordsystem on grid.
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_basis)
    rfm.reference_metric()

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return func
```

<a id='adm_cart_to_bssn_cart'></a>

## Step 2.c: `ADM_Cart_to_BSSN_Cart()`: Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis \[Back to [top](#toc)\]
$$\label{adm_cart_to_bssn_cart}$$


```python
def Cfunction_ADM_Cart_to_BSSN_Cart(include_T4UU=False):
    includes = []

    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    desc = "Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis"
    c_type = "static void"
    name = "ADM_Cart_to_BSSN_Cart"
    params = """paramstruct *restrict params, const REAL xCart[3], const ADM_Cart_basis_struct *restrict ADM_Cart_basis,
                                  BSSN_Cart_basis_struct *restrict BSSN_Cart_basis"""

    # Extract desired rfm basis from reference_metric::CoordSystem
    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    # Set CoordSystem to Cartesian
    par.set_parval_from_str("reference_metric::CoordSystem", "Cartesian")
    rfm.reference_metric()

    gammaCartDD = ixp.declarerank2("ADM_Cart_basis->gammaDD", "sym01")
    KCartDD     = ixp.declarerank2("ADM_Cart_basis->KDD", "sym01")

    import BSSN.BSSN_in_terms_of_ADM as BitoA
    BitoA.trK_AbarDD_aDD(gammaCartDD, KCartDD)
    BitoA.gammabarDD_hDD(gammaCartDD)
    BitoA.cf_from_gammaDD(gammaCartDD)

    body = r"""
  // *In the Cartesian basis*, convert ADM quantities gammaDD & KDD
  //   into BSSN gammabarDD, AbarDD, cf, and trK.
  BSSN_Cart_basis->alpha = ADM_Cart_basis->alpha;
  BSSN_Cart_basis->betaU0 = ADM_Cart_basis->betaU0;
  BSSN_Cart_basis->betaU1 = ADM_Cart_basis->betaU1;
  BSSN_Cart_basis->betaU2 = ADM_Cart_basis->betaU2;
  BSSN_Cart_basis->BU0 = ADM_Cart_basis->BU0;
  BSSN_Cart_basis->BU1 = ADM_Cart_basis->BU1;
  BSSN_Cart_basis->BU2 = ADM_Cart_basis->BU2;
"""
    list_of_output_exprs    = [BitoA.cf, BitoA.trK]
    list_of_output_varnames = ["BSSN_Cart_basis->cf", "BSSN_Cart_basis->trK"]
    for i in range(3):
        for j in range(i, 3):
            list_of_output_exprs += [BitoA.gammabarDD[i][j]]
            list_of_output_varnames += ["BSSN_Cart_basis->gammabarDD" + str(i) + str(j)]
            list_of_output_exprs += [BitoA.AbarDD[i][j]]
            list_of_output_varnames += ["BSSN_Cart_basis->AbarDD" + str(i) + str(j)]
    if include_T4UU:
        T4CartUU = ixp.declarerank2("ADM_Cart_basis->T4UU", "sym01", DIM=4)
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += ["BSSN_Cart_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))
    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore reference metric globals to desired reference metric.
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_basis)
    rfm.reference_metric()

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return func
```

<a id='bssn_cart_to_rescaled_bssn_rfm'></a>

## Step 2.d: `BSSN_Cart_to_rescaled_BSSN_rfm()`: Convert Cartesian-basis BSSN vectors/tensors *except* $\lambda^i$, to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities. \[Back to [top](#toc)\]
$$\label{bssn_cart_to_rescaled_bssn_rfm}$$

By the time this function is called, all BSSN tensors and vectors are in the Cartesian coordinate basis $x^i_{\rm Cart} = (x,y,z)$, but we need them in the curvilinear coordinate basis $x^i_{\rm rfm}=$`(xx0,xx1,xx2)` set by the `"reference_metric::CoordSystem"` variable. Empirically speaking, it is far easier to write `(x(xx0,xx1,xx2),y(xx0,xx1,xx2),z(xx0,xx1,xx2))` than the inverse, so we will compute the Jacobian matrix

$$
{\rm Jac\_dUCart\_dDrfmUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm rfm}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUrfm\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm rfm}}{\partial x^j_{\rm Cart}},
$$

using NRPy+'s `generic_matrix_inverter3x3()` function. In terms of these, the transformation of BSSN tensors from Spherical to `"reference_metric::CoordSystem"` coordinates may be written:

\begin{align}
\beta^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Cart}} \beta^\ell_{\rm Cart}\\
B^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Cart}} B^\ell_{\rm Cart}\\
\bar{\gamma}^{\rm rfm}_{ij} &= 
\frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm rfm}}
\frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm rfm}} \bar{\gamma}^{\rm Cart}_{\ell m}\\
\bar{A}^{\rm rfm}_{ij} &= 
\frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm rfm}}
\frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm rfm}} \bar{A}^{\rm Cart}_{\ell m}
\end{align}

The above basis transforms are included in functions `basis_transform_tensorDD_from_Cartesian_to_rfmbasis()` and `basis_transform_vectorU_from_Cartesian_to_rfmbasis()` in `reference_metric.py`, and we use them below.

After the basis transform has been performed, we perform tensor rescalings to compute the evolved variables $h_{ij}$, $a_{ij}$, $\text{vet}^i$, and $\text{bet}^i$:

$\bar{\gamma}_{ij}$ is rescaled $h_{ij}$ according to the prescription described in the [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

$$
h_{ij} = (\bar{\gamma}_{ij} - \hat{\gamma}_{ij})/\text{ReDD[i][j]}.
$$

Further $\bar{A}_{ij}$, $\beta^i$, $B^i$ are rescaled to $a_{ij}$, $\text{vet}^i$, and $\text{bet}^i$ respectively via the standard formulas (found in [the covariant BSSN formulation tutorial](Tutorial-BSSN_formulation.ipynb); also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

\begin{align}
a_{ij} &= \bar{A}_{ij}/\text{ReDD[i][j]} \\
\text{vet}^i &= \beta^i/\text{ReU[i]} \\
\text{bet}^i &= B^i/\text{ReU[i]}.
\end{align}

Finally, functions depending on the stress-energy tensor $T^{\mu\nu}$ all assume it to be *unrescaled*. Thus we do not rescale $T^{\mu\nu}$.


```python
def Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(rel_path_to_Cparams=os.path.join("."), include_T4UU=False):
    includes = []

    desc = r"""Convert Cartesian-basis BSSN vectors/tensors *except* lambda^i,
to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities"""
    c_type = "static void"
    name = "BSSN_Cart_to_rescaled_BSSN_rfm"
    params = """paramstruct *restrict params, const REAL xCart[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis"""

    body = r"""
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
"""

    # Define the input variables:
    gammabarCartDD = ixp.declarerank2("BSSN_Cart_basis->gammabarDD", "sym01")
    AbarCartDD     = ixp.declarerank2("BSSN_Cart_basis->AbarDD", "sym01")
    betaCartU = ixp.declarerank1("BSSN_Cart_basis->betaU")
    BCartU    = ixp.declarerank1("BSSN_Cart_basis->BU")

    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammabarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, gammabarCartDD)
    AbarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, AbarCartDD)
    betaU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, betaCartU)
    BU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, BCartU)

    # Next rescale:
    vetU = ixp.zerorank1()
    betU = ixp.zerorank1()
    hDD  = ixp.zerorank2()
    aDD  = ixp.zerorank2()
    for i in range(3):
        vetU[i] = betaU[i] / rfm.ReU[i]
        betU[i] =    BU[i] / rfm.ReU[i]
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]

    if include_T4UU:
        T4CartUU = ixp.declarerank2("BSSN_Cart_basis->T4UU", "sym01", DIM=4)
        T4UU = rfm.basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(Jac_dUrfm_dDCartUD, T4CartUU)

    alpha, cf, trK = sp.symbols('BSSN_Cart_basis->alpha BSSN_Cart_basis->cf BSSN_Cart_basis->trK', real=True)

    list_of_output_exprs    = [alpha, cf, trK]
    list_of_output_varnames = ["rescaled_BSSN_rfm_basis->alpha",
                               "rescaled_BSSN_rfm_basis->cf",
                               "rescaled_BSSN_rfm_basis->trK"]
    for i in range(3):
        list_of_output_exprs += [vetU[i]]
        list_of_output_varnames += ["rescaled_BSSN_rfm_basis->vetU" + str(i)]
        list_of_output_exprs += [betU[i]]
        list_of_output_varnames += ["rescaled_BSSN_rfm_basis->betU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [hDD[i][j]]
            list_of_output_varnames += ["rescaled_BSSN_rfm_basis->hDD" + str(i) + str(j)]
            list_of_output_exprs += [aDD[i][j]]
            list_of_output_varnames += ["rescaled_BSSN_rfm_basis->aDD" + str(i) + str(j)]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                # T4UU IS ASSUMED NOT RESCALED; RESCALINGS ARE HANDLED WITHIN BSSN RHSs, etc.
                list_of_output_exprs += [T4UU[mu][nu]]
                list_of_output_varnames += ["rescaled_BSSN_rfm_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True, rel_path_to_Cparams=rel_path_to_Cparams)
    return func
```

<a id='lambda'></a>

## Step 2.e: `initial_data_lambdaU_grid_interior()`: Compute $\lambda^i$ from finite-difference derivatives of rescaled metric quantities \[Back to [top](#toc)\]
$$\label{lambda}$$

We compute $\bar{\Lambda}^i$ (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)), from finite-difference derivatives of rescaled metric quantities $h_{ij}$:

$$
\bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).
$$

The [reference_metric.py](../edit/reference_metric.py) module provides us with analytic expressions for $\hat{\Gamma}^i_{jk}$, so here we need only compute finite-difference expressions for $\bar{\Gamma}^i_{jk}$, based on the values for $h_{ij}$ provided in the initial data. Once $\bar{\Lambda}^i$ has been computed, we apply the usual rescaling procedure:
$$
\lambda^i = \bar{\Lambda}^i/\text{ReU[i]},
$$
and then output the result to a C file using the NRPy+ finite-difference C output routine.


```python
# initial_data_lambdaU_grid_interior() computes lambdaU from
#  finite-difference derivatives of rescaled metric quantities
def Cfunction_initial_data_lambdaU_grid_interior():
    includes = []
    c_type = "static void"

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Compute lambdaU in " + output_Coord + " coordinates"
    name = "initial_data_lambdaU_grid_interior"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"""
    # Step 7: Compute $\bar{\Lambda}^i$ from finite-difference derivatives of rescaled metric quantities

    # We will need all BSSN gridfunctions to be defined, as well as
    #     expressions for gammabarDD_dD in terms of exact derivatives of
    #     the rescaling matrix and finite-difference derivatives of
    #     hDD's. This functionality is provided by BSSN.BSSN_unrescaled_and_barred_vars,
    #     which we call here to overwrite above definitions of gammabarDD,gammabarUU, etc.
    Bq.gammabar__inverse_and_derivs()  # Provides gammabarUU and GammabarUDD
    gammabarUU    = Bq.gammabarUU
    GammabarUDD   = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU0"), rhs=lambdaU[0]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU1"), rhs=lambdaU[1]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU2"), rhs=lambdaU[2])]
    body = fin.FD_outputC("returnstring", lambdaU_expressions,
                           params="outCverbose=False,includebraces=False,preindent=0")
    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        loopopts="InteriorPoints,Read_xxs",
        enableCparameters=True)
    return func
```

<a id='register_driver_function'></a>

## Step 2.f: Putting it all together: Register `initial_data_reader__convert_ADM_..._to_BSSN()` C function \[Back to [top](#toc)\]
$$\label{register_driver_function}$$

As discussed above, the function `initial_data_reader__convert_ADM_..._to_BSSN()` is our core driver function, generally called directly by `main()` to read/construct initial data in terms of ADM quantities in the spherical/Cartesian basis, and convert these quantities to rescaled BSSN variables. As such, it incorporates all the above functions just described.


```python
def add_to_Cfunction_dict_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(addl_includes=None,
                                                                               rel_path_to_Cparams=os.path.join("."),
                                                                               input_Coord="Spherical",
                                                                               include_T4UU=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    if par.parval_from_str("finite_difference::enable_FD_functions"):
        includes += ["finite_difference_functions.h"]
    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            print("Error: addl_includes must be a list.")
            sys.exit(1)
        includes += addl_includes

    def T4UU_prettyprint():
        return r"""
  REAL T4UU00,T4UU01,T4UU02,T4UU03;
  REAL        T4UU11,T4UU12,T4UU13;
  REAL               T4UU22,T4UU23;
  REAL                      T4UU33;
"""
    prefunc = """
// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL gammaDD00,gammaDD01,gammaDD02,gammaDD11,gammaDD12,gammaDD22;
  REAL KDD00,KDD01,KDD02,KDD11,KDD12,KDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} ADM_Cart_basis_struct;\n"
    ##############
    prefunc += """
// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL cf, trK;
  REAL gammabarDD00,gammabarDD01,gammabarDD02,gammabarDD11,gammabarDD12,gammabarDD22;
  REAL AbarDD00,AbarDD01,AbarDD02,AbarDD11,AbarDD12,AbarDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} BSSN_Cart_basis_struct;\n"
    ##############
    prefunc += """
// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0,vetU1,vetU2, betU0,betU1,betU2;
  REAL cf, trK;
  REAL hDD00,hDD01,hDD02,hDD11,hDD12,hDD22;
  REAL aDD00,aDD01,aDD02,aDD11,aDD12,aDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} rescaled_BSSN_rfm_basis_struct;\n"
    ##############
    ##############
    prefunc += Cfunction_ADM_SphorCart_to_Cart(input_Coord=input_Coord, include_T4UU=include_T4UU)
    prefunc += Cfunction_ADM_Cart_to_BSSN_Cart(                         include_T4UU=include_T4UU)
    prefunc += Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(rel_path_to_Cparams=rel_path_to_Cparams,
                                                        include_T4UU=include_T4UU)
    prefunc += Cfunction_initial_data_lambdaU_grid_interior()

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Read in ADM initial data in the " + input_Coord + " basis, and convert to BSSN data in " + output_Coord + " coordinates"
    c_type = "void"
    name = "initial_data_reader__convert_ADM_" + input_Coord + "_to_BSSN"
    params = """griddata_struct *restrict griddata, ID_persist_struct *restrict ID_persist,
                                                             void ID_function(const paramstruct *params, const REAL xCart[3],
                                                                              const ID_persist_struct *restrict ID_persist,
                                                                              initial_data_struct *restrict initial_data)"""

    body = r"""
  const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0,0,Nxx_plus_2NGHOSTS0, i1,0,Nxx_plus_2NGHOSTS1, i2,0,Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];  xx_to_Cart(&griddata->params, griddata->xx, i0,i1,i2, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(&griddata->params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(&griddata->params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(&griddata->params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(&griddata->params, xCart, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3S(i0,i1,i2);

    // Output data to gridfunctions
"""
    gf_list = ["alpha", "trK", "cf"]
    for i in range(3):
        gf_list += ["vetU"+str(i), "betU"+str(i)]
        for j in range(i, 3):
            gf_list += ["hDD"+str(i)+str(j), "aDD"+str(i)+str(j)]
    for gf in sorted(gf_list):
        body += "    griddata->gridfuncs.y_n_gfs[IDX4ptS("+gf.upper()+"GF, idx3)] = rescaled_BSSN_rfm_basis."+gf+";\n"
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                gf = "T4UU" + str(mu) + str(nu)
                body += "    griddata->gridfuncs.auxevol_gfs[IDX4ptS("+gf.upper()+"GF, idx3)] = rescaled_BSSN_rfm_basis."+gf+";\n"
    body += """
  } // END LOOP over all gridpoints on given grid

  initial_data_lambdaU_grid_interior(&griddata->params, griddata->xx, griddata->gridfuncs.y_n_gfs);
"""

    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return pickle_NRPy_env()
```

<a id='nbd'></a>

# Step 3: `register_NRPy_basic_defines()`: Register `ID_data_struct` within `NRPy_basic_defines.h` \[Back to [top](#toc)\]
$$\label{nbd}$$

Other than its core use as a means to store ADM input quantities, `ID_data_struct` is designed to be extensible. For example, it may be used to store e.g., pseudospectral coefficients for TwoPunctures, initial data gridfunctions from NRPyElliptic, pointers to TOV 1D data from the TOV solver, etc.


```python
# Other than its core use as a means to store ADM input quantities,
# `initial_data_struct` is designed to be extensible. For example, it may be
# used to store e.g., pseudospectral coefficients for TwoPunctures,
# initial data gridfunctions from NRPyElliptic, pointers to TOV 1D data
# from the TOV solver, etc.
def register_NRPy_basic_defines(ID_persist_struct_contents_str="", include_T4UU=False):
    Nbd = r"""typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
"""
    if include_T4UU:
        Nbd += """
  REAL T4SphorCartUU00,T4SphorCartUU01,T4SphorCartUU02,T4SphorCartUU03;
  REAL                 T4SphorCartUU11,T4SphorCartUU12,T4SphorCartUU13;
  REAL                                 T4SphorCartUU22,T4SphorCartUU23;
  REAL                                                 T4SphorCartUU33;
"""
    Nbd += """
} initial_data_struct;
"""

    Nbd += "typedef struct __ID_persist_struct__ {\n"
    Nbd += ID_persist_struct_contents_str + "\n"
    Nbd += "} ID_persist_struct;\n"
    outC_NRPy_basic_defines_h_dict["BSSN_initial_data"] = Nbd
```

<a id='code_validation'></a>

# Step 4: Code Validation against  `BSSN.ADM_Initial_Data_Reader__BSSN_Converter` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify that functions within

1. this tutorial and 
2. the NRPy+ [BSSN.ADM_Initial_Data_Reader__BSSN_Converter](../edit/BSSN/ADM_Initial_Data_Reader__BSSN_Converter.py) module

are *identical*. This notebook serves as documentation for the Python module, and the Python module is meant to be called from external routines. Thus this validation test ensures the documentation remains consistent with the Python module.


```python
import BSSN.ADM_Initial_Data_Reader__BSSN_Converter as AID

funclist = [(add_to_Cfunction_dict_exact_ADM_ID_function, AID.add_to_Cfunction_dict_exact_ADM_ID_function),
            (Cfunction_ADM_SphorCart_to_Cart, AID.Cfunction_ADM_SphorCart_to_Cart),
            (Cfunction_ADM_Cart_to_BSSN_Cart, AID.Cfunction_ADM_Cart_to_BSSN_Cart),
            (Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm, AID.Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm),
            (Cfunction_initial_data_lambdaU_grid_interior, AID.Cfunction_initial_data_lambdaU_grid_interior),
            (add_to_Cfunction_dict_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN, AID.add_to_Cfunction_dict_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN),
            (register_NRPy_basic_defines, AID.register_NRPy_basic_defines)
           ]

import inspect
for func in funclist:
    # https://stackoverflow.com/questions/20059011/check-if-two-python-functions-are-equal
    if inspect.getsource(func[0]) != inspect.getsource(func[1]):
        with open(func[0].__name__ + "_Jupyter_notebook_version.c", "w") as file:
            file.write(inspect.getsource(func[0]))
        with open(func[1].__name__ + "_Python_module_version.c", "w") as file:
            file.write(inspect.getsource(func[1]))
        print("ERROR: function " + func[0].__name__ + " is not the same as the Ccodegen_library version!")
        print(" For more info, try this:")
        print("diff " + func[0].__name__ + "_Jupyter_notebook_version.c" + " " + func[1].__name__ + "_Python_module_version.c")
        sys.exit(1)

print("PASS! ALL FUNCTIONS ARE IDENTICAL")
```

    PASS! ALL FUNCTIONS ARE IDENTICAL


<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

Once the following code finishes running, the generated PDF may be found at the following location within the directory you have the NRPy+ tutorial saved:
[Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.pdf](Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.pdf)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data_Reader__BSSN_Converter")
```

    Created Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.tex, and compiled
        LaTeX file to PDF file Tutorial-
        ADM_Initial_Data_Reader__BSSN_Converter.pdf

