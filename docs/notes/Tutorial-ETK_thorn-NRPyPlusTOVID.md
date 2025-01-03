<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# `NRPyPlusTOVID`: An Einstein Toolkit Thorn for Piecewise-Polytrope TOV neutron star initial data

## Author: Zach Etienne and Leo Werneck

## This notebook introduces the ETK module `NRPyPlusTOVID`, used for generating initial data for equilibrium neutron stars. The module employs NRPy+ for conversion from SymPy expressions to C-code kernel. It outlines methods for transforming this initial data to conform to the Toolkit's standards, including conversion to spherical and Cartesian ADMBase and HydroBase variables.

**Notebook Status:** <font color='orange'><b> Partially Validated </b></font>

**Validation Notes:** NRPy+ TOV initial data generation module validated against [Josh Faber's TOV initial data solver](https://ccrg.rit.edu/~jfaber/BNSID/TOV/), as described in the [NRPy+ implementation notes of the TOV solution for piecewise-polytrope neutron stars](Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb).

### NRPy+ Source Code for this module: [TOV/TOV_Solver.py](../edit/TOV/TOV_Solver.py) [\[tutorial\]](Tutorial-Tutorial-ADM_Initial_Data-TOV.ipynb) Constructs numerical solution to TOV equations for neutron stars with piecewise polytrope equations of state

## Introduction:
In this part of the tutorial, we will construct an Einstein Toolkit (ETK) thorn (module) that will set up [TOV initial data](https://en.wikipedia.org/wiki/Tolman–Oppenheimer–Volkoff_equation) for an equilibrium neutron star. As documented in the [Piecewise Polytrope NRPy+ tutorial](Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb), piecewise-polytrope equations of state are supported, which closely approximate realistic nuclear equations of state appropriate for neutron star matter. In the [Tutorial-Tutorial-ADM_Initial_Data-TOV](Tutorial-Tutorial-ADM_Initial_Data-TOV.ipynb) tutorial notebook, we used NRPy+ to construct the SymPy expressions for these initial data. 

We will construct this thorn in two steps.

1. Call on NRPy+ to convert the SymPy expressions for the initial data into one C-code kernel.
1. Write the C code and linkages to the Einstein Toolkit infrastructure (i.e., the .ccl files) to complete this Einstein Toolkit module.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): **Call on NRPy+ to generate the TOV solution given a piecewise-polytrope equation of state; output the data to a text file**
1. [Step 2](#initial_data): **Converting TOV initial data so that it can be used by the Einstein Toolkit**
    1. [Step 2.a](#initial_data__interpolation): Interpolate the TOV data file as needed
    1. [Step 2.b](#initial_data__tov_to_adm_sph): Converting the TOV variables to ADM variables in Spherical coordinates
    1. [Step 2.c](#initial_data__convert_adm_sph_to_admbase): Convert Spherical ADM quantities to `ADMBase` (Cartesian) variables $\left\{\alpha,\beta^i,\gamma_{ij},K_{ij}\right\}$
    1. [Step 2.d](#initial_data__convert_to_hydrobase): Convert TOV solution quantities to `HydroBase` variables $\left\{P,\rho_{\rm baryonic},\epsilon,v_{\rm (n)}^i\right\}$
1. [Step 3](#einstein): **Interfacing with the Einstein Toolkit**
    1. [Step 3.a](#einstein_c): Constructing the Einstein Toolkit C-code calling functions that include the C code kernels
    1. [Step 3.b](#einstein_ccl): CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure
    1. [Step 3.c](#einstein_list): Add the C code to the Einstein Toolkit compilation list
1. [Step 4](#latex_pdf_output): **Output this notebook to $\LaTeX$-formatted PDF**

<a id='initializenrpy'></a>

# Step 1: Call on NRPy+ to generate the TOV solution given a piecewise-polytrope equation of state; output the data to a text file \[Back to [top](#toc)\]
$$\label{initializenrpy}$$




```python
# Step 1: Import needed core NRPy+ modules
from outputC import lhrh         # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import shutil, os                # Standard Python modules for multiplatform OS-level functions

# Create directory for NRPyPlusTOVID thorn & subdirectories in case they don't exist.
outrootdir = "NRPyPlusTOVID/"
cmd.mkdir(os.path.join(outrootdir))
outdir = os.path.join(outrootdir,"src") # Main C code output directory
cmd.mkdir(outdir)

# Step 1.a: This is an Einstein Toolkit (ETK) thorn. Here we
#           tell NRPy+ that gridfunction memory access will
#           therefore be in the "ETK" style.
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

# Step 1.b: NRPyPlusTOVID uses Cartesian coordinates, so
#           we specify the reference metric to be Cartesian here:
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# ADJUST THIS PARAMETER IF DESIRED.
# "Single" = Single Polytrope
# "APR4"   = APR4 Piecewise Polytrope
# "SLY"    = SLy  Piecewise Polytrope
# .---------------------------------------------.
# | For all available names please look in the  |
# | TOV/Piecewise_Polytrope__dict.py NRPy+ file |
# .---------------------------------------------.
# vvvvvvvvvvvvvvvv
EOSname = "Single"
# EOSname = "SLy"
# EOSname = "APR4"
# ^^^^^^^^^^^^^^^^

# Import our TOV solver, which supports both single
# and piecewise polytropic EOSs
import TOV.TOV_Solver as TOV
import TOV.Polytropic_EOSs as poly

if EOSname=="Single":
    # Set neos = 1 (single polytrope)
    neos = 1

    # Set rho_poly_tab (not needed for a single polytrope)
    rho_poly_tab = []

    # Set Gamma_poly_tab
    Gamma_poly_tab = [2.0]

    # Set K_poly_tab0
    K_poly_tab0 = 1. # ZACH NOTES: CHANGED FROM 100.

    rhob_central = 0.129285309 # M/R_Schw = 1.468770268913230e-01

    # Set the eos quantities
    eos = poly.set_up_EOS_parameters__complete_set_of_input_variables(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab0)
    import time
    start = time.time()
    TOV.TOV_Solver(eos,
                   outfile="outputTOVpolytrope.txt",
                   rho_baryon_central=rhob_central,
                   verbose = True)
    print("Single Polytrope TOV solution generated in: "+str(time.time()-start)+" s")
    print("Initial data file: outputTOVpolytrope.txt")
else:
    # Set up the EOS parameters
    eos = poly.set_up_EOS_parameters__Read_et_al_input_variables(EOSname)

    # Set up the initial condition for the pressure by
    # selecting a central baryon density
#     rhob_central = 2.0 # M/R_Schw = 3.303692404611947e-01
#     rhob_central = 1.0 # M/R_Schw = 2.051637558540178e-01
    rhob_central = 0.8 # M/R_Schw = 1.470662481999595e-01

    # Solve the TOV equations given our EOS and central density
    import time
    start = time.time()
    outfilename = "outputTOVpolytrope-"+EOSname+".txt"
    TOV.TOV_Solver(eos,outfile=outfilename,rho_baryon_central=rhob_central,verbose=True)
    print("PPEOS "+EOSname+" TOV solution generated in: "+str(time.time()-start)+" s")
    print("Initial data file: "+outfilename)
```

    1256 1256 1256 1256 1256 1256
    Just generated a TOV star with
    * M        = 1.405031497682765e-01 ,
    * R_Schw   = 9.566039886703138e-01 ,
    * R_iso    = 8.100079558158075e-01 ,
    * M/R_Schw = 1.468770268913230e-01 
    
    Single Polytrope TOV solution generated in: 0.19237208366394043 s
    Initial data file: outputTOVpolytrope.txt


<a id='initial_data'></a>

# Step 2: Converting TOV initial data so that it can be used by the Einstein Toolkit \[Back to [top](#toc)\]
$$\label{initial_data}$$

Main driver function:

* Looping over all gridpoints:
    * Read in `const CCTK_REAL rr = r[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)];`
    * **Given this radius call interpolation driver to get all the base TOV quantities**
    * **Convert TOV spacetime quantities to ADM quantities in *spherical* basis**
    * Call the Cartesian ADMBase converter
    * Call the HydroBase converter

<a id='initial_data__interpolation'></a>

## Step 2.a: Interpolate the TOV data file as needed \[Back to [top](#toc)\]
$$\label{initial_data__interpolation}$$

We start by interpolating the TOV data file to the gridpoints used by ETK, using the [tov_interp.h](../edit/TOV/tov_interp.h) file, which using Lagrange polynomial interpolation (for more details on the usage of this interpolation file, please look at the [start-to-finish TOV initial data tutorial notebook](Tutorial-Start_to_Finish-BSSNCurvilinear-TOV_initial_data.ipynb)).

Keep in mind that the TOV data file just written stored $\left(r,\rho(r),\rho_{\text{baryonic}}(r),P(r),M(r),e^{\nu(r)}\right)$, where $\rho(r)$ is the total mass-energy density (cf. $\rho_{\text{baryonic}}$).


```python
shutil.copy(os.path.join("TOV","tov_interp.h"),outdir)

with open(os.path.join(outdir,"interpolate_TOV_solution_to_point.h"), "w") as file:
    file.write("""
/* Load the TOV_interpolate_1D() function */
#include "tov_interp.h"

/* This function returns the TOV quantities at point rr
 * by interpolating the data in the TOV initial data file.
 */
void interpolate_TOV_solution_to_point(const CCTK_REAL rr, ID_inputs other_inputs,
                           CCTK_REAL *exp_4phi, CCTK_REAL *expnu,
                           CCTK_REAL *Pressure, CCTK_REAL *rho_baryon, CCTK_REAL *rho__total_energy_density) {

  /* The mass valus is not used, but we have to
   * store it in this dummy variable because the
   * initial data file contains it.
   */
  CCTK_REAL M;

  /* Perform the interpolation, returning:
   *  - rho__total_energy_density
   *  - rho_baryon
   *  - Pressure
   *  - Mass (dummy variable, unused)
   *  - exp(nu)
   *  - exp(4phi)
   */
  TOV_interpolate_1D(rr,other_inputs.Rbar,other_inputs.Rbar_idx,other_inputs.interp_stencil_size,
                     other_inputs.numlines_in_file,
                     other_inputs.r_Schw_arr,other_inputs.rho_arr,other_inputs.rho_baryon_arr,other_inputs.P_arr,other_inputs.M_arr,
                     other_inputs.expnu_arr,other_inputs.exp4phi_arr,other_inputs.rbar_arr,
                     rho__total_energy_density,rho_baryon,Pressure,&M,expnu,exp_4phi);

}\n""")
```

<a id='initial_data__tov_to_adm_sph'></a>

## Step 2.b: Converting the TOV variables to ADM variables in Spherical coordinates \[Back to [top](#toc)\]
$$\label{initial_data__tov_to_adm_sph}$$

Now we perform the interpolation of the TOV quantities to ADM quantities in spherical coordinates, using (see [the TOV initial data tutorial notebook](Tutorial-ADM_Initial_Data-TOV.ipynb) for more details):


\begin{equation}
\boxed{
\begin{aligned}
\alpha &= e^{\nu(\bar{r})/2} \\
\beta^k &= 0 \\
\gamma_{\bar{r}\bar{r}} &= e^{4\phi}\\
\gamma_{\theta\theta} &= e^{4\phi} \bar{r}^2 \\
\gamma_{\phi\phi} &= e^{4\phi} \bar{r}^2 \sin^2 \theta \\
\end{aligned}
}
\end{equation}


```python
with open(os.path.join(outdir,"convert_TOV_spacetime_vars_to_ADM_vars.h"), "w") as file:
    file.write("""
/* This function converts TOV quantities into
 * ADM quantities in Spherical coordinates.
 */
void convert_TOV_spacetime_vars_to_ADM_vars( const CCTK_REAL rr, const CCTK_REAL th,
                           const CCTK_REAL IDexp_4phi, const CCTK_REAL IDexpnu,
                           CCTK_REAL *IDalpha,
                           CCTK_REAL *IDgammaDD00, CCTK_REAL *IDgammaDD01, CCTK_REAL *IDgammaDD02,
                           CCTK_REAL *IDgammaDD11, CCTK_REAL *IDgammaDD12, CCTK_REAL *IDgammaDD22) {

  /***************************************************************
   * Convert TOV quantities to ADM quantities in Spherical basis *
   ***************************************************************
   *
   * First we convert the lapse function:
   * .------------------.
   * | alpha = e^(nu/2) |
   * .------------------.
   */
   *IDalpha = sqrt(IDexpnu);

  /* Next we convert the metric function:
   * .----------------------------------------.
   * | gamma_{00} = e^{4phi}                  |
   * .----------------------------------------.
   * | gamma_{11} = e^{4phi} r^2              |
   * .----------------------------------------.
   * | gamma_{22} = e^{4phi} r^2 sin^2(theta) |
   * .----------------------------------------.
   * | All other components are zero.         |
   * .----------------------------------------.
   */
   *IDgammaDD00 = IDexp_4phi;
   *IDgammaDD11 = IDexp_4phi * rr * rr;
   *IDgammaDD22 = IDexp_4phi * rr * rr * sin(th) * sin(th);
   *IDgammaDD01 = 0.0;
   *IDgammaDD02 = 0.0;
   *IDgammaDD12 = 0.0;

}\n""")
```

<a id='initial_data__convert_adm_sph_to_admbase'></a>

## Step 2.c: Convert Spherical ADM quantities to `ADMBase` (Cartesian) variables $\left\{\alpha,\beta^i,\gamma_{ij},K_{ij}\right\}$ \[Back to [top](#toc)\]
$$\label{initial_data__convert_adm_sph_to_admbase}$$

The [TOV line element](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) in *Schwarzschild coordinates* is written (in the $-+++$ form):
$$
ds^2 = - c^2 e^\nu dt^2 + \left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 + r^2 d\Omega^2.
$$

In *isotropic coordinates* with $G=c=1$ (i.e., the initial coordinate slicing and units we prefer to use), the ($-+++$ form) line element is written:
$$
ds^2 = - e^{\nu} dt^2 + e^{4\phi} \left(d\bar{r}^2 + \bar{r}^2 d\Omega^2\right),
$$
where $\phi$ here is the *conformal factor*.

The ADM 3+1 line element for this diagonal metric in isotropic spherical coordinates is given by:
$$
ds^2 = (-\alpha^2 + \beta_k \beta^k) dt^2 + \gamma_{\bar{r}\bar{r}} d\bar{r}^2 + \gamma_{\theta\theta} d\theta^2+ \gamma_{\phi\phi} d\phi^2,
$$

from which we can immediately read off the ADM quantities:
\begin{align}
\alpha &= e^{\nu(\bar{r})/2} \\
\beta^k &= 0 \\
\gamma_{\bar{r}\bar{r}} &= e^{4\phi}\\
\gamma_{\theta\theta} &= e^{4\phi} \bar{r}^2 \\
\gamma_{\phi\phi} &= e^{4\phi} \bar{r}^2 \sin^2 \theta \\
\end{align}


```python
thismodule   = __name__

IDalpha   = par.Cparameters("REAL", thismodule, "IDalpha", 1e300) # IDalpha must be set in C
IDbetaU   = ixp.zerorank1() # beta^i is zero
IDgammaDD = ixp.zerorank2()
for i in range(3):
    for j in range(i,3):
        IDgammaDD[i][j] = par.Cparameters("REAL", thismodule, "IDgammaDD"+str(i)+str(j), 1e300) # IDgammaDD must be set in C
        IDgammaDD[j][i] = IDgammaDD[i][j]
IDKDD     = ixp.zerorank2() # K_{ij} is zero
```

As this ETK module expects Cartesian coordinates, and the TOV solution above is in the spherical basis, we next perform the Jacobian transformations necessary to convert into the Cartesian basis:

All ADM tensors and vectors are in the Spherical coordinate basis $x^i_{\rm Sph} = (r,\theta,\phi)$, but we need them in the Cartesian coordinate basis $x^i_{\rm Cart}=$`(xx0,xx1,xx2)` set by the `"reference_metric::CoordSystem"` variable. Empirically speaking, it is far easier to write `(x(xx0,xx1,xx2),y(xx0,xx1, xx2),z(xx0,xx1,xx2))` than the inverse, so we will compute the Jacobian matrix

$$
{\rm Jac\_dUSph\_dDrfmUD[i][j]} = \frac{\partial x^i_{\rm Sph}}{\partial x^j_{\rm Cart}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUrfm\_dDSphUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm Sph}},
$$

using NRPy+'s `generic_matrix_inverter3x3()` function. 

In terms of these, the transformation of ADM tensors from Spherical to `"reference_metric::CoordSystem==Cartesian"` coordinates may be written:

\begin{align}
\gamma^{\rm Cart}_{ij} &= 
\frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm Sph}}
\frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm Sph}} \gamma^{\rm Sph}_{\ell m}
\end{align}

Since $\beta^i=K_{ij}=0$ in this case, and $\alpha$ is not a tensor, only the above Jacobian transformation need be performed:


```python
# Transform initial data to our coordinate system:
# First compute Jacobian and its inverse
drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff(rfm.xxSph[0],rfm.xx[1]), sp.diff(rfm.xxSph[0],rfm.xx[2])],
                                       [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                       [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv()

# Declare as gridfunctions the final quantities we will output for the initial data
alpha   = gri.register_gridfunctions("EVOL","alpha")
betaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")
gammaDD = ixp.register_gridfunctions_for_single_rank2("EVOL","gammaDD","sym01")
KDD     = ixp.register_gridfunctions_for_single_rank2("EVOL","KDD","sym01")

alpha = IDalpha # No Jacobian necessary!
betaU = IDbetaU # Because beta^i = 0
KDD   = IDKDD   # Because K_{ij} = 0

for i in range(3):
    for j in range(3):
        # Matrices are stored in row, column format, so (i,j) <-> (row,column)
        gammaDD[i][j] = 0
        for k in range(3):
            for l in range(3):
                gammaDD[i][j] += drrefmetric__dx_0UDmatrix[(k,i)]*drrefmetric__dx_0UDmatrix[(l,j)]*IDgammaDD[k][l]

# -={ Spacetime quantities: Generate C code from expressions and output to file }=-
ADMQuantities_to_print = [\
                          lhrh(lhs=gri.gfaccess("out_gfs","alpha"),rhs=alpha),\
                          lhrh(lhs=gri.gfaccess("out_gfs","betaU0"),rhs=betaU[0]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","betaU1"),rhs=betaU[1]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","betaU2"),rhs=betaU[2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD00"),rhs=gammaDD[0][0]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD01"),rhs=gammaDD[0][1]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD02"),rhs=gammaDD[0][2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD11"),rhs=gammaDD[1][1]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD12"),rhs=gammaDD[1][2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","gammaDD22"),rhs=gammaDD[2][2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD00"),rhs=KDD[0][0]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD01"),rhs=KDD[0][1]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD02"),rhs=KDD[0][2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD11"),rhs=KDD[1][1]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD12"),rhs=KDD[1][2]),\
                          lhrh(lhs=gri.gfaccess("out_gfs","KDD22"),rhs=KDD[2][2])
                          ]

with open(os.path.join(outdir,"ADMQuantities.h"),"w") as file:
    ADMQuantities_CcodeKernel = fin.FD_outputC("returnstring",ADMQuantities_to_print,
                                               params="outCverbose=False,includebraces=False,preindent=1")
    file.write("""
static inline
void ADMQuantities(const cGH* restrict const cctkGH, const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,

                   const CCTK_REAL *restrict xx0GF,const CCTK_REAL *restrict xx1GF,const CCTK_REAL *restrict xx2GF,

                   const CCTK_REAL IDalpha,
                   const CCTK_REAL IDgammaDD00,const CCTK_REAL IDgammaDD01, const CCTK_REAL IDgammaDD02,
                   const CCTK_REAL IDgammaDD11,const CCTK_REAL IDgammaDD12, const CCTK_REAL IDgammaDD22,

                   CCTK_REAL *alphaGF,CCTK_REAL *betaU0GF,CCTK_REAL *betaU1GF,CCTK_REAL *betaU2GF,
                   CCTK_REAL *gammaDD00GF, CCTK_REAL *gammaDD01GF, CCTK_REAL *gammaDD02GF,
                   CCTK_REAL *gammaDD11GF, CCTK_REAL *gammaDD12GF, CCTK_REAL *gammaDD22GF,
                   CCTK_REAL *KDD00GF, CCTK_REAL *KDD01GF, CCTK_REAL *KDD02GF,
                   CCTK_REAL *KDD11GF, CCTK_REAL *KDD12GF, CCTK_REAL *KDD22GF) {
    const CCTK_REAL xx0 = xx0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];
    const CCTK_REAL xx1 = xx1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];
    const CCTK_REAL xx2 = xx2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];

"""+ADMQuantities_CcodeKernel+"""
}
    """)
```

<a id='initial_data__convert_to_hydrobase'></a>

## Step 2.d: Convert TOV solution quantities to `HydroBase` variables $\left\{P,\rho_{\rm baryonic},\epsilon,v_{\rm (n)}^i\right\}$ \[Back to [top](#toc)\]
$$\label{initial_data__convert_to_hydrobase}$$


The TOV solver outputs pressure $P$, the *total* energy density $\rho$, and the baryonic density $\rho_{\rm baryonic}$ as a function of the stellar radius (in isotropic coordinates by default). 

Then, the `HydroBase` quantities $\rho^{\rm HB}_{\rm baryonic}$, internal energy $\epsilon^{\rm HB}$, and pressure $P^{\rm HB}$ are given in terms of these variables via

\begin{align}
P^{\rm HB} &= P; \\
\rho^{\rm HB}_{\rm baryonic} &= \rho_{\rm baryonic}, \\
\rho &= \rho_{\rm baryonic} \left(1 + \epsilon_{\rm cold}\right) \\
\implies \epsilon_{\rm cold} &= \frac{\rho}{\rho_{\rm baryonic}} - 1\\
\epsilon^{\rm HB} &= \epsilon_{\rm cold}, \\
\end{align}
[the NRPy+ piecewise polytrope tutorial notebook](Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb#rhob_from_pcold). Note that $\rho_{\rm baryonic}$ will be floored to a nonzero atmosphere value, so that computing $\epsilon$ will never involve a division by zero.

The TOV star is motionless, with all spatial components of the 4-velocity $u^i=0$ and (as seen above) zero shift $\beta^i$. Thus the Valencia 3-velocity (i.e., the 3-velocity normal to the spatial slice) $v_{\rm (n)}^i$ is given by

$$
v_{\rm (n)}^{i,{\rm HB}} = 0
$$


```python
IDValencia3velocityU = ixp.zerorank1() # Valencia 3-velocity is zero
IDPressure     = par.Cparameters("REAL", thismodule, "IDPressure", 1e300) # IDPressure must be set in C
IDrho_baryonic = par.Cparameters("REAL", thismodule, "IDrho_baryonic", 1e300) # IDrho_baryonic must be set in C
IDrho__total_energy_density = par.Cparameters("REAL", thismodule, "IDrho__total_energy_density", 1e300) # IDrho__total_energy_density must be set in C

# Declare as gridfunctions the final quantities we will output for the initial data
Valencia3velocityU = ixp.register_gridfunctions_for_single_rank1("EVOL","Valencia3velocityU")
Pressure, rho_baryonic, epsilon = gri.register_gridfunctions("EVOL",["Pressure", "rho_baryonic", "epsilon"])

Valencia3velocityU = IDValencia3velocityU # Because all components of Valencia3velocityU are *zero*
Pressure     = IDPressure
rho_baryonic = IDrho_baryonic
epsilon      = IDrho__total_energy_density / IDrho_baryonic - sp.sympify(1)

# -={ Spacetime quantities: Generate C code from expressions and output to file }=-
HydroQuantities_to_print = [\
                            lhrh(lhs=gri.gfaccess("out_gfs","Pressure"),rhs=Pressure),\
                            lhrh(lhs=gri.gfaccess("out_gfs","rho_baryonic"),rhs=rho_baryonic),\
                            lhrh(lhs=gri.gfaccess("out_gfs","epsilon"),rhs=epsilon),\
                            lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU0"),rhs=Valencia3velocityU[0]),\
                            lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU1"),rhs=Valencia3velocityU[1]),\
                            lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU2"),rhs=Valencia3velocityU[2])
                           ]

with open(os.path.join(outdir,"HydroQuantities.h"),"w") as file:
    HydroQuantities_CcodeKernel = fin.FD_outputC("returnstring",HydroQuantities_to_print,
                                               params="outCverbose=False,includebraces=False,preindent=2")
    file.write("""
static inline
void HydroQuantities(const cGH* restrict const cctkGH, const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,

                   const CCTK_REAL IDPressure, const CCTK_REAL IDrho_baryonic,
                   const CCTK_REAL IDrho__total_energy_density,

                   CCTK_REAL *PressureGF,CCTK_REAL *rho_baryonicGF,
                   CCTK_REAL *epsilonGF,
                   CCTK_REAL *Valencia3velocityU0GF,
                   CCTK_REAL *Valencia3velocityU1GF,
                   CCTK_REAL *Valencia3velocityU2GF) {
    DECLARE_CCTK_PARAMETERS;
    if(IDrho__total_energy_density <= 0 || IDrho_baryonic <= 0 || IDPressure <= 0) {
        rho_baryonicGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = rho_atmosphere;
        PressureGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)]     = K_atmosphere*pow(rho_atmosphere,Gamma_atmosphere);
        epsilonGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)]      = 0;
        Valencia3velocityU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
        Valencia3velocityU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
        Valencia3velocityU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
    } else {
"""+HydroQuantities_CcodeKernel+"""
        // Apply pressure depletion.
        PressureGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] *= (1.0 - Pressure_depletion_factor);
    }
}
    """)
```

<a id='einstein'></a>

# Step 3: Interfacing with the Einstein Toolkit \[Back to [top](#toc)\]
$$\label{einstein}$$


<a id='einstein_c'></a>

## Step 3.a: Constructing the Einstein Toolkit C-code calling functions that include the C code kernels \[Back to [top](#toc)\]
$$\label{einstein_c}$$

We will write another C file with the functions we need here.


```python
with open(os.path.join(outdir,"InitialData.c"), "w") as file:
    file.write("""
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

// Declare initial data input struct:
//          stores data from initial data solver,
//          so they can be put on the numerical grid.
typedef struct __ID_inputs {
    CCTK_REAL Rbar;
    int Rbar_idx;
    int interp_stencil_size;
    int numlines_in_file;
    CCTK_REAL *r_Schw_arr,*rho_arr,*rho_baryon_arr,*P_arr,*M_arr,*expnu_arr,*exp4phi_arr,*rbar_arr;
} ID_inputs;


#include "ADMQuantities.h"
#include "HydroQuantities.h"
#include "interpolate_TOV_solution_to_point.h"
#include "convert_TOV_spacetime_vars_to_ADM_vars.h"

// Alias for "vel" vector gridfunction:
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

void read_TOV_input_data_from_file(ID_inputs *TOV_in) {

    DECLARE_CCTK_PARAMETERS;

    // Step 1: Set up TOV initial data
    // Step 1.a: Read TOV initial data from data file
    // Open the data file:
    char filename[100];
    sprintf(filename,"%s",TOV_filename); // TOV_filename is a CCTK_PARAMETER
    FILE *in1Dpolytrope = fopen(filename, "r");
    if (in1Dpolytrope == NULL) {
        fprintf(stderr,"ERROR: could not open file %s\\n",filename);
        exit(1);
    }
    // Count the number of lines in the data file:
    int numlines_in_file = count_num_lines_in_file(in1Dpolytrope);
    // Allocate space for all data arrays:
    CCTK_REAL *r_Schw_arr     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rho_arr        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rho_baryon_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *P_arr          = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *M_arr          = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *expnu_arr      = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *exp4phi_arr    = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rbar_arr       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);

    // Read from the data file, filling in arrays.
    // read_datafile__set_arrays() may be found in TOV/tov_interp.h
    if(read_datafile__set_arrays(in1Dpolytrope, r_Schw_arr,rho_arr,rho_baryon_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr) == 1) {
        fprintf(stderr,"ERROR WHEN READING FILE %s!\\n",filename);
        exit(1);
    }
    fclose(in1Dpolytrope);
    REAL Rbar = -100;
    int Rbar_idx = -100;
    for(int i=1;i<numlines_in_file;i++) {
        if(rho_arr[i-1]>0 && rho_arr[i]==0) { Rbar = rbar_arr[i-1]; Rbar_idx = i-1; }
    }
    if(Rbar<0) {
        fprintf(stderr,"Error: could not find rbar=Rbar from data file.\\n");
        exit(1);
    }

    TOV_in->Rbar = Rbar;
    TOV_in->Rbar_idx = Rbar_idx;

    const int interp_stencil_size = 12;
    TOV_in->interp_stencil_size = interp_stencil_size;
    TOV_in->numlines_in_file = numlines_in_file;

    TOV_in->r_Schw_arr     = r_Schw_arr;
    TOV_in->rho_arr        = rho_arr;
    TOV_in->rho_baryon_arr = rho_baryon_arr;
    TOV_in->P_arr          = P_arr;
    TOV_in->M_arr          = M_arr;
    TOV_in->expnu_arr      = expnu_arr;
    TOV_in->exp4phi_arr    = exp4phi_arr;
    TOV_in->rbar_arr       = rbar_arr;
    /* END TOV INPUT ROUTINE */
}

void NRPyPlusTOVID_ET_InitialData(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    ID_inputs TOV_in;
    read_TOV_input_data_from_file(&TOV_in);

#pragma omp parallel for
  for(CCTK_INT i2=0;i2<cctk_lsh[2];i2++) for(CCTK_INT i1=0;i1<cctk_lsh[1];i1++) for(CCTK_INT i0=0;i0<cctk_lsh[0];i0++) {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i0,i1,i2);
        CCTK_REAL rr = r[idx];
        CCTK_REAL th = acos(z[idx]/rr);

        CCTK_REAL IDexp_4phi,IDnu,IDPressure,IDrho_baryonic,IDrho__total_energy_density;
        interpolate_TOV_solution_to_point(rr, TOV_in, &IDexp_4phi,&IDnu,
                                          &IDPressure,&IDrho_baryonic,&IDrho__total_energy_density);

        CCTK_REAL IDalpha,IDgammaDD00,IDgammaDD01,IDgammaDD02,IDgammaDD11,IDgammaDD12,IDgammaDD22;
        convert_TOV_spacetime_vars_to_ADM_vars(rr, th, IDexp_4phi,IDnu,
            &IDalpha,&IDgammaDD00,&IDgammaDD01,&IDgammaDD02,&IDgammaDD11,&IDgammaDD12,&IDgammaDD22);

        HydroQuantities(cctkGH, i0,i1,i2,
                        IDPressure,IDrho_baryonic,IDrho__total_energy_density,
                        press,rho,eps,velx,vely,velz);

        ADMQuantities(cctkGH,i0,i1,i2,
                      x,y,z,
                      IDalpha,IDgammaDD00,IDgammaDD01,IDgammaDD02,IDgammaDD11,IDgammaDD12,IDgammaDD22,
                      alp,betax,betay,betaz,
                      gxx,gxy,gxz,gyy,gyz,gzz,
                      kxx,kxy,kxz,kyy,kyz,kzz);
  }

  free(TOV_in.r_Schw_arr);
  free(TOV_in.rho_arr);
  free(TOV_in.rho_baryon_arr);
  free(TOV_in.P_arr);
  free(TOV_in.M_arr);
  free(TOV_in.expnu_arr);
  free(TOV_in.exp4phi_arr);
  free(TOV_in.rbar_arr);
}
""")
```

<a id='einstein_ccl'></a>

## Step 3.b: CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{einstein_ccl}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn:

1. `interface.ccl}`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns. Specifically, this file governs the interaction between this thorn and others; more information can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 
With "implements", we give our thorn its unique name. By "inheriting" other thorns, we tell the Toolkit that we will rely on variables that exist and are declared "public" within those functions.


```python
%%writefile $outrootdir/interface.ccl
implements: NRPyPlusTOVID
inherits: admbase grid hydrobase
```

    Overwriting NRPyPlusTOVID//interface.ccl


2. `param.ccl`: specifies free parameters within the thorn, enabling them to be set at runtime. It is required to provide allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3).


```python
%%writefile $outrootdir/param.ccl
shares: grid
shares: ADMBase
USES CCTK_INT lapse_timelevels
USES CCTK_INT shift_timelevels
USES CCTK_INT metric_timelevels

USES KEYWORD metric_type

EXTENDS KEYWORD initial_data
{
  "NRPyPlusTOVID" :: "Initial data from NRPyPlusTOVID solution"
}
EXTENDS KEYWORD initial_lapse
{
  "NRPyPlusTOVID" :: "Initial lapse from NRPyPlusTOVID solution"
}
EXTENDS KEYWORD initial_shift
{
  "NRPyPlusTOVID" :: "Initial shift from NRPyPlusTOVID solution"
}
EXTENDS KEYWORD initial_dtlapse
{
  "NRPyPlusTOVID" :: "Initial dtlapse from NRPyPlusTOVID solution"
}
EXTENDS KEYWORD initial_dtshift
{
  "NRPyPlusTOVID" :: "Initial dtshift from NRPyPlusTOVID solution"
}

shares: HydroBase
EXTENDS KEYWORD initial_hydro
{
  "NRPyPlusTOVID" :: "Initial GRHD data from NRPyPlusTOVID solution"
}

#["r_in","r_at_max_density","a","M"] A_b, kappa, gamma
restricted:
CCTK_STRING TOV_filename "Which interpolator should I use"
{
 ".+" :: "Any nonempty string"
} "outputTOVpolytrope.txt"

restricted:
CCTK_REAL rho_atmosphere "Atmosphere baryonic density"
{
 0:* :: "Physical values"
-1   :: "forbidden value to make sure it is explicitly set in the parfile"
} -1

restricted:
CCTK_REAL K_atmosphere "Polytropic K to be used with the EOS corresponding to rho_atmosphere"
{
 0:* :: "Physical values"
-1   :: "forbidden value to make sure it is explicitly set in the parfile"
} -1

restricted:
CCTK_REAL Gamma_atmosphere "Polytropic Gamma to be used with the EOS corresponding to rho_atmosphere"
{
 0:* :: "Physical values"
-1   :: "forbidden value to make sure it is explicitly set in the parfile"
} -1

restricted:
CCTK_REAL Pressure_depletion_factor "Pressure depletion factor = Pdf: P => (1-Pdf)*P"
{
 0:* :: "Greater than or equal to zero, where zero is no depletion and default."
} 0.0
```

    Overwriting NRPyPlusTOVID//param.ccl


3. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions. $\text{schedule.ccl}$'s official documentation may be found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). 

We specify here the standardized ETK "scheduling bins" in which we want each of our thorn's functions to run.


```python
%%writefile $outrootdir/schedule.ccl
STORAGE: ADMBase::metric[metric_timelevels], ADMBase::curv[metric_timelevels], ADMBase::lapse[lapse_timelevels], ADMBase::shift[shift_timelevels]

schedule NRPyPlusTOVID_ET_InitialData IN HydroBase_Initial
{
  LANG: C
  READS: grid::x(Everywhere)
  READS: grid::y(Everywhere)
  READS: grid::y(Everywhere)
  WRITES: admbase::alp(Everywhere)
  WRITES: admbase::betax(Everywhere)
  WRITES: admbase::betay(Everywhere)
  WRITES: admbase::betaz(Everywhere)
  WRITES: admbase::kxx(Everywhere)
  WRITES: admbase::kxy(Everywhere)
  WRITES: admbase::kxz(Everywhere)
  WRITES: admbase::kyy(Everywhere)
  WRITES: admbase::kyz(Everywhere)
  WRITES: admbase::kzz(Everywhere)
  WRITES: admbase::gxx(Everywhere)
  WRITES: admbase::gxy(Everywhere)
  WRITES: admbase::gxz(Everywhere)
  WRITES: admbase::gyy(Everywhere)
  WRITES: admbase::gyz(Everywhere)
  WRITES: admbase::gzz(Everywhere)
  WRITES: hydrobase::vel[0](Everywhere)
  WRITES: hydrobase::vel[1](Everywhere)
  WRITES: hydrobase::vel[2](Everywhere)
  WRITES: hydrobase::rho(Everywhere)
  WRITES: hydrobase::eps(Everywhere)
  WRITES: hydrobase::press(Everywhere)
} "Set up general relativistic hydrodynamic (GRHD) fields for TOV initial data"
```

    Overwriting NRPyPlusTOVID//schedule.ccl


<a id='einstein_list'></a>

## Step 3.c: Add the C code to the Einstein Toolkit compilation list \[Back to [top](#toc)\]
$$\label{einstein_list}$$

We will also need `make.code.defn`, which indicates the list of files that need to be compiled. This thorn only has the one C file to compile.


```python
%%writefile $outdir/make.code.defn
SRCS = InitialData.c
```

    Overwriting NRPyPlusTOVID/src/make.code.defn


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-NRPyPlusTOVID.pdf](Tutorial-ETK_thorn-NRPyPlusTOVID.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-NRPyPlusTOVID")
```

    Created Tutorial-ETK_thorn-NRPyPlusTOVID.tex, and compiled LaTeX file to
        PDF file Tutorial-ETK_thorn-NRPyPlusTOVID.pdf

