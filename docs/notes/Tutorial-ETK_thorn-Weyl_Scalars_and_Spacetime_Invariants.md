<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Weyl Scalars and Invariants: An Introduction to Einstein Toolkit Diagnostic Thorns

## Author: Patrick Nelson & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This tutorial describes the creation of a diagnostic thorn for the Einstein Toolkit (ETK) that computes Weyl scalars and invariants, derived from SymPy expressions in the NRPy+ tutorial. Steps involve converting SymPy expressions into C-code kernels and building ETK infrastructure (interface.ccl, param.ccl, and schedule.ccl files). The result is a versatile module for efficient numerical simulations.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** Numerical results from this module have been confirmed to agree with the trusted WeylScal4 Einstein Toolkit thorn to roundoff error.

### NRPy+ Source Code for this module:  
* [WeylScal4NRPy/WeylScalars_Cartesian.py](../edit/WeylScal4NRPy/WeylScalars_Cartesian.py)
* [WeylScal4NRPy/WeylScalarInvariants_Cartesian.py](../edit/WeylScal4NRPy/WeylScalarInvariants_Cartesian.py)

which are fully documented in the NRPy+ [Tutorial-WeylScalars-Cartesian](Tutorial-WeylScalars-Cartesian.ipynb) module on using NRPy+ to construct the Weyl scalars and invariants as SymPy expressions.

## Introduction:
In the [previous tutorial notebook](Tutorial-WeylScalars-Cartesian.ipynb), we constructed within SymPy full expressions for the real and imaginary components of all five Weyl scalars $\psi_0$, $\psi_1$, $\psi_2$, $\psi_3$, and $\psi_4$ as well as the Weyl invariants. So that we can easily access these expressions, we have ported the Python code needed to generate the Weyl scalar SymPy expressions to [WeylScal4NRPy/WeylScalars_Cartesian.py](../edit/WeylScal4NRPy/WeylScalars_Cartesian.py), and the Weyl invariant SymPy expressions to [WeylScal4NRPy/WeylScalarInvariants_Cartesian.py](../edit/WeylScal4NRPy/WeylScalarInvariants_Cartesian.py).

Here we will work through the steps necessary to construct an Einstein Toolkit diagnostic thorn (module), starting from these SymPy expressions, which computes these expressions using ADMBase gridfunctions as input. This tutorial is in two steps:

1. Call on NRPy+ to convert the SymPy expressions for the Weyl Scalars and associated Invariants into one C-code kernel for each.
1. Write the C code and build up the needed Einstein Toolkit infrastructure (i.e., the .ccl files).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 1](#nrpy): Call on NRPy+ to convert the SymPy expressions for the Weyl scalars and associated invariants into one C-code kernel for each
1. [Step 2](#etk): Interfacing with the Einstein Toolkit
    1. [Step 2.a](#etkc): Constructing the Einstein Toolkit C-code calling functions that include the C code kernels
    1. [Step 2.b](#cclfiles): CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure
    1. [Step 2.c](#etk_list): Add the C file to Einstein Toolkit compilation list
1. [Step 3](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='nrpy'></a>

# Step 1: Call on NRPy+ to convert the SymPy expressions for the Weyl scalars and associated invariants into one C-code kernel for each \[Back to [top](#toc)\]
$$\label{nrpy}$$

<font color='red'><b>WARNING</b></font>: It takes some time to generate the CSE-optimized C code kernels for these quantities, especially the Weyl scalars... expect 5 minutes on a modern computer.


```python
from outputC import lhrh, outCfunction  # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import os                        # Standard Python modules for multiplatform OS-level functions

# Since we are writing an Einstein Toolkit thorn, we must set our memory access style to "ETK".
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

import WeylScal4NRPy.WeylScalars_Cartesian as weyl
par.set_parval_from_str("output_scalars","all_psis_and_invariants")

weyl.WeylScalars_Cartesian()

output_scalars = par.parval_from_str("output_scalars")

!mkdir WeylScal4NRPy     2>/dev/null # 2>/dev/null: Don't throw an error or warning if the directory already exists.
!mkdir WeylScal4NRPy/src 2>/dev/null # 2>/dev/null: Don't throw an error or warning if the directory already exists.

scalars_lhrh = [lhrh(lhs=gri.gfaccess("out_gfs","psi4r"),rhs=weyl.psi4r),
                lhrh(lhs=gri.gfaccess("out_gfs","psi4i"),rhs=weyl.psi4i)]

if output_scalars in ('all_psis', 'all_psis_and_invariants'):
    scalars_lhrh = [
                    lhrh(lhs=gri.gfaccess("out_gfs","psi4r"),rhs=weyl.psi4r),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi4i"),rhs=weyl.psi4i),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi3r"),rhs=weyl.psi3r),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi3i"),rhs=weyl.psi3i),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi2r"),rhs=weyl.psi2r),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi2i"),rhs=weyl.psi2i),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi1r"),rhs=weyl.psi1r),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi1i"),rhs=weyl.psi1i),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi0r"),rhs=weyl.psi0r),
                    lhrh(lhs=gri.gfaccess("out_gfs","psi0i"),rhs=weyl.psi0i)
                   ]

# Generating the CSE is the slowest
# operation in this notebook, and much of the CSE
# time is spent sorting CSE expressions. Disabling
# this sorting makes the C codegen 3-4x faster,
# but the tradeoff is that every time this is
# run, the CSE patterns will be different
# (though they should result in mathematically
# *identical* expressions). You can expect
# roundoff-level differences as a result.
psis_CcodeKernel = fin.FD_outputC("returnstring",scalars_lhrh,params="outCverbose=False,CSE_sorting=none")

cctk_params_readin_str = r"""
  DECLARE_CCTK_PARAMETERS;
  if(cctk_nghostzones[0] != cctk_nghostzones[1] ||
     cctk_nghostzones[0] != cctk_nghostzones[2]) {
    CCTK_ERROR("cctk_nghostzones[i] must be the same in all directions i");
  }
  const CCTK_INT NGHOSTS CCTK_ATTRIBUTE_UNUSED = cctk_nghostzones[0];
  const CCTK_INT Nxx0 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[0] - 2*NGHOSTS;
  const CCTK_INT Nxx1 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[1] - 2*NGHOSTS;
  const CCTK_INT Nxx2 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[2] - 2*NGHOSTS;
  const CCTK_INT Nxx_plus_2NGHOSTS0 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[0];
  const CCTK_INT Nxx_plus_2NGHOSTS1 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[1];
  const CCTK_INT Nxx_plus_2NGHOSTS2 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[2];
"""

desc = "Calculate the Weyl Scalars"
name = "calc_psis"
outCfunction(
    outfile  = os.path.join("WeylScal4NRPy","src",name+".h"), desc=desc, name=name,
    params   ="""const cGH* restrict const cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,
               const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,
               const CCTK_REAL *xGF,const CCTK_REAL *yGF,const CCTK_REAL *zGF,
               const CCTK_REAL *gammaDD00GF,const CCTK_REAL *gammaDD01GF,const CCTK_REAL *gammaDD02GF,const CCTK_REAL *gammaDD11GF,const CCTK_REAL *gammaDD12GF,const CCTK_REAL *gammaDD22GF,
               const CCTK_REAL     *kDD00GF,const CCTK_REAL     *kDD01GF,const CCTK_REAL     *kDD02GF,const CCTK_REAL     *kDD11GF,const CCTK_REAL     *kDD12GF,const CCTK_REAL     *kDD22GF,
               CCTK_REAL *psi4rGF,CCTK_REAL *psi4iGF,
               CCTK_REAL *psi3rGF,CCTK_REAL *psi3iGF,
               CCTK_REAL *psi2rGF,CCTK_REAL *psi2iGF,
               CCTK_REAL *psi1rGF,CCTK_REAL *psi1iGF,
               CCTK_REAL *psi0rGF,CCTK_REAL *psi0iGF""",
    preloop=cctk_params_readin_str,
    body     = psis_CcodeKernel,
    loopopts ="InteriorPoints",
    enableCparameters=False)

# Reset the registered gridfunctions list.
gri.glb_gridfcs_list = []
#par.set_parval_from_str("WeylScal4NRPy.WeylScalars_Cartesian::output_scalars","all_psis_and_invariants")
output_scalars = par.parval_from_str("output_scalars")

import WeylScal4NRPy.WeylScalarInvariants_Cartesian as invar
invar.WeylScalarInvariants_Cartesian()
invars_lhrh = [
               lhrh(lhs=gri.gfaccess("out_gfs","curvIr"),rhs=invar.curvIr),
               lhrh(lhs=gri.gfaccess("out_gfs","curvIi"),rhs=invar.curvIi),
               lhrh(lhs=gri.gfaccess("out_gfs","curvJr"),rhs=invar.curvJr),
               lhrh(lhs=gri.gfaccess("out_gfs","curvJi"),rhs=invar.curvJi),
               lhrh(lhs=gri.gfaccess("out_gfs","J1curv"),rhs=invar.J1curv),
               lhrh(lhs=gri.gfaccess("out_gfs","J2curv"),rhs=invar.J2curv),
               lhrh(lhs=gri.gfaccess("out_gfs","J3curv"),rhs=invar.J3curv),
               lhrh(lhs=gri.gfaccess("out_gfs","J4curv"),rhs=invar.J4curv)
              ]

invars_CcodeKernel = fin.FD_outputC("returnstring",invars_lhrh,
                    params="outCverbose=False,CSE_sorting=none")# Generating the CSE is the slowest
                                                                # operation in this notebook, and much of the CSE
                                                                # time is spent sorting CSE expressions. Disabling
                                                                # this sorting makes the C codegen 3-4x faster,
                                                                # but the tradeoff is that every time this is
                                                                # run, the CSE patterns will be different
                                                                # (though they should result in mathematically
                                                                # *identical* expressions). You can expect
                                                                # roundoff-level differences as a result.

desc = "Calculate the Weyl Invariants"
name = "calc_invars"
outCfunction(
    outfile  = os.path.join("WeylScal4NRPy","src",name+".h"), desc=desc, name=name,
    params   ="""const cGH* restrict const cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,
                 const CCTK_REAL *psi4rGF,const CCTK_REAL *psi4iGF,
                 const CCTK_REAL *psi3rGF,const CCTK_REAL *psi3iGF,
                 const CCTK_REAL *psi2rGF,const CCTK_REAL *psi2iGF,
                 const CCTK_REAL *psi1rGF,const CCTK_REAL *psi1iGF,
                 const CCTK_REAL *psi0rGF,const CCTK_REAL *psi0iGF,
                 CCTK_REAL *curvIrGF,CCTK_REAL *curvIiGF,
                 CCTK_REAL *curvJrGF,CCTK_REAL *curvJiGF,
                 CCTK_REAL *J1curvGF,CCTK_REAL *J2curvGF,
                 CCTK_REAL *J3curvGF,CCTK_REAL *J4curvGF""",
    preloop=cctk_params_readin_str,
    body     = invars_CcodeKernel,
    loopopts ="InteriorPoints",
    enableCparameters=False)

```

    Output C function calc_psis() to file WeylScal4NRPy/src/calc_psis.h
    Output C function calc_invars() to file WeylScal4NRPy/src/calc_invars.h


<a id='etk'></a>

# Step 2: Interfacing with the Einstein Toolkit \[Back to [top](#toc)\]
$$\label{etk}$$


<a id='etkc'></a>

## Step 2.a: Constructing the Einstein Toolkit calling functions that include the C code kernels \[Back to [top](#toc)\]
$$\label{etkc}$$

Now that we have generated the C code kernels (`WeylScal4NRPy_psis.h` and `WeylScal4NRPy_invars.h`) express the Weyl scalars and invariants as CSE-optimized finite-difference expressions, we next need to write the C code functions that incorporate these kernels and are called by the Einstein Toolkit scheduler.


```python
%%writefile WeylScal4NRPy/src/WeylScal4NRPy.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "calc_psis.h"
#include "calc_invars.h"

extern void weylscal4_mainfunction(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if(cctk_iteration % WeylScal4NRPy_calc_every != 0) { return; }

  const CCTK_REAL invdx0 = 1.0 / (CCTK_DELTA_SPACE(0));
  const CCTK_REAL invdx1 = 1.0 / (CCTK_DELTA_SPACE(1));
  const CCTK_REAL invdx2 = 1.0 / (CCTK_DELTA_SPACE(2));

  /* Now, to calculate psi4: */
  calc_psis(cctkGH,cctk_lsh,cctk_nghostzones,
            invdx0,invdx1,invdx2,
            x,y,z,
            gxx,gxy,gxz,gyy,gyz,gzz,
            kxx,kxy,kxz,kyy,kyz,kzz,
            psi4r,psi4i,
            psi3r,psi3i,
            psi2r,psi2i,
            psi1r,psi1i,
            psi0r,psi0i);

  if (CCTK_EQUALS(output_scalars, "all_psis_and_invariants")) {
    calc_invars(cctkGH,cctk_lsh,cctk_nghostzones,
      	        psi4r,psi4i,
                psi3r,psi3i,
                psi2r,psi2i,
                psi1r,psi1i,
                psi0r,psi0i,
                NRPycurvIr,NRPycurvIi,
                NRPycurvJr,NRPycurvJi,
                NRPyJ1curv,NRPyJ2curv,
                NRPyJ3curv,NRPyJ4curv);
    }
}

```

    Overwriting WeylScal4NRPy/src/WeylScal4NRPy.c


<a id='cclfiles'></a>

## Step 2.b: CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{cclfiles}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn:

1. `interface.ccl`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns.
1. `param.ccl`: specifies free parameters within the thorn.
1. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions.

Let's start with `interface.ccl`. The [official Einstein Toolkit (Cactus) documentation](http://einsteintoolkit.org/usersguide/UsersGuide.html) defines what must/should be included in an `interface.ccl` file [**here**](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 


```python
%%writefile WeylScal4NRPy/interface.ccl

# With "implements", we give our thorn its unique name.
implements: WeylScal4NRPy

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits:   admbase Boundary Grid methodoflines

# Tell the Toolkit that we want the various Weyl scalars
#    and invariants to be visible to other thorns by using
#    the keyword "public". Note that declaring these
#    gridfunctions *does not* allocate memory for them;
#    that is done by the schedule.ccl file.
public:
CCTK_REAL NRPyPsi4_group type=GF timelevels=3 tags='tensortypealias="Scalar" tensorweight=0 tensorparity=1'
{
  psi4r, psi4i
} "Psi4_group"

public:
CCTK_REAL NRPyPsi3210_group type=GF timelevels=3 tags='tensortypealias="Scalar" tensorweight=0 tensorparity=1'
{
  psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i
} "Psi3210_group"

public:
CCTK_REAL NRPyInvars_group type=GF timelevels=3 tags='tensortypealias="Scalar" tensorweight=0 tensorparity=1'
{
  NRPycurvIr,NRPycurvIi,NRPycurvJr,NRPycurvJi,NRPyJ1curv,NRPyJ2curv,NRPyJ3curv,NRPyJ4curv
} "NRPyInvars_group"
```

    Overwriting WeylScal4NRPy/interface.ccl


We will now write the file `param.ccl`. This file allows the listed parameters to be set at runtime. We also give allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3). 

The first parameter specifies how many time levels need to be stored. Generally when using the ETK's adaptive-mesh refinement (AMR) driver [Carpet](https://carpetcode.org/), three timelevels are needed so that the diagnostic quantities can be properly interpolated and defined across refinement boundaries. 

The second parameter determines how often we will calculate $\psi_4$, and the third parameter indicates whether just $\psi_4$, all Weyl scalars, or all Weyl scalars and invariants are going to be output. The third parameter is currently specified entirely within NRPy+, so by this point, it is *not* a free parameter. Thus it is not quite correct to include it in this list of *free* parameters (FIXME).


```python
%%writefile WeylScal4NRPy/param.ccl

restricted:
CCTK_INT timelevels "Number of active timelevels" STEERABLE=RECOVER
{
  0:3 :: ""
} 3

restricted:
CCTK_INT WeylScal4NRPy_calc_every "WeylScal4_psi4_calc_Nth_calc_every" STEERABLE=ALWAYS
{
  *:* :: ""
} 1

private:
CCTK_KEYWORD output_scalars "Whether to output all Weyl scalars, just psi4, or all scalars and invariants"
{
  "all_psis" :: ""
  "all_psis_and_invariants" :: ""
} "all_psis"

```

    Overwriting WeylScal4NRPy/param.ccl


Finally, we will write the file `schedule.ccl`; its official documentation is found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). This file dictates when the various parts of the thorn will be run. We first assign storage for both the real and imaginary components of $\psi_4$, and then specify that we want our code run in  the `MoL_PseudoEvolution` schedule group (consistent with the original `WeylScal4` Einstein Toolkit thorn), after the ADM variables are set. At this step, we declare that we will be writing code in C. We also specify the gridfunctions that we wish to read in from memory--in our case, we need all the components of $K_{ij}$ (the spatial extrinsic curvature) and $\gamma_{ij}$ (the physical [as opposed to conformal] 3-metric), in addition to the coordinate values. Note that the ETK adopts the widely-used convention that components of $\gamma_{ij}$ are prefixed in the code with $\text{g}$ and not $\gamma$.


```python
%%writefile WeylScal4NRPy/schedule.ccl

STORAGE: NRPyPsi4_group[timelevels]
if (CCTK_EQUALS(output_scalars, "all_psis_and_invariants") || CCTK_EQUALS(output_scalars, "all_psis"))
{
 STORAGE: NRPyPsi3210_group[timelevels]
}
if (CCTK_EQUALS(output_scalars, "all_psis_and_invariants"))
{
 STORAGE: NRPyInvars_group[timelevels]
}

schedule group WeylScal4NRPy_group in MoL_PseudoEvolution after ADMBase_SetADMVars
{
} "Schedule WeylScal4NRPy group"

schedule weylscal4_mainfunction in WeylScal4NRPy_group
{
  LANG: C
   READS: admbase::kxx(Everywhere)
   READS: admbase::kxy(Everywhere)
   READS: admbase::kxz(Everywhere)
   READS: admbase::kyy(Everywhere)
   READS: admbase::kyz(Everywhere)
   READS: admbase::kzz(Everywhere)
   READS: admbase::gxx(Everywhere)
   READS: admbase::gxy(Everywhere)
   READS: admbase::gxz(Everywhere)
   READS: admbase::gyy(Everywhere)
   READS: admbase::gyz(Everywhere)
   READS: admbase::gzz(Everywhere)
   READS: grid::x(Everywhere)
   READS: grid::y(Everywhere)
   READS: grid::z(Everywhere)
   WRITES: WeylScal4NRPy::psi4i(Interior)
   WRITES: WeylScal4NRPy::psi4r(Interior)
   WRITES: WeylScal4NRPy::psi3i(Interior)
   WRITES: WeylScal4NRPy::psi3r(Interior)
   WRITES: WeylScal4NRPy::psi2i(Interior)
   WRITES: WeylScal4NRPy::psi2r(Interior)
   WRITES: WeylScal4NRPy::psi1i(Interior)
   WRITES: WeylScal4NRPy::psi1r(Interior)
   WRITES: WeylScal4NRPy::psi0i(Interior)
   WRITES: WeylScal4NRPy::psi0r(Interior)
   WRITES: WeylScal4NRPy::NRPycurvIi(Interior)
   WRITES: WeylScal4NRPy::NRPycurvIr(Interior)
   WRITES: WeylScal4NRPy::NRPycurvJi(Interior)
   WRITES: WeylScal4NRPy::NRPycurvJr(Interior)
   WRITES: WeylScal4NRPy::NRPyJ1curv(Interior)
   WRITES: WeylScal4NRPy::NRPyJ2curv(Interior)
   WRITES: WeylScal4NRPy::NRPyJ3curv(Interior)
   WRITES: WeylScal4NRPy::NRPyJ4curv(Interior)

} "Call WeylScal4NRPy main function"

```

    Overwriting WeylScal4NRPy/schedule.ccl


<a id='etk_list'></a>

## Step 2.c: Tell the Einstein Toolkit to compile the C code \[Back to [top](#toc)\]
$$\label{etk_list}$$

The `make.code.defn` lists the source files that need to be compiled. Naturally, this thorn has only the one C file $-$ written above $-$ to compile:


```python
%%writefile WeylScal4NRPy/src/make.code.defn

SRCS = WeylScal4NRPy.c
```

    Overwriting WeylScal4NRPy/src/make.code.defn


<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-Weyl_Scalars_and_Spacetime_Invariants.pdf](Tutorial-ETK_thorn-Weyl_Scalars_and_Spacetime_Invariants.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-Weyl_Scalars_and_Spacetime_Invariants")
```

    Created Tutorial-ETK_thorn-Weyl_Scalars_and_Spacetime_Invariants.tex, and
        compiled LaTeX file to PDF file Tutorial-ETK_thorn-
        Weyl_Scalars_and_Spacetime_Invariants.pdf

