<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# sbPoynETNRPy: An Einstein Toolkit Thorn for Computing the 4-Velocity Time-Component $u^0$, the Magnetic Field Measured by a Comoving Observer $b^{\mu}$, and the Poynting Vector $S^i$

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook describes how the `sbPoynETNRPy` thorn uses SymPy expressions and NRPy+ for generating C-code kernels, computing 4-velocity time-component $u^0$, comoving observer's magnetic field $b^\mu$, and Poynting vector $S_i$. It elaborates on integration with the Einstein Toolkit infrastructure, the building of specific ccl files, and the use of NRPy+'s outputC for relevant expression generation.
<br>

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated against the hand-written smallbPoynET in WVUThorns_diagnostics (a trusted code), which itself is based on expressions in IllinoisGRMHD... which was validated against the original GRMHD code of the Illinois NR group.

## Introduction:
In the [previous tutorial notebook](Tutorial-u0_smallb_Poynting-Cartesian.ipynb), we constructed within SymPy full expressions for the 4-velocity time-component $u^0$, the magnetic field (measured by a comoving observer) $b^{\mu}$, and the Poynting vector $S^i$.

Here we will work through the steps necessary to construct an Einstein Toolkit diagnostic thorn (module) that uses ADMBase and HydroBase variables as input into the NRPy+-generated SymPy expressions for $b^{\mu}$, $b^2$, and the Poynting Vector $S^i$, outputting to gridfunctions `smallb4U[]`, `smallb2etk` (the "etk" suffix must be appended because base gridfunction names ending in numbers are not allowed in NRPy+), and `SPoyn[]`, respectively. 

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This tutorial is organized as follows

1. [Step 1](#initializenrpy): Call on NRPy+ to convert the SymPy expressions for $b^{\mu}$, $b^2$, and the Poynting Vector $S^i$ to C code kernels
1. [Step 2](#etk): Build up the needed Einstein Toolkit infrastructure to implement the NRPy+-generated C code kernels
    1. [Step 2.a](#etkc): Write the C code functions called by the Einstein Toolkit scheduler that incorporate the above ".h" files
    1. [Step 2.b](#cclfiles): CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure
    1. [Step 2.c](#etksys): Inform the Einstein Toolkit build system of the C code
1. [Step 3](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a> 

# Step 1: Call on NRPy+ to convert the SymPy expressions for $b^{\mu}$, $b^2$, and the Poynting Vector $S^i$ to C code kernels \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1a: import all needed modules from NRPy+:
import indexedexp as ixp        # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from outputC import outputC     # NRPy+: Basic C code output functionality
import grid as gri              # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par  # NRPy+: parameter interface
```


```python
# Step 1b: Initialize  parameters (stub; there are none for this module)
thismodule = __name__
```

We will disable verbose output in the NRPy+ outputC function. This is an important step in this case because our final expressions are very large. Verbose output, when enabled, will print (in comments) the input SymPy expressions to the top of the file *without* CSE, resulting here in an *enormous* output file.

We will also declare the additional gridfunctions we need for this thorn:

**Inputs from ADMBase:**
* the physical metric $\gamma_{ij}$
* the spacetime gauge quantities $\alpha$ and $\beta^i$

**Inputs from HydroBase:**
* the Valencia 3-velocity $v^i_{(n)}$
* the densitized magnetic field of a normal observer $\tilde{B}^i$

**Output gridfunctions:**
* the magnetic field as observed in a frame comoving with the plasma $b^\mu$ (`smallb4U[]}`)
* twice the magnetic pressure $2 P_{\rm mag} = b_\mu b^\mu = b^2$ (`smallb2etk`)
* the Poynting vector $S^i$ (`SPoyn[]`)


```python
# Step 1c: Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

# Step 1d: declare the additional gridfunctions (i.e., functions whose values are declared
#          at every grid point, either inside or outside of our SymPy expressions) needed
#         for this thorn

# INPUT GRIDFUNCTIONS:
gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01") # The AUX or EVOL designation is *not*
                                                                                # used in diagnostic modules.
betaU = ixp.register_gridfunctions_for_single_rank1("AUX","betaU") # The AUX or EVOL designation is *not*
                                                                   # used in diagnostic modules.
alpha = gri.register_gridfunctions("AUX","alpha") # The AUX or EVOL designation is *not*
                                                  # used in diagnostic modules.
ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU") # The AUX or EVOL designation is *not*
                                                                             # used in diagnostic modules.
BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU") # The AUX or EVOL designation is *not*
                                                             # used in diagnostic modules.

# OUTPUT GRIDFUNCTIONS:
smallb4U = ixp.register_gridfunctions_for_single_rank1("AUX","smallb4U",DIM=4) # The AUX or EVOL designation is *not*
                                                                               # used in diagnostic modules.
smallb2etk = gri.register_gridfunctions("AUX","smallb2etk") # The AUX or EVOL designation is *not*
                                                            # used in diagnostic modules.
PoynSU = ixp.register_gridfunctions_for_single_rank1("AUX","PoynSU") # The AUX or EVOL designation is *not*
                                                                     # used in diagnostic modules.

# Step 1f: Call the NRPy+ module to set up the SymPy expressions for the output, as well as the C code for computing u^0
import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0etc
u0etc.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)

# Step 1g: Set the gridfunction memory access type to "ETK":
par.set_parval_from_str("GridFuncMemAccess","ETK")

# Step 1h: Make output directories:
!mkdir sbPoynETNRPy     2>/dev/null # 2>/dev/null: Don't throw an error or warning if the directory already exists.
!mkdir sbPoynETNRPy/src 2>/dev/null # 2>/dev/null: Don't throw an error or warning if the directory already exists.

# Step 1i: Output routine for computing u0:
with open("sbPoynETNRPy/src/u0.h", "w") as file:
    file.write(str(u0etc.computeu0_Cfunction))
    print("Wrote to file \""+file.name+"\"")

# Step 1j: Use NRPy+'s outputC to convert the SymPy expressions for smallb4U, smallb2etk, and PoynSU to C code:
#outputC([u0etc.smallb4U[0],u0etc.smallb4U[1],u0etc.smallb4U[2],u0etc.smallb4U[3],u0etc.smallb2etk,
outputC([u0etc.smallb4U[0],u0etc.smallb4U[1],u0etc.smallb4U[2],u0etc.smallb4U[3],u0etc.smallb2etk,
         u0etc.PoynSU[0],u0etc.PoynSU[1],u0etc.PoynSU[2]],
         [gri.gfaccess("","smallb4U0"),gri.gfaccess("","smallb4U1"),gri.gfaccess("","smallb4U2"),gri.gfaccess("","smallb4U3"),
          gri.gfaccess("","smallb2etk"),
          gri.gfaccess("","PoynSU0"),gri.gfaccess("","PoynSU1"),gri.gfaccess("","PoynSU2")],
       filename="sbPoynETNRPy/src/smallb4U_smallb2etk_PoynSU.h",
       params="outCverbose=False") # <- Force outCverbose=False for this
                                   #    module to avoid gigantic C file filled with the
                                   #    non-CSE expressions for the Weyl scalars.
```

    Wrote to file "sbPoynETNRPy/src/u0.h"
    Wrote to file "sbPoynETNRPy/src/smallb4U_smallb2etk_PoynSU.h"


<a id='etk'></a>

# Step 2: Build up the needed Einstein Toolkit infrastructure to implement the NRPy+-generated C code kernels \[Back to [top](#toc)\]
$$\label{etk}$$


<a id='etkc'></a> 

## Step 2.a: Write the C code functions called by the Einstein Toolkit scheduler that incorporate the above ".h" files \[Back to [top](#toc)\]
$$\label{etkc}$$


```python
%%writefile sbPoynETNRPy/src/sbPoynETNRPy.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void sbPoynETNRPy_lowlevel(const cGH* restrict const cctkGH,const int *cctk_lsh,
                           const CCTK_REAL *gammaDD00GF,const CCTK_REAL *gammaDD01GF,const CCTK_REAL *gammaDD02GF,
                           const CCTK_REAL *gammaDD11GF,const CCTK_REAL *gammaDD12GF,const CCTK_REAL *gammaDD22GF,
                           const CCTK_REAL *alphaGF,
                           const CCTK_REAL *betaU0GF,const CCTK_REAL *betaU1GF,const CCTK_REAL *betaU2GF,
                           const CCTK_REAL *vel,const CCTK_REAL *Bvec,

                           CCTK_REAL *smallb4U0GF,CCTK_REAL *smallb4U1GF,CCTK_REAL *smallb4U2GF,CCTK_REAL *smallb4U3GF,
                           CCTK_REAL *smallb2etkGF,
                           CCTK_REAL *PoynSU0GF,CCTK_REAL *PoynSU1GF,CCTK_REAL *PoynSU2GF) {

    DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
    for(int i2=0;i2<cctk_lsh[2];i2++) for(int i1=0;i1<cctk_lsh[1];i1++) for(int i0=0;i0<cctk_lsh[0];i0++) {
        const CCTK_REAL gammaDD00 = gammaDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL gammaDD01 = gammaDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL gammaDD02 = gammaDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL gammaDD11 = gammaDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL gammaDD12 = gammaDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL gammaDD22 = gammaDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];

        const CCTK_REAL alpha = alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];

        const CCTK_REAL betaU0 = betaU0GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL betaU1 = betaU1GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
        const CCTK_REAL betaU2 = betaU2GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];

        // Valencia 3-velocity may be adjusted due to the velocity ceiling.
        CCTK_REAL ValenciavU0 = vel[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 0)];
        CCTK_REAL ValenciavU1 = vel[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 1)];
        CCTK_REAL ValenciavU2 = vel[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 2)];

        const CCTK_REAL BU0 = Bvec[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 0)];
        const CCTK_REAL BU1 = Bvec[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 1)];
        const CCTK_REAL BU2 = Bvec[CCTK_GFINDEX4D(cctkGH, i0,i1,i2, 2)];

        CCTK_REAL u0;
#include "u0.h"
#include "smallb4U_smallb2etk_PoynSU.h"
    }
}

extern void sbPoynETNRPy(CCTK_ARGUMENTS) {

    DECLARE_CCTK_PARAMETERS;
    DECLARE_CCTK_ARGUMENTS;

    if(sbPoynETNRPy_calc_every<=0 || cctk_iteration%sbPoynETNRPy_calc_every!=0) { return; }

    /* Calculate smallb4U[], smallb2etk, and PoynSU[]: */
    sbPoynETNRPy_lowlevel(cctkGH,cctk_lsh,
                          gxx,gxy,gxz,gyy,gyz,gzz,
                          alp,
                          betax,betay,betaz,
                          vel,Bvec,

                          smallb4U0,smallb4U1,smallb4U2,smallb4U3,
                          smallb4_sq,
                          PoynSU0,PoynSU1,PoynSU2);
}
```

    Overwriting sbPoynETNRPy/src/sbPoynETNRPy.c


<a id='cclfiles'></a> 

## Step 2.b: CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{cclfiles}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn:

1. `interface.ccl`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns.
1. `param.ccl`: specifies free parameters within the thorn.
1. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions.

Let's start with `interface.ccl`. The [official Einstein Toolkit (Cactus) documentation](http://einsteintoolkit.org/usersguide/UsersGuide.html) defines what must/should be included in an `interface.ccl` file [**here**](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 


```python
%%writefile sbPoynETNRPy/interface.ccl

# With "implements", we give our thorn its unique name.
implements: sbPoynETNRPy

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits:   ADMBase Boundary Grid HydroBase MethodofLines

# Tell the Toolkit that we want the various Weyl scalars
#    and invariants to be visible to other thorns by using
#    the keyword "public". Note that declaring these
#    gridfunctions *does not* allocate memory for them;
#    that is done by the schedule.ccl file.
public:
CCTK_REAL smallb4U_group type=GF timelevels=3
{
  smallb4U0,smallb4U1,smallb4U2,smallb4U3
} "smallb4U 4-vector"

public:
CCTK_REAL smallb4_sq_group type=GF timelevels=3
{
  smallb4_sq
} "smallb^{mu} squared == twice the magnetic pressure"

public:
CCTK_REAL PoynSU_group type=GF timelevels=3
{
  PoynSU0,PoynSU1,PoynSU2
} "Poynting 3-vector"
```

    Overwriting sbPoynETNRPy/interface.ccl


We will now write the file `param.ccl`. This file allows the listed parameters to be set at runtime. We also give allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3). 

The first parameter specifies how many time levels need to be stored. Generally when using the ETK's adaptive-mesh refinement (AMR) driver [Carpet](https://carpetcode.org/), three time levels are needed so that the diagnostic quantities can be properly interpolated and defined across refinement boundaries. 

The second parameter determines how often we will calculate $b^\mu$, $b^2$, and $S^i$.

The third parameter sets the maximum allowed Lorentz factor when computing $u^0$ (i.e., $\Gamma_{\rm max}$, as defined in the [previous tutorial notebook](Tutorial-u0_smallb_Poynting-Cartesian.ipynb)).


```python
%%writefile sbPoynETNRPy/param.ccl

shares: HydroBase
USES CCTK_INT timelevels

restricted:
CCTK_INT timelevels "Number of active timelevels" STEERABLE=RECOVER
{
  0:3 :: ""
} 3

restricted:
CCTK_INT sbPoynETNRPy_calc_every "Compute these quantities every sbPoynETNRPy_calc_every iterations." STEERABLE=ALWAYS
{
  *:* :: ""
} 1

restricted:
CCTK_REAL GAMMA_SPEED_LIMIT "Maximum Lorentz factor."
{
 1:* :: "Positive > 1, though you'll likely have troubles in GRMHD far above 10, or far above 2000 in GRFFE."
} 10.0

```

    Overwriting sbPoynETNRPy/param.ccl


Finally, we will write the file `schedule.ccl`; its official documentation is found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). 

This file registers the function we wish to call, `sbPoynETNRPy`, with the Einstein Toolkit scheduler.


```python
%%writefile sbPoynETNRPy/schedule.ccl

STORAGE: smallb4U_group[timelevels]
STORAGE: smallb4_sq_group[timelevels]
STORAGE: PoynSU_group[timelevels]

schedule group sbPoynETNRPy_group in MoL_PseudoEvolution after ADMBase_SetADMVars
{
} "Schedule sbPoynETNRPy group"

schedule sbPoynETNRPy in sbPoynETNRPy_group
{
    LANG: C
    READS: admbase::gxx(Everywhere)
    READS: admbase::gxy(Everywhere)
    READS: admbase::gxz(Everywhere)
    READS: admbase::gyy(Everywhere)
    READS: admbase::gyz(Everywhere)
    READS: admbase::gzz(Everywhere)

    READS: admbase::alpha(Everywhere)

    READS: admbase::betax(Everywhere)
    READS: admbase::betay(Everywhere)
    READS: admbase::betaz(Everywhere)

    READS: HydroBase::vel(Everywhere)
    READS: HydroBase::Bvec(Everywhere)

    WRITES: sbPoynETNRPy::smallb4U0(Everywhere)
    WRITES: sbPoynETNRPy::smallb4U1(Everywhere)
    WRITES: sbPoynETNRPy::smallb4U2(Everywhere)
    WRITES: sbPoynETNRPy::smallb4U3(Everywhere)

    WRITES: sbPoynETNRPy::smallb4_sq(Everywhere)

    WRITES: sbPoynETNRPy::PoynSU0(Everywhere)
    WRITES: sbPoynETNRPy::PoynSU1(Everywhere)
    WRITES: sbPoynETNRPy::PoynSU2(Everywhere)
} "Call sbPoynETNRPy main function, to compute $b^mu$, $b^2$, and $S^i$"
```

    Overwriting sbPoynETNRPy/schedule.ccl


<a id='etksys'></a> 

## Step 2.c: Inform the Einstein Toolkit build system of the C code \[Back to [top](#toc)\]
$$\label{etksys}$$

The `make.code.defn` lists the source files that need to be compiled. Naturally, this thorn has only the one C file $-$ written above $-$ to compile:


```python
%%writefile sbPoynETNRPy/src/make.code.defn

SRCS = sbPoynETNRPy.c
```

    Overwriting sbPoynETNRPy/src/make.code.defn


<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-u0_smallb_Poynting.pdf](Tutorial-ETK_thorn-u0_smallb_Poynting.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-u0_smallb_Poynting")
```

    Created Tutorial-ETK_thorn-u0_smallb_Poynting.tex, and compiled LaTeX file
        to PDF file Tutorial-ETK_thorn-u0_smallb_Poynting.pdf

