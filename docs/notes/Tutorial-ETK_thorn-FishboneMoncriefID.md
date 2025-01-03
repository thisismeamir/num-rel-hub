<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# `FishboneMoncriefID`: An Einstein Toolkit Initial Data Thorn for Fishbone-Moncrief initial data

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook demonstrates the construction of an Einstein Toolkit thorn, converting SymPy expressions for Fishbone-Moncrief initial data into a C-code kernel with NRPy+. It highlights the necessity of adjusting polytropic constants due to rescaling of the disk's maximum density to unity, preserving the polytropic equation of state. The tutorial further delineates the creation of the thorn's interface with the larger Einstein Toolkit infrastructure, including 'ccl' file formation.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** Agrees with trusted Fishbone-Moncrief initial data module in HARM3D. Also generates results in agreement with the trusted version sent to Event Horizon Telescope (EHT) GRMHD code comparison project collaborators. This thorn was used for the [IllinoisGRMHD](http://illinoisgrmhd.net) contribution to the [EHT GRMHD code comparison project](https://arxiv.org/abs/1904.04923).

### NRPy+ Source Code for this module: [FishboneMoncriefID/FishboneMoncriefID.py](../edit/FishboneMoncriefID/FishboneMoncriefID.py) [\[tutorial\]](Tutorial-FishboneMoncriefID.ipynb) Constructs SymPy expressions for [Fishbone-Moncrief initial data](Tutorial-FishboneMoncriefID.ipynb)

## Introduction:
In this part of the tutorial, we will construct an Einstein Toolkit (ETK) thorn (module) that will set up Fishbone-Moncrief initial data. In the [Tutorial-FishboneMoncriefID](Tutorial-FishboneMoncriefID.ipynb) tutorial notebook, we used NRPy+ to construct the SymPy expressions for Fishbone-Moncrief initial data. 

We will construct this thorn in two steps.

1. Call on NRPy+ to convert the SymPy expressions for the initial data into one C-code kernel.
1. Write the C code and linkages to the Einstein Toolkit infrastructure (i.e., the .ccl files) to complete this Einstein Toolkit module.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Call on NRPy+ to convert the SymPy expression for the Fishbone-Moncrief initial data into a C-code kernel
1. [Step 2](#einstein): Interfacing with the Einstein Toolkit
    1. [Step 2.a](#einstein_c): Constructing the Einstein Toolkit C-code calling functions that include the C code kernels
    1. [Step 2.b](#einstein_ccl): CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure
    1. [Step 2.c](#einstein_list): Add the C code to the Einstein Toolkit compilation list
1. [Step 3](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Call on NRPy+ to convert the SymPy expression for the Fishbone-Moncrief initial data into a C-code kernel \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

After importing the core modules, we will set `GridFuncMemAccess` to `ETK`. SymPy expressions for Fishbone-Moncrief initial data are written inside [FishboneMoncriefID/FishboneMoncriefID.py](../edit/FishboneMoncriefID/FishboneMoncriefID.py), and we simply import them for use here.


```python
# Step 1: Call on NRPy+ to convert the SymPy expression for the
#         Fishbone-Moncrief initial data into a C-code kernel

# Step 1a: Import needed NRPy+ core modules:
from outputC import lhrh,outputC # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import os                        # Standard Python modules for multiplatform OS-level functions
import FishboneMoncriefID.FishboneMoncriefID as fmid # Stores closed-form SymPy expressions for F-M initial data.

# Step 1b: This is an Einstein Toolkit (ETK) thorn. Here we
#          tell NRPy+ that gridfunction memory access will
#          therefore be in the "ETK" style.
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

# Step 1c: Within the ETK, the 3D gridfunctions x, y, and z store the
#         Cartesian grid coordinates. Setting the gri.xx[] arrays
#         to point to these gridfunctions forces NRPy+ to treat
#         the Cartesian coordinate gridfunctions properly --
#         reading them from memory as needed.
xcoord,ycoord,zcoord = gri.register_gridfunctions("AUX",["xcoord","ycoord","zcoord"])
gri.xx[0] = xcoord
gri.xx[1] = ycoord
gri.xx[2] = zcoord

# Step 1d: Call the FishboneMoncriefID() function from within the
#          FishboneMoncriefID/FishboneMoncriefID.py module. This
#          sets all the ID gridfunctions.
fmid.FishboneMoncriefID()
Valencia3velocityU = ixp.register_gridfunctions_for_single_rank1("EVOL","Valencia3velocityU")

# -={ Spacetime quantities: Generate C code from expressions and output to file }=-
KerrSchild_to_print = [\
                     lhrh(lhs=gri.gfaccess("out_gfs","alpha"),rhs=fmid.IDalpha),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU0"),rhs=fmid.IDbetaU[0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU1"),rhs=fmid.IDbetaU[1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","betaU2"),rhs=fmid.IDbetaU[2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD00"),rhs=fmid.IDgammaDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD01"),rhs=fmid.IDgammaDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD02"),rhs=fmid.IDgammaDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD11"),rhs=fmid.IDgammaDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD12"),rhs=fmid.IDgammaDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","gammaDD22"),rhs=fmid.IDgammaDD[2][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD00"),rhs=fmid.IDKDD[0][0]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD01"),rhs=fmid.IDKDD[0][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD02"),rhs=fmid.IDKDD[0][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD11"),rhs=fmid.IDKDD[1][1]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD12"),rhs=fmid.IDKDD[1][2]),\
                     lhrh(lhs=gri.gfaccess("out_gfs","KDD22"),rhs=fmid.IDKDD[2][2]),\
                     ]
# Force outCverbose=False for this module to avoid gigantic C files
# filled with the non-CSE expressions.
KerrSchild_CcodeKernel = fin.FD_outputC("returnstring",KerrSchild_to_print,params="outCverbose=False")

# -={ GRMHD quantities: Generate C code from expressions and output to file }=-
FMdisk_GRHD_rho_initial_to_print = [lhrh(lhs=gri.gfaccess("out_gfs","rho_initial"),rhs=fmid.rho_initial)]
FMdisk_GRHD_rho_initial_CcodeKernel = fin.FD_outputC("returnstring",FMdisk_GRHD_rho_initial_to_print)

FMdisk_GRHD_velocities_to_print = [\
                                 lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU0"),rhs=fmid.IDValencia3velocityU[0]),\
                                 lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU1"),rhs=fmid.IDValencia3velocityU[1]),\
                                 lhrh(lhs=gri.gfaccess("out_gfs","Valencia3velocityU2"),rhs=fmid.IDValencia3velocityU[2]),\
                                 ]
FMdisk_GRHD_velocities_CcodeKernel = fin.FD_outputC("returnstring",FMdisk_GRHD_velocities_to_print)

# Step 1f: Create directories for the thorn if they don't exist.
Ccodesdir = "FishboneMoncriefID"
cmd.mkdir(Ccodesdir)
cmd.mkdir(os.path.join(Ccodesdir,"src"))

# Step 1g: Write the C code kernel to file.
with open(os.path.join(Ccodesdir,"src","KerrSchild.h"), "w") as file:
    file.write(str(KerrSchild_CcodeKernel.replace("time","cctk_time")))

with open(os.path.join(Ccodesdir,"src","FMdisk_GRHD_velocities.h"), "w") as file:
    file.write(str(FMdisk_GRHD_velocities_CcodeKernel.replace("time","cctk_time")))

with open(os.path.join(Ccodesdir,"src","FMdisk_GRHD_rho_initial.h"), "w") as file:
    file.write(str(FMdisk_GRHD_rho_initial_CcodeKernel.replace("time","cctk_time")))

hm1string = outputC(fmid.hm1,"hm1",filename="returnstring")
with open(os.path.join(Ccodesdir,"src","FMdisk_GRHD_hm1.h"), "w") as file:
    file.write(str(hm1string))
```

<a id='einstein'></a>

# Step 2: Interfacing with the Einstein Toolkit \[Back to [top](#toc)\]
$$\label{einstein}$$


<a id='einstein_c'></a>

## Step 2.a: Constructing the Einstein Toolkit C-code calling functions that include the C code kernels \[Back to [top](#toc)\]
$$\label{einstein_c}$$

Here we construct `InitialData.c`, which contains C driver functions that pull in the necessary NRPy+ C-code kernels.

First we set up driver routines to specify the Kerr-Schild metric and the Fishbone-Moncrief disk velocity at a given gridpoint.


```python
%%writefile $Ccodesdir/src/InitialData.c
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h> // Needed for rand()

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

// Alias for "vel" vector gridfunction:
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

void FishboneMoncrief_KerrSchild(const cGH* restrict const cctkGH,const CCTK_INT *cctk_lsh,
                                 const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,
                                 const CCTK_REAL *xcoordGF,const CCTK_REAL *ycoordGF,const CCTK_REAL *zcoordGF,
                                 CCTK_REAL *alphaGF,CCTK_REAL *betaU0GF,CCTK_REAL *betaU1GF,CCTK_REAL *betaU2GF,
                                 CCTK_REAL *gammaDD00GF,CCTK_REAL *gammaDD01GF,CCTK_REAL *gammaDD02GF,CCTK_REAL *gammaDD11GF,CCTK_REAL *gammaDD12GF,CCTK_REAL *gammaDD22GF,
                                 CCTK_REAL     *KDD00GF,CCTK_REAL     *KDD01GF,CCTK_REAL     *KDD02GF,CCTK_REAL     *KDD11GF,CCTK_REAL     *KDD12GF,CCTK_REAL     *KDD22GF)
{

  DECLARE_CCTK_PARAMETERS

#include "KerrSchild.h"

}

void FishboneMoncrief_FMdisk_GRHD_velocities(const cGH* restrict const cctkGH,const CCTK_INT *cctk_lsh,
                                             const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,
                                             const CCTK_REAL *xcoordGF,const CCTK_REAL *ycoordGF,const CCTK_REAL *zcoordGF,
                                             CCTK_REAL *Valencia3velocityU0GF, CCTK_REAL *Valencia3velocityU1GF, CCTK_REAL *Valencia3velocityU2GF)
{

  DECLARE_CCTK_PARAMETERS

#include "FMdisk_GRHD_velocities.h"

}
```

    Overwriting FishboneMoncriefID/src/InitialData.c


Next, we set up the driver function for setting all metric and hydrodynamical fields $\rho,P,\epsilon,v^i$.

**Important**: Suppose the Fishbone-Moncrief initial data yield a density $\rho(r,\theta)$ (which is valid for all Fishbone-Moncrief disks centered at the origin, $r=0$, as F-M disks are axisymmetric). Then the disk will have pressure
$$
P = \kappa \rho^\Gamma.
$$

Since the disk is not self-gravitating, we are allowed to rescale the maximum density in the disk to be one in code units; i.e., $\rho_{\rm max}=1$. This may be incompatible with the initial choice of polytropic constant $\kappa$, as rescaling the density results in a rescaling of pressure $P$, as follows.

When we rescale $\rho$ so that the maximum density in the disk is one, we make the following transformation:
$$
\rho \to \rho' = \frac{\rho}{\rho_{\rm max}}.
$$
Since pressure has units of $\rho c^2$, and we use $G=c=1$ units, pressure must therefore be rescaled by the same factor:
\begin{align}
P \to P' &= \frac{P}{\rho_{\rm max}} \\
&= \frac{\kappa \rho^\Gamma}{\rho_{\rm max}} \\
&= \kappa \frac{\rho^\Gamma}{\rho_{\rm max}} \\
&= \kappa \frac{(\rho' \rho_{\rm max})^\Gamma}{\rho_{\rm max}} \\
&= \kappa \rho_{\rm max}^{\Gamma-1} (\rho')^\Gamma \\
&= \kappa' (\rho')^\Gamma
\end{align}

Thus the polytropic equation of state is still valid, but only if 
$$
\kappa' = \kappa \rho_{\rm max}^{\Gamma-1} = \frac{P_{\rm max}}{\rho_{\rm max}}.
$$
As e.g., `IllinoisGRMHD` requires that the initial $P'$ be given as a polytropic equation of state, with $P'_{\rm cold} = \kappa' (\rho')^\Gamma$, $\kappa'$ must be input into the `FishboneMoncriefID` (and `IllinoisGRMHD`) thorns instead of $\kappa$. If this does not happen, the code will error out, providing the correct value for $\kappa'$ that must be set in the parameter file.


```python
%%writefile -a $Ccodesdir/src/InitialData.c

void FishboneMoncrief_ET_GRHD_initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Fishbone-Moncrief Disk Initial data.");
  CCTK_VINFO("Using input parameters of\n a = %e,\n M = %e,\nr_in = %e,\nr_at_max_density = %e\nkappa = %e\ngamma = %e",a,M,r_in,r_at_max_density,kappa,gamma);

  // First compute maximum pressure and density
  CCTK_REAL P_max, rho_max;
  {
    CCTK_REAL hm1;
    CCTK_REAL xcoord = r_at_max_density;
    CCTK_REAL ycoord = 0.0;
    CCTK_REAL zcoord = 0.0;
    {
#include "FMdisk_GRHD_hm1.h"
    }
    rho_max = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) );
    P_max   = kappa * pow(rho_max, gamma);
  }

  // We enforce units such that rho_max = 1.0; if these units are not obeyed, then
  //    we error out. If we did not error out, then the value of kappa used in all
  //    EOS routines would need to be changed, and generally these appear as
  //    read-only parameters.
  if(fabs(P_max/rho_max - kappa) > 1e-8) {
    printf("Error: To ensure that P = kappa*rho^Gamma, where rho_max = 1.0,\n");
    printf("       you must set (in your parfile) the polytropic constant kappa = P_max/rho_max = %.15e\n\n",P_max/rho_max);
    printf(" Needed values for kappa, for common values of Gamma:\n");
    printf(" For Gamma =4/3, use kappa=K_initial=K_poly = 4.249572342020724e-03 to ensure rho_max = 1.0\n");
    printf(" For Gamma =5/3, use kappa=K_initial=K_poly = 6.799315747233158e-03 to ensure rho_max = 1.0\n");
    printf(" For Gamma = 2,  use kappa=K_initial=K_poly = 8.499144684041449e-03 to ensure rho_max = 1.0\n");
    exit(1);
  }

#pragma omp parallel for
  for(CCTK_INT k=0;k<cctk_lsh[2];k++) for(CCTK_INT j=0;j<cctk_lsh[1];j++) for(CCTK_INT i=0;i<cctk_lsh[0];i++) {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL xcoord = x[idx];
        CCTK_REAL ycoord = y[idx];
        CCTK_REAL zcoord = z[idx];
        CCTK_REAL rr = r[idx];

        FishboneMoncrief_KerrSchild(cctkGH,cctk_lsh,
                                    i,j,k,
                                    x,y,z,
                                    alp,betax,betay,betaz,
                                    gxx,gxy,gxz,gyy,gyz,gzz,
                                    kxx,kxy,kxz,kyy,kyz,kzz);

        CCTK_REAL hm1;
        bool set_to_atmosphere=false;
        if(rr > r_in) {
          {
#include "FMdisk_GRHD_hm1.h"
          }
          if(hm1 > 0) {
            rho[idx] = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) ) / rho_max;
            press[idx] = kappa*pow(rho[idx], gamma);
            // P = (\Gamma - 1) rho epsilon
            eps[idx] = press[idx] / (rho[idx] * (gamma - 1.0));
            FishboneMoncrief_FMdisk_GRHD_velocities(cctkGH,cctk_lsh,
                                                    i,j,k,
                                                    x,y,z,
                                                    velx,vely,velz);
          } else {
            set_to_atmosphere=true;
          }
        } else {
          set_to_atmosphere=true;
        }
        // Outside the disk? Set to atmosphere all hydrodynamic variables!
        if(set_to_atmosphere) {
          // Choose an atmosphere such that
          //   rho =       1e-5 * r^(-3/2), and
          //   P   = k rho^gamma
          // Add 1e-100 or 1e-300 to rr or rho to avoid divisions by zero.
          rho[idx] = 1e-5 * pow(rr + 1e-100,-3.0/2.0);
          press[idx] = kappa*pow(rho[idx], gamma);
          eps[idx] = press[idx] / ((rho[idx] + 1e-300) * (gamma - 1.0));
          w_lorentz[idx] = 1.0;
          velx[idx] = 0.0;
          vely[idx] = 0.0;
          velz[idx] = 0.0;
        }
      }

  CCTK_INT final_idx = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,cctk_lsh[1]-1,cctk_lsh[2]-1);
  CCTK_VINFO("=====   OUTPUTS   =====");
  CCTK_VINFO("betai: %e %e %e \ngij: %e %e %e %e %e %e \nKij: %e %e %e %e %e %e\nalp: %e\n",betax[final_idx],betay[final_idx],betaz[final_idx],gxx[final_idx],gxy[final_idx],gxz[final_idx],gyy[final_idx],gyz[final_idx],gzz[final_idx],kxx[final_idx],kxy[final_idx],kxz[final_idx],kyy[final_idx],kyz[final_idx],kzz[final_idx],alp[final_idx]);
  CCTK_VINFO("rho: %.15e\nPressure: %.15e\nvx: %.15e\nvy: %.15e\nvz: %.15e",rho[final_idx],press[final_idx],velx[final_idx],vely[final_idx],velz[final_idx]);
}

void FishboneMoncrief_ET_GRHD_initial__perturb_pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for(CCTK_INT k=0;k<cctk_lsh[2];k++) for(CCTK_INT j=0;j<cctk_lsh[1];j++) for(CCTK_INT i=0;i<cctk_lsh[0];i++) {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
        // Generate random number in range [0,1),
        // snippet courtesy http://daviddeley.com/random/crandom.htm
        CCTK_REAL random_number_between_0_and_1 = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        CCTK_REAL random_number_between_min_and_max = random_min + (random_max - random_min)*random_number_between_0_and_1;
        press[idx] = press[idx]*(1.0 + random_number_between_min_and_max);
        // Add 1e-300 to rho to avoid division by zero when density is zero.
        eps[idx] = press[idx] / ((rho[idx] + 1e-300) * (gamma - 1.0));
      }
}
```

    Appending to FishboneMoncriefID/src/InitialData.c


<a id='einstein_ccl'></a>

## Step 2.b: CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{einstein_ccl}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn:

1. `interface.ccl}`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns. Specifically, this file governs the interaction between this thorn and others; more information can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 
With "implements", we give our thorn its unique name. By "inheriting" other thorns, we tell the Toolkit that we will rely on variables that exist and are declared "public" within those functions.


```python
%%writefile $Ccodesdir/interface.ccl
implements: FishboneMoncriefID
inherits: admbase grid hydrobase
```

    Overwriting FishboneMoncriefID/interface.ccl


2. `param.ccl`: specifies free parameters within the thorn, enabling them to be set at runtime. It is required to provide allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3).


```python
%%writefile $Ccodesdir/param.ccl
shares: grid
shares: ADMBase
USES CCTK_INT lapse_timelevels
USES CCTK_INT shift_timelevels
USES CCTK_INT metric_timelevels

USES KEYWORD metric_type

EXTENDS KEYWORD initial_data
{
  "FishboneMoncriefID" :: "Initial data from FishboneMoncriefID solution"
}
EXTENDS KEYWORD initial_lapse
{
  "FishboneMoncriefID" :: "Initial lapse from FishboneMoncriefID solution"
}
EXTENDS KEYWORD initial_shift
{
  "FishboneMoncriefID" :: "Initial shift from FishboneMoncriefID solution"
}
EXTENDS KEYWORD initial_dtlapse
{
  "FishboneMoncriefID" :: "Initial dtlapse from FishboneMoncriefID solution"
}
EXTENDS KEYWORD initial_dtshift
{
  "FishboneMoncriefID" :: "Initial dtshift from FishboneMoncriefID solution"
}

shares: HydroBase
EXTENDS KEYWORD initial_hydro
{
  "FishboneMoncriefID" :: "Initial GRHD data from FishboneMoncriefID solution"
}

#["r_in","r_at_max_density","a","M"] A_b, kappa, gamma
restricted:
CCTK_REAL r_in "Fixes the inner edge of the disk"
{
 0.0:* :: "Must be positive"
} 6.0

restricted:
CCTK_REAL r_at_max_density "Radius at maximum disk density. Needs to be > r_in"
{
 0.0:* :: "Must be positive"
} 12.0

restricted:
CCTK_REAL a "The spin parameter of the black hole"
{
 0:1.0 :: "Positive values, up to 1. Negative disallowed, as certain roots are chosen in the hydro fields setup. Check those before enabling negative spins!"
} 0.9375

restricted:
CCTK_REAL M "Kerr-Schild BH mass. Probably should always set M=1."
{
 0.0:* :: "Must be positive"
} 1.0

restricted:
CCTK_REAL A_b "Scaling factor for the vector potential"
{
 *:* :: ""
} 1.0

restricted:
CCTK_REAL kappa "Equation of state: P = kappa * rho^gamma"
{
 0.0:* :: "Positive values"
} 1.0e-3

restricted:
CCTK_REAL gamma "Equation of state: P = kappa * rho^gamma"
{
 0.0:* :: "Positive values"
} 1.3333333333333333333333333333

##################################
# PRESSURE PERTURBATION PARAMETERS
private:
CCTK_REAL random_min "Floor value of random perturbation to initial pressure, where perturbed pressure = pressure*(1.0 + (random_min + (random_max-random_min)*RAND[0,1)))"
{
  *:*      :: "Any value"
} -0.02

private:
CCTK_REAL random_max "Ceiling value of random perturbation to initial pressure, where perturbed pressure = pressure*(1.0 + (random_min + (random_max-random_min)*RAND[0,1)))"
{
  *:*      :: "Any value"
} 0.02


```

    Overwriting FishboneMoncriefID/param.ccl


3. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions. $\text{schedule.ccl}$'s official documentation may be found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). 

We specify here the standardized ETK "scheduling bins" in which we want each of our thorn's functions to run.


```python
%%writefile $Ccodesdir/schedule.ccl
STORAGE: ADMBase::metric[metric_timelevels], ADMBase::curv[metric_timelevels], ADMBase::lapse[lapse_timelevels], ADMBase::shift[shift_timelevels]

schedule FishboneMoncrief_ET_GRHD_initial IN HydroBase_Initial
{
  LANG: C
  READS: grid::x(Everywhere)
  READS: grid::y(Everywhere)
  READS: grid::z(Everywhere)
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
  WRITES: hydrobase::vel(Everywhere) # Note that vel is a vector gridfunction.
  WRITES: hydrobase::rho(Everywhere)
  WRITES: hydrobase::eps(Everywhere)
  WRITES: hydrobase::press(Everywhere)
} "Set up general relativistic hydrodynamic (GRHD) fields for Fishbone-Moncrief disk"

schedule FishboneMoncrief_ET_GRHD_initial__perturb_pressure IN CCTK_INITIAL AFTER Seed_Magnetic_Fields BEFORE IllinoisGRMHD_ID_Converter
{
    LANG: C
} "Add random perturbation to initial pressure, after seed magnetic fields have been set up (in case we'd like the seed magnetic fields to depend on the pristine pressures)"

```

    Overwriting FishboneMoncriefID/schedule.ccl


<a id='einstein_list'></a>

## Step 2.c: Add the C code to the Einstein Toolkit compilation list \[Back to [top](#toc)\]
$$\label{einstein_list}$$

We will also need `make.code.defn`, which indicates the list of files that need to be compiled. This thorn only has the one C file to compile.


```python
%%writefile $Ccodesdir/src/make.code.defn
SRCS = InitialData.c
```

    Overwriting FishboneMoncriefID/src/make.code.defn


<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-FishboneMoncriefID.pdf](Tutorial-ETK_thorn-FishboneMoncriefID.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-FishboneMoncriefID")
```

    Created Tutorial-ETK_thorn-FishboneMoncriefID.tex, and compiled LaTeX file
        to PDF file Tutorial-ETK_thorn-FishboneMoncriefID.pdf

