<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# `interp_sph_grids_ETK`: An Einstein Toolkit Module for Interpolation to Spherical Grids

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This module is designed to interpolate arbitrary quantities on [Einstein Toolkit](https://einsteintoolkit.org/) Adaptive-Mesh Refinement (AMR) grids (using the [Carpet](https://carpetcode.org/) AMR infrastructure) to numerical grids with spherical sampling. The tutorial elaborates on core ETK C routines for interpolation, NRPy+ output for required gridfunctions, and the necessary CCL files for interfacing with the Toolkit.

**Notebook Status:** <font color='red'><b> In progress </b></font>

**Validation Notes:** This module has not yet undergone validation testing.

## Introduction:
Given some set of $N$ quantities $\mathbf{Q}=\{Q_0,Q_1,Q_2,...,Q_{N-2},Q_{N-1}\}$, this module performs the following for each $Q_i$:

1. Evaluate $Q_i$ at all gridpoints that are not ghost zones. Sometimes $Q_i$ is computed using finite difference derivatives, so this is necessary.
1. Call upon Carpet's interpolation and interprocessor synchronization functions to fill in $Q_i$ at all ghost zones, *except* at the outer boundary. We do not generally trust $Q_i$ at the outer boundary due to errors associated with the approximate outer boundary conditions. 
1. At this point, $Q_i$ is set at all gridpoints except ghost zones at the outer boundary. Interpolate $Q_i$ to the spherical grids, **maintaining the Cartesian basis for all vectors and tensors**, and append the result to a file.

This tutorial notebook takes a three-part structure. First, all the needed core Einstein Toolkit (ETK) C routines for interpolation are presented. Second, NRPy+ is used to output gridfunctions needed on the spherical grids. Third, the needed files for interfacing this module with the rest of the Einstein Toolkit (ccl files) are specified.

<a id='toc'></a>

# Table of Contents: 
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#etkmodule): Setting up the Core C Code for the Einstein Toolkit Module
    1. [Step 1.a](#etk_interp): Low-Level ETK Interpolation Function
    1. [Step 1.b](#sphericalgridnotes): Setting up the Spherical Grids
    1. [Step 1.c](#fileformat): Outputting to File
    1. [Step 1.d](#maininterpolator) The Main Interpolation Driver Function
1. [Step 2](#nrpy): Use NRPy+ C Output to Set All Output Gridfunctions
    1. [Step 2.a](#nrpy_list_of_funcs_interp): Set up NRPy-based `list_of_functions_to_interpolate.h` 
        1. [Step 2.a.i](#nrpygrmhd): GRMHD quantities (***IN PROGRESS***)
        1. [Step 2.a.ii](#nrpy4metric): Compute all 10 components of the 4-metric $g_{\mu\nu}$
        1. [Step 2.a.iii](#nrpy4christoffels): Compute all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$
    1. [Step 2.b](#nrpy_c_callingfunction): C code calling function for the NRPy+ C output
    1. [Step 2.c](#nrpygetgfinterporder): `get_gf_name()` and `get_interp_order()` functions
    1. [Step 2.d](#nrpy_interp_counter): C Code for Initializing and incrementing `InterpCounter`
1. [Step 3](#cclfiles): Interfacing with the rest of the Einstein Toolkit; Setting up CCL files
    1. [Step 3.a](#makecodedefn): `make.code.defn`
    1. [Step 3.b](#interfaceccl): `interface.ccl`
    1. [Step 3.c](#paramccl): `param.ccl`
    1. [Step 3.d](#scheduleccl): `schedule.ccl`
1. [Step 4](#readingoutputfile): Python Script for Reading the Output File
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='etkmodule'></a>

# Step 1: Setting up the Core C Code for the Einstein Toolkit Module \[Back to [top](#toc)\]
$$\label{etkmodule}$$

First we set up the output directories for the ETK module:


```python
import cmdline_helper as cmd  # NRPy+: Multi-platform Python command-line interface
import shutil, os             # Standard Python modules for multiplatform OS-level functions, benchmarking

# Create C code output directory:
Ccodesdir = "interp_sph_grids_ETK"
# First remove C code output directory and all subdirectories if they exist
# Courtesy https://stackoverflow.com/questions/303200/how-do-i-remove-delete-a-folder-that-is-not-empty
shutil.rmtree(Ccodesdir, ignore_errors=True)
# Then create a fresh directory
cmd.mkdir(Ccodesdir)
cmd.mkdir(os.path.join(Ccodesdir,"src/"))
```

<a id='etk_interp'></a>

## Step 1.a: Low-Level ETK Interpolation Function \[Back to [top](#toc)\]
$$\label{etk_interp}$$

We start by writing the low-level interpolation function **`Interpolate_to_sph_grid()`**, which  to file. 

**`Interpolate_to_sph_grid()`** takes as input

* **cctkGH**: Information about the underlying Cactus/Carpet grid hierarchy.
* **interp_num_points**: Number of destination interpolation points
* **point_x_temp, point_y_temp, point_z_temp**: Cartesian $(x,y,z)$ location for each of the **interp_num_points** interpolation points.
* **input_array_names[1]**: List of input gridfunction names to interpolate. We will do this only one gridfunction at a time, for gridfunction $Q_i$, as described above.

**`Interpolate_to_sph_grid()`** outputs:

* **output_f[1]**: The gridfunction **input_array_names[1]** interpolated to the set of **interp_num_points** specified in the input.


```python
%%writefile $Ccodesdir/src/Interpolate_to_sph_grid.h

void Interpolate_to_sph_grid(cGH *cctkGH,CCTK_INT interp_num_points, CCTK_INT interp_order,
                             CCTK_REAL *point_x_temp,CCTK_REAL *point_y_temp,CCTK_REAL *point_z_temp,
                             const CCTK_STRING input_array_names[1], CCTK_REAL *output_f[1]) {
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr;

  const CCTK_INT NUM_INPUT_ARRAYS=1;
  const CCTK_INT NUM_OUTPUT_ARRAYS=1;

  CCTK_STRING coord_system = "cart3d";

  // Set up handles
  const CCTK_INT coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);
  }

  const CCTK_INT operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  char interp_order_string[10];
  snprintf(interp_order_string, 10, "order=%d", interp_order);
  CCTK_STRING interpolator_pars = interp_order_string;
  CCTK_INT param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  }

  CCTK_INT operand_indices[NUM_INPUT_ARRAYS]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {
    operand_indices[i] = i;
  }
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");


  CCTK_INT operation_codes[NUM_INPUT_ARRAYS];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS  ; i++) {
    operation_codes[i] = 0;
  }
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operation_codes, "operation_codes");

  const void* interp_coords[3]
    = { (const void *) point_x_temp,
        (const void *) point_y_temp,
        (const void *) point_z_temp };

  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "COULD NOT FIND VARIABLE '%s'.",
        input_array_names[i]);
      exit(1);
    }
  }

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS ; i++) {
    output_array_types[i] = CCTK_VARIABLE_REAL;
  }

  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) output_f[0] };

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               3, // number of dimensions
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_num_points,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS, // Number of output arrays
                               output_array_types,
                               output_arrays);
  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    exit(1);
  }
}
```

    Writing interp_sph_grids_ETK/src/Interpolate_to_sph_grid.h


<a id='sphericalgridnotes'></a>

## Step 1.b: Setting up the Spherical Grids \[Back to [top](#toc)\]
$$\label{sphericalgridnotes}$$

+ By default, we set logarithmic radial coordinates: $r(x_{0,i}) = R_0 + e^{x_{0,i}}$, where

    + $x_{0,i} = x_{0, \mathrm{beg}} + \left(i+\frac{1}{2}\right) \Delta x_0$
    + $x_{0, {\mathrm{beg}}} = \log\left( R_{\mathrm{in}} - R_0 \right)$
    + $\Delta x_0 = \frac{1}{N_0}\log\left(\frac{R_\mathrm{out} - R_0}{R_\mathrm{in} - R_0}\right)$


+ As for the polar angle $\theta$, there are two options:
    + **Option 1**: 
  $$ \theta(x_{1,j})  \, = \, \theta_c \, + \, \left( \pi - 2 \theta_c \right) x_{1,j} \, + \, \xi \, \sin\left(2 \pi x_{1,j} \right),$$
      where
        + $x_{1,j} = x_{1, \mathrm{beg}} + \left(j+\frac{1}{2}\right) \Delta x_1$
        + $\Delta x_1 = \frac{1}{N_1}$

    + **Option 2**: 
  $$ \theta(x_{1,j}) = \frac{\pi}{2} \left[  1  + \left(1-\xi \right) \left(2 x_{1,j} - 1 \right) + \left( \xi - \frac{2 \theta_c}{\pi} \right) \left( 2 x_{1,j} - 1 \right)^n \right],$$
      where
        + $n$ is odd
        + $x_{1,j} = x_{1, \mathrm{beg}} + \left(j+\frac{1}{2}\right) \Delta x_1$
        + $\Delta x_1 = \frac{1}{N_1}$


+ The azimuthal angle $\phi$ is uniform, so that $\phi(x_{2,k}) = x_{2,k}$:
  
    + $x_{2,k} \in [0,2\pi]$
    + $x_{2,k} = x_{2, \mathrm{beg}} + \left(k+\frac{1}{2}\right)\Delta x_{2}$
    + $\Delta x_{2} = \frac{ 2 \pi }{N_2}$


```python
%%writefile $Ccodesdir/src/Set_up_interp_points_on_sph_grid.h

void sph_grid_Interpolate_many_pts__set_interp_pts(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dx0 = log( (Rout - R0) / (Rin - R0) ) / ((CCTK_REAL)N0);
  CCTK_REAL dx1 =      1.0 / ((CCTK_REAL)N1);
  CCTK_REAL dx2 = 2.0*M_PI / ((CCTK_REAL)N2);
  CCTK_REAL x0_beg = log( Rin - R0 );
  CCTK_INT which_pt = 0;
  for(CCTK_INT k=0;k<N2;k++) for(CCTK_INT j=0;j<N1;j++) for(CCTK_INT i=0;i<N0;i++) {
    CCTK_REAL x0_i = x0_beg + ((CCTK_REAL)i + 0.5)*dx0;
    CCTK_REAL rr = R0 + exp(x0_i);

    CCTK_REAL x1_j = x1_beg + ((CCTK_REAL)j + 0.5)*dx1;
    CCTK_REAL th = -1e300;
    if(theta_option == 1) {
       th = th_c + (M_PI - 2.0*th_c)*x1_j + xi*sin(2.0*M_PI*x1_j);
    } else if (theta_option == 2) {
       th = M_PI/2.0 * ( 1.0 + (1.0 - xi)*(2.0*x1_j - 1.0) + (xi - 2.0*th_c/M_PI)*pow(2.0*x1_j - 1.0 ,th_n) );
    } else {
       printf("Error: theta_option = %d NOT SUPPORTED.",theta_option);
       exit(1);
    }

    CCTK_REAL x2_k = x2_beg + ((CCTK_REAL)k + 0.5)*dx2;
    CCTK_REAL ph = x2_k;

    points_x[which_pt] = rr*sin(th)*cos(ph);
    points_y[which_pt] = rr*sin(th)*sin(ph);
    points_z[which_pt] = rr*cos(th);
    which_pt++;
  }
}
```

    Writing interp_sph_grids_ETK/src/Set_up_interp_points_on_sph_grid.h


<a id='fileformat'></a>

## Step 1.c: Outputting to File (File format notes) \[Back to [top](#toc)\]
$$\label{fileformat}$$

Since they take almost no space relative to the data chunks, we attach the entire metadata to each interpolated function that is output:


```python
%%writefile $Ccodesdir/src/output_to_file.h

#include "define_NumInterpFunctions.h"

void output_to_file(CCTK_ARGUMENTS,char gf_name[100],int *order,CCTK_REAL *output_f[1]) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char filename[100];
  sprintf (filename, "%s/interp_sph_grids.dat", out_dir);
  FILE *file;
  if(*InterpCounter == 1) {
    file = fopen (filename,"w");
  } else {
    file = fopen (filename,"a+");
  }
  if (! file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "interp_sph_grid__ET_thorn: Cannot open output file '%s'", filename);
    exit(1);
  }

  fwrite(gf_name, 100*sizeof(char), 1, file);
  fwrite(order, sizeof(CCTK_INT), 1, file);

  fwrite(&N0, sizeof(CCTK_INT), 1, file);
  fwrite(&R0, sizeof(CCTK_REAL), 1, file);
  fwrite(&Rin, sizeof(CCTK_REAL), 1, file);
  fwrite(&Rout, sizeof(CCTK_REAL), 1, file);

  fwrite(&N1, sizeof(CCTK_INT), 1, file);
  fwrite(&x1_beg, sizeof(CCTK_REAL), 1, file);
  fwrite(&theta_option, sizeof(CCTK_INT), 1, file);
  fwrite(&th_c, sizeof(CCTK_REAL), 1, file);
  fwrite(&xi, sizeof(CCTK_REAL), 1, file);
  fwrite(&th_n, sizeof(CCTK_INT), 1, file);

  fwrite(&N2, sizeof(CCTK_INT), 1, file);
  fwrite(&x2_beg, sizeof(CCTK_REAL), 1, file);

  CCTK_REAL magic_number = 1.130814081305130e-21;
  fwrite(&magic_number, sizeof(CCTK_REAL), 1, file);
  fwrite(&cctk_iteration, sizeof(CCTK_INT), 1, file);
  fwrite(&cctk_time, sizeof(CCTK_REAL), 1, file);
  for(CCTK_INT i=0;i<1;i++) {
    fwrite(output_f[i], sizeof(CCTK_REAL)*N0*N1*N2, 1, file);
  }

  fclose(file);
}
```

    Writing interp_sph_grids_ETK/src/output_to_file.h


<a id='maininterpolator'></a>

## Step 1.d: The Main Interpolation Driver Function \[Back to [top](#toc)\]
$$\label{maininterpolator}$$

The **`Interpolate_to_sph_grid_main_function()`** function calls the above functions as follows:

1. **`sph_grid_Interpolate_many_pts__set_interp_pts()`**: First set up the spherical grids
1. **`Interpolate_to_sph_grid()`**: Output


```python
%%writefile $Ccodesdir/src/main_function.cc

// Include needed ETK & C library header files:
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// Needed for dealing with Cactus/ETK infrastructure
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
// Needed for low-level interpolation functions
#include "util_Table.h"
#include "util_String.h"

// Include locally-defined C++ functions:
#include "Set_up_interp_points_on_sph_grid.h"
#include "Interpolate_to_sph_grid.h"
#include "output_to_file.h"
#include "get_interp_order.h"
#include "get_gf_name.h"

void Interpolate_to_sph_grid_main_function(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Perform interpolation only at iteration == interp_out_iteration:
  if(cctk_iteration != interp_out_iteration) return;

  // Set up spherically sampled interpolation grid arrays points_x,points_y,points_z:
  sph_grid_Interpolate_many_pts__set_interp_pts(CCTK_PASS_CTOC);

  // Set up output array:
  CCTK_REAL *output_f[1];
  output_f[0] = output_interped;
  // The name of the input gridfunction is always "interp_sph_grids_ETK::interped_gf":
  const CCTK_STRING input_array_names[1] = { "interp_sph_grids_ETK::interped_gf" };

  // Perform interpolation!
  int order = get_interp_order(*InterpCounter);
  char gf_name[100];
  get_gf_name(*InterpCounter,gf_name);
  printf("Interpolating\033[1m %s \033[0m... using interpolation order = %d\n",gf_name,order);
  Interpolate_to_sph_grid(cctkGH, N0*N1*N2, order,
                             points_x,points_y,points_z, input_array_names, output_f);

  if(CCTK_MyProc(cctkGH)==0) {
    for(int i=0;i<N0*N1*N2;i++) {
        if(output_f[0][i] > 1e20) {
            printf("BAD POINT: %s %d %e %e %e %e\n",gf_name,i,points_x[i],points_y[i],points_z[i], output_f[0][i]);
//            exit(0);
        }
    }
    output_to_file(CCTK_PASS_CTOC,gf_name,&order,output_f);
    printf("Interpolate_to_sph_grid_main_function(): Just output to file at iteration %d\n",cctk_iteration);
  } else {
    printf("Interpolate_to_sph_grid_main_function(): Process !=0 waiting for file output at iteration %d\n",cctk_iteration);
  }
}
```

    Writing interp_sph_grids_ETK/src/main_function.cc


<a id='nrpy'></a>

# Step 2: Use NRPy+ C Output to Set All Output Gridfunctions \[Back to [top](#toc)\]
$$\label{nrpy}$$


```python
# Step 2 P0: Import needed NRPy+ modules/parameters
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
from outputC import lhrh         # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import loop as lp                # NRPy+: Generate C code loops

par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

from collections import namedtuple
gf_interp = namedtuple('gf_interp', 'description order')
gf_interp_list = []
gf_interp_list.append(gf_interp("dummy -- used because this is a 1-offset array",order="-100"))

interped_gf = gri.register_gridfunctions("AUX","interped_gf")

def interp_fileout(which_InterpCounter, expression, filename):
    kernel = fin.FD_outputC("returnstring",lhrh(lhs=gri.gfaccess("out_gfs","interped_gf"),rhs=expression),"outCverbose=False")
    output_type="a"
    if which_InterpCounter == 1:
        output_type="w"

    with open(filename, output_type) as file:
        file.write("if(*InterpCounter == "+str(which_InterpCounter)+") {\n")
        file.write(lp.loop(["i2","i1","i0"],
                           ["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],\
                           ["cctk_lsh[2]-cctk_nghostzones[2]",
                            "cctk_lsh[1]-cctk_nghostzones[1]",
                            "cctk_lsh[0]-cctk_nghostzones[0]"],\
                           ["1","1","1"],\
                           ["#pragma omp parallel for","",""],"   ",kernel))
        file.write("}\n")
    # If successful, return incremented which_InterpCounter:
    return which_InterpCounter+1
```

<a id='nrpy_list_of_funcs_interp'></a>

## Step 2.a:  Set up NRPy-based `list_of_functions_to_interpolate.h` \[Back to [top](#toc)\]
$$\label{nrpy_list_of_funcs_interp}$$

First specify NRPy+ output file and initialize `which_InterpCounter`, which keeps track of the number of interpolated functions on the grid


```python
# Step 2.a:  Set up NRPy-based "list_of_functions_to_interpolate.h"

NRPyoutfilename = os.path.join(Ccodesdir,"src","list_of_functions_to_interpolate.h")

which_InterpCounter = 1
```

<a id='nrpygrmhd'></a>

### Step 2.a.i: GRMHD quantities (*IN PROGRESS*) \[Back to [top](#toc)\]
$$\label{nrpygrmhd}$$

These include

* $\rho_b$, the baryonic density (i.e., the HydroBase variable $\verb|rho|$)
* $P$, the total gas pressure (i.e., the HydroBase variable $\verb|press|$)
* $\Gamma v_{(n)}^i$, the Valencia 3-velocity times the Lorentz factor (i.e., the HydroBase 3-gridfuntion $\verb|vel|$, multiplied by the Lorentz factor). This definition of velocity has the advantage that after interpolation, it will not violate $u^\mu u_\mu = -1$. In terms of the IllinoisGRMHD 3-velocity $v^i = u^i / u^0$, the Valencia 3-velocity is given by (Eq. 11 of [Etienne *et al*](https://arxiv.org/pdf/1501.07276.pdf)):

$$
v_{(n)}^i = \frac{1}{\alpha} \left(v^i + \beta^i\right).
$$
Further, $\Gamma = \alpha u^0$ is given by (as shown [here](Tutorial-u0_smallb_Poynting-Cartesian.ipynb)):
$$
\Gamma = \alpha u^0 = \sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}.
$$
Therefore, $\Gamma v_{(n)}^i$ is given by
$$
\Gamma v_{(n)}^i = \frac{1}{\alpha} \left(v^i + \beta^i\right) \sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}.
$$

* $A_i$, the *unstaggered* magnetic vector potential.
* $B^i$, the *unstaggered* magnetic field vector (output only for validation purposes).


```python
# Step 2.a.i: GRMHD quantities (*IN PROGRESS*)
# INPUT GRIDFUNCTIONS: The AUX or EVOL designation is *not* used in diagnostic modules.
gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01")
betaU = ixp.register_gridfunctions_for_single_rank1("AUX","betaU")
alpha = gri.register_gridfunctions("AUX","alpha")

DIM=3

gf_interp_list.append(gf_interp("IGM density primitive",order="1"))
rho_b       = gri.register_gridfunctions("AUX","rho_b")
interp_expr = rho_b
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("IGM pressure primitive",order="1"))
P = gri.register_gridfunctions("AUX","P")
interp_expr = P
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

Next we implement:
$$
v_{(n)}^i = \frac{1}{\alpha} \left(v^i + \beta^i\right),
$$
and
$$
\Gamma v_{(n)}^i = \sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}} v_{(n)}^i.
$$


```python
IGMvU = ixp.register_gridfunctions_for_single_rank1("AUX","IGMvU")
Valenciav = ixp.zerorank1()
for i in range(DIM):
    Valenciav[i] = 1/alpha * (IGMvU[i] + betaU[i])
v_dot_v = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        v_dot_v += gammaDD[i][j]*Valenciav[i]*Valenciav[j]

Gamma_times_ValenciavU = ixp.zerorank1()
for i in range(DIM):
    Gamma_times_ValenciavU[i] = sp.sqrt(1/(1 - v_dot_v))*Valenciav[i]
    gf_interp_list.append(gf_interp("Lorentz factor, times Valencia vU"+str(i),order="1"))
    interp_expr = Gamma_times_ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# For testing:
# gf_interp_list.append(gf_interp("Lorentz factor",order="1"))
# interp_expr = v_dot_v
# which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("Valencia vU"+str(i),order="1"))
#     interp_expr = Valenciav[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU")
for i in range(DIM):
    gf_interp_list.append(gf_interp("IGM magnetic field component B"+str(i),order="1"))
    interp_expr = BU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy4metric'></a>

### Step 2.a.ii: Compute all 10 components of the 4-metric $g_{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{nrpy4metric}$$

We are given $\gamma_{ij}$, $\alpha$, and $\beta^i$ from ADMBase, and the 4-metric is given in terms of these quantities as
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$


```python
# Step 2.a.ii: Compute all 10 components of the 4-metric g_{\mu\nu}
# Eq. 2.121 in B&S
betaD = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        betaD[i] += gammaDD[i][j]*betaU[j]

# Now compute the beta contraction.
beta2 = sp.sympify(0)
for i in range(DIM):
    beta2 += betaU[i]*betaD[i]

# Eq. 2.122 in B&S
g4DD = ixp.zerorank2(DIM=4)
g4DD[0][0] = -alpha**2 + beta2
for i in range(DIM):
    g4DD[i+1][0] = g4DD[0][i+1] = betaD[i]
for i in range(DIM):
    for j in range(DIM):
        g4DD[i+1][j+1] = gammaDD[i][j]

for mu in range(4):
    for nu in range(mu,4):
        gf_interp_list.append(gf_interp("4-metric component g4DD"+str(mu)+str(nu),order="4"))
        interp_expr = g4DD[mu][nu]
        which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy4christoffels'></a>

### Step 2.a.iii: Compute all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$ \[Back to [top](#toc)\]
$$\label{nrpy4christoffels}$$

By definition,
$$
\Gamma^{\mu}_{\nu\delta} = \frac{1}{2} g^{\mu\eta} \left(g_{\eta\nu,\delta} + g_{\eta\delta,\nu} - g_{\nu\delta,\eta}  \right)
$$

Recall that $g_{\mu\nu}$ is given from $\gamma_{ij}$, $\alpha$, and $\beta^i$ via
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$

The derivatives $g_{\mu\nu,\eta}$ are then computed in terms of finite-difference derivatives of the input ADM gridfunctions $\gamma_{ij}$, $\alpha$, and $\beta^i$, **assuming that the 4-metric is static, so that $\partial_t g_{\mu\nu}=0$ for all $\mu$ and $\nu$**.

To compute $g^{\mu\nu}$, we use the standard formula (Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf)):
$$
g^{\mu\nu} = \begin{pmatrix} 
-\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
\frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
\end{pmatrix},
$$
where $\gamma^{ij}$ is given by the inverse of $\gamma_{ij}$.


```python
# Step 2.a.iii: Compute all 40 4-Christoffels \Gamma^{\mu}_{\nu\delta}

betaDdD = ixp.zerorank2()
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
            betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

# Eq. 2.122 in B&S
g4DDdD = ixp.zerorank3(DIM=4)
alpha_dD   = ixp.declarerank1("alpha_dD")
for i in range(DIM):
    # Recall that g4DD[0][0] = -alpha^2 + betaU[i]*betaD[i]
    g4DDdD[0][0][i+1] += -2*alpha*alpha_dD[i]
    for j in range(DIM):
        g4DDdD[0][0][i+1] += betaU_dD[j][i]*betaD[j] + betaU[j]*betaDdD[j][i]

for i in range(DIM):
    for j in range(DIM):
        # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
        g4DDdD[i+1][0][j+1] = g4DDdD[0][i+1][j+1] = betaDdD[i][j]
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that g4DD[i][j] = gammaDD[i][j]
            g4DDdD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]

gammaUU, dummyDET = ixp.symm_matrix_inverter3x3(gammaDD)

g4UU = ixp.zerorank2(DIM=4)
g4UU[0][0] = -1 / alpha**2
for i in range(DIM):
    g4UU[0][i+1] = g4UU[i+1][0] = betaU[i]/alpha**2
for i in range(DIM):
    for j in range(DIM):
        g4UU[i+1][j+1] = gammaUU[i][j] - betaU[i]*betaU[j]/alpha**2
```

Again, we are to compute:
$$
\Gamma^{\mu}_{\nu\delta} = \frac{1}{2} g^{\mu\eta} \left(g_{\eta\nu,\delta} + g_{\eta\delta,\nu} - g_{\nu\delta,\eta}  \right)
$$


```python
Gamma4UDD = ixp.zerorank3(DIM=4)
for mu in range(4):
    for nu in range(4):
        for delta in range(4):
            for eta in range(4):
                Gamma4UDD[mu][nu][delta] += sp.Rational(1,2)*g4UU[mu][eta]*\
                (g4DDdD[eta][nu][delta] + g4DDdD[eta][delta][nu] - g4DDdD[nu][delta][eta])

# Now output the 4-Christoffels to file:
for mu in range(4):
    for nu in range(4):
        for delta in range(nu,4):
            gf_interp_list.append(gf_interp("4-Christoffel GammaUDD"+str(mu)+str(nu)+str(delta),order="4"))
            interp_expr = Gamma4UDD[mu][nu][delta]
            which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy_c_callingfunction'></a>

## Step 2.b: C code calling function for the NRPy+ C output \[Back to [top](#toc)\]
$$\label{nrpy_c_callingfunction}$$

In the above blocks, we wrote and appended to a file `list_of_functions_to_interpolate.h`. Here we write the calling function for this C code.


```python
%%writefile $Ccodesdir/src/construct_function_to_interpolate__store_to_interped_gf.cc
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Set the gridfunction interped_gf, according to the interpolation counter variable interp_counter.
//    For example, we might interpolate "IllinoisGRMHD::rho_b" if interp_counter==0. The following
//    function takes care of these
void list_of_functions_to_interpolate(cGH *cctkGH,const CCTK_INT *cctk_lsh,const CCTK_INT *cctk_nghostzones,
                                     const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,
                                     const CCTK_INT *InterpCounter,
                                     const CCTK_REAL *rho_bGF,const CCTK_REAL *PGF,
                                     const CCTK_REAL *IGMvU0GF,const CCTK_REAL *IGMvU1GF,const CCTK_REAL *IGMvU2GF,
                                     const CCTK_REAL *BU0GF,const CCTK_REAL *BU1GF,const CCTK_REAL *BU2GF,
                                     const CCTK_REAL *gammaDD00GF,const CCTK_REAL *gammaDD01GF,const CCTK_REAL *gammaDD02GF,
                                     const CCTK_REAL *gammaDD11GF,const CCTK_REAL *gammaDD12GF,const CCTK_REAL *gammaDD22GF,
                                     const CCTK_REAL *betaU0GF,const CCTK_REAL *betaU1GF,const CCTK_REAL *betaU2GF,
                                     const CCTK_REAL *alphaGF,   CCTK_REAL *interped_gfGF) {
#include "list_of_functions_to_interpolate.h"
}

void construct_function_to_interpolate__store_to_interped_gf(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL invdx0 = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdx1 = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdx2 = 1.0 / CCTK_DELTA_SPACE(2);
  list_of_functions_to_interpolate(cctkGH,cctk_lsh,cctk_nghostzones,invdx0,invdx1,invdx2,
                                   InterpCounter,
                                   rho_b,P,
                                   vx,vy,vz,
                                   Bx,By,Bz,
                                   gxx,gxy,gxz,gyy,gyz,gzz,
                                   betax,betay,betaz,alp, interped_gf);
// interped_gf will be interpolated across AMR boundaries, meaning that
//    it must be prolongated. Only gridfunctions with 3 timelevels stored
//    may be prolongated (provided time_interpolation_order is set to the
//    usual value of 2). We should only call this interpolation routine
//    at iterations in which all gridfunctions are on the same timelevel
//    (usually a power of 2), which will ensure that the following
//    "filling of the timelevels" is completely correct.
#pragma omp parallel for
    for(int i=0;i<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];i++) {
        interped_gf_p[i]   = interped_gf[i];
        interped_gf_p_p[i] = interped_gf[i];
    }
}
```

    Writing interp_sph_grids_ETK/src/construct_function_to_interpolate__store_to_interped_gf.cc


<a id='nrpygetgfinterporder'></a>

## Step 2.c: `get_gf_name()` and `get_interp_order()` functions \[Back to [top](#toc)\]
$$\label{nrpygetgfinterporder}$$


```python
with open(os.path.join(Ccodesdir,"src","get_gf_name.h"), "w") as file:
    file.write("void get_gf_name(const int InterpCounter,char gf_name[100]) {\n")
    for i in range(1,which_InterpCounter):
        file.write("    if(InterpCounter=="+str(i)+") { snprintf(gf_name,100,\""+gf_interp_list[i].description+"\"); return; }\n")
    file.write("    printf(\"Error. InterpCounter = %d unsupported. I should not be here.\\n\",InterpCounter); exit(1);\n")
    file.write("}\n")

with open(os.path.join(Ccodesdir,"src","get_interp_order.h"), "w") as file:
    file.write("int get_interp_order(const int InterpCounter) {\n")
    for i in range(1,which_InterpCounter):
        file.write("    if(InterpCounter=="+str(i)+") return "+gf_interp_list[i].order+";\n")
    file.write("    return -1;\n}\n")
```

<a id='nrpy_interp_counter'></a>

## Step 2.d: C Code for Initializing and incrementing `InterpCounter` \[Back to [top](#toc)\]
$$\label{nrpy_interp_counter}$$

The gridfunctions are interpolated one at a time based on the current value of the index quantity `InterpCounter`. Here we write the C code needed for initializing and incrementing this variable.


```python
with open(os.path.join(Ccodesdir,"src","define_NumInterpFunctions.h"), "w") as file:
    file.write("#define NumInterpFunctions "+str(which_InterpCounter)+"\n")
```


```python
%%writefile $Ccodesdir/src/interp_counter.cc
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "define_NumInterpFunctions.h"

void SphGrid_InitializeInterpCounterToZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  *InterpCounter = 0;

  if(verbose==2) printf("interp_sph_grids_ETK: Just set InterpCounter to %d\n",*InterpCounter);
}

void SphGrid_InitializeInterpCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration == interp_out_iteration) {
    *InterpCounter = 1;
    if(verbose==2) printf("interp_sph_grids_ETK: Just set InterpCounter to %d ; ready to start looping over interpolated gridfunctions!\n",
                          *InterpCounter);
  }
}

// This function increments InterpCounter if we are at the interp_out_iteration until
// it hits NumInterpFunctions. At this iteration, InterpCounter is set to zero, which
// exits the loop.
void SphGrid_IncrementInterpCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    if(*InterpCounter == NumInterpFunctions-1) {
        *InterpCounter = 0;
        if(verbose==2) printf("interp_sph_grids_ETK: Finished! Just zeroed InterpCounter.\n");
    } else {
        (*InterpCounter)++;
        if(verbose==2) printf("interp_sph_grids_ETK: Just incremented InterpCounter to %d of %d\n",*InterpCounter,NumInterpFunctions-1);
    }
}
```

    Writing interp_sph_grids_ETK/src/interp_counter.cc


<a id='cclfiles'></a>

# Step 3: Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{cclfiles}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn:

1. `interface.ccl`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns.
1. `param.ccl`: specifies free parameters within the thorn.
1. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions.

<a id='makecodedefn'></a>

## Step 3.a: `make.code.defn` \[Back to [top](#toc)\]
$$\label{makecodedefn}$$

Before writing the "ccl" files, we first add Einstein Toolkit's equivalent of a Makefile, the `make.code.defn` file:


```python
%%writefile $Ccodesdir/src/make.code.defn
# Main make.code.defn file for thorn interp_sph_grids_ETK

# Source files in this directory
SRCS =  main_function.cc interp_counter.cc construct_function_to_interpolate__store_to_interped_gf.cc
```

    Writing interp_sph_grids_ETK/src/make.code.defn


<a id='interfaceccl'></a>

## Step 3.b: `interface.ccl` \[Back to [top](#toc)\]
$$\label{interfaceccl}$$

Let's now write `interface.ccl`. The [official Einstein Toolkit (Cactus) documentation](http://einsteintoolkit.org/usersguide/UsersGuide.html) defines what must/should be included in an `interface.ccl` file [**here**](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 


```python
%%writefile $Ccodesdir/interface.ccl

# With "implements", we give our thorn its unique name.
implements: interp_sph_grids_ETK

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits:   admbase IllinoisGRMHD Grid

# Tell the Toolkit that we want "interped_gf" and "InterpCounter"
#    and invariants to NOT be visible to other thorns, by using
#    the keyword "private". Note that declaring these
#    gridfunctions here *does not* allocate memory for them;
#    that is done by the schedule.ccl file.
private:
CCTK_REAL interpolation_gf type=GF timelevels=3 tags='Checkpoint="no"'
{
  interped_gf
} "Gridfunction containing output from interpolation."

int InterpCounterVar type = SCALAR tags='checkpoint="no"'
{
  InterpCounter
} "Counter that keeps track of which function we are interpolating."

CCTK_REAL interp_pointcoords_and_output_arrays TYPE=ARRAY DISTRIB=CONSTANT DIM=1 SIZE=N0*N1*N2 tags='checkpoint="no"'
{
  points_x,points_y,points_z,
  output_interped
}
```

    Writing interp_sph_grids_ETK/interface.ccl


<a id='paramccl'></a>

## Step 3.c: `param.ccl` \[Back to [top](#toc)\]
$$\label{paramccl}$$

We will now write the file `param.ccl`. This file allows the listed parameters to be set at runtime. We also give allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3). 


```python
%%writefile $Ccodesdir/param.ccl

# Output the interpolated data to the IO::out_dir directory:
shares: IO
USES STRING out_dir

restricted:

########################################
# BASIC THORN STEERING PARAMETERS
CCTK_INT interp_out_iteration "Which iteration to interpolate to spherical grids?" STEERABLE=ALWAYS
{
  0:* :: ""
} 960000

## Interpolator information
CCTK_STRING interpolator_name "Which interpolator to use?" STEERABLE=ALWAYS
{
  ".+" :: "Any nonempty string; an unsupported value will throw an error."
} "Lagrange polynomial interpolation"

CCTK_INT verbose "Set verbosity level: 1=useful info; 2=moderately annoying (though useful for debugging)" STEERABLE=ALWAYS
{
  0:2 :: "0 = no output; 1=useful info; 2=moderately annoying (though useful for debugging)"
} 2
########################################
# SPHERICAL COORDINATE SYSTEM PARAMETERS
CCTK_INT N0 "Number of points in r direction" STEERABLE=ALWAYS
{
  0:* :: ""
} 96

CCTK_INT N1 "Number of points in theta direction" STEERABLE=ALWAYS
{
  0:* :: ""
} 96

CCTK_INT N2 "Number of points in phi direction" STEERABLE=ALWAYS
{
  0:* :: ""
} 96

##########
# Cartesian position of center of spherical grid (usually center of BH) -- CURRENTLY UNSUPPORTED!
CCTK_REAL x_center "x-position of center." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0

CCTK_REAL y_center "y-position of center." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0

CCTK_REAL z_center "z-position of center." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0

##########
# Radial parameters:
CCTK_REAL R0 "Radial offset: r(x0) = R_0 + exp(x0). Probably should keep it set to zero." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0

CCTK_REAL Rin "x0 offset: x0 = log(Rin-R0) + (i + 0.5)Dx0." STEERABLE=ALWAYS
{
  0:* :: ""
} 1.08986052555408

CCTK_REAL Rout "Dx0 = log( (Rout-R0) / (Rin-R0) )/N0" STEERABLE=ALWAYS
{
  0:* :: ""
} 80.0

##########
# Theta parameters:
CCTK_REAL x1_beg "x1 offset: x1 = x1_beg + (j + 0.5)Dx1. Probably should keep it set to zero." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0

CCTK_INT theta_option "Which prescription for theta should be used? 1 or 2?" STEERABLE=ALWAYS
{
  1:2 :: ""
} 1

CCTK_REAL th_c "theta_c: Angular cutout size for theta = 0 and pi" STEERABLE=ALWAYS
{
  0:* :: ""
} 0.053407075111026485 # 0.017*pi

CCTK_REAL xi "Amplitude of nonlinear part of the theta distribution." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.25

CCTK_INT th_n "Power of nonlinear part of theta distribution. Only for theta_option=2" STEERABLE=ALWAYS
{
  0:* :: ""
} 9

##########
# Phi parameters:
CCTK_REAL x2_beg "x2 offset: x2 = x2_beg + (k + 0.5)Dx2. Probably should keep it set to zero." STEERABLE=ALWAYS
{
  0:* :: ""
} 0.0
########################################
```

    Writing interp_sph_grids_ETK/param.ccl


<a id='scheduleccl'></a>

## Step 3.d: `schedule.ccl` \[Back to [top](#toc)\]
$$\label{scheduleccl}$$

Finally, we will write the file `schedule.ccl`; its official documentation is found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). 

This file declares storage for variables declared in the `interface.ccl`file and specifies when the various parts of the thorn will be run:


```python
%%writefile $Ccodesdir/schedule.ccl

STORAGE: interpolation_gf[3]
STORAGE: InterpCounterVar
STORAGE: interp_pointcoords_and_output_arrays

#############################
SCHEDULE SphGrid_InitializeInterpCounterToZero AT CCTK_INITIAL
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable to zero"

SCHEDULE SphGrid_InitializeInterpCounterToZero AT CCTK_POST_RECOVER_VARIABLES
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable to zero"

SCHEDULE SphGrid_InitializeInterpCounter before SphGrid_InterpGroup AT CCTK_ANALYSIS
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable"
##################

SCHEDULE GROUP SphGrid_InterpGroup AT CCTK_ANALYSIS BEFORE CarpetLib_printtimestats BEFORE CarpetLib_printmemstats AFTER Convert_to_HydroBase WHILE interp_sph_grids_ETK::InterpCounter
{
} "Perform all spherical interpolations. This group is only actually scheduled at cctk_iteration==interp_out_iteration."

SCHEDULE construct_function_to_interpolate__store_to_interped_gf in SphGrid_InterpGroup before DoSum
{
  STORAGE: interpolation_gf[3],InterpCounterVar,interp_pointcoords_and_output_arrays
  OPTIONS: GLOBAL,LOOP-LOCAL
  SYNC: interpolation_gf
  LANG: C
} "Construct the function to interpolate"

SCHEDULE Interpolate_to_sph_grid_main_function in SphGrid_InterpGroup after construct_function_to_interpolate__store_to_interped_gf
{
  OPTIONS: GLOBAL
  LANG: C
} "Perform interpolation and output result to file."
#######
SCHEDULE SphGrid_IncrementInterpCounter in SphGrid_InterpGroup after Interpolate_to_sph_grid_main_function
{
  LANG: C
  OPTIONS: GLOBAL
} "Increment InterpCounter variable, or set to zero once loop is complete."
##################
```

    Writing interp_sph_grids_ETK/schedule.ccl


<a id='readingoutputfile'></a>

# Step 4: Python Script for Reading the Output File \[Back to [top](#toc)\]
$$\label{readingoutputfile}$$

Here is a Python code for reading the output file generated by this thorn. It is based on a collection of Python scripts written by Bernard Kelly, available [here](https://bitbucket.org/zach_etienne/nrpy/src/master/mhd_diagnostics/). 

After generating the output file `interp_sph_grids.dat` using the Einstein Toolkit thorn above, this script will read in all the data. Processing can then be done by straightforward modification of this script. Save the script as "Interp_Sph_ReadIn.py", and run it using the command

**`python Interp_Sph_ReadIn.py interp_sph_grids.dat 58 outfile`**

Currently the last parameter "outfile" is required but not used.

```python
"""
interp_sph_grids.dat File Reader. Compatible with Python 2.7+ and 3.6+ at least.

Zachariah B. Etienne

Based on Python scripts written by Bernard Kelly:
https://bitbucket.org/zach_etienne/nrpy/src/master/mhd_diagnostics/

Find the latest version of this reader at the bottom of this Jupyter notebook:
https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids.ipynb

Usage instructions:

From the command-line, run via:
python Interp_Sph_ReadIn.py interp_sph_grids.dat [number of gridfunctions (58 or so)] [outfile]

Currently the last parameter "outfile" is required but not actually used.
"""
import numpy as np
import struct
import sys
import argparse

parser = argparse.ArgumentParser(description='Read file.')
parser.add_argument("datafile", help="main data file")
parser.add_argument("number_of_gridfunctions", help="number of gridfunctions")

parser.add_argument("outfileroot", help="root of output file names")

args = parser.parse_args()

datafile = args.datafile
outfileroot = args.outfileroot
number_of_gridfunctions = int(args.number_of_gridfunctions)

print("reading from "+str(datafile))

"""
read_char_array():
Reads a character array of size="size"
from a file (with file handle = "filehandle")
and returns the character array as a proper 
Python string.
"""
def read_char_array(filehandle,size):
    reached_end_of_string = False
    chartmp = struct.unpack(str(size)+'s', filehandle.read(size))[0]

    #https://docs.python.org/3/library/codecs.html#codecs.decode
    char_array_orig = chartmp.decode('utf-8',errors='ignore')

    char_array = ""
    for i in range(len(char_array_orig)):
        char = char_array_orig[i]
        # C strings end in '\0', which in Python-ese is '\x00'.
        #   As characters read after the end of the string will
        #   generally be gibberish, we no longer append 
        #   to the output string after '\0' is reached.
        if   sys.version_info[0]==3 and bytes(char.encode('utf-8')) == b'\x00':
            reached_end_of_string = True
        elif sys.version_info[0]==2 and char ==  '\x00':
            reached_end_of_string = True

        if reached_end_of_string == False:
            char_array += char
        else:
            pass # Continue until we've read 'size' bytes
    return char_array

"""
read_header()
Reads the header from a file.
"""
def read_header(filehandle):
    # This function makes extensive use of Python's struct.unpack
    # https://docs.python.org/3/library/struct.html
    # First store gridfunction name and interpolation order used:
    # fwrite(gf_name, 100*sizeof(char), 1, file);
    gf_name = read_char_array(filehandle,100)
    # fwrite(order, sizeof(CCTK_INT), 1, file);
    order = struct.unpack('i',filehandle.read(4))[0]

    # Then the radial grid parameters:
    # fwrite( & N0, sizeof(CCTK_INT), 1, file);
    N0 =    struct.unpack('i',filehandle.read(4))[0]
    # fwrite( & R0, sizeof(CCTK_REAL), 1, file);
    R0 =    struct.unpack('d',filehandle.read(8))[0]
    # fwrite( & Rin, sizeof(CCTK_REAL), 1, file);
    Rin =   struct.unpack('d',filehandle.read(8))[0]
    # fwrite( & Rout, sizeof(CCTK_REAL), 1, file);
    Rout =  struct.unpack('d',filehandle.read(8))[0]

    # Then the grid parameters related to the theta coordinate:
    # fwrite( & N1, sizeof(CCTK_INT), 1, file);
    N1           = struct.unpack('i', filehandle.read(4))[0]
    # fwrite( & x1_beg, sizeof(CCTK_REAL), 1, file);
    x1_beg       = struct.unpack('d', filehandle.read(8))[0]
    # fwrite( & theta_option, sizeof(CCTK_INT), 1, file);
    theta_option = struct.unpack('i', filehandle.read(4))[0]
    # fwrite( & th_c, sizeof(CCTK_REAL), 1, file);
    th_c         = struct.unpack('d', filehandle.read(8))[0]
    # fwrite( & xi, sizeof(CCTK_REAL), 1, file);
    xi           = struct.unpack('d', filehandle.read(8))[0]
    # fwrite( & th_n, sizeof(CCTK_INT), 1, file);
    th_n         = struct.unpack('i', filehandle.read(4))[0]

    # Then the grid parameters related to the phi coordinate:
    # fwrite( & N2, sizeof(CCTK_INT), 1, file);
    N2     = struct.unpack('i', filehandle.read(4))[0]
    # fwrite( & x2_beg, sizeof(CCTK_REAL), 1, file);
    x2_beg = struct.unpack('d', filehandle.read(8))[0]

    magic_number_check = 1.130814081305130e-21
    # fwrite( & magic_number, sizeof(CCTK_REAL), 1, file);
    magic_number = struct.unpack('d', filehandle.read(8))[0]
    if magic_number != magic_number_check:
        print("Error: Possible file corruption: Magic number mismatch. Found magic number = "+str(magic_number)+" . Expected "+str(magic_number_check))
        exit(1)
    # fwrite( & cctk_iteration, sizeof(CCTK_INT), 1, file);
    cctk_iteration = struct.unpack('i', filehandle.read(4))[0]
    # fwrite( & cctk_time, sizeof(CCTK_REAL), 1, file);
    cctk_time      = struct.unpack('d', filehandle.read(8))[0]

    return gf_name,order,N0,R0,Rin,Rout,N1,x1_beg,theta_option,th_c,xi,th_n,N2,x2_beg,cctk_iteration,cctk_time

# Now open the file and read all the data
with open(datafile,"rb") as f:
    # Main loop over all gridfunctions
    for i in range(number_of_gridfunctions):
        # Data are output in chunks, one gridfunction at a time, with metadata
        #    for each gridfunction stored at the top of each chunk
        # First read in the metadata:
        gf_name, order, N0, R0, Rin, Rout, N1, x1_beg, theta_option, th_c, xi, th_n, N2, x2_beg, cctk_iteration, cctk_time = read_header(f)
        print("\n\nReading gridfunction "+gf_name)
        data_chunk_size = N0*N1*N2*8 # 8 bytes per double-precision number
        # Next read in the full gridfunction data
        bytechunk = f.read(data_chunk_size)
        # Process the data using NumPy's frombuffer() function:
        #   https://docs.scipy.org/doc/numpy/reference/generated/numpy.frombuffer.html
        buffer_res = np.frombuffer(bytechunk)
        # Reshape the data into a 3D NumPy array:
        #   https://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
        this_data = buffer_res.reshape(N0,N1,N2)

        # Sanity check: Make sure the output in the "middle" of the grid looks reasonable.
        ii = int(N0/2)
        jj = int(N1/2)
        kk = int(N2/2)
        with open("output-gf"+str(i)+".txt","w") as file:
            for ii in range(N0):
                for kk in range(N2):
                    r  = ii*1.0/N0
                    th = (jj*1.0)*np.pi/N1
                    ph = (kk*1.0)*2.0*np.pi/N2
                    xx = r*np.sin(th)*np.cos(ph)
                    yy = r*np.sin(th)*np.sin(ph)
                    zz = r*np.cos(th)
                    file.write(str(xx)+" "+str(yy)+" "+str(zz)+" "+str(this_data[kk,jj,ii])+"\n")

```

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids.pdf](Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids")
```

    Created Tutorial-ETK_thorn-Interpolation_to_Spherical_Grids.tex, and
        compiled LaTeX file to PDF file Tutorial-ETK_thorn-
        Interpolation_to_Spherical_Grids.pdf

