<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# `interp_arbgrid_MO_ETK`: An Einstein Toolkit module for interpolation to arbitrary grids, at multiple interpolation orders, in a Cartesian basis. 

## (Includes notes on transformations to other coordinate bases.)

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This module is designed to interpolate arbitrary quantities on [Einstein Toolkit](https://einsteintoolkit.org/) Adaptive-Mesh Refinement (AMR) grids (using the [Carpet](https://carpetcode.org/) AMR infrastructure) to numerical grids with arbitrary sampling. The tutorial elaborates on core ETK C routines for interpolation, NRPy+ output for required gridfunctions, and the necessary CCL files for interfacing with the Toolkit.

**Validation Status:** <font color='red'><b> In progress</b></font>

**Validation Notes:** This module is currently undergoing validation testing.


## Introduction: 

Given some set of $N$ quantities $\mathbf{Q}=\{Q_0,Q_1,Q_2,...,Q_{N-2},Q_{N-1}\}$, this module performs the following for each $Q_i$:

1. Evaluate $Q_i$ at all gridpoints that are not ghost zones. Sometimes $Q_i$ is computed using finite difference derivatives, so this is necessary.
1. Call upon Carpet's interpolation and interprocessor synchronization functions to fill in $Q_i$ at all ghost zones, *except* at the outer boundary. We do not generally trust $Q_i$ at the outer boundary due to errors associated with the approximate outer boundary conditions. 
1. At this point, $Q_i$ is set at all gridpoints except ghost zones at the outer boundary. Interpolate $Q_i$ to the desired output grids, **maintaining the Cartesian basis for all vectors and tensors**, and append the result to a file.

This tutorial notebook takes a three-part structure. First, all the needed core Einstein Toolkit (ETK) C routines for interpolation are presented. Second, NRPy+ is used to output gridfunctions needed on the output grids. Third, the needed files for interfacing this module with the rest of the Einstein Toolkit (ccl files) are specified.

<a id='toc'></a>

# Table of Contents: 
$$\label{toc}$$

1. [Step 1](#etkmodule): Setting up the Core C Code for the Einstein Toolkit Module
    1. [Step 1.a](#etk_interp): Low-Level Einstein Toolkit Interpolation Function
    1. [Step 1.b](#fileformat): Outputting to File
    1. [Step 1.c](#maininterpolator): The Main Interpolator Driver Function
    1. [Step 1.d](#standalonerandompoints): Standalone C code to output random points data 
1. [Step 2](#nrpy): Using NRPy+ to Generate C Code for Needed Gridfunctions
    1. [Step 2.a](#nrpy_list_of_funcs_interp): Set up NRPy-based `list_of_functions_to_interpolate.h`
        1. [Step 2.a.i](#nrpygrmhd): All GRMHD quantities, except vector potential $A_i$
        1. [Step 2.a.ii](#unstaggera): Unstagger $A_i$ and add to "list of functions to interpolate"
        1. [Step 2.a.iii](#nrpy4metric): Compute all 10 components of the 4-metric $g_{\mu\nu}$
        1. [Step 2.a.iv](#nrpy4christoffels_cartesian):Compute all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$
        1. [Step 2.a.v](#nrpy4christoffels_spherical):  Notes on computing all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$ in the Spherical basis
        1. [Step 2.a.vi](#nrpybasisxform): Notes on basis transforming all Cartesian basis quantities to spherical
        1. [Step 2.a.vii](#psi4andfriends): Output Weyl scalars $\psi_0$ through $\psi_4$, as well as Weyl invariants $J$ and $I$, from the `WeylScal4` ETK thorn
        1. [Step 2.a.viii](#constraints_gfs): Hamiltonian and momentum constraints
        1. [Step 2.a.ix](#nuc_eos_gfs): Nuclear (tabulated) equation of state quantities
    1. [Step 2.b](#nrpy_c_calling_function): C code calling function for the NRPy+ C output
    1. [Step 2.c](#nrpygetgfname): The `get_gf_name()` function
    1. [Step 2.d](#nrpy_interp_counter): C Code for Initializing and incrementing `InterpCounter`
    1. [Step 2.e](#validationagainstfm): Validation of interpolated data against exact Fishbone-Moncrief data
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
Ccodesdir = "interp_arbgrid_MO_ETK"
# First remove C code output directory and all subdirectories if they exist
# Courtesy https://stackoverflow.com/questions/303200/how-do-i-remove-delete-a-folder-that-is-not-empty
shutil.rmtree(Ccodesdir, ignore_errors=True)
# Then create a fresh directory
cmd.mkdir(Ccodesdir)
cmd.mkdir(os.path.join(Ccodesdir,"src/"))
cmd.mkdir(os.path.join(Ccodesdir,"src","standalone/"))
```

<a id='etk_interp'></a>

## Step 1.a: Low-Level ETK Interpolation Function \[Back to [top](#toc)\]
$$\label{etk_interp}$$

We start by writing the low-level interpolation function **`Interpolate_to_dest_grid()`**, which  to file. 

**`Interpolate_to_dest_grid()`** takes as input

* **cctkGH**: Information about the underlying Cactus/Carpet grid hierarchy.
* **interp_num_points**: Number of destination interpolation points
* **point_x_temp, point_y_temp, point_z_temp**: Cartesian $(x,y,z)$ location for each of the **interp_num_points** interpolation points.
* **input_array_names[1]**: List of input gridfunction names to interpolate. We will do this only one gridfunction at a time, for gridfunction $Q_i$, as described above.

**`Interpolate_to_dest_grid()`** outputs:

* **output_f[1]**: The gridfunction **input_array_names[1]** interpolated to the set of **interp_num_points** specified in the input.


```python
%%writefile $Ccodesdir/src/Interpolate_to_dest_grid.h

void Interpolate_to_dest_grid(const cGH *restrict cctkGH,
                              const CCTK_INT interp_num_points,
                              const CCTK_INT interp_order,
                              char interpolator_name[100],
                              const CCTK_REAL *restrict point_x_temp,
                              const CCTK_REAL *restrict point_y_temp,
                              const CCTK_REAL *restrict point_z_temp,
                              const CCTK_STRING input_array_names[1],
                              CCTK_REAL *output_f[1]) {
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

  const void* interp_coords[3] = { (const void *) point_x_temp,
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

  void* output_arrays[NUM_OUTPUT_ARRAYS] = { (void *) output_f[0] };

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

    Writing interp_arbgrid_MO_ETK/src/Interpolate_to_dest_grid.h



```python
%%writefile $Ccodesdir/src/interpolate_set_of_points_in_file.h

#define ALLOCATE_2D_GENERIC(type,array,ni,nj) type **array=(type **)malloc(ni * sizeof(type *)); \
  for(int cc = 0; cc < ni; cc++)                 array[cc]=(type * )malloc(nj * sizeof(type));
#define FREE_2D_GENERIC(type,array,ni,nj) for(int cc = 0; cc < ni;cc++) free((void *)array[cc]); \
  /**/                                                                  free((void *)array);
#include "output_to_file.h"

// Calls the above function and output_to_file().
void interpolate_set_of_points_in_file(CCTK_ARGUMENTS,char filename_basename[100],char gf_name[100],char interpolator_name[100],int num_interp_orders,int *interp_orders_list) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS; // Needed for x_center,y_center,z_center
  // Set up output array:
  // The name of the input gridfunction is always "interp_arbgrid_MO_ETK::interped_gf":
  const CCTK_STRING input_array_names[1] = { "interp_arbgrid_MO_ETK::interped_gf" };
  CCTK_REAL *points_dest_grid_x,*points_dest_grid_y,*points_dest_grid_z; // Coordinates of points of destination grid
  CCTK_REAL **output_f; // Output to be written to dataset, will be filled with NaNs at points out of bounds
  // For benchmarking purposes:
  time_t start_timer,end_timer;
  time(&start_timer); // Resolution of one second...
  CCTK_REAL time_in_seconds;

  int num_dest_grid_points;
  if(CCTK_MyProc(cctkGH)==0) {
    // Step 1: Read list of desired interpolation destination points from file:

    // Step 1.a: Read integer at top of file indicating number of points.
    int num_dest_grid_pointsx,num_dest_grid_pointsy,num_dest_grid_pointsz;
    char pointsx_filename[100]; snprintf(pointsx_filename,100,"%s-x.dat",filename_basename);
    printf("Reading list of x data points from file %s...\n",pointsx_filename);
    FILE *pointsx_file = fopen(pointsx_filename, "rb");
    if(!pointsx_file) { printf("Error: Unable to open %s\n",pointsx_filename); exit(1); }
    fread(&num_dest_grid_pointsx, sizeof(int), 1, pointsx_file);

    char pointsy_filename[100]; snprintf(pointsy_filename,100,"%s-y.dat",filename_basename);
    printf("Reading list of y data points from file %s...\n",pointsy_filename);
    FILE *pointsy_file = fopen(pointsy_filename, "rb");
    if(!pointsy_file) { printf("Error: Unable to open %s\n",pointsy_filename); exit(1); }
    fread(&num_dest_grid_pointsy, sizeof(int), 1, pointsy_file);

    char pointsz_filename[100]; snprintf(pointsz_filename,100,"%s-z.dat",filename_basename);
    printf("Reading list of z data points from file %s...\n",pointsz_filename);
    FILE *pointsz_file = fopen(pointsz_filename, "rb");
    if(!pointsz_file) { printf("Error: Unable to open %s\n",pointsz_filename); exit(1); }
    fread(&num_dest_grid_pointsz, sizeof(int), 1, pointsz_file);

    // Step 1.a.i: Sanity check: make sure that num_dest_grid_pointsx == num_dest_grid_pointsy == num_dest_grid_pointsz
    if(num_dest_grid_pointsx != num_dest_grid_pointsy || num_dest_grid_pointsy != num_dest_grid_pointsz) {
      printf("Error: Failed sanity check. Number of interpolation points different in %s-{x,y,z}.dat data files!\n",
             filename_basename);
      exit(1);
    } else {
      // If sanity check passes:
      num_dest_grid_points = num_dest_grid_pointsx;
    } // END sanity check

      // Step 1.b: Allocate memory for destination grids and interpolation output
    if(num_dest_grid_points <= 0 || num_dest_grid_points > 2000000000) {
      printf("Error: Failed sanity check. Number of interpolation points was found to be: %d",num_dest_grid_points);
      exit(1);
    } // END sanity check
    points_dest_grid_x = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
    points_dest_grid_y = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
    points_dest_grid_z = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
    output_f = (CCTK_REAL **)malloc(1 * sizeof(CCTK_REAL *));
    for(int cc = 0; cc < 1; cc++) output_f[cc]=(CCTK_REAL *)malloc(num_dest_grid_points * sizeof(CCTK_REAL));

    // Step 1.c: Store cell-centered points to allocated memory.
    fread(points_dest_grid_x, sizeof(CCTK_REAL), num_dest_grid_points, pointsx_file);
    fread(points_dest_grid_y, sizeof(CCTK_REAL), num_dest_grid_points, pointsy_file);
    fread(points_dest_grid_z, sizeof(CCTK_REAL), num_dest_grid_points, pointsz_file);
    int magic_numberx; fread(&magic_numberx, sizeof(int), 1, pointsx_file);
    int magic_numbery; fread(&magic_numbery, sizeof(int), 1, pointsy_file);
    int magic_numberz; fread(&magic_numberz, sizeof(int), 1, pointsz_file);
    int correct_magicnum = -349289480;
    if(magic_numberx != correct_magicnum || magic_numbery != correct_magicnum || magic_numberz != correct_magicnum) {
      printf("Error: Failed sanity check. Magic numbers in x,y,z data files were: %d %d %d, respectively, but should have been: %d",
             magic_numberx,magic_numbery,magic_numberz,correct_magicnum);
      exit(1);
    }
    fclose(pointsx_file);
    fclose(pointsy_file);
    fclose(pointsz_file);
    time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);
    printf("Finished in %e seconds.\n",time_in_seconds);

    // Step 1.d: Apply offset to x,y,z coordinates to ensure they are centered on (x_center,y_center,z_center)
    //           For example, if a black hole is situated at (x,y,z) = (1,2,3), then we set
    //           (x_center,y_center,z_center) = (1,2,3) in our ETK parameter file (i.e., with extension .par)
    //           and if we desire a point at (x_dest,y_dest,z_dest) = (0,0,0) on the *destination* grid,
    //           this will correspond to point (x_src,y_src,z_src) = (1,2,3) = (x_center,y_center,z_center)
    //           on the source grid. Thus the translation between source and destination grids is given by
    //           (x_src,y_src,z_src) = (x_dest+x_center, y_dest+y_center, z_dest+z_center),
    //           where (x_src,y_src,z_src) = (points_dest_grid_x[i],points_dest_grid_y[i],points_dest_grid_z[i]) for point i.
    for(int point=0;point<num_dest_grid_points;point++) {
      points_dest_grid_x[point] += x_center;
      points_dest_grid_y[point] += y_center;
      points_dest_grid_z[point] += z_center;
    }

    // Step 1.e: Look for the points that are out of bounds and set their coordinates to out_of_bounds_interp_xyz.
    //           At the end, we will replace the interpolation output at these points with NaNs.
    printf("Looking for points that are out of bounds.\n");
    int num_out_of_bounds_points = 0;
    for(int i=0;i<num_dest_grid_points;i++){
      if(fabs(points_dest_grid_x[i])>out_of_bounds_interp_xyz ||
         fabs(points_dest_grid_y[i])>out_of_bounds_interp_xyz ||
         fabs(points_dest_grid_z[i])>out_of_bounds_interp_xyz) {
        points_dest_grid_x[i] = out_of_bounds_interp_xyz;
        points_dest_grid_y[i] = out_of_bounds_interp_xyz;
        points_dest_grid_z[i] = out_of_bounds_interp_xyz;
        num_out_of_bounds_points++;
      }
    }
    time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);
    printf("Finished in %e seconds, found %i points out of bounds.\n", time_in_seconds, num_out_of_bounds_points);

  } // END if(CCTK_MyProc(cctkGH)==0)


  // Step 1.f: Looping over interp order as desired, interpolate to destination points & output to file
  for(int order_i=0; order_i<num_interp_orders; order_i++) {
    int order = interp_orders_list[order_i];
    printf("Interpolating\033[1m %s \033[0m... using interpolation order = %d\n",gf_name,order);
    printf("and %s \n",interpolator_name);
    if(CCTK_MyProc(cctkGH)==0) {
      Interpolate_to_dest_grid(cctkGH, num_dest_grid_points, order, interpolator_name,
                               points_dest_grid_x,points_dest_grid_y,points_dest_grid_z, input_array_names, output_f);

      // Step 1.f.i: Sanity check -- check for bad point:
#pragma omp parallel for
      for(int i=0;i<num_dest_grid_points;i++) {
        if(output_f[0][i] > 1e20) {
          printf("BAD POINT: %s %d %e %e %e %e\n",gf_name,i,points_dest_grid_x[i],points_dest_grid_y[i],points_dest_grid_z[i], output_f[0][i]);
          exit(1);
        }
      }
      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);
      printf("Finished in %d seconds. Next: Filling cells out of bounds with NaNs\n",time_in_seconds);
      // Step 1.f.ii: Filling cells out of bounds with NaNs:
      for(int i=0;i<num_dest_grid_points;i++){
        if(fabs(points_dest_grid_x[i])==out_of_bounds_interp_xyz ||
           fabs(points_dest_grid_y[i])==out_of_bounds_interp_xyz ||
           fabs(points_dest_grid_z[i])==out_of_bounds_interp_xyz) {
          output_f[0][i] = 0.0 / 0.0; // ie NaN
        }
      }
      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);
      printf("Finished in %e seconds. Next: Interpolate_to_dest_grid_main_function(): Outputting to file at iteration %d\n",time_in_seconds,cctk_iteration);
      output_to_file(CCTK_PASS_CTOC,gf_name,&order,&num_dest_grid_points,output_f);
      time(&end_timer); time_in_seconds = difftime(end_timer,start_timer); time(&start_timer);
      printf("Finished in %e seconds. Interpolate_to_dest_grid_main_function(): Finished output to file at iteration %d\n",time_in_seconds,cctk_iteration);
    } else {
      // On all MPI processes that are nonzero, only call the interpolation function
      //    to ensure the MPI calls from the actual interpolation (driven by proc==0) are seen.
      // Setting num_dest_grid_points to zero results in a segfault on certain (ahem, Frontera)
      //    systems. So we set num_dest_grid_points = 1 and interpolate to the origin
      //    only for MPI processes that are nonzero, leaving the heavy lifting to MPI process 0.
      num_dest_grid_points = 1;
      points_dest_grid_x = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
      points_dest_grid_y = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
      points_dest_grid_z = (CCTK_REAL  *)malloc(sizeof(CCTK_REAL)*num_dest_grid_points);
      output_f = (CCTK_REAL **)malloc(1 * sizeof(CCTK_REAL *));
      for(int cc = 0; cc < 1; cc++) output_f[cc]=(CCTK_REAL *)malloc(num_dest_grid_points * sizeof(CCTK_REAL));

      Interpolate_to_dest_grid(cctkGH, num_dest_grid_points, order, interpolator_name,
                               points_dest_grid_x,points_dest_grid_y,points_dest_grid_z, input_array_names, output_f);
      output_f[0][0] = 0.0;
    } // END if(CCTK_MyProc(cctkGH)==0)
  } // END for(int order_i=0; order_i<num_interp_orders; order_i++)

  // Step 1.g: Free memory for destination grids and interpolation output
  free(points_dest_grid_x);
  free(points_dest_grid_y);
  free(points_dest_grid_z);
  FREE_2D_GENERIC(CCTK_REAL,output_f,1,num_dest_grid_points);
} // END function
#undef ALLOCATE_2D_GENERIC
#undef FREE_2D_GENERIC
```

    Writing interp_arbgrid_MO_ETK/src/interpolate_set_of_points_in_file.h


<a id='fileformat'></a>

## Step 1.b: Outputting to File (File format notes) \[Back to [top](#toc)\]
$$\label{fileformat}$$

Since they take almost no space relative to the data chunks, we attach the entire metadata to each interpolated function that is output:


```python
%%writefile $Ccodesdir/src/output_to_file.h

#include "define_NumInterpFunctions.h"

// output_to_file() starts order and InterpCounter both with the value 1
void output_to_file(CCTK_ARGUMENTS,
                    char gf_name[100],
                    CCTK_INT *restrict order,
                    CCTK_INT *restrict num_interp_points,
                    CCTK_REAL *restrict output_f[1]) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char filename[100];
  sprintf (filename, "%s/interp_dest_grids_MO.dat", out_dir);
  FILE *file;
  if(*InterpCounter == 1 && *order==1) {
    file = fopen (filename,"w");
    printf("WRITING to file %s\n",filename);
    // Write EOS information to the beginning of the file
    char eos_info[256];
    if( enable_nuc_eos ) {
      // If using nuclear EOS, give the table name
      sprintf(eos_info,"Nuclear (tabulated) EOS from file %s",nuceos_table_name);
    }
    else {
      // If using hybrid EOS, then give Gamma_th and the number of polytropic pieces
      sprintf(eos_info,"Hybrid EOS with Gamma_th = %.15e and %d polytropic pieces",Gamma_th,neos);
    }
    fwrite(eos_info, 256*sizeof(char), 1, file);
  }
  else {
    file = fopen (filename,"a+");
    printf("Appending to file %s\n",filename);
  }
  if (! file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "interp_dest_grid__ET_thorn: Cannot open output file '%s'", filename);
    exit(1);
  }

  fwrite(gf_name, 100*sizeof(char), 1, file);
  fwrite(order, sizeof(CCTK_INT), 1, file);
  fwrite(num_interp_points, sizeof(CCTK_INT),1,file);

  CCTK_REAL magic_number = 1.130814081305130e-21;
  fwrite(&magic_number, sizeof(CCTK_REAL), 1, file);
  fwrite(&cctk_iteration, sizeof(CCTK_INT), 1, file);
  fwrite(&cctk_time, sizeof(CCTK_REAL), 1, file);
  for(CCTK_INT i=0;i<1;i++) {
    fwrite(output_f[i], sizeof(CCTK_REAL)*(*num_interp_points), 1, file);
  }

  fclose(file);
}

```

    Writing interp_arbgrid_MO_ETK/src/output_to_file.h


<a id='maininterpolator'></a>

## Step 1.c: The Main Interpolation Driver Function \[Back to [top](#toc)\]
$$\label{maininterpolator}$$

The **`Interpolate_to_dest_grid_main_function()`** function calls the above functions as follows:

1. **`Interpolate_to_dest_grid()`** ([Above section](#etk_interp)): Interpolates to destination grid and calls
    1. **`output_to_file()`** ([Above section](#fileformat)): Outputs information about interpolation, as well as interpolation result, to file


```python
%%writefile $Ccodesdir/src/main_function.cc


// Include needed ETK & C library header files:
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> // for benchmarking
// Needed for dealing with Cactus/ETK infrastructure
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
// Needed for low-level interpolation functions
#include "util_Table.h"
#include "util_String.h"

// Include locally-defined C++ functions:
#include "Interpolate_to_dest_grid.h"
#include "get_gf_name.h"
#include "interpolate_set_of_points_in_file.h"

void Interpolate_to_dest_grid_main_function(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Perform interpolation only at iteration == interp_out_iteration:
  if(cctk_iteration != interp_out_iteration) return;
  // Perform interpolation!
  // Process zero (CCTK_MyProc(cctkGH)==0) is responsible for directing the interpolation.
  //    All other processes must see the cctk_InterpGridArrays() within Interpolate_to_dest_grid(),
  //    so that the MPI calls work properly, but these nonzero processes can call
  //    Interpolate_to_dest_grid() with number of interpolated points set to zero, and
  //    without needing a malloc().
  char gf_name[100]; get_gf_name(*InterpCounter,gf_name);
  char filename_basename[100];
  char interpolator_name[100];

  // The following "if" statement is destination-code dependent.
  // In our case, since we are interpolating variables for harm3d,
  // we interpolate the vector potential to the corners of each cell.
  // Every other quantity is interpolated to the center of each cell.
 if(strncmp(gf_name,"Unstaggered",11) == 0){
    sprintf(filename_basename,"corner_points");
  }
  else{
    sprintf(filename_basename,"cell_centered_points");
  }

  int num_interp_orders,*interp_orders_list;
  // 4-metric, 4-Christoffels and A_mu only output Hermite interpolation order==3
  // so we secure continuity of first derivative.
  if( (strncmp(gf_name,"4-",2) == 0) ||
      (strncmp(gf_name,"Unstaggered",11) == 0) ||
     (strncmp(gf_name,"Constraints",11) == 0)) {
    num_interp_orders = 1;
    sprintf(interpolator_name, "Hermite polynomial interpolation");
    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);
    interp_orders_list[0] = 3;
  } else {
    num_interp_orders = 3;
    sprintf(interpolator_name, "Lagrange polynomial interpolation");
    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);
    int count = 0; for(int order=1;order<=4;order*=2) { interp_orders_list[count] = order; count++; }
  }

  interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name,interpolator_name,num_interp_orders,interp_orders_list);
  free(interp_orders_list);

  // Now perform interpolation of 4-metric on
  //   faces (i-1/2,j,k), (i,j-1/2,k), (i,j,k-1/2) and corners (i-1/2,j-1/2,k-1/2)
  if(strncmp(gf_name,"4-metric",8) == 0) {
    num_interp_orders = 1;
    sprintf(interpolator_name, "Hermite polynomial interpolation");
    interp_orders_list = (int *)malloc(sizeof(int)*num_interp_orders);
    interp_orders_list[0] = 3;

    char gf_name_new[100];

    sprintf(filename_basename,"faceim_points");
    snprintf(gf_name_new,100,"faceim (i-1/2,j,k): %s",gf_name);
    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);

    sprintf(filename_basename,"facejm_points");
    snprintf(gf_name_new,100,"facejm (i,j-1/2,k): %s",gf_name);
    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);

    sprintf(filename_basename,"facekm_points");
    snprintf(gf_name_new,100,"facekm (i,j,k-1/2): %s",gf_name);
    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);

    sprintf(filename_basename,"corner_points");
    snprintf(gf_name_new,100,"cornr (i-1/2,j-1/2,k-1/2): %s",gf_name);
    interpolate_set_of_points_in_file(CCTK_PASS_CTOC,filename_basename,gf_name_new,interpolator_name,num_interp_orders,interp_orders_list);
  } // END if(strncmp(gf_name,"4-metric",8) == 0)
} // END function

```

    Writing interp_arbgrid_MO_ETK/src/main_function.cc


<a id='standalonerandompoints'></a>

## Step 1.d: Standalone C code to output random points data \[Back to [top](#toc)\]
$$\label{standalonerandompoints}$$


```python
%%writefile $Ccodesdir/src/standalone/standalone_C_code_genpoints.c

// Part P1: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

const double xyzmin = -1000.0;
const double xyzmax =  1000.0;

void write_to_xyz_files(int num_interp_points, char filename_basename[100]) {
  char filenamex[100],filenamey[100],filenamez[100];
  snprintf(filenamex,100,"%s-x.dat",filename_basename);
  snprintf(filenamey,100,"%s-y.dat",filename_basename);
  snprintf(filenamez,100,"%s-z.dat",filename_basename);
  FILE *filex = fopen(filenamex,"wb");
  FILE *filey = fopen(filenamey,"wb");
  FILE *filez = fopen(filenamez,"wb");

  // Write file headers:
  fwrite(&num_interp_points, sizeof(int), 1, filex);
  fwrite(&num_interp_points, sizeof(int), 1, filey);
  fwrite(&num_interp_points, sizeof(int), 1, filez);

  // Write guts of file:
  for(int ii=0;ii<num_interp_points;ii++) {
    double rngx = xyzmin + (xyzmax - xyzmin)*drand48(); // drand48() returns between 0.0 & 1.0
    double rngy = xyzmin + (xyzmax - xyzmin)*drand48();
    double rngz = xyzmin + (xyzmax - xyzmin)*drand48();
    fwrite(&rngx, sizeof(double), 1, filex);
    fwrite(&rngy, sizeof(double), 1, filey);
    fwrite(&rngz, sizeof(double), 1, filez);
  }

  // Write magic number as file footers:
  int magic_number = -349289480;
  fwrite(&magic_number, sizeof(int), 1, filex);
  fwrite(&magic_number, sizeof(int), 1, filey);
  fwrite(&magic_number, sizeof(int), 1, filez);

  // Close files.
  fclose(filex);
  fclose(filey);
  fclose(filez);
}

int main(int argc, const char *argv[]) {

  // Step 0a: Read command-line input, error out if nonconformant
  if(argc != 2 || atoi(argv[1]) < 1) {
    printf("Error: Expected one command-line argument: ./standalone_C_code_genpoints [num_interp_points],\n");
    exit(1);
  }

  const int num_interp_points = atoi(argv[1]);

  char filename_basename[100];
  sprintf(filename_basename,"cell_centered_points");
  write_to_xyz_files(num_interp_points, filename_basename);
  sprintf(filename_basename,"faceim_points");
  write_to_xyz_files(num_interp_points, filename_basename);
  sprintf(filename_basename,"facejm_points");
  write_to_xyz_files(num_interp_points, filename_basename);
  sprintf(filename_basename,"facekm_points");
  write_to_xyz_files(num_interp_points, filename_basename);
  sprintf(filename_basename,"corner_points");
  write_to_xyz_files(num_interp_points, filename_basename);

  return 0;
}
```

    Writing interp_arbgrid_MO_ETK/src/standalone/standalone_C_code_genpoints.c


<a id='nrpy'></a>

# Step 2: Use NRPy+ C Output to Set All Output Gridfunctions \[Back to [top](#toc)\]
$$\label{nrpy}$$



```python
# Step 1: Import needed NRPy+ parameters
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
from outputC import lhrh         # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import loop as lp                # NRPy+: Generate C code loops

par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

from collections import namedtuple
gf_interp = namedtuple('gf_interp', 'gf_description')
gf_interp_list = []
gf_interp_list.append(gf_interp("dummy -- used because this is a 1-offset array"))

interped_gf = gri.register_gridfunctions("AUX","interped_gf")

def interp_fileout(which_InterpCounter, expression, filename):
    kernel = fin.FD_outputC("returnstring",lhrh(lhs=gri.gfaccess("out_gfs","interped_gf"),rhs=expression),"outCverbose=False")
    output_type="a"
    if which_InterpCounter == 1:
        # Write the file header, which includes #define's for GAMMA_SPEED_LIMIT and TINYDOUBLE:
        with open(filename, "w") as file:
            file.write("// Parameters needed to ensure velocity computations are robust:\n")
            file.write("#define GAMMA_SPEED_LIMIT 20\n")
            file.write("#define TINYDOUBLE        1e-100\n\n")
    compute_xx0xx1xx2 = ""
    if "SPHERICAL" in gf_interp_list[which_InterpCounter].gf_description:
        compute_xx0xx1xx2 = """
// ONLY NEEDED/USED IF CONVERTING TO SPHERICAL BASIS:
const double Cartx = x[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - x_center;
const double Carty = y[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - y_center;
const double Cartz = z[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] - z_center;

const double xx0 = sqrt(Cartx*Cartx + Carty*Carty + Cartz*Cartz);
const double xx1 = acos(Cartz/xx0);
const double xx2 = atan2(Carty,Cartx);\n
"""
    with open(filename, output_type) as file:
        if( "temperature" in str(interp_expr) or
            "Y_e" in str(interp_expr) or
            "eps" in str(interp_expr) or
            "entropy"in str(interp_expr) ):
            file.write("if(*InterpCounter == "+str(which_InterpCounter)+" && enable_nuc_eos) {\n")
        else:
            file.write("if(*InterpCounter == "+str(which_InterpCounter)+") {\n")
        file.write("    // Interpolating: "+gf_interp_list[which_InterpCounter].gf_description+"\n")
        file.write(lp.loop(["i2","i1","i0"],
                           ["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],\
                           ["cctk_lsh[2]-cctk_nghostzones[2]",
                            "cctk_lsh[1]-cctk_nghostzones[1]",
                            "cctk_lsh[0]-cctk_nghostzones[0]"],\
                           ["1","1","1"],\
                           ["#pragma omp parallel for","",""],"   ",
                           compute_xx0xx1xx2+kernel))
        file.write("}\n")
    # If successful, return incremented which_InterpCounter:
    return which_InterpCounter+1
```

<a id='nrpy_list_of_funcs_interp'></a>

## Step 2.a: Set up NRPy-based `list_of_functions_to_interpolate.h` \[Back to [top](#toc)\]
$$\label{nrpy_list_of_funcs_interp}$$

First specify NRPy+ output file and initialize `which_InterpCounter`, which keeps track of the number of interpolated functions on the grid


```python
NRPyoutfilename = os.path.join(Ccodesdir,"src","list_of_functions_to_interpolate.h")

which_InterpCounter = 1
```

<a id='nrpygrmhd'></a>

### Step 2.a.i: GRMHD quantities \[Back to [top](#toc)\]
$$\label{nrpygrmhd}$$

These include
* $\rho_b$, the baryonic density (e.g. $\text{HydroBase::rho}$ or $\text{IllinoisGRMHD::rho_b}$)
* $P$, the total gas pressure (e.g. $\text{HydroBase::press}$ or $\text{IllinoisGRMHD::P}$)
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
# INPUT GRIDFUNCTIONS: The AUX or EVOL designation is *not* used in diagnostic modules.

gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01")
alpha = gri.register_gridfunctions("AUX","alpha")
betaU = ixp.register_gridfunctions_for_single_rank1("AUX","betaU")

# Add a constant beta offset, to account for linear
#      (i.e., constant velocity) coordinate drift.
# Note that beta_offsetU's are set in param.ccl.
# As beta_offsetU is constant in space, it has no
#      impact on betaU_dD's.
beta_offsetU0,beta_offsetU1,beta_offsetU2 = par.Cparameters("REAL","modulenamedoesntmatter",
                                                            ["beta_offsetU0","beta_offsetU1","beta_offsetU2"],
                                                            [0.0,0.0,0.0])
betaU[0] += beta_offsetU0
betaU[1] += beta_offsetU1
betaU[2] += beta_offsetU2

# Tensors are given in Cartesian basis:
# Derivatives of metric
gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
alpha_dD   = ixp.declarerank1("alpha_dD")

DIM=3

IGMvU = ixp.register_gridfunctions_for_single_rank1("AUX","IGMvU")
BU    = ixp.register_gridfunctions_for_single_rank1("AUX","BU")

gf_interp_list.append(gf_interp("IGM density primitive"))
rho_b       = gri.register_gridfunctions("AUX","rho_b")
interp_expr = rho_b
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("IGM pressure primitive"))
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

We use expressions for $v_{(n)}^i$, $\Gamma v_{(n)}^i$, and $v^i$ as implemented [in the GRHD equations notebook](Tutorial-GRHD_Equations-Cartesian.ipynb#convertvtou) and corresponding [GRHD.equations Python module](../edit/GRHD/equations.py). These expressions enforce a speed limit on $\Gamma$ to ensure e.g., the denominator within the radical in the above expression for $\Gamma v_{(n)}^i$ is never negative, which would result in `NaN`s. 


```python
import GRHD.equations as Ge
Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, IGMvU)

# ValenciavU = ixp.zerorank1()
# for i in range(DIM):
#     ValenciavU[i] = 1/alpha * (IGMvU[i] + betaU[i])
# v_dot_v = sp.sympify(0)
# for i in range(DIM):
#     for j in range(DIM):
#         v_dot_v += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]
# u4Uzero = sp.sqrt(1/(1 - v_dot_v))/alpha # u^0 = LorentzGamma/alpha

u4Uzero = Ge.u4U_ito_vU[0]
gf_interp_list.append(gf_interp("u^0: zero (time) component of 4-velocity"))
interp_expr = u4Uzero
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

Gamma_times_ValenciavU = ixp.zerorank1()
Gamma      = alpha*Ge.u4U_ito_vU[0]
ValenciavU = ixp.zerorank1()
for i in range(DIM):
    ValenciavU[i] = (1/alpha * (Ge.rescaledvU[i] + betaU[i])) # simplify?

for i in range(DIM):
    Gamma_times_ValenciavU[i] = Gamma*ValenciavU[i]
```

Next we'll basis transform $W^i = \Gamma v_{\rm n}^i$ from Cartesian to the spherical basis:

Within [`reference_metric.py`](../edit/reference_metric.py), the `compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()` function defines Jacobians relative to the center of the source (reference metric) grid, at a point $x^j_{\rm src}=$(`xx0,xx1,xx2`)${}_{\rm src}$ on the source grid:
$$
{\rm Jac\_dUCart\_dDsrcUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm src}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUsrc\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm src}}{\partial x^j_{\rm Cart}},
$$

using NRPy+'s `generic_matrix_inverter3x3()` function. 

In terms of these, the transformation of $W^i$ from Cartesian coordinates to `"reference_metric::CoordSystem=Spherical"`may be written:

\begin{align}
W^i_{\rm Sph} &= \frac{\partial x^i_{\rm Sph}}{\partial x^\ell_{\rm Cart}} W^\ell_{\rm Cart}
\end{align}


```python
import reference_metric as rfm    # NRPy+: Reference metric support

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

# Step 2.a: Construct Jacobian & Inverse Jacobians:
Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

Sph_basis_Gamma_times_ValenciavU = ixp.zerorank1()
for i in range(DIM):
    for l in range(DIM):
        Sph_basis_Gamma_times_ValenciavU[i] += Jac_dUrfm_dDCartUD[i][l] * Gamma_times_ValenciavU[l]

for i in range(DIM):
    gf_interp_list.append(gf_interp("*SPHERICAL BASIS* Lorentz factor, times Valencia vU"+str(i)))
    interp_expr = Sph_basis_Gamma_times_ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```


```python
for i in range(DIM):
    gf_interp_list.append(gf_interp("(speed-limited) Valencia 3-velocity vU"+str(i)))
    interp_expr = ValenciavU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

for i in range(DIM):
    if i==0:
        gf_interp_list.append(gf_interp("Local grid resolution dx=dy=dz"))
        invdx0 = sp.symbols('invdx0', real=True)
        interp_expr = 1/invdx0
    else:
        gf_interp_list.append(gf_interp("(speed-limited) IGM 3-velocity vU"+str(i)+" = u^i divided by u^0"))
        interp_expr = Ge.u4U_ito_vU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# For testing:
# gf_interp_list.append(gf_interp("Lorentz factor"))
# interp_expr = v_dot_v
# which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("Valencia vU"+str(i)))
#     interp_expr = Valenciav[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

for i in range(DIM):
    gf_interp_list.append(gf_interp("IGM magnetic field component B"+str(i)))
    interp_expr = BU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='unstaggera'></a>

### Step 2.a.ii: Unstagger $A_i$ and add to "list of functions to interpolate" \[Back to [top](#toc)\]
$$\label{unstaggera}$$

First generate the C code needed to unstagger the A-fields.


```python
%%writefile $Ccodesdir/src/unstagger_A_fields.cc

// Include needed ETK & C library header files:
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// Needed for dealing with Cactus/ETK infrastructure
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void unstagger_A_fields(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Set Ai_unstaggered = Ai and exit the function if A fields are unstaggered already.
    if(A_fields_are_staggered == 0) {
#pragma omp parallel for
        for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
            int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
            Ax_unstaggered[index] = Ax[index];
            Ay_unstaggered[index] = Ay[index];
            Az_unstaggered[index] = Az[index];
        }
        return;
    }
    printf("Unstaggering A fields on grid with dx = %e!\n",CCTK_DELTA_SPACE(0));
    // If A fields are staggered (IllinoisGRMHD-style), then unstagger them:
    // First unstagger A_x, which is defined at (i, j+1/2, k+1/2). Unstaggering
    //   is as simple as A_x(i,j,k) = 1/4 * (A_x(i,j-1/2,k-1/2)+A_x(i,j-1/2,k+1/2)+A_x(i,j+1/2,k-1/2)+A_x(i,j+1/2,k+1/2))
#pragma omp parallel for
    for(int k=1;k<cctk_lsh[2];k++) for(int j=1;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        Ax_unstaggered[index] = 0.25*(Ax[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Ax[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] +
                                      Ax[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] + Ax[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
    }
#pragma omp parallel for
    for(int k=1;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=1;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        Ay_unstaggered[index] = 0.25*(Ay[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Ay[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] +
                                      Ay[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] + Ay[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
    }
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=1;j<cctk_lsh[1];j++) for(int i=1;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        Az_unstaggered[index] = 0.25*(Az[CCTK_GFINDEX3D(cctkGH,i,j,k)]     + Az[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] +
                                      Az[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] + Az[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
    }
}
```

    Writing interp_arbgrid_MO_ETK/src/unstagger_A_fields.cc


Next we instruct NRPy+ to interpolate the unstaggered gridfunctions.


```python
Ax_unstaggered = gri.register_gridfunctions("AUX","Ax_unstaggered")
Ay_unstaggered = gri.register_gridfunctions("AUX","Ay_unstaggered")
Az_unstaggered = gri.register_gridfunctions("AUX","Az_unstaggered")

gf_interp_list.append(gf_interp("Unstaggered vector potential component Ax"))
interp_expr = Ax_unstaggered
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Unstaggered vector potential component Ay"))
interp_expr = Ay_unstaggered
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Unstaggered vector potential component Az"))
interp_expr = Az_unstaggered
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy4metric'></a>

### Step 2.a.iii: Compute all 10 components of the 4-metric $g_{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{nrpy4metric}$$

We are given $\gamma_{ij}$, $\alpha$, and $\beta^i$ from ADMBase, and the 4-metric is given in terms of these quantities as
$$
g_{\mu\nu} = \begin{pmatrix} 
-\alpha^2 + \beta^k \beta_k & \beta_i \\
\beta_j & \gamma_{ij}
\end{pmatrix}.
$$


```python
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
        gf_interp_list.append(gf_interp("4-metric component g4DD"+str(mu)+str(nu)))
        interp_expr = g4DD[mu][nu]
        which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy4christoffels_cartesian'></a>

### Step 2.a.iv: Compute all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$ in Cartesian coordinates \[Back to [top](#toc)\]
$$\label{nrpy4christoffels_cartesian}$$

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
betaDdD = ixp.zerorank2()

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
            betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

# Eq. 2.122 in B&S
g4DDdD = ixp.zerorank3(DIM=4)
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
            gf_interp_list.append(gf_interp("4-Christoffel GammaUDD"+str(mu)+str(nu)+str(delta)))
            interp_expr = Gamma4UDD[mu][nu][delta]
            which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy4christoffels_spherical'></a>

### Step 2.a.v: Notes on computing all 40 4-Christoffels $\Gamma^{\mu}_{\nu\delta}$ in the Spherical basis \[Back to [top](#toc)\]
$$\label{nrpy4christoffels_spherical}$$

As explained in [Eq. 3.15 of Carroll's lecture notes on GR](https://ned.ipac.caltech.edu/level5/March01/Carroll3/Carroll3.html), while connection coefficients (a.k.a. Christoffel symbols) are not tensors, differences in connection coefficients are tensors.

Thus we may define

$$
\Delta^\mu_{\nu\delta} = \Gamma^\mu_{\nu\delta} - \hat{\Gamma}^\mu_{\nu\delta},
$$

where for example $\Gamma^\mu_{\nu\delta}$ is the connection related to the curved spacetime 4-metric in some basis and $\hat{\Gamma}^\mu_{\nu\delta}$ is the connection related to the flat spacetime 4-metric in the same basis.

We are given the 4-metric data in Cartesian coordinates, for which $\hat{\Gamma}^\mu_{\nu\delta}=0$. The basis transform to spherical coordinates is then straightforward:

\begin{align}
\Delta^\mu_{\text{Sph}\ \nu\delta} &= \Gamma^\mu_{\text{Sph}\ \nu\delta} - \hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta} \\
&= \frac{\partial x^\mu_{\rm Sph}}{\partial x^\alpha_{\rm Cart}}
\frac{\partial x^\beta_{\rm Cart}}{\partial x^\nu_{\rm Sph}}
\frac{\partial x^\gamma_{\rm Cart}}{\partial x^\delta_{\rm Sph}} \Delta^\alpha_{\text{Cart}\ \beta\gamma} \\
&= \frac{\partial x^\mu_{\rm Sph}}{\partial x^\alpha_{\rm Cart}}
\frac{\partial x^\beta_{\rm Cart}}{\partial x^\nu_{\rm Sph}}
\frac{\partial x^\gamma_{\rm Cart}}{\partial x^\delta_{\rm Sph}} \Gamma^\alpha_{\text{Cart}\ \beta\gamma} \\
\implies \Gamma^\mu_{\text{Sph}\ \nu\delta} &= \frac{\partial x^\mu_{\rm Sph}}{\partial x^\alpha_{\rm Cart}}
\frac{\partial x^\beta_{\rm Cart}}{\partial x^\nu_{\rm Sph}}
\frac{\partial x^\gamma_{\rm Cart}}{\partial x^\delta_{\rm Sph}} \Gamma^\alpha_{\text{Cart}\ \beta\gamma} +
\hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta}
\end{align}

**Define $\hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta}$.**

By definition,
$$
\hat{\Gamma}^{\mu}_{\nu\delta} = \frac{1}{2} \hat{g}^{\mu\eta} \left(\hat{g}_{\eta\nu,\delta} + \hat{g}_{\eta\delta,\nu} - \hat{g}_{\nu\delta,\eta}  \right).
$$

In static spherical coordinates, $\hat{g}_{\nu\delta}$ is given by 

$$
\hat{g}_{\mu\nu} = \begin{pmatrix} 
-1 & 0 \\
0 & \hat{\gamma}_{ij}
\end{pmatrix},
$$
so the inverse is easy to compute:
$$
\hat{g}^{\mu\nu} = \begin{pmatrix} 
-1 & 0 \\
0 & 1/\hat{\gamma}_{ij}
\end{pmatrix}.
$$
Here is the NRPy+ code implementation of $\hat{g}_{\mu\nu}$, $\hat{g}^{\mu\nu}$, and $\hat{g}_{\eta\nu,\delta}$:


```python
# import reference_metric as rfm
# # Set the desired *output* coordinate system to Spherical:
# #par.set_parval_from_str("reference_metric::CoordSystem","NobleSphericalThetaOptionOne")
# par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
# print("Calling reference_metric()...")
# rfm.reference_metric()
# print("Just finished calling reference_metric()...")


# g4hatDD = ixp.zerorank2(DIM=4)
# g4hatUU = ixp.zerorank2(DIM=4)
# g4hatDD[0][0] = sp.sympify(-1)
# g4hatUU[0][0] = sp.sympify(-1)
# for j in range(3):
#     g4hatDD[j+1][j+1] = rfm.ghatDD[j][j]
#     g4hatUU[j+1][j+1] = 1/rfm.ghatDD[j][j]
# g4hatDDdD = ixp.zerorank3(DIM=4)
# for eta in range(4):
#     for nu in range(4):
#         for j in range(3): # Time derivatives are all zero, so g4hatDDdD[eta][nu][0] = 0 (as initialized).
#             g4hatDDdD[eta][nu][j+1] = sp.diff(g4hatDD[eta][nu],rfm.xx[j])
```

Next we compute the 4-Christoffels $\hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta}$.


```python
# Gamma4hatSphUDD = ixp.zerorank3(DIM=4)
# for mu in range(4):
#     for nu in range(4):
#         for delta in range(4):
#             for eta in range(4):
#                 Gamma4hatSphUDD[mu][nu][delta] += sp.Rational(1,2)*g4hatUU[mu][eta]* \
#                 ( g4hatDDdD[eta][nu][delta] + g4hatDDdD[eta][delta][nu] - g4hatDDdD[nu][delta][eta] )

# # Here are the results, cf. Eq 18 of https://arxiv.org/pdf/1211.6632.pdf
# sp.pretty_print(Gamma4hatSphUDD)
```

Finally, compute $\Gamma^\mu_{\text{Sph}\ \nu\delta}$. Recall from above that
\begin{align}
\Gamma^\mu_{\text{Sph}\ \nu\delta} &= \frac{\partial x^\mu_{\rm Sph}}{\partial x^\alpha_{\rm Cart}}
\frac{\partial x^\beta_{\rm Cart}}{\partial x^\nu_{\rm Sph}}
\frac{\partial x^\gamma_{\rm Cart}}{\partial x^\delta_{\rm Sph}} \Gamma^\alpha_{\text{Cart}\ \beta\gamma} +
\hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta}
\end{align}

<a id='nrpybasisxform'></a>

### Step 2.a.vi:  Notes on basis transforming all Cartesian basis quantities to spherical \[Back to [top](#toc)\]
$$\label{nrpybasisxform}$$

All tensors and vectors are in the Cartesian coordinate basis $x^i_{\rm Cart} = (x,y,z)$, but we need them in the curvilinear coordinate basis $x^i_{\rm rfm}$=`(xx0,xx1,xx2)`=$(r,\theta,\phi)$ set by NRPy+'s `"reference_metric::CoordSystem"` variable (we'll set this parameter to `"Spherical"`). 

Empirically speaking, it is usually easier to write `(x(xx0,xx1,xx2),y(xx0,xx1,xx2),z(xx0,xx1,xx2))` than the inverse, so we will compute the Jacobian matrix

$$
{\rm Jac\_dUSph\_dDrfmUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm rfm}},
$$

via exact differentiation (courtesy SymPy), and the inverse Jacobian
$$
{\rm Jac\_dUrfm\_dDSphUD[i][j]} = \frac{\partial x^i_{\rm rfm}}{\partial x^j_{\rm Cart}},
$$

using NRPy+'s `generic\_matrix\_inverter3x3()` function. In terms of these, the transformation of vectors and rank-2 fully covariant tensors from Cartesian to `"reference_metric::CoordSystem"` (Spherical) coordinates may be written:

\begin{align}
g^{\rm rfm}_{\mu\nu} &= 
\frac{\partial x^{\alpha}_{\rm Cart}}{\partial x^{\mu}_{\rm rfm}}
\frac{\partial x^{\beta}_{\rm Cart}}{\partial x^{\nu}_{\rm rfm}} g^{\rm Cart}_{\alpha \beta} \\
\Gamma^\mu_{\text{Sph}\ \nu\delta} &= \frac{\partial x^\mu_{\rm Sph}}{\partial x^\alpha_{\rm Cart}}
\frac{\partial x^\beta_{\rm Cart}}{\partial x^\nu_{\rm Sph}}
\frac{\partial x^\gamma_{\rm Cart}}{\partial x^\delta_{\rm Sph}} \Gamma^\alpha_{\text{Cart}\ \beta\gamma} +
\hat{\Gamma}^\mu_{\text{Sph}\ \nu\delta}
\end{align}


```python
# Jac_dUCart_dDrfmUD = ixp.zerorank2()
# for i in range(DIM):
#     for j in range(DIM):
#         Jac_dUCart_dDrfmUD[i][j] = sp.simplify(sp.diff(rfm.xx_to_Cart[i],rfm.xx[j]))

# Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)
# Jac4_dUCart_dDrfmUD = ixp.zerorank2(DIM=4)
# Jac4_dUrfm_dDCartUD = ixp.zerorank2(DIM=4)
# for alp in range(4):
#     for bet in range(4):
#         if alp==0 or bet==0:
#             Jac4_dUCart_dDrfmUD[alp][bet] = sp.sympify(1) # Time components unchanged
#             Jac4_dUrfm_dDCartUD[alp][bet] = sp.sympify(1) # Time components unchanged
#         else:
#             Jac4_dUCart_dDrfmUD[alp][bet] = sp.simplify(Jac_dUCart_dDrfmUD[alp-1][bet-1])
#             Jac4_dUrfm_dDCartUD[alp][bet] = sp.simplify(Jac_dUrfm_dDCartUD[alp-1][bet-1])

# Gamma4SphUDD = ixp.zerorank3(DIM=4)
# for mu in range(4):
#     for nu in range(4):
#         for delt in range(4):
#             Gamma4SphUDD[mu][nu][delt] = Gamma4hatSphUDD[mu][nu][delt]
#             for alp in range(4):
#                 for bet in range(4):
#                     for gam in range(4):
#                         Gamma4SphUDD[mu][nu][delt] += \
#                 Jac4_dUrfm_dDCartUD[mu][alp]*Jac4_dUCart_dDrfmUD[bet][nu]*Jac4_dUCart_dDrfmUD[gam][delt] * \
#                           Gamma4UDD[alp][bet][gam]

# # Now output the Spherical 4-Christoffels to file:
# for mu in range(4):
#     for nu in range(4):
#         for delt in range(nu,4):
#             gf_interp_list.append(gf_interp("4-Christoffel component in SPHERICAL BASIS: GammaSphUDD"+str(mu)+str(nu)+str(delt)))
#             interp_expr = Gamma4SphUDD[mu][nu][delt]
#             which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

Output the 4-metric in Spherical coordinates to file:

\begin{align}
g^{\rm rfm}_{\mu\nu} &= 
\frac{\partial x^{\alpha}_{\rm Cart}}{\partial x^{\mu}_{\rm rfm}}
\frac{\partial x^{\beta}_{\rm Cart}}{\partial x^{\nu}_{\rm rfm}} g^{\rm Cart}_{\alpha \beta}
\end{align}


```python
# g4SphDD = ixp.zerorank2(DIM=4)
# for mu in range(4):
#     for nu in range(4):
#         for alp in range(4):
#             for bet in range(4):
#                 g4SphDD[mu][nu] += Jac4_dUCart_dDrfmUD[alp][mu]*Jac4_dUCart_dDrfmUD[bet][nu]*g4DD[alp][bet]

# for mu in range(4):
#     for nu in range(mu,4):
#         gf_interp_list.append(gf_interp("4-metric component in SPHERICAL BASIS: g4SphDD"+str(mu)+str(nu)))
#         interp_expr = g4SphDD[mu][nu]
#         which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

Next output the various GRMHD 3-vectors in the Spherical basis

\begin{align}
v^i_{\rm Sph} &= 
\frac{\partial x^{i}_{\rm Sph}}{\partial x^{j}_{\rm Cart}} v^j_{\rm Cart}
\end{align}


```python
# IGMvSphU                  = ixp.zerorank1()
# ValenciavSphU             = ixp.zerorank1()
# Gamma_times_ValenciavSphU = ixp.zerorank1()
# BSphU                     = ixp.zerorank1()

# for i in range(DIM):
#     for j in range(DIM):
#         IGMvSphU[i]                  += Jac_dUrfm_dDCartUD[i][j] * IGMvU[j]
#         ValenciavSphU[i]             += Jac_dUrfm_dDCartUD[i][j] * ValenciavU[j]
#         Gamma_times_ValenciavSphU[i] += Jac_dUrfm_dDCartUD[i][j] * Gamma_times_ValenciavU[j]
#         BSphU[i]                     += Jac_dUrfm_dDCartUD[i][j] * BU[j]

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("IGM 3-velocity vU"+str(i)+" = u^i/u^0 in SPHERICAL BASIS"))
#     interp_expr = IGMvSphU[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("Valencia 3-velocity vU"+str(i)+" in SPHERICAL BASIS"))
#     interp_expr = ValenciavSphU[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("Lorentz factor, times Valencia vU"+str(i)+" in SPHERICAL BASIS"))
#     interp_expr = Gamma_times_ValenciavSphU[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

# for i in range(DIM):
#     gf_interp_list.append(gf_interp("IGM magnetic field component B"+str(i)+" in SPHERICAL BASIS"))
#     interp_expr = BSphU[i]
#     which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

 <a id='psi4andfriends'></a>

### Step 2.a.vii: ): Output Weyl scalars $\psi_0$ through $\psi_4$, as well as Weyl invariants $J$ and $I$, from the `WeylScal4` ETK thorn \[Back to [top](#toc)\]
$$\label{psi4andfriends}$$


```python
Weylgfs = ["Psi0r","Psi0i","Psi1r","Psi1i","Psi2r","Psi2i","Psi3r","Psi3i","Psi4r","Psi4i",
                                    "curvIr","curvIi","curvJr","curvJi"]
Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,curvIr,curvIi,curvJr,curvJi = \
    gri.register_gridfunctions("AUX",Weylgfs)
count = 0
for gf in [Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,curvIr,curvIi,curvJr,curvJi]:
    gf_interp_list.append(gf_interp("4-Weyl scalar or invariant "+Weylgfs[count]))
    interp_expr = gf
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
    count = count + 1
```

<a id='constraints_gfs'></a>

### Step 2.a.viii: Hamiltonian and momentum constraints \[Back to [top](#toc)\]
$$\label{constraints_gfs}$$

We now setup the code to output the Hamiltonian and momentum constraint gridfunctions, that is:

* $\mathcal{H}$ (e.g., $\text{ML_BSSN::H}$ or $\text{Baikal::HGF}$)
* $\mathcal{M}^{0}$ (e.g., $\text{ML_BSSN::M1}$ or $\text{Baikal::MU0GF}$)
* $\mathcal{M}^{1}$ (e.g., $\text{ML_BSSN::M2}$ or $\text{Baikal::MU1GF}$)
* $\mathcal{M}^{2}$ (e.g., $\text{ML_BSSN::M3}$ or $\text{Baikal::MU2GF}$)


```python
gf_interp_list.append(gf_interp("Constraints - Hamiltonian"))
H = gri.register_gridfunctions("AUX","InterpH")
interp_expr = H
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

MU = ixp.register_gridfunctions_for_single_rank1("AUX","InterpMU")

# Sph_basis_MU = ixp.zerorank1()
# for i in range(DIM):
#     for l in range(DIM):
#         Sph_basis_MU[i] += Jac_dUrfm_dDCartUD[i][l] * MU[l]

for i in range(DIM):
    gf_interp_list.append(gf_interp("Constraints - Momentum M^"+chr(ord('x')+i)))
#     interp_expr = Sph_basis_MU[i]
    interp_expr = MU[i]
    which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nuc_eos_gfs'></a>

### Step 2.a.ix: Nuclear (tabulated) equation of state quantities \[Back to [top](#toc)\]
$$\label{nuc_eos_gfs}$$

When using nuclear (tabulated) equations of state, additional gridfunctions are available for output. These are:

* $T$, the temperature (e.g. $\text{HydroBase::temperature}$ or $\text{IllinoisGRMHD::igm_temperature}$)
* $Y_{e}$, the electron fraction (e.g. $\text{HydroBase::Y_e}$ or $\text{IllinoisGRMHD::igm_Ye}$)
* $\epsilon$, the specific internal energy (e.g. $\text{HydroBase::eps}$ or $\text{IllinoisGRMHD::igm_eps}$)
* $S$, the entropy (e.g. $\text{HydroBase::entropy}$ or $\text{IllinoisGRMHD::igm_entropy}$)


```python
gf_interp_list.append(gf_interp("Temperature primitive"))
temperature = gri.register_gridfunctions("AUX","temperature")
interp_expr = temperature
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Electron fraction primitive"))
Ye = gri.register_gridfunctions("AUX","Y_e")
interp_expr = Ye
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Specific internal energy primitive"))
eps = gri.register_gridfunctions("AUX","eps")
interp_expr = eps
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)

gf_interp_list.append(gf_interp("Entropy primitive"))
entropy = gri.register_gridfunctions("AUX","entropy")
interp_expr = entropy
which_InterpCounter = interp_fileout(which_InterpCounter,interp_expr,NRPyoutfilename)
```

<a id='nrpy_c_calling_function'></a>

## Step 2.b: C code calling function for the NRPy+ C output \[Back to [top](#toc)\]
$$\label{nrpy_c_calling_function}$$

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
void list_of_functions_to_interpolate(cGH *restrict cctkGH,
                                      const CCTK_INT *restrict cctk_lsh,
                                      const CCTK_INT *restrict cctk_nghostzones,
                                      const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
                                      const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,
                                      const CCTK_INT *restrict InterpCounter,
                                      const CCTK_REAL *restrict rho_bGF,
                                      const CCTK_REAL *restrict PGF,
                                      const CCTK_REAL *restrict temperatureGF,
                                      const CCTK_REAL *restrict Y_eGF,
                                      const CCTK_REAL *restrict epsGF,
                                      const CCTK_REAL *restrict entropyGF,
                                      const CCTK_REAL *restrict IGMvU0GF,const CCTK_REAL *restrict IGMvU1GF,const CCTK_REAL *restrict IGMvU2GF,
                                      const CCTK_REAL *restrict BU0GF,const CCTK_REAL *restrict BU1GF,const CCTK_REAL *restrict BU2GF,
                                      const CCTK_REAL *restrict gammaDD00GF,const CCTK_REAL *restrict gammaDD01GF,const CCTK_REAL *restrict gammaDD02GF,
                                      const CCTK_REAL *restrict gammaDD11GF,const CCTK_REAL *restrict gammaDD12GF,const CCTK_REAL *restrict gammaDD22GF,
                                      const CCTK_REAL *restrict betaU0GF,const CCTK_REAL *restrict betaU1GF,const CCTK_REAL *restrict betaU2GF,
                                      const CCTK_REAL *restrict alphaGF,
                                      CCTK_REAL *restrict interped_gfGF,
                                      CCTK_REAL *restrict Ax_unstaggeredGF,CCTK_REAL *restrict Ay_unstaggeredGF,CCTK_REAL *restrict Az_unstaggeredGF,
                                      const CCTK_REAL *restrict Psi0rGF,const CCTK_REAL *restrict Psi0iGF,
                                      const CCTK_REAL *restrict Psi1rGF,const CCTK_REAL *restrict Psi1iGF,
                                      const CCTK_REAL *restrict Psi2rGF,const CCTK_REAL *restrict Psi2iGF,
                                      const CCTK_REAL *restrict Psi3rGF,const CCTK_REAL *restrict Psi3iGF,
                                      const CCTK_REAL *restrict Psi4rGF,const CCTK_REAL *restrict Psi4iGF,
                                      const CCTK_REAL *restrict curvIrGF,const CCTK_REAL *restrict curvIiGF,
                                      const CCTK_REAL *restrict curvJrGF,const CCTK_REAL *restrict curvJiGF,
                                      const CCTK_REAL *restrict InterpHGF,
                                      const CCTK_REAL *restrict InterpMU0GF,
                                      const CCTK_REAL *restrict InterpMU1GF,
                                      const CCTK_REAL *restrict InterpMU2GF) {
  DECLARE_CCTK_PARAMETERS;
#include "list_of_functions_to_interpolate.h"
}

void construct_function_to_interpolate__store_to_interped_gf(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  printf("Called construct_function_to_interpolate__store_to_interped_gf() on grid with dx = %e!\n",CCTK_DELTA_SPACE(0));
  const CCTK_REAL invdx0 = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdx1 = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdx2 = 1.0 / CCTK_DELTA_SPACE(2);

  // .------------------------------------.
  // | Constraint equations gridfunctions |
  // .------------------------------------.
  // We get the pointer using the CCTK_VarDataPtr() function because
  // if the gridfunction does not exist in this version of IGM then
  // we will get compilation errors.
  int timelevel = 0;

  // Hamiltonian constraint gridfunction
  CCTK_REAL *InterpHGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,HVarString));
  if( !InterpHGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Hamiltonian constraint - Couldn't get data pointer of input array variable '%s'", HVarString);

  // Momentum constraint gridfunction (x-direction)
  CCTK_REAL *InterpMU0GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M0VarString));
  if( !InterpMU0GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (x-direction) - Couldn't get data pointer of input array variable '%s'", M0VarString);

  // Momentum constraint gridfunction (y-direction)
  CCTK_REAL *InterpMU1GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M1VarString));
  if( !InterpMU1GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (y-direction) - Couldn't get data pointer of input array variable '%s'", M1VarString);

  // Momentum constraint gridfunction (z-direction)
  CCTK_REAL *InterpMU2GF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,M2VarString));
  if( !InterpMU2GF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Momentum constraint (z-direction) - Couldn't get data pointer of input array variable '%s'", M2VarString);

  // .-------------------------------------------------------.
  // | Additional hydro quantities for nuclear/tabulated EOS |
  // .-------------------------------------------------------.
  CCTK_REAL *epsGF, *temperatureGF, *Y_eGF, *entropyGF;
  if( enable_nuc_eos ) {
    // Temperature gridfunction
    temperatureGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,temperatureVarString));
    if( !temperatureGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "temperature - Couldn't get data pointer of input array variable '%s'", temperatureVarString);

    // Electron fraction gridfunction
    Y_eGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,Y_eVarString));
    if( !Y_eGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "Y_e - Couldn't get data pointer of input array variable '%s'", Y_eVarString);

    // Specific internal energy gridfunction
    epsGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,epsVarString));
    if( !epsGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "eps - Couldn't get data pointer of input array variable '%s'", epsVarString);

    // Entropy gridfunction
    entropyGF = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel,entropyVarString));
    if( !entropyGF ) CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, "entropy - Couldn't get data pointer of input array variable '%s'", entropyVarString);
  }
  else {
    // If we are not using nuclear EOS, then set the gridfunction pointers to NULL.
    temperatureGF = NULL;
    Y_eGF         = NULL;
    epsGF         = NULL;
    entropyGF     = NULL;
  }

  list_of_functions_to_interpolate(cctkGH,cctk_lsh,cctk_nghostzones,
                                   x,y,z,
                                   invdx0,invdx1,invdx2,
                                   InterpCounter,
                                   rho_b,P,temperatureGF,Y_eGF,epsGF,entropyGF,
                                   vx,vy,vz,
                                   Bx,By,Bz,
                                   gxx,gxy,gxz,gyy,gyz,gzz,
                                   betax,betay,betaz,alp, interped_gf,
                                   Ax_unstaggered,Ay_unstaggered,Az_unstaggered,
                                   Psi0r,Psi0i,Psi1r,Psi1i,Psi2r,Psi2i,Psi3r,Psi3i,Psi4r,Psi4i,
                                   curvIr,curvIi,curvJr,curvJi,
                                   InterpHGF,InterpMU0GF,InterpMU1GF,InterpMU2GF);

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

    Writing interp_arbgrid_MO_ETK/src/construct_function_to_interpolate__store_to_interped_gf.cc


<a id='nrpygetgfname'></a>

## Step 2.c: The `get_gf_name()` function \[Back to [top](#toc)\]

$$\label{nrpygetgfname}$$


```python
with open(os.path.join(Ccodesdir,"src","get_gf_name.h"), "w") as file:
    file.write("void get_gf_name(const int InterpCounter,char gf_name[100]) {\n")
    for i in range(1,which_InterpCounter):
        file.write("    if(InterpCounter=="+str(i)+") { snprintf(gf_name,100,\""+gf_interp_list[i].gf_description+"\"); return; }\n")
    file.write("    printf(\"Error. InterpCounter = %d unsupported. I should not be here.\\n\",InterpCounter); exit(1);\n")
    file.write("}\n")
```

<a id='nrpy_interp_counter'></a>

## Step 2.d: C Code for Initializing and incrementing `InterpCounter` \[Back to [top](#toc)\]
$$\label{nrpy_interp_counter}$$

The gridfunctions are interpolated one at a time based on the current value of the index quantity `InterpCounter`. Here we write the C code needed for initializing and incrementing this variable.


```python
with open(os.path.join(Ccodesdir,"src","define_NumInterpFunctions.h"), "w") as file:
    file.write("#define NumInterpFunctions ("+str(which_InterpCounter)+"-4*(1-enable_nuc_eos))\n")
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

void ArbGrid_InitializeInterpCounterToZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  *InterpCounter = 0;

  if(verbose==2) printf("interp_arbgrid_MO_ETK: Just set InterpCounter to %d\n",*InterpCounter);
}

void ArbGrid_InitializeInterpCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration == interp_out_iteration) {
    *InterpCounter = 1;
    if(verbose==2) printf("interp_arbgrid_MO_ETK: Just set InterpCounter to %d ; ready to start looping over interpolated gridfunctions!\n",
                          *InterpCounter);
  }
}

// This function increments InterpCounter if we are at the interp_out_iteration until
// it hits NumInterpFunctions. At this iteration, InterpCounter is set to zero, which
// exits the loop.
void ArbGrid_IncrementInterpCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(*InterpCounter == NumInterpFunctions-1) {
    *InterpCounter = 0;
    if(verbose==2) printf("interp_arbgrid_MO_ETK: Finished! Just zeroed InterpCounter.\n");
  } else {
    (*InterpCounter)++;
    if(verbose==2) printf("interp_arbgrid_MO_ETK: Just incremented InterpCounter to %d of %d\n",*InterpCounter,NumInterpFunctions-1);
  }
}
```

    Writing interp_arbgrid_MO_ETK/src/interp_counter.cc


<a id='validationagainstfm'></a>

# Step 2.e: Validation of interpolated data against exact Fishbone-Moncrief data \[Back to [top](#toc)\]
$$\label{validationagainstfm}$$




```python
# # Step 1c: Call the FishboneMoncriefID() function from within the
# #          FishboneMoncriefID/FishboneMoncriefID.py module.
# import FishboneMoncriefID.FishboneMoncriefID as fmid

# old_glb_gridfcs_list = gri.glb_gridfcs_list
# # Step 1: Set up the Fishbone-Moncrief initial data. This sets all the ID gridfunctions.
# gri.glb_gridfcs_list = [] # Reset list of gridfunctions
# fmid.FishboneMoncriefID("Spherical")
# gammaDD = ixp.zerorank2()

# DIM = 3
# for i in range(DIM):
#     for j in range(DIM):
#         if i<=j:
#             gammaDD[i][j] = fmid.IDgammaDD[i][j]
#         else:
#             gammaDD[i][j] = fmid.IDgammaDD[j][i]

# # gamma_{ij} v^i_{(n)} v^j_{(n)}
# Gammacontraction = sp.sympify(0)
# for i in range(DIM):
#     for j in range(DIM):
#         Gammacontraction += gammaDD[i][j] * fmid.IDValencia3velocityU[i] * fmid.IDValencia3velocityU[j]

# Gammafactor = sp.sqrt(1 / (1 - Gammacontraction))

# # -={ F-M quantities: Generate C code from expressions and output to file }=-
# FishboneMoncrief_to_print = [\
#                      lhrh(lhs="Gammafactor",rhs=Gammafactor),\
#                      lhrh(lhs="Gamma_times_ValenciavU0",rhs=Gammafactor*fmid.IDValencia3velocityU[0]),\
#                      lhrh(lhs="Gamma_times_ValenciavU1",rhs=Gammafactor*fmid.IDValencia3velocityU[1]),\
#                      lhrh(lhs="Gamma_times_ValenciavU2",rhs=Gammafactor*fmid.IDValencia3velocityU[2]),\
#                      ]
# fin.FD_outputC(os.path.join(Ccodesdir,"src","FM_Gamma__Gamma_times_Valenciavs_sphbasis.h"),FishboneMoncrief_to_print,
#                params="outCverbose=False,CSE_enable=True")

# # Restore old gridfunctions list:
# gri.glb_gridfcs_list = old_glb_gridfcs_list
```


```python
# %%writefile $Ccodesdir/src/FM_validation.cc

# #include <assert.h>
# #include <stdio.h>
# #include <stdlib.h>
# #include <string.h>
# #include <math.h>
# #include <ctype.h>
# // Needed for dealing with Cactus/ETK infrastructure
# #include "cctk.h"
# #include "cctk_Arguments.h"
# #include "cctk_Parameters.h"
# // Needed for low-level interpolation functions
# #include "util_Table.h"
# #include "util_String.h"

# // C++ function prototypes:
# extern void Interpolate_to_dest_grid(const cGH *cctkGH,const CCTK_INT interp_num_points, const CCTK_INT interp_order,
#                                      const CCTK_REAL *point_x_temp,const CCTK_REAL *point_y_temp,const CCTK_REAL *point_z_temp,
#                                      const CCTK_STRING input_array_names[1], CCTK_REAL *output_f[1]);
# extern void get_gf_name(const int InterpCounter,char gf_name[100]);

# #define FREE_2D_GENERIC(type,array,ni,nj) for(int cc = 0; cc < ni;cc++) free((void *)array[cc]); \
# /**/                                                                    free((void *)array);

# void FM_validation(CCTK_ARGUMENTS)
# {
#     DECLARE_CCTK_ARGUMENTS;
#     DECLARE_CCTK_PARAMETERS;

#     const CCTK_INT sph_Nr  = 3200;
#     const CCTK_INT sph_Nth = 1;
#     const CCTK_INT sph_Nph = 160;
#     const CCTK_REAL sph_rmin = 0.1;
#     const CCTK_REAL sph_rmax = 50.0;
#     const CCTK_REAL sph_thmin = M_PI/2.0;
#     const CCTK_REAL sph_thmax = M_PI/2.0;
#     const CCTK_REAL sph_phmin = 0;
#     const CCTK_REAL sph_phmax = 2.0*M_PI;

#     const CCTK_INT num_interp_points = sph_Nr*sph_Nth*sph_Nph;

#   // STEP 1: IF GAMMA*VALENCIA,PROCEED. IF RHO_B, OUTPUT FM FOR VEL DATA. ELSE RETURN.
#   // Perform interpolation only at iteration == interp_out_iteration:
#   if(cctk_iteration != interp_out_iteration) return;
#   char gf_name[100]; get_gf_name(*InterpCounter,gf_name);

# //  if(strncmp(gf_name,"Lorentz factor, times Valencia",30) == 0) {
#   if(0 == 0) {
#       // Perform interpolation!
#       // Process zero (CCTK_MyProc(cctkGH)==0) is responsible for directing the interpolation.
#       //    All other processes must see the cctk_InterpGridArrays() within Interpolate_to_dest_grid(),
#       //    so that the MPI calls work properly, but these nonzero processes can call
#       //    Interpolate_to_dest_grid() with number of interpolated points set to zero, and
#       //    without needing a malloc().

#       if(CCTK_MyProc(cctkGH)==0) {
#         CCTK_REAL *points_x,*points_y,*points_z,**output_f;
#         // The name of the input gridfunction is always "interp_arbgrid_MO_ETK::interped_gf":
#         const CCTK_STRING input_array_names[1] = { "interp_arbgrid_MO_ETK::interped_gf" };

#         // STEP 1: Construct list of desired interpolation destination points:
#         points_x = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_interp_points);
#         points_y = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_interp_points);
#         points_z = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*num_interp_points);
#         output_f = (CCTK_REAL **)malloc(1 * sizeof(CCTK_REAL *));
#         for(int cc = 0; cc < 1; cc++) output_f[cc]=(CCTK_REAL *)malloc(num_interp_points * sizeof(CCTK_REAL));

#         // STEP 2: ALLOCATE INTERPOLATION ARRAYS SET INTERPOLATION POINT ARRAYS
#         CCTK_INT pointcount = 0;
#         for(int ir=0;ir<sph_Nr;ir++) for(int ith=0;ith<sph_Nth;ith++) for(int iph=0;iph<sph_Nph;iph++) {
#             const CCTK_REAL r  = sph_rmin  + (CCTK_REAL)ir /((CCTK_REAL)sph_Nr ) * (sph_rmax - sph_rmin);
#             const CCTK_REAL th = sph_thmin + (CCTK_REAL)ith/((CCTK_REAL)sph_Nth) * (sph_thmax - sph_thmin);
#             const CCTK_REAL ph = sph_phmin + (CCTK_REAL)iph/((CCTK_REAL)sph_Nph) * (sph_phmax - sph_phmin);
#             points_x[pointcount] = r*sin(th)*cos(ph);
#             points_y[pointcount] = r*sin(th)*sin(ph);
#             points_z[pointcount] = r*cos(th);
#             pointcount++;
#         } // END for(int ir=0;ir<sph_Nr;ir++) for...

#         // STEP 3: Looping over interp order as desired, interpolate to destination points & output to file
#         for(int order=1;order<=4;order*=2) {
#             printf("ASCII FM Validation: %d pts; Interpolating\033[1m %s \033[0m... using interpolation order = %d\n",num_interp_points,gf_name,order);
#       //Interpolate_to_dest_grid(cGH *cctkGH,CCTK_INT interp_num_points, CCTK_INT interp_order,
#       //                        CCTK_REAL *point_x_temp,CCTK_REAL *point_y_temp,CCTK_REAL *point_z_temp,
#       //                        const CCTK_STRING input_array_names[1], CCTK_REAL *output_f[1])
#             Interpolate_to_dest_grid(cctkGH, num_interp_points, order,
#                                      points_x,points_y,points_z, input_array_names, output_f);

#             // Step 1.d.i: Sanity check -- check for bad point:
#         #pragma omp parallel for
#             for(int i=0;i<num_interp_points;i++) {
#                 if(output_f[0][i] > 1e20) {
#                     printf("BAD POINT: %s %d %e %e %e %e\n",gf_name,i,points_x[i],points_y[i],points_z[i], output_f[0][i]);
#                     exit(1);
#                 } // END if(output_f[0][i] > 1e20)
#             } // END for(int i=0;i<num_interp_points;i++)
#             char filename[500];
#             sprintf (filename, "%s/validation_points-%s-order%d.asc", out_dir,gf_name,order);
#             FILE *file;
#             file = fopen (filename,"w");
#             printf("WRITING to file %s\n",filename);
#             if (! file) {
#                 CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
#                             "interp_dest_grid__ET_thorn: Cannot open ASCII output file '%s'", filename);
#                 exit(1);
#             } // END if (! file)
#             pointcount = 0;
#             for(int ir=0;ir<sph_Nr;ir++) for(int ith=0;ith<sph_Nth;ith++) for(int iph=0;iph<sph_Nph;iph++) {
#                 const CCTK_REAL xx = points_x[pointcount];
#                 const CCTK_REAL yy = points_y[pointcount];
#                 const CCTK_REAL zz = points_z[pointcount];
#                 fprintf(file,"%e %e %e %e\n",xx,yy,zz,output_f[0][pointcount]);
#                 pointcount++;
#             }
#             fclose(file);
#         } // END for(int order=1;order<=4;order*=2)

#         // STEP 3: FREE THE MALLOCs for destination grids and interpolation output
#         free(points_x);
#         free(points_y);
#         free(points_z);
#         FREE_2D_GENERIC(CCTK_REAL,output_f,1,num_interp_points);

#       } else if(CCTK_MyProc(cctkGH)!=0) {
#         // On all MPI processes that are nonzero, only call the interpolation function

#         CCTK_REAL *points_x,*points_y,*points_z,**output_f;
#         // The name of the input gridfunction is always "interp_arbgrid_MO_ETK::interped_gf":
#         const CCTK_STRING input_array_names[1] = { "interp_arbgrid_MO_ETK::interped_gf" };

#         //    to ensure the MPI calls from the actual interpolation (driven by proc==0) are seen.
#         for(int order=1;order<=4;order*=2) {
#             Interpolate_to_dest_grid(cctkGH, 0, order,points_x,points_y,points_z, input_array_names, output_f);
#         } // END for(int order=1;order<=4;order*=2)
#       } // END if(CCTK_MyProc(cctkGH)...)

#   }
#     if(strncmp(gf_name,"IGM density primitive",21) == 0 && CCTK_MyProc(cctkGH)==0) {
#       char filename[500];
#       sprintf (filename, "%s/FMvalidation_points.asc", out_dir);
#       FILE *file;
#       file = fopen (filename,"w");
#       printf("WRITING to file %s\n",filename);

#       for(int ir=0;ir<sph_Nr;ir++) for(int ith=0;ith<sph_Nth;ith++) for(int iph=0;iph<sph_Nph;iph++) {
#         const CCTK_REAL xx0 = sph_rmin  + (CCTK_REAL)ir /((CCTK_REAL)sph_Nr ) * (sph_rmax - sph_rmin);
#         const CCTK_REAL xx1 = sph_thmin + (CCTK_REAL)ith/((CCTK_REAL)sph_Nth) * (sph_thmax - sph_thmin);
#         const CCTK_REAL xx2 = sph_phmin + (CCTK_REAL)iph/((CCTK_REAL)sph_Nph) * (sph_phmax - sph_phmin);
#         const CCTK_REAL xx = xx0*sin(xx1)*cos(xx2);
#         const CCTK_REAL yy = xx0*sin(xx1)*sin(xx2);
#         const CCTK_REAL zz = xx0*cos(xx1);
#         if(xx0 < r_in) {
#             fprintf(file,"%e %e %e %e %e %e %e\n",xx,yy,zz,
#                     0.0,0.0,0.0,1.0);
#         } else {
#             CCTK_REAL Gammafactor,Gamma_times_ValenciavU0,Gamma_times_ValenciavU1,Gamma_times_ValenciavU2;
#     #include "FM_Gamma__Gamma_times_Valenciavs_sphbasis.h"
#             fprintf(file,"%e %e %e %e %e %e %e\n",xx,yy,zz,
#                     Gamma_times_ValenciavU0,Gamma_times_ValenciavU1,Gamma_times_ValenciavU2,Gammafactor);
#         } // END if(xx0 < r_in)
#       } // END for(int ir=0;ir<sph_Nr;ir++) for...
#       fclose(file);
#       return;
#   } // END if(strncmp(gf_name,"Lorentz factor, times Valencia",30) == 0)
# } // END function
# #undef FREE_2D_GENERIC
```

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
# Main make.code.defn file for thorn interp_arbgrid_MO_ETK

# Source files in this directory
SRCS =  main_function.cc unstagger_A_fields.cc interp_counter.cc \
        construct_function_to_interpolate__store_to_interped_gf.cc # FM_validation.cc # <- For FishboneMoncriefID validation
```

    Writing interp_arbgrid_MO_ETK/src/make.code.defn


<a id='interfaceccl'></a>

## Step 3.b: `interface.ccl` \[Back to [top](#toc)\]
$$\label{interfaceccl}$$

Let's now write `interface.ccl`. The [official Einstein Toolkit (Cactus) documentation](http://einsteintoolkit.org/usersguide/UsersGuide.html) defines what must/should be included in an `interface.ccl` file [**here**](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-178000D2.2). 


```python
%%writefile $Ccodesdir/interface.ccl

# With "implements", we give our thorn its unique name.
implements: interp_arbgrid_MO_ETK

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits:   admbase IllinoisGRMHD Grid
inherits:   WeylScal4 # Needed for Weyl scalars psi4, psi3, psi..., and Weyl invariants I & J.
# For FM ID comparisons:
# inherits:   FishboneMoncriefID

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

CCTK_REAL unstaggered_A_fields type=GF timelevels=3 tags='Checkpoint="no"'
{
  Ax_unstaggered,Ay_unstaggered,Az_unstaggered
} "Unstaggered A-field components."


int InterpCounterVar type = SCALAR tags='checkpoint="no"'
{
  InterpCounter
} "Counter that keeps track of which function we are interpolating."
```

    Writing interp_arbgrid_MO_ETK/interface.ccl


<a id='paramccl'></a>

## Step 3.c: `param.ccl` \[Back to [top](#toc)\]
$$\label{paramccl}$$

We will now write the file `param.ccl`. This file allows the listed parameters to be set at runtime. We also give allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-183000D2.3). 


```python
%%writefile $Ccodesdir/param.ccl

# Output the interpolated data to the IO::out_dir directory:
shares: IO
USES STRING out_dir

# These parameters are used to output Hybrid EOS information
shares: IllinoisGRMHD
USES CCTK_INT neos
USES CCTK_REAL Gamma_th

# These parameters are used to output nuclear (tabulated) EOS information
shares: EOS_Omni
USES CCTK_STRING nuceos_table_name

# For FM ID comparisons:
# shares: FishboneMoncriefID
# USES KEYWORD M
# USES KEYWORD a
# USES KEYWORD r_in
# USES KEYWORD r_at_max_density

restricted:

########################################
# BASIC THORN STEERING PARAMETERS
CCTK_INT interp_out_iteration "Which iteration to interpolate to destination grids?" STEERABLE=ALWAYS
{
  0:* :: ""
} 960000

## Interpolator information
CCTK_INT verbose "Set verbosity level: 1=useful info; 2=moderately annoying (though useful for debugging)" STEERABLE=ALWAYS
{
  0:2 :: "0 = no output; 1=useful info; 2=moderately annoying (though useful for debugging)"
} 2

CCTK_INT A_fields_are_staggered "Are A fields staggered? 1 = yes; 0 = no. Default to yes." STEERABLE=ALWAYS
{
  0:1 :: ""
} 1

##########
# Cartesian position of center of output grid (usually center of BH).
CCTK_REAL x_center "x-position of center." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

CCTK_REAL y_center "y-position of center." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

CCTK_REAL z_center "z-position of center." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

##########
# Shift offset:
CCTK_REAL beta_offsetU0 "Offset to betax, to account for coordinate drift in x direction." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

CCTK_REAL beta_offsetU1 "Offset to betay, to account for coordinate drift in y direction." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

CCTK_REAL beta_offsetU2 "Offset to betaz, to account for coordinate drift in z direction." STEERABLE=ALWAYS
{
  *:* :: ""
} 0.0

CCTK_REAL out_of_bounds_interp_xyz "Do not interpolate points with fabs(xyz) > out_of_bounds_interp_xyz, where xyz are centered at xyz_center (usually center of BH). Fill dataset with NaN instead." STEERABLE=ALWAYS
{
  0:* :: "Any positive number"
} 1E100

###########
# Nuclear EOS options
CCTK_INT enable_nuc_eos "Use nuclear (tabulated/microphysics) equation of state? 1 = yes; 0 = no. Default to no." STEERABLE=ALWAYS
{
  0:1 :: ""
} 0

CCTK_STRING temperatureVarString "Temperature GF name. Defaults to IllinoisGRMHD's." STEERABLE=ALWAYS
{
  "IllinoisGRMHD::igm_temperature" :: "IllinoisGRMHD's temperature gridfunction name"
  "HydroBase::temperature"         :: "HydroBase's temperature gridfunction name"
  ".+"                             :: "Or use you can use your own thorn's temperature gridfunction name"
} "IllinoisGRMHD::igm_temperature"

CCTK_STRING Y_eVarString "Electron fraction GF name. Defaults to IllinoisGRMHD's." STEERABLE=ALWAYS
{
  "IllinoisGRMHD::igm_Ye" :: "IllinoisGRMHD's electron fraction gridfunction name"
  "HydroBase::Y_e"        :: "HydroBase's  electron fraction gridfunction name"
  ".+"                    :: "Or use you can use your own thorn's  electron fraction gridfunction name"
} "IllinoisGRMHD::igm_Ye"

CCTK_STRING epsVarString "Specific internal energy GF name. Defaults to IllinoisGRMHD's." STEERABLE=ALWAYS
{
  "IllinoisGRMHD::igm_eps" :: "IllinoisGRMHD's specific internal energy gridfunction name"
  "HydroBase::eps"         :: "HydroBase's specific internal energy gridfunction name"
  ".+"                     :: "Or use you can use your own thorn's specific internal energy gridfunction name"
} "IllinoisGRMHD::igm_eps"

CCTK_STRING entropyVarString "Entropy GF name. Defaults to IllinoisGRMHD's." STEERABLE=ALWAYS
{
  "IllinoisGRMHD::igm_entropy" :: "IllinoisGRMHD's entropy gridfunction name"
  "HydroBase::entropy"         :: "HydroBase's entropy gridfunction name"
  ".+"                         :: "Or use you can use your own thorn's entropy gridfunction name"
} "IllinoisGRMHD::igm_entropy"

CCTK_STRING HVarString "Hamiltonian constraint GF name. Defaults to Baika's." STEERABLE=ALWAYS
{
  "Baikal::HGF" :: "Baikal's Hamiltonian constraint gridfunction name"
  "ML_BSSN::H"  :: "ML_BSSN's Hamiltonian constraint gridfunction name"
  ".+"                         :: "Or use you can use your own thorn's entropy gridfunction name"
} "Baikal::HGF"

CCTK_STRING M0VarString "Momentum constraint (x-direction) Defaults to Baikal's." STEERABLE=ALWAYS
{
  "Baikal::MU0GF" :: "Baikal's Momentum constraint (x-direction) gridfunction name"
  "ML_BSSN::M1"   :: "ML_BSSN's Momentum constraint (x-direction) gridfunction name"
  ".+"                         :: "Or use you can use your own thorn's entropy gridfunction name"
} "Baikal::MU0GF"

CCTK_STRING M1VarString "Momentum constraint (y-direction) GF name. Defaults to Baikal's." STEERABLE=ALWAYS
{
  "Baikal::MU1GF" :: "Baikal's Momentum constraint (y-direction) gridfunction name"
  "ML_BSSN::M2"   :: "ML_BSSN's Momentum constraint (y-direction) gridfunction name"
  ".+"                         :: "Or use you can use your own thorn's entropy gridfunction name"
} "Baikal::MU1GF"

CCTK_STRING M2VarString "Momentum constraint (z-direction) GF name. Defaults to Baikal's." STEERABLE=ALWAYS
{
  "Baikal::MU2GF" :: "Baikal's Momentum constraint (z-direction) gridfunction name"
  "ML_BSSN::M3"   :: "ML_BSSN's Momentum constraint (z-direction) gridfunction name"
  ".+"                         :: "Or use you can use your own thorn's entropy gridfunction name"
} "Baikal::MU2GF"
```

    Writing interp_arbgrid_MO_ETK/param.ccl


<a id='scheduleccl'></a>

## Step 3.d: `schedule.ccl` \[Back to [top](#toc)\]
$$\label{scheduleccl}$$

Finally, we will write the file `schedule.ccl`; its official documentation is found [here](http://einsteintoolkit.org/usersguide/UsersGuidech12.html#x17-186000D2.4). 

This file declares storage for variables declared in the `interface.ccl` file and specifies when the various parts of the thorn will be run:


```python
%%writefile $Ccodesdir/schedule.ccl

STORAGE: interpolation_gf[3]
STORAGE: unstaggered_A_fields[3]
STORAGE: InterpCounterVar
# STORAGE: interp_pointcoords_and_output_arrays

#############################
SCHEDULE ArbGrid_InitializeInterpCounterToZero AT CCTK_INITIAL
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable to zero"

SCHEDULE ArbGrid_InitializeInterpCounterToZero AT CCTK_POST_RECOVER_VARIABLES
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable to zero"

SCHEDULE ArbGrid_InitializeInterpCounter before ArbGrid_InterpGroup AT CCTK_ANALYSIS
{
  LANG: C
  OPTIONS: GLOBAL
} "Initialize InterpCounter variable"
##################

SCHEDULE GROUP ArbGrid_InterpGroup AT CCTK_ANALYSIS BEFORE CarpetLib_printtimestats BEFORE CarpetLib_printmemstats AFTER Convert_to_HydroBase WHILE interp_arbgrid_MO_ETK::InterpCounter
{
} "Perform all interpolations. This group is only actually scheduled at cctk_iteration==interp_out_iteration."

SCHEDULE unstagger_A_fields in ArbGrid_InterpGroup before construct_function_to_interpolate__store_to_interped_gf
{
  STORAGE: unstaggered_A_fields[3]
  OPTIONS: GLOBAL,LOOP-LOCAL
  SYNC:    unstaggered_A_fields
  LANG: C
} "Unstagger A fields."

SCHEDULE construct_function_to_interpolate__store_to_interped_gf in ArbGrid_InterpGroup before DoSum
{
  STORAGE: interpolation_gf[3],InterpCounterVar
  OPTIONS: GLOBAL,LOOP-LOCAL
  SYNC:    interpolation_gf
  LANG: C
} "Construct the function to interpolate"

SCHEDULE Interpolate_to_dest_grid_main_function in ArbGrid_InterpGroup after construct_function_to_interpolate__store_to_interped_gf
{
  OPTIONS: GLOBAL
  LANG: C
} "Perform interpolation and output result to file."

# For FishboneMoncriefID validation only.
# SCHEDULE FM_validation in ArbGrid_InterpGroup after Interpolate_to_dest_grid_main_function
# {
#   OPTIONS: GLOBAL
#   LANG: C
# } "Perform interpolation and output result to 2D ASCII file."

#######
SCHEDULE ArbGrid_IncrementInterpCounter in ArbGrid_InterpGroup after Interpolate_to_dest_grid_main_function
{
  LANG: C
  OPTIONS: GLOBAL
} "Increment InterpCounter variable, or set to zero once loop is complete."
##################
```

    Writing interp_arbgrid_MO_ETK/schedule.ccl


<a id='readingoutputfile'></a>

# Step 4: Python Script for Reading the Output File \[Back to [top](#toc)\]
$$\label{readingoutputfile}$$

Here is a Python code for reading the output file generated by this thorn. It is based on a collection of Python scripts written by Bernard Kelly, available [here](https://bitbucket.org/zach_etienne/nrpy/src/master/mhd_diagnostics/). 

After generating the output file `interp_arbgrid_MO_ETK.dat` using the Einstein Toolkit thorn above, this script will read in all the data. Processing can then be done by straightforward modification of this script. Save the script as "Interp_Arb_ReadIn.py", and run it using the command

**`python Interp_Arb_ReadIn.py interp_arbgrid_MO_ETK.dat 58 outfile`**

Currently the last parameter "outfile" is required but not used.

```python
"""
interp_arbgrid_MO_ETK.dat File Reader. Compatible with Python 2.7+ and 3.6+ at least.

Zachariah B. Etienne

Based on Python scripts written by Bernard Kelly:
https://bitbucket.org/zach_etienne/nrpy/src/master/mhd_diagnostics/

Find the latest version of this reader at the bottom of this Jupyter notebook:
https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-ETK_thorn-Interpolation_to_Arbitrary_Grids_multi_order.ipynb

Usage instructions:

From the command-line, run via:
python Interp_Arb_ReadIn.py interp_arbgrid_MO_ETK.dat [number of gridfunctions (58 or so)] [outfile]

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
        #   generally be gibberish, we no inter append 
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

    # Then the number of interpolation points (stored as an int)
    num_interp_points = struct.unpack('i',filehandle.read(4))[0]

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

    return gf_name,order,num_interp_points,cctk_iteration,cctk_time

# Now open the file and read all the data
with open(datafile,"rb") as f:
    # Main loop over all gridfunctions
    for i in range(number_of_gridfunctions):
        # Data are output in chunks, one gridfunction at a time, with metadata
        #    for each gridfunction stored at the top of each chunk
        # First read in the metadata:
        gf_name, order, num_interp_points, cctk_iteration, cctk_time = read_header(f)
        print("\nReading gridfunction "+gf_name+", stored at interp order = "+str(order))
        data_chunk_size = num_interp_points*8 # 8 bytes per double-precision number
        # Next read in the full gridfunction data
        bytechunk = f.read(data_chunk_size)
        # Process the data using NumPy's frombuffer() function:
        #   https://docs.scipy.org/doc/numpy/reference/generated/numpy.frombuffer.html
        buffer_res = np.frombuffer(bytechunk)
        # Reshape the data into a 3D NumPy array:
        #   https://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
        #this_data = buffer_res.reshape(N0,N1,N2)

        # Sanity check: Output data at all points:
        with open("output-gf"+str(i)+".txt","w") as file:
            for ii in range(num_interp_points):
                file.write(str(ii) + "\t" + str(buffer_res[ii])+"\n")
```

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-Interpolation_to_Arbitrary_Grids_multi_order.pdf](Tutorial-ETK_thorn-Interpolation_to_Arbitrary_Grids_multi_order.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-Interpolation_to_Arbitrary_Grids_multi_order")
```

    Created Tutorial-ETK_thorn-
        Interpolation_to_Arbitrary_Grids_multi_order.tex, and compiled LaTeX
        file to PDF file Tutorial-ETK_thorn-
        Interpolation_to_Arbitrary_Grids_multi_order.pdf

