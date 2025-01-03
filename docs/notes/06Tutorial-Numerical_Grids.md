<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Numerical Grids in NRPy+

## Author: Zach Etienne

### NRPy+ Source Code for this module: [grid.py](../edit/grid.py)

Solving partial differential equations on the computer with finite difference techniques requires that we sample our solutions to these equations on numerical grids. We call all sampled functions that are stored on our numerical grid *gridfunctions*. NRPy+'s grid module adds the capability of registering gridfunctions in NRPy, setting basic parameters of a numerical grid, and providing functions to other modules regarding reading and writing of gridfunctions to memory in the C code.

Parameters in this module include:
* **grid::DIM** -- the dimension of the grid (e.g., a 3D numerical grid will have grid::DIM=3).
* **grid::Nx\[DIM\]** -- an integer array yielding the size of the grid in each direction. E.g., in Cartesian coordinates Nx\[0\] will be set to the number of grid points in the x direction. *This is a parameter that is set at C runtime, not in NRPy+; NRPy+ simply generates the declaration of this parameter in the C code.*
* **grid::MemAllocStyle** -- how the gridfunction is allocated in memory in the C code. This is used when generating the C code to ensure that gridfunctions are read as sequentially in memory as possible, to avoid [cache misses](https://en.wikipedia.org/wiki/CPU_cache#CACHE-MISS). There are currently two MemAllocStyles supported:
    * If the following loop accesses the grid function in a contiguous manner in memory, we set "grid::MemAllocStyle=012":
        * for(int i0=0;i0<Nx\[0\];i0++) for(int i1=0;i1<Nx\[1\];i1++) for(int i2=0;i2<Nx\[2\];i2++) {
    * Alternatively, the "grid::MemAllocStyle=210" is the reverse:
        * for(int i2=0;i2<Nx\[2\];i2++) for(int i1=0;i1<Nx\[1\];i1++) for(int i0=0;i0<Nx\[0\];i0++) {
* **grid::GridFuncMemAccess** -- specifies how gridfunction data is accessed from memory. For example,
    * In the Einstein Toolkit (grid::GridFuncMemAccess="ETK"), the datum from gridfunction dummy at point (i0,j0,k0) is accessed via "dummyGF\[CCTK_GFINDEX3D(cctkGH,i0,j0,k0)\]".
    * In SENR (grid::GridFuncMemAccess="SENRlike"), the datum from gridfunction dummy in gridfunction array gfarry at point (i0,j0,k0) is accessed via "gfarray\[IDX4D(DUMMYGF,i0,j0,k0)\]".
    * *Special note*: NRPy+ is code agnostic; additional types can be easily be added by modifying the function gfaccess() in the NRPy+ grid module. It should only require about 5 lines of code.

Functions in this module include:
* **gfaccess(gfarrayname = "",varname = "",ijklstring = "")**: given a gridfunction array name, a variable name, and a string indicating the coordinates of the point in memory, return the string to access the gridfunction in memory at this data point. See grid::GridFuncMemAccess parameter above for a complete description.
* **register_gridfunctions(gf_type,gf_names)**: returns either a single gridfunction or list of gridfunctions as SymPy variables.
    * gf_type can be set to either "EVOL" or "AUX".
        * Setting to "EVOL" denotes a grid function that will need to be stepped forward in time within the C code's timestepping algorithm.
        * Setting to "AUX" denotes any other grid function.
    * gf_names can be a single gridfunction or a set of gridfunctions, all with the same gf_type.
* **variable_type(var)**: first searches the list of registered gridfunctions and parameters for gridfunctions or C parameters; then outputs "gridfunction" or "Cparameter" if one or the other is found, respectively. Otherwise it will output "other"
* **gridfunction_lists()**:
    + The function returns three lists, corresponding to the names (strings) of the evolved, auxiliary, and "auxevol" gridfunction names respectively. Each is interpreted by NRPy+ as follows...
        + "Evolved" gridfunctions in NRPy+ will be automatically handled with the [MoLtimestepping](Tutorial-Method_of_Lines-C_Code_Generation.ipynb) module, including memory allocation (`RK_Allocate_Memory.h`), MoL updates (`RK_MoL.h`), and memory deallocation (`RK_Free_Memory.h`).
        + "Auxevol" gridfunctions in NRPy+ must be manually allocated and deallocated by the user. These gridfunctions provide additional storage for gridfunctions needed within the MoL step.
        + "Auxiliary" gridfunctions in NRPy+, soon to be ***deprecated*** and replaced by the "diagnostic" gridfunction type, are used for diagnostic purposes only, at the end of each timestep. As only a single MoL timelevel is needed to launch the next timestep, "auxiliary" gridfunctions use memory allocated for "evolved" gridfunctions at a different timelevel and can be reached via the `*diagnostic_gfs` pointer. "auxiliary"/"diagnostic" gridfunctions therefore cannot be used with "Euler" timestepping, and one cannot use more "auxiliary"/"diagnostic" gridfunctions than the number of "evolved" gridfunctions.
    + The aliases `#define`'d in "outdir/gridfunction_defines.h" are meant to be human-friendly, so that accessing each gridfunction in C code can be done by its human-friendly alias (e.g., test_gfs\[IDX4(VVGF,i0,i1,i2\] instead of the less-friendly test_gfs\[IDX4(1,i0,i1,i2\]). 
    * Example: if we register with NRPy+ only two gridfunctions uu and vv, which are evolved quantities (i.e., represent the solution of the PDEs we are solving, and are registered with gftype == "EVOL"), then the first returned list (all gridfunctions  registered as EVOL) will be \["uu","vv"\], and the second (all gridfunctions registered as AUX) will be the empty list: \[\]. Also, this function will create a file with the following content:

```C
/* This file is automatically generated by NRPy+. Do not edit. */
/* EVOLVED VARIABLES: */
#define NUM_EVOL_GFS 2
#define UUGF 0
#define VVGF 1

/* AUXILIARY VARIABLES: */
#define NUM_AUX_GFS 0

/* AUXEVOL VARIABLES: */
#define NUM_AUXEVOL_GFS 0
```

Let's now register a gridfunction called "phi" in the "in_gfs" gridfunction array, then print C code that specifies how to access the gridfunction from the single point in memory (i0,i1,i2). Next, we will demonstrate that phi is a normal SymPy variable, and finally, we will confirm that its type is "gridfunction" using variable_type():


```python
import grid as gri
import NRPy_param_funcs as par
import cmdline_helper as cmd

# Register gridfunction phi, as type "AUX".
# WARNING: register_gridfunctions can only be run once on a given gridfunction;
#          you'll need to reset the Jupyter kernel before trying again:
phi = gri.register_gridfunctions("AUX","phi")

print("Here's how to access the gridfunction phi, in grid array in_gfs, at point (i0,i1,i2):")
print("SENR-like memory access: "+gri.gfaccess("in_gfs","phi","i0,i1,i2"))

par.set_paramsvals_value("grid::GridFuncMemAccess = ETK")
print("ETK memory access: "+gri.gfaccess("in_gfs","phi","i0,i1,i2"))

# Note that phi can now be used as a usual SymPy variable:
print("\n\nTo demonstrate that phi is just a regular SymPy variable, we will now print its square:")
print(phi**2)

print("\n\n\"phi\" is of type \""+ gri.variable_type(phi)+"\"")

evolved_variables_list,auxiliary_variables_list, auxevol_variables_list = \
        gri.gridfunction_lists()[0:3]
print("\n\n")
print("Here is the list of registered evolved variables: ",evolved_variables_list)
print("... and here is the list of registered auxiliary variables: ",auxiliary_variables_list)
print("... and here is the list of registered auxevol variables: ",auxevol_variables_list)
print("\n\n ... and here is the output of gridfunction_defines():")
print(gri.gridfunction_defines())
```

    Here's how to access the gridfunction phi, in grid array in_gfs, at point (i0,i1,i2):
    SENR-like memory access: aux_gfs[IDX4S(PHIGF, i0,i1,i2)]
    ETK memory access: phiGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)]
    
    
    To demonstrate that phi is just a regular SymPy variable, we will now print its square:
    phi**2
    
    
    "phi" is of type "gridfunction"
    
    
    
    Here is the list of registered evolved variables:  []
    ... and here is the list of registered auxiliary variables:  ['phi']
    ... and here is the list of registered auxevol variables:  []
    
    
     ... and here is the output of gridfunction_defines():
    // EVOLVED VARIABLES:
    #define NUM_EVOL_GFS 0
    
    
    // AUXILIARY VARIABLES:
    #define NUM_AUX_GFS 1
    #define PHIGF	0
    
    
    // AUXEVOL VARIABLES:
    #define NUM_AUXEVOL_GFS 0
    
    
    // EXTERNAL VARIABLES:
    #define NUM_EXTERNAL_GFS 0
    

