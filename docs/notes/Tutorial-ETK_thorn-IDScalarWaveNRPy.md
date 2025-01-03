<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# IDScalarWaveNRPy: An Einstein Toolkit Initial Data Thorn for the Scalar Wave Equation

## Author: Terrence Pierre Jacques & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## The notebook outlines the construction of an Einstein Toolkit (ETK) thorn for setting up initial data for the scalar wave initial value problem. The ETK module employs NRPy+ to transform SymPy expressions for initial data into a C-code kernel, which is then integrated with the ETK infrastructure. The tutorial demonstrates the conversion of SymPy expressions into C-code, interfacing with the ETK, constructing C-code calling functions, defining module interactions with ETK, and the addition of C-code to the ETK compilation list.

[comment]: <> (Notebook Status and Validation Notes: TODO)

### NRPy+ Source Code for this module: [ScalarWave/InitialData.py](../edit/ScalarWave/InitialData.py) [\[**tutorial**\]](Tutorial-ScalarWave.ipynb) Contructs the SymPy expressions for spherical gaussian and plane-wave initial data

## Introduction:
In this part of the tutorial, we will construct an Einstein Toolkit (ETK) thorn (module) that will set up *initial data* for the scalar wave initial value problem. In a [previous tutorial notebook](Tutorial-ScalarWave.ipynb), we used NRPy+ to construct the SymPy expressions for either spherical Gaussian or plane-wave initial data. This thorn is largely based on and should function similarly to the $\text{IDScalarWaveC}$ thorn included in the Einstein Toolkit (ETK) $\text{CactusWave}$ arrangement.

We will construct this thorn in two steps.

1. Call on NRPy+ to convert the SymPy expressions for the initial data into one C-code kernel.
1. Write the C code and linkages to the Einstein Toolkit infrastructure (i.e., the .ccl files) to complete this Einstein Toolkit module.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$ 

This notebook is organized as follows

1. [Step 1](#initializenrpy): Call on NRPy+ to convert the SymPy expression for the scalar wave initial data into a C-code kernel
1. [Step 2](#einstein): Interfacing with the Einstein Toolkit
    1. [Step 2.a](#einstein_c): Constructing the Einstein Toolkit C-code calling functions that include the C code kernels
    1. [Step 2.b](#einstein_ccl): CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure
    1. [Step 2.c](#einstein_list): Add the C code to the Einstein Toolkit compilation list
1. [Step 3](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules \[Back to [top](#toc)\]

$$\label{initializenrpy}$$


```python
# Step 1: Import needed core NRPy+ modules
from outputC import lhrh         # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd  # NRPy+: Multi-platform Python command-line interface
import os, sys                   # Standard Python modules for multiplatform OS-level functions

# Step 1a: Create directories for the thorn if they don't exist.
# Create directory for WaveToyNRPy thorn & subdirectories in case they don't exist.
outrootdir = "IDScalarWaveNRPy/"
cmd.mkdir(os.path.join(outrootdir))
outdir = os.path.join(outrootdir,"src") # Main C code output directory
cmd.mkdir(outdir)

# Step 1b: This is an Einstein Toolkit (ETK) thorn. Here we
#          tell NRPy+ that gridfunction memory access will
#          therefore be in the "ETK" style.
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")
```

<a id='einstein'></a>

# Step 2: Interfacing with the Einstein Toolkit \[Back to [top](#toc)\]
$$\label{einstein}$$


<a id='einstein_c'></a>

## Step 2.a: Constructing the Einstein Toolkit C-code calling functions that include the C code kernels \[Back to [top](#toc)\]
$$\label{einstein_c}$$

Using SymPy, we construct the exact expressions for all scalar wave initial data currently supported in NRPy, documented in [Tutorial-ScalarWave.ipynb](Tutorial-ScalarWave.ipynb). We write the generated C codes into different C files, corresponding to the type of initial data they may want to choose at run time. Note that the code below can be easily extensible to include other types of initial data.


```python
# Step 1c: Call the InitialData() function from within the
#          ScalarWave/InitialData.py module.
import ScalarWave.InitialData as swid

_time = par.Cparameters("REAL", __name__, "time", 0.0)

# Step 1e: Call the InitialData() function to set up initial data.
#         Options include:
#    "PlaneWave": monochromatic (single frequency/wavelength) plane wave
#    "SphericalGaussian": spherically symmetric Gaussian, with default stdev=3
ID_options = ["PlaneWave", "SphericalGaussian"]
for ID in ID_options:
    gri.glb_gridfcs_list = []

# Within the ETK, the 3D gridfunctions x, y, and z store the
# Cartesian grid coordinates. Setting the gri.xx[] arrays
# to point to these gridfunctions forces NRPy+ to treat
# the Cartesian coordinate gridfunctions properly --
# reading them from memory as needed.
    x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])
    rfm.xx[0] = x
    rfm.xx[1] = y
    rfm.xx[2] = z
    swid.InitialData(WaveType=ID,
                     default_sigma=0.25,
                     default_k0=1.0,
                     default_k1=0.,
                     default_k2=0.)

    # Step 1f: Register uu and vv gridfunctions so they can be written to by NRPy.
    uu,vv = gri.register_gridfunctions("EVOL",["uu","vv"])

    # Step 1g: Set the uu and vv gridfunctions to the uu_ID & vv_ID variables
    #         defined by InitialData_PlaneWave().
    uu = swid.uu_ID
    vv = swid.vv_ID

    # Step 1h: Create the C code output kernel.
    ScalarWave_ID_SymbExpressions = [\
                            lhrh(lhs=gri.gfaccess("out_gfs","uu"),rhs=uu),\
                            lhrh(lhs=gri.gfaccess("out_gfs","vv"),rhs=vv),]

    ScalarWave_ID_CcodeKernel = fin.FD_outputC("returnstring",ScalarWave_ID_SymbExpressions)

    ScalarWave_ID_looped = lp.loop(["i2","i1","i0"],["0","0","0"],["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],\
                                   ["1","1","1"],["#pragma omp parallel for","",""],"",\
                                   ScalarWave_ID_CcodeKernel.replace("time","cctk_time"))

    #  Write the C code kernel to file.
    with open(os.path.join(outdir,"ScalarWave_"+ID+"ID.h"), "w") as file:
        file.write(str(ScalarWave_ID_looped))

```

<a id='einstein_ccl'></a>

## Step 2. b: CCL files - Define how this module interacts and interfaces with the larger Einstein Toolkit infrastructure \[Back to [top](#toc)\]
$$\label{einstein_ccl}$$

Writing a module ("thorn") within the Einstein Toolkit requires that three "ccl" files be constructed, all in the root directory of the thorn.

1. `interface.ccl`: defines the gridfunction groups needed, and provides keywords denoting what this thorn provides and what it should inherit from other thorns. Specifically, this file governs the interaction between this thorn and others; more information can be found in the [official Einstein Toolkit documentation](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-179000D2.2). 
With "implements", we give our thorn its unique name. By "inheriting" other thorns, we tell the Toolkit that we will rely on variables that exist and are declared "public" within those functions.


```python
evol_gfs_list    = []
for i in range(len(gri.glb_gridfcs_list)):
    if gri.glb_gridfcs_list[i].gftype == "EVOL":
        evol_gfs_list.append(   gri.glb_gridfcs_list[i].name+"GF")

# NRPy+'s finite-difference code generator assumes gridfunctions
#    are alphabetized; not sorting may result in unnecessary
#    cache misses.
evol_gfs_list.sort()

with open(os.path.join(outrootdir,"interface.ccl"), "w") as file:
    file.write("""
# With "implements", we give our thorn its unique name.
implements: IDScalarWaveNRPy

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits: WaveToyNRPy grid
""")
```

2. `param.ccl`: specifies free parameters within the thorn, enabling them to be set at runtime. It is required to provide allowed ranges and default values for each parameter. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-184000D2.3).


```python
def keep_param__return_type(paramtuple):
    keep_param = True # We'll not set some parameters in param.ccl;
                      #   e.g., those that should be #define'd like M_PI.
    typestring = ""
    # Separate thorns within the ETK take care of grid/coordinate parameters;
    #   thus we ignore NRPy+ grid/coordinate parameters:
    if paramtuple.module == "grid" or paramtuple.module == "reference_metric" or paramtuple.parname == "wavespeed":
        keep_param = False

    partype = paramtuple.type
    if partype == "bool":
        typestring += "BOOLEAN "
    elif partype == "REAL":
        if paramtuple.defaultval != 1e300: # 1e300 is a magic value indicating that the C parameter should be mutable
            typestring += "CCTK_REAL "
        else:
            keep_param = False
    elif partype == "int":
        typestring += "CCTK_INT "
    elif partype == "#define":
        keep_param = False
    elif partype == "char":
        # FIXME: char/string parameter types should in principle be supported
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    else:
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    return keep_param, typestring

paramccl_str="""
# This param.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.

shares: grid

USES KEYWORD type

shares: WaveToyNRPy

USES REAL wavespeed

restricted:
CCTK_KEYWORD initial_data "Type of initial data"
{"""

for ID in ID_options:
    paramccl_str +='''
  "'''+ID+'''"      :: "'''+ID+'"'

paramccl_str +='''
} "'''+ID+'''"

'''
paramccl_str +="""
restricted:

"""

for i in range(len(par.glb_Cparams_list)):
    # keep_param is a boolean indicating whether we should accept or reject
    #    the parameter. singleparstring will contain the string indicating
    #    the variable type.
    keep_param, singleparstring = keep_param__return_type(par.glb_Cparams_list[i])

    if keep_param:
        parname = par.glb_Cparams_list[i].parname
        partype = par.glb_Cparams_list[i].type
        singleparstring += parname + " \""+ parname +" (see NRPy+ for parameter definition)\"\n"
        singleparstring += "{\n"
        if partype != "bool":
            singleparstring += " *:* :: \"All values accepted. NRPy+ does not restrict the allowed ranges of parameters yet.\"\n"
        singleparstring += "} "+str(par.glb_Cparams_list[i].defaultval)+"\n\n"

        paramccl_str += singleparstring
with open(os.path.join(outrootdir,"param.ccl"), "w") as file:
    file.write(paramccl_str)
```

3. `schedule.ccl`: allocates storage for gridfunctions, defines how the thorn's functions should be scheduled in a broader simulation, and specifies the regions of memory written to or read from gridfunctions. More information on this file's syntax can be found in the [official Einstein Toolkit documentation](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-187000D2.4). 

We specify here the standardized ETK "scheduling bins" in which we want each of our thorn's functions to run.


```python
with open(os.path.join(outrootdir,"schedule.ccl"), "w") as file:
    file.write("""
# This schedule.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.

if (CCTK_EQUALS (initial_data, "PlaneWave"))
{
    schedule IDScalarWaveNRPy_param_check at CCTK_PARAMCHECK
    {
      LANG: C
      OPTIONS: global
    } "Check sanity of parameters"
}

schedule IDScalarWaveNRPy_InitialData at CCTK_INITIAL as WaveToy_InitialData
{
 STORAGE: WaveToyNRPy::scalar_fields[3]
  LANG:          C
} "Initial data for 3D wave equation"
""")
```

<a id='einstein_list'></a>

## Step 2.c: Add the C code to the Einstein Toolkit compilation list \[Back to [top](#toc)\]
$$\label{einstein_list}$$

We will also need `make.code.defn`, which indicates the list of files that need to be compiled. This thorn only has the one C file to compile.


```python
make_code_defn_list = []
def append_to_make_code_defn_list(filename):
    if filename not in make_code_defn_list:
        make_code_defn_list.append(filename)
    return os.path.join(outdir,filename)
```


```python
with open(append_to_make_code_defn_list("InitialData.c"),"w") as file:
    file.write("""

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void IDScalarWaveNRPy_param_check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (kk0 == 0 && kk1 == 0 && kk2 == 0) {
     CCTK_WARN(0,"kk0==kk1==kk2==0: Zero wave vector cannot be normalized. Set one of the kk's to be != 0.");
  }
}

void IDScalarWaveNRPy_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const CCTK_REAL *xGF = x;
  const CCTK_REAL *yGF = y;
  const CCTK_REAL *zGF = z;

  if (CCTK_EQUALS (initial_data, "PlaneWave")) {
      #include "ScalarWave_PlaneWaveID.h"
  } else if (CCTK_EQUALS (initial_data, "SphericalGaussian")) {
      #include "ScalarWave_SphericalGaussianID.h"
  }
}
""")
```


```python
with open(os.path.join(outdir,"make.code.defn"), "w") as file:
    file.write("""
# Main make.code.defn file for thorn WaveToyNRPy

# Source files in this directory
SRCS =""")
    filestring = ""
    for i in range(len(make_code_defn_list)):
        filestring += "      "+make_code_defn_list[i]
        if i != len(make_code_defn_list)-1:
            filestring += " \\\n"
        else:
            filestring += "\n"
    file.write(filestring)
```

<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ETK_thorn-IDScalarWaveNRPy.pdf](Tutorial-ETK_thorn-IDScalarWaveNRPy.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ETK_thorn-IDScalarWaveNRPy")
```

    [WARNING] Duplicate link reference '[comment]' at line 17 column 1
    Created Tutorial-ETK_thorn-IDScalarWaveNRPy.tex, and compiled LaTeX file to
        PDF file Tutorial-ETK_thorn-IDScalarWaveNRPy.pdf

