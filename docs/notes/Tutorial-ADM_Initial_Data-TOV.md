<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# [TOV](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) Initial Data

## Authors: Phil Chang, Zach Etienne, & Leo Werneck
### Formatting improvements courtesy Brandon Clark

## This notebook establishes initial data for a [TOV](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) star in spherical, isotropic coordinates. The notebook presents a module for setting up ADM metric quantities and the stress-energy tensor $T^{\mu\nu}$ from interpolated TOV data and demonstrates reliability through validation against the Python `TOV.TOV_Ccodegen_library` module. Detailed steps for reading and interpolating TOV data files, while avoiding the Gibbs phenomenon at the stellar surface, are outlined.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This module has been validated to exhibit convergence to zero of the Hamiltonian constraint violation at the expected order to the exact solution (see [start-to-finish TOV module](Tutorial-Start_to_Finish-BSSNCurvilinear-TOV_initial_data.ipynb) for a full test). Note that convergence at the surface of the star is lower order due to the sharp drop to zero in $T^{\mu\nu}$.

### NRPy+ Source Code for this module: [TOV/TOV_Solver.py](../edit/TOV/TOV_Solver.py)

[comment]: <> (Introduction: TODO)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#tov): The TOV Equations
1. [Step 3](#code_validation): Code Validation against `TOV.TOV_Solver` NRPy+ module
1. [Step 4](#c_code_generation): C-code routines (library) for reading TOV data files
    1. [Step 4.a](#id_persist): Declare `ID_persist`, the C `struct` data type that will store data read from file
    1. [Step 4.b](#read_data_file): `TOV_read_data_file_set_ID_persist()`: Read TOV data file and store data to the `ID_persist` struct
    1. [Step 4.c](#interp_data_file): `TOV_interpolate_1D()`: Interpolate TOV data to any desired distance from the center of the star
    1. [Step 4.d](#tov_id_function): `TOV_ID_function()`: Driver routine for setting up ADM metric quantities and $T^{\mu\nu}$ from interpolated TOV data
        1. [Step 4.d.i](#adm_ito_tov): ADM metric quantities, written in terms of the TOV solution
        1. [Step 4.d.ii](#tmunu_ito_tov): Stress-energy tensor $T^{\mu\nu}$, written in terms of the TOV solution
        1. [Step 4.d.iii](#tov_id_function_itself): The `TOV_ID_function()` itself
1. [Step 5](#ccodegen_library_validation): Validate C code generation against `TOV.TOV_Ccodegen_library` Python module 
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$


```python
# Step 1: Import needed Python/NRPy+ modules
import numpy as np                  # NumPy: A numerical methods module for Python
import scipy.integrate as si        # SciPy: Python module for mathematics, science, and engineering applications
import math, sys                    # Standard Python modules for math; multiplatform OS-level functions
import TOV.Polytropic_EOSs as ppeos # NRPy+: Piecewise polytrope equation of state support
from outputC import outputC,add_to_Cfunction_dict  # NRPy+: Core C code output module
import sympy as sp                  # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par      # NRPy+: Parameter interface
import indexedexp as ixp            # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm      # NRPy+: Reference metric support
```

<a id='tov'></a>

# Step 2: The TOV equations \[Back to [top](#toc)\]
$$\label{tov}$$

The [TOV line element](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) in terms of the *Schwarzschild coordinate* $r$ is written (in the $-+++$ form):
$$
ds^2 = - c^2 e^\nu dt^2 + \left(1 - \frac{2Gm}{rc^2}\right)^{-1} dr^2 + r^2 d\Omega^2,
$$
where $m(r)$ is the mass-energy enclosed at a given $r$, and is equal to the total star's mass outside the stellar radius $r=R$.

In terms of the *isotropic coordinate* $\bar{r}$ with $G=c=1$ (i.e., the coordinate system and units we'd prefer to use), the ($-+++$ form) line element is written:
$$
ds^2 = - e^{\nu} dt^2 + e^{4\phi} \left(d\bar{r}^2 + \bar{r}^2 d\Omega^2\right),
$$
where $\phi$ here is the *conformal factor*.

Setting components of the above line element equal to one another, we get (in $G=c=1$ units):

\begin{align}
r^2 &= e^{4\phi} \bar{r}^2 \implies e^{4\phi} = \frac{r^2}{\bar{r}^2} \\
\left(1 - \frac{2m}{r}\right)^{-1} dr^2 &= e^{4\phi} d\bar{r}^2 \\
\implies \frac{d\bar{r}(r)}{dr} &= \left(1 - \frac{2m}{r} \right)^{-1/2} \frac{\bar{r}(r)}{r}.
\end{align}

The TOV equations provide radial ODEs for the pressure and $\nu$ (from [the Wikipedia article on the TOV solution](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation)):

\begin{align}
\frac{dP}{dr} &= - \frac{1}{r} \left( \frac{\rho + P}{2} \right) \left(\frac{2 m}{r} + 8 \pi r^2 P\right) \left(1 - \frac{2 m}{r}\right)^{-1} \\
\frac{d \nu}{d r} &= \frac{1}{r}\left(1 - \frac{2 m}{r}\right)^{-1} \left(\frac{2 m}{r} + 8 \pi r^2 P\right) \\
\end{align}

Assuming a polytropic equation of state, which relates the pressure $P$ to the baryonic rest-mass density $\rho_B$,

$$
P(\rho_B) = K \rho_B^\Gamma,
$$
the specific internal energy will be given by
$$
\epsilon = \frac{P}{\rho_B (\Gamma - 1)},
$$

so the total mass-energy density $\rho$ is given by
$$
\rho = \rho_B (1 + \epsilon).
$$

Given this, the mass-energy $m(r)$ density is the solution to the ODE:
$$
\frac{dm(r)}{dr} = 4\pi r^2 \rho(r)
$$

Thus the full set of ODEs that need to be solved is given by

$$
\boxed{
\begin{matrix}
\frac{dP}{dr} &=& - \frac{1}{r} \left( \frac{\rho + P}{2} \right) \left(\frac{2 m}{r} + 8 \pi r^2 P\right) \left(1 - \frac{2 m}{r}\right)^{-1} \\
\frac{d \nu}{d r} &=& \frac{1}{r}\left(1 - \frac{2 m}{r}\right)^{-1} \left(\frac{2 m}{r} + 8 \pi r^2 P\right) \\
\frac{m(r)}{dr} &=& 4\pi r^2 \rho(r) \\
\frac{d\bar{r}(r)}{dr} &=& \left(1 - \frac{2m}{r} \right)^{-1/2} \frac{\bar{r}(r)}{r}
\end{matrix}
}\ .
$$

The following code solves these equations, and was largely written by Phil Chang.


```python
# Step 2: The TOV equations

## TOV SOLVER FOR SINGLE AND PIECEWISE POLYTROPES
## Authors: Phil Chang, Zachariah B. Etienne, Leo Werneck

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-Start_to_Finish-BSSNCurvilinear-TOV_initial_data.ipynb


# Inputs:
# * Output data file name
# * rho_baryon_central, the central density of the TOV star.
# * n, the polytropic equation of state index. n=1 models cold, degenerate neutron star matter.
# * K_Polytrope, the polytropic constant.
# * Verbose output toggle (default = True)

# Output: An initial data file (default file name = "outputTOVpolytrope.txt") that well
#         samples the (spherically symmetric) solution both inside and outside the star.
#         It is up to the initial data module to perform the 1D interpolation to generate
#         the solution at arbitrary radius. The file has the following columns:
# Column 1: Schwarzschild radius
# Column 2: rho(r), *total* mass-energy density (as opposed to baryonic rest-mass density)
# Column 3: P(r), Pressure
# Column 4: m(r), mass enclosed
# Column 5: e^{nu(r)}, g_{tt}(r)
# Column 6: e^{4 phi(r)}, conformal factor g_{rr}(r)
# Column 7: rbar(r), Isotropic radius

# rbar refers to the isotropic radius, and
# R_Schw refers to the Schwarzschild radius

def TOV_Solver(eos,
               outfile = "outputTOVpolytrope.txt",
               rho_baryon_central = 0.129285,
               verbose = True,
               return_M_and_RSchw = False,
               accuracy = "medium",
               integrator_type = "default",
               no_output_File = False,
               pressure_renormalization=1): # reset the pressure to stellar oscillations studies

    def TOV_rhs(r_Schw, y) :
    # In \tilde units
    #
        P    = y[0]
        m    = y[1]
        # nu   = y[2] # nu is not needed as input into TOV_rhs
        rbar = y[3]

        # Compute rho_b and eps_cold, to be used below
        # to compute rho_(total)
        rho_baryon, eps_cold = ppeos.Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold(eos,P)

#         with open("rhob_P_cold_and_eps_cold.dat","a+") as file:
#             file.write(str(r_Schw).format("%.15e")+"  "+str(rho_baryon).format("%.15e")+"  "+str(P).format("%.15e")+"  "+str(eps_cold).format("%.15e")+"\n")

        # Compute rho, the *total* mass-energy density:
        # .------------------------------.
        # | rho = (1 + eps)*rho_(baryon) |
        # .------------------------------.
        # with eps = eps_cold, for the initial data.
        rho = (1.0 + eps_cold)*rho_baryon

#         m = 4*math.pi/3. * rho*r_Schw**3
        if( r_Schw < 1e-4 or m <= 0.):
            # From https://github.com/natj/tov/blob/master/tov.py#L33:
            # dPdr = -cgs.G*(eden + P/cgs.c**2)*(m + 4.0*pi*r**3*P/cgs.c**2)
            # dPdr = dPdr/(r*(r - 2.0*cgs.G*m/cgs.c**2))
            dPdrSchw = -(rho + P)*(4.*math.pi/3.*r_Schw*rho + 4.*math.pi*r_Schw*P)/(1.-8.*math.pi*rho*r_Schw*r_Schw)
            drbardrSchw = 1./(1. - 8.*math.pi*rho*r_Schw*r_Schw)**0.5
        else:
            dPdrSchw = -(rho + P)*(m + 4.*math.pi*r_Schw**3*P)/(r_Schw*r_Schw*(1.-2.*m/r_Schw))
            drbardrSchw = 1./(1. - 2.*m/r_Schw)**0.5*rbar/r_Schw

        dmdrSchw  =  4.*math.pi*r_Schw*r_Schw*rho
        dnudrSchw = -2./(P + rho)*dPdrSchw
        return [dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw]

    def integrateStar( eos, P ):

        if accuracy == "medium":
            min_step_size = 1e-5
            max_step_size = 1e-2
            integrator    = 'dop853'
        elif accuracy == "low":
            min_step_size = 1e-3
            max_step_size = 1e-1
            integrator    = 'dopri5'
        elif accuracy == "verylow":
            min_step_size = 1e-1
            max_step_size = 5e-1
            integrator    = 'dopri5'
        elif accuracy == "high":
            min_step_size = 1e-5
            max_step_size = 1e-5
            integrator    = 'dop853'
        elif accuracy == "veryhigh":
            min_step_size = 1e-7
            max_step_size = 1e-6
            integrator    = 'dop853'
        else:
            print("Unknown accuracy option: "+str(accuracy))

        if integrator_type == "default":
            pass
        else:
            integrator = integrator_type

        integrator = si.ode(TOV_rhs).set_integrator(integrator)#,rtol=1e-4,atol=1e-4)
#         integrator = si.ode(TOV_rhs).set_integrator('dopri5',rtol=1e-4)
        y0 = [P, 0., 0., 0.]
        r_Schw = 0.
        integrator.set_initial_value(y0,r_Schw)
        dr_Schw = min_step_size
        P = y0[0]

        PArr      = []
        r_SchwArr = []
        mArr      = []
        nuArr     = []
        rbarArr   = []

        while integrator.successful() and P > 1e-19*y0[0] :
            P, m, nu, rbar = integrator.integrate(r_Schw + dr_Schw)
            # Update the value of r_Schw to the latest integrated value
            r_Schw += dr_Schw

            dPdrSchw, dmdrSchw, _dnudrSchw, _drbardrSchw = TOV_rhs( r_Schw, [P,m,nu,rbar])
            dr_Schw = 0.1*min(abs(P/dPdrSchw), abs(m/dmdrSchw))
            dr_Schw = min(dr_Schw, max_step_size)
            PArr.append(P)
            r_SchwArr.append(r_Schw)
            mArr.append(m)
            nuArr.append(nu)
            rbarArr.append(rbar)

        M = mArr[-1]
        R_Schw = r_SchwArr[-1]

        if no_output_File == True:
            return R_Schw, M

        # Apply integration constant to ensure rbar is continuous across TOV surface
        for ii in range(len(rbarArr)):
            rbarArr[ii] *= 0.5*(np.sqrt(R_Schw*(R_Schw - 2.0*M)) + R_Schw - M) / rbarArr[-1]

        nuArr_np = np.array(nuArr)
        # Rescale solution to nu so that it satisfies BC: exp(nu(R))=exp(nutilde-nu(r=R)) * (1 - 2m(R)/R)
        #   Thus, nu(R) = (nutilde - nu(r=R)) + log(1 - 2*m(R)/R)
        nuArr_np = nuArr_np - nuArr_np[-1] + math.log(1.-2.*mArr[-1]/r_SchwArr[-1])

        r_SchwArrExtend_np = 10.**(np.arange(0.01,5.0,0.01))*r_SchwArr[-1]

        r_SchwArr.extend(r_SchwArrExtend_np)
        mArr.extend(r_SchwArrExtend_np*0. + M)
        PArr.extend(r_SchwArrExtend_np*0.)
        exp2phiArr_np = np.append( np.exp(nuArr_np), 1. - 2.*M/r_SchwArrExtend_np)
        nuArr.extend(np.log(1. - 2.*M/r_SchwArrExtend_np))
        rbarArr.extend( 0.5*(np.sqrt(r_SchwArrExtend_np**2 - 2.*M*r_SchwArrExtend_np) + r_SchwArrExtend_np - M) )

        #phiArr_np = np.append( np.exp(nuArr_np), 1. - 2.*M/r_SchwArrExtend_np)

        # Appending to a Python array does what one would reasonably expect.
        #   Appending to a numpy array allocates space for a new array with size+1,
        #   then copies the data over... over and over... super inefficient.
        r_SchwArr_np     = np.array(r_SchwArr)
        PArr_np          = np.array(PArr)
        rho_baryonArr_np = np.array(PArr) # This is just to initialize the array
        for j in range(len(PArr_np)):
            # Compute rho_b from P
            rho_baryonArr_np[j] = ppeos.Polytrope_EOS__compute_rhob_from_P_cold(eos,PArr_np[j])

        mArr_np               = np.array(mArr)
        rbarArr_np            = np.array(rbarArr)
        confFactor_exp4phi_np = (r_SchwArr_np/rbarArr_np)**2

        # Compute the *total* mass-energy density (as opposed to the *baryonic* mass density)
        rhoArr_np = []
        for i in range(len(PArr)):
            rho_baryon, eps_cold = ppeos.Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold(eos,PArr[i])
            rho = (1.0 + eps_cold ) * rho_baryon
            rhoArr_np.append(rho)

        if verbose:
            print(len(r_SchwArr_np),len(rhoArr_np),len(rho_baryonArr_np),len(PArr_np),len(mArr_np),len(exp2phiArr_np))

        PArr_np *= pressure_renormalization # set for pressure renormalization studies

        # Special thanks to Leonardo Werneck for pointing out this issue with zip()
        if sys.version_info[0] < 3:
            np.savetxt(outfile, zip(r_SchwArr_np,rhoArr_np,rho_baryonArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np),
                       fmt="%.15e")
        else:
            np.savetxt(outfile, list(zip(r_SchwArr_np,rhoArr_np,rho_baryonArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np)),
                       fmt="%.15e")

        return R_Schw, M

    # Set initial condition from rho_baryon_central
    P_initial_condition = ppeos.Polytrope_EOS__compute_P_cold_from_rhob(eos, rho_baryon_central)

    # Integrate the initial condition
    R_Schw_TOV, M_TOV = integrateStar(eos, P_initial_condition)
    if verbose:
        print("Just generated a TOV star with R_Schw = %.15e , M = %.15e , M/R_Schw = %.15e ." %(R_Schw_TOV,M_TOV,(M_TOV / R_Schw_TOV)))

    if return_M_and_RSchw:
        return M_TOV, R_Schw_TOV

############################
# Single polytrope example #
############################
# Set neos = 1 (single polytrope)
neos = 1

# Set rho_poly_tab (not needed for a single polytrope)
rho_poly_tab = []

# Set Gamma_poly_tab
Gamma_poly_tab = [2.0]

# Set K_poly_tab0
K_poly_tab0 = 1. # ZACH NOTES: CHANGED FROM 100.

# Set the eos quantities
eos = ppeos.set_up_EOS_parameters__complete_set_of_input_variables(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab0)

# Set initial condition (Pressure computed from central density)
rho_baryon_central  = 0.129285

M_TOV, R_Schw_TOV = TOV_Solver(eos,outfile="outputTOVpolytrope.txt",rho_baryon_central=0.129285,verbose = True,
                               return_M_and_RSchw=True,accuracy="medium",integrator_type="default",
                               pressure_renormalization=1.0)
```

    1256 1256 1256 1256 1256 1256
    Just generated a TOV star with R_Schw = 9.566044579232513e-01 , M = 1.405030336771405e-01 , M/R_Schw = 1.468768334847266e-01 .


<a id='code_validation'></a>

# Step 3: Code Validation \[Back to [top](#toc)\]
$$\label{code_validation}$$

Here, as a code validation check, we verify agreement in the SymPy expressions for these TOV initial data between

1. this tutorial and 
2. the NRPy+ [TOV.TOV_Solver](../edit/TOV/TOV_Solver.py) module.


```python
# Step 3: Code Validation against TOV.TOV_Solver module

import filecmp

import TOV.TOV_Solver as TOV
TOV.TOV_Solver(eos,
               outfile="outputTOVpolytrope-validation.txt",
               rho_baryon_central=0.129285,
               verbose = True,
               accuracy="medium",
               integrator_type="default",
               no_output_File = False)

if filecmp.cmp('outputTOVpolytrope.txt',
               'outputTOVpolytrope-validation.txt') == False:
    print("ERROR: TOV initial data test FAILED!")
    sys.exit(1)
else:
    print("TOV initial data test PASSED.")
```

    1256 1256 1256 1256 1256 1256
    Just generated a TOV star with
    * M        = 1.405030336771405e-01 ,
    * R_Schw   = 9.566044579232513e-01 ,
    * R_iso    = 8.100085557410308e-01 ,
    * M/R_Schw = 1.468768334847266e-01 
    
    TOV initial data test PASSED.


<a id='c_code_generation'></a>

# Step 4: C-code routines (library) for reading TOV data files \[Back to [top](#toc)\]
$$\label{c_code_generation}$$

<a id='id_persist'></a>

## Step 4.a: Declare `ID_persist`, the C `struct` data type that will store data read from file \[Back to [top](#toc)\]
$$\label{id_persist}$$

The TOV data file stores various physical quantities describing the star's pressure, density, and gravitational fields, *all as a function of distance from the center of the star*, `rbar`. `Rbar` is the radius of the star.

The `ID_persist` struct takes the following form

```c
  REAL Rbar;      // rbar corresponding to the outermost radius at which density is nonzero
  int Rbar_idx;   // Index (line of data file) corresponding to Rbar
  int interp_stencil_size;  // Lagrange polynomial interpolation stencil size
  int numlines_in_file;     // Number of radii stored in the file
  // The following arrays store stellar information at all numlines_in_file radii:
  REAL *restrict r_Schw_arr; // Stellar radial coordinate in units of Schwarzschild radius
  REAL *restrict rho_arr;    // Mass-energy density
  REAL *restrict rho_baryon_arr; // Baryonic mass density
  REAL *restrict P_arr;  // Pressure
  REAL *restrict M_arr;  // Integrated rest mass
  REAL *restrict expnu_arr;    // Metric quantity
  REAL *restrict exp4phi_arr;  // Metric quantity
  REAL *restrict rbar_arr;  // Rbar coordinate
```


```python
def ID_persist_str():
    return r"""
  REAL Rbar;      // rbar corresponding to the outermost radius at which density is nonzero
  int Rbar_idx;   // Index (line of data file) corresponding to Rbar
  int interp_stencil_size;  // Lagrange polynomial interpolation stencil size
  int numlines_in_file;     // Number of radii stored in the file
  // The following arrays store stellar information at all numlines_in_file radii:
  REAL *restrict r_Schw_arr; // Stellar radial coordinate in units of Schwarzschild radius
  REAL *restrict rho_arr;    // Mass-energy density
  REAL *restrict rho_baryon_arr; // Baryonic mass density
  REAL *restrict P_arr;  // Pressure
  REAL *restrict M_arr;  // Integrated rest mass
  REAL *restrict expnu_arr;    // Metric quantity
  REAL *restrict exp4phi_arr;  // Metric quantity
  REAL *restrict rbar_arr;  // Rbar coordinate
"""
```

<a id='read_data_file'></a>

## Step 4.b: `TOV_read_data_file_set_ID_persist()`: Read TOV data file and store data to the `ID_persist` struct \[Back to [top](#toc)\]
$$\label{read_data_file}$$

The following function also sets the `interp_stencil_size` parameter, indicating the total size of the Lagrange polynomial interpolation stencil. The default of 12 reflects a quite high interpolation order, corresponding to an 11th-order polynomial being fit through the data at various radii to estimate stellar quantities at a single desired (arbitrary) distance from the center of the star.

As hydrodynamic data (like density and pressure) at the stellar radius sharply drop to zero, the interpolation algorithm offsets the center of the stencil so that the interpolation never crosses the stellar surface. This avoids the [Gibbs phenomenon](https://en.wikipedia.org/wiki/Gibbs_phenomenon), ensuring super high fidelity of the output.


```python
def add_to_Cfunction_dict_TOV_read_data_file_set_ID_persist(interp_stencil_size=12):
    includes = ["NRPy_basic_defines.h"]
    desc = "Returns the number of lines in a TOV data file."
    c_type = "void"
    name = "TOV_read_data_file_set_ID_persist"
    params = "const char *input_filename, ID_persist_struct *ID_persist"
    body = r"""
  char filename[100];
  snprintf(filename, 100, input_filename);
  FILE *TOV_solution_datafile = fopen(filename, "r");
  if(TOV_solution_datafile == NULL) {
    fprintf(stderr,"ERROR: could not open TOV solution data file %s\n",filename);
    exit(1);
  }

  // First set the interpolation stencil size. Note that the interpolation is set up to
  //   avoid interpolations across the star's surface. Hence we can use a super high
  //   order interpolant.
  ID_persist->interp_stencil_size = """+str(interp_stencil_size)+r""";

  int numlines_in_file = 0;
  {
    char * line = NULL;

    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, TOV_solution_datafile)) != -1) {
      numlines_in_file++;
    }
    rewind(TOV_solution_datafile);

    free(line);
  }
  ID_persist->numlines_in_file = numlines_in_file;

  // Now that numlines_in_file is set, we can now allocate memory for all arrays.
  {
    ID_persist->r_Schw_arr     = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rho_arr        = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rho_baryon_arr = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->P_arr          = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->M_arr          = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->expnu_arr      = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->exp4phi_arr    = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rbar_arr       = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
  }

  {
    char * line = NULL;

    size_t len = 0;
    ssize_t read;

    int which_line = 0;
    while ((read = getline(&line, &len, TOV_solution_datafile)) != -1) {
      // Define the line delimiters (i.e., the stuff that goes between the data on a given
      //     line of data.  Here, we define both spaces " " and tabs "\t" as data delimiters.
      const char delimiters[] = " \t";

      // Now we define "token", a pointer to the first column of data
      char *token;

      // Each successive time we call strtok(NULL,blah), we read in a new column of data from
      //     the originally defined character array, as pointed to by token.

      token=strtok(line, delimiters); if(token==NULL) { fprintf(stderr, "Error reading %s\n", filename); exit(1); }
      ID_persist->r_Schw_arr[which_line]     = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rho_arr[which_line]        = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rho_baryon_arr[which_line] = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->P_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->M_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->expnu_arr[which_line]      = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->exp4phi_arr[which_line]    = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rbar_arr[which_line]       = strtod(token, NULL);

      which_line++;
    }
    free(line);

    fclose(TOV_solution_datafile);
  }

  {
    // Finally set Rbar and Rbar_idx
    ID_persist->Rbar     = -100.0;
    ID_persist->Rbar_idx = -100;
    for(int i=1;i<numlines_in_file;i++) {
      if(ID_persist->rho_arr[i-1] > 0  &&  ID_persist->rho_arr[i] == 0) {
        ID_persist->Rbar = ID_persist->rbar_arr[i-1];
        ID_persist->Rbar_idx = i-1;
      }
    }
    if(ID_persist->Rbar < 0) {
      fprintf(stderr,"Error: could not find rbar=Rbar (i.e., the surface of the star) from data file.\n");
      exit(1);
    }
  }
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
```

<a id='interp_data_file'></a>

## Step 4.c: `TOV_interpolate_1D()`: Interpolate TOV data to any desired distance from the center of the star \[Back to [top](#toc)\]
$$\label{interp_data_file}$$

The TOV solution stored in the data file samples only a very limited number of points inside and outside the star. However we'll need data at arbitrary radii between these points. `TOV_interpolate_1D()` provides this service.

`TOV_interpolate_1D()` applies 1-dimensional Lagrange polynomial interpolation with stencil size `ID_persist->interp_stencil_size` (default 12) to obtain TOV data between sampled points in the ordered data file (each point corresponds to a specific radius).

Once a desired output radius `rrbar` is chosen, `TOV_interpolate_1D()` calls the included bisection algorithm `bisection_idx_finder()` to find the closest sample point with radius `rbar`, to desired output point `rrbar`. If the stencil crosses the star's surface, it then adjusts the stencil to avoid the [Gibbs phenomenon](https://en.wikipedia.org/wiki/Gibbs_phenomenon) associated with Lagrange interpolating data with a kink.


```python
def add_to_Cfunction_dict_TOV_interpolate_1D():
    includes=["NRPy_basic_defines.h"]
    prefunc = r"""
// Find interpolation index using Bisection root-finding algorithm:
static inline int bisection_idx_finder(const REAL rrbar, const int numlines_in_file, const REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  REAL y1 = rrbar-rbar_arr[x1];
  REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 >= 0) {
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
    exit(1);
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
    if(y_midpoint*y1 < 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
      if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) return x1;
      // Otherwiser return x2:
      return x2;
    }
  }
  fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
  exit(1);
}
"""
    desc = """Read a TOV solution from data file and perform
1D interpolation of the solution to a desired radius.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
    c_type="void"
    name="TOV_interpolate_1D"
    params="""REAL rrbar,const ID_persist_struct *ID_persist,
                        REAL *restrict rho,REAL *restrict rho_baryon,REAL *restrict P,REAL *restrict M,REAL *restrict expnu,REAL *restrict exp4phi"""
    body = r"""
  // First unpack ID_persist struct (this should be pretty quick relative to the rest of the routine):
  const REAL Rbar               = ID_persist->Rbar;
  const int Rbar_idx            = ID_persist->Rbar_idx;
  const int interp_stencil_size = ID_persist->interp_stencil_size;
  const int numlines_in_file    = ID_persist->numlines_in_file;
  const REAL *restrict rbar_arr       = ID_persist->rbar_arr;
  const REAL *restrict r_Schw_arr     = ID_persist->r_Schw_arr;
  const REAL *restrict rho_arr        = ID_persist->rho_arr;
  const REAL *restrict rho_baryon_arr = ID_persist->rho_baryon_arr;
  const REAL *restrict P_arr          = ID_persist->P_arr;
  const REAL *restrict M_arr          = ID_persist->M_arr;
  const REAL *restrict expnu_arr      = ID_persist->expnu_arr;
  const REAL *restrict exp4phi_arr    = ID_persist->exp4phi_arr;

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rrbar,numlines_in_file,rbar_arr);


  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

  // -= Do not allow the interpolation stencil to cross the star's surface =-
  // max index is when idxmin + (interp_stencil_size-1) = Rbar_idx
  //  -> idxmin at most can be Rbar_idx - interp_stencil_size + 1
  if(rrbar < Rbar) {
    idxmin = MIN(idxmin,Rbar_idx - interp_stencil_size + 1);
  } else {
    idxmin = MAX(idxmin,Rbar_idx+1);
    idxmin = MIN(idxmin,numlines_in_file - interp_stencil_size + 1);
  }
  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  REAL rbar_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rbar_sample[i-idxmin] = rbar_arr[i];
  }
  REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *rho = 0.0;
  *rho_baryon = 0.0;
  *P = 0.0;
  *M = 0.0;
  *expnu = 0.0;
  *exp4phi = 0.0;

  REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw      += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho        += l_i_of_r[i-idxmin] * rho_arr[i];
    *rho_baryon += l_i_of_r[i-idxmin] * rho_baryon_arr[i];
    *P          += l_i_of_r[i-idxmin] * P_arr[i];
    *M          += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu      += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi    += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

  if(rrbar > Rbar) {
    *rho        = 0;
    *rho_baryon = 0;
    *P          = 0;
    *M          = M_arr[Rbar_idx+1];
    *expnu      = 1. - 2.*(*M) / r_Schw;
    *exp4phi    = pow(r_Schw / rrbar,2.0);
  }
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
```

<a id='tov_id_function'></a>

## Step 4.d: `TOV_ID_function()`: Driver routine for setting up ADM metric quantities and $T^{\mu\nu}$ from interpolated TOV data \[Back to [top](#toc)\]
$$\label{tov_id_function}$$

`TOV_ID_function()` is a function with standard input and output parameters, providing the necessary linkage to the NRPy+ BSSN initial data infrastructure. This infrastructure takes as input ADM metric data and $T^{\mu\nu}$ in the spherical or Cartesian basis. TOV data, of course, are provided in the spherical basis.

In other words, `TOV_ID_function()` takes as input the desired Cartesian point $(x,y,z)$ at which we desire initial data, and calls the `TOV_interpolate_1D()` routine to get raw data at that point. `TOV_ID_function()` converts the interpolated quantities to ADM metric and $T^{\mu\nu}$, stored in the standardized `initial_data` struct.

Conversions between data within the TOV solution output file and ADM/$T^{\mu\nu}$ are as follows:

<a id='adm_ito_tov'></a>

### Step 4.d.i: ADM metric quantities, written in terms of the TOV solution \[Back to [top](#toc)\]
$$\label{adm_ito_tov}$$

The [TOV line element](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) in *Schwarzschild coordinates* is written (in the $-+++$ form):
$$
ds^2 = - c^2 e^\nu dt^2 + \left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 + r^2 d\Omega^2.
$$

In *isotropic coordinates* with $G=c=1$ (i.e., the coordinate system we'd prefer to use), the ($-+++$ form) line element is written:
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
def ADM_quantities_ito_TOV_soln(rbar, theta, expnu, exp4phi):
    # in TOV ID, betaU=BU=KDD=0
    alpha = sp.sqrt(expnu)

    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()

    gammaDD = ixp.zerorank2(DIM=3)
    gammaDD[0][0] = exp4phi
    gammaDD[1][1] = exp4phi * rbar ** 2
    gammaDD[2][2] = exp4phi * rbar ** 2 * sp.sin(theta) ** 2

    KDD = ixp.zerorank2()

    return alpha, betaU, BU, gammaDD, KDD
```

<a id='tmunu_ito_tov'></a>

### Step 4.d.ii: Stress-energy tensor $T^{\mu\nu}$, written in terms of the TOV solution \[Back to [top](#toc)\]
$$\label{tmunu_ito_tov}$$


We will also need the stress-energy tensor $T^{\mu\nu}$. [As discussed here](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation), the stress-energy tensor is diagonal:

\begin{align}
T^t_t &= -\rho \\
T^i_j &= P \delta^i_j \\
\text{All other components of }T^\mu_\nu &= 0.
\end{align}

Since $\beta^i=0$ the inverse metric expression simplifies to (Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf)):
$$
g^{\mu\nu} = \begin{pmatrix} 
-\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
\frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
\end{pmatrix} =
\begin{pmatrix} 
-\frac{1}{\alpha^2} & 0 \\
0 & \gamma^{ij}
\end{pmatrix},
$$

and since the 3-metric is diagonal we get

\begin{align}
\gamma^{\bar{r}\bar{r}} &= e^{-4\phi}\\
\gamma^{\theta\theta} &= e^{-4\phi}\frac{1}{\bar{r}^2} \\
\gamma^{\phi\phi} &= e^{-4\phi}\frac{1}{\bar{r}^2 \sin^2 \theta}.
\end{align}

Thus raising $T^\mu_\nu$ yields a diagonal $T^{\mu\nu}$

\begin{align}
T^{tt} &= -g^{tt} \rho = \frac{1}{\alpha^2} \rho = e^{-\nu(\bar{r})} \rho \\
T^{\bar{r}\bar{r}} &= g^{\bar{r}\bar{r}} P = \frac{1}{e^{4 \phi}} P \\
T^{\theta\theta} &= g^{\theta\theta} P = \frac{1}{e^{4 \phi}\bar{r}^2} P\\
T^{\phi\phi} &= g^{\phi\phi} P = \frac{1}{e^{4\phi}\bar{r}^2 \sin^2 \theta} P 
\end{align}


```python
def T4UU_ito_TOV_soln(rbar, theta, rho, P, expnu, exp4phi):

    T4UU = ixp.zerorank2(DIM=4)
    T4UU[0][0] = rho/expnu
    T4UU[1][1] = P/exp4phi
    T4UU[2][2] = P/(exp4phi*rbar**2)
    T4UU[3][3] = P/(exp4phi*rbar**2*sp.sin(theta)**2)

    return T4UU
```

<a id='tov_id_function_itself'></a>

### Step 4.d.iii: The `TOV_ID_function()` itself \[Back to [top](#toc)\]
$$\label{tov_id_function_itself}$$


```python
def add_to_Cfunction_dict_TOV_ID_function():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """This function sets TOV initial data at a single point xCart[3] = (x,y,z),
outputting to initial_data_struct *restrict initial_data"""
    c_type = "void"
    name = "TOV_ID_function"
    params = "const paramstruct *params, const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"
    body  = """  // First set r(=rbar), theta, phi in terms of xCart[3]:
  const REAL Cartx=xCart[0], Carty=xCart[1], Cartz=xCart[2];"""
    desired_rfm_coord = par.parval_from_str("reference_metric::CoordSystem")
    # "Spherical" is the native basis for TOV initial data.
    par.set_parval_from_str("reference_metric::CoordSystem", "Spherical")
    rfm.reference_metric()
    body += r"""
  REAL rbar, theta, phi;
  {
""" + outputC(rfm.Cart_to_xx[:3], ["rbar", "theta", "phi"], filename="returnstring",
              params="outCverbose=False,includebraces=False,preindent=2")
    body += r"""
  }

  // Next set gamma_{ij} in spherical basis
  REAL rho,rho_baryon,P,M,expnu,exp4phi;
  TOV_interpolate_1D(rbar,ID_persist, &rho,&rho_baryon,&P,&M,&expnu,&exp4phi);
"""
    rbar, theta, rho, P, expnu, exp4phi = sp.symbols('rbar theta rho P expnu exp4phi', real=True)
    alpha, betaU, BU, gammaDD, KDD = ADM_quantities_ito_TOV_soln(rbar, theta, expnu, exp4phi)
    T4UU = T4UU_ito_TOV_soln(rbar, theta, rho, P, expnu, exp4phi)

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
    for mu in range(4):
        for nu in range(mu, 4):
            list_of_output_exprs += [T4UU[mu][nu]]
            list_of_output_varnames += ["initial_data->T4SphorCartUU" + str(mu) + str(nu)]

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
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
```

<a id='ccodegen_library_validation'></a>

# Step 5: Validate C code generation against `TOV.TOV_Ccodegen_library` Python module \[Back to [top](#toc)\]
$$\label{ccodegen_library_validation}$$

To ensure this documentation remains maintained in lockstep with development of the separate `TOV.TOV_Ccodegen_library` Python module, we register all the C functions both here and from the Python module and compare them. If they are found to be *identical*, the test passes. Otherwise it fails and the notebook will exit with an error.


```python
import TOV.TOV_Ccodegen_library as TCL
import sys

funclist = [("ID_persist_str", ID_persist_str, TCL.ID_persist_str),
            ("add_to_Cfunction_dict_TOV_read_data_file_set_ID_persist", add_to_Cfunction_dict_TOV_read_data_file_set_ID_persist, TCL.add_to_Cfunction_dict_TOV_read_data_file_set_ID_persist),
            ("add_to_Cfunction_dict_TOV_interpolate_1D", add_to_Cfunction_dict_TOV_interpolate_1D, TCL.add_to_Cfunction_dict_TOV_interpolate_1D),
            ("ADM_quantities_ito_TOV_soln", ADM_quantities_ito_TOV_soln, TCL.ADM_quantities_ito_TOV_soln),
            ("T4UU_ito_TOV_soln", T4UU_ito_TOV_soln, TCL.T4UU_ito_TOV_soln),
            ("ADM_quantities_ito_TOV_soln", ADM_quantities_ito_TOV_soln, TCL.ADM_quantities_ito_TOV_soln),
            ("add_to_Cfunction_dict_TOV_ID_function", add_to_Cfunction_dict_TOV_ID_function, TCL.add_to_Cfunction_dict_TOV_ID_function)
           ]

if sys.version_info.major >= 3:
    import inspect
    for func in funclist:
        # https://stackoverflow.com/questions/20059011/check-if-two-python-functions-are-equal
        if inspect.getsource(func[1]) != inspect.getsource(func[2]):
            print("inspect.getsource(func[1]):")
            print(inspect.getsource(func[1]))
            print("inspect.getsource(func[2]):")
            print(inspect.getsource(func[2]))
            print("ERROR: function " + func[0] + " is not the same as the Ccodegen_library version!")
            sys.exit(1)

    print("PASS! ALL FUNCTIONS ARE IDENTICAL")
else:
    print("SORRY CANNOT CHECK FUNCTION IDENTITY WITH PYTHON 2. PLEASE update your Python installation.")
```

    PASS! ALL FUNCTIONS ARE IDENTICAL


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-ADM_Initial_Data-TOV](Tutorial-ADM_Initial_Data-TOV.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-ADM_Initial_Data-TOV")
```

    Created Tutorial-ADM_Initial_Data-TOV.tex, and compiled LaTeX file to PDF
        file Tutorial-ADM_Initial_Data-TOV.pdf

