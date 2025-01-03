<img height="60px" src="https://blackholesathome.net/images/NRPyPlusLogo.png" align="center" hspace="20px" vspace="5px">
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

## The NRPy+ Tutorial: An Introduction to Python-Based Code Generation for Numerical Relativity... and Beyond!

### Lead author: [Zachariah B. Etienne](https://etienneresearch.com) $\leftarrow$ Please feel free to email comments, revisions, or errata! 

***If you are unfamiliar with using Jupyter Notebooks, first review the official [Jupyter Notebook Basics Guide](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Notebook%20Basics.html).***

### PART 1: Basic Functionality of NRPy+, a First Application
##### NRPy+ Basics
+ [NRPy+: Introduction & Motivation](https://nrpyplus.net) (NRPy+ home page)
+ [NRPy+: 10-minute Overview](Tutorial-NRPyPlus_10_Minute_Overview.ipynb) $\leftarrow$ **Start here**
+ [Solving the scalar wave equation with `NumPy`](Tutorial-Solving_the_Scalar_Wave_Equation_with_NumPy.ipynb) $\leftarrow$ **Then go here** (fully worked example of a hyperbolic PDE solver, written in pure Python)
+ [Basic C Code Output, NRPy+'s Parameter Interface](Tutorial-Coutput__Parameter_Interface.ipynb)
+ [`cmdline_helper`: Multi-platform command-line helper functions](Tutorial-cmdline_helper.ipynb) (*Courtesy Brandon Clark*)
+ [Numerical Grids](Tutorial-Numerical_Grids.ipynb)
+ [Indexed Expressions (e.g., tensors, pseudotensors, etc.)](Tutorial-Indexed_Expressions.ipynb)
+ [Loop Generation](Tutorial-Loop_Generation_Cache_Blocking.ipynb)
+ [Finite Difference Derivatives](Tutorial-Finite_Difference_Derivatives.ipynb)
    + Instructional notebook: [How NRPy+ Computes Finite Difference Derivative Coefficients](Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs.ipynb)
    + **Start-to-Finish Example**: [Finite-Difference Playground: A Complete C Code for Validating NRPy+-Based Finite Differences](Tutorial-Start_to_Finish-Finite_Difference_Playground.ipynb)
+ [NRPy+ SymPy LaTeX Interface](Tutorial-SymPy_LaTeX_Interface.ipynb) (*Courtesy Ken Sible*)
    + **Start-to-Finish Example**: [LaTeX Interface: BSSN (Cartesian)](Tutorial-LaTeX_Interface_Example-BSSN_Cartesian.ipynb)
+ Method of Lines for PDEs: Step PDEs forward in time using ODE methods
    + Solving ODEs using explicit Runge Kutta methods (*Courtesy Brandon Clark*)
        + [The Family of explicit Runge-Kutta methods and their Butcher tables](Tutorial-RK_Butcher_Table_Dictionary.ipynb)
        + [Validating Runge Kutta Butcher tables using truncated Taylor series](Tutorial-RK_Butcher_Table_Validation.ipynb)
    + [Generating C Code to implement Method of Lines timestepping with explicit Runge Kutta-like methods](Tutorial-Method_of_Lines-C_Code_Generation.ipynb) (*Courtesy Brandon Clark*)
+ Convenient mathematical operations
    + [Representing `min(a,b)` and `max(a,b)` without `if()` statements; defining piecewise functions](Tutorial-Min_Max_and_Piecewise_Expressions.ipynb) (*Courtesy Patrick Nelson*)
    + [Symbolic Tensor Rotation using Quaternions](Tutorial-Symbolic_Tensor_Rotation.ipynb) (*Courtesy Ken Sible*)
+ Contributing to NRPy+
    + [The NRPy+ Tutorial Style Guide](Tutorial-Template_Style_Guide.ipynb) (*Courtesy Brandon Clark*)
    + [Adding Unit Tests](Tutorial-UnitTesting.ipynb) (*Courtesy Kevin Lituchy*)

### PART 2: Basic Physics Applications 
##### Using NRPy+ to Numerically Solve the Scalar Wave Equation in Cartesian Coordinates
+ Application: [The Scalar **Wave Equation** in Cartesian Coordinates](Tutorial-ScalarWave.ipynb)
    + **Start-to-Finish Example**: [Numerically Solving the Scalar Wave Equation: A Complete C Code](Tutorial-Start_to_Finish-ScalarWave.ipynb)
    + Solving the Wave Equation with the <font color='green'>**Einstein Toolkit**</font> (*Courtesy Patrick Nelson & Terrence Pierre Jacques*)
        + [<font color='green'>**IDScalarWaveNRPy**</font>: Initial data for the scalar wave equation](Tutorial-ETK_thorn-IDScalarWaveNRPy.ipynb) 
        + [<font color='green'>**WaveToyNRPy**</font>: Solving the scalar wave equation, using the method of lines](Tutorial-ETK_thorn-WaveToyNRPy.ipynb)

##### Diagnostic Notebooks: Gravitational Wave Extraction in Cartesian coordinates
+ Application: [All Weyl scalars and invariants in Cartesian Coordinates](Tutorial-WeylScalarsInvariants-Cartesian.ipynb) (*Courtesy Patrick Nelson*)
    + [<font color='green'>**WeylScal4NRPy**</font>: An <font color='green'>**Einstein Toolkit**</font> Diagnostic Thorn](Tutorial-ETK_thorn-Weyl_Scalars_and_Spacetime_Invariants.ipynb) (*Courtesy Patrick Nelson*)


##### Solving the Effective-One-Body Equations of Motion

+ Application: [SEOBNR: The Spinning-Effective-One-Body-Numerical-Relativity Hamiltonian, version 3](Tutorial-SEOBNR_v3_Hamiltonian.ipynb)
    + Solving the SEOBNR Hamiltonian equations of motion (<font color='red'>**in progress**</font>)
        + [Initial data: Setting the initial trajectory](in_progress/Tutorial-SEOBNR_Initial_Conditions.ipynb)
        + [SymPy-generated exact derivatives of the SEOBNR Hamiltonian](in_progress/Tutorial-Spinning_Effective_One_Body_Numerical_Relativity_Hamiltonian-Cartesian.ipynb)
    + [The SEOBNRv4P Hamiltonian](Tutorial-SEOBNR_v4P_Hamiltonian.ipynb)

##### NRPyPN: Validated Post-Newtonian Expressions for Input into Wolfram Mathematica, SymPy, or Highly-Optimized C Codes

+ [NRPyPN Main Menu](NRPyPN/NRPyPN.ipynb) $\leftarrow$ includes NRPyPN Table of Contents and a quick interface for setting up low-eccentricity (up to 3.5 PN order) momentum parameters for binary black hole initial data.

### PART 3: Solving PDEs in Curvilinear Coordinate Systems
+ [Moving beyond Cartesian Grids: Reference Metrics](Tutorial-Reference_Metric.ipynb)
+ **Start-to-Finish Example**: [Implementation of Curvilinear Boundary Conditions, Including for Tensorial Quantities](Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb)
+ Application: [**The Scalar Wave Equation** in Curvilinear Coordinates, using a Reference Metric](Tutorial-ScalarWaveCurvilinear.ipynb)
    + **Start-to-Finish Example**: [Numerically Solving the Scalar Wave Equation in Curvilinear Coordinates: A Complete C Code](Tutorial-Start_to_Finish-ScalarWaveCurvilinear.ipynb)

### PART 4: Numerical Relativity $-$ BSSN in Curvilinear Coordinates

+ [**Overview: Covariant BSSN formulation of general relativity in curvilinear coordinates**](Tutorial-BSSN_formulation.ipynb)
    + [Construction of useful BSSN quantities](Tutorial-BSSN_quantities.ipynb)
    + [BSSN time-evolution equations](Tutorial-BSSN_time_evolution-BSSN_RHSs.ipynb)
    + [Time-evolution equations for BSSN gauge quantities $\alpha$ and $\beta^i$](Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb)
    + [Hamiltonian and momentum constraint equations](Tutorial-BSSN_constraints.ipynb)
    + [Enforcing the conformal 3-metric $\det{\bar{\gamma}_{ij}}=\det{\hat{\gamma}_{ij}}$ constraint](Tutorial-BSSN_enforcing_determinant_gammabar_equals_gammahat_constraint.ipynb)
    + [Writing quantities of ADM formalism in terms of BSSN quantities](Tutorial-ADM_in_terms_of_BSSN.ipynb)
    + [Basis transformations of BSSN variables](Tutorial-BSSN-basis_transforms.ipynb)
+ **Initial data notebooks**. Initial data are set in terms of standard [ADM formalism](https://en.wikipedia.org/wiki/ADM_formalism) spacetime quantities.
    + [Non-Spinning ("static trumpet") black hole initial data](Tutorial-ADM_Initial_Data-StaticTrumpet.ipynb) (*Courtesy Terrence Pierre Jacques & Ian Ruchlin*)
    + [Spinning UIUC black hole initial data](Tutorial-ADM_Initial_Data-UIUC_BlackHole.ipynb) (*Courtesy Terrence Pierre Jacques & Ian Ruchlin*)
    + [Spinning Shifted Kerr-Schild black hole initial data](Tutorial-ADM_Initial_Data-ShiftedKerrSchild.ipynb) (*Courtesy George Vopal*)
    + [Brill-Lindquist initial data: Two-black-holes released from rest](Tutorial-ADM_Initial_Data-Brill-Lindquist.ipynb)
    + [Black hole accretion disk initial data (Fishbone-Moncrief)](Tutorial-FishboneMoncriefID.ipynb)
        + [<font color='green'>**FishboneMoncriefID**</font>: Setting up Fishbone-Moncrief disk initial data within the <font color='green'>**Einstein Toolkit**</font>, using HydroBase variables as input](Tutorial-ETK_thorn-FishboneMoncriefID.ipynb)
    + [Neutron Star initial data: The Tolman-Oppenheimer-Volkoff (TOV) solution](Tutorial-ADM_Initial_Data-TOV.ipynb) (*Courtesy Phil Chang*)
        + [Implementation of Single and Piecewise Polytropic EOSs](Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb) (*Courtesy Leo Werneck*)
+ **BSSN initial data reader**
    + [Initial Data Reader/Converter: ADM Spherical/Cartesian to BSSN-with-a-reference-metric](Tutorial-ADM_Initial_Data_Reader__BSSN_Converter.ipynb) (Use this module for initial data conversion if the initial data are known *exactly*. The BSSN quantity $\lambda^i$ will be computed exactly using SymPy from given ADM quantities.)
        + [**Start-to-Finish *exact* initial data validation notebook**](Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb): Confirms all exact initial data types listed above satisfy Einstein's equations of general relativity. (*Courtesy Brandon Clark & George Vopal*)
        + [**Start-to-Finish *numerical* initial data validation notebook**: The TOV solution](Tutorial-Start_to_Finish-BSSNCurvilinear-TOV_Initial_Data.ipynb): Neutron star initial data, confirms numerical errors converge to zero at expected order. (TOV initial data are generated via [the *numerical* solution of a system of ODEs](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation), thus are known only numerically)
+ **Diagnostic curvilinear BSSN modules**
    + [The gravitational wave Weyl scalar $\psi_4$, in arbitrary curvilinear coordinates](Tutorial-Psi4.ipynb)
        + [Constructing the quasi-Kinnersley tetrad for $\psi_4$](Tutorial-Psi4_tetrads.ipynb)
            + [Start-to-Finish validation of above expressions in Cartesian coordinates, against Patrick Nelson's Weyl Scalars & Invariants notebook](BSSN/Psi4Cartesianvalidation/Tutorial-Psi4-Cartesian_validation.ipynb)
        + [Spin-weighted spherical harmonics](Tutorial-SpinWeighted_Spherical_Harmonics.ipynb) (*Courtesy Brandon Clark*)

+ **Start-to-Finish curvilinear BSSN simulation examples**:
    + [<font color='purple'>**Colliding black holes!**</font>](Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide.ipynb)
    + [The "Hydro without Hydro" test: evolving the spacetime fields of TOV star with $T^{\mu\nu}$ assumed static](Tutorial-Start_to_Finish-BSSNCurvilinear-Neutron_Star-Hydro_without_Hydro.ipynb)

### PART 5: Numerical Relativity $-$ General Relativistic Hydrodynamics (GRHD), Force-Free Electrodynamics (GRFFE), & Magnetohydrodynamics (GRMHD)

+ [The equations of general relativistic hydrodynamics (**GRHD**), in Cartesian coordinates](Tutorial-GRHD_Equations-Cartesian.ipynb)
+ [The equations of general relativistic, force-free electrodynamics (**GRFFE**), in Cartesian coordinates](Tutorial-GRFFE_Equations-Cartesian.ipynb)
+ [The equations of general relativistic magnetohydrodynamics (**GRMHD**), in Cartesian coordinates](Tutorial-GRMHD_Equations-Cartesian.ipynb)
    + [Derivation of the GRMHD evolution equations from $T^{\mu\nu}$ in our formalism](Tutorial-Derivation_of_GRMHD_Evolution_Equations.ipynb)
+ [**`IllinoisGRMHD`** ](IllinoisGRMHD/doc/) with piecewise-polytrope/hybrid equation of state support (*Courtesy Leo Werneck*): <font color='red'>***In progress***</font>
    + Diagnostic notebook: [<font color='green'>**sbPoynETNRPy**</font>: Evaluating $b^\mu$ and $S^i$ in the <font color='green'>**Einstein Toolkit**</font>, using HydroBase variables as input](Tutorial-ETK_thorn-u0_smallb_Poynting.ipynb)
         + Tutorial notebook: [Computing the 4-Velocity Time-Component $u^0$, the Magnetic Field Measured by a Comoving Observer $b^{\mu}$, and the Poynting Vector $S^i$](Tutorial-u0_smallb_Poynting-Cartesian.ipynb)


```python

```


```python

```
