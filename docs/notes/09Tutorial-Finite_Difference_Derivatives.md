<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# NRPy+'s Finite Difference Interface

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This tutorial examines numerical derivatives via finite difference methods using the NRPy+ module. It outlines FD techniques, which involve creating an Nth-degree polynomial for a function on a uniform grid and approximating its derivatives. Key functions, `compute_fdcoeffs_fdstencl()` and `FD_outputC()`, which determine finite difference coefficients and generate respective C code, are detailed.

### NRPy+ Source Code for this module: [finite_difference.py](../edit/finite_difference.py)


<a id='toc'></a>

# Table of Contents 
$$\label{toc}$$ 

This notebook is organized as follows

1. [Preliminaries](#fdd): Introduction to Finite Difference Derivatives
1. [Step 1](#fdmodule): The finite_difference NRPy+ module
    1. [Step 1.a](#fdcoeffs_func): The `compute_fdcoeffs_fdstencl()` function
        1. [Step 1.a.i](#exercise): Exercise: Using `compute_fdcoeffs_fdstencl()`
    1. [Step 1.b](#fdoutputc): The  `FD_outputC()` function
1. [Step 2](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF    

<a id='fdd'></a>

# Preliminaries: Introduction to Finite Difference Derivatives \[Back to [top](#toc)\]
$$\label{fdd}$$

Suppose we have a *uniform* numerical grid in one dimension; say, the Cartesian $x$ direction. Since the grid is uniform, the spacing between successive grid points is $\Delta x$, and the position of the $i$th point is given by

$$x_i = x_0 + i \Delta x.$$

Then, given a function $u(x)$ on this uniform grid, we will adopt the notation

$$u(x_i) = u_i.$$

We wish to approximate derivatives of $u_i$ at some nearby point (in this tutorial, we will consider derivatives at one of the sampled points $x_i$) using [finite difference](https://en.wikipedia.org/wiki/Finite_difference). (FD) techniques. 

FD techniques are usually constructed as follows:
* First, find the unique $N$th-degree polynomial that passes through $N+1$ sampled points of our function $u$ in the neighborhood of where we wish to find the derivative. 
* Then, provided $u$ is smooth and properly-sampled, the $n$th derivative of the polynomial (where $n\le N-1$; *Exercise: Justify this inequality*) is approximately equal to the $n$th derivative of $u$. We call this the **$n$th-order finite difference derivative of $u$**. 
* So long as the function $u$ is smooth and properly sampled, the relative error between the exact and the finite difference derivative $u^{(n)}$ will generally decrease as the polynomial degree or sampling density increases.

The $n$th finite difference derivative of $u(x)$ at $x=x_i$ can then be written in the form
$$u^{(n)}(x_i)_{\text{FD}} = \sum_{j=0}^{N} u_j a_j,$$
where the $a_j$'s are known as *finite difference coefficients*. So long as the $N$th-degree polynomial that passes through the $N+1$ points is unique, the corresponding set of $a_j$'s is unique as well.

There are multiple ways to compute the finite difference coefficients $a_j$, including solving for the $N$th-degree polynomial that passes through the function at the sampled points. However, the most popular and most straightforward way involves Taylor series expansions about sampled points near the point where we wish to evaluate the derivative.

**Recommended: Learn more about the algorithm NRPy+ adopts to automatically compute finite difference derivatives: ([How NRPy+ Computes Finite Difference Coefficients](Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs.ipynb))**


<a id='fdmodule'></a>

# Step 1: The finite_difference NRPy+ module \[Back to [top](#toc)\]
$$\label{fdmodule}$$

The finite_difference NRPy+ module contains one parameter:

* **FD_CENTDERIVS_ORDER**: An integer indicating the requested finite difference accuracy order (not the order of the derivative) , where FD_CENTDERIVS_ORDER = [the size of the finite difference stencil in each direction, plus one].

The finite_difference NRPy+ module contains two core functions: `compute_fdcoeffs_fdstencl()` and `FD_outputC()`. The first is a low-level function normally called only by `FD_outputC()`, which computes and outputs finite difference coefficients and the numerical grid indices (stencil) corresponding to each coefficient:

<a id='fdcoeffs_func'></a>

## Step 1.a:  The `compute_fdcoeffs_fdstencl()` function \[Back to [top](#toc)\]
$$\label{fdcoeffs_func}$$

**compute_fdcoeffs_fdstencl(derivstring,FDORDER=-1)**:
* Output nonzero finite difference coefficients and corresponding numerical stencil as lists, using as inputs:
    * **derivstring**: indicates the precise type and direction derivative desired, with  options:
        * **centered derivatives**, where the center of the finite difference stencil corresponds to the point where the derivative is desired;
            * For a first-order derivative, set derivstring to "D"+"dirn", where "dirn" is an integer denoting direction. For a second-order derivative, set derivstring to "DD"+"dirn1"+"dirn2", where "dirn1" and "dirn2" are integers denoting the direction of each derivative. Currently, only $1 \le N \le 2$ is supported (extension to higher-order derivatives is straightforward). Examples in 3D Cartesian coordinates (x,y,z):
                * the derivative operator $\partial_x^2$ corresponds to derivstring = "DD00,"
                * the derivative operator $\partial_x \partial_y$ corresponds to derivstring = "DD01,"
                * the derivative operator $\partial_z$ corresponds to derivstring = "D2."
        * **up—or downwinded—derivatives**, where the center of the finite difference stencil is *one gridpoint* up or down from where the derivative is requested;
            * Set derivstring to "upD"+"dirn" or "dnD"+"dirn", where "dirn" is an integer denoting direction. Example in 3D Cartesian coordinates (x,y,z):
                * the upwinded derivative operator $\partial_x$ corresponds to derivstring = "dupD0."
        * and **Kreiss-Oliger dissipation derivatives**, where the center of the finite difference stencil corresponds to the point where the dissipation will be applied.
            * Set derivstring to "dKOD"+"dirn", where "dirn" is an integer denoting direction. Example in 3D Cartesian coordinates (x,y,z):
                * the Kreiss-Oliger derivative operator $\partial_z^\text{KO}$ corresponds to derivstring = "dKOD2."
    * **FDORDER**: an *optional* parameter that, if set to a positive even integer, overrides FD_CENTDERIVS_ORDER.

Within NRPy+, `compute_fdcoeffs_fdstencl()` is only called from `FD_outputC()`. Regardless, this function provides a nice interface for evaluating finite difference coefficients, as shown below:


```python
# Import the finite difference module
import finite_difference as fin  # NRPy+: Finite difference C code generation module

fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl("dDD00")
print(fdcoeffs)
print(fdstencl)
```

    [-1/12, 4/3, -5/2, 4/3, -1/12]
    [[-2, 0, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0], [2, 0, 0, 0]]


Interpreting the output, notice first that $\texttt{fdstencl}$ is a list of coordinate indices, where up to 4 dimension indices are supported (higher dimensions are possible and can be straightforwardly added, though be warned about [The Curse of Dimensionality](https://en.wikipedia.org/wiki/Curse_of_dimensionality)).

Thus NRPy+ found that for some function $u$, the fourth-order accurate finite difference operator at point $x_{i0}$ is given by

$$[\partial_{x}^{2} u]^\text{FD4}_{i0} = \frac{1}{\Delta x^{2}} \left[ -\frac{1}{12} \left(u_{i0-2,i1,i2,i3} + u_{i0+2,i1,i2,i3}\right) - \frac{5}{2}u_{i0,i1,i2,i3} + \frac{4}{3}\left(u_{i0-1,i1,i2,i3} + u_{i0+1,i1,i2,i3}\right)\right]$$

Notice also that multiplying by the appropriate power of $\frac{1}{\Delta x}$ term is up to the user of this function.

In addition, if the gridfunction $u$ exists on a grid that is less than four (spatial) dimensions, it is up to the user to truncate the additional index information.

<a id='exercise'></a>

### Step 1.a.i: Exercise: Using `compute_fdcoeffs_fdstencl()` \[Back to [top](#toc)\]
$$\label{exercise}$$

Using `compute_fdcoeffs_fdstencl()` write the necessary loops to output the finite difference coefficient tables in the Wikipedia article on [finite difference coefficients](https://en.wikipedia.org/wiki/Finite_difference_coefficients), for first and second centered derivatives (i.e., up to $\partial_i^2$)  up to eighth-order accuracy. [Solution, courtesy Brandon Clark](Tutorial-Finite_Difference_Derivatives-FDtable_soln.ipynb).

<a id='fdoutputc'></a>

## Step 1.b: The  `FD_outputC()` function \[Back to [top](#toc)\]
$$\label{fdoutputc}$$

**FD_outputC(filename,sympyexpr_list)**: C code generator for finite-difference expressions.

C codes that evaluate expressions with finite difference derivatives on numerical grids generally consist of three  components, all existing within a loop over "interior" gridpoints; at a given gridpoint, the code must do the following things. 

1. Read gridfunctions from memory at all points needed to evaluate the finite difference derivatives or the gridfunctions themselves.
2. Perform arithmetic, including computation of finite difference stencils.
3. Write the output from the arithmetic to other gridfunctions.

To minimize cache misses and maximize potential compiler optimizations, it is generally recommended to segregate the above three steps. `FD_outputC()` first analyzes the input expressions, searching for derivatives of gridfunctions. The search is very easy, as NRPy+ requires a very specific syntax for derivatives: 
* `gf_dD0` denotes the first derivative of gridfunction "gf" in direction zero.
* `gf_dupD0` denotes the upwinded first derivative of gridfunction "gf" in direction zero.
* `gf_ddnD0` denotes the downwinded first derivative of gridfunction "gf" in direction zero.
* `gf_dKOD2` denotes the Kreiss-Oliger dissipation operator of gridfunction "gf" in direction two.
Each time `FD_outputC()` finds a derivative (including references to the gridfunction directly \["zeroth"-order derivatives\]) in this way, it calls `compute_fdcoeffs_fdstencl()` to record the specific locations in memory from which the underlying gridfunction must be read to evaluate the appropriate finite difference derivative.

`FD_outputC()` then orders this list of points for all gridfunctions and points in memory, optimizing memory reads based on how the gridfunctions are stored in memory (set via parameter MemAllocStyle in the NRPy+ grid module). It then completes step 1. 

For step 2, `FD_outputC()` exports all of the finite difference expressions, as well as the original expressions input into the function, to `outputC()` to generate the optimized C code. Step 3 follows trivially from just being careful with the bookkeeping in the above steps.

`FD_outputC()` takes two arguments:
* **filename**: set to "stdout" to print to screen. Otherwise, specify a filename.
* **sympyexpr_list**: a single named tuple or list of named tuples of type "lhrh", where the lhrh type refers to the simple structure:
    * **lhrh(left-hand side of equation, right-hand side of the equation)**

Time for an example: let's compute 
$$
\texttt{output} = \text{phi_dDD00} = \partial_x^2 \phi(x,t),
$$
where $\phi$ is a function of space and time, though we only store its spatial values at a given time (*a la* the [Method of Lines](https://reference.wolfram.com/language/tutorial/NDSolveMethodOfLines.html), described & implemented in next the [Scalar Wave Equation module](Tutorial-Start_to_Finish-ScalarWave.ipynb)). 

As detailed above, the suffix $\text{_dDD00}$ tells NRPy+ to construct the second finite difference derivative of gridfunction $\texttt{phi}$ with respect to coordinate $xx0$ (in this case $xx0$ is simply the Cartesian coordinate $x$). Here is the NRPy+ implementation:


```python
from outputC import lhrh         # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import finite_difference as fin  # NRPy+: Finite difference C code generation module

# Set the spatial dimension to 1
par.set_paramsvals_value("grid::DIM = 1")

# Register the input gridfunction "phi" and the gridfunction to which data are output, "output":
phi, output = gri.register_gridfunctions("AUX",["phi","output"])

# Declare phi_dDD as a rank-2 indexed expression: phi_dDD[i][j] = \partial_i \partial_j phi
phi_dDD = ixp.declarerank2("phi_dDD","nosym")

# Set output to \partial_0^2 phi
output = phi_dDD[0][0]

# Output to the screen the core C code for evaluating the finite difference derivative
fin.FD_outputC("stdout",lhrh(lhs=gri.gfaccess("out_gf","output"),rhs=output))
```

    {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
       */
      /*
       *  Original SymPy expression:
       *  "const double phi_dDD00 = invdx0**2*(-5*phi/2 + 4*phi_i0m1/3 - phi_i0m2/12 + 4*phi_i0p1/3 - phi_i0p2/12)"
       */
      const double phi_i0m2 = aux_gfs[IDX2S(PHIGF, i0-2)];
      const double phi_i0m1 = aux_gfs[IDX2S(PHIGF, i0-1)];
      const double phi = aux_gfs[IDX2S(PHIGF, i0)];
      const double phi_i0p1 = aux_gfs[IDX2S(PHIGF, i0+1)];
      const double phi_i0p2 = aux_gfs[IDX2S(PHIGF, i0+2)];
      const double FDPart1_Rational_5_2 = 5.0/2.0;
      const double FDPart1_Rational_1_12 = 1.0/12.0;
      const double FDPart1_Rational_4_3 = 4.0/3.0;
      const double phi_dDD00 = ((invdx0)*(invdx0))*(FDPart1_Rational_1_12*(-phi_i0m2 - phi_i0p2) + FDPart1_Rational_4_3*(phi_i0m1 + phi_i0p1) - FDPart1_Rational_5_2*phi);
      /*
       * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
       */
      /*
       *  Original SymPy expression:
       *  "aux_gfs[IDX2S(OUTPUTGF, i0)] = phi_dDD00"
       */
      aux_gfs[IDX2S(OUTPUTGF, i0)] = phi_dDD00;
    }


Some important points about the above code follow.
* The gridfunction PHIGF samples some function $\phi(x)$ at discrete uniform points in $x$, labeled $x_i$ at all points $i\in [0,N]$, so that 
$$\phi(x_i) = \phi_{i}=\text{in_gfs[IDX2(PHIGF, i)]}.$$ 
* For a *uniformly* sampled function with constant grid spacing (sample rate) $\Delta x$, $x_i$ is defined as $x_i = x_0 + i \Delta x$.
* The variable $\texttt{invdx0}$ must be defined by the user in terms of the uniform gridspacing $\Delta x$ as $\texttt{invdx0} = \frac{1}{\Delta x}$. 
     * *Aside*: why do we choose to multiply by $1/\Delta x$ instead of dividing the expression by $\Delta x$, which would seem much more straightforward? 
         * *Answer*: as discussed in the [first part of the tutorial](Tutorial-Coutput__Parameter_Interface.ipynb), division of floating-point numbers on modern CPUs is far more expensive than multiplication, usually by a factor of ~3 or more.

<a id='latex_pdf_output'></a>

# Step 2: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Finite_Difference_Derivatives.pdf](Tutorial-Finite_Difference_Derivatives.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Finite_Difference_Derivatives")
```

    Created Tutorial-Finite_Difference_Derivatives.tex, and compiled LaTeX file
        to PDF file Tutorial-Finite_Difference_Derivatives.pdf

