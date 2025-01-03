# NRPy+ 10-Minute Overview

## Author: Zach Etienne

## This notebook presents NRPy+, a Python framework that converts intricate SymPy expressions into optimized C/C++ code for use in numerical relativity. It emphasizes the computation of 3-Christoffel symbols $\Gamma^{i}_{jk}$ and numerical derivative designations within expressions. 
[comment]: <> (added this abstract in case. omit of it's unnecessary.)

<a id='toc'></a>

# Table of Contents

This notebook is organized as follows

1. [Part 1](#why): Why NRPy+?
1. [Part 2](#christoffel_symbols): Constructing 3-Christoffels $\Gamma^i_{jk}$ as symbolic expressions in terms of 3-metric $\gamma_{ij}$ and derivatives
1. [Part 3](#outputc): `outputC()` example: Output $\Gamma^i_{jk}$ expressions as optimized C code, assuming derivatives already specified
1. [Part 4](#fd_outputc): `FD_outputC()` example: Specify numerical derivatives within $\Gamma^i_{jk}$ expressions as finite differences
1. [Part 5](#what_next): What next? Navigating the NRPy+ tutorial
1. [Part 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF    
# Part 1: Why NRPy+? \[Back to [top](#toc)\]
A core problem faced by numerical relativity is not that the techniques we use are particularly complex, but that *the equations numerical relativists solve are very complex*.

[Einstein notation](https://en.wikipedia.org/wiki/Einstein_notation) certainly makes the complexity more manageable, and is the most common way to perform tensorial mathematics in numerical relativity.

NRPy+ is built upon the idea that equations written in Einstein notation *are a form of code*, and when this code is expressed in NRPy+'s native (Python) language

* rank-0 tensors (scalars) are *variables*
* rank-1 tensors are *lists*
* rank-2 tensors are *lists of lists*
* ... and so forth.

Further, implied tensor summations are simply *loops*.

NRPy+ combines the above ideas with the incredibly powerful [SymPy](https://www.sympy.org/) computer algebra package (think: Mathematica for Python, but fully free & open source), with a custom-built code generation infrastructure for converting complex SymPy expressions directly into highly-optimized C/C++ code.

Importantly, NRPy+ builds on the idea of *learning by example*. The core NRPy+ repository contains more than 100 pedagogical and well-formatted Jupyter notebook tutorials, covering topics of core relevance to the field of numerical relativity. About 25 of these tutorials generate *complete C codes* capable of e.g., evolving the scalar wave equation, Maxwell's equations, and Einstein's equations of general relativity -- all in a variety of coordinate systems. All ~100 Jupyter notebooks are linked to from [the main NRPy+ tutorial table of contents.](NRPyPlus_Tutorial.ipynb)

This 10-minute overview however is designed to introduce only the very basic features of NRPy+, with a core focus on the idea that *NRPy+ can be used to benefit any numerical relativity code*.
# Part 2: Constructing 3-Christoffels $\Gamma^i_{jk}$ as symbolic expressions in terms of 3-metric $\gamma_{ij}$ and derivatives 

**Problem statement**: Given a three-metric $\gamma_{ij}$, construct all 18 independent Christoffel symbols $\Gamma^i_{jk}$, which involves first derivatives of the metric. Assume that $\gamma_{ij}$ *and its derivatives* are given numerically, requiring the derivatives to be defined as symbols.

In NRPy+ we adopt a rigid syntax for tensors and indexed expressions involving Python lists, so that for example

* $\gamma_{ij}=$ `gammaDD[i][j]`
* $\gamma_{ij,k}=$ `gammaDD_dD[i][j][k]`

Christoffel symbols (of the first kind) are defined as ([source](https://en.wikipedia.org/wiki/Christoffel_symbols#Christoffel_symbols_of_the_first_kind)):

$$
\begin{align}
\Gamma_{ij}^k &= \frac{1}{2} \gamma^{kl}\left(\gamma_{jl,i} + \gamma_{il,j} - \gamma_{ij,l}\right)\\
\end{align}
$$
So first we'll define $\gamma_{ij}$ and its inverse using NRPy+ functions


```python
import sympy as sp        # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# gammaDD is a rank-2 indexed expression symmetric in indexes 0 and 1
gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01", DIM=3)

# gammaUU is the inverse of gammaDD, and
# gammaDET (unused) is the determinant of gammaDD
gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)  # inverts a 3x3 symmetric matrix
```

Below we interact a bit with these generated expressions, confirming they 

*Exercise to students*: Verify that $\gamma_{ij} \gamma^{ij} = 3$ using the provided data structures & Python loops. (*Hint 1*: Scroll down a couple of cells to see how the Christoffel symbols are implemented.  *Hint 2*: You will need to use SymPy's `simplify()` function to simplify the expression obtained.)


```python
print("Check that gammaDD[0][1] = gammaDD[1][0]:")
print("gammaDD[0][1] = ", gammaDD[0][1])
print("gammaDD[1][0] = ", gammaDD[1][0])

print("\nOutput gammaUU[1][0] = ", gammaUU[1][0])
print("\nCheck that gammaUU[1][0] - gammaUU[0][1] = 0:")
print("gammaUU[1][0] - gammaUU[0][1] = ", gammaUU[1][0] - gammaUU[0][1])
```

    Check that gammaDD[0][1] = gammaDD[1][0]:
    gammaDD[0][1] =  gammaDD01
    gammaDD[1][0] =  gammaDD01
    
    Output gammaUU[1][0] =  (-gammaDD01*gammaDD22 + gammaDD02*gammaDD12)/(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*gammaDD12**2 - gammaDD01**2*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - gammaDD02**2*gammaDD11)
    
    Check that gammaUU[1][0] - gammaUU[0][1] = 0:
    gammaUU[1][0] - gammaUU[0][1] =  0


Define Christoffel symbols in terms of the inverse metric and metric first derivatives:
$$
\Gamma_{ij}^k = \frac{1}{2} \gamma^{kl}\left(\gamma_{jl,i} + \gamma_{il,j} - \gamma_{ij,l}\right)
$$


```python
# First define symbolic expressions for metric derivatives
gammaDD_dD = ixp.declarerank3("gammaDD_dD", symmetry="sym01", DIM=3)

# Initialize GammaUDD (3-Christoffel) to zero
GammaUDD = ixp.zerorank3(DIM=3)
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                GammaUDD[k][i][j] += sp.Rational(1, 2) * gammaUU[k][l] * (
                    gammaDD_dD[j][l][i] + gammaDD_dD[i][l][j] - gammaDD_dD[i][j][l])
```

Now let's confirm that $\Gamma^{k}_{ij} = \Gamma^k_{ji}$:


```python
for i in range(3):
    for j in range(3):
        for k in range(3):
            print(GammaUDD[i][j][k] - GammaUDD[i][k][j])
```

# Part 3:  `outputC()` example: Output $\Gamma^i_{jk}$ expressions as optimized C code, assuming derivatives already specified 

At the core of NRPy+ is the ability to convert SymPy expressions to highly optimized C code.

**Problem statement**: Output all 18 unique Christoffel symbols with 3 different levels of optimization supported by NRPy+'s core C output routine `outputC()`.

First we store all 18 unique Christoffel symbols, as well as their desired variable names in C, to Python lists.

```python
symbols_list = []
varname_list = []
for i in range(3):
    for j in range(3):
        for k in range(j, 3):
            symbols_list += [GammaUDD[i][j][k]]
            varname_list += ["ChristoffelUDD" + str(i) + str(j) + str(k)]
```

Next we input these lists into NRPy+'s C/C++ code generation function `outputC()`, at three different levels of optimization.

**Optimization Level 0 (don't ever do it this way)**: Compute each Christoffel symbol independently.

```python
import outputC as outC  # NRPy+: Core C code output module
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=False,outCverbose=False")
```

Notice in the above code that many expressions are recomputed time and time again. This is *incredibly inefficient*, and generally compilers won't optimize this properly. So in optimization level 1, we use common subexpression elimination.

**Optimization Level 1**: Use [common subexpression elimination (CSE)](https://en.wikipedia.org/wiki/Common_subexpression_elimination) to group common subexpressions.

```python
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=True,outCverbose=False")
```

**Optimization Level 2**: Use CSE and take advantage of [single instruction, multiple data (SIMD) macros](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data). NRPy+ translates these macros into assembler-level-optimized code for (x86_64) CPUs. Use case: looping over data on a numerical grid; can evaluate expressions at up to 8 gridpoints simultaneously *per CPU core*. Note that `FusedMulSubSIMD(a,b,c)` and `FusedMulAddSIMD(a,b,c)` perform *two* floating-point operations per cycle (on modern CPUs), whereas for other arithmetic operations at most one FP operation can be performed per cycle.


```python
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=True,enable_SIMD=True,outCverbose=False")
```

# Part 4:  `FD_outputC()` example: Specify numerical derivatives within $\Gamma^i_{jk}$ expressions as finite differences 

To emphasize the infrastructure-agnostic nature of NRPy+, in the above C codes we assumed that the derivatives were already computed (e.g., using finite-difference derivatives, discontinuous Galerkin methods, pseudospectral methods, finite-element methods, etc.)

In this part, we demonstrate NRPy+'s `FD_outputC()` code, which basically prepends the above C/C++ codes with the code needed for computing arbitrary-order finite-difference derivatives.

To do this, NRPy+ makes the standard assumption that the underlying grid is *uniform* and that derivatives are taken with respect to functions stored at each point on the numerical grid. Appropriately, we call such functions "gridfunctions".

In the following, we assume that the input, $\gamma_{ij}$, is stored at each gridpoint. Also, we wish to store each independent component of $\Gamma^i_{jk}$ at each gridpoint as output. The below code cell also introduces the NRPy+ parameter interface (accessed via `NRPy_param_funcs.py`) to set the finite-differencing order to 6.

We'll use Optimization Level 1, to make the code easier to read; change `enable_SIMD=False` in the below `FD_outputC()` function call to `enable_SIMD=True` to see the Optimization Level 2 version.


```python
import grid as gri               # NRPy+: Functions having to do with numerical grids
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: parameter interface

# First specify the finite-differencing order to be 6th order (all even orders > 0 supported!)
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 6)

# Next register gridfunctions for Christoffel symbols
gri.register_gridfunctions("AUX", varname_list, rank=3, is_indexed=True, DIM=3)

# Then register gamma_{ij} as gridfunctions
ixp.register_gridfunctions_for_single_rank2("EVOL", "gammaDD", "sym01", DIM=3)

# Finally output the C code using FD_outputC()
# Step 1: FD_outputC() requires that left-hand side and right-hand side of each
#     expression be specified in a named-tuple "lhrh" defined in outputC.py
outputC_lhrh = []
for i, varname in enumerate(varname_list):
    outputC_lhrh += [outC.lhrh(lhs=gri.gfaccess("out_gf", varname), rhs=symbols_list[i])]
# Step 2: call FD_outputC.
fin.FD_outputC("stdout", outputC_lhrh, params="CSE_enable=True,enable_SIMD=False,outCverbose=False")
```

    {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
       */
      const double gammaDD00_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2-3)];
      const double gammaDD00_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2-2)];
      const double gammaDD00_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2-1)];
      const double gammaDD00_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1-3,i2)];
      const double gammaDD00_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1-2,i2)];
      const double gammaDD00_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1-1,i2)];
      const double gammaDD00_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0-3,i1,i2)];
      const double gammaDD00_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0-2,i1,i2)];
      const double gammaDD00_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0-1,i1,i2)];
      const double gammaDD00 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2)];
      const double gammaDD00_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0+1,i1,i2)];
      const double gammaDD00_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0+2,i1,i2)];
      const double gammaDD00_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0+3,i1,i2)];
      const double gammaDD00_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1+1,i2)];
      const double gammaDD00_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1+2,i2)];
      const double gammaDD00_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1+3,i2)];
      const double gammaDD00_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2+1)];
      const double gammaDD00_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2+2)];
      const double gammaDD00_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD00GF, i0,i1,i2+3)];
      const double gammaDD01_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2-3)];
      const double gammaDD01_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2-2)];
      const double gammaDD01_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2-1)];
      const double gammaDD01_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1-3,i2)];
      const double gammaDD01_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1-2,i2)];
      const double gammaDD01_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1-1,i2)];
      const double gammaDD01_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0-3,i1,i2)];
      const double gammaDD01_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0-2,i1,i2)];
      const double gammaDD01_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0-1,i1,i2)];
      const double gammaDD01 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2)];
      const double gammaDD01_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0+1,i1,i2)];
      const double gammaDD01_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0+2,i1,i2)];
      const double gammaDD01_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0+3,i1,i2)];
      const double gammaDD01_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1+1,i2)];
      const double gammaDD01_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1+2,i2)];
      const double gammaDD01_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1+3,i2)];
      const double gammaDD01_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2+1)];
      const double gammaDD01_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2+2)];
      const double gammaDD01_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD01GF, i0,i1,i2+3)];
      const double gammaDD02_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2-3)];
      const double gammaDD02_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2-2)];
      const double gammaDD02_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2-1)];
      const double gammaDD02_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1-3,i2)];
      const double gammaDD02_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1-2,i2)];
      const double gammaDD02_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1-1,i2)];
      const double gammaDD02_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0-3,i1,i2)];
      const double gammaDD02_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0-2,i1,i2)];
      const double gammaDD02_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0-1,i1,i2)];
      const double gammaDD02 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2)];
      const double gammaDD02_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0+1,i1,i2)];
      const double gammaDD02_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0+2,i1,i2)];
      const double gammaDD02_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0+3,i1,i2)];
      const double gammaDD02_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1+1,i2)];
      const double gammaDD02_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1+2,i2)];
      const double gammaDD02_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1+3,i2)];
      const double gammaDD02_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2+1)];
      const double gammaDD02_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2+2)];
      const double gammaDD02_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD02GF, i0,i1,i2+3)];
      const double gammaDD11_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2-3)];
      const double gammaDD11_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2-2)];
      const double gammaDD11_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2-1)];
      const double gammaDD11_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1-3,i2)];
      const double gammaDD11_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1-2,i2)];
      const double gammaDD11_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1-1,i2)];
      const double gammaDD11_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0-3,i1,i2)];
      const double gammaDD11_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0-2,i1,i2)];
      const double gammaDD11_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0-1,i1,i2)];
      const double gammaDD11 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2)];
      const double gammaDD11_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0+1,i1,i2)];
      const double gammaDD11_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0+2,i1,i2)];
      const double gammaDD11_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0+3,i1,i2)];
      const double gammaDD11_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1+1,i2)];
      const double gammaDD11_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1+2,i2)];
      const double gammaDD11_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1+3,i2)];
      const double gammaDD11_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2+1)];
      const double gammaDD11_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2+2)];
      const double gammaDD11_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD11GF, i0,i1,i2+3)];
      const double gammaDD12_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2-3)];
      const double gammaDD12_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2-2)];
      const double gammaDD12_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2-1)];
      const double gammaDD12_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1-3,i2)];
      const double gammaDD12_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1-2,i2)];
      const double gammaDD12_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1-1,i2)];
      const double gammaDD12_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0-3,i1,i2)];
      const double gammaDD12_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0-2,i1,i2)];
      const double gammaDD12_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0-1,i1,i2)];
      const double gammaDD12 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2)];
      const double gammaDD12_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0+1,i1,i2)];
      const double gammaDD12_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0+2,i1,i2)];
      const double gammaDD12_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0+3,i1,i2)];
      const double gammaDD12_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1+1,i2)];
      const double gammaDD12_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1+2,i2)];
      const double gammaDD12_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1+3,i2)];
      const double gammaDD12_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2+1)];
      const double gammaDD12_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2+2)];
      const double gammaDD12_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD12GF, i0,i1,i2+3)];
      const double gammaDD22_i0_i1_i2m3 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2-3)];
      const double gammaDD22_i0_i1_i2m2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2-2)];
      const double gammaDD22_i0_i1_i2m1 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2-1)];
      const double gammaDD22_i0_i1m3_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1-3,i2)];
      const double gammaDD22_i0_i1m2_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1-2,i2)];
      const double gammaDD22_i0_i1m1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1-1,i2)];
      const double gammaDD22_i0m3_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0-3,i1,i2)];
      const double gammaDD22_i0m2_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0-2,i1,i2)];
      const double gammaDD22_i0m1_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0-1,i1,i2)];
      const double gammaDD22 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2)];
      const double gammaDD22_i0p1_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0+1,i1,i2)];
      const double gammaDD22_i0p2_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0+2,i1,i2)];
      const double gammaDD22_i0p3_i1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0+3,i1,i2)];
      const double gammaDD22_i0_i1p1_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1+1,i2)];
      const double gammaDD22_i0_i1p2_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1+2,i2)];
      const double gammaDD22_i0_i1p3_i2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1+3,i2)];
      const double gammaDD22_i0_i1_i2p1 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2+1)];
      const double gammaDD22_i0_i1_i2p2 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2+2)];
      const double gammaDD22_i0_i1_i2p3 = in_gfs[IDX4S(GAMMADD22GF, i0,i1,i2+3)];
      const double FDPart1_Rational_3_4 = 3.0/4.0;
      const double FDPart1_Rational_3_20 = 3.0/20.0;
      const double FDPart1_Rational_1_60 = 1.0/60.0;
      const double gammaDD_dD000 = invdx0*(FDPart1_Rational_1_60*(-gammaDD00_i0m3_i1_i2 + gammaDD00_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD00_i0m2_i1_i2 - gammaDD00_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD00_i0m1_i1_i2 + gammaDD00_i0p1_i1_i2));
      const double gammaDD_dD001 = invdx1*(FDPart1_Rational_1_60*(-gammaDD00_i0_i1m3_i2 + gammaDD00_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD00_i0_i1m2_i2 - gammaDD00_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD00_i0_i1m1_i2 + gammaDD00_i0_i1p1_i2));
      const double gammaDD_dD002 = invdx2*(FDPart1_Rational_1_60*(-gammaDD00_i0_i1_i2m3 + gammaDD00_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD00_i0_i1_i2m2 - gammaDD00_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD00_i0_i1_i2m1 + gammaDD00_i0_i1_i2p1));
      const double gammaDD_dD010 = invdx0*(FDPart1_Rational_1_60*(-gammaDD01_i0m3_i1_i2 + gammaDD01_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD01_i0m2_i1_i2 - gammaDD01_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD01_i0m1_i1_i2 + gammaDD01_i0p1_i1_i2));
      const double gammaDD_dD011 = invdx1*(FDPart1_Rational_1_60*(-gammaDD01_i0_i1m3_i2 + gammaDD01_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD01_i0_i1m2_i2 - gammaDD01_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD01_i0_i1m1_i2 + gammaDD01_i0_i1p1_i2));
      const double gammaDD_dD012 = invdx2*(FDPart1_Rational_1_60*(-gammaDD01_i0_i1_i2m3 + gammaDD01_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD01_i0_i1_i2m2 - gammaDD01_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD01_i0_i1_i2m1 + gammaDD01_i0_i1_i2p1));
      const double gammaDD_dD020 = invdx0*(FDPart1_Rational_1_60*(-gammaDD02_i0m3_i1_i2 + gammaDD02_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD02_i0m2_i1_i2 - gammaDD02_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD02_i0m1_i1_i2 + gammaDD02_i0p1_i1_i2));
      const double gammaDD_dD021 = invdx1*(FDPart1_Rational_1_60*(-gammaDD02_i0_i1m3_i2 + gammaDD02_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD02_i0_i1m2_i2 - gammaDD02_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD02_i0_i1m1_i2 + gammaDD02_i0_i1p1_i2));
      const double gammaDD_dD022 = invdx2*(FDPart1_Rational_1_60*(-gammaDD02_i0_i1_i2m3 + gammaDD02_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD02_i0_i1_i2m2 - gammaDD02_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD02_i0_i1_i2m1 + gammaDD02_i0_i1_i2p1));
      const double gammaDD_dD110 = invdx0*(FDPart1_Rational_1_60*(-gammaDD11_i0m3_i1_i2 + gammaDD11_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD11_i0m2_i1_i2 - gammaDD11_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD11_i0m1_i1_i2 + gammaDD11_i0p1_i1_i2));
      const double gammaDD_dD111 = invdx1*(FDPart1_Rational_1_60*(-gammaDD11_i0_i1m3_i2 + gammaDD11_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD11_i0_i1m2_i2 - gammaDD11_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD11_i0_i1m1_i2 + gammaDD11_i0_i1p1_i2));
      const double gammaDD_dD112 = invdx2*(FDPart1_Rational_1_60*(-gammaDD11_i0_i1_i2m3 + gammaDD11_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD11_i0_i1_i2m2 - gammaDD11_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD11_i0_i1_i2m1 + gammaDD11_i0_i1_i2p1));
      const double gammaDD_dD120 = invdx0*(FDPart1_Rational_1_60*(-gammaDD12_i0m3_i1_i2 + gammaDD12_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD12_i0m2_i1_i2 - gammaDD12_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD12_i0m1_i1_i2 + gammaDD12_i0p1_i1_i2));
      const double gammaDD_dD121 = invdx1*(FDPart1_Rational_1_60*(-gammaDD12_i0_i1m3_i2 + gammaDD12_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD12_i0_i1m2_i2 - gammaDD12_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD12_i0_i1m1_i2 + gammaDD12_i0_i1p1_i2));
      const double gammaDD_dD122 = invdx2*(FDPart1_Rational_1_60*(-gammaDD12_i0_i1_i2m3 + gammaDD12_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD12_i0_i1_i2m2 - gammaDD12_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD12_i0_i1_i2m1 + gammaDD12_i0_i1_i2p1));
      const double gammaDD_dD220 = invdx0*(FDPart1_Rational_1_60*(-gammaDD22_i0m3_i1_i2 + gammaDD22_i0p3_i1_i2) + FDPart1_Rational_3_20*(gammaDD22_i0m2_i1_i2 - gammaDD22_i0p2_i1_i2) + FDPart1_Rational_3_4*(-gammaDD22_i0m1_i1_i2 + gammaDD22_i0p1_i1_i2));
      const double gammaDD_dD221 = invdx1*(FDPart1_Rational_1_60*(-gammaDD22_i0_i1m3_i2 + gammaDD22_i0_i1p3_i2) + FDPart1_Rational_3_20*(gammaDD22_i0_i1m2_i2 - gammaDD22_i0_i1p2_i2) + FDPart1_Rational_3_4*(-gammaDD22_i0_i1m1_i2 + gammaDD22_i0_i1p1_i2));
      const double gammaDD_dD222 = invdx2*(FDPart1_Rational_1_60*(-gammaDD22_i0_i1_i2m3 + gammaDD22_i0_i1_i2p3) + FDPart1_Rational_3_20*(gammaDD22_i0_i1_i2m2 - gammaDD22_i0_i1_i2p2) + FDPart1_Rational_3_4*(-gammaDD22_i0_i1_i2m1 + gammaDD22_i0_i1_i2p1));
      /*
       * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
       */
      const double FDPart35 = -gammaDD_dD001 + 2*gammaDD_dD010;
      const double FDPart37 = -gammaDD_dD002 + 2*gammaDD_dD020;
      const double FDPart39 = -gammaDD_dD012 + gammaDD_dD021 + gammaDD_dD120;
      const double FDPart310 = gammaDD_dD012 - gammaDD_dD021 + gammaDD_dD120;
      const double FDPart311 = -gammaDD_dD112 + 2*gammaDD_dD121;
      const double FDPart312 = 2*gammaDD_dD011 - gammaDD_dD110;
      const double FDPart313 = gammaDD_dD012 + gammaDD_dD021 - gammaDD_dD120;
      const double FDPart314 = 2*gammaDD_dD122 - gammaDD_dD221;
      const double FDPart315 = 2*gammaDD_dD022 - gammaDD_dD220;
      const double FDPart33 = (1.0/2.0)/(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - ((gammaDD02)*(gammaDD02))*gammaDD11);
      const double FDPart34 = FDPart33*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12)));
      const double FDPart36 = FDPart33*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12);
      const double FDPart38 = FDPart33*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11);
      const double FDPart316 = FDPart33*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02);
      const double FDPart317 = FDPart33*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02)));
      const double FDPart318 = FDPart33*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01)));
      aux_gfs[IDX4S(CHRISTOFFELUDD000GF, i0, i1, i2)] = FDPart34*gammaDD_dD000 + FDPart35*FDPart36 + FDPart37*FDPart38;
      aux_gfs[IDX4S(CHRISTOFFELUDD001GF, i0, i1, i2)] = FDPart34*gammaDD_dD001 + FDPart36*gammaDD_dD110 + FDPart38*FDPart39;
      aux_gfs[IDX4S(CHRISTOFFELUDD002GF, i0, i1, i2)] = FDPart310*FDPart36 + FDPart34*gammaDD_dD002 + FDPart38*gammaDD_dD220;
      aux_gfs[IDX4S(CHRISTOFFELUDD011GF, i0, i1, i2)] = FDPart311*FDPart38 + FDPart312*FDPart34 + FDPart36*gammaDD_dD111;
      aux_gfs[IDX4S(CHRISTOFFELUDD012GF, i0, i1, i2)] = FDPart313*FDPart34 + FDPart36*gammaDD_dD112 + FDPart38*gammaDD_dD221;
      aux_gfs[IDX4S(CHRISTOFFELUDD022GF, i0, i1, i2)] = FDPart314*FDPart36 + FDPart315*FDPart34 + FDPart38*gammaDD_dD222;
      aux_gfs[IDX4S(CHRISTOFFELUDD100GF, i0, i1, i2)] = FDPart316*FDPart37 + FDPart317*FDPart35 + FDPart36*gammaDD_dD000;
      aux_gfs[IDX4S(CHRISTOFFELUDD101GF, i0, i1, i2)] = FDPart316*FDPart39 + FDPart317*gammaDD_dD110 + FDPart36*gammaDD_dD001;
      aux_gfs[IDX4S(CHRISTOFFELUDD102GF, i0, i1, i2)] = FDPart310*FDPart317 + FDPart316*gammaDD_dD220 + FDPart36*gammaDD_dD002;
      aux_gfs[IDX4S(CHRISTOFFELUDD111GF, i0, i1, i2)] = FDPart311*FDPart316 + FDPart312*FDPart36 + FDPart317*gammaDD_dD111;
      aux_gfs[IDX4S(CHRISTOFFELUDD112GF, i0, i1, i2)] = FDPart313*FDPart36 + FDPart316*gammaDD_dD221 + FDPart317*gammaDD_dD112;
      aux_gfs[IDX4S(CHRISTOFFELUDD122GF, i0, i1, i2)] = FDPart314*FDPart317 + FDPart315*FDPart36 + FDPart316*gammaDD_dD222;
      aux_gfs[IDX4S(CHRISTOFFELUDD200GF, i0, i1, i2)] = FDPart316*FDPart35 + FDPart318*FDPart37 + FDPart38*gammaDD_dD000;
      aux_gfs[IDX4S(CHRISTOFFELUDD201GF, i0, i1, i2)] = FDPart316*gammaDD_dD110 + FDPart318*FDPart39 + FDPart38*gammaDD_dD001;
      aux_gfs[IDX4S(CHRISTOFFELUDD202GF, i0, i1, i2)] = FDPart310*FDPart316 + FDPart318*gammaDD_dD220 + FDPart38*gammaDD_dD002;
      aux_gfs[IDX4S(CHRISTOFFELUDD211GF, i0, i1, i2)] = FDPart311*FDPart318 + FDPart312*FDPart38 + FDPart316*gammaDD_dD111;
      aux_gfs[IDX4S(CHRISTOFFELUDD212GF, i0, i1, i2)] = FDPart313*FDPart38 + FDPart316*gammaDD_dD112 + FDPart318*gammaDD_dD221;
      aux_gfs[IDX4S(CHRISTOFFELUDD222GF, i0, i1, i2)] = FDPart314*FDPart316 + FDPart315*FDPart38 + FDPart318*gammaDD_dD222;
    }


<a id='what_next'></a>

# Part 5: What next? Navigating the NRPy+ tutorial \[Back to [top](#toc)\]
$$\label{what_next}$$

As mentioned previously, NRPy+ is meant to be "learn by example". To that end, there are more than 100 fully documented Jupyter notebooks covering a wide variety of topics relevant to numerical relativity research and optimized algorithms for solving PDEs numerically.

So the answer to the question "*What next?*" is, naturally, "*Whatever you like!*" To continue the journey, check out [**the main NRPy+ tutorial table of contents**.](NRPyPlus_Tutorial.ipynb).

<a id='latex_pdf_output'></a>

# Part 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-NRPyPlus_10_Minute_Overview.pdf](Tutorial-NRPyPlus_10_Minute_Overview.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-NRPyPlus_10_Minute_Overview")
```

    Traceback (most recent call last):
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/traitlets.py", line 632, in get
        value = obj._trait_values[self.name]
                ~~~~~~~~~~~~~~~~~^^^^^^^^^^^
    KeyError: 'template_paths'
    
    During handling of the above exception, another exception occurred:
    
    Traceback (most recent call last):
      File "/home/kid-a/Documents/projects/num-rel-hub/env/bin/jupyter-nbconvert", line 8, in <module>
        sys.exit(main())
                 ^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/jupyter_core/application.py", line 283, in launch_instance
        super().launch_instance(argv=argv, **kwargs)
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/config/application.py", line 1075, in launch_instance
        app.start()
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/nbconvertapp.py", line 420, in start
        self.convert_notebooks()
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/nbconvertapp.py", line 586, in convert_notebooks
        self.exporter = cls(config=self.config)
                        ^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/templateexporter.py", line 350, in __init__
        super().__init__(config=config, **kw)
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/exporter.py", line 123, in __init__
        self._init_preprocessors()
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/templateexporter.py", line 535, in _init_preprocessors
        conf = self._get_conf()
               ^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/templateexporter.py", line 553, in _get_conf
        for path in map(Path, self.template_paths):
                              ^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/traitlets.py", line 687, in __get__
        return t.cast(G, self.get(obj, cls))  # the G should encode the Optional
                         ^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/traitlets.py", line 635, in get
        default = obj.trait_defaults(self.name)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/traitlets.py", line 1897, in trait_defaults
        return t.cast(Sentinel, self._get_trait_default_generator(names[0])(self))
                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/traitlets/traitlets.py", line 1241, in __call__
        return self.func(*args, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/templateexporter.py", line 564, in _template_paths
        template_names = self.get_template_names()
                         ^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/kid-a/Documents/projects/num-rel-hub/env/lib/python3.12/site-packages/nbconvert/exporters/templateexporter.py", line 652, in get_template_names
        raise ValueError(msg)
    ValueError: No template sub-directory with name '../nbconvert_latex_settings' found in the following paths:
    	/home/kid-a/Documents/projects/num-rel-hub/env/share/jupyter
    	/home/kid-a/.local/share/jupyter
    	/usr/local/share/jupyter
    	/usr/share/jupyter
    Created Tutorial-NRPyPlus_10_Minute_Overview.tex, and compiled LaTeX file
        to PDF file Tutorial-NRPyPlus_10_Minute_Overview.pdf



```python

```


```python

```


```python

```
