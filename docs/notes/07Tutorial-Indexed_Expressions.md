<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Indexed Expressions: Representing and manipulating tensors, pseudotensors, etc. in NRPy+

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook demonstrates the use of NRPy+ in symbolic tensor computations and numerical code generation. It introduces techniques for manipulating indexed expressions with varying ranks, demonstrating practical applications such as tensor contractions and the conversion of C code with SIMD support.

**Notebook Status:** <font color='red'><b> Not Validated </b></font>

### NRPy+ Source Code for this module: [indexedexp.py](../edit/indexedexp.py)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#idx1): Rank-1 Indexed Expressions
    1. [Step 2.a](#dot): Performing a Dot Product
1. [Step 3](#idx2): Rank-2 and Higher Indexed Expressions 
    1. [Step 3.a](#con): Creating C Code for the contraction variable 
    1. [Step 3.b](#simd): Enable SIMD support
1. [Step 4](#exc): Exercise
1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from NRPy+ for dealing with indexed expressions and ouputting C code. 


```python
# The NRPy_param_funcs module sets up global structures that manage free parameters within NRPy+
import NRPy_param_funcs as par   # NRPy+: Parameter interface
# The indexedexp module defines various functions for defining and managing indexed quantities like tensors and pseudotensors
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
# The grid module defines various parameters related to a numerical grid or the dimensionality of indexed expressions
# For example, it declares the parameter DIM, which specifies the dimensionality of the indexed expression
import grid as gri               # NRPy+: Functions having to do with numerical grids
from outputC import outputC      # NRPy+: Basic C code output functionality
```

<a id='idx1'></a>

# Step 2: Rank-1 Indexed Expressions \[Back to [top](#toc)\]
$$\label{idx1}$$

Indexed expressions of rank 1 are stored as [Python lists](https://www.tutorialspoint.com/python/python_lists.htm). 

There are two standard ways to declare indexed expressions:
+ **Initialize indexed expression to zero:** 
    + **zerorank1(DIM=-1)** $\leftarrow$ As we will see below, initializing to zero is useful if the indexed expression depends entirely on some other indexed or non-indexed expressions.
        + **DIM** is an *optional* parameter that, if set to -1, will default to the dimension as set in the **grid** module: `par.parval_from_str("grid::DIM")`. Otherwise, the rank-1 indexed expression will have dimension **DIM**.
+ **Initialize indexed expression symbolically:** 
    + **declarerank1(symbol, DIM=-1)**. 
        + As in **`zerorank1()`**, **DIM** is an *optional* parameter that, if set to -1, will default to the dimension as set in the **grid** module: `par.parval_from_str("grid::DIM")`. Otherwise, the rank-1 indexed expression will have dimension **DIM**.

`zerorank1()` and `declarerank1()` are both wrapper functions for the more general function `declare_indexedexp()`.
+ **declare_indexedexp(rank, symbol=None, symmetry=None, dimension=None)**.
    + The following are optional parameters: **symbol**, **symmetry**, and **dimension**. If **symbol** is not specified, then `declare_indexedexp()` will initialize an indexed expression to zero. If **symmetry** is not specified or has value "nosym", then an indexed expression will not be symmetrized, which has no relevance for an indexed expression of rank 1. If **dimension** is not specified or has value -1, then **dimension** will default to the dimension as set in the **grid** module: `par.parval_from_str("grid::DIM")`.

For example, the 3-vector $\beta^i$ (upper index denotes contravariant) can be initialized to zero as follows:


```python
# Declare rank-1 contravariant ("U") vector
betaU = ixp.zerorank1()

# Print the result. It's a list of zeros!
print(betaU)
```

    [0, 0, 0]


Next set $\beta^i = \sum_{j=0}^i j = \{0,1,3\}$


```python
# Get the dimension we just set, so we know how many indices to loop over
DIM = par.parval_from_str("grid::DIM")

for i in range(DIM): # sum i from 0 to DIM-1, inclusive
    for j in range(i+1): # sum j from 0 to i, inclusive
        betaU[i] += j

print("The 3-vector betaU is now set to: "+str(betaU))
```

    The 3-vector betaU is now set to: [0, 1, 3]


Alternatively, the 3-vector $\beta^i$ can be initialized **symbolically** as follows:


```python
# Set the dimension to 3
par.set_parval_from_str("grid::DIM",3)

# Declare rank-1 contravariant ("U") vector
betaU = ixp.declarerank1("betaU")

# Print the result. It's a list!
print(betaU)
```

    [betaU0, betaU1, betaU2]


Declaring $\beta^i$ symbolically is standard in case `betaU0`, `betaU1`, and `betaU2` are defined elsewhere (e.g., read in from main memory as a gridfunction).

As can be seen, NRPy+'s standard naming convention for indexed rank-1 expressions is 
+ **\[base variable name\]+\["U" for contravariant (up index) or "D" for covariant (down index)\]**

*Caution*: after declaring the vector, `betaU0`, `betaU1`, and `betaU2` can only be accessed or manipulated through list access; i.e., via `betaU[0]`, `betaU[1]`, and `betaU[2]`, respectively. Attempts to access `betaU0` directly will fail.

Knowing this, let's multiply `betaU1` by 2:


```python
betaU[1] *= 2
print("The 3-vector betaU is now set to "+str(betaU))
print("The component betaU[1] is now set to "+str(betaU[1]))
```

    The 3-vector betaU is now set to [betaU0, 2*betaU1, betaU2]
    The component betaU[1] is now set to 2*betaU1


<a id='dot'></a>

## Step 2.a: Performing a Dot Product \[Back to [top](#toc)\]
$$\label{dot}$$

Next, let's declare the variable $\beta_j$ and perform the dot product $\beta^i \beta_i$:


```python
# First set betaU back to its initial value
betaU = ixp.declarerank1("betaU")

# Declare beta_j:
betaD = ixp.declarerank1("betaD")

# Get the dimension we just set, so we know how many indices to loop over
DIM = par.parval_from_str("grid::DIM")

# Initialize dot product to zero
dotprod = 0

# Perform dot product beta^i beta_i
for i in range(DIM):
    dotprod += betaU[i]*betaD[i]

# Print result!
print(dotprod)
```

    betaD0*betaU0 + betaD1*betaU1 + betaD2*betaU2


<a id='idx2'></a>

# Step 3: Rank-2 and Higher Indexed Expressions \[Back to [top](#toc)\]
$$\label{idx2}$$

Moving to higher ranks, rank-2 indexed expressions are stored as lists of lists, rank-3 indexed expressions as lists of lists of lists, etc. For example:

+ the covariant rank-2 tensor $g_{ij}$ is declared as `gDD[i][j]` in NRPy+, so that e.g., `gDD[0][2]` is stored with name `gDD02` and
+ the rank-2 tensor $T^{\mu}{}_{\nu}$ is declared as `TUD[m][n]` in NRPy+ (index names are of course arbitrary).

*Caveat*: note that it is currently up to the user to determine whether the combination of indexed expressions makes sense; NRPy+ does not track whether up and down indices are written consistently.

NRPy+ supports symmetries in indexed expressions (above rank 1), so that if $h_{ij} = h_{ji}$, then declaring `hDD[i][j]` to be symmetric in NRPy+ will result in both `hDD[0][2]` and `hDD[2][0]` mapping to the *single* SymPy variable `hDD02`.

To see how this works in NRPy+, let's define in NRPy+ a symmetric, rank-2 tensor $h_{ij}$ in three dimensions, and then compute the contraction, which should be given by $$con = h^{ij}h_{ij} = h_{00} h^{00} + h_{11} h^{11} + h_{22} h^{22} + 2 (h_{01} h^{01} + h_{02} h^{02} + h_{12} h^{12}).$$


```python
# Get the dimension we just set (should be set to 3).
DIM = par.parval_from_str("grid::DIM")

# Declare h_{ij}=hDD[i][j] and h^{ij}=hUU[i][j]
hUU = ixp.declarerank2("hUU","sym01")
hDD = ixp.declarerank2("hDD","sym01")

# Perform sum h^{ij} h_{ij}, initializing contraction result to zero
con = 0
for i in range(DIM):
    for j in range(DIM):
        con += hUU[i][j]*hDD[i][j]

# Print result
print(con)
```

    hDD00*hUU00 + 2*hDD01*hUU01 + 2*hDD02*hUU02 + hDD11*hUU11 + 2*hDD12*hUU12 + hDD22*hUU22


<a id='con'></a>

## Step 3.a: Creating C Code for the contraction variable $\text{con}$ \[Back to [top](#toc)\]
$$\label{con}$$

Next let's create the C code for the contraction variable $\text{con}$, without CSE (common subexpression elimination)


```python
outputC(con,"con")
```

    /*
     *  Original SymPy expression:
     *  "con = hDD00*hUU00 + 2*hDD01*hUU01 + 2*hDD02*hUU02 + hDD11*hUU11 + 2*hDD12*hUU12 + hDD22*hUU22"
     */
    {
      con = hDD00*hUU00 + 2*hDD01*hUU01 + 2*hDD02*hUU02 + hDD11*hUU11 + 2*hDD12*hUU12 + hDD22*hUU22;
    }
    


<a id='simd'></a>

## Step 3.b: Enable SIMD support \[Back to [top](#toc)\]
$$\label{simd}$$

Finally, let's see how it looks with SIMD support enabled


```python
outputC(con,"con",params="enable_SIMD=True")
```

    /*
     *  Original SymPy expression:
     *  "con = hDD00*hUU00 + 2*hDD01*hUU01 + 2*hDD02*hUU02 + hDD11*hUU11 + 2*hDD12*hUU12 + hDD22*hUU22"
     */
    {
      const double tmp_Integer_2 = 2.0;
      const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);
    
      con = FusedMulAddSIMD(hDD22, hUU22, FusedMulAddSIMD(_Integer_2, MulSIMD(hDD01, hUU01), FusedMulAddSIMD(_Integer_2, MulSIMD(hDD02, hUU02), FusedMulAddSIMD(_Integer_2, MulSIMD(hDD12, hUU12), FusedMulAddSIMD(hDD00, hUU00, MulSIMD(hDD11, hUU11))))));
    }
    


<a id='exc'></a>

# Step 4: Exercise \[Back to [top](#toc)\]
$$\label{exc}$$

Setting $\beta^i$ via the declarerank1(), write the NRPy+ code required to generate the needed C code for the lowering operator: $g_{ij} \beta^i$, and set the result to C variables `betaD0out`, `betaD1out`, and `betaD2out` [solution](Tutorial-Indexed_Expressions_soln.ipynb). *Hint: You will want to use the `zerorank1()` function*

**To complete this exercise, you must first reset all variables in the notebook:**


```python
# *Uncomment* the below %reset command and then press <Shift>+<Enter>.
#    Respond with "y" in the dialog box to reset all variables.

# %reset
```

**Write your solution below:**


```python

```

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Indexed_Expressions.pdf](Tutorial-Indexed_Expressions.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Indexed_Expressions")
```

    Created Tutorial-Indexed_Expressions.tex, and compiled LaTeX file to PDF
        file Tutorial-Indexed_Expressions.pdf

