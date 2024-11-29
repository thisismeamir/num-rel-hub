# Overview
In this document we would examine a python package called SERN/NRPy+, for numerical relativity and beyond. This document is the basics of this package which is going through its first notebook tutorial (NRPy+ 10-Minute Overview)

# Introduction
Einstein notation certainly makes relativistic equations look much simpler than before, NRPy+ assumes that using this notation, each geometrical object (scalar, vector, tensor, ...) is a code. For NRPy+:
1. A variable represents a rank-0 tensor,
2. A list represents a rank-1 tensor (vector),
3. A list of lists represents a rank-2 tensor,
4. A list of lists of lists ...
NRPy+ combines this idea with SymPy a python library for symbolical computation (computer algebra system), but this is not the end. In fact, NRPy+ has a pipeline to transform these complex symbolical expressions into robust C/C++ code.

# Using NRPy+
## Constructing $\!^{(3)}\Gamma^i_{jk}$ 
We began our journey in NRPy+ by constructing the $\!^{(3)}\Gamma^i_{jk}$ as symbolic expressions in terms of 3-metric $\gamma_{ij}$ and its derivatives. 

$$
\begin{align}  
\Gamma_{ij}^k &= \frac{1}{2} \gamma^{kl}\left(\gamma_{jl,i} + \gamma_{il,j} - \gamma_{ij,l}\right)\\  
\end{align}
$$
The above equation is how a normal theoretical relativist would try to construct the values. For NRPy+ the task becomes:

```python
import sympy as sp        # SymPy: The Python computer algebra package upon which NRPy+ depends  
import indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support  
```

After importing we would introduce our metric:

```python
# gammaDD is a rank-2 indexed expression symmetric in indexes 0 and 1  
gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01", DIM=3)  
```

note that we have introduced a symbol `"gammaDD"`, symmetry `"sym01"` and the dimension `DIM=3`, and now we would calculate $\gamma^{ij}$ as well as its determinant $\gamma$: 

```python
# gammaUU is the inverse of gammaDD, and  
# gammaDET (unused) is the determinant of gammaDD  
gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)  # inverts a 3x3 symmetric matrix
```

Now we can start defining the Christoffel symbols:

```python
# First define symbolic expressions for metric derivatives  
gammaDD_dD = ixp.declarerank3("gammaDD_dD", symmetry="sym01", DIM=3)  
```

this introduces an abstract symbolic expression for the derivative of metric. It is abstract because we have not yet introduced the components of neither this nor the metric.

```python
# Initialize GammaUDD (3-Christoffel) to zero  
GammaUDD = ixp.zerorank3(DIM=3)  
for i in range(3):  
    for j in range(3):  
        for k in range(3):  
            for l in range(3):  
                GammaUDD[k][i][j] += sp.Rational(1, 2) * gammaUU[k][l] * (  
                    gammaDD_dD[j][l][i] + gammaDD_dD[i][l][j] - gammaDD_dD[i][j][l])
```

Here we initialized $\Gamma$ as a rank-three tensor with dimensions of three as well. Then we iterated through each component and set each as the expression (the Christoffel definition told us).

As we know the lower indices of Christoffel symbols are symmetric, which means:

$$
\Gamma^i_{jk} = \Gamma^i_{kj}
$$
we can check that out here as well:

```python
for i in range(3):  
    for j in range(3):  
        for k in range(3):  
            print(GammaUDD[i][j][k] - GammaUDD[i][k][j])
```

which would result zero for all the combinations. Now we can convert this code into a highly optimized C code using:
1. Storing all the symbols and variables:
2. Optimizing

```python
symbols_list = []  
varname_list = []  
for i in range(3):  
    for j in range(3):  
        for k in range(j, 3):  
            symbols_list += [GammaUDD[i][j][k]]  
            varname_list += ["ChristoffelUDD" + str(i) + str(j) + str(k)]
```

then

```python
import outputC as outC  # NRPy+: Core C code output module  
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=False,outCverbose=False")
```

Notice in the above code that many expressions are recomputed time and time again. This is *incredibly inefficient*, and generally compilers won't optimize this properly. So in optimization level 1, we use common subexpression elimination.  
  
**Optimization Level 1**: Use [common subexpression elimination (CSE)](https://en.wikipedia.org/wiki/Common_subexpression_elimination) to group common subexpressions

```python
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=True,outCverbose=False")
```

**Optimization Level 2**: Use CSE and take advantage of [single instruction, multiple data (SIMD) macros](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data). NRPy+ translates these macros into assembler-level-optimized code for (x86_64) CPUs. Use case: looping over data on a numerical grid; can evaluate expressions at up to 8 gridpoints simultaneously *per CPU core*.

```python
outC.outputC(symbols_list, varname_list, filename="stdout", params="CSE_enable=True,enable_SIMD=True,outCverbose=False")
```
