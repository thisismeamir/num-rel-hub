<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Loop Generation and Cache Blocking

## Author: Ken Sible

## This notebook presents the generation of C-based loops with arbitrary dimensions and cache blocking, using NRPy+. The `loop.py` module produces single or nested loops with specified options and supports cache blocking to improve spatial locality, as demonstrated using examples of matrix-vector multiplication and others, thereby enhancing computational efficiency.
[comment]: <> (**Notebook Status:** <font color='red'><b> Not Validated </b></font>)

### NRPy+ Source Code for this module:
1. [loop.py](../edit/loop.py); [\[**tutorial**\]](Tutorial-Loop_Generation_Cache_Blocking.ipynb) The loop.py script will generate a single or nested loop of arbitrary dimension in C, and has support for [cache blocking](https://software.intel.com/en-us/articles/how-to-use-loop-blocking-to-optimize-memory-use-on-32-bit-intel-architecture).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

1. [Step 1](#generation): Demonstrate Loop Generation
    1. [Step 1.a](#simple): Generate a simple, three-dimensional loop with specified options
    1. [Step 1.b](#arbitrary): Generate a nested loop of arbitrary dimension with exact parameters
1. [Step 2](#exercise): Exercise (Loop Generation)
1. [Step 3](#tiling): Demonstrate Cache Blocking
1. [Step 4](#latex_pdf_output): $\LaTeX$ PDF Output

<a id='generation'></a>

# Step 1: Demonstrate Loop Generation \[Back to [top](#toc)\]
$$\label{generation}$$

In the following section, we demonstrate single and nested loop generation in C using NRPy+.

In [loop.py](../edit/loop.py), the following functions are implemented for loop generation:

- ```loop(idx_var, lower_bound, upper_bound, increment, pragma, padding="  ", interior="", tile_size="")```
    + ```idx_var```: index variable for the loop (```idxvar[0]```$\Rightarrow$ outermost loop, ```idxvar[N - 1]```$\Rightarrow$ innermost loop)
    + ```lower_bound```:&nbsp;&nbsp;&nbsp;lower bound for ```idxvar``` or ```idxvar[i]```
    + ```upper_bound```:&nbsp;&nbsp;&nbsp;upper bound for ```idxvar``` or ```idxvar[i]```
    + ```increment```:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;increment for ```idxvar``` or ```idxvar[i]```
    + ```pragma```: OpenMP pragma for ```idxvar``` or ```idxvar[i]```
    + ```padding```: (*optional*) padding before a line (tab stop)
    + ```interior```: &nbsp;&nbsp;(*optional*) interior of the loop
    + ```tile_size```: &nbsp;&nbsp;(*optional*) tile size for cache blocking
- ```simple_loop(options, interior)```
    + ```options```: options for loop generation
    + ```interior```: interior of the loop


```python
from loop import loop, simple_loop # Import NRPy+ module for loop generation
```

<a id='simple'></a>

## Step 1.a:  `simple_loop()` \[Back to [top](#toc)\]
$$\label{simple}$$

The `simple_loop()` function will generate a simple loop in C (for use inside of a function) with specified options.


```python
# 'AllPoints': loop over all points on a numerical grid, including ghost zones
print(simple_loop('AllPoints', '// <INTERIOR>'))
```

      #pragma omp parallel for
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
            // <INTERIOR>
          } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
      } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    



```python
# 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
print(simple_loop('InteriorPoints', '// <INTERIOR>'))
```

      #pragma omp parallel for
      for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
        for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
          for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
            // <INTERIOR>
          } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
        } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
      } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    



```python
# 'Read_xxs': read the xx[3][:] 1D coordinate arrays, as some interior dependency exists
print(simple_loop('AllPoints Read_xxs', '// <INTERIOR>'))
#HERE
```

      #pragma omp parallel for
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          const REAL xx1 = xx[1][i1];
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
            const REAL xx0 = xx[0][i0];
            // <INTERIOR>
          } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
      } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    



```python
# 'enable_rfm_precompute': enable pre-computation of reference metric
print(simple_loop('AllPoints enable_rfm_precompute', '// <INTERIOR>'))
```

      #pragma omp parallel for
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        #include "rfm_files/rfm_struct__read2.h"
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          #include "rfm_files/rfm_struct__read1.h"
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
            #include "rfm_files/rfm_struct__read0.h"
            // <INTERIOR>
          } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
      } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    



```python
# 'enable_SIMD': enable SIMD support (https://en.wikipedia.org/wiki/SIMD)
print(simple_loop('AllPoints enable_rfm_precompute enable_SIMD', '// <INTERIOR>'))
```

      #pragma omp parallel for
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        #include "rfm_files/rfm_struct__SIMD_outer_read2.h"
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          #include "rfm_files/rfm_struct__SIMD_outer_read1.h"
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0 += SIMD_width) {
            #include "rfm_files/rfm_struct__SIMD_inner_read0.h"
            // <INTERIOR>
          } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0 += SIMD_width)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
      } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    



```python
# 'DisableOpenMP': disable OpenMP parallelization (https://en.wikipedia.org/wiki/OpenMP)
print(simple_loop('AllPoints DisableOpenMP', '// <INTERIOR>'))
```

      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
            // <INTERIOR>
          } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
      } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    


<a id='arbitrary'></a>

## Step 1.b:  `loop()` \[Back to [top](#toc)\]
$$\label{arbitrary}$$

The `loop()` function will generate a nested loop of arbitrary dimension in C with exact parameters (i.e. more control than `simple_loop()`).


```python
# Generate a one-dimensional loop over i from 0 to N, stepping by 1, with specified loop body
print(loop('i', '0', 'N', '1', '', interior='// <INTERIOR>'))
```

      for (int i = 0; i < N; i++) {
        // <INTERIOR>
      } // END LOOP: for (int i = 0; i < N; i++)
    



```python
# Generate a two-dimensional loop over i from 0 to Nx, stepping by 1,
#   and j from 0 to Ny, stepping by 1, with specified loop body
print(loop(['i', 'j'], ['0', '0'], ['Nx', 'Ny'], ['1', '1'], ['', ''], interior='// <INTERIOR>'))
```

      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          // <INTERIOR>
        } // END LOOP: for (int j = 0; j < Ny; j++)
      } // END LOOP: for (int i = 0; i < Nx; i++)
    


<a id='exercise'></a>

# Step 2: Exercise (Loop Generation) \[Back to [top](#toc)\]
$$\label{exercise}$$

**Goal:** Reproduce the following output using the loop generation infrastructure in NRPy+ ([solution](Tutorial-Loop_Generation_Cache_Blocking_soln.ipynb))

```
for (int n = 0; n < (Nt - 1); n++) {
    u[n][0] = u[n][Nx] = 0;
    for (int k = 1; k < (Nx - 1); k++) {
        u[n + 1][k] = u[n][k] + r*(u[n][k + 1] - 2*u[n][k] + u[n][k - 1]);
    } // END LOOP: for (int k = 1; k < (Nx - 1); k++)
    for (int k = 0; k < Nx; k++) {
        u[n][k] = u[n + 1][k];
    } // END LOOP: for (int k = 0; k < Nx; k++)
} // END LOOP: for (int n = 0; n < (Nt - 1); n++)
```


```python
# Write Solution Here
```

<a id='tiling'></a>

# Step 3: Demonstrate Cache Blocking \[Back to [top](#toc)\]
$$\label{tiling}$$

In the following section, we demonstrate cache blocking (loop tiling) using NRPy+. The advantage of cache blocking is improved [spatial locality](https://en.wikipedia.org/wiki/Locality_of_reference) by blocking or tiling a data structure to fit inside the cache. We minimize the number of cache misses that occur (reduce memory bandwidth pressure) by reusing the subset of our data structure that was cached rather than accessing main memory ([source](https://software.intel.com/en-us/articles/cache-blocking-techniques)). The following example of matrix-vector multiplication will demonstrate cache blocking using NRPy+.


```python
print('// Untiled Loop')
print(loop(['i', 'j'], ['0', '0'], ['N', 'N'], ['1', '1'], ['', ''], interior='c[i] += a[i][j] * b[j];'))

print('// Tiled Loop\n#define MIN(x, y) (((x) < (y)) ? (x) : (y))')
print(loop(['i', 'j'], ['0', '0'], ['N', 'N'], ['1', '1'], ['', ''], interior='c[i] += a[i][j] * b[j];', tile_size=['2', '2']))
```

    // Untiled Loop
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          c[i] += a[i][j] * b[j];
        } // END LOOP: for (int j = 0; j < N; j++)
      } // END LOOP: for (int i = 0; i < N; i++)
    
    // Tiled Loop
    #define MIN(x, y) (((x) < (y)) ? (x) : (y))
      for (int iB = 0; iB < N; iB += 2) {
        for (int jB = 0; jB < N; jB += 2) {
          for (int i = iB; i < MIN(N, iB + 2); i++) {
            for (int j = jB; j < MIN(N, jB + 2); j++) {
              c[i] += a[i][j] * b[j];
            } // END LOOP: for (int j = jB; j < MIN(N, jB + 2); j++)
          } // END LOOP: for (int i = iB; i < MIN(N, iB + 2); i++)
        } // END LOOP: for (int jB = 0; jB < N; jB += 2)
      } // END LOOP: for (int iB = 0; iB < N; iB += 2)
    


Order of Memory Access (Untiled): `(a[0][0], b[0]), (a[0][1], b[1]), (a[0][2], b[2]), (a[0][3], b[3]), ...`

Order of Memory Access (Tiled): &nbsp;&nbsp;&nbsp;&nbsp;`(a[0][0], b[0]), (a[0][1], b[1]), (a[1][0], b[0]), (a[1][1], b[1]), ...`

We should remark that cache blocking might not be a valid loop optimization whenever the order of memory access affects the resulting output. However, that potential issue does not occur in our example with matrix-vector multiplication. Moreover, the block or tile size will depend on the CPU architecture, and hence experimentation is required to determine the optimal size.

<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Loop_Generation_Cache_Blocking.pdf](Tutorial-Loop_Generation_Cache_Blocking.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Loop_Generation_Cache_Blocking")
```

    Created Tutorial-Loop_Generation_Cache_Blocking.tex, and compiled LaTeX
        file to PDF file Tutorial-Loop_Generation_Cache_Blocking.pdf

