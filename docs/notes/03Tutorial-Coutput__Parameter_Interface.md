<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# C Output and Parameter Interfaces

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

## Exploring C output and parameter interfaces in NRPy+, this notebook initializes core Python/NRPy+ modules, performs common subexpression elimination (CSE), and generates C code. It further delves into the NRPy+ parameter interface and demonstrates how Single Instruction, Multiple Data (SIMD) paradigms can optimize NRPy+ generated C code.

### Required reading if you are unfamiliar with programming or [computer algebra systems](https://en.wikipedia.org/wiki/Computer_algebra_system). Otherwise, use for reference; you should be able to pick up the syntax as you follow the tutorial.
+ **[Python Tutorial](https://docs.python.org/3/tutorial/index.html)**
+ **[SymPy Tutorial](http://docs.sympy.org/latest/tutorial/intro.html)**

### NRPy+ Source Code for this module:  
* [outputC.py](../edit/outputC.py)
* [NRPy_param_funcs.py](../edit/NRPy_param_funcs.py)
* [SIMD.py](../edit/SIMD.py)

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

The module is organized as follows:

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#sympy_ccode): Common Subexpression Elimination (CSE)
1. [Step 3](#coutput): **Let's generate some C code!** NRPy+'s core C code output routine, `Coutput()`
    1. [Step 3.a](#cfunction): **Wrap it up!** NRPy+'s C function wrapper routine, `outCfunction()`
1. [Step 4](#param): **Oh, the features you'll see!** Parameters in NRPy+
    1. [Step 4.a](#param_func): `NRPy_param_funcs`: The NRPy+ Parameter Interface
1. [Step 5](#simd): **Warp speed!** SIMD (Single Instruction, Multiple Data) in NRPy+-Generated C Code
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from Python/NRPy+ for dealing with parameter interfaces and outputting C code.   


```python
# Step 1: Initialize core Python/NRPy+ modules
from outputC import outputC,outCfunction # NRPy+: Core C code output module
import NRPy_param_funcs as par           # NRPy+: parameter interface
import sympy as sp                       # SymPy: The Python computer algebra package upon which NRPy+ depends
```

<a id='sympy_ccode'></a>

# Step 2: Common Subexpression Elimination (CSE) \[Back to [top](#toc)\]
$$\label{sympy_ccode}$$

Let's begin with a simple [SymPy](http://www.sympy.org/) worksheet that makes use of SymPy's built in C code generator function, [ccode](http://docs.sympy.org/dev/modules/utilities/codegen.html)(), to evaluate the expression $x = b^2 \sin (2a) + \frac{c}{\sin (2a)}$.


```python
# Step 2: Common Subexpression Elimination

# Declare some variables, using SymPy's symbols() function
a,b,c = sp.symbols("a b c")

# Set x = b^2*sin(2*a) + c/sin(2*a).
x = b**2*sp.sin(2*a) + c/(sp.sin(2*a))

# Convert the expression into C code
sp.ccode(x)
```




    'pow(b, 2)*sin(2*a) + c/sin(2*a)'



Computation of this expression in C requires 3 multiplications, one division, two sin() function calls, and one addition. Multiplications, additions, and subtractions typically require one clock cycle per SIMD element on a modern CPU, while divisions can require ~3x longer, and transcendental functions ~20x longer than additions or multiplications (See, e.g., [this page](https://software.intel.com/sites/landingpage/IntrinsicsGuide/#techs=AVX&expand=118), [this page](http://www.agner.org/optimize/microarchitecture.pdf), or [this page](http://nicolas.limare.net/pro/notes/2014/12/16_math_speed/) for more details). 

One goal in generating C codes involving mathematical expressions in NRPy+ is to minimize the number of floating point operations, and SymPy provides a means to do this, known as [common subexpression elimination](https://en.wikipedia.org/wiki/Common_subexpression_elimination), or CSE.

CSE algorithms search for common patterns within expressions and declare them as new variables, so they need not be computed again. To call SymPy's CSE algorithm, we need only pass the expression to [sp.cse()](http://docs.sympy.org/latest/modules/simplify/simplify.html#sympy.simplify.cse_main.cse):


```python
print(sp.cse(x))
```

    ([(x0, sin(2*a))], [b**2*x0 + c/x0])


As you can see, SymPy returned a list with two elements. The first element, $(\texttt{x0, sin(2*a)})$, indicates that a new variable $\texttt{x0}$ should be set to $\texttt{sin(2*a)}$. The second element yields the expression for our original expression $x$ in terms of the original variables, as well as the new variable $\texttt{x0}$. 

$$\texttt{x0} = \sin(2*a)$$ is the common subexpression, so that the final expression $x$ is given by $$x = pow(b,2)*\texttt{x0} + c/\texttt{x0}.$$

Thus, at the cost of a new variable assignment, SymPy's CSE has decreased the computational cost by one multiplication and one sin() function call.

NRPy+ makes full use of SymPy's CSE algorithm in generating optimized C codes, and in addition automatically adjusts expressions like `pow(x,2)` into `((x)*(x))`.

*Caveat: In order for a CSE to function optimally, it needs to know something about the cost of basic mathematical operations versus the cost of declaring a new variable. SymPy's CSE algorithm does not make any assumptions about cost, instead opting to declare new variables any time a common pattern is found more than once. The degree to which this is suboptimal is unclear.*

<a id='coutput'></a>

# Step 3: **Let's generate some C code!** NRPy+'s core C code output routine, `Coutput()` \[Back to [top](#toc)\]
$$\label{coutput}$$

NRPy+'s `outputC()` function provides the core of NRPy+ functionality. It builds upon SymPy's `ccode()` and `cse()` functions and adds the ability to generate [SIMD](https://en.wikipedia.org/wiki/SIMD) [compiler intrinsics](https://software.intel.com/sites/landingpage/IntrinsicsGuide/) for modern Intel and AMD-based CPUs. 

As `outputC()` is at the heart of NRPy+, it will be useful to understand how it is called:

```python
outputC(sympyexpr, output_varname_str, filename = "stdout", params = "", prestring = "", poststring = "")
```

`outputC()` requires at least two arguments: 
+ **sympyexpr** is a SymPy expression or a list of SymPy expressions
+ **output_varname_str** is the variable name to assign the SymPy expression, or alternatively the list of variable names to assign the SymPy expressions. If a list is provided, it must be the same length as the list of SymPy expressions.

Additional, optional arguments to `outputC()` include
+ **filename** (third argument; defaults to "stdout" if unspecified): 
     + "stdout" = print to the screen
     + "filename.c" = output to filename.c
     + "returnstring" = return C output as a string. I.e., call 
         + string = outputC(sympyexpr, output_varname_str, filename = "returnstring")
         + ... and then manipulate the string directly.
+ **params** (fourth argument; defaults to "" if unspecified): A comma-separated list of tunable parameters. For example: *params="preindent=1,includebraces=False,declareoutputvars=False,SIMD_debug=True"* Parameters can be listed in any order, and repeats are allowed; the final repeated value will be the value that is set. List of parameters:
    + *preindent*: (integer, defaults to 0) The number of tab stops to add to C code output
    + *includebraces*: (True or False, defaults to True) Wrap the C output expression in curly braces?
    + *declareoutputvars*: (True or False, defaults to False) Prepend the output variable with the variable type, thus defining the output variable.
    + *outCfileaccess*: ("w" or "a", defaults to "w") Write ("w") or append ("a") to the C output file.
    + *outCverbose*: (True or False, defaults to True) Output a comment block displaying the input SymPy expressions, if set to True.
    + *CSE_enable*: (True or False, defaults to True) If set to True, common-subexpression elimination (CSE) will be used.
    + *CSE_varprefix*: (Any string without spaces, defaults to "tmp") Prefix each temporary variable in the CSE with this string.
    + *enable_SIMD*: (True or False, defaults to False) If set to True, C code output exclusively uses SIMD compiler intrinsics, which must be linked to the actual intrinsics for a given compiler/SIMD library through C macros.
    + *SIMD_debug*: (True or False, defaults to False) Verify for each expression that the SIMD output matches the input (non-SIMD) expression.
+ **prestring** (fifth argument; defaults to "" -- empty string -- if unspecified): Preface C code output with prestring.
+ **poststring** (sixth argument): Same as prestring, but places poststring at the end of the C code output.

Notice that by default, CSE is enabled (fourth function argument). Thus if we call outputC with two arguments, NRPy+ will process the expression through SymPy's CSE:


```python
# Step 3: NRPy+'s C code output routine, `Coutput()`

# Declare some variables, using SymPy's symbols() function
a,b,c = sp.symbols("a b c")

# Set x = b^2*sin(2*a) + c/sin(2*a).
x = b**2*sp.sin(2*a) + c/(sp.sin(2*a))

outputC(x,"x")
```

    /*
     *  Original SymPy expression:
     *  "x = b**2*sin(2*a) + c/sin(2*a)"
     */
    {
      const double tmp0 = sin(2*a);
      x = ((b)*(b))*tmp0 + c/tmp0;
    }
    


<a id='cfunction'></a>

## Step 3.a: **Wrap it up!** NRPy+'s C function wrapper routine, `outCfunction()` \[Back to [top](#toc)\]
$$\label{cfunction}$$


```python
# Declare some variables, using SymPy's symbols() function
a,b,c = sp.symbols("a b c")

# Set x = b^2*sin(2*a) + c/sin(2*a).
x = b**2*sp.sin(2*a) + c/(sp.sin(2*a))

desc="Output x(a,b,c) = b^2*sin(2*a) + c/sin(2*a)"
name="output_x_of_a_b_c"
string = outCfunction(
    outfile  = "returnstring", desc=desc, name=name,
    params   = "const double a,const double b,const double c, double *x",
    body     = outputC(x,"*x",filename="returnstring",params="includebraces=False,preindent=1"),
    enableCparameters=False)
print(string)
```

    /*
     * Output x(a,b,c) = b^2*sin(2*a) + c/sin(2*a)
     */
    void output_x_of_a_b_c(const double a,const double b,const double c, double *x) {
    
      /*
       *  Original SymPy expression:
       *  "*x = b**2*sin(2*a) + c/sin(2*a)"
       */
      const double tmp0 = sin(2*a);
      *x = ((b)*(b))*tmp0 + c/tmp0;
    }
    


## <a id='param'></a>

# Step 4: Oh, the features you'll see! Parameters in NRPy+ \[Back to [top](#toc)\]
$$\label{param}$$

*TL;DR: When adding new features to NRPy+ or to modules that use NRPy+, it is strongly recommended to take advantage of NRPy+'s parameter interface.*

As documented above, NRPy+'s `outputC()` routine accepts up to six inputs. Suppose we have a project that makes use of NRPy+ to generate *multiple C codes* for a project. It is reasonable to expect that these six inputs might vary from one C code to the next in the same project. (For example, sometimes a C code will be sufficiently simple that CSE only acts to obfuscate.) Thus we include these six inputs as part of the function call.

Suppose we wanted to add another feature to `outputC()` that is universal to our project. If `outputC()`'s behavior were only steerable with inputs into the function call, then the number of inputs will balloon with the number of features, making the entire NRPy+ codebase far less manageable. To address this problem, while at the same time making the modules and functions within NRPy+ more easily extensible, we have introduced a parameter interface.

<a id='param_func'></a>

## Step 4.a: NRPy_param_funcs: The NRPy+ Parameter Interface \[Back to [top](#toc)\]
$$\label{param_func}$$

The **`NRPy_param_funcs`** module manages the parameter interface in NRPy+, and parameter information is stored in two global data structures defined within this module:

* glb_params_list\[\]: The list of registered parameters. Each item in the list is a [named tuple](https://docs.python.org/2/library/collections.html#collections.namedtuple) of the type `glb_param`, where the type is defined
    * glb_param(\[parameter type\],\[module name\],\[parameter name\],\[default value\]).
* glb_paramsvals_list\[\]: The list of parameter values. When a new glb_param is appended to the glb_params_list\[\], the corresponding element in glb_paramsvals_list\[\] is set to the default value. This value can be overwritten with
    * parameter files or parameter file overrides *when running NRPy+ in command-line mode*, as follows:
        * **python nrpy.py \[PARAMETER FILE\] \[PARAMETER FILE OVERRIDES\]**), or
    * with set_paramsvals_value("modulename::variablename = \[value\]") *when running NRPy+ in interactive mode*.

**Example**: Suppose you write a new module, *mymodule* (in "mymodule.py") that depends on NRPy+, which contains a free parameter $n$, which is an integer (we set integers in NRPy+ as *type="int"*). A reasonable default value is $n=2$. To register this parameter with NRPy+, set the following at the top of your mymodule.py:


```python
# Step 4.a: NRPy_param_funcs: The NRPy+ Parameter Interface

par.initialize_param(par.glb_param(type="int", module="mymodule", parname="n", defaultval=2))
```

At any time, you can find the parameter's value via the `par.parval_from_str()` function, which accepts a string in one of two formats: "`variablename`" or "`modulename::variablename`". *Warning*: If more than one module sets the parameter with variable name `"n"`, `par.parval_from_str("n")` will produce an error.


```python
print(par.parval_from_str("n"))
print(par.parval_from_str("mymodule::n"))
```

    2
    2


Next, let's overwrite the default parameter value of `"mymodule::n"` to be 4 instead:


```python
par.set_paramsvals_value("mymodule::n = 4")
print(par.parval_from_str("mymodule::n"))
```

    4


**Warning**: Setting NRPy+ parameters via direct calls to `par.set_paramsvals_value("modulename::variablename")`, when in non-interactive mode (e.g., running in a Jupyter or iPython notebook) is *strongly* discouraged, and in the future may result in an error message.

<a id='simd'></a>

# Step 5: Warp speed! SIMD (Single Instruction, Multiple Data) in NRPy+-Generated C Code \[Back to [top](#toc)\]
$$\label{simd}$$

Taking advantage of a CPU's SIMD instruction set can yield very nice performance boosts, but only when the CPU can be used to process a large data set that can be performed in parallel. It enables the computation of multiple parts of the data set at once. 

For example, given the expression 
$$\texttt{double x = a*b},$$ 
where $\texttt{double}$ precision variables $\texttt{a}$ and $\texttt{b}$ vary at each point on a computational grid, AVX compiler intrinsics will enable the multiplication computation at *four* grid points *each clock cycle*, *on each CPU core*. Therefore, without these intrinsic, the computation might take four times longer. Compilers can sometimes be smart enough to "vectorize" the loops over data, but when the mathematical expressions become too complex (e.g., in the context of numerically solving Einstein's equations of general relativity), the compiler will simply give up and refuse to enable SIMD vectorization.

As SIMD intrinsics can differ from one CPU to another, and even between compilers, NRPy+ outputs generic C macros for common arithmetic operations and transcendental functions. In this way, the C code's Makefile can decide the most optimal SIMD intrinsics for the given CPU's instruction set and compiler. For example, most modern CPUs support [AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions), and a majority support up to [AVX2](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2), while some support up to [AVX512](https://en.wikipedia.org/wiki/AVX-512) instruction sets. For a full list of compiler intrinsics, see the [official Intel SIMD intrinsics documentation](https://software.intel.com/sites/landingpage/IntrinsicsGuide/).

To see how this works, let's return to our NRPy+ `outputC()` CSE example above, but this time enabling SIMD intrinsics:


```python
# Step 5: Taking Advantage of SIMD (Single Instruction, Multiple Data) in NRPy+-Generated C Code

# Declare some variables, using SymPy's symbols() function
a,b,c = sp.symbols("a b c")

# Set x = b^2*sin(2*a) + c/sin(2*a).
x = b**2*sp.sin(2*a) + c/(sp.sin(2*a))

outputC(x,"x",params="enable_SIMD=True")
```

    /*
     *  Original SymPy expression:
     *  "x = b**2*sin(2*a) + c/sin(2*a)"
     */
    {
      const double tmp_Integer_1 = 1.0;
      const REAL_SIMD_ARRAY _Integer_1 = ConstSIMD(tmp_Integer_1);
    
      const double tmp_Integer_2 = 2.0;
      const REAL_SIMD_ARRAY _Integer_2 = ConstSIMD(tmp_Integer_2);
    
      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);
    
      const REAL_SIMD_ARRAY tmp_0 = SinSIMD(MulSIMD(_Integer_2, a));
      x = FusedMulAddSIMD(tmp_0, MulSIMD(b, b), DivSIMD(c, tmp_0));
    }
    


The above SIMD code does the following.
* First it fills a constant SIMD array of type `REAL_SIMD_ARRAY `with the integer 2 to the double-precision 2.0. The larger C code in which the above-generated code will be embedded should automatically `#define REAL_SIMD_ARRAY` to e.g., _m256d or _m512d for AVX or AVX512, respectively. In other words, AVX intrinsics will need to set 4 double-precision variables in `REAL_SIMD_ARRAY` to 2.0, and AVX-512 intrinsics will need to set 8.
* Then it changes all arithmetic operations to be in the form of SIMD "functions", which are in fact #define'd in the larger C code as compiler intrinsics. 

FusedMulAddSIMD(a,b,c) performs a fused-multiply-add operation (i.e., `FusedMulAddSIMD(a,b,c)`=$a*b+c$), which can be performed on many CPUs nowadays (with FMA or AVX-512 instruction support) with a *single clock cycle*, at nearly the same expense as a single addition or multiplication, 

Note that it is assumed that the SIMD code exists within a suitable set of nested loops, in which the innermost loop increments every 4 in the case of AVX double precision or 8 in the case of AVX-512 double precision.

As an additional note, NRPy+'s SIMD routines are aware that the C `pow(x,y)` function is exceedingly expensive when $|\texttt{y}|$ is a small integer. It will automatically convert such expressions into either multiplications of x or one-over multiplications of x, as follows (notice there are no calls to `PowSIMD()` intrinsics!):


```python
# Declare some variables, using SymPy's symbols() function
a,b,c = sp.symbols("a b c")

# Set x = b^2*sin(2*a) + c/sin(2*a).
x = b**2 + a**(-3) + c*a**(sp.Rational(1,2))

outputC(x,"x", params="enable_SIMD=True")
```

    /*
     *  Original SymPy expression:
     *  "x = sqrt(a)*c + b**2 + a**(-3)"
     */
    {
      const double tmp_Integer_1 = 1.0;
      const REAL_SIMD_ARRAY _Integer_1 = ConstSIMD(tmp_Integer_1);
    
      x = FusedMulAddSIMD(b, b, FusedMulAddSIMD(c, SqrtSIMD(a), DivSIMD(_Integer_1, MulSIMD(MulSIMD(a, a), a))));
    }
    


For those who would like to maximize fused-multiply-adds (FMAs) and fused-multiply-subtracts (FMSs), NRPy+ has more advanced pattern matching, which can be enabled via the `params="SIMD_find_more_FMAsFMSs=True"` option. **Note that finding more FMAs and FMSs may actually degrade performance, and the default behavior is found to be optimal on x86_64 CPUs.** In the below example, notice that the more advanced pattern matching finds another FMA:


```python
print("// SIMD_find_more_FMAsFMSs=True:\n// searches for more FMAs/FMSs, which has been found to degrade performance on some CPUs:")
outputC(x,"x", params="enable_SIMD=True,SIMD_find_more_FMAsFMSs=True")
```

    // SIMD_find_more_FMAsFMSs=True:
    // searches for more FMAs/FMSs, which has been found to degrade performance on some CPUs:
    /*
     *  Original SymPy expression:
     *  "x = sqrt(a)*c + b**2 + a**(-3)"
     */
    {
      const double tmp_Integer_1 = 1.0;
      const REAL_SIMD_ARRAY _Integer_1 = ConstSIMD(tmp_Integer_1);
    
      x = FusedMulAddSIMD(b, b, FusedMulAddSIMD(c, SqrtSIMD(a), DivSIMD(_Integer_1, MulSIMD(MulSIMD(a, a), a))));
    }
    


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Coutput__Parameter_Interface.pdf](Tutorial-Coutput__Parameter_Interface.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Coutput__Parameter_Interface")
```

    Created Tutorial-Coutput__Parameter_Interface.tex, and compiled LaTeX file
        to PDF file Tutorial-Coutput__Parameter_Interface.pdf
