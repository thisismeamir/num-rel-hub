# Overview

NRPy+ provides a tool to make efficient (in terms of computational cycles) C code that can be optimized based on common CPU architectures. Another thing that NRPy+ suggests is that although compilers would try their best to optimize codes, sometimes in complex expressions (like the ones that we are going to see in Numerical Relativity), they just give up! NRPy+ tries to also optimize in those situations (obviously because the package's whole meaning is in those situations).

# CSE
NRPy+ uses sympy (pythons symbolic computation module) to enhance its C code production. 

The core idea is to use `sp.cse(expr)` to find repeating expressions (or any means to make the expression computationally efficient) and then replace the expression with a new one. as an example: 

$$
\sin(x +2)+ x +2
$$
would turn into, a one time computation of x+2 and then substituting it into the following expression

$$
\sin(\alpha)+\alpha
$$
Where $\alpha$ is the expression $x+2$ but computed. For this example we just used one addition instead of two for computing $x+2$ both inside the sine and as the other expression.

# NRPy+ `outputC()`

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

Notice that by default, CSE is enabled (fourth function argument). Thus if we call outputC with two arguments, NRPy+ will process the expression through SymPy's CSE.

# NRPy+'s C function wrapper routine, `outCfunction()`

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

# NRPy_param_funcs: The NRPy+ Parameter Interface 


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