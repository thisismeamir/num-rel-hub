<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Adding Unit Tests to the NRPy+ Unit Testing Infrastructure

## Author: Kevin Lituchy

## This module outlines the process of establishing and interpreting unit tests within NRPy+'s framework, facilitating rigorous checks to prevent unnoticed discrepancies and ensuring code reliability across different versions.

## Introduction:
The goal of this module is to give the user an overview/understanding of NRPy+'s Unit Testing framework, which will give the user enough information to begin creating their own unit tests. We will begin by giving an overview of the important prerequisite knowledge to make the most out of unit tests. Next, we give an explanation for the user interaction within the unit testing framework; this will give the user the ability to create tests for themselves. Then we give the user some insight into interpreting the output from their unit tests. Finally, a full example using a test module will be run through in full, both with and without errors, to give the user a realistic sense of what unit testing entails.

For in-depth explanations of all subfunctions (not user-interactable), see the [UnitTesting Function Reference](./UnitTesting/UnitTesting_Function_Reference.ipynb). This may not be essential to get unit tests up and running, but it will be invaluable information if the user ever wants to make modifications to the unit testing code.


<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 1](#motivation): Motivation and Prerequisite Knowledge
    1. [Step 1.a](#dicts): Dictionaries
    1. [Step 1.b](#logging): Logging
1. [Step 2](#interaction): User Interaction
    1. [Step 2.a](#testfile): Test File
    1. [Step 2.b](#trustedvaluesdict): trusted_values_dict
    1. [Step 2.c](#bash): Bash Script
1. [Step 3](#output): Interpreting output
1. [Step 4](#checklist): Checklist
    1. [Step 4.a](#4a): Directory Creation
    1. [Step 4.b](#4b): Test File Creation
    1. [Step 4.c](#4c): Input Parameters
    1. [Step 4.d](#4d): Fill-in bash script and run

1. [Step 5](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF
    


<a id='motivation'></a>

# Step 1: Motivation and Prerequisite Knowledge \[Back to [top](#toc)\]
$$\label{motivation}$$

What is the purpose of unit testing, and why should you do it? To begin
thinking about that, consider what subtleties can occur within your code
that is almost unnoticeable to the eye, but wind up giving you an
incorrect result. You could make a small optimization, and not notice
any change in your result. However, maybe the optimization you made only
works on Python 3 and not Python 2, or it changes a value by some tiny
amount--too small to be noticeable at a simple glance, but enough to
make a difference in succeeding calculations.

This is where unit testing comes in. By initially calculating values for
the globals of your modules in a **trusted** version of your code and
storing those values in a dictionary, you can then easily check if
something stopped working correctly by comparing your newly calculated
values to the ones you've stored. On the frontend, there are four
concepts essential to understand to get your unit tests up and running:
`trusted_values_dict`, `create_test`, your testing module (which
will simply be referred to as `test_file`), and a bash script (which
will simply be referred to as `bash_script`). There is also some
important prerequisite knowledge that may be helpful to grasp before
beginning your testing. There are many functions at play in the backend
as well, all of which are described in the Function Reference. Mastery of these functions may not be
essential to get your tests up-and-running, but some basic understanding
of them with undoubtedly help with debugging.

An important caveat is that the unit testing does not test the
**correctness** of your code or your variables. The unit tests act as a
protective measure to ensure that nothing was broken between versions of
your code; it gets its values by running _your_ code, so if something
starts out incorrect, it will be stored as incorrect in the system.


<a id='dicts'></a>

## Step 1.a: Dictionaries \[Back to [top](#toc)\]
$$\label{dicts}$$

Dictionaries are used throughout the unit testing infrastructure. The user must create simple dictionaries to pass to our testing functions. If you know nothing about dictionaries, we recommend [this](https://www.w3schools.com/python/python_dictionaries.asp) article; it will get you up to speed for simple dictionary creation.

<a id='logging'></a>

## Step 1.b: Logging \[Back to [top](#toc)\]
$$\label{logging}$$

Logging is a Python module that allows the user to specify their desired level of output by modifying a parameter, rather than having to use if-statements and print-statements. We allow the user to change the level of output through a parameter `logging_level`, in which we support the following levels.

`ERROR`: only print when an error occurs.

`INFO`: print general information about test beginning, completion, major function calls, etc., as well as everything above. Recommended. 

`DEBUG`: print maximum amount of information--every comparison, as well as everything above.

A good way to think of these logging levels is that `INFO` is the default, `ERROR` is similar to a non-verbose mode, and `DEBUG` is similar to a verbose mode.

<a id='interaction'></a>

# Step 2: User Interaction \[Back to [top](#toc)\]
$$\label{interaction}$$

Within the module the user is intending to test, a directory named `tests` should be created. This will house the test file for the given module and its associated `trusted_values_dict`. For example, if I intend to test `BSSN`, I will create a new directory `BSSN/tests`. Within the `tests` directory, the user should create a file called `test_(module).py`--or `test_BSSN.py` for the given example. 


<a id='testfile'></a>

## Step 2.a: Test File \[Back to [top](#toc)\]
$$\label{testfile}$$

The test file is how the user inputs their module, functions, and global information to the testing suite. For the purpose of consistency, we've created a skeleton for the test file (found [here](../edit/UnitTesting/test_skeleton.py)) that contains all information the user must specify. The user should change the name of the function to something relevant to their test. However, note that the name of the function must begin with `test_` in order for the bash script to successfully run the test -- this is the default naming scheme for most test suites/software. Inside the function, multiple fields are required to be filled out by the user; these fields are `module`, `module_name`, and `function_and_global_dict`. Below the function, there is some code that begins with `if __name__ == '__main__':`. The user can ignore this code as it does backend work and makes sure to pass the proper information for the test.

`module` is a string representing the module to be tested. 

`module_name` is a string representing the name of the module.

`function_and_global_dict` is a dictionary whose keys are string representations of functions that the user would like to be called on `module` and whose values are lists of string representations of globals that can be acquired by running their respective functions on `module`.

Example:

```
def test_BrillLindquist():

    module = 'BSSN.BrillLindquist'

    module_name = 'bl'

    function_and_global_dict = {'BrillLindquist()': 
                                ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']}
                                
    create_test(module, module_name, function_and_global_dict)
```

In most cases, this simple structure is enough to do exactly what the user wants. Sometimes, however, there is other information that needs to be passed into the test--this is where optional arguments come in.

The tests can take two optional arguments, `logging_level` and `initialization_string_dict`.

`logging_level` follows the same scheme as described [above](#logging). 

`initialization_string_dict` is a dictionary whose keys are functions that **must** also be in `function_and_global_dict` and whose values are strings containing well-formed Python code. The strings are executed as Python code before their respective function is called on the module. The purpose of this argument is to allow the user to do any necessary NRPy+ setup before they call their function.

Example:

```
def test_quantities():

    module = 'BSSN.BSSN_quantities'

    module_name = 'BSSN_quantities'

    function_and_global_dict = {'BSSN_basic_tensors()': ['gammabarDD', 'AbarDD', 'LambdabarU', 'betaU', 'BU']}
                                
    logging_level = 'DEBUG'

    initialization_string = '''
import reference_metric as rfm
rfm.reference_metric()
rfm.ref_metric__hatted_quantities()
'''

    initialization_string_dict = {'BSSN_basic_tensors()': initialization_string}
                                
    create_test(module, module_name, function_and_global_dict, logging_level=logging_level,          
                initialization_string_dict=initialization_string_dict)
```

An important thing to note is that even though `initialization_string` looks odd with its indentation, this is necessary for Python to interpret it correctly. If it was indented, Python would think you were trying to indent that code when it shouldn't be, and an error will occur.

A question you may be wondering is why we need to create a new dictionary for the initialization string, instead of just passing it as its own argument. This is because the testing suite can accept multiple function calls, each with its own associated global list, in one function. It then naturally follows that we need `initailization_string_dict` to allow each function call to have its own code that runs before its function call. In the following example, the function `BSSN_basic_tensors()` has an initialization string, but the function `declare_BSSN_gridfunctions_if_not_declared_already()` doesn't. You can also clearly see they each have their own associated globals.

Example:

```
def test_quantities():

    module = 'BSSN.BSSN_quantities'

    module_name = 'BSSN_quantities'

    function_and_global_dict = {'declare_BSSN_gridfunctions_if_not_declared_already()': 
                                ['hDD', 'aDD', 'lambdaU', 'vetU', 'betU', 'trK', 'cf', 'alpha'], 
                                
                                'BSSN_basic_tensors()': ['gammabarDD', 'AbarDD', 'LambdabarU', 'betaU', 'BU']}
                                
    logging_level = 'DEBUG'

    initialization_string = '''
import reference_metric as rfm
rfm.reference_metric()
rfm.ref_metric__hatted_quantities()
'''

    initialization_string_dict = {'BSSN_basic_tensors()': initialization_string}
                                
    create_test(module, module_name, function_and_global_dict, logging_level=logging_level,          
                initialization_string_dict=initialization_string_dict)
```

Lastly, within a single test file, you can define multiple test functions. It's as simple as defining a new function whose name starts with `test_` in the file and making sure to fill out the necessary fields.

<a id='trustedvaluesdict'></a>

## Step 2.b: trusted_values_dict \[Back to [top](#toc)\]
$$\label{trustedvaluesdict}$$

At this point, it should be understood that our test suite will compare trusted values of your variables to newly calculated values to ensure that no variables were unknowingly modified. The `trusted_values_dict` acts as the means of storing the trusted value for each variable with the purpose of future comparison. A new `trusted_values_dict` is created by default when a test file is run for the first time--it's visible in `tests/trusted_values_dict.py`. Note that if you run your code but can't see the file, refresh your IDE--it's there, sometimes IDE's just get confused when you create a file within Python. The default structure of all `trusted_value_dict` files is as follows.

```
from mpmath import mpf, mp, mpc
from UnitTesting.standard_constants import precision

mp.dps = precision
trusted_values_dict = {}

```

The proper code to copy into this file will be printed to the console when a test is run. The test suite will also automatically write its calculated globals' values for a given function to this file in the proper format--make sure to check that things seem correct, though! Remember that the `trusted_values_dict` stores **trusted**, not necessarily **correct**, values for each global.

<a id='bash'></a>

## Step 2.c: Bash Script \[Back to [top](#toc)\]
$$\label{bash}$$

In order to successfully run all the user's unit tests and properly integrate testing with TravisCI, we use a bash script as the 'hub' of all the tests to be run. This makes it easy for the user to comment out tests they don't want to run, add new tests to be automatically run with one line, etc. 

We offer a skeleton file, [`run_NRPy_UnitTests`](../edit/UnitTesting/run_NRPy_UnitTests.sh), which contains all the proper code to be easily run with minimum user interaction. All the user must do is call the `add_test` function on the test file they'd like to be run underneath the `TODO` comment. There are many examples in the file that show exactly how to create a new test. Then to add more tests, simply go to the next line and add another test. It's as simple as that! 

To run the bash script, open up a terminal, type in the path of the bash script, and then pick the Python interpreter to run the code--for example, `./UnitTesting/run_NRPy_UnitTests.sh python` or `./UnitTesting/run_NRPy_UnitTests.sh python3`.


There's an additional field in the bash script called `rerun_if_fail`. It is a boolean that, if true, will automatically rerun the tests that failed with their `logging_level` set to `DEBUG`. This gives the user a plethora of debugging information that should make it much easier to figure out the issue. We'd recommend enabling it if there are only a couple of modules that failed, as there is a very large amount of information printed for each failing module. However, it is an invaluable resource for figuring out a bug in your code, so keep it in mind when tests are failing. 


<a id='output'></a>

# Step 3: Interpreting Output \[Back to [top](#toc)\]
$$\label{output}$$

Once a user's tests are fully set up, they need to be able to interpret the output of their tests; doing this allows the user to easily figure out what went wrong, why, and how to fix it. The amount of output for a given module is of course dependent on its logging level. For the purposes of this tutorial, we will assume that `logging_level` is set to `INFO`.

While running a test, the output is printed to the console that tells the user what is occurring at what point in time.

Example of successful test run console output:

```
Testing test_u0_smallb_Poynting__Cartesian...

INFO:root: Creating file /home/kevin/virtpypy/nrpyunittesting/u0_smallb_Poynting__Cartesian/tests/u0sbPoyn__compute_u0_smallb_Poynting__Cartesian__test.py...
INFO:root: ...Success: File created.

INFO:root: Currently working on function compute_u0_smallb_Poynting__Cartesian() in module u0sbPoyn...

INFO:root: Importing trusted_values_dict...
INFO:root: ...Success: Imported trusted_values_dict.

INFO:root: Calling evaluate_globals...
INFO:root: ...Success: evaluate_globals ran without errors.

INFO:root: Calling cse_simplify_and_evaluate_sympy_expressions...
INFO:root: ...Success: cse_simplify_and_evaluate_sympy_expressions ran without errors.

INFO:root: Calling calc_error...
INFO:root: ...Success: calc_error ran without errors.

.
----------------------------------------------------------------------
Ran 1 test in 3.550s

OK
INFO:root: Test for function compute_u0_smallb_Poynting__Cartesian() in module u0sbPoyn passed! Deleting test file...
INFO:root: ...Deletion successful. Test complete.

----------------------------------------------------------------------
----------------------------------------------------------------------

All tests passed!

```

Step-by-step, how do we interpret this output? 

The first line tells the user what is being called--this is the function that they define in their current test file being run; this is why a descriptive function name is important.

Next, a file is created in the same directory as the current test file that actually runs the tests--a new file is created to ensure that a clean Python environment is used for each test.

Each function test within the test file is then called successively -- in this example, that is `compute_u0_smallb_Poynting__Cartesian()` in module `u0sbPoyn`. 

The function's associated `trusted_values_dict` is imported.

Each global in the function's global list is then evaluated to a SymPy expression in `evaluate_globals`.

Each global's SymPy expression is then evaluated into a numerical value in `cse_simplify_and_evaluate_sympy_expressions`.

Each global is compared to its trusted value in `calc_error`.

In this example, since no values differed, the test passed--there is nothing the user has to do.

The purpose of giving the user this output is to make it as easy as possible, when something fails, to figure out why and how to fix it. Say the user wasn't given the above output--instead, after `INFO:root: Calling calc_error...` is printed, an error is printed and the program quits--then the user knows the error occurred somewhere in `calc_error`, and it's easy to figure out what. If this output wasn't given, it would be extremely difficult to bugfix.

Now let's consider an example output when the trusted values for a couple of globals differ from the newly calculated values:

```
Testing function test_u0_smallb_Poynting__Cartesian...

INFO:root: Creating file /home/kevin/virtpypy/nrpyunittesting/u0_smallb_Poynting__Cartesian/tests/u0sbPoyn__compute_u0_smallb_Poynting__Cartesian__test.py...
INFO:root: ...Success: File created.

INFO:root: Currently working on function compute_u0_smallb_Poynting__Cartesian() in module u0sbPoyn...

INFO:root: Importing trusted_values_dict...
INFO:root: ...Success: Imported trusted_values_dict.

INFO:root: Calling evaluate_globals...
INFO:root: ...Success: evaluate_globals ran without errors.

INFO:root: Calling cse_simplify_and_evaluate_sympy_expressions...
INFO:root: ...Success: cse_simplify_and_evaluate_sympy_expressions ran without errors.

INFO:root: Calling calc_error...
ERROR:root:

Variable(s) ['g4DD[0][0]', 'g4DD[0][1]'] in module u0sbPoyn failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for 
u0sbPoyn__compute_u0_smallb_Poynting__Cartesian__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict. 
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: 2019-08-14
# Reason for changing values: TODO
trusted_values_dict['u0sbPoyn__compute_u0_smallb_Poynting__Cartesian__globals'] = {'g4DD[0][0]': mpf('1.42770464273047624140299713523'), 'g4DD[0][1]': mpf('0.813388473397507463814385913308'), 'g4DD[0][2]': mpf('0.652706348793296836714132090803'), 'g4DD[0][3]': mpf('1.22429414375154980405074869244'), 'g4DD[1][0]': mpf('0.813388473397507463814385913308'), 'g4DD[1][1]': mpf('0.657497767033916602485987823457'), 'g4DD[1][2]': mpf('0.057738705167452830657737194997'), 'g4DD[1][3]': mpf('0.391026617743468030141684721457'), 'g4DD[2][0]': mpf('0.652706348793296836714132090803'), 'g4DD[2][1]': mpf('0.057738705167452830657737194997'), 'g4DD[2][2]': mpf('0.142350778742078798444481435581'), 'g4DD[2][3]': mpf('0.723120760610660329170684690325'), 'g4DD[3][0]': mpf('1.22429414375154980405074869244'), 'g4DD[3][1]': mpf('0.391026617743468030141684721457'), 'g4DD[3][2]': mpf('0.723120760610660329170684690325'), 'g4DD[3][3]': mpf('0.919283767179900235255729512573'), 'g4UU[0][0]': mpf('-3.03008926847944211197781568781'), 'g4UU[0][1]': mpf('2.25487680174746330097911429618'), 'g4UU[0][2]': mpf('0.883964088219310673292773829627'), 'g4UU[0][3]': mpf('2.38097417962378184338842987037'), 'g4UU[1][0]': mpf('2.25487680174746330097911429618'), 'g4UU[1][1]': mpf('-0.109478866681308257858746913533'), 'g4UU[1][2]': mpf('-1.57673807708112257475456728614'), 'g4UU[1][3]': mpf('-1.71617440778374167644096697698'), 'g4UU[2][0]': mpf('0.883964088219310673292773829627'), 'g4UU[2][1]': mpf('-1.57673807708112257475456728614'), 'g4UU[2][2]': mpf('-2.06437348264308527564608946534'), 'g4UU[2][3]': mpf('1.11728919891515454886216209391'), 'g4UU[3][0]': mpf('2.38097417962378184338842987037'), 'g4UU[3][1]': mpf('-1.71617440778374167644096697698'), 'g4UU[3][2]': mpf('1.11728919891515454886216209391'), 'g4UU[3][3]': mpf('-2.23203972375107587678882970577'), 'PoynSU[0]': mpf('0.103073801363157111172177901697'), 'PoynSU[1]': mpf('0.11100316917740755485837448786'), 'PoynSU[2]': mpf('-0.00451075406485067218999293829888'), 'smallb2etk': mpf('0.164454779456120937541853683919'), 'smallb4D[0]': mpf('0.567950228622914592095169687713'), 'smallb4D[1]': mpf('0.286535540626704686153523625219'), 'smallb4D[2]': mpf('0.10714698030450909234705495631'), 'smallb4D[3]': mpf('0.455828728852934996540499932291'), 'smallb4U[0]': mpf('0.105192967134481035810628308481'), 'smallb4U[1]': mpf('0.298063886336868595407751154048'), 'smallb4U[2]': mpf('0.338357239072152954142217077157'), 'smallb4U[3]': mpf('-0.0371837983175520496394928051176'), 'u0': mpf('0.751914772923022001194226504595'), 'uBcontraction': mpf('0.214221928967111307128784385185'), 'uD[0]': mpf('0.216251888123560253383388291707'), 'uD[1]': mpf('0.167535113620266400428280039145'), 'uD[2]': mpf('0.232536332826570514618343904792'), 'uU[0]': mpf('-0.364066660324468905733189036042'), 'uU[1]': mpf('-0.0378849772494775056256865716175'), 'uU[2]': mpf('-0.476480636313712572229280970632')}

#####

.
----------------------------------------------------------------------
Ran 1 test in 3.481s

OK
ERROR:root: Test for function compute_u0_smallb_Poynting__Cartesian() in module u0sbPoyn failed! Please examine test file.

----------------------------------------------------------------------
----------------------------------------------------------------------

Tests failed!

Failures:

u0_smallb_Poynting__Cartesian/tests/test_u0_smallb_Poynting__Cartesian.py: ['test_u0_smallb_Poynting__Cartesian'] 

----------------------------------------------------------------------


```

This seems like a lot to take in, but it's not too difficult to understand once fairly well acquainted with the output. The beginning is identical to the output from the successful test run, up until `calc_error` is called. This gives the user some insight that there was an error during `calc_error`, which gives an indication that at least one global had a different calculated and trusted value.

The next line confirms this suspicion: globals `g4DD[0][0]` and `g4DD[0][1]` had differing trusted and calculated values. This gives the user the **exact** information they're looking for. `g4DD` is a rank-2 tensor as per NRPy+ naming convention, and the user can very easily see the indices that failed. If the user expected this--say they found a bug in their code that generated this global--they can then copy the new `trusted_values_dict` entry into their `trusted_values_dict` and comment out/delete the old entry. Then by re-running the test, there should no longer be an error--the trusted and calculated values should be the same.

The next output tells the user that the test failed, and to examine the test file. This is necessary if there's an unexpected failure--it will help the user figure out why it occurred.

Finally, an additional output is given that tells the user all the test files and their respective functions that failed. This may seem repetitive, but there's a good reason for it. Say 20 modules were being tested, 10 of which had failures. Then to figure out what failed, the user would have to scroll through all the output and keep a mental note of what failed. By printing everything that failed at the very end, the user gets instant insight into what failed; this may help the user figure out why.


<a id='checklist'></a>

# Step 4: Checklist for adding a new test \[Back to [top](#toc)\]
$$\label{checklist}$$

<a id='4a'></a>

## Step 4.a: Directory Creation \[Back to [top](#toc)\]
$$\label{4a}$$

Create a `tests` directory in the directory of the module being tested.

<a id='4b'></a>

## Step 4.b: Test File Creation  \[Back to [top](#toc)\]
$$\label{4b}$$

Create a `test_file.py` in the `tests` directory based off [UnitTesting/test_skeleton.py](../edit/UnitTesting/test_skeleton.py).

<a id='4c'></a>

## Step 4c: Input Parameters \[Back to [top](#toc)\]
$$\label{4c}$$

Change the name of `test_your_module()` to whatever you're testing, making sure the name starts with `test_`. Fill in the following paremeters in `test_file.py`: `module`, `module_name`, `function_and_global_dict`, and if need be, `logging_level` and `initialization_string`.

<a id='4d'></a>

## Step 4d: Fill-in bash script and run \[Back to [top](#toc)\]
$$\label{4d}$$

Use the `add_test` function in the [run_NRPy_UnitTests.sh](../edit/UnitTesting/run_NRPy_UnitTests.sh) to create your test below the `TODO` line. Run the bash script to run your test!

<a id='latex_pdf_output'></a>

# Step 5: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename [Tutorial-UnitTesting.pdf](Tutorial-UnitTesting.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-UnitTesting")
```

    Created Tutorial-UnitTesting.tex, and compiled LaTeX file to PDF file
        Tutorial-UnitTesting.pdf

