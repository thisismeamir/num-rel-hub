<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# `cmdline_helper.py`: Multi-platform command-line helper functions
## Authors: Brandon Clark & Zach Etienne

## The notebook presents [cmdline_helper.py](../edit/cmdline_helper.py), a multi-platform Python module for interacting with the command-line within the NRPy+ tutorial framework. It provides functions for checking if an executable exists, compiling and executing C code, executing input strings, and file manipulations such as deletion or directory creation. Its functionality ensures compatibility across Linux, Windows, and Mac OS command-line interfaces.

**Notebook Status:** <font color=orange><b> Self-Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). **Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**

### NRPy+ Python Module associated with this notebook: [cmdline_helper.py](../edit/cmdline_helper.py)

## Introduction:
Throughout the NRPy+ tutorial, there are a handful of modules that require interaction with the command line, to compile C code, manipulate files, execute code, etc. This module serves as a reference for Python functions that exist in [cmdline_helper](../edit/cmdline_helper.py), which is designed to be compatible with Linux, Windows, and Mac OS command line interfaces.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
1. [Step 2](#functions): The Functions
    1. [Step 2.a](#checkexec): **`check_executable_exists()`**
    1. [Step 2.b](#compile): **`C_compile()`** and **`new_C_compile()`**
    1. [Step 2.c](#execute): **`Execute()`**
    1. [Step 2.d](#output): **`Execute_input_string()`**
    1. [Step 2.e](#delete):  **`delete_existing_files()`** & **`mkdir()`**
1. [Step 3](#code_validation): Code Validation against `cmdline_helper.py`NRPy+ module
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize core NRPy+ modules \[Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing the necessary Python modules.


```python
%%writefile cmdline_helper-validation.py
# As documented in the NRPy+ tutorial notebook
# Tutorial-cmdline_helper.ipynb, this Python script
# provides a multi-platform means to run executables,
# remove files, and compile code.

# Basic functions:
# check_executable_exists(): Check to see whether an executable exists.
#                            Error out or return False if it does not exist;
#                            return True if executable exists in PATH.
# C_compile(): Compile C code using gcc.
# Execute(): Execute generated executable file, using taskset
#            if available. Calls Execute_input_string() to
#            redirect output from stdout & stderr to desired
#            destinations.
# Execute_input_string(): Executes an input string and redirects
#            output from stdout & stderr to desired destinations.
# delete_existing_files(file_or_wildcard):
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS

# Authors: Brandon Clark
#          Zach Etienne
#          zachetie **at** gmail **dot* com
#          Kevin Lituchy

import io, os, shlex, subprocess, sys, time, multiprocessing, getpass, platform, glob
```

    Overwriting cmdline_helper-validation.py


<a id='functions'></a>

# Step 2: The Functions \[Back to [top](#toc)\]
$$\label{functions}$$

<a id='checkexec'></a>

## Step 2.a: `check_executable_exists()` \[Back to [top](#toc)\]
$$\label{checkexec}$$

`check_executable_exists()` takes the required string `exec_name` (i.e., the name of the executable) as its first input. Its second input is the optional boolean `error_if_not_found`, which defaults to `True` (so that it exits with an error if the executable is not found).

`check_executable_exists()` returns `True` if the executable exists, `False` if the executable does not exist, and `error_if_not_found` is set to `False`.


```python
%%writefile -a cmdline_helper-validation.py

# check_executable_exists(): Check to see whether an executable exists.
#                            Error out or return False if it does not exist;
#                            return True if executable exists in PATH.
def check_executable_exists(exec_name, error_if_not_found=True):
    cmd = "where" if os.name == "nt" else "which"
    try:
        subprocess.check_output([cmd, exec_name])
    except subprocess.CalledProcessError:
        if error_if_not_found:
            print("Sorry, cannot execute the command: " + exec_name)
            sys.exit(1)
        else:
            return False
    return True
```

    Appending to cmdline_helper-validation.py


<a id='compile'></a>

## Step 2.b: **`C_compile()`** and **`new_C_compile()`** \[Back to [top](#toc)\]
$$\label{compile}$$

### `C_compile()`

The `C_compile()` function takes the following inputs as **strings**
* Path name to the generated C_file, `"main_C_output_path"`, and
* Name of the executable playground file, `"main_C_output_file"`.

The `C_compile()` function first checks for a ***gcc compiler***, which is a must when compiling C code within the NRPy+ tutorial. The function then removes any existing executable file. After that, the function constructs a script `compile_string` to run the compilation based on the function inputs and the operating system (OS) in use.

Finally, it runs the actual compilation, by passing the compilation script `compile_string` on to the `Execute_input_string()` function, see [Step 2.d](#output).


```python
%%writefile -a cmdline_helper-validation.py

# C_compile(): Write a function to compile the Main C code into an executable file
def C_compile(main_C_output_path, main_C_output_file, compile_mode="optimized", custom_compile_string="", additional_libraries=""):
    print("Compiling executable...")
    # Step 1: Check for gcc compiler
    check_executable_exists("gcc")

    if additional_libraries != "":
        additional_libraries = " " + additional_libraries

    # Step 2: Delete existing version of executable
    if os.name == "nt":
        main_C_output_file += ".exe"
    delete_existing_files(main_C_output_file)

    # Step 3: Compile the executable
    if compile_mode=="safe":
        compile_string = "gcc -std=gnu99 -O2 -g -fopenmp "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            # Step 3.A: Maybe gcc is actually clang in disguise (as in MacOS)?!
            #           https://stackoverflow.com/questions/33357029/using-openmp-with-clang
            print("Most safe failed. Probably on MacOS. Replacing -fopenmp with -fopenmp=libomp:")
            compile_string = "gcc -std=gnu99 -O2 -fopenmp=libomp "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="icc":
        check_executable_exists("icc")
        compile_string = "icc -std=gnu99 -O2 -xHost -qopenmp -unroll "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="custom":
        Execute_input_string(custom_compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="optimized":
        compile_string = "gcc -std=gnu99 -Ofast -fopenmp -march=native -funroll-loops "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            # Step 3.A: Revert to more compatible gcc compile option
            print("Most optimized compilation failed. Removing -march=native:")
            compile_string = "gcc -std=gnu99 -Ofast -fopenmp -funroll-loops "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        if not os.path.isfile(main_C_output_file):
            # Step 3.B: Revert to maximally compatible gcc compile option
            print("Next-to-most optimized compilation failed. Moving to maximally-compatible gcc compile option:")
            compile_string = "gcc -std=gnu99 -O2 "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        # Step 3.C: If there are still missing components within the compiler, say compilation failed
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="emscripten":
        compile_string = "emcc -std=gnu99 -s -O3 -march=native -funroll-loops -s ALLOW_MEMORY_GROWTH=1 "\
            +str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm "+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        # If there are still missing components within the compiler, say compilation failed
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed.")
            sys.exit(1)
    else:
        print("Sorry, compile_mode = \""+compile_mode+"\" unsupported.")
        sys.exit(1)

    print("Finished compilation.")
```

    Appending to cmdline_helper-validation.py


### `new_C_compile()`

The `new_C_compile()` function first constructs a `Makefile` from all functions registered to `outputC`'s `outC_function_dict`, and then attempts to compile the code using a parallel `make`. If that fails (e.g., due to the `make` command not being found or optimizations not being supported), it instead constructs a script that builds the code in serial with most compiler optimizations disabled.

`new_C_compile()` takes the following inputs:
* `Ccodesrootdir`: Path to the generated C code, `Makefile`, and executable.
* `exec_name`: File name of executable.
* `uses_free_parameters_h=False` (optional): Will cause the compilation of `main.c` to depend on `free_parameters.h`, such that if the latter is updated `main.c` will be recompiled due to `make` being run.
* `compiler_opt_option="fast"` (optional): Optimization option for compilation. Choose from `fast`, `fastdebug`, or `debug`. The latter two will enable the `-g` flag.
* `addl_CFLAGS=None` (optional): A list of additional compiler flags. E.g., `[-funroll-loops,-ffast-math]`.
* `addl_libraries=None` (optional): A list of additional libraries against which to link.
* `mkdir_Ccodesrootdir=True` (optional): Attempt to make the `Ccodesrootdir` directory if it doesn't exist. If it does exist, this has no effect.
* `CC="gcc"` (optional): Choose the compiler. `"gcc"` should be used for the C compiler and `"g++"` for C++. FIXME: Other compilers may throw warnings due to default compilation flags being incompatible.
* `attempt=1` (optional, **do not touch**): An internal flag used to recursively call this function again in case the `Makefile` build fails (runs a shell script in serial with debug options disabled, as a backup).


```python
%%writefile -a cmdline_helper-validation.py


from outputC import construct_Makefile_from_outC_function_dict
def new_C_compile(Ccodesrootdir, exec_name, uses_free_parameters_h=False,
                  compiler_opt_option="fast", addl_CFLAGS=None,
                  addl_libraries=None, mkdir_Ccodesrootdir=True, CC="gcc", attempt=1):
    check_executable_exists("gcc")
    use_make = check_executable_exists("make", error_if_not_found=False)

    construct_Makefile_from_outC_function_dict(Ccodesrootdir, exec_name, uses_free_parameters_h,
                                               compiler_opt_option, addl_CFLAGS,
                                               addl_libraries, mkdir_Ccodesrootdir, use_make, CC=CC)
    orig_working_directory = os.getcwd()
    os.chdir(Ccodesrootdir)
    if use_make:
        Execute_input_string("make -j" + str(int(multiprocessing.cpu_count()) + 2), os.devnull)
    else:
        Execute_input_string(os.path.join("./", "backup_script_nomake.sh"))
    os.chdir(orig_working_directory)

    if not os.path.isfile(os.path.join(Ccodesrootdir, exec_name)) and attempt == 1:
        print("Optimized compilation FAILED. Removing optimizations (including OpenMP) and retrying with debug enabled...")
        # First clean up object files.
        filelist = glob.glob(os.path.join(Ccodesrootdir, "*.o"))
        for file in filelist:
            os.remove(file)
        # Then retry compilation (recursion)
        new_C_compile(Ccodesrootdir, exec_name, uses_free_parameters_h,
                      compiler_opt_option="debug", addl_CFLAGS=addl_CFLAGS,
                      addl_libraries=addl_libraries, mkdir_Ccodesrootdir=mkdir_Ccodesrootdir, CC=CC, attempt=2)
    if not os.path.isfile(os.path.join(Ccodesrootdir, exec_name)) and attempt == 2:
        print("Sorry, compilation failed")
        sys.exit(1)
    print("Finished compilation.")

```

    Appending to cmdline_helper-validation.py


<a id='execute'></a>

## Step 2.c: `Execute()` \[Back to [top](#toc)\]
$$\label{execute}$$

The `Execute()` function takes the following inputs as **strings**
* Name of the executable file, `main_C_output_file`,
* **(Optional):** Any necessary arguments associated with the executable file output, `executable_output_arguments`, and 
* **(Optional):** Name of a file to store output during execution, `executable_output_file_name`.

The `Execute()` function first removes any existing output files. It then begins to construct the script `execute_string` in order to execute the executable file that has been generated by the `C_compile()` function. `execute_string` is built based on the function inputs and the operating system (OS) in use.

Finally, it runs the actual execution, by passing the execution script `execute_string` on to the `Execute_input_string()` function, see [Step 2.d](#output).


```python
%%writefile -a cmdline_helper-validation.py


# Execute(): Execute generated executable file, using taskset
#            if available. Calls Execute_input_string() to
#            redirect output from stdout & stderr to desired
#            destinations.
def Execute(executable, executable_output_arguments="", file_to_redirect_stdout=os.devnull, verbose=True):
    # Step 1: Delete old version of executable file
    if file_to_redirect_stdout != os.devnull:
        delete_existing_files(file_to_redirect_stdout)

    # Step 2: Build the script for executing the desired executable
    execute_string = ""
    # When in Windows...
    # https://stackoverflow.com/questions/1325581/how-do-i-check-if-im-running-on-windows-in-python
    if os.name == "nt":
        # ... do as the Windows do
        # https://stackoverflow.com/questions/49018413/filenotfounderror-subprocess-popendir-windows-7
        execute_prefix = "cmd /c " # Run with cmd /c executable [options] on Windows
    else:
        execute_prefix = "./"      # Run with ./executable [options] on Linux & Mac
    taskset_exists = check_executable_exists("taskset", error_if_not_found=False)
    if taskset_exists:
        execute_string += "taskset -c 0"
        on_4900hs = False
        if platform.system() == "Linux" and \
                "AMD Ryzen 9 4900HS" in str(subprocess.check_output("cat /proc/cpuinfo", shell=True)):
            on_4900hs = True
        if not on_4900hs and getpass.getuser() != "jovyan": # on mybinder, username is jovyan, and taskset -c 0 is the fastest option.
            # If not on mybinder and taskset exists:
            has_HT_cores = False  # Does CPU have hyperthreading cores?
            if platform.processor() != '': # If processor string returns null, then assume CPU does not support hyperthreading.
                                           # This will yield correct behavior on ARM (e.g., cell phone) CPUs.
                has_HT_cores=True
            if has_HT_cores == True:
                # NOTE: You will observe a speed-up by using only *PHYSICAL* (as opposed to logical/hyperthreading) cores:
                N_cores_to_use = int(multiprocessing.cpu_count()/2) # To account for hyperthreading cores
            else:
                N_cores_to_use = int(multiprocessing.cpu_count()) # Use all cores if none are hyperthreading cores.
                                                                  # This will happen on ARM (e.g., cellphone) CPUs
            for i in range(N_cores_to_use-1):
                execute_string += ","+str(i+1)
        if on_4900hs:
            execute_string = "taskset -c 1,3,5,7,9,11,13,15"
        execute_string += " "
    execute_string += execute_prefix+executable+" "+executable_output_arguments

    # Step 3: Execute the desired executable
    Execute_input_string(execute_string, file_to_redirect_stdout, verbose)
```

    Appending to cmdline_helper-validation.py


<a id='output'></a>

## Step 2.d:  `Execute_input_string()` \[Back to [top](#toc)\]
$$\label{output}$$

The `Execute_input_string()` function takes the following inputs as strings
* The script to be executed, `input_string`, and
* An output file name for any needed redirects, `executable_output_file_name`. 

The `Execute_input_string()` executes a script, outputting `stderr` to the screen and redirecting any additional outputs from the executable to the specified `executable_output_file_name`. 



```python
%%writefile -a cmdline_helper-validation.py

# Execute_input_string(): Executes an input string and redirects
#            output from stdout & stderr to desired destinations.
def Execute_input_string(input_string, file_to_redirect_stdout=os.devnull, verbose=True):

    if verbose:
        print("(EXEC): Executing `"+input_string+"`...")
    start = time.time()
    # https://docs.python.org/3/library/subprocess.html
    if os.name != 'nt':
        args = shlex.split(input_string)
    else:
        args = input_string

    # https://stackoverflow.com/questions/18421757/live-output-from-subprocess-command
    filename = "tmp.txt"
    with io.open(filename, 'w') as writer, io.open(filename, 'rb', buffering=-1) as reader, io.open(file_to_redirect_stdout, 'wb') as rdirect:
        process = subprocess.Popen(args, stdout=rdirect, stderr=writer)
        while process.poll() is None:
            # https://stackoverflow.com/questions/21689365/python-3-typeerror-must-be-str-not-bytes-with-sys-stdout-write/21689447
            sys.stdout.write(reader.read().decode('utf-8'))
            time.sleep(0.2)
        # Read the remaining
        sys.stdout.write(reader.read().decode('utf-8'))
    delete_existing_files(filename)
    end = time.time()
    if verbose:
        print("(BENCH): Finished executing in "+'{:#.2f}'.format(round(end-start, 2))+" seconds.")
```

    Appending to cmdline_helper-validation.py


<a id='delete'></a>

## Step 2.e:  `delete_existing_files()` & `mkdir()` \[Back to [top](#toc)\]
$$\label{delete}$$

The `delete_existing_files()` function takes a string, `file_or_wildcard`, as input.

`delete_existing_files()` deletes any existing files that match the pattern given by `file_or_wildcard`. Deleting files is important when running the same code multiple times, ensuring that you're not reusing old data from a previous run, or seeing the same plot from a previous output. 

The `mkdir()` function makes a directory if it does not yet exist. It passes the input string "newpath" through `os.path.join()` to ensure that forward slashes are replaced by backslashes in Windows environments.


```python
%%writefile -a cmdline_helper-validation.py

# delete_existing_files(file_or_wildcard):
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS
def delete_existing_files(file_or_wildcard):
    delete_string = ""
    if os.name == "nt":
        delete_string += "del " + file_or_wildcard
    else:
        delete_string += "rm -f " + file_or_wildcard
    os.system(delete_string)

# https://stackoverflow.com/questions/1274405/how-to-create-new-folder
def mkdir(newpath):
    if not os.path.exists(os.path.join(newpath)):
        os.makedirs(os.path.join(newpath))
```

    Appending to cmdline_helper-validation.py



```python
%%writefile -a cmdline_helper-validation.py

# TO BE RUN ONLY FROM nrpytutorial or nrpytutorial/subdir/
def output_Jupyter_notebook_to_LaTeXed_PDF(notebookname, verbose=True):
    in_nrpytutorial_rootdir = os.getcwd().split("/")[-1] == "nrpytutorial"
    if sys.version_info[0] == 3:
        location_of_template_file = "."
        if not in_nrpytutorial_rootdir:
            location_of_template_file = ".."
        Execute_input_string(r"jupyter nbconvert --to latex --template="
                             +os.path.join(location_of_template_file, "nbconvert_latex_settings")
                             +r" --log-level='WARN' "+notebookname+".ipynb",verbose=False)
    else:
        Execute_input_string(r"jupyter nbconvert --to latex --log-level='WARN' "+notebookname+".ipynb",verbose=False)
    for _i in range(3):  # _i is an unused variable.
        Execute_input_string(r"pdflatex -interaction=batchmode "+notebookname+".tex",verbose=False)
    delete_existing_files(notebookname+".out "+notebookname+".aux "+notebookname+".log")
    if verbose:
        import textwrap
        wrapper = textwrap.TextWrapper(initial_indent="",subsequent_indent="    ",width=75)
        print(wrapper.fill("Created "+notebookname+".tex, and compiled LaTeX file to PDF file "+notebookname+".pdf"))
```

    Appending to cmdline_helper-validation.py


<a id='code_validation'></a>

# Step 3: Code Validation against `cmdline_helper.py` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

To validate the code in this tutorial we check for agreement between the files

1. `cmdline_helper-validation.py` (written in this tutorial) and
1. the NRPy+ [cmdline_helper.py](../edit/cmdline_helper.py) module



```python
import difflib
import sys
def trim_lines(lines):
    new_lines=[]
    for line in lines:
        x = line.rstrip()
        if x != '':
             new_lines += [x]
    return new_lines

print("Printing difference between original cmdline_helper.py and this code, cmdline_helper-validation.py.")
# Open the files to compare
with open("cmdline_helper.py") as file1, open("cmdline_helper-validation.py") as file2:
    # Read the lines of each file
    file1_lines = trim_lines(file1.readlines())
    file2_lines = trim_lines(file2.readlines())
    num_diffs = 0
    for line in difflib.unified_diff(file1_lines, file2_lines, fromfile="cmdline_helper.py", tofile="cmdline_helper-validation.py"):
        sys.stdout.writelines(line)
        num_diffs = num_diffs + 1
    if num_diffs == 0:
        print("No difference. TEST PASSED!")
    else:
        print("ERROR: Disagreement found with .py file. See differences above.")
        sys.exit(1)
```

    Printing difference between original cmdline_helper.py and this code, cmdline_helper-validation.py.



    ---------------------------------------------------------------------------

    FileNotFoundError                         Traceback (most recent call last)

    Cell In[1], line 13
         11 print("Printing difference between original cmdline_helper.py and this code, cmdline_helper-validation.py.")
         12 # Open the files to compare
    ---> 13 with open("cmdline_helper.py") as file1, open("cmdline_helper-validation.py") as file2:
         14     # Read the lines of each file
         15     file1_lines = trim_lines(file1.readlines())
         16     file2_lines = trim_lines(file2.readlines())


    File /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/IPython/core/interactiveshell.py:284, in _modified_open(file, *args, **kwargs)
        277 if file in {0, 1, 2}:
        278     raise ValueError(
        279         f"IPython won't let you open fd={file} by default "
        280         "as it is likely to crash IPython. If you know what you are doing, "
        281         "you can use builtins' open."
        282     )
    --> 284 return io_open(file, *args, **kwargs)


    FileNotFoundError: [Errno 2] No such file or directory: 'cmdline_helper-validation.py'


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-cmdline_helper.pdf](Tutorial-cmdline_helper.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-cmdline_helper")
```

    Created Tutorial-cmdline_helper.tex, and compiled LaTeX file to PDF file
        Tutorial-cmdline_helper.pdf

