<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Calculating Spin-Weighted Spherical Harmonics
## Authors: Zach Etienne & Brandon Clark

## This notebook presents a Python code, designed to calculate spin-weighted spherical harmonics utilizing Sympy and the Goldberg function [Goldberg et al. (1967)](https://aip.scitation.org/doi/10.1063/1.1705135). The implementation is verified against a recognized Mathematica notebook and related NRPy+ module, and is then output as a C code.

**Notebook Status:** <font color='green'><b> Validated </b></font>

**Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). In addition, its results have been validated against a [trusted Mathematica notebook](https://demonstrations.wolfram.com/versions/source.jsp?id=SpinWeightedSphericalHarmonics&version=0012).

### NRPy+ Source Code for this module: [SpinWeight_minus2_SphHarmonics/SpinWeight_minus2_SphHarmonics.py](../edit/SpinWeight_minus2_SphHarmonics/SpinWeight_minus2_SphHarmonics.py)

## Introduction:
This tutorial notebook defines a Python function for computing spin-weighted spherical harmonics using Sympy. Spin-weight $s=-2$ spherical harmonics are the natural basis for decomposing gravitational wave data.

The tutorial contains code necessary to validate the resulting expressions assuming $s=-2$ against a trusted Mathematica notebook (validated for all $(\ell,m)$ up to $\ell=8$. Finally it outputs a C code capable of computing $_{-2}Y_{\ell m} (\theta, \phi)$ for all $(\ell,m)$ for $\ell=0$ up to `maximum_l`.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows:

1. [Step 1](#initializenrpy): Initialize needed Python/NRPy+ modules
1. [Step 2](#gbf): Defining the Goldberg function
1. [Step 3](#math_code_validation): Code Validation against Mathematica script
1. [Step 4](#ccode): Generate C-code function for computing s=-2 spin-weighted spherical harmonics, using NRPy+
1. [Step 5](#code_validation): Code Validation against SpinWeight_minus2_SphHarmonics/SpinWeight_minus2_SphHarmonics NRPy+ module
1. [Step 6](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='initializenrpy'></a>

# Step 1: Initialize needed Python/NRPy+ modules [Back to [top](#toc)\]
$$\label{initializenrpy}$$

Let's start by importing all the needed modules from NRPy+:


```python
# Step 1: Initialize needed Python/NRPy+ modules
from outputC import outputC   # NRPy+: Core C code output module
import sympy as sp            # SymPy: The Python computer algebra package upon which NRPy+ depends
import os, sys                # Standard Python modules for multiplatform OS-level functions

# Step 1.a: Set maximum l to which we will validate the spin-weighted spherical harmonics with s=-2:
maximum_l = 4 # Note that we have validated against Mathematica up to and including l=8 -- perfect agreement.
```

<a id='gbf'></a>

# Step 2: Defining the Goldberg function [Back to [top](#toc)\]
$$\label{gbf}$$

One way to calculate the spin-weighted spherical harmonics is using the following formula
from [Goldberg et al. (1967)](https://aip.scitation.org/doi/10.1063/1.1705135):

$$ _sY_{\ell m} (\theta, \phi) = \left(-1\right)^m \sqrt{ \frac{(\ell+m)! (\ell-m)! (2\ell+1)} {4\pi (\ell+s)! (\ell-s)!} } \sin^{2\ell} \left( \frac{\theta}{2} \right) \times\sum_{r=0}^{\ell-s} {\ell-s \choose r} {\ell+s \choose r+s-m} \left(-1\right)^{\ell-r-s} e^{i m \phi} \cot^{2r+s-m} \left( \frac{\theta} {2} \right)$$


```python
# Step 2: Defining the Goldberg function

# Step 2.a: Declare SymPy symbols:
th, ph = sp.symbols('th ph',real=True)

# Step 2.b: Define the Goldberg formula for spin-weighted spherical harmonics
#           (https://aip.scitation.org/doi/10.1063/1.1705135);
#           referenced & described in Wikipedia Spin-weighted spherical harmonics article:
#           https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244
def Y(s, l, m, th, ph, GenerateMathematicaCode=False):
    Sum = 0
    for r in range(l-s + 1):
        if GenerateMathematicaCode == True:
            # Mathematica needs expression to be in terms of cotangent, so that code validation below
            #    yields identity with existing Mathematica notebook on spin-weighted spherical harmonics.
            Sum +=  sp.binomial(l-s, r)*sp.binomial(l+s, r+s-m)*(-1)**(l-r-s)*sp.exp(sp.I*m*ph)*sp.cot(th/2)**(2*r+s-m)
        else:
            # SymPy C code generation cannot handle the cotangent function, so define cot(th/2) as 1/tan(th/2):
            Sum +=  sp.binomial(l-s, r)*sp.binomial(l+s, r+s-m)*(-1)**(l-r-s)*sp.exp(sp.I*m*ph)/sp.tan(th/2)**(2*r+s-m)

    return (-1)**m*sp.simplify(sp.sqrt(sp.factorial(l+m)*sp.factorial(l-m)*(2*l+1)/(4*sp.pi*sp.factorial(l+s)*sp.factorial(l-s)))*sp.sin(th/2)**(2*l)*Sum)
```

<a id='math_code_validation'></a>

# Step 3: Code Validation against Mathematica script [Back to [top](#toc)\]
$$\label{math_code_validation}$$

To validate the code we wish to compare it with an existent [Mathematica notebook](https://demonstrations.wolfram.com/versions/source.jsp?id=SpinWeightedSphericalHarmonics&version=0012). We will validate the code using a spin-value of $s=-2$ and $\ell = 8,7,6,5,4,3,2,1,0$ while leaving $m$, $\theta$, and $\phi$ unknown. 


```python
# Step 3: Code Validation against Mathematica notebook:
#         https://demonstrations.wolfram.com/versions/source.jsp?id=SpinWeightedSphericalHarmonics&version=0012

# # For the l=0 case m=0, otherwise there is a divide-by-zero in the Y() function above.
# print("FullSimplify[Y[-2, 0, 0, th, ph]-"+str(sp.mathematica_code(sp.simplify(Y(-2, 0, 0, th, ph,GenerateMathematicaCode=True))))+"] \n") # Agrees with Mathematica notebook for l = 0

# # Check the other cases
# for l in range(1,maximum_l+1): # Agrees with Mathematica notebook for  l = 1, 2, 4, 5, 6, 7, 8;
#     print("FullSimplify[Y[-2, "+str(l)+", m, th, ph]-("+
#           str(sp.mathematica_code(sp.simplify(Y(-2, l, m, th, ph, GenerateMathematicaCode=True)))).replace("binomial","Binomial").replace("factorial","Factorial")+")] \n")
```

<a id='ccode'></a>

# Step 4: Generate C-code function for computing s=-2 spin-weighted spherical harmonics, using NRPy+ [Back to [top](#toc)\]
$$\label{ccode}$$


```python
# Step 4: Generating C Code function for computing
#         s=-2 spin-weighted spherical harmonics,
#         using NRPy+'s outputC() function.

outCparams = "preindent=3,outCfileaccess=a,outCverbose=False,includebraces=True"

with open(os.path.join("SpinWeight_minus2_SphHarmonics","SpinWeight_minus2_SphHarmonics.h"), "w") as file:
    file.write("""
void SpinWeight_minus2_SphHarmonics(const int l, const int m, const REAL th, const REAL ph,
                                   REAL *reYlmswm2_l_m, REAL *imYlmswm2_l_m) {
if(l<0 || l>"""+str(maximum_l)+""" || m<-l || m>+l) {
    printf("ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,"""+str(maximum_l)+"""] and only m=[-l,+l] is defined.\\n");
    printf("       You chose l=%d and m=%d, which is out of these bounds.\\n",l,m);
    exit(1);
}\n""")

    file.write("switch(l) {\n")
    for l in range(maximum_l+1): # Output values up to and including l=8.
        file.write("    case "+str(l)+":\n")
        file.write("        switch(m) {\n")
        for m in range(-l,l+1):
            file.write("            case "+str(m)+":\n")
            Y_m2_lm = Y(-2, l, m, th, ph)
            Cstring = outputC([sp.re(Y_m2_lm),sp.im(Y_m2_lm)],["*reYlmswm2_l_m","*imYlmswm2_l_m"],
                              "returnstring",outCparams)
            file.write(Cstring)
            file.write("                  return;\n")
        file.write("        }  /* End switch(m) */\n")
    file.write("    } /* End switch(l) */\n")
    file.write("} /* End function SpinWeight_minus2_SphHarmonics() */\n")
```

<a id='code_validation'></a>

# [Step 5](#code_validation): Code Validation against `SpinWeight_minus2_SphHarmonics.SpinWeight_minus2_SphHarmonics` NRPy+ module \[Back to [top](#toc)\]
$$\label{code_validation}$$

As additional validation, we verify agreement in the SymPy expressions for the spin-weight -2 spherical harmonics expressions between
1. this tutorial and 
2. the NRPy+ [`SpinWeight_minus2_SphHarmonics.SpinWeight_minus2_SphHarmonics`](../edit/SpinWeight_minus2_SphHarmonics/SpinWeight_minus2_SphHarmonics.py) module.


```python
import SpinWeight_minus2_SphHarmonics.SpinWeight_minus2_SphHarmonics as swm2
swm2.SpinWeight_minus2_SphHarmonics(maximum_l=4,filename=os.path.join("SpinWeight_minus2_SphHarmonics","SpinWeight_minus2_SphHarmonics-NRPymodule.h"))

print("\n\n### BEGIN VALIDATION TESTS ###")
import filecmp
fileprefix = os.path.join("SpinWeight_minus2_SphHarmonics","SpinWeight_minus2_SphHarmonics")
if filecmp.cmp(fileprefix+"-NRPymodule.h",fileprefix+".h") == False:
    print("VALIDATION TEST FAILED ON file: "+fileprefix+".h"+".")
    sys.exit(1)

print("VALIDATION TEST PASSED on file: "+fileprefix+".h")
print("### END VALIDATION TESTS ###")
```

    
    
    ### BEGIN VALIDATION TESTS ###
    VALIDATION TEST PASSED on file: SpinWeight_minus2_SphHarmonics/SpinWeight_minus2_SphHarmonics.h
    ### END VALIDATION TESTS ###


<a id='latex_pdf_output'></a>

# Step 6: Output this notebook to $\LaTeX$-formatted PDF \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-SpinWeighted_Spherical_Harmonics.pdf](Tutorial-SpinWeighted_Spherical_Harmonics.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-SpinWeighted_Spherical_Harmonics")
```

    Created Tutorial-SpinWeighted_Spherical_Harmonics.tex, and compiled LaTeX
        file to PDF file Tutorial-SpinWeighted_Spherical_Harmonics.pdf

