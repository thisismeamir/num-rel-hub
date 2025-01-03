**Google Tag Manager Integration**
=====================================

### Overview of Google Tag Manager Integration

In this section, we will explore the integration of Google Tag Manager into a web page.

### Theory Review

#### Introduction to Google Tag Manager

Google Tag Manager is a free tool provided by Google that allows marketers and developers to manage their website's tags (tracking codes) in one place. It simplifies the process of adding and managing tracking codes on your website, without requiring IT support or coding expertise.

```python
# No code needed for this example
```

### Code Implementation


```html
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>
```

### Theory Review

#### How Google Tag Manager Works

Google Tag Manager allows you to add and manage tracking codes on your website through a user-friendly interface. Here's how it works:

1.  **Install the Google Tag Manager container**: Add the Google Tag Manager script to your website by copying and pasting the code into your HTML file.
2.  **Create tags**: Create custom tags using the Google Tag Manager interface, which allows you to choose from a variety of tracking codes (e.g., Google Analytics).
3.  **Trigger tags**: Define triggers that determine when each tag should fire.

### Output


```html
The code above adds and initializes Google Tag Manager on your website.
```

Note: This code snippet must be placed in the `<head>` section of your HTML file to work correctly.

$$\text{Google Tag Manager} \rightarrow \text{add tracking codes} \rightarrow \text{analyze user behavior}$$

### Example Use Case

Suppose we want to track page views on our website using Google Analytics. We can create a custom tag in Google Tag Manager and trigger it when the page loads.

```python
# Example code for creating a custom tag in Google Tag Manager
tag_name = "Page View Tracker"
trigger_type = "Page Load"
```

This will allow us to track page views on our website using Google Analytics.

### Benefits of Google Tag Manager

Google Tag Manager offers several benefits, including:

*   **Easy-to-use interface**: The Google**Equations of General Relativistic Hydrodynamics (GRHD)**
===========================================================

### Overview of GRHD

In this section, we will explore the equations of General Relativistic Hydrodynamics (GRHD).

### Theory Review

#### Introduction to GRHD

General Relativistic Hydrodynamics is a theoretical framework that combines Einstein's theory of general relativity with hydrodynamics. It is used to study the behavior of fluids in extreme environments, such as neutron stars and black holes.

```python
# No code needed for this example
```

### Equations of GRHD

The equations of GRHD can be written as:

$$\nabla_\mu ( \rho u^\mu ) = 0$$

and

$$\nabla_\nu ( T^{\mu\nu} ) = 0$$

where $\rho$ is the energy density, $u^\mu$ is the four-velocity, and $T^{\mu\nu}$ is the stress-energy tensor.

### Code Implementation


```python
import sympy as sp

# Define the variables
x, y, z = sp.symbols('x y z')
t = sp.symbols('t')

# Define the energy density
rho = x**2 + y**2 + z**2 + t**2

# Define the four-velocity
u0 = 1 / (1 - x**2)
ux = -y / u0
uy = x / u0
uz = 0
ut = 0

# Define the stress-energy tensor
T00 = rho * u0**2
T01 = rho * ux * ut
T02 = rho * uy * ut
T03 = rho * uz * ut
T10 = T01
T11 = (rho + p) * ux**2
T12 = (rho + p) * ux * uy
T13 = (rho + p) * ux * uz
T20 = T02
T21 = T12
T22 = (rho + p) * uy**2
T23 = (rho + p) * uy * uz
T30 = T03
T31 = T13
T32 = T23
T33 = (rho + p) * uz**2

# Print the stress-energy tensor
print(T00)
```

### Theory Review

#### How GRHD Works

GRHD**Authors: Zach Etienne & Patrick Nelson**
==========================================

### Overview of the Authors

In this section, we will explore the authors of a particular document or project.

### Theory Review

#### Introduction to Authors

An author is a person who creates written content, such as a book, article, or research paper. In the context of academic and professional work, authors are often recognized for their contributions to their field.

```python
# No code needed for this example
```

### Code Implementation


```markdown
**Authors:**
Zach Etienne & Patrick Nelson
```

### Theory Review

#### Benefits of Recognizing Authors

Recognizing authors is important because it:

*   **Acknowledges their work**: By crediting the author, we acknowledge the time and effort they put into creating the content.
*   **Provides accountability**: Knowing who created the content can help to hold them accountable for any errors or inaccuracies.

### Output


```markdown
The output will be:
**Authors:**
Zach Etienne & Patrick Nelson
```

Note that this is a basic example of recognizing authors. In real-world applications, you may need to include additional information, such as their affiliations or contact details.

### Example Use Case

Suppose we are creating a research paper and want to acknowledge the contributions of our team members. We can use the following code:

```markdown
**Authors:**
Zach Etienne
Patrick Nelson
John Doe
Jane Smith
```

This will ensure that all contributors are recognized for their work.

### Benefits of Citing Authors

Citing authors has several benefits, including:

*   **Academic integrity**: By citing authors, we demonstrate our commitment to academic integrity and respect for the work of others.
*   **Transparency**: Citing authors provides transparency and clarity about who contributed to the content.**General Relativistic Hydrodynamics (GRHD) Notebook**
=====================================================

### Overview of GRHD Notebook

This notebook documents and constructs quantities useful for building symbolic expressions for the equations of general relativistic hydrodynamics.

### Theory Review

#### Introduction to General Relativistic Hydrodynamics

General Relativistic Hydrodynamics is a theoretical framework that combines Einstein's theory of general relativity with hydrodynamics. It is used to study the behavior of fluids in extreme environments, such as neutron stars and black holes.

```python
# Import necessary libraries
import sympy as sp
```

### Code Implementation


```markdown
**Notebook Status:** <font color='orange'><b> Self-Validated </b></font>

**Validation Notes:**
This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#validation-notes).

**Quantities Used in this Notebook:**

*   `rho`: energy density
*   `u0`: four-velocity (time component)
*   `ux`, `uy`, `uz`: four-velocity components (space components)

```python
# Define the variables
x, y, z = sp.symbols('x y z')
t = sp.symbols('t')

# Define the energy density
rho = x**2 + y**2 + z**2 + t**2

# Define the four-velocity components
u0 = 1 / (1 - x**2)
ux = -y / u0
uy = x / u0
uz = 0
```

### Theory Review

#### Benefits of Using SymPy for GRHD

Using SymPy to build symbolic expressions for GRHD has several benefits, including:

*   **Easy to read and write**: SymPy expressions are easy to read and write, making it easier to understand and modify the code.
*   **Fast computation**: SymPy can perform fast computations on large expressions.

### Output


```markdown
The output will be:
`rho`
`u0`
`ux`
`uy`
`uz`

```

Note that this is a basic example of using SymPy for GRHD. In real-world applications, you may need to use more complex expressions and functions.

### Example Use Case

Suppose we want to build a symbolic expression for the stress-energy tensor in GRHD. We can use the following code:

```python
# Define**Code Validation**
=====================

### Overview of Code Validation

In this section, we will explore the process of validating code.

### Theory Review

#### Introduction to Code Validation

Code validation is an essential step in ensuring that our code works correctly and produces accurate results. It involves testing our code against expected outputs or behaviors.

```python
# No code needed for this example
```

### Code Implementation


```markdown
**Code Validation**
-------------------

*   **Unit tests**: Unit tests are small, focused tests that check the correctness of a single function or method.
*   **Integration tests**: Integration tests verify that multiple components work together as expected.

**Additional validation tests may have been performed, but are as yet, undocumented. (TODO)**
```

### Theory Review

#### Types of Code Validation

There are several types of code validation, including:

*   **Unit testing**: Unit testing involves writing small, focused tests to verify the correctness of individual functions or methods.
*   **Integration testing**: Integration testing checks that multiple components work together as expected.

### Output


```markdown
The output will be:
`# Code Validation`

`## Unit tests`
`### Integration tests`

`## Additional validation tests may have been performed, but are as yet, undocumented. (TODO)`
```

Note that this is a basic example of code validation. In real-world applications, you may need to perform more extensive testing and validation.

### Example Use Case

Suppose we want to validate our code for calculating the stress-energy tensor in general relativity. We can write unit tests to verify the correctness of individual functions and methods:

```python
# Define a function to calculate the stress-energy tensor
def calculate_stress_energy_tensor(rho, u0, ux, uy, uz):
    # Calculate the stress-energy tensor components
    T00 = rho * u0**2
    T01 = rho * ux * u0
    T02 = rho * uy * u0
    T03 = rho * uz * u0
    
    return T00, T01, T02, T03

# Write unit tests to verify the correctness of the function
def test_calculate_stress_energy_tensor():
    # Test case 1: Verify that the stress-energy tensor components are calculated correctly for a specific input
    rho = 2.0
    u0 = 1.0
    ux = -1.0
    uy = 1**General Relativistic Hydrodynamics Equations**
=====================================================

### Overview of General Relativistic Hydrodynamics Equations

In this section, we will explore the equations of general relativistic hydrodynamics in conservative form.

### Theory Review

#### Introduction to General Relativistic Hydrodynamics

General Relativistic Hydrodynamics is a theoretical framework that combines Einstein's theory of general relativity with hydrodynamics. It is used to study the behavior of fluids in extreme environments, such as neutron stars and black holes.

```python
# No code needed for this example
```

### Equations of General Relativistic Hydrodynamics

The equations of general relativistic hydrodynamics can be written as:

$$\partial_t \rho_* + \partial_j \left(\rho_* v^j\right) = 0$$

$$\partial_t \tilde{\tau} + \partial_j \left(\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j \right) = s$$

$$\partial_t \tilde{S}_i + \partial_j \left(\alpha \sqrt{\gamma} T^j{}_i \right) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i},$$

where we assume $T^{\mu\nu}$ is the stress-energy tensor of a perfect fluid:

$$
T^{\mu\nu} = \rho_0 h u^{\mu} u^{\nu} + P g^{\mu\nu},
$$

the $s$ source term is given in terms of ADM quantities via

$$
s = \alpha \sqrt{\gamma}\left[\left(T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \right)K_{ij}
- \left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha \right],
$$

and 

$$
v^j = \frac{u^j}{u^0} \\
\rho_* = \alpha\sqrt{\gamma} \rho_0 u^0 \\
h = 1 + \epsilon + \frac{P}{\rho_0}.
$$

### Code Implementation


```python
import sympy as sp**A Note on Notation**
=====================

### Overview of Notation Conventions

In this section, we will explore the notation conventions used in NRPy+.

### Theory Review

#### Introduction to Notation Conventions

NRPy+ uses a specific set of notation conventions to represent mathematical expressions. These conventions are essential for understanding and working with NRPy+ code.

```python
# No code needed for this example
```

### Greek Indices vs. Latin Indices

In NRPy+, Greek indices refer to four-dimensional quantities, while Latin indices refer to three-dimensional quantities.

*   **Greek Indices**: Represent four-dimensional quantities where the zeroth component indicates the temporal (time) component.
*   **Latin Indices**: Represent three-dimensional quantities, which can be counterintuitive since Python lists start indexing from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first spatial direction.

### Code Implementation


```python
# Import necessary libraries
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')
T4EMUU = ixp.zerorank2(DIM=4)

# Iterate over Greek indices (mu, nu)
for mu in range(4):
    for nu in range(4):
        # Calculate first term of b^2 u^\mu u^\nu
        T4EMUU[mu, nu] += sp.symbols('b')**2 * sp.symbols('u_mu') * sp.symbols('u_nu')
```

### Theory Review

#### Why Use Greek Indices for Four-Dimensional Quantities?

Greek indices are used to represent four-dimensional quantities because they allow for a more natural and concise representation of mathematical expressions. In particular, the zeroth component of Greek indices always corresponds to the temporal (time) component.

$$
b^2 u^\mu u^\nu = \left( b^2 u^0 u^0 + b^2 u^1 u^1 + b^2 u^2 u^2 + b^2 u^3 u^3 \right)
$$

### Output


```python
The output will be:
`b**2*u_mu*u_nu`
```

Note that this is a basic example of using Greek indices in NRPy+. In real-world applications, you may need to work with more complex mathematical expressions.

### Example Use Case

Suppose we want to calculate the stress-energy**Mixed Greek and Latin Indices**
================================

### Overview of Mixed Greek and Latin Indices

In this section, we will explore how to handle mixed Greek and Latin indices.

### Theory Review

#### Introduction to Mixed Greek and Latin Indices

When working with mixed Greek and Latin indices, it is essential to understand the conventions used in NRPy+. Specifically, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one.

```python
# No code needed for this example
```

### Handling Mixed Greek and Latin Indices

To handle mixed Greek and Latin indices, we need to follow specific rules:

*   **Greek Index in a Three-Vector**: When a Greek index appears in a three-vector, it should be decremented.
*   **Latin Index in a Four-Vector**: When a Latin index appears in a four-vector, it should be incremented.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')
gammaDD = ixp.zerorank2(DIM=3)
betaU = ixp.zerorank1(DIM=3)

# Calculate $\beta_i = \gamma_{ij} \beta^j$
for i in range(3):
    for j in range(3):
        gammaDD[i][j] += sp.symbols('gamma_ij') * betaU[j]
```

### Theory Review

#### Why Offset One Set of Indices by One?

When working with mixed Greek and Latin indices, it is necessary to offset one set of indices by one because Python lists start indexing from 0. By following this convention, we can ensure that our mathematical expressions are represented correctly.

$$
\beta_i = \gamma_{ij} \beta^j = \left( \gamma_{00} \beta^0 + \gamma_{01} \beta^1 + \gamma_{02} \beta^2 + \gamma_{03} \beta^3 \right)
$$

### Output


```python
The output will be:
`beta_i`
```

Note that this is a basic example of handling mixed Greek and Latin indices. In real-world applications, you may need to work with more complex mathematical expressions.

### Example Use Case

Suppose we want to calculate the expression $\frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g**Expression: $$\alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2$$**
====================================================================

### Overview of the Expression

In this section, we will explore the expression $\alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2$.

### Theory Review

#### Introduction to the Expression

The expression $$\alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2$$ represents a term in the stress-energy tensor, where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the metric, $T^{\mu \nu}_{\rm EM}$ is the electromagnetic stress-energy tensor, and $\partial_i g_{\mu \nu}$ is the derivative of the metric with respect to the $i^{th}$ coordinate.

```python
# No code needed for this example
```

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')
alpsqrtgam = sp.symbols('alpha sqrt(gamma)')
T4EMUU = ixp.zerorank2(DIM=4)
g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)

# Calculate expression
for i in range(3):
    for mu in range(4):
        for nu in range(4):
            S_tilde_rhsD[i] += alpsqrtgam * T4EMUU[mu][nu] * g4DD_zerotimederiv_dD[mu][nu][i+1] / 2
```

### Theory Review

#### Why Use the Expression?

The expression $\alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2$ is used to represent a term in the stress-energy tensor. This term represents the interaction between the electromagnetic field and the metric.

$$
T^{\mu \nu}_{\rm EM} = \left( F^{\mu \rho} F^\nu_\rho - \frac{1}{4} g^{\mu \nu} F_{\**Table of Contents**
=====================

### Overview of the Notebook

This notebook is organized into several sections, each containing a set of functions that construct various families of quantities.

$$\label{toc}$$

### Functions for Constructing Quantities

Each section below corresponds to a function that constructs a specific family of quantities.


#### Step 1: Constructing Basic Quantities

In this step, we will define the basic quantities used throughout the notebook. These include:

*   `alpha`
*   `betaU`
*   `gammaDD`
*   `g4DD`

```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Define basic quantities
alpha = ixp.symbols('alpha')
betaU = ixp.zerorank1(DIM=3)
gammaDD = ixp.zerorank2(DIM=3)
g4DD = ixp.zerorank4(DIM=4)
```

#### Step 2: Constructing Stress-Energy Tensor

In this step, we will define the stress-energy tensor using the basic quantities defined in the previous section.

```python
# Define stress-energy tensor
T4EMUU = ixp.zerorank2(DIM=4)

# Calculate stress-energy tensor components
for mu in range(4):
    for nu in range(4):
        T4EMUU[mu][nu] += alpha**2 * g4DD[mu][nu]
```

#### Step 3: Constructing Metric Derivatives

In this step, we will define the derivatives of the metric using the basic quantities defined in the previous sections.

```python
# Define metric derivatives
g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)

# Calculate metric derivative components
for mu in range(4):
    for nu in range(4):
        for i in range(3):
            g4DD_zerotimederiv_dD[mu][nu][i+1] += alpha * g4DD[mu][nu]
```

### Theory Review

#### Why Organize the Notebook in this Way?

The notebook is organized into several sections to make it easier to follow and understand. Each section corresponds to a specific function that constructs a family of quantities.

$$
\text{Step 1} \rightarrow \text{Basic Quantities}
$$**Importing Modules**
=====================

### Overview of Module Imports

In this section, we will discuss the necessary imports for NRPy+ and Python.

### Theory Review

#### Introduction to Module Imports

Module imports are essential in programming as they allow us to use pre-existing code libraries to perform complex tasks. In the context of NRPy+, we need to import modules that provide functionalities for numerical relativity.

```python
# No code needed for this example
```

### Code Implementation


```python
import sympy as sp
import numpy as np

from nrpy import *

# Define variables
ixp = sp.IndexedBase('ixp')

# Import necessary NRPy+ modules
nrpy_modules = ['nrpy']
for module in nrpy_modules:
    exec(f'import {module}')

# Import additional Python modules
python_modules = ['numpy', 'scipy']
for module in python_modules:
    exec(f'import {module}')
```

### Theory Review

#### Why Use Module Imports?

Module imports are used to reuse code written by others. This allows us to focus on the main task at hand without having to reinvent the wheel.

*   **Modularity**: Modular code is easier to understand, maintain, and modify.
*   **Reusability**: Imported modules can be reused in other projects or contexts.

### Output


```python
The output will be:
`sp`
`np`
```

Note that this is a basic example of importing NRPy+ and Python modules. In real-world applications, you may need to import additional modules depending on the specific requirements of your project.

### Example Use Case

Suppose we want to use the `numpy` module for numerical computations:

```python
# Import numpy module
import numpy as np

# Perform a simple numerical computation
x = np.array([1, 2, 3])
y = x**2
print(y)
```

This will output: `[1 4 9]`**Defining the Stress-Energy Tensor**
=====================================

### Overview of the Stress-Energy Tensor

In this section, we will define the stress-energy tensor $T^{\mu\nu}$ and its corresponding components.

$$
T^{\mu\nu} = \left( F^{\mu\rho} F^\nu_\rho - \frac{1}{4} g^{\mu\nu} F_{\rho\sigma} F^{\rho\sigma} \right)
$$

### Code Implementation


```python
import sympy as sp
import numpy as np

# Define variables
ixp = sp.IndexedBase('ixp')

# Define the stress-energy tensor components
T4UU = ixp.zerorank2(DIM=4)

for mu in range(4):
    for nu in range(4):
        T4UU[mu][nu] += (F_mu_rho * F_nu_rho - 0.25 * g_mu_nu * F_rho_sigma * F_rho_sigma)
```

### Theory Review

#### What is the Stress-Energy Tensor?

The stress-energy tensor $T^{\mu\nu}$ represents the distribution of energy and momentum in a system.

*   **Energy-Momentum**: The stress-energy tensor encodes information about the energy and momentum of a system.
*   **Conservation Laws**: The stress-energy tensor satisfies conservation laws, ensuring that energy and momentum are conserved in a system.

### Output


```python
The output will be:
`T^{\mu\nu}`
```

Note that this is a basic example of defining the stress-energy tensor. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the energy density $T^{00}$:

```python
# Calculate the energy density
energy_density = T4UU[0][0]
print(energy_density)
```

This will output: `[expression]`**Writing Conservative Variables**
================================

### Overview of Conservative Variables

In this section, we will write the conservative variables in terms of the primitive variables.

### Theory Review

#### What are Conservative Variables?

Conservative variables are a set of quantities that describe the state of a system in a way that is invariant under certain transformations. In the context of fluid dynamics and numerical relativity, conservative variables can be used to describe the evolution of a system over time.

*   **Invariance**: Conservative variables are invariant under certain transformations, making them useful for describing systems that undergo changes due to external factors.
*   **Evolution**: The use of conservative variables allows us to study the evolution of a system over time in a more robust and accurate way.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute sqrt(gamma) DET
sqrtgamDET = ixp.symbols('sqrt(gamma)') * sp.det(ixp.zerorank2(DIM=3))

# Compute rho_star
rho_star = ixp.symbols('rho_star')
for i in range(3):
    rho_star += ixp.symbols('u_i')**2

# Compute tau_tilde
tau_tilde = ixp.symbols('tau_tilde')

# Compute S_tildeD
S_tildeD = ixp.zerorank1(DIM=3)
for i in range(3):
    S_tildeD[i] += ixp.symbols('u_i') * rho_star
```

### Theory Review

#### Why Write Conservative Variables?

Writing conservative variables allows us to describe the state of a system in a more robust and accurate way. This is particularly useful when studying systems that undergo changes due to external factors.

*   **Accuracy**: Writing conservative variables ensures that our description of a system is accurate and takes into account all relevant physical processes.
*   **Robustness**: The use of conservative variables makes our analysis more robust and less prone to errors.

### Output


```python
The output will be:
`sqrt(gamma) DET`
`rho_star`
`tau_tilde`
```

Note that this is a basic example of writing conservative variables. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the energy density rho:

```python**Defining Fluxes for GRHD Equations**
=====================================

### Overview of GRHD Equations

In this section, we will define the fluxes for the General Relativity Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

The GRHD equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Define fluxes
fluxes = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        fluxes[mu][nu] += (rho_star * u_nu) + (tau_tilde * u_mu)

# Compute conservative variables
cons_vars = ixp.zerorank1(DIM=4)
for i in range(3):
    cons_vars[i] += rho_star * u_i

# Compute fluxes for GRHD equations
flux_grhd = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        flux_grhd[mu][nu] += (rho_star * u_nu) + (tau_tilde * u_mu)
```

### Theory Review

#### What are the Fluxes?

The fluxes are a set of quantities that describe the transport of energy and momentum across the boundaries of a system.

*   **Energy Transport**: The fluxes capture the transport of energy due to the motion of the fluid.
*   **Momentum Transport**: The fluxes also describe the transport of momentum, which is essential for studying systems in strong-field gravity.

### Output


```python
The output will be:
`flux_grhd`
```

Note that this is a basic example of defining fluxes for GRHD equations. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the energy flux:

```python
**Defining $\rho_*$ Flux Term**
=============================

### Overview of $\rho_*$ Flux Term

In this section, we will define the $\rho_*$ flux term for the General Relativity Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

The GRHD equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute vU from u4U without speed limit
vU = ixp.zerorank1(DIM=3)
for i in range(3):
    vU[i] += (u_0 * u_i) / u_0**2

# Compute rho* flux term
rho_star_fluxU = ixp.zerorank1(DIM=3)
for i in range(3):
    rho_star_fluxU[i] += rho_star * vU[i]
```

### Theory Review

#### What is the $\rho_*$ Flux Term?

The $\rho_*$ flux term captures the transport of fluid mass across the boundaries of a system.

*   **Fluid Mass Transport**: The $\rho_*$ flux term describes how fluid mass moves from one region to another.
*   **Conservation Law**: The $\rho_*$ flux term is an essential component of the GRHD equations, which describe the conservation of fluid mass and energy.

### Output


```python
The output will be:
`rho_star_fluxU`
```

Note that this is a basic example of defining the $\rho_*$ flux term. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the fluid mass transport:

```python
# Calculate fluid mass transport
fluid_mass_transport = rho_star_fluxU[0]
print(fluid_mass_transport)
```

This will output: `[expression]`**Defining $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms**
=====================================================

### Overview of $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms

In this section, we will define the $\tilde{\tau}$ and $\tilde{S}_i$ flux terms for the General Relativity Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

The GRHD equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute tau_tilde flux term
tau_tilde_fluxU = ixp.zerorank1(DIM=3)
for i in range(3):
    tau_tilde_fluxU[i] += (tau_tilde * u_i)

# Compute tilde_S_i flux term
tilde_S_i_fluxUD = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        tilde_S_i_fluxUD[mu][nu] += (alpha * sqrt(gamma) * T_nu_mu)
```

### Theory Review

#### What are the $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms?

The $\tilde{\tau}$ and $\tilde{S}_i$ flux terms capture the transport of energy and momentum across the boundaries of a system.

*   **Energy Transport**: The $\tilde{\tau}$ flux term describes how energy moves from one region to another.
*   **Momentum Transport**: The $\tilde{S}_i$ flux terms describe how momentum moves in different directions.
*   **Conservation Law**: The $\tilde{\tau}$ and $\tilde{S}_i$ flux terms are essential components of the GRHD equations, which describe the conservation of energy and momentum.

### Output


```python
The output will be:
`tau_tilde_fluxU**Defining Source Terms for GRHD Equations**
=====================================

### Overview of Source Terms

In this section, we will define the source terms that appear on the right-hand sides (RHSs) of the General Relativity Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

The GRHD equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute source term for density evolution equation
source_term_density = ixp.symbols('source_term_density')
for i in range(3):
    source_term_density += (alpha * sqrt(gamma) * rho_star)

# Compute source term for momentum evolution equations
source_terms_momentum = ixp.zerorank1(DIM=3)
for i in range(3):
    source_terms_momentum[i] += (alpha * sqrt(gamma) * tau_tilde)
```

### Theory Review

#### What are the Source Terms?

The source terms capture the effects of various physical processes on the evolution of the fluid.

*   **Density Evolution Equation**: The source term for density evolution equation represents the rate at which density changes due to various physical processes.
*   **Momentum Evolution Equations**: The source terms for momentum evolution equations represent the rates at which momentum changes due to various physical processes.
*   **Curved Spacetime Effects**: The source terms take into account the effects of curved spacetime on the fluid, such as gravitational waves and frame-dragging.

### Output


```python
The output will be:
`source_term_density`
`source_terms_momentum`

```

Note that this is a basic example of defining source terms for GRHD equations. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the rate at which density changes due to various physical processes:

```python
# Calculate rate of change of density
rate**Defining $s$ Source Term**
==========================

### Overview of $s$ Source Term

In this section, we will define the $s$ source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to GRHD Equations

The General Relativity Hydrodynamics (GRHD) equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute s source term
s_source_term = ixp.symbols('s')
for i in range(3):
    s_source_term += (alpha * sqrt(gamma) * T0i)
```

### Theory Review

#### What is the $s$ Source Term?

The $s$ source term captures the effects of various physical processes on the evolution of the fluid.

*   **Energy Momentum Transport**: The $s$ source term represents the rate at which energy and momentum are transported across the boundaries of a system.
*   **Curved Spacetime Effects**: The $s$ source term takes into account the effects of curved spacetime on the fluid, such as gravitational waves and frame-dragging.

### Output


```python
The output will be:
`s_source_term`

```

Note that this is a basic example of defining the $s$ source term. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the rate at which energy and momentum are transported across the boundaries of a system:

```python
# Calculate rate of energy-momentum transport
rate_energy_momentum_transport = s_source_term
print(rate_energy_momentum_transport)
```

This will output: `[expression]`**Defining Source Term on RHS of $\tilde{S}_i$ Equation**
=====================================================

### Overview of $\tilde{S}_i$ Equation

In this section, we will define the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to GRHD Equations

The General Relativity Hydrodynamics (GRHD) equations describe the evolution of a fluid in a curved spacetime. The equations consist of a set of conservation laws that govern the behavior of the fluid.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute source term on RHS of tilde_S_i equation
source_term_tilde_S_i = ixp.symbols('source_term_tilde_S_i')
for mu in range(4):
    for nu in range(4):
        source_term_tilde_S_i += (alpha * sqrt(gamma) * T_nu_mu)
```

### Theory Review

#### What is the Source Term on RHS of $\tilde{S}_i$ Equation?

The source term on RHS of $\tilde{S}_i$ equation captures the effects of various physical processes on the evolution of the fluid.

*   **Energy Momentum Transport**: The source term represents the rate at which energy and momentum are transported across the boundaries of a system.
*   **Curved Spacetime Effects**: The source term takes into account the effects of curved spacetime on the fluid, such as gravitational waves and frame-dragging.

### Output


```python
The output will be:
`source_term_tilde_S_i`

```

Note that this is a basic example of defining the source term on RHS of $\tilde{S}_i$ equation. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the rate at which energy and momentum are transported across the boundaries of a system:

```python
# Calculate rate of energy-momentum transport
rate_energy_momentum_transport =**Computing $g_{\mu\nu,i}$ in Terms of ADM Quantities and Their Derivatives**
=====================================================================

### Overview of Computing $g_{\mu\nu,i}$

In this section, we will compute the derivatives of the four-metric $g_{\mu\nu}$ with respect to the coordinate $\mathbf{x}^i$, denoted as $g_{\mu\nu,i}$.

### Theory Review

#### Introduction to ADM Quantities and Their Derivatives

The Arnowitt-Deser-Misner (ADM) formalism is a mathematical framework used to describe the evolution of general relativity. The ADM quantities are defined in terms of the three-metric $\gamma_{ij}$ and its derivatives with respect to time.

*   **Three-Metric**: The three-metric $\gamma_{ij}$ represents the metric on the spatial hypersurface.
*   **Derivatives of Three-Metric**: The derivatives of the three-metric with respect to time are used to compute the ADM quantities.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute g4DD_zerotimederiv_dD
g4DD_zerotimederiv_dD = ixp.symbols('g4DD_zerotimederiv_dD')
for mu in range(4):
    for nu in range(4):
        g4DD_zerotimederiv_dD[mu][nu] += (1/2) * (gamma_inv[mu, lambda] * gamma_inv[nu, sigma] * gamma_lambda_sigma)
```

### Theory Review

#### What is the Expression for $g_{\mu\nu,i}$?

The expression for $g_{\mu\nu,i}$ can be obtained by taking the derivative of the four-metric with respect to the coordinate $\mathbf{x}^i$.

*   **Derivative of Four-Metric**: The derivative of the four-metric with respect to the coordinate $\mathbf{x}^i$ is given by $g_{\mu\nu,i}$.
*   **ADM Quantities and Their Derivatives**: The expression for $g_{\mu\nu,i}$ involves the ADM quantities and their derivatives with respect to time.

### Output


```python
The output will be:
`g4DD_zerotimederiv_d**Computing the Source Term of the $\tilde{S}_i$ Equation**
=====================================================

### Overview of Computing the Source Term of the $\tilde{S}_i$ Equation

In this section, we will compute the source term of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to the $\tilde{S}_i$ Equation

The $\tilde{S}_i$ equation is a component of the General Relativity Hydrodynamics (GRHD) equations, which describe the evolution of a fluid in a curved spacetime.

*   **Conservation Laws**: The GRHD equations are based on the conservation of energy and momentum.
*   **Curved Spacetime**: The GRHD equations take into account the curvature of spacetime, which is essential for describing systems in strong-field gravity.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Compute S_tilde_source_termD
S_tilde_source_termD = ixp.symbols('S_tilde_source_termD')
for mu in range(4):
    for nu in range(4):
        S_tilde_source_termD[mu][nu] += (1/2) * (alpha * sp.sqrt(gamma) * T_mu_nu)
```

### Theory Review

#### What is the Expression for the Source Term of the $\tilde{S}_i$ Equation?

The source term of the $\tilde{S}_i$ equation can be computed using the formula:

$$\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}$$

*   **ADM Quantities and Their Derivatives**: The source term involves the ADM quantities and their derivatives with respect to time.
*   **Curved Spacetime Effects**: The source term takes into account the effects of curved spacetime on the fluid, such as gravitational waves and frame-dragging.

### Output


```python
The output will be:
`S_tilde_source_termD`

```

Note that this is a basic example of computing the source term of the $\tilde{S}_i$ equation. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to calculate the**Conversion of $v^i$ to $u^\mu$**
================================

### Overview of Conversion of $v^i$ to $u^\mu$

In this section, we will discuss the conversion of the fluid velocity $v^i$ to the four-velocity $u^\mu$. This is an important step in converting the primitive variables to conservative variables.

### Theory Review

#### Introduction to Fluid Velocity and Four-Velocity

*   **Fluid Velocity**: The fluid velocity $v^i$ represents the velocity of the fluid with respect to a coordinate frame.
*   **Four-Velocity**: The four-velocity $u^\mu$ is a four-vector that represents the velocity of an object in spacetime.

### Code Implementation


```python
import sympy as sp

# Define variables
ixp = sp.IndexedBase('ixp')

# Convert vU to u4U with rescaling (Valencia method)
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit():
    u4U = ixp.symbols('u4U')
    for mu in range(4):
        u4U[mu] += vU * (1 - (vU**2) / c**2)**(-0.5)
    
    return u4U

# Convert vU to u4U with rescaling (Valencia method)
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit():
    u4U = ixp.symbols('u4U')
    for mu in range(4):
        u4U[mu] += vU * (1 - (vU**2) / c**2)**(-0.5)
    
    return u4U

# Compute u4U
u4U = u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit()
```

### Theory Review

#### What is the Expression for $u^\mu$?

The expression for $u^\mu$ can be obtained by applying the Valencia method, which involves rescaling the fluid velocity $v^i$.

*   **Rescaling**: The Valencia method introduces a rescaling factor that depends on the speed of light $c$ and the fluid velocity $v^i$.
*   **Four-Velocity Expression**: The resulting expression for $u^\mu$**Declaring Variables and Constructing GRHD Equations**
=====================================================

### Overview of Declaring Variables and Constructing GRHD Equations

In this section, we will declare the ADM and hydrodynamical input variables, and construct the General Relativity Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to ADM and Hydrodynamical Input Variables

*   **ADM Quantities**: The Arnowitt-Deser-Misner (ADM) quantities describe the properties of the spacetime.
*   **Hydrodynamical Quantities**: The hydrodynamical quantities describe the properties of the fluid.

### Code Implementation


```python
import sympy as sp

# Declare ADM and hydrodynamical input variables
ADM_vars = {
    'gamma': ixp.symbols('gamma'),
    'alpha': ixp.symbols('alpha'),
    'beta_i': ixp.symbols('beta_i')
}

hydro_vars = {
    'rho': ixp.symbols('rho'),
    'u_i': ixp.symbols('u_i'),
    'P': ixp.symbols('P')
}

# Construct GRHD equations
def construct_GRHD_equations(ADM_vars, hydro_vars):
    # Density evolution equation
    density_evolution = sp.Eq(ixp.Derivative(rho, ADM_vars['t']), -ixp.Derivative(alpha * rho, ADM_vars['t']))
    
    # Momentum evolution equations
    momentum_evolution = [sp.Eq(ixp.Derivative(u_i, ADM_vars['t']), -ixp.Derivative(alpha * u_i, ADM_vars['t'])) for i in range(3)]
    
    return density_evolution, momentum_evolution

# Construct GRHD equations
density_evolution, momentum_evolution = construct_GRHD_equations(ADM_vars, hydro_vars)
```

### Theory Review

#### What are the ADM and Hydrodynamical Input Variables?

The ADM quantities describe the properties of the spacetime, while the hydrodynamical quantities describe the properties of the fluid.

*   **ADM Quantities**: The ADM quantities include the metric $\gamma$, the lapse function $\alpha$, and the shift vector $\beta_i$.
*   **Hydrodynamical Quantities**: The hydrodynamical quantities include the density $\rho$, the velocity $u_i$, and the pressure $P$.

### Output


```python
**Code Validation Against GRHD Equations**
=====================================

### Overview of Code Validation Against GRHD Equations

In this section, we will validate our code by comparing it with the existing `GRHD.equations` module in the NRPy+ framework.

### Theory Review

#### Introduction to GRHD Equations and NRPy+

*   **GRHD Equations**: The General Relativity Hydrodynamics (GRHD) equations describe the evolution of a fluid in a curved spacetime.
*   **NRPy+ Framework**: The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.

### Code Implementation


```python
import numpy as np
import sympy as sp

# Load GRHD equations from NRPy+
from nrpy import load_grhd_equations

# Define variables
ixp = sp.IndexedBase('ixp')

# Get GRHD equations from NRPy+
GRHD_equations = load_grhd_equations()

# Compare our code with the GRHD equations from NRPy+
for equation in GRHD_equations:
    if not sp.simplify(equation) == sp.simplify(compute_GRHD_equations()):
        print(f"Error: {equation} does not match the corresponding equation in NRPy+")
```

### Theory Review

#### What is Code Validation?

Code validation involves checking that our code produces the same results as a known and trusted implementation.

*   **Known Implementation**: In this case, we are comparing our code with the existing `GRHD.equations` module in the NRPy+ framework.
*   **Trustworthy Results**: If our code produces the same results as the known implementation, then it is likely to be correct.

### Output


```python
The output will be:
```

Note that this is a basic example of code validation. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to check if our implementation of the GRHD equations matches the existing implementation in NRPy+:

```python
# Check if our implementation of the GRHD equations matches the existing implementation in NRPy+
if sp.simplify(compute_GRHD_equations()) == sp.simplify(GRHD_equations):
    print("Our code is correct!")
else:
    print("Error: Our code does not match the expected result.")
```

This will output `Our code**Outputting This Notebook to a LaTeX-formatted PDF File**
=====================================================

### Overview of Outputting This Notebook to a LaTeX-formatted PDF File

In this section, we will output the current notebook to a $\LaTeX$-formatted PDF file.

### Theory Review

#### Introduction to LaTeX Formatting

*   **$\LaTeX$ Formatting**: The $\LaTeX$ formatting language is used to create mathematical documents.
*   **PDF Output**: We can use tools like Sphinx or Jupyter Book to output our notebook as a $\LaTeX$-formatted PDF file.

### Code Implementation


```python
import nbconvert

# Create a new HTMLExporter object
html_exporter = nbconvert.HTMLExporter()

# Convert the notebook to an HTML string
html_string, resources = html_exporter.from_notebook_node(nb)

# Save the HTML string as a LaTeX-formatted PDF file
with open('output.pdf', 'w') as f:
    f.write(html_string)
```

### Theory Review

#### What is $\LaTeX$ Formatting?

$\LaTeX$ formatting is a typesetting system that allows us to create mathematical documents.

*   **Mathematical Documents**: We can use $\LaTeX$ to create documents with mathematical equations, formulas, and other special characters.
*   **PDF Output**: By using tools like Sphinx or Jupyter Book, we can output our notebook as a $\LaTeX$-formatted PDF file.

### Output


```python
The output will be:
```

Note that this is a basic example of outputting a notebook to a $\LaTeX$-formatted PDF file. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to output our notebook as a $\LaTeX$-formatted PDF file:

```python
# Output the notebook as a LaTeX-formatted PDF file
nbconvert(nb, 'output.pdf', {'pdf_exporter': 'sphinx'})
```

This will output `output.pdf`**Step 1: Importing Needed Modules**
=====================================

### Overview of Importing Needed Modules

In this step, we will import the necessary NRPy+ and Python modules required for our project.

### Theory Review

#### Introduction to NRPy+

*   **NRPy+**: The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.
*   **Modules**: We need to import various modules from NRPy+, including `nrpy` and `utils`.

### Code Implementation


```python
import os
import sys

# Import necessary NRPy+ modules
from nrpy import *

# Import necessary Python modules
import numpy as np
```

### Theory Review

#### What are the Necessary Modules?

We need to import various modules from NRPy+, including `nrpy` and `utils`.

*   **NRPy+:** We need to import the `nrpy` module, which contains functions for numerical relativity.
*   **numpy:** We also need to import the `numpy` module, which is a Python library for efficient numerical computation.

### Output


```python
The output will be:
```

Note that this is a basic example of importing necessary modules. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to import the `nrpy` module:

```python
# Import the nrpy module
import nrpy as nr
```

This will import the `nrpy` module and make its functions available for use in our code.**Table of Contents**
=====================

### Overview of Table of Contents

This document serves as a guide for the NRPy+ tutorial. The table of contents below provides an overview of the topics covered in this tutorial.

### Theory Review

#### Introduction to NRPy+

*   **NRPy+:** The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.
*   **Modules:** We need to import various modules from NRPy+, including `nrpy` and `utils`.

### Code Implementation


```python
import os
import sys

# Import necessary NRPy+ modules
from nrpy import *

# Import necessary Python modules
import numpy as np
```

### Theory Review

#### What are the Necessary Modules?

We need to import various modules from NRPy+, including `nrpy` and `utils`.

*   **NRPy+:** We need to import the `nrpy` module, which contains functions for numerical relativity.
*   **numpy:** We also need to import the `numpy` module, which is a Python library for efficient numerical computation.

### Output


```python
The output will be:
```

Note that this is a basic example of importing necessary modules. In real-world applications, you may need to consider additional components or modifications depending on the specific requirements of your project.

### Example Use Case

Suppose we want to import the `nrpy` module:


```python
# Import the nrpy module
import nrpy as nr
```

This will import the `nrpy` module and make its functions available for use in our code.

## Step 1: Import needed NRPy+ & Python modules

[Back to [top](\#toc)\]**Step 1: Importing Core NRPy+ Modules**
=====================================

### Overview of Importing Core NRPy+ Modules

In this step, we will import the necessary core modules from the NRPy+ framework.

### Code Implementation


```python
from outputC import nrpyAbs
```

### Theory Review

#### Introduction to NRPy+

*   **NRPy+:** The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.
*   **Core Modules:** The core modules are the fundamental building blocks of the NRPy+ framework.

### Code Explanation


```python
from outputC import nrpyAbs
```

This code imports the `nrpyAbs` module from the `outputC` package. This module contains various functions and classes for numerical relativity, including:

*   **Numerical Integration:** Functions for performing numerical integration of differential equations.
*   **Root Finding:** Functions for finding roots of polynomials and other mathematical expressions.

### Theory Review

#### What are the Core Modules?

The core modules in NRPy+ include:

*   **nrpyAbs:** This module contains functions for numerical integration, root finding, and other mathematical operations.
*   **nrpyFuncs:** This module contains functions for creating and manipulating numerical functions.
*   **nrpyMath:** This module contains functions for performing mathematical operations.

### Example Use Cases

Suppose we want to use the `nrpyAbs` module to perform numerical integration:


```python
import nrpyAbs as na

# Define a function to integrate
def integrand(x):
    return x**2 + 1

# Perform numerical integration
result = na.integrate(integrand, 0, 1)

print(result)
```

This will output the result of the numerical integration.

## Step 2: Importing Additional NRPy+ Modules


```python
from outputC import nrpyFuncs
```

This code imports the `nrpyFuncs` module from the `outputC` package. This module contains functions for creating and manipulating numerical functions.

## Step 3: Importing Core Python Modules


```python
import numpy as np
```

This code imports the NumPy library, which is a powerful tool for efficient numerical computation in Python.**NRPy+: Core C Code Output Module**
=====================================

### Overview of the NRPy+ Core C Code Output Module

In this section, we will discuss the core C code output module of the NRPy+ framework.

### Theory Review

#### Introduction to NRPy+

*   **NRPy+:** The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.
*   **Core Modules:** The core modules are the fundamental building blocks of the NRPy+ framework.

### Code Implementation


```python
import NRPy_param_funcs as par
```

### Theory Review

#### What is the Core C Code Output Module?

The core C code output module is responsible for generating the C code that implements the numerical relativity algorithms in the NRPy+ framework.

*   **C Code Generation:** The core C code output module generates the C code that implements the numerical relativity algorithms.
*   **Algorithm Implementation:** The generated C code implements the numerical relativity algorithms, including numerical integration and root finding.

### Code Explanation


```python
import NRPy_param_funcs as par
```

This code imports the `NRPy_param_funcs` module, which is a collection of functions for setting up parameters in the NRPy+ framework. The `par` variable is an alias for this module.

### Theory Review

#### How Does the Core C Code Output Module Work?

The core C code output module works by:

*   **Parsing Input Parameters:** Parsing the input parameters and setting them up in the NRPy+ framework.
*   **Generating C Code:** Generating the C code that implements the numerical relativity algorithms.
*   **Compiling C Code:** Compiling the generated C code into a library or executable.

### Example Use Cases

Suppose we want to use the core C code output module to generate the C code for a simple numerical integration algorithm:


```python
import NRPy_param_funcs as par

# Set up input parameters
par.set_parval_from_inputargs('NRPy_step_size', 1e-6)
par.set_parval_from_inputargs('NRPy_max_steps', 1000)

# Generate the C code for the numerical integration algorithm
C_code = par.generate_C_code()

print(C_code)
```

This will output the generated C code.

## Step 2: Setting Up Parameters


```python
par.set_parval_from_inputargs('NRPy_step_size', 1e-6)
**NRPy+: Parameter Interface**
=============================

### Overview of the NRPy+ Parameter Interface

The NRPy+ parameter interface is a set of functions for setting up parameters in the NRPy+ framework.

### Theory Review

#### Introduction to Parameters in NRPy+

*   **Parameters:** Parameters are values that need to be specified by the user before running an NRPy+ simulation.
*   **NRPy+:** The Numerical Relativity in Python (NRPy+) framework is a software tool for numerical relativity.

### Code Implementation


```python
import sympy as sp
```

### Theory Review

#### What are Parameters in NRPy+?

Parameters in NRPy+ are values that need to be specified by the user before running an NRPy+ simulation. These parameters can include:

*   **Step Size:** The step size for numerical integration.
*   **Max Steps:** The maximum number of steps for numerical integration.

### Code Explanation


```python
import sympy as sp
```

This code imports the `sympy` library, which is a Python library for symbolic mathematics.

### Theory Review

#### How are Parameters Set Up in NRPy+?

Parameters in NRPy+ are set up using the following functions:

*   **`set_parval_from_inputargs()`:** This function sets the value of a parameter based on an input argument.
*   **`get_parval()`:** This function returns the current value of a parameter.

### Example Use Cases

Suppose we want to use the NRPy+ parameter interface to set up a step size for numerical integration:


```python
import sympy as sp
from NRPy_param_funcs import par

# Set up input parameters
par.set_parval_from_inputargs('NRPy_step_size', 1e-6)

# Get the current value of the step size
step_size = par.get_parval('NRPy_step_size')

print(step_size)
```

This will output the current value of the step size.

## Step 2: Setting Up Parameters with `set_parval_from_inputargs()`

```python
par.set_parval_from_inputargs('NRPy_max_steps', 1000)
```

This code sets up an input parameter for the maximum number of steps.

## Step 3: Getting the Current Value of a Parameter


```python
step_size = par.get_parval('NRPy_step_size')
```

This code gets the current value of the step**SymPy: A Computer Algebra Package Used by NRPy+**
=====================================================

### Overview of SymPy

SymPy is a Python library for symbolic mathematics. It is used by the NRPy+ framework to perform various calculations.

### Theory Review

#### Introduction to SymPy

*   **SymPy:** SymPy is a computer algebra package that can be used for symbolic mathematics.
*   **NRPy+:** The Numerical Relativity in Python (NRPy+) framework uses SymPy for its calculations.

### Code Implementation


```python
import indexedexp as ixp
```

### Theory Review

#### What is IndexedExp?

IndexedExp is a module that provides an interface to SymPy's indexed expressions. It is used by NRPy+ to perform various calculations.

*   **Indexed Expressions:** Indexed expressions are used to represent mathematical quantities with multiple indices.
*   **SymPy:** The `indexedexp` module uses SymPy to perform symbolic mathematics.

### Code Explanation


```python
import indexedexp as ixp
```

This code imports the `indexedexp` module, which provides an interface to SymPy's indexed expressions.

### Theory Review

#### How is IndexedExp Used by NRPy+?

The `indexedexp` module is used by NRPy+ to perform various calculations. It allows users to represent mathematical quantities with multiple indices and perform symbolic mathematics on them.

*   **Multiple Indices:** The `indexedexp` module provides a way to represent mathematical quantities with multiple indices.
*   **Symbolic Mathematics:** The module uses SymPy's capabilities for symbolic mathematics to perform various calculations.

### Example Use Cases

Suppose we want to use the `indexedexp` module to calculate the derivative of a function with respect to one of its variables:


```python
import indexedexp as ixp
from sympy import symbols, diff

# Define the function
f = x**2 + 3*x - 4

# Calculate the derivative of f with respect to x
derivative_f_x = diff(f, x)

print(derivative_f_x)
```

This will output the derivative of `f` with respect to `x`.

## Step 2: Creating an Indexed Expression


```python
import indexedexp as ixp

# Define the coordinates
x, y, z = symbols('x y z')

# Create an indexed expression for the coordinate x
xi = ixp.declare('xi', x)

print(xi**NRPy+: Support for Symbolic Indexed Expressions**
=====================================================

### Overview of Symbolic Indexed Expression Support in NRPy+

In this section, we will discuss the symbolic indexed expression (e.g., tensors, vectors, etc.) support provided by the NRPy+ framework.

### Theory Review

#### Introduction to Symbolic Indexed Expressions

*   **Symbolic Indexed Expressions:** These are mathematical expressions that can be manipulated symbolically using computer algebra systems like SymPy.
*   **Tensors and Vectors:** Tensors and vectors are types of symbolic indexed expressions used in numerical relativity.

### Code Implementation


```python
import nrpy
from nrpy import *
```

### Theory Review

#### What are Symbolic Indexed Expressions?

Symbolic indexed expressions are mathematical expressions that can be manipulated symbolically using computer algebra systems like SymPy. They are used to represent tensors and vectors in numerical relativity.

*   **Tensors:** Tensors are multi-dimensional arrays with indices that transform under changes of coordinate system.
*   **Vectors:** Vectors are one-dimensional arrays with an index that transforms under changes of coordinate system.

### Code Explanation


```python
import nrpy
from nrpy import *
```

This code imports the `nrpy` module and all its submodules, which provide a set of functions for numerical relativity. The `*` symbol is used to import all functions from the `nrpy` module.

### Theory Review

#### How are Symbolic Indexed Expressions Used in NRPy+?

The symbolic indexed expression support provided by NRPy+ allows users to represent tensors and vectors as mathematical expressions that can be manipulated using computer algebra systems like SymPy. This enables users to perform complex calculations involving tensors and vectors, such as tensor contractions and vector dot products.

*   **Tensor Contractions:** Tensor contractions are operations on tensors that involve the contraction of indices.
*   **Vector Dot Products:** Vector dot products are operations on vectors that involve the product of two or more vectors.

### Example Use Cases

Suppose we want to use the symbolic indexed expression support provided by NRPy+ to calculate the stress-energy tensor:


```python
import nrpy
from nrpy import *

# Define the coordinates
x, y, z = symbols('x y z')

# Create a symbolic indexed expression for the stress-energy tensor
T00 = x**2 + 3*x - 4

print(T00)
```

This**Step 2: Defining the Stress-Energy Tensor**
=============================================

### Overview of Defining the Stress-Energy Tensor

In this step, we will define the stress-energy tensor $T^{\mu\nu}$ and $T^\mu{}_\nu$.

### Theory Review

#### Introduction to the Stress-Energy Tensor

*   **Stress-Energy Tensor:** The stress-energy tensor is a mathematical object that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The stress-energy tensor is defined as:

$$
T^{\mu\nu} = \rho u^\mu u^\nu + p g^{\mu\nu}
$$

where $\rho$ is the density, $u^\mu$ is the four-velocity, and $p$ is the pressure.

### Code Implementation


```python
import sympy as sp

# Define the coordinates
x = sp.symbols('x')

# Define the stress-energy tensor T^{\mu\nu}
T00 = x**2 + 3*x - 4
T01 = x + 1
T10 = x + 1
T11 = x**2 - 3

print("T00 =", T00)
print("T01 =", T01)
print("T10 =", T10)
print("T11 =", T11)
```

### Theory Review

#### What are the Components of the Stress-Energy Tensor?

The stress-energy tensor has four components: $T^{\mu\nu}$. The components are defined as:

*   **$T^{00}$:** This component represents the energy density.
*   **$T^{01}$ and $T^{10}$:** These components represent the momentum fluxes.
*   **$T^{11}$:** This component represents the pressure.

### Code Explanation


```python
import sympy as sp

# Define the coordinates
x = sp.symbols('x')

# Define the stress-energy tensor T^{\mu\nu}
T00 = x**2 + 3*x - 4
T01 = x + 1
T10 = x + 1
T11 = x**2 - 3
```

This code defines the components of the stress-energy tensor using SymPy. The `x` variable is a symbol that represents the spatial coordinate.

### Theory Review

#### How are the Components of the Stress-Energy Tensor Related?

The components**Recall of the Stress-Energy Tensor**
=====================================

### Overview of the Stress-Energy Tensor

In this section, we will recall the definition of the stress-energy tensor $T^{\mu\nu}$ and its relation to the energy density $\rho_0$ and pressure $P$.

### Theory Review

#### Introduction to the Stress-Energy Tensor

*   **Stress-Energy Tensor:** The stress-energy tensor is a mathematical object that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The stress-energy tensor is defined as:

$$
T^{\mu\nu} = \rho_0 h u^{\mu} u^{\nu} + P g^{\mu\nu}
$$

where $\rho_0$ is the rest mass density, $h = 1 + \epsilon + \frac{P}{\rho_0}$, and $g^{\mu\nu}$ is the metric tensor.

### Code Implementation


```python
import sympy as sp

# Define the symbols
x, y, z, t = sp.symbols('x y z t')
rho_0 = sp.symbols('rho_0')
P = sp.symbols('P')
epsilon = sp.symbols('epsilon')

# Define the stress-energy tensor T^{\mu\nu}
T00 = rho_0 * (1 + epsilon + P/rho_0)
T01 = 0
T10 = 0
T11 = rho_0 * (1 + epsilon + P/rho_0)

print("T00 =", T00)
print("T01 =", T01)
print("T10 =", T10)
print("T11 =", T11)
```

### Theory Review

#### What are the Components of the Stress-Energy Tensor?

The stress-energy tensor has four components: $T^{\mu\nu}$. The components are defined as:

*   **$T^{00}$:** This component represents the energy density.
*   **$T^{01}$ and $T^{10}$:** These components represent the momentum fluxes.
*   **$T^{11}$:** This component represents the pressure.

### Code Explanation


```python
import sympy as sp

# Define the symbols
x, y, z, t = sp.symbols('x y z t')
rho_0 = sp.symbols('rho_0')
P = sp.s**Defining Enthalpy (h)**
=========================

### Overview of Defining Enthalpy (h)

In this step, we will define the enthalpy ($h$) as a function of rest mass density ($\rho_b$), pressure ($P$), and specific internal energy ($\epsilon$).

### Theory Review

#### Introduction to Enthalpy (h)

*   **Enthalpy:** The enthalpy is a thermodynamic property that represents the total energy per unit rest mass.
*   **Definition:** The enthalpy is defined as:

$$
h = 1 + \epsilon + \frac{P}{\rho_b}
$$

where $\rho_b$ is the rest mass density, $P$ is the pressure, and $\epsilon$ is the specific internal energy.

### Code Implementation


```python
def compute_enthalpy(rho_b,P,epsilon):
    global h
    h = 1 + epsilon + P/rho_b
```

This code defines a function `compute_enthalpy` that takes three arguments: `rho_b`, `P`, and `epsilon`. The function returns the enthalpy (`h`) as defined above.

### Theory Review

#### Why is Enthalpy Important?

The enthalpy is an important concept in thermodynamics, particularly in the context of numerical relativity. It allows us to describe the total energy per unit rest mass, which is essential for simulations involving relativistic fluids.

*   **Relativistic Fluid Dynamics:** The enthalpy plays a crucial role in relativistic fluid dynamics, where it is used to describe the behavior of fluids under extreme conditions.
*   **Numerical Simulations:** In numerical simulations, the enthalpy is often used as a conserved quantity to track the evolution of the system.

### Code Explanation


```python
def compute_enthalpy(rho_b,P,epsilon):
    global h
    h = 1 + epsilon + P/rho_b
```

This code defines a function `compute_enthalpy` that takes three arguments: `rho_b`, `P`, and `epsilon`. The function uses the `global` keyword to modify the value of `h` in the current scope.

### Example Use Case


```python
# Define the inputs
rho_b = 1.0  # Rest mass density
P = 0.5      # Pressure
epsilon = 0.2  #**Computing the Stress-Energy Tensor $T^{mu nu}$**
=====================================================

### Overview of Computing the Stress-Energy Tensor $T^{mu nu}$

In this step, we will define a function `compute_T4UU` that computes the stress-energy tensor $T^{mu nu}$ using the inputs provided.

### Theory Review

#### Introduction to the Stress-Energy Tensor $T^{mu nu}$

*   **Stress-Energy Tensor:** The stress-energy tensor $T^{mu nu}$ is a 4-dimensional tensor that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The stress-energy tensor is defined as:

$$
T^{\mu\nu} = \rho_0 h u^\mu u^\nu + P g^{\mu\nu}
$$

where $\rho_0$ is the rest mass density, $h$ is the enthalpy, $u^\mu$ is the four-velocity, and $g^{\mu\nu}$ is the metric tensor.

### Code Implementation


```python
def compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U):
    global T4UU

    compute_enthalpy(rho_b,P,epsilon)
    
    # Compute h
    h = 1 + epsilon + P/rho_b
    
    # Initialize the stress-energy tensor
    T4UU = gammaDD * u4U**2 + P * gammaDD
    
    return T4UU
```

This code defines a function `compute_T4UU` that takes seven inputs: `gammaDD`, `betaU`, `alpha`, `rho_b`, `P`, `epsilon`, and `u4U`. The function uses the `global` keyword to modify the value of `T4UU` in the current scope.

### Theory Review

#### How is the Stress-Energy Tensor $T^{mu nu}$ Related to Other Quantities?

The stress-energy tensor $T^{mu nu}$ is related to other quantities such as the rest mass density $\rho_0$, enthalpy $h$, four-velocity $u^\mu$, and metric tensor $g^{\mu\nu}$. The relationship between these quantities can be expressed using the following equations:

*   **Energy Density:** The energy density is given by:

$$
\rho = \rho_0 h u^0 u^0 + P g**Defining $g^{mu nu}$ in Terms of ADM Quantities**
=====================================================

### Overview of Defining $g^{mu nu}$ in Terms of ADM Quantities

In this step, we will define the metric tensor $g^{mu nu}$ in terms of the ADM quantities using the `AB4m` module.

### Theory Review

#### Introduction to ADM Quantities

*   **ADM Quantities:** The Arnowitt-Deser-Misner (ADM) quantities are a set of variables used to describe the spacetime geometry.
*   **Definition:** The ADM quantities include:
	+ $\gamma_{ij}$: The 3-metric
	+ $\beta^i$: The shift vector
	+ $\alpha$: The lapse function

### Code Implementation


```python
import BSSN.ADMBSSN_tofrom_4metric as AB4m

# Define the ADM quantities
ADM_quantities = AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
```

This code imports the `AB4m` module and uses its `g4UU_ito_BSSN_or_ADM` function to convert the 4-metric $g^{\mu nu}$ into ADM quantities.

### Theory Review

#### How is $g^{mu nu}$ Related to Other Quantities?

The metric tensor $g^{mu nu}$ is related to other quantities such as the 3-metric $\gamma_{ij}$, shift vector $\beta^i$, and lapse function $\alpha$. The relationship between these quantities can be expressed using the following equations:

*   **Metric Tensor:** The metric tensor is given by:

$$
g^{\mu\nu} = \left( \begin{array}{cccc}
\alpha^{-2} & -\alpha^{-2}\beta_i \\
-\alpha^{-2}\beta_j & \gamma_{ij}
\end{array} \right)
$$

where $\alpha$ is the lapse function, $\beta^i$ is the shift vector, and $\gamma_{ij}$ is the 3-metric.

### Code Explanation


```python
import BSSN.ADMBSSN_tofrom_4metric as AB4m

# Define the ADM quantities
ADM_quantities = AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
**Computing the Stress-Energy Tensor $T^{mu nu}$**
=====================================================

### Overview of Computing the Stress-Energy Tensor $T^{mu nu}$

In this final step, we will compute the stress-energy tensor $T^{mu nu}$ using the previously defined quantities.

### Theory Review

#### Introduction to the Stress-Energy Tensor $T^{mu nu}$

*   **Stress-Energy Tensor:** The stress-energy tensor $T^{mu nu}$ is a 4-dimensional tensor that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The stress-energy tensor is defined as:

$$
T^{\mu\nu} = \rho_0 h u^\mu u^\nu + P g^{\mu\nu}
$$

where $\rho_0$ is the rest mass density, $h$ is the enthalpy, $u^\mu$ is the four-velocity, and $g^{\mu\nu}$ is the metric tensor.

### Code Implementation


```python
# Initialize T4UU as a 4x4 matrix
T4UU = ixp.zerorank2(DIM=4)

# Compute each component of T4UU
for mu in range(4):
    for nu in range(4):
        # Component T4UU[mu][nu] is given by the sum of two terms
        T4UU[mu][nu] = rho_b * h * u4U[mu]*u4U[nu] + P*AB4m.g4UU[mu][nu]

# The computed stress-energy tensor T4UU is now available for use
```

This code initializes the $T^{mu nu}$ matrix as a 4x4 matrix using `ixp.zerorank2(DIM=4)`. It then computes each component of the matrix by iterating over the indices $mu$ and $nu$, using the previously defined quantities to compute the components.

### Theory Review

#### How is the Stress-Energy Tensor $T^{mu nu}$ Related to Other Quantities?

The stress-energy tensor $T^{mu nu}$ is related to other quantities such as the rest mass density $\rho_0$, enthalpy $h$, four-velocity $u^\mu$, and metric tensor $g^{\mu\nu}$. The relationship between these quantities can be expressed using the following equations:

*   **Energy Density:****Computing the Stress-Energy Tensor $T^{mu}_{nu}$**
=====================================================

### Overview of Computing the Stress-Energy Tensor $T^{mu}_{nu}$

In this step, we will define a function `compute_T4UD` to compute the stress-energy tensor $T^{mu}_{nu}$.

### Theory Review

#### Introduction to the Stress-Energy Tensor $T^{mu}_{nu}$

*   **Stress-Energy Tensor:** The stress-energy tensor is a mathematical object that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The stress-energy tensor is defined as:

$$
T^{\mu}_{\nu} = T^{\mu\delta} g_{\delta \nu}
$$

where $g_{\delta \nu}$ is the metric tensor.

### Code Implementation


```python
def compute_T4UD(gammaDD,betaU,alpha, T4UU):
    global T4UD
    
    # Compute T^{mu}_{nu} by contracting T^{mu nu} with g_{nu delta}
    T4UD = ixp.zerorank2(DIM=4)
    
    for mu in range(4):
        for nu in range(4):
            T4UD[mu][nu] = T4UU[mu][nu]
```

This code defines a function `compute_T4UD` that takes four arguments: `gammaDD`, `betaU`, `alpha`, and `T4UU`. The function uses the `global` keyword to modify the value of `T4UD` in the current scope.

### Theory Review

#### How is $T^{mu}_{nu}$ Related to Other Quantities?

The stress-energy tensor $T^{mu}_{nu}$ is related to other quantities such as the metric tensor $g_{\delta \nu}$. The relationship between these quantities can be expressed using the following equations:

*   **Contracted Stress-Energy Tensor:** The contracted stress-energy tensor is given by:

$$
T^{\mu}_{\nu} = T^{\mu\delta} g_{\delta \nu}
$$

where $g_{\delta \nu}$ is the metric tensor.

### Code Explanation


```python
def compute_T4UD(gammaDD,betaU,alpha, T4UU):
    global T4UD
    
    # Compute T^{mu}_{nu} by contracting T^{mu nu} with g**Computing the Contracted Stress-Energy Tensor $T^{\mu}_{\nu}$**
=================================================================

### Overview of Computing the Contracted Stress-Energy Tensor $T^{\mu}_{\nu}$

In this step, we will compute the contracted stress-energy tensor $T^{\mu}_{\nu}$ using the previously computed stress-energy tensor $T^{mu nu}$.

### Theory Review

#### Introduction to the Contracted Stress-Energy Tensor $T^{\mu}_{\nu}$

*   **Contracted Stress-Energy Tensor:** The contracted stress-energy tensor is a mathematical object that describes the distribution of energy and momentum in spacetime.
*   **Definition:** The contracted stress-energy tensor is defined as:

$$
T^{\mu}_{\nu} = T^{mu \delta} g_{\delta \nu}
$$

where $g_{\delta \nu}$ is the metric tensor.

### Code Implementation


```python
# Compute the contracted stress-energy tensor T^mu_nu
T_mu_nu = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        T_mu_nu[mu][nu] = T4UU[mu][nu]
```

This code computes the contracted stress-energy tensor $T^{\mu}_{\nu}$ by contracting the previously computed stress-energy tensor $T^{mu nu}$ with the metric tensor $g_{\delta \nu}$.

### Theory Review

#### How is the Contracted Stress-Energy Tensor $T^{\mu}_{\nu}$ Related to Other Quantities?

The contracted stress-energy tensor $T^{\mu}_{\nu}$ is related to other quantities such as the metric tensor $g_{\delta \nu}$. The relationship between these quantities can be expressed using the following equations:

*   **Contracted Stress-Energy Tensor:** The contracted stress-energy tensor is given by:

$$
T^{\mu}_{\nu} = T^{mu \delta} g_{\delta \nu}
$$

where $g_{\delta \nu}$ is the metric tensor.

### Code Explanation


```python
# Compute the contracted stress-energy tensor T^mu_nu
T_mu_nu = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        T_mu_nu[mu][nu] =**Converting to Conserved Quantities**
=====================================

### Overview of Converting to Conserved Quantities

In this section, we will convert the primitive variables to conserved quantities using the relationship between them.

### Theory Review

#### Introduction to Primitive and Conserved Quantities

*   **Primitive Variables:** The primitive variables are the input variables that describe the system's initial conditions.
*   **Conserved Quantities:** The conserved quantities are the output variables that describe the system's evolution over time.
*   **Relationship between Primitive and Conserved Quantities:** The relationship between primitive and conserved quantities is given by:

$$
T^{mu \nu} = T^{\prime mu \nu}
$$

where $T^{\prime mu \nu}$ are the conserved quantities.

### Code Implementation


```python
import BSSN.ADMBSSN_tofrom_4metric as AB4m

# Convert to ADM quantities
AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

# Initialize T^mu_nu
T4UD = ixp.zerorank2(DIM=4)
for mu in range(4):
    for nu in range(4):
        for delta in range(4):
            # Calculate the contraction of T^{\prime mu \nu} with g_{delta \nu}
            T4UD[mu][nu] += T4UU[mu][delta]*AB4m.g4DD[delta][nu]
```

This code converts the primitive variables to conserved quantities using the relationship between them.

### Theory Review

#### How is the Contraction of $T^{\prime mu \nu}$ with $g_{\delta \nu}$ Related to Other Quantities?

The contraction of $T^{\prime mu \nu}$ with $g_{\delta \nu}$ is related to other quantities such as the metric tensor $g_{\delta \nu}$. The relationship between these quantities can be expressed using the following equations:

*   **Contracted Stress-Energy Tensor:** The contracted stress-energy tensor is given by:

$$
T^{\mu}_{\nu} = T^{\prime mu \nu}
$$

where $T^{\prime mu \nu}$ are the conserved quantities.

### Code Explanation


```python
import BSSN.ADMBSSN_tofrom_4metric as**Writing Conservative Variables in Terms of Primitive Variables**
=============================================================

### Overview of Writing Conservative Variables in Terms of Primitive Variables

In this step, we will write the conservative variables in terms of the primitive variables.

### Theory Review

#### Introduction to Conservative and Primitive Variables

*   **Conservative Variables:** The conservative variables are a set of quantities that describe the system's evolution over time.
*   **Primitive Variables:** The primitive variables are the input variables that describe the system's initial conditions.
*   **Relationship between Conservative and Primitive Variables:** The relationship between conservative and primitive variables is given by:

$$
\left(\begin{array}{c} D \\ S_j \\ J \end{array}\right) = \left(\begin{array}{ccc} 1 & 0 & 0 \\ 0 & \alpha^{-2}\beta_i & 0 \\ 0 & 0 & 1/\sqrt{\gamma} \end{array}\right) \left(\begin{array}{c} \rho_0 \\ S^i \\ J^s \end{array}\right)
$$

where $\rho_0$ is the rest mass density, $S^i$ are the momentum variables, and $J^s$ is the matter flux.

### Code Implementation


```python
# Define the primitive variables
rho_0 = 1.0  # Rest mass density
S_i = 2.0    # Momentum variables
J_s = 3.0    # Matter flux

# Compute the conservative variables using the relationship between them
D = rho_0
S_j = alpha**-2 * beta_i * S_i
J = 1 / np.sqrt(gamma) * J_s
```

This code computes the conservative variables in terms of the primitive variables.

### Theory Review

#### How is the Relationship between Conservative and Primitive Variables Related to Other Quantities?

The relationship between conservative and primitive variables is related to other quantities such as the lapse function $\alpha$, shift vector $\beta_i$, and 3-metric $\gamma_{ij}$. The relationship between these quantities can be expressed using the following equations:

*   **Lapse Function:** The lapse function is given by:

$$
\alpha = \frac{\partial t}{\partial \tau}
$$

where $t$ is the time coordinate and $\tau$ is the proper time.
*   **Shift Vector:** The shift vector is given**Converting Primitive Variables to Conservative Variables**
===========================================================

### Overview of Converting Primitive Variables to Conservative Variables

In this section, we will convert the primitive variables to conservative variables using the relationships between them.

### Theory Review

#### Introduction to Primitive and Conservative Variables

*   **Primitive Variables:** The primitive variables are the input variables that describe the system's initial conditions.
*   **Conservative Variables:** The conservative variables are a set of quantities that describe the system's evolution over time.
*   **Relationship between Primitive and Conservative Variables:** The relationship between primitive and conservative variables is given by:

$$
\left(\begin{array}{c} \rho_* \\ \tilde{\tau} \\ \tilde{S}_i \end{array}\right) = \left(\begin{array}{ccc} \alpha\sqrt{\gamma} & 0 & 0 \\ -\rho_* & \alpha^2\sqrt{\gamma} & 0 \\ 0 & \alpha\sqrt{\gamma} & 0 \end{array}\right) \left(\begin{array}{c} \rho_0 u^0 \\ T^{00} \\ T^0{}_i \end{array}\right)
$$

where $\rho_0$ is the rest mass density, $u^0$ is the zeroth component of the 4-velocity, and $T^{00}$ and $T^0{}_i$ are components of the stress-energy tensor.

### Code Implementation


```python
import numpy as np

# Define the primitive variables
rho_0 = 1.0  # Rest mass density
u0 = 2.0     # Zeroth component of the 4-velocity
T00 = 3.0    # Component of the stress-energy tensor
Ti = 4.0     # Components of the stress-energy tensor

# Define the conservative variables using the relationship between them
rho_star = alpha * np.sqrt(gammaDET) * rho_0 * u0
tau_tilde = alpha**2 * np.sqrt(gammaDET) * T00 - rho_star
S_i_tilde = alpha * np.sqrt(gammaDET) * Ti
```

This code computes the conservative variables in terms of the primitive variables.

### Theory Review

#### How is the Relationship between Primitive and Conservative Variables Related to Other Quantities?

The relationship between primitive and conservative variables is related to other quantities**Computing Conservative Variables in Terms of Primitive Variables**
================================================================

### Overview of Computing Conservative Variables in Terms of Primitive Variables

In this step, we will compute the conservative variables in terms of the primitive variables.

### Theory Review

#### Introduction to Conservative and Primitive Variables

*   **Conservative Variables:** The conservative variables are a set of quantities that describe the system's evolution over time.
*   **Primitive Variables:** The primitive variables are the input variables that describe the system's initial conditions.
*   **Relationship between Conservative and Primitive Variables:** The relationship between conservative and primitive variables is given by:

$$
\left(\begin{array}{c} \rho_* \\ \tilde{\tau} \\ \tilde{S}_i \end{array}\right) = \left(\begin{array}{ccc} \alpha\sqrt{\gamma} & 0 & 0 \\ -\rho_* & \alpha^2\sqrt{\gamma} & 0 \\ 0 & \alpha\sqrt{\gamma} & 0 \end{array}\right) \left(\begin{array}{c} \rho_0 u^0 \\ T^{00} \\ T^0{}_i \end{array}\right)
$$

where $\rho_0$ is the rest mass density, $u^0$ is the zeroth component of the 4-velocity, and $T^{00}$ and $T^0{}_i$ are components of the stress-energy tensor.

### Code Implementation


```python
import sympy as sp
import ixp

# Define the function to compute sqrt(gammaDET)
def compute_sqrtgammaDET(gammaDD):
    global sqrtgammaDET
    
    # Compute gammaUU and gammaDET using ixp.symm_matrix_inverter3x3
    _gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    
    # Compute sqrt(gammaDET) using sympy.sqrt
    sqrtgammaDET = sp.sqrt(gammaDET)

# Define the function to compute rho_star
def compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U):
    global rho_star
    
    # Compute rho_* using the formula: rho_* = alpha * sqrt(gammaDET) * rho_0 * u^0
    rho_star = alpha * sqrtgammaDET * rho_b * u4U[0]
```

This code computes the conservative variables in**Computing Flux Variables**
==========================

### Overview of Computing Flux Variables

In this section, we will compute the flux variables in terms of the primitive variables.

### Theory Review

#### Introduction to Flux Variables

*   **Flux Variables:** The flux variables are a set of quantities that describe the system's evolution over time.
*   **Definition:** The flux variables are defined as:

$$
\left(\begin{array}{c} D \\ S_j \\ J \end{array}\right) = \left(\begin{array}{ccc} 1 & 0 & 0 \\ 0 & \alpha^{-2}\beta_i & 0 \\ 0 & 0 & 1/\sqrt{\gamma} \end{array}\right) \left(\begin{array}{c} \rho_* \\ S^i \\ J^s \end{array}\right)
$$

where $\rho_*$ is the conservative density, $S^i$ are the momentum variables, and $J^s$ is the matter flux.

### Code Implementation


```python
# Define the function to compute rho_star
def compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U):
    global rho_star
    
    # Compute rho_* using the formula: rho_* = alpha * sqrt(gammaDET) * rho_0 * u^0
    rho_star = alpha * sqrtgammaDET * rho_b * u4U[0]

# Define the function to compute tau_tilde
def compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star):
    global tau_tilde
    
    # Compute tau_tilde using the formula: tau_tilde = alpha^2 * sqrt(gammaDET) * T^{00} - rho_*
    tau_tilde = alpha**2*sqrtgammaDET*T4UU[0][0] - rho_star

# Define the function to compute S_tildeD
def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    global S_tildeD
    
    # Initialize S_tildeD as a 3-element array
    S_tildeD = ixp.zerorank1(DIM=3)
    
    # Compute each component of S_tildeD using the formula: S_j = alpha * sqrt(gammaDET) * T^0{}_j
    for i in range(3):
        S_tildeD**Computing Flux Variables for GRHD Equations**
=============================================

### Overview of Computing Flux Variables for GRHD Equations

In this step, we will define the flux variables for the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The GRHD equations are defined as:

$$
\partial_t \rho_* + \partial_i S_j = 0
$$

where $\rho_*$ is the conservative density, $S_j$ are the momentum variables, and $\partial_t$ and $\partial_i$ are the partial derivatives with respect to time and space.

### Code Implementation


```python
# Define the function to compute D (conservative density)
def compute_D(alpha, sqrtgammaDET, rho_b, u4U):
    global D
    
    # Compute D using the formula: D = \rho_*
    D = alpha * sqrtgammaDET * rho_b * u4U[0]

# Define the function to compute S_j (momentum variables)
def compute_S_j(alpha, sqrtgammaDET, betaU, T4UD):
    global S_j
    
    # Initialize S_j as a 3-element array
    S_j = ixp.zerorank1(DIM=3)
    
    # Compute each component of S_j using the formula: S_j = \alpha^2 * sqrt(gammaDET) * \beta_i * T^{0i}
    for i in range(3):
        S_j[i] = alpha**2 * sqrtgammaDET * betaU[i] * T4UD[0][i+1]

# Define the function to compute J (matter flux)
def compute_J(alpha, sqrtgammaDET, u4U, S_tildeD):
    global J
    
    # Compute J using the formula: J = 1 / \sqrt{\gamma} * S^s
    J = alpha**-2 * sqrtgammaDET * u4U[0] * S_tildeD[0]
```

This code defines the flux variables for the GRHD equations.

### Theory Review

#### How are the Flux Variables Related to Other Quantities?

The flux variables are related to other quantities such as the conservative density $\**Computing Flux Terms for GRHD Equations**
=========================================

### Overview of Computing Flux Terms for GRHD Equations

In this section, we will compute the flux terms for the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Flux Terms

*   **GRHD Flux Terms:** The GRHD flux terms are a set of quantities that describe the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The GRHD flux terms are defined as:

$$
\left(\begin{array}{c} D \\ S_j \\ J \end{array}\right) = \left(\begin{array}{ccc} 1 & 0 & 0 \\ 0 & \alpha^{-2}\beta_i & 0 \\ 0 & 0 & 1/\sqrt{\gamma} \end{array}\right) \left(\begin{array}{c} \rho_* \\ S^i \\ J^s \end{array}\right)
$$

where $\rho_*$ is the conservative density, $S^i$ are the momentum variables, and $J^s$ is the matter flux.

### Code Implementation


```python
# Define the function to compute the D-flux term
def compute_D_flux_term(alpha, sqrtgammaDET, rho_b, u4U):
    global D_flux_term
    
    # Compute the D-flux term using the formula: D-flux = \rho_*
    D_flux_term = alpha * sqrtgammaDET * rho_b * u4U[0]

# Define the function to compute the S_j-flux terms
def compute_S_j_flux_terms(alpha, betaU, T4UD):
    global S_j_flux_terms
    
    # Initialize S_j_flux_terms as a 3-element array
    S_j_flux_terms = ixp.zerorank1(DIM=3)
    
    # Compute each component of the S_j-flux terms using the formula: S_j-flux = \alpha^2 * sqrt(gammaDET) * \beta_i * T^{0i}
    for i in range(3):
        S_j_flux_terms[i] = alpha**2 * sqrtgammaDET * betaU[i] * T4UD[0][i+1]

# Define the function to compute the J-flux term
def compute_J_flux_term(alpha, u4U, S_t**Computing $\rho_*$ Flux Term**
================================

### Overview of Computing $\rho_*$ Flux Term

In this step, we will compute the $\rho_*$ flux term for the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to $\rho_*$ Flux Term

*   **$\rho_*$ Flux Term:** The $\rho_*$ flux term is a quantity that describes the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The $\rho_*$ flux term is defined as:

$$
\left(\begin{array}{c} D \\ S_j \\ J \end{array}\right) = \left(\begin{array}{ccc} 1 & 0 & 0 \\ 0 & \alpha^{-2}\beta_i & 0 \\ 0 & 0 & 1/\sqrt{\gamma} \end{array}\right) \left(\begin{array}{c} \rho_* \\ S^i \\ J^s \end{array}\right)
$$

where $\rho_*$ is the conservative density, $S^i$ are the momentum variables, and $J^s$ is the matter flux.

### Code Implementation


```python
# Define the function to compute the D-flux term
def compute_D_flux_term(alpha, sqrtgammaDET, rho_b, u4U):
    global D_flux_term
    
    # Compute the D-flux term using the formula: D-flux = \rho_*
    D_flux_term = alpha * sqrtgammaDET * rho_b * u4U[0]
```

This code computes the $\rho_*$ flux term.

### Theory Review

#### How is the $\rho_*$ Flux Term Related to Other Quantities?

The $\rho_*$ flux term is related to other quantities such as the conservative density $\rho_*$, the momentum variables $S^i$, and the matter flux $J^s$. The relationship between these quantities can be expressed using the following equations:

*   **Conservative Density:** The conservative density is given by:

$$
\rho_* = \alpha\sqrt{\gamma} \rho_0 u^0
$$

where $\rho_0$ is the rest mass density and $u^0$ is the zeroth component of the 4-velocity.
*   **Momentum Variables:** The momentum variables are given by:

$$**Computing Flux Term for Conservative Density**
=============================================

### Overview of Computing Flux Term for Conservative Density

In this section, we will compute the flux term for the conservative density.

### Theory Review

#### Introduction to Flux Term for Conservative Density

*   **Flux Term:** The flux term is a quantity that describes the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The flux term is defined as:

$$
\partial_t \rho_* + \partial_j (\rho_* v^j) = 0
$$

where $\rho_*$ is the conservative density and $v^j$ is the velocity.

### Code Implementation


```python
# Define the function to compute the velocity
def compute_velocity(u4U):
    global velocity
    
    # Compute the velocity using the formula: v^j = u^j/u^0
    for i in range(3):
        velocity[i] = u4U[i+1] / u4U[0]

# Define the function to compute the flux term
def compute_flux_term(alpha, sqrtgammaDET, rho_b, velocity):
    global flux_term
    
    # Compute the flux term using the formula: \rho_* v^j
    for i in range(3):
        flux_term[i] = alpha * sqrtgammaDET * rho_b * velocity[i]
```

This code computes the flux term for the conservative density.

### Theory Review

#### How is the Flux Term Related to Other Quantities?

The flux term is related to other quantities such as the conservative density $\rho_*$, the velocity $v^j$, and the metric tensor. The relationship between these quantities can be expressed using the following equations:

*   **Conservative Density:** The conservative density is given by:

$$
\rho_* = \alpha\sqrt{\gamma} \rho_0 u^0
$$

where $\rho_0$ is the rest mass density and $u^0$ is the zeroth component of the 4-velocity.
*   **Velocity:** The velocity is given by:

$$
v^j = u^j/u^0
$$**Computing Flux Terms for GRHD Equations**
==========================================

### Overview of Computing Flux Terms for GRHD Equations

In this section, we will define the flux terms for the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The GRHD equations are defined as:

$$
\partial_t \rho_* + \partial_j (\rho_* v^j) = 0
$$

where $\rho_*$ is the conservative density and $v^j$ is the velocity.

### Code Implementation


```python
import numpy as np

# Define the function to compute the flux term for conservative density
def compute_rho_star_flux_term(alpha, sqrtgammaDET, rho_b, u4U):
    global rho_star_flux_term
    
    # Compute the flux term using the formula: \rho_* v^j
    rho_star_flux_term = alpha * sqrtgammaDET * rho_b * (u4U[1] / u4U[0], u4U[2] / u4U[0], u4U[3] / u4U[0])

# Define the function to compute the flux term for momentum density
def compute_S_j_flux_term(alpha, betaU, T4UD):
    global S_j_flux_term
    
    # Compute the flux term using the formula: \alpha^2 \beta_i T^{0i}
    S_j_flux_term = alpha**2 * np.array([betaU[1] * T4UD[0][1], betaU[2] * T4UD[0][2], betaU[3] * T4UD[0][3]])

# Define the function to compute the flux term for energy density
def compute_J_flux_term(alpha, u4U, S_tildeD):
    global J_flux_term
    
    # Compute the flux term using the formula: 1 / \sqrt{\gamma} (\tilde{S}_s)
    J_flux_term = (alpha**-2 * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * u4U[0] * S_tildeD[0])
```

This**Computing $v^i$ Components from $u^\mu$**
==========================================

### Overview of Computing $v^i$ Components from $u^\mu$

In this section, we will compute the components of the 3-velocity $v^i$ from the 4-velocity $u^\mu$. This is necessary for computing the flux term $\rho_* v^j$.

### Theory Review

#### Introduction to 3-Velocity and 4-Velocity

*   **4-Velocity:** The 4-velocity is a vector that describes the velocity of an object in spacetime. It is given by $u^\mu = (\gamma, \gamma \beta^i)$ where $\gamma$ is the Lorentz factor and $\beta^i$ are the components of the 3-velocity.
*   **3-Velocity:** The 3-velocity is a vector that describes the velocity of an object in space. It is given by $v^i = \beta^i / \gamma$.

### Code Implementation


```python
# Define the function to compute vU from u4U
def compute_vU_from_u4U__no_speed_limit(u4U):
    global vU
    
    # Compute each component of vU using the formula: v^i = u^i / u^0
    for i in range(3):
        vU[i] = u4U[i+1] / u4U[0]
```

This code computes the components of the 3-velocity $v^i$ from the 4-velocity $u^\mu$.**Computing 3-Velocity Components**
==================================

### Overview of Computing 3-Velocity Components

In this section, we will compute the components of the 3-velocity $v^i$ from the 4-velocity $u^\mu$. This is necessary for computing the flux term $\rho_* v^j$.

### Theory Review

#### Introduction to 3-Velocity and 4-Velocity

*   **4-Velocity:** The 4-velocity is a vector that describes the velocity of an object in spacetime. It is given by $u^\mu = (\gamma, \gamma \beta^i)$ where $\gamma$ is the Lorentz factor and $\beta^i$ are the components of the 3-velocity.
*   **3-Velocity:** The 3-velocity is a vector that describes the velocity of an object in space. It is given by $v^i = \beta^i / \gamma$.

### Code Implementation


```python
# Define the function to compute vU from u4U
def compute_vU_from_u4U__no_speed_limit(u4U):
    global vU
    
    # Initialize vU as a 3-element array
    vU = ixp.zerorank1(DIM=3)
    
    # Compute each component of vU using the formula: v^i = u^i / u^0
    for j in range(3):
        vU[j] = u4U[j+1]/u4U[0]
```

This code computes the components of the 3-velocity $v^i$ from the 4-velocity $u^\mu$.

### Theory Review

#### Mathematical Derivation of 3-Velocity Components

The 3-velocity components can be computed using the following mathematical formula:

$$
v^i = \frac{u^i}{u^0}
$$

where $u^i$ are the spatial components of the 4-velocity and $u^0$ is the time component of the 4-velocity.

### Example Use Case

This code can be used to compute the 3-velocity components for a given 4-velocity. For example, if we have the 4-velocity $u^\mu = (2, 1, 1, 1)$, we can use this code to compute the 3-velocity components**Computing $\rho_* v^i$ Flux**
================================

### Overview of Computing $\rho_* v^i$ Flux

In this section, we will compute the flux term $\rho_* v^i$, which is a component of the GRHD equations.

### Theory Review

#### Introduction to $\rho_* v^i$ Flux

*   **$\rho_* v^i$ Flux:** The $\rho_* v^i$ flux is a quantity that represents the change in the conservative density over time due to the velocity of the fluid.
*   **Definition:** The $\rho_* v^i$ flux is given by:

$$
\rho_* v^i = \rho_* u^i
$$

where $\rho_*$ is the conservative density and $u^i$ are the components of the 3-velocity.

### Code Implementation


```python
# Define the function to compute rho_star_fluxU
def compute_rho_star_fluxU(vU, rho_star):
    global rho_star_fluxU
    
    # Initialize rho_star_fluxU as a 3-element array
    rho_star_fluxU = ixp.zerorank1(DIM=3)
    
    # Compute each component of rho_star_fluxU using the formula: rho_* v^i = rho_* u^i
    for j in range(3):
        rho_star_fluxU[j] = rho_star * vU[j]
```

This code computes the $\rho_* v^i$ flux.

### Theory Review

#### Mathematical Derivation of $\rho_* v^i$ Flux

The $\rho_* v^i$ flux can be computed using the following mathematical formula:

$$
\rho_* v^i = \rho_* u^i
$$

where $\rho_*$ is the conservative density and $u^i$ are the components of the 3-velocity.

### Example Use Case

This code can be used to compute the $\rho_* v^i$ flux for a given conservative density and velocity. For example, if we have the conservative density $\rho_* = 1$ and the velocity $v^i = (1, 2, 3)$, we can use this code to compute the $\rho_* v^i$ flux.

### Output

The output of this code is a 3-element array representing the $\rho_* v^i$ flux. The elements of the array are computed using the formula**Computing $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms**
==========================================================

### Overview of Computing $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms

In this section, we will compute the flux terms for the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The GRHD equations are defined as:

$$
\partial_t \rho_* + \partial_j (\rho_* v^j) = 0
$$

where $\rho_*$ is the conservative density and $v^j$ is the velocity.

### Code Implementation


```python
# Define the function to compute tau_tilde flux term
def compute_tau_tilde_flux_term(alpha, T4UU):
    global tau_tilde_flux_term
    
    # Compute the flux term using the formula: \alpha^2 \sqrt{\gamma} T^{00}
    tau_tilde_flux_term = alpha**2 * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * T4UU[0][0]

# Define the function to compute S_i flux terms
def compute_S_i_flux_terms(alpha, betaU, T4UD):
    global S_i_flux_terms
    
    # Initialize S_i_flux_terms as a 3-element array
    S_i_flux_terms = ixp.zerorank1(DIM=3)
    
    # Compute each component of the S_i flux terms using the formula: \alpha \sqrt{\gamma} T^0{}_i
    for i in range(3):
        S_i_flux_terms[i] = alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * betaU[i] * T4UD[0][i+1]
```

This code computes the flux terms for the GRHD equations.

### Theory Review

#### Mathematical Derivation of $\tilde{\tau}$ and $\tilde{S}_i$ Flux Terms

The $\tilde{\tau}$ and $\tilde{S}_i$ flux terms can be computed using the following mathematical formulas:

*  **Computing Terms for GRHD Equations**
=====================================

### Overview of Computing Terms for GRHD Equations

In this section, we will compute the terms that go inside the $\partial_j$'s on the left-hand side of the GRHD equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Definition:** The GRHD equations are defined as:

$$
\partial_t \tilde{\tau} + \partial_j (\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j) = s \\
\partial_t \tilde{S}_i + \partial_j (\alpha \sqrt{\gamma} T^j{}_i) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $\tilde{\tau}$ and $\tilde{S}_i$ are the internal energy density and momentum density, respectively.

### Code Implementation


```python
# Define the function to compute alpha squared gamma T^{0j}
def compute_alpha2_gamma_T_0j(alpha, sqrtgammaDET, T4UD):
    global alpha2_gamma_T_0j
    
    # Compute the term using the formula: \alpha^2 \sqrt{\gamma} T^{0j}
    for j in range(3):
        alpha2_gamma_T_0j[j] = alpha**2 * np.sqrt(sqrtgammaDET) * T4UD[0][j+1]

# Define the function to compute rho_* v^j
def compute_rho_star_v_j(rho_b, u4U):
    global rho_star_v_j
    
    # Compute the term using the formula: \rho_* v^j
    for j in range(3):
        rho_star_v_j[j] = rho_b * u4U[j+1]

# Define the function to compute alpha gamma T^{j}_i
def compute_alpha_gamma_T_j_i(alpha, betaU, T4UD):
    global alpha_gamma_T_j_i
    
    # Compute the term using the formula: \alpha \sqrt{\gamma} T^j{}_i
    for i in range(3):
        alpha_gamma_T_j_i[i]**Computing $\tilde{\tau}$ Flux**
================================

### Overview of Computing $\tilde{\tau}$ Flux

In this section, we will compute the flux term for the internal energy density, $\tilde{\tau}$. This is a component of the GRHD equations.

### Theory Review

#### Introduction to Internal Energy Density and Flux

*   **Internal Energy Density:** The internal energy density, $\tilde{\tau}$, represents the total energy content of the fluid per unit volume.
*   **Flux Term:** The flux term for $\tilde{\tau}$ is given by:

$$
\partial_t \tilde{\tau} + \partial_j (\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j) = s
$$

where $s$ is the source term.

### Code Implementation


```python
# Define the function to compute tau_tilde_fluxU
def compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU,T4UU, rho_star):
    global tau_tilde_fluxU
    
    # Initialize tau_tilde_fluxU as a 3-element array
    tau_tilde_fluxU = ixp.zerorank1(DIM=3)
    
    # Compute each component of the tau_tilde_fluxU using the formula: 
    # alpha^2 \sqrt{\gamma} T^{0j} - rho_* v^j
    for j in range(3):
        tau_tilde_fluxU[j] = alpha**2*sqrtgammaDET*T4UU[0][j+1] - rho_star*vU[j]
```

This code computes the flux term for the internal energy density.

### Theory Review

#### Mathematical Derivation of $\tilde{\tau}$ Flux Term

The $\tilde{\tau}$ flux term can be computed using the following mathematical formula:

$$
\partial_j (\alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j) = \alpha^2 \sqrt{\gamma} T^{0j} - \rho_* v^j
$$

where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric, $T^{0j}$ is the 0-j component of the stress-energy tensor, and $v^j$ is the j-component of the velocity.

### Example Use**Computing $\tilde{S}_i$ Flux**
==============================

### Overview of Computing $\tilde{S}_i$ Flux

In this section, we will compute the flux term for the momentum density, $\tilde{S}_i$. This is a component of the GRHD equations.

### Theory Review

#### Introduction to Momentum Density and Flux

*   **Momentum Density:** The momentum density, $\tilde{S}_i$, represents the total momentum content of the fluid per unit volume.
*   **Flux Term:** The flux term for $\tilde{S}_i$ is given by:

$$
\partial_t \tilde{S}_i + \partial_j (\alpha \sqrt{\gamma} T^j{}_i) = \frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $g_{\mu\nu}$ is the metric tensor.

### Code Implementation


```python
# Define the function to compute S_tilde_fluxUD
def compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD):
    global S_tilde_fluxUD
    
    # Initialize S_tilde_fluxUD as a 3x3 matrix
    S_tilde_fluxUD = ixp.zerorank2(DIM=3)
    
    # Compute each component of the S_tilde_fluxUD using the formula: 
    # \alpha \sqrt{\gamma} T^j{}_i
    for j in range(3):
        for i in range(3):
            S_tilde_fluxUD[j][i] = alpha*sqrtgammaDET*T4UD[j+1][i+1]
```

This code computes the flux term for the momentum density.

### Theory Review

#### Mathematical Derivation of $\tilde{S}_i$ Flux Term

The $\tilde{S}_i$ flux term can be computed using the following mathematical formula:

$$
\partial_j (\alpha \sqrt{\gamma} T^j{}_i) = \alpha \sqrt{\gamma} T^j{}_i
$$

where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric, $T^j{}_i$ is the j-i component of the stress-energy tensor.

### Example Use Case

This code can be used to**Computing Source Terms**
=========================

### Overview of Computing Source Terms

In this section, we will compute the source terms that appear on the right-hand sides (RHS) of the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations and Source Terms

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Source Terms:** The source terms appear on the RHS of the GRHD equations and represent the effects of external forces or sources on the fluid.

### Code Implementation


```python
# Define the function to compute s (source term for internal energy density)
def compute_s(alpha, betaU, S_tildeD, T4UD):
    global s
    
    # Compute the source term using the formula: s = alpha \sqrt{\gamma} (u^i S_{,i})
    s = alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * u4U[1] * S_tildeD[1]

# Define the function to compute source term for momentum density
def compute_source_term_S_i(alpha, betaU, T4UD):
    global source_term_S_i
    
    # Compute the source term using the formula: 
    # \frac{1}{2} alpha\sqrt{\gamma} (T^{\mu\nu} g_{\mu\nu,i})
    for i in range(3):
        source_term_S_i[i] = 0.5 * alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * T4UD[0][i+1]
```

This code computes the source terms that appear on the RHS of the GRHD equations.

### Theory Review

#### Mathematical Derivation of Source Terms

The source terms can be computed using the following mathematical formulas:

*   **Source Term for Internal Energy Density:** The source term for internal energy density is given by:

$$
s = \alpha \sqrt{\gamma} (u^i S_{,i})
$$

where $u^i$ is the i-component of the 4-velocity and $S_{,i}$ is the i-component of the momentum density.

*   **Source Term**Source Terms for GRHD Equations**
==================================

### Overview of Source Terms for GRHD Equations

In this section, we will discuss the source terms that appear on the right-hand sides (RHS) of the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations and Source Terms

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the evolution of relativistic fluids in curved spacetime.
*   **Source Terms:** The source terms appear on the RHS of the GRHD equations and represent the effects of external forces or sources on the fluid.

### Code Implementation


```python
# Define the function to compute s (source term for internal energy density)
def compute_s(alpha, betaU, S_tildeD, T4UD):
    global s
    
    # Compute the source term using the formula: 
    # alpha \sqrt{\gamma} (u^i S_{,i})
    s = alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * u4U[1] * S_tildeD[1]

# Define the function to compute source term for momentum density
def compute_source_term_S_i(alpha, betaU, T4UD):
    global source_term_S_i
    
    # Compute the source term using the formula: 
    # 0.5 alpha\sqrt{\gamma} (T^{\mu\nu} g_{\mu\nu,i})
    for i in range(3):
        source_term_S_i[i] = 0.5 * alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * T4UD[0][i+1]
```

This code computes the source terms that appear on the RHS of the GRHD equations.

### Theory Review

#### Mathematical Derivation of Source Terms

The source terms can be computed using the following mathematical formulas:

*   **Source Term for Internal Energy Density:** The source term for internal energy density is given by:

$$
s = \alpha \sqrt{\gamma} (u^i S_{,i})
$$

where $u^i$ is the i-component of the 4-velocity and $S_{,i}$ is the i-component of the momentum density.

**Computing Source Term for Internal Energy Density**
=====================================================

### Overview of Computing Source Term for Internal Energy Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to $\tilde{\tau}$ Equation and Source Term

*   **$\tilde{\tau}$ Equation:** The $\tilde{\tau}$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.
*   **Source Term:** The source term appears on the RHS of the $\tilde{\tau}$ equation and represents the effects of external forces or sources on the internal energy density.

### Code Implementation


```python
# Define the function to compute s (source term for internal energy density)
def compute_s(alpha, betaU, S_tildeD, T4UD):
    global s
    
    # Compute the source term using the formula: 
    # alpha \sqrt{\gamma} (u^i S_{,i})
    s = alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * u4U[1] * S_tildeD[1]
```

This code computes the source term that appears on the RHS of the $\tilde{\tau}$ equation.

### Theory Review

#### Mathematical Derivation of Source Term for Internal Energy Density

The source term can be computed using the following mathematical formula:

$$
s = \alpha \sqrt{\gamma} (u^i S_{,i})
$$

where $u^i$ is the i-component of the 4-velocity and $S_{,i}$ is the i-component of the momentum density.

### Example Use Case

This code can be used to compute the source term for internal energy density in a variety of astrophysical simulations.**Computing Source Term for Internal Energy Density**
=====================================================

### Overview of Computing Source Term for Internal Energy Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to $\tilde{\tau}$ Equation and Source Term

*   **$\tilde{\tau}$ Equation:** The $\tilde{\tau}$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.
*   **Source Term:** The source term appears on the RHS of the $\tilde{\tau}$ equation and represents the effects of external forces or sources on the internal energy density.

### Mathematical Derivation of Source Term

The source term can be computed using the following mathematical formula:

$$
s = \alpha \sqrt{\gamma} (u^i S_{,i})
$$

where $u^i$ is the i-component of the 4-velocity and $S_{,i}$ is the i-component of the momentum density.

### Code Implementation


```python
# Define the function to compute s_source_term
def compute_s_source_term(KDD, betaU, alpha, sqrtgammaDET, alpha_dD, T4UU):
    global s_source_term
    
    # Compute the source term using the formula: 
    # \alpha \sqrt{\gamma} (u^i S_{,i})
    s_source_term = sp.sympify(0)
    
    # Define Term 1
    term_1 = KDD[0][1]*betaU[1]*betaU[2] + 2*KDD[0][2]*betaU[2] + KDD[2][2]
    
    # Define Term 2
    term_2 = - (T4UU[0][1]*betaU[1] + T4UU[0][2])
    
    # Compute the source term
    s_source_term += alpha*sqrtgammaDET*(term_1*alpha_dD[2][2] + term_2*sp.diff(alpha,1))
```

This code computes the source term that appears on the RHS of the $\tilde{\tau}$ equation.

### Example Use Case

This code can be used to compute the source term for internal energy density in a variety of astrophysical simulations.**Computing Term 1 of Source Term**
=====================================

### Overview of Computing Term 1 of Source Term

In this section, we will compute the first term of the source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to Term 1

*   **Term 1:** The first term of the source term is given by:

$$
\underbrace{T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \left(K_{ij}\right)}_{\text{Term 1}}
$$

where $T^{00}$, $T^{0i}$, and $T^{ij}$ are the components of the stress-energy tensor.

### Code Implementation


```python
# Compute Term 1
for i in range(3):
    for j in range(3):
        s_source_term += (T4UU[0][0]*betaU[i]*betaU[j] + 2*T4UU[0][i+1]*betaU[j] + T4UU[i+1][j+1])*KDD[i][j]
```

This code computes the first term of the source term.

### Mathematical Derivation of Term 1

The first term can be computed using the following mathematical formula:

$$
T^{00}\beta^i\beta^j + 2 T^{0i}\beta^j + T^{ij} \left(K_{ij}\right)
$$

where $T^{00}$, $T^{0i}$, and $T^{ij}$ are the components of the stress-energy tensor.

### Example Use Case

This code can be used to compute the first term of the source term for internal energy density in a variety of astrophysical simulations.**Computing Term 2 of Source Term**
=====================================

### Overview of Computing Term 2 of Source Term

In this section, we will compute the second term of the source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to Term 2

*   **Term 2:** The second term of the source term is given by:

$$
\underbrace{-\left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha}_{\text{Term 2}}
$$

where $T^{00}$ and $T^{0i}$ are the components of the stress-energy tensor.

### Code Implementation


```python
# Compute Term 2
for i in range(3):
    s_source_term += -(T4UU[0][0]*betaU[i] + T4UU[0][i+1])*alpha_dD[i]
```

This code computes the second term of the source term.

### Mathematical Derivation of Term 2

The second term can be computed using the following mathematical formula:

$$
-\left(T^{00}\beta^i + T^{0i} \right)\partial_i\alpha
$$

where $T^{00}$ and $T^{0i}$ are the components of the stress-energy tensor.

### Example Use Case

This code can be used to compute the second term of the source term for internal energy density in a variety of astrophysical simulations.**Computing Term 3 of Source Term**
=====================================

### Overview of Computing Term 3 of Source Term

In this section, we will compute the third term of the source term that appears on the right-hand side (RHS) of the $\tilde{\tau}$ equation.

### Theory Review

#### Introduction to Term 3

*   **Term 3:** The third term of the source term is given by:

$$
\underbrace{\alpha \sqrt{\gamma}}_{\text{Term 3}}
$$

where $\alpha$ is the lapse function and $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric.

### Code Implementation


```python
# Compute Term 3
s_source_term *= alpha*sqrtgammaDET
```

This code computes the third term of the source term by multiplying the previous result by $\alpha \sqrt{\gamma}$.

### Mathematical Derivation of Term 3

The third term can be computed using the following mathematical formula:

$$
\alpha \sqrt{\gamma}
$$

where $\alpha$ is the lapse function and $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric.

### Example Use Case

This code can be used to compute the final value of the source term for internal energy density in a variety of astrophysical simulations.

### Theory Review

#### Properties of Term 3

*   **Multiplicative Factor:** The third term is multiplied by $\alpha \sqrt{\gamma}$, which represents the effects of time dilation and spatial curvature on the internal energy density.
*   **Importance in GRHD Equations:** The source term plays a crucial role in the GRHD equations, as it determines the evolution of the internal energy density over time.

### Code Implementation


```python
# Define the function to compute the source term
def compute_s_source_term(KDD, betaU, alpha, sqrtgammaDET, alpha_dD, T4UU):
    global s_source_term
    
    # Compute Term 1
    for i in range(3):
        for j in range(3):
            s_source_term += (T4UU[0][0]*betaU[i]*betaU[j] + 2*T4UU[0][i+1]*betaU[j] + T4UU[i+1][j+1])*KDD[i][j]
    
    # Compute Term 2
    for i in range**Computing Source Term for Momentum Density**
=============================================

### Overview of Computing Source Term for Momentum Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to $\tilde{S}_i$ Equation and Source Term

*   **$\tilde{S}_i$ Equation:** The $\tilde{S}_i$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.
*   **Source Term:** The source term appears on the RHS of the $\tilde{S}_i$ equation and represents the effects of external forces or sources on the momentum density.

### Code Implementation


```python
# Define the function to compute the source term for momentum density
def compute_source_term_S_i(alpha, betaU, KDD, T4UD):
    global source_term_S_i
    
    # Compute the source term using the formula: 
    # 0.5 alpha\sqrt{\gamma} (T^{\mu\nu} g_{\mu\nu,i})
    for i in range(3):
        source_term_S_i[i] = 0.5 * alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * T4UD[0][i+1]
```

This code computes the source term that appears on the RHS of the $\tilde{S}_i$ equation.

### Mathematical Derivation of Source Term for Momentum Density

The source term can be computed using the following mathematical formula:

$$
\frac{1}{2} \alpha \sqrt{\gamma} (T^{\mu\nu} g_{\mu\nu,i})
$$

where $T^{\mu\nu}$ is the stress-energy tensor and $g_{\mu\nu,i}$ is the i-component of the metric tensor.

### Example Use Case

This code can be used to compute the source term for momentum density in a variety of astrophysical simulations.**Computing Source Term for Momentum Density**
=============================================

### Overview of Computing Source Term for Momentum Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to $\tilde{S}_i$ Equation and Source Term

*   **$\tilde{S}_i$ Equation:** The $\tilde{S}_i$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.
*   **Source Term:** The source term appears on the RHS of the $\tilde{S}_i$ equation and represents the effects of external forces or sources on the momentum density.

### Mathematical Derivation of Source Term for Momentum Density

The source term can be computed using the following mathematical formula:

$$
\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $T^{\mu\nu}$ is the stress-energy tensor and $g_{\mu\nu,i}$ is the i-component of the metric tensor.

### Code Implementation


```python
# Define the function to compute the source term for momentum density
def compute_source_term_S_i(alpha, betaU, KDD, T4UD):
    global source_term_S_i
    
    # Compute the source term using the formula: 
    # 0.5 alpha\sqrt{\gamma} (T^{\mu\nu} g_{\mu\nu,i})
    for i in range(3):
        source_term_S_i[i] = 0.5 * alpha * np.sqrt(1.0/np.array([1.0, 1.0, 1.0])) * T4UD[0][i+1]
```

This code computes the source term that appears on the RHS of the $\tilde{S}_i$ equation.

### Mathematical Derivation of Metric Derivatives

The metric derivatives $g_{\mu\nu,i}$ can be computed using the following mathematical formula:

$$
g_{\mu\nu,i} = \frac{\partial g_{\mu\nu}}{\partial x^i}
$$

where $x^i$ is the i-component of the coordinate vector.

### Example Use Case

This code can be used to compute the source term**Computing Metric Derivatives**
=============================

### Overview of Computing Metric Derivatives

In this section, we will compute the metric derivatives $g_{\mu\nu,i}$ in terms of the ADM quantities and their derivatives.

### Theory Review

#### Introduction to Metric Derivatives

*   **Metric Derivatives:** The metric derivatives $g_{\mu\nu,i}$ represent the partial derivatives of the metric tensor with respect to the coordinate vector.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of Metric Derivatives

The metric derivatives can be computed using the following mathematical formula:

$$
g_{\mu\nu,i} = \frac{\partial g_{\mu\nu}}{\partial x^i}
$$

where $x^i$ is the i-component of the coordinate vector.

### Code Implementation


```python
# Define the function to compute metric derivatives
def compute_metric_derivatives(gDD, xD):
    global metric_derivatives
    
    # Compute the metric derivatives using the formula: 
    # g_{\mu\nu,i} = \frac{\partial g_{\mu\nu}}{\partial x^i}
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                metric_derivatives[mu][nu][i] = sp.diff(gDD[mu][nu], xD[i])
```

This code computes the metric derivatives in terms of the ADM quantities and their derivatives.

### Example Use Case

This code can be used to compute the metric derivatives for a variety of spacetime geometries.**Computing Metric Derivatives**
=============================

### Overview of Computing Metric Derivatives

In this section, we will compute the metric derivatives $g_{\mu\nu,k}$ in terms of the ADM quantities and their derivatives.

### Theory Review

#### Introduction to Metric Derivatives

*   **Metric Derivatives:** The metric derivatives $g_{\mu\nu,k}$ represent the partial derivatives of the 4-metric with respect to the coordinate vector.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of Metric Derivatives

The metric derivatives can be computed using the following mathematical formula:

$$
g_{\mu\nu,k} = \begin{pmatrix}
-2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k} & \beta_{i,k} \\
\beta_{j,k} & \gamma_{ij,k}
\end{pmatrix},
$$

where $\beta_i = \gamma_{ij} \beta^j$.

### Code Implementation


```python
# Define the function to compute metric derivatives
def compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha, gammaDD_dD, betaU_dD, alpha_dD):
    global g4DD_zerotimederiv_dD
    
    # Compute the metric derivatives using the formula: 
    # g_{\mu\nu,k} = \begin{pmatrix}
    # -2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k} & \beta_{i,k} \\
    # \beta_{j,k} & \gamma_{ij,k}
    # \end{pmatrix}
    for k in range(3):
        g4DD_zerotimederiv_dD[0][0][k] = -2 * alpha * alpha_dD[k]
        g4DD_zerotimederiv_dD[0][1][k] = betaU_dD[k] + betaU[1]*betaU_dD[0]
        g4DD_zerotimederiv_dD[1][0][k] = betaU**Computing Beta Derivatives**
=============================

### Overview of Computing Beta Derivatives

In this section, we will compute the beta derivatives using Equation 2.121 from B&S.

### Theory Review

#### Introduction to Beta Derivatives

*   **Beta Derivatives:** The beta derivatives represent the partial derivatives of the shift vector with respect to the coordinate vector.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of Beta Derivatives

The beta derivatives can be computed using the following mathematical formula:

$$
\beta_i = \gamma_{ij} \beta^j
$$

where $\gamma_{ij}$ is the spatial metric and $\beta^j$ is the j-component of the shift vector.

### Code Implementation


```python
# Define the function to compute beta derivatives
def compute_beta_derivatives(gammaDD, betaU):
    global betaD
    
    # Compute the beta derivatives using the formula: 
    # beta_i = gamma_{ij} beta^j
    for i in range(3):
        betaD[i] = 0.0
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]
```

This code computes the beta derivatives using Equation 2.121 from B&S.

### Mathematical Derivation of Beta Derivative Second Order Terms

The second-order terms of the beta derivative can be computed using the following mathematical formula:

$$
\partial_k \beta_i = \gamma_{ij,k} \beta^j + \gamma_{ij} \beta^j_{,k}
$$

where $\gamma_{ij,k}$ is the k-component of the spatial metric derivative.

### Code Implementation


```python
# Define the function to compute beta derivative second-order terms
def compute_beta_derivative_second_order_terms(gammaDD_dD, betaU_dD):
    global betaDdD
    
    # Compute the beta derivative second-order terms using the formula: 
    # partial_k beta_i = gamma_{ij,k} beta^j + gamma_{ij} beta^j_{,k}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                betaDdD[i][j][k] = 0.0
               **Computing Beta Derivative Second-Order Terms**
=============================================

### Overview of Computing Beta Derivative Second-Order Terms

In this section, we will compute the beta derivative second-order terms using Equation 2.121 from B&S.

### Theory Review

#### Introduction to Beta Derivative Second-Order Terms

*   **Beta Derivative Second-Order Terms:** The beta derivative second-order terms represent the partial derivatives of the shift vector with respect to the coordinate vector, up to the second order.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of Beta Derivative Second-Order Terms

The beta derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k \beta_i = \gamma_{ij,k} \beta^j + \gamma_{ij} \beta^j_{,k}
$$

where $\gamma_{ij,k}$ is the k-component of the spatial metric derivative.

### Code Implementation


```python
# Define the function to compute beta derivative second-order terms
def compute_beta_derivative_second_order_terms(gammaDD_dD, gammaDD, betaU_dD, betaU):
    global betaDdD
    
    # Compute the beta derivative second-order terms using the formula: 
    # partial_k beta_i = gamma_{ij,k} beta^j + gamma_{ij} beta^j_{,k}
    for i in range(3):
        for k in range(3):
            betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j]
            betaDdD[i][k] += gammaDD[i][j]*betaU_dD[j][k]
```

This code computes the beta derivative second-order terms using Equation 2.121 from B&S.

### Mathematical Derivation of Beta Derivative Second-Order Terms

The mathematical derivation of the beta derivative second-order terms is based on the following:

*   **Beta Derivative:** The beta derivative is given by $\beta_i = \gamma_{ij} \beta^j$.
*   **Spatial Metric Derivative:** The spatial metric derivative is given by $\gamma_{ij,k}$.

### Example Use Case

This code can be used to compute the beta derivative second-order terms for a variety of spacetime geometries.**Computing 4-Metric Derivative Second-Order Terms**
==================================================

### Overview of Computing 4-Metric Derivative Second-Order Terms

In this section, we will compute the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Theory Review

#### Introduction to 4-Metric Derivative Second-Order Terms

*   **4-Metric Derivative Second-Order Terms:** The 4-metric derivative second-order terms represent the partial derivatives of the 4-metric with respect to the coordinate vector, up to the second order.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of 4-Metric Derivative Second-Order Terms

The 4-metric derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k g_{\mu\nu} = \begin{pmatrix}
-2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k} & \beta_{i,k} \\
\beta_{j,k} & \gamma_{ij,k}
\end{pmatrix},
$$

where $\alpha$ is the lapse function, $\beta^j$ is the j-component of the shift vector, and $\gamma_{ij}$ is the spatial metric.

### Code Implementation


```python
# Define the function to compute 4-metric derivative second-order terms
def compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha):
    global g4DD_zerotimederiv_dD
    
    # Compute the 4-metric derivative second-order terms using the formula: 
    # partial_k g_{\mu\nu} = \begin{pmatrix}
    # -2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k} & \beta_{i,k} \\
    # \beta_{j,k} & \gamma_{ij,k}
    # \end{pmatrix},
    g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)
    
    for k in range(3):
        for mu in range(4):
           **Computing 4-Metric Derivative Second-Order Terms**
==================================================

### Overview of Computing 4-Metric Derivative Second-Order Terms

In this section, we will compute the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Theory Review

#### Introduction to 4-Metric Derivative Second-Order Terms

*   **4-Metric Derivative Second-Order Terms:** The 4-metric derivative second-order terms represent the partial derivatives of the 4-metric with respect to the coordinate vector, up to the second order.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of 4-Metric Derivative Second-Order Terms

The 4-metric derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k g_{00} = -2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k}
$$

where $\alpha$ is the lapse function and $\beta^j$ is the j-component of the shift vector.

### Code Implementation


```python
# Define the function to compute 4-metric derivative second-order terms
def compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha):
    global g4DD_zerotimederiv_dD
    
    # Compute the 4-metric derivative second-order terms using the formula: 
    # partial_k g_{00} = -2 \alpha\alpha_{,k} + \beta^j_{,k} \beta_j + \beta^j \beta_{j,k}
    g4DD_zerotimederiv_dD[0][0] = -alpha**2 + betaU[j]*betaD[j]
    
    for k in range(3):
        g4DD_zerotimederiv_dD[0][0][k+1] += -2*alpha*alpha_dD[k]
        
        # Compute the second-order terms
        for j in range(3):
            g4DD_zerotimederiv_dD[0][0][k+1] += betaU_dD[j][k]*betaD[j] + betaU[j]*betaDdD[j][**Computing 4-Metric Derivative Second-Order Terms**
==================================================

### Overview of Computing 4-Metric Derivative Second-Order Terms

In this section, we will compute the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Theory Review

#### Introduction to 4-Metric Derivative Second-Order Terms

*   **4-Metric Derivative Second-Order Terms:** The 4-metric derivative second-order terms represent the partial derivatives of the 4-metric with respect to the coordinate vector, up to the second order.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of 4-Metric Derivative Second-Order Terms

The 4-metric derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k g_{0i} = \beta_{i,k}
$$

where $\beta_i$ is the i-component of the shift vector.

### Code Implementation


```python
# Define the function to compute 4-metric derivative second-order terms
def compute_g4DD_zerotimederiv_dD(gammaDD, betaU):
    global g4DD_zerotimederiv_dD
    
    # Compute the 4-metric derivative second-order terms using the formula: 
    # partial_k g_{0i} = beta_{i,k}
    for i in range(3):
        g4DD[i][0] = g4DD[0][i] = betaD[i]
        
        # Compute the second-order terms
        for k in range(3):
            g4DD_zerotimederiv_dD[i+1][0][k+1] = g4DD_zerotimederiv_dD[0][i+1][k+1] = betaDdD[i][k]

    # Compute the remaining second-order terms
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ```
This code computes the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Mathematical Derivation of Spatial Metric Derivative Second-Order Terms

The spatial metric derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k \gamma_{ij**Computing 4-Metric Derivative Second-Order Terms**
==================================================

### Overview of Computing 4-Metric Derivative Second-Order Terms

In this section, we will compute the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Theory Review

#### Introduction to 4-Metric Derivative Second-Order Terms

*   **4-Metric Derivative Second-Order Terms:** The 4-metric derivative second-order terms represent the partial derivatives of the 4-metric with respect to the coordinate vector, up to the second order.
*   **ADM Quantities:** The ADM quantities are used to describe the spacetime geometry in terms of the lapse function, shift vector, and metric tensor.

### Mathematical Derivation of 4-Metric Derivative Second-Order Terms

The 4-metric derivative second-order terms can be computed using the following mathematical formula:

$$
\partial_k g_{ij} = \gamma_{ij,k}
$$

where $\gamma_{ij}$ is the spatial metric and $\gamma_{ij,k}$ is the k-component of the spatial metric derivative.

### Code Implementation


```python
# Define the function to compute 4-metric derivative second-order terms
def compute_g4DD_zerotimederiv_dD(gammaDD, gammaDD_dD):
    global g4DD_zerotimederiv_dD
    
    # Compute the 4-metric derivative second-order terms using the formula: 
    # partial_k g_{ij} = gamma_{ij,k}
    for i in range(3):
        for j in range(3):
            g4DD[i][j] = gammaDD[i][j]
            
            # Compute the second-order terms
            for k in range(3):
                g4DD_zerotimederiv_dD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]
```

This code computes the 4-metric derivative second-order terms using Equation 2.122 from B&S.

### Example Use Case

This code can be used to compute the 4-metric derivative second-order terms for a variety of spacetime geometries.**Computing Source Term for Momentum Density**
=============================================

### Overview of Computing Source Term for Momentum Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to Source Term

*   **Source Term:** The source term represents the effects of external forces or sources on the momentum density.
*   **$\tilde{S}_i$ Equation:** The $\tilde{S}_i$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.

### Mathematical Derivation of Source Term

The source term can be computed using the following mathematical formula:

$$
\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric, $T^{\mu\nu}$ is the stress-energy tensor, and $g_{\mu\nu,i}$ is the i-component of the 4-metric derivative.

### Code Implementation


```python
# Define the function to compute source term for momentum density
def compute_source_term_S_i(alpha, betaU, KDD, T4UD):
    global g4DD
    
    # Compute the source term using the formula: 
    # source term = 0.5 * alpha * sqrt(gamma) * T^mu nu * g_{mu nu,i}
    for i in range(3):
        source_term = 0.5 * alpha * math.sqrt(KDD[i][i]) * T4UD[0][i]
        
        # Add the source term to the momentum density
        S_i += source_term
```

This code computes the source term that appears on the RHS of the $\tilde{S}_i$ equation.

### Example Use Case

This code can be used to compute the source term for a variety of spacetime geometries.**Computing Source Term for Momentum Density**
=============================================

### Overview of Computing Source Term for Momentum Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to Source Term

*   **Source Term:** The source term represents the effects of external forces or sources on the momentum density.
*   **$\tilde{S}_i$ Equation:** The $\tilde{S}_i$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.

### Mathematical Derivation of Source Term

The source term can be computed using the following mathematical formula:

$$
\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric, $T^{\mu\nu}$ is the stress-energy tensor, and $g_{\mu\nu,i}$ is the i-component of the 4-metric derivative.

### Code Implementation


```python
# Define the function to compute source term for momentum density
def compute_source_term_S_i(alpha, betaU, KDD, T4UD):
    global g4DD
    
    # Compute the source term using the formula: 
    # source term = 0.5 * alpha * sqrt(gamma) * T^mu nu * g_{mu nu,i}
    for i in range(3):
        source_term = 0.5 * alpha * math.sqrt(KDD[i][i]) * T4UD[0][i]
        
        # Add the source term to the momentum density
        S_i += source_term
    
    # Return the source term
    return source_term
```

This code computes the source term that appears on the RHS of the $\tilde{S}_i$ equation.

### Example Use Case

This code can be used to compute the source term for a variety of spacetime geometries.**Computing Source Term for Momentum Density**
=============================================

### Overview of Computing Source Term for Momentum Density

In this section, we will compute the source term that appears on the right-hand side (RHS) of the $\tilde{S}_i$ equation.

### Theory Review

#### Introduction to Source Term

*   **Source Term:** The source term represents the effects of external forces or sources on the momentum density.
*   **$\tilde{S}_i$ Equation:** The $\tilde{S}_i$ equation is a component of the General Relativistic Hydrodynamics (GRHD) equations.

### Mathematical Derivation of Source Term

The source term can be computed using the following mathematical formula:

$$
\frac{1}{2} \alpha\sqrt{\gamma} T^{\mu\nu} g_{\mu\nu,i}
$$

where $\alpha$ is the lapse function, $\sqrt{\gamma}$ is the square root of the determinant of the spatial metric, $T^{\mu\nu}$ is the stress-energy tensor, and $g_{\mu\nu,i}$ is the i-component of the 4-metric derivative.

### Code Implementation


```python
# Define the function to compute source term for momentum density
def compute_S_tilde_source_termD(alpha, sqrtgammaDET, g4DD_zerotimederiv_dD, T4UU):
    global S_tilde_source_termD
    
    # Compute the source term using the formula: 
    # source term = 0.5 * alpha * sqrt(gamma) * T^mu nu * g_{mu nu,i}
    S_tilde_source_termD = ixp.zerorank1(DIM=3)
    
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                S_tilde_source_termD[i] += sp.Rational(1,2)*alpha*sqrtgammaDET*T4UU[mu][nu]*g4DD_zerotimederiv_dD[mu][nu][i+1]
```

This code computes the source term that appears on the RHS of the $\tilde{S}_i$ equation.

### Example Use Case

This code can be used to compute the source term for a variety of spacetime geometries.**Conversion of Velocity Components**
=====================================

### Overview of Conversion Process

In this section, we will convert the velocity components from the spatial basis ($v^i$) to the spacetime basis ($u^\mu$).

### Theory Review

#### Introduction to Velocity Components

*   **Velocity Components:** The velocity components represent the rate of change of the position with respect to time.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the spatial basis to the spacetime basis is given by:

$$
u^\mu = \frac{\partial x^\mu}{\partial t} + v^i \gamma_{i}^\mu
$$

where $x^\mu$ is the position in spacetime, $\gamma_{i}^\mu$ is the spatial metric in the spacetime basis, and $v^i$ is the velocity component in the spatial basis.

### Code Implementation


```python
# Define the function to convert velocity components from spatial basis to spacetime basis
def convert_v_to_u(v, gammaDD):
    global u
    
    # Compute the conversion using the formula: 
    # u^\mu = partial x^\mu / partial t + v^i * gamma_{i}^\mu
    for mu in range(4):
        for i in range(3):
            u[mu] += v[i]*gammaDD[i][mu]
```

This code computes the velocity components in the spacetime basis using the conversion formula.

### Example Use Case

This code can be used to convert the velocity components from the spatial basis to the spacetime basis for a variety of spacetime geometries.**Computing 4-Velocity Components**
=====================================

### Overview of Computing 4-Velocity Components

In this section, we will compute the 4-velocity components $u^\mu$ from the given velocity components in the spatial basis $v^i$.

### Theory Review

#### Introduction to 4-Velocity Components

*   **4-Velocity Components:** The 4-velocity components represent the rate of change of the position with respect to time in spacetime.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the spatial basis to the spacetime basis is given by:

$$
u^\mu = \frac{\partial x^\mu}{\partial t} + v^i \gamma_{i}^\mu
$$

where $x^\mu$ is the position in spacetime, $\gamma_{i}^\mu$ is the spatial metric in the spacetime basis, and $v^i$ is the velocity component in the spatial basis.

### Code Implementation


```python
# Define the function to compute 4-velocity components from spatial basis to spacetime basis
def compute_u_from_v(v, betaU, gammaDD):
    global u
    
    # Compute the conversion using the formula: 
    # u^\mu = partial x^\mu / partial t + v^i * gamma_{i}^\mu
    for mu in range(4):
        for i in range(3):
            u[mu] += v[i]*gammaDD[i][mu]
```

This code computes the 4-velocity components $u^\mu$ from the given velocity components in the spatial basis.

### Algorithm for Computing $u^0$

The algorithm for computing $u^0$ involves several steps:

1. Choose a maximum Lorentz factor $\Gamma_{\rm max}$.
2. Compute $R=\gamma_{ij}v^i_{(n)}v^j_{(n)}=1 - \frac{1}{\Gamma^2}$.
3. If $R \le 1 - \frac{1}{\Gamma_{\rm max}^2}$, then skip the next step.
4. Otherwise if $R > 1 - \frac{1}{\Gamma_{\rm max}^2}$ then adjust $**Computing 4-Velocity Components**
=====================================

### Overview of Computing 4-Velocity Components

In this section, we will compute the 4-velocity components $u^\mu$ from the given Valencia 3-velocity $v_{(n)}^i$.

### Theory Review

#### Introduction to 4-Velocity Components

*   **4-Velocity Components:** The 4-velocity components represent the rate of change of the position with respect to time in spacetime.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the Valencia 3-velocity $v_{(n)}^i$ to the 4-velocity $u^\mu$ is given by:

$$
\alpha v_{(n)}^i = \frac{u^i}{u^0} + \beta^i
$$

where $\alpha$ is the lapse function, $\beta^i$ is the shift vector, and $v_{(n)}^i$ is the Valencia 3-velocity.

### Code Implementation


```python
# Define the function to compute 4-velocity components from Valencia 3-velocity
def convert_valencia_to_u(v_n, alpha, betaU):
    global u
    
    # Compute the conversion using the formula: 
    # u^i = alpha u^0 * (v_{(n)}^i - beta^i)
    for mu in range(4):
        u[mu] = 0
    for i in range(3):
        u[i+1] += alpha*(v_n[i] - betaU[i])
```

This code computes the 4-velocity components $u^\mu$ from the given Valencia 3-velocity $v_{(n)}^i$.

### Speed Limiter

To prevent numerical errors, a speed limiter is applied to ensure that the Lorentz factor $\Gamma = \frac{1}{\sqrt{1-R}}$ does not become too large. The speed limiter sets an upper limit on the Lorentz factor, preventing it from exceeding a certain value.

### Code Implementation


```python
# Define the function to apply speed limiter
def apply_speed_limiter(R, gamma_max):
    global R_star
    
    # Compute the speed limiter**Computing 4-Velocity Components**
=====================================

### Overview of Computing 4-Velocity Components

In this section, we will compute the 4-velocity components $u^\mu$ from the given Valencia 3-velocity $v_{(n)}^i$.

### Theory Review

#### Introduction to 4-Velocity Components

*   **4-Velocity Components:** The 4-velocity components represent the rate of change of the position with respect to time in spacetime.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the Valencia 3-velocity $v_{(n)}^i$ to the 4-velocity $u^\mu$ is given by:

$$
\alpha v_{(n)}^i = \frac{u^i}{u^0} + \beta^i
$$

where $\alpha$ is the lapse function, $\beta^i$ is the shift vector, and $v_{(n)}^i$ is the Valencia 3-velocity.

### Code Implementation


```python
# Define the function to compute 4-velocity components from Valencia 3-velocity
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU):
    global rescaledValenciavU
    
    # Compute the conversion using the formula: 
    # u^i = alpha u^0 * (v_{(n)}^i - beta^i)
    
    # Apply speed limiter
    R_star = 1 - (1 / GAMMA_SPEED_LIMIT**2)
    rescaledValenciavU = ValenciavU.copy()
    
    for i in range(3):
        if R_star < gammaDD[i][i]:
            rescaledValenciavU[i] *= math.sqrt(R_star / gammaDD[i][i])
    
    # Compute the 4-velocity components
    u = [0]*4
    for mu in range(4):
        for i in range(3):
            u[mu] += alpha * (rescaledValenciavU[i] - betaU[i]) * gammaDD[i][mu]
```

This code computes the 4-velocity components $u**Computing 4-Velocity Components**
=====================================

### Overview of Computing 4-Velocity Components

In this section, we will compute the 4-velocity components $u^\mu$ from the given Valencia 3-velocity $v_{(n)}^i$.

### Theory Review

#### Introduction to 4-Velocity Components

*   **4-Velocity Components:** The 4-velocity components represent the rate of change of the position with respect to time in spacetime.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the Valencia 3-velocity $v_{(n)}^i$ to the 4-velocity $u^\mu$ is given by:

$$
\alpha v_{(n)}^i = \frac{u^i}{u^0} + \beta^i
$$

where $\alpha$ is the lapse function, $\beta^i$ is the shift vector, and $v_{(n)}^i$ is the Valencia 3-velocity.

### Code Implementation


```python
# Define the function to compute 4-velocity components from Valencia 3-velocity
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU):
    global rescaledValenciavU
    
    # Compute the conversion using the formula: 
    # u^i = alpha u^0 * (v_{(n)}^i - beta^i)
    
    # Apply speed limiter
    R_star = 1 - (1 / GAMMA_SPEED_LIMIT**2)
    rescaledValenciavU = ValenciavU.copy()
    
    for i in range(3):
        if R_star < gammaDD[i][i]:
            rescaledValenciavU[i] *= math.sqrt(R_star / gammaDD[i][i])
    
    # Compute the 4-velocity components
    u = [0]*4
    for mu in range(4):
        for i in range(3):
            u[mu] += alpha * (rescaledValenciavU[i] - betaU[i]) * gammaDD[i][mu]
```

This code computes the 4-velocity components $u**Computing 4-Velocity Components**
=====================================

### Overview of Computing 4-Velocity Components

In this section, we will compute the 4-velocity components $u^\mu$ from the given Valencia 3-velocity $v_{(n)}^i$.

### Theory Review

#### Introduction to 4-Velocity Components

*   **4-Velocity Components:** The 4-velocity components represent the rate of change of the position with respect to time in spacetime.
*   **Spatial Basis vs. Spacetime Basis:** The velocity components can be represented in either the spatial basis or the spacetime basis.

### Mathematical Derivation of Conversion Formula

The conversion formula from the Valencia 3-velocity $v_{(n)}^i$ to the 4-velocity $u^\mu$ is given by:

$$
\alpha v_{(n)}^i = \frac{u^i}{u^0} + \beta^i
$$

where $\alpha$ is the lapse function, $\beta^i$ is the shift vector, and $v_{(n)}^i$ is the Valencia 3-velocity.

### Code Implementation


```python
# Define the function to compute 4-velocity components from Valencia 3-velocity
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU):
    global u4U_ito_ValenciavU, rescaledValenciavU
    
    # Compute the conversion using the formula: 
    # u^i = alpha u^0 * (v_{(n)}^i - beta^i)
    
    # Apply speed limiter
    R_star = 1 - (1 / GAMMA_SPEED_LIMIT**2)
    rescaledValenciavU = ValenciavU.copy()
    
    for i in range(3):
        if R_star < gammaDD[i][i]:
            rescaledValenciavU[i] *= math.sqrt(R_star / gammaDD[i][i])
    
    # Compute the 4-velocity components
    u = [0]*4
    for mu in range(4):
        for i in range(3):
            u[mu] += alpha * (rescaledValenciavU[i] - betaU[i]) * gammaDD[i][mu]
    
    #**Computing the Lorentz Factor**
==============================

### Overview of Computing the Lorentz Factor

In this section, we will compute the Lorentz factor $R$ from the given Valencia 3-velocity components $v^i$.

### Theory Review

#### Introduction to the Lorentz Factor

*   **Lorentz Factor:** The Lorentz factor is a dimensionless quantity that measures the amount of time dilation caused by special relativity.
*   **Computing the Lorentz Factor:** The Lorentz factor can be computed using the formula: $$ R = \gamma_{ij} v^i v^j $$
	+   where $\gamma_{ij}$ is the 3-metric, and $v^i$ is the Valencia 3-velocity component.

### Code Implementation


```python
# Define the function to compute the Lorentz factor
def compute_lorentz_factor(gammaDD, ValenciavU):
    global R
    
    # Compute the Lorentz factor using the formula: 
    # R = gamma_{ij} v^i v^j
    R = sp.sympify(0)
    
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]
```

This code computes the Lorentz factor $R$ from the given Valencia 3-velocity components.

### Theory Review

#### Introduction to SymPy

*   **SymPy:** SymPy is a Python library for symbolic mathematics.
	+   It can be used to perform computations involving mathematical expressions and equations.

### Code Implementation


```python
# Import the SymPy library
import sympy as sp
    
# Define the Lorentz factor using SymPy
R = sp.sympify(0)
```

This code uses SymPy to define the Lorentz factor $R$ as a symbolic expression.**Defining a Global Variable**
=============================

### Overview of Defining a Global Variable

In this section, we will define a global variable `GAMMA_SPEED_LIMIT` using the `par.Cparameters()` function from the `param` module.

### Theory Review

#### Introduction to Global Variables

*   **Global Variables:** Global variables are variables that can be accessed from any part of a program.
	+   They can be used to store values that need to be shared between different parts of a program.

#### Using `par.Cparameters()` Function

*   **`par.Cparameters()` Function:** The `par.Cparameters()` function is used to define global parameters in the C code.
	+   It takes four arguments: the type of parameter, the module name, the default value, and the actual value from the input file.

### Code Implementation


```python
# Import the param module
import par
    
# Define a global variable using par.Cparameters()
GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)
```

This code defines a global variable `GAMMA_SPEED_LIMIT` with a default value of 10.0. The actual value will be overwritten in the main C code.

### Theory Review

#### Importance of Default Value

*   **Default Value:** The default value is not terribly important here, since we can overwrite it in the main C code.
	+   However, having a default value allows us to provide a fallback value if the actual value is not specified.

#### Type and Module Name

*   **Type:** The type of parameter is set to "REAL", indicating that it is a real-valued parameter.
*   **Module Name:** The module name is set to `thismodule`, which is a global variable defined earlier.

### Code Implementation


```python
# Define the module name
thismodule = "GRHD"
```

This code defines the module name as "GRHD", which will be used to specify the module where the parameter is defined.**Defining a Global Variable**
=============================

### Overview of Defining a Global Variable

In this section, we will define a global variable `GAMMA_SPEED_LIMIT` using the `par.Cparameters()` function from the `param` module.

### Theory Review

#### Introduction to Global Variables

*   **Global Variables:** Global variables are variables that can be accessed from any part of a program.
	+   They can be used to store values that need to be shared between different parts of a program.

#### Using `par.Cparameters()` Function

*   **`par.Cparameters()` Function:** The `par.Cparameters()` function is used to define global parameters in the C code.
	+   It takes four arguments: the type of parameter, the module name, the default value, and the actual value from the input file.

### Code Implementation


```python
# Import the param module
import par
    
# Define a global variable using par.Cparameters()
GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)
```

This code defines a global variable `GAMMA_SPEED_LIMIT` with a default value of 10.0.

### Theory Review

#### Importance of Default Value

*   **Default Value:** The default value is used as a fallback if the actual value is not specified in the input file.
	+   It provides a starting point for the simulation, and can be adjusted later based on the specific requirements of the problem.

#### Default Value Based On

*   **`par.Cparameters()` Function:** The `par.Cparameters()` function allows us to specify a default value based on other parameters or values in the code.
	+   This provides flexibility and makes it easier to adjust the simulation settings without having to modify the code itself.

### Code Implementation


```python
# Define the module name
thismodule = "GRHD"
```

This code defines the module name as "GRHD", which will be used to specify the module where the parameter is defined.**Importing Modules**
=====================

### Overview of Importing Modules

In this section, we will import the necessary modules for the IllinoisGRMHD simulation.

### Theory Review

#### Introduction to Modules

*   **Modules:** Modules are pre-written code libraries that can be imported into a program to perform specific tasks.
	+   They provide a way to reuse code and make it easier to write complex programs by breaking them down into smaller, more manageable pieces.

#### Importing the IllinoisGRMHD Module

*   **`import` Statement:** The `import` statement is used to import modules into a program.
	+   It allows us to use functions, variables, and other code from the module within our own code.

### Code Implementation


```python
# Import the param module
import par
    
# Define the module name
thismodule = "GRHD"
```

This code imports the `param` module, which is used to define global parameters in the C code. The module name is also defined as "GRHD", which will be used to specify the module where the parameter is defined.

### Theory Review

#### Importance of Module Importation

*   **Module Importation:** Importing modules provides a way to reuse code and make it easier to write complex programs.
	+   It allows us to focus on writing our own code, rather than having to implement every function or feature from scratch.

#### IllinoisGRMHD Module

*   **IllinoisGRMHD Module:** The IllinoisGRMHD module is a pre-written code library that provides functions and variables for performing the IllinoisGRMHD simulation.
	+   It contains the necessary code for setting up the simulation, including defining global parameters, initializing variables, and performing calculations.**Defining Global Variables**
=============================

### Overview of Defining Global Variables

In this section, we will define two global variables: `Rmax` and `GiRaFFE`.

### Theory Review

#### Introduction to Global Variables

*   **Global Variables:** Global variables are variables that can be accessed from any part of a program.
	+   They can be used to store values that need to be shared between different parts of a program.

### Code Implementation


```python
# Define global variable GiRaFFE with default value 2000.0
GiRaFFE = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 2000.0)
```

This code defines the `GiRaFFE` global variable with a default value of 2000.0.

### Theory Review

#### Importance of Default Values

*   **Default Values:** Default values are used as fallbacks if actual values are not specified in the input file.
	+   They provide a starting point for simulations and can be adjusted later based on specific requirements.

### Code Implementation


```python
# Define global variable Rmax using Lorentz factor calculation
Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
```

This code defines the `Rmax` global variable using a calculation involving the Lorentz factor (`GAMMA_SPEED_LIMIT`) and default value for `GiRaFFE`.

### Theory Review

#### Lorentz Factor Calculation

*   **Lorentz Factor:** The Lorentz factor is a dimensionless quantity that measures time dilation in special relativity.
	+   In this case, it's used to calculate the maximum allowed velocity (`Rmax`) based on the `GAMMA_SPEED_LIMIT`.

### Mathematics


$$ R_{\text{max}} = 1 - \frac{1}{(\gamma \cdot c)^2} $$
```python
# Define global variable GAMMA_SPEED_LIMIT using par.Cparameters()
GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)
```

This code defines the `GAMMA_SPEED_LIMIT` global variable with a default value of 10.0.**Calculating the Lorentz Factor**
=====================================

### Overview of Calculating the Lorentz Factor

In this section, we will calculate the Lorentz factor (`R`), and then use it to determine `Rstar`.

### Theory Review

#### Introduction to the Lorentz Factor

*   **Lorentz Factor:** The Lorentz factor is a dimensionless quantity that measures time dilation in special relativity.
	+   It's calculated using the formula: $$ R = \gamma_{ij} v^i v^j $$
	+   Where $\gamma_{ij}$ is the 3-metric, and $v^i$ is the Valencia 3-velocity component.

### Code Implementation


```python
# Calculate the Lorentz factor (R) using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
R = sp.sympify(0)
for i in range(3):
    for j in range(3):
        R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

# Calculate the minimum of Rmax and R
Rstar = min(Rmax, R)
```

This code calculates the Lorentz factor (`R`) using the formula: $$ R = \gamma_{ij} v^i v^j $$
Then it uses the `min()` function to determine `Rstar`, which is the minimum of `Rmax` and `R`.

### Theory Review

#### Importance of Minimizing the Lorentz Factor

*   **Minimizing the Lorentz Factor:** The Lorentz factor (`R`) is used to calculate `Rstar`.
	+   By minimizing `R`, we ensure that the simulation remains physically valid.

### Code Implementation


```python
# Calculate Rmax using the formula: Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
```

This code calculates `Rmax` using the formula: $$ R_{\text{max}} = 1 - \frac{1}{(\gamma \cdot c)^2} $$
Where $\gamma$ is the Lorentz factor, and $c$ is the speed of light.

### Mathematics


$$ R_{\text{**Calculating the Lorentz Factor**
=====================================

### Overview of Calculating the Lorentz Factor

In this section, we will calculate the Lorentz factor (`R`) and determine `Rstar`.

### Theory Review

#### Introduction to the Lorentz Factor

*   **Lorentz Factor:** The Lorentz factor is a dimensionless quantity that measures time dilation in special relativity.
	+   It's calculated using the formula: $$ R = \gamma_{ij} v^i v^j $$
	+   Where $\gamma_{ij}$ is the 3-metric, and $v^i$ is the Valencia 3-velocity component.

### Code Implementation


```python
# Calculate the Lorentz factor (R) using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
R = sp.sympify(0)
for i in range(3):
    for j in range(3):
        R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

# Calculate the minimum of Rmax and R
if R < Rmax:
    # If R is less than Rmax, then Rstar = 0.5*(Rmax+R-Rmax+R) = R
    Rstar = R
```

This code calculates the Lorentz factor (`R`) using the formula: $$ R = \gamma_{ij} v^i v^j $$
Then it checks if `R` is less than `Rmax`. If it is, then `Rstar` is set to `R`.

### Theory Review

#### Importance of Minimizing the Lorentz Factor

*   **Minimizing the Lorentz Factor:** The Lorentz factor (`R`) is used to calculate `Rstar`.
	+   By minimizing `R`, we ensure that the simulation remains physically valid.

### Code Implementation


```python
# Calculate Rmax using the formula: Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
```

This code calculates `Rmax` using the formula: $$ R_{\text{max}} = 1 - \frac{1}{(\gamma \cdot c)^2} $**Calculating the Lorentz Factor**
=====================================

### Overview of Calculating the Lorentz Factor

In this section, we will calculate the Lorentz factor (`R`) and determine `Rstar`.

### Theory Review

#### Introduction to the Lorentz Factor

*   **Lorentz Factor:** The Lorentz factor is a dimensionless quantity that measures time dilation in special relativity.
	+   It's calculated using the formula: $$ R = \gamma_{ij} v^i v^j $$
	+   Where $\gamma_{ij}$ is the 3-metric, and $v^i$ is the Valencia 3-velocity component.

### Code Implementation


```python
# Calculate the Lorentz factor (R) using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
R = sp.sympify(0)
for i in range(3):
    for j in range(3):
        R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

# Calculate the minimum of Rmax and R
if R >= Rmax:
    # If R is greater than or equal to Rmax, then Rstar = 0.5*(Rmax+R-Rmax+R) = Rmax
    Rstar = Rmax
```

This code calculates the Lorentz factor (`R`) using the formula: $$ R = \gamma_{ij} v^i v^j $$
Then it checks if `R` is greater than or equal to `Rmax`. If it is, then `Rstar` is set to `Rmax`.

### Theory Review

#### Importance of Minimizing the Lorentz Factor

*   **Minimizing the Lorentz Factor:** The Lorentz factor (`R`) is used to calculate `Rstar`.
	+   By minimizing `R`, we ensure that the simulation remains physically valid.

### Code Implementation


```python
# Calculate Rmax using the formula: Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
Rmax = 1 - 1 / (GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)
```

This code calculates `Rmax` using the formula: $$ R_{\text{max}} = 1 - \frac{1}{**Avoiding Division by Zero**
=============================

### Overview of Avoiding Division by Zero

In this section, we will discuss how to avoid division by zero in the calculation of `Rstar`.

### Theory Review

#### Introduction to Division by Zero

*   **Division by Zero:** When dividing a number by zero, an error occurs.
	+   This can happen when calculating `Rstar` if `R` is equal to `Rmax`.
*   **Why it happens:** The division by zero occurs because we are trying to calculate the difference between two identical numbers (`Rmax - R`).
	+   When these numbers are identical, their difference is zero.

### Code Implementation


```python
# Calculate Rstar using the formula: Rstar = sp.Rational(1, 2) * (Rmax + R - nrpyAbs(Rmax - R))
Rstar = sp.Rational(1, 2) * (Rmax + R - nrpyAbs(Rmax - R))

# Add TINYDOUBLE to R below to avoid a 0/0
R += TINYDOUBLE
```

This code calculates `Rstar` using the formula: $$ \text{Rstar} = \frac{1}{2} (\text{Rmax} + \text{R} - |\text{Rmax} - \text{R}|) $$
Then it adds a small value (`TINYDOUBLE`) to `R` to avoid division by zero.

### Theory Review

#### Importance of Avoiding Division by Zero

*   **Avoiding Division by Zero:** By adding a small value to `R`, we ensure that the simulation remains physically valid.
	+   This is because division by zero can cause errors in the calculation of `Rstar`.

### Code Implementation


```python
# Define TINYDOUBLE as a small positive number
TINYDOUBLE = 1e-100

# Add TINYDOUBLE to R below to avoid a 0/0
R += TINYDOUBLE
```

This code defines `TINYDOUBLE` as a small positive number (1e-100). This value is added to `R` to avoid division by zero.

### Mathematics


$$ \text{Rstar} = \frac{1}{2} (\text{Rmax} + \text{R} - |\text{Rmax} - \text{R}|) $$**Valencia 3-Velocity Components**
=====================================

### Overview of Valencia 3-Velocity Components

In this section, we will discuss the Valencia 3-velocity components.

### Theory Review

#### Introduction to Valencia 3-Velocity Components

*   **Valencia 3-Velocity Components:** The Valencia 3-velocity components are a set of three velocity components that describe the motion of an object in three-dimensional space.
	+   They are used to calculate the Lorentz factor, which is a dimensionless quantity that measures time dilation in special relativity.

### Code Implementation


```python
# Define ValenciavU as a list of three Valencia 3-velocity components
ValenciavU = [0.0, 0.0, 0.0]
```

This code defines `ValenciavU` as a list of three Valencia 3-velocity components, all of which are initialized to zero.

### Theory Review

#### Importance of Valencia 3-Velocity Components

*   **Importance:** The Valencia 3-velocity components are used to calculate the Lorentz factor, which is a fundamental quantity in special relativity.
	+   They must be accurately calculated in order to obtain accurate results in simulations and calculations.

### Mathematics


$$ \text{ValenciavU} = (\text{v}_x, \text{v}_y, \text{v}_z) $$

*   **Valencia 3-Velocity Components:** The Valencia 3-velocity components are defined as:
	+   $$ \text{v}_x = \frac{\text{dx}}{\text{dt}} $$
	+   $$ \text{v}_y = \frac{\text{dy}}{\text{dt}} $$
	+   $$ \text{v}_z = \frac{\text{dz}}{\text{dt}} $$
*   **Units:** The Valencia 3-velocity components have units of length/time.

### Code Implementation


```python
# Define the Lorentz factor calculation function
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD**Tiny Doubles**
================

### Overview of Tiny Doubles

In this section, we will discuss the concept of tiny doubles.

### Theory Review

#### Introduction to Tiny Doubles

*   **Tiny Doubles:** Tiny doubles are a type of floating-point number that is used in numerical computations.
	+   They have a small value, typically on the order of 1e-100 or smaller.

### Code Implementation


```python
# Define TINYDOUBLE as a tiny double value
TINYDOUBLE = 1e-100

# Use TINYDOUBLE to avoid division by zero
R += TINYDOUBLE
```

This code defines `TINYDOUBLE` as a tiny double value and uses it to avoid division by zero in the calculation of `R`.

### Theory Review

#### Importance of Tiny Doubles

*   **Importance:** Tiny doubles are used to prevent numerical instability in calculations.
	+   They can be used to add a small amount to a value, preventing division by zero or other numerical errors.

### Mathematics


$$ \text{TINYDOUBLE} = 1e-100 $$

*   **Definition:** The tiny double `TINYDOUBLE` is defined as a small positive number (1e-100).

### Code Implementation


```python
# Define the TINYDOUBLE value as a constant
TINYDOUBLE = 1.0e-100

# Use TINYDOUBLE to avoid division by zero
R += TINYDOUBLE
```

This code defines `TINYDOUBLE` as a constant and uses it to add a small amount to `R`, preventing division by zero.

### Theory Review

#### Avoiding Division by Zero

*   **Importance:** Avoiding division by zero is crucial in numerical computations.
	+   Tiny doubles can be used to prevent this error from occurring.**Till the End of Time**
==========================

### Overview of Till the End of Time

In this section, we will discuss the concept of "till the end of time" and its relation to numerical computations.

### Theory Review

#### Introduction to Numerical Computations

*   **Numerical Computations:** Numerical computations are used to solve mathematical problems by using approximations.
	+   They are an essential tool in various fields, including physics, engineering, and computer science.

### Code Implementation


```python
# Define a function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define a function to calculate the Valencia 3-velocity components
def valencia_velocity(components):
    # Calculate the Valencia 3-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines two functions: `lorentz_factor` and `valencia_velocity`. The first function calculates the Lorentz factor using the formula: $$ R = \gamma_{ij} v^i v^j $$, while the second function calculates the Valencia 3-velocity components using the formula: $$ v_i = \frac{d x_i}{d t} $$.

### Theory Review

#### Importance of Numerical Computations

*   **Importance:** Numerical computations are essential in various fields, including physics, engineering, and computer science.
	+   They provide a way to solve mathematical problems by using approximations.

### Mathematics


$$ R = \gamma_{ij} v^i v^j $$

*   **Lorentz Factor:** The Lorentz factor is calculated using the formula: $$ R = \gamma_{ij} v^i v^j $$
	+   Where $\gamma_{ij}$ is the 3-metric, and $v^i**TINYDOUBLE Parameter**
==========================

### Overview of TINYDOUBLE Parameter

In this section, we will discuss the `TINYDOUBLE` parameter and its relation to numerical computations.

### Theory Review

#### Introduction to Numerical Computations

*   **Numerical Computations:** Numerical computations are used to solve mathematical problems by using approximations.
	+   They are an essential tool in various fields, including physics, engineering, and computer science.

### Code Implementation


```python
# Import the par module
import par

# Define TINYDOUBLE as a parameter
TINYDOUBLE = par.Cparameters("REAL", "GAMMA_SPEED_LIMIT", 1e-100)
```

This code defines `TINYDOUBLE` as a parameter using the `par.Cparameters()` function. The parameter is set to a small value of 1e-100.

### Theory Review

#### Importance of TINYDOUBLE Parameter

*   **Importance:** The `TINYDOUBLE` parameter is used to prevent numerical instability in calculations.
	+   It provides a way to add a small amount to values, preventing division by zero or other numerical errors.

### Mathematics


$$ \text{TINYDOUBLE} = 1e-100 $$

*   **Definition:** The `TINYDOUBLE` parameter is defined as a small positive number (1e-100).

### Code Implementation


```python
# Define the function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define the function to calculate the Valencia 3-velocity components
def valencia_velocity(components):
    # Calculate the Valencia 3-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines two functions: `lorentz_factor` and `valencia_velocity`. The first function calculates the Lorentz factor**Defining TINYDOUBLE**
=======================

### Overview of Defining TINYDOUBLE

In this section, we will discuss the definition of `TINYDOUBLE` as a parameter.

### Theory Review

#### Introduction to Parameters

*   **Parameters:** Parameters are values that are used in calculations and can be adjusted or modified.
	+   They provide a way to customize the behavior of algorithms and simulations.

### Code Implementation


```python
# Import the par module
import par

# Define TINYDOUBLE as a parameter
TINYDOUBLE = par.Cparameters("REAL", thismodule, "TINYDOUBLE", 1e-100)
```

This code defines `TINYDOUBLE` as a parameter using the `par.Cparameters()` function. The parameter is set to a small value of 1e-100.

### Theory Review

#### Importance of TINYDOUBLE Parameter

*   **Importance:** The `TINYDOUBLE` parameter is used to prevent numerical instability in calculations.
	+   It provides a way to add a small amount to values, preventing division by zero or other numerical errors.

### Mathematics


$$ \text{TINYDOUBLE} = 1e-100 $$

*   **Definition:** The `TINYDOUBLE` parameter is defined as a small positive number (1e-100).

### Code Implementation


```python
# Define the function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define the function to calculate the Valencia 3-velocity components
def valencia_velocity(components):
    # Calculate the Valencia 3-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines two functions: `lorentz_factor` and `valencia_velocity`. The first function calculates the Lorentz factor**Rescaled Valencia 3-Velocity**
=====================================

### Overview of Rescaling the Valencia 3-Velocity

In this section, we will discuss the concept of rescaling the Valencia 3-velocity.

### Theory Review

#### Introduction to Valencia 3-Velocity

*   **Valencia 3-Velocity:** The Valencia 3-velocity is a vector quantity that represents the velocity of an object in three-dimensional space.
	+   It's used to calculate various physical quantities, such as the Lorentz factor.

### Code Implementation


```python
# Define the function to rescale the Valencia 3-velocity
def rescaled_valencia_velocity(ValenciavU):
    # Rescale the Valencia 3-velocity using the formula: v' = GAMMA_SPEED_LIMIT * v / sqrt(GAMMA_SPEED_LIMIT^2 + v^2)
    v_prime = [0.0, 0.0, 0.0]
    for i in range(3):
        v_prime[i] = GAMMA_SPEED_LIMIT * ValenciavU[i] / math.sqrt(GAMMA_SPEED_LIMIT**2 + ValenciavU[i]**2)

    return v_prime
```

This code defines a function `rescaled_valencia_velocity` that takes the Valencia 3-velocity as input and returns its rescaled version.

### Theory Review

#### Rescaling the Valencia 3-Velocity

*   **Rescaling:** Rescaling is used to limit the speed of an object to a maximum value, known as the speed limit.
	+   In this case, we're using the Lorentz factor (GAMMA_SPEED_LIMIT) to rescale the Valencia 3-velocity.

### Mathematics


$$ v' = \frac{GAMMA\_SPEED\_LIMIT v}{\sqrt{GAMMA\_SPEED\_LIMIT^2 + v^2}} $$

*   **Rescaled Velocity:** The rescaled velocity is calculated using the formula above.
	+   Where `v` is the original Valencia 3-velocity, and `GAMMA_SPEED_LIMIT` is the Lorentz factor.

### Code Implementation


```python
# Define the function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * Valenciav**Rescaling the Valencia 3-Velocity**
=====================================

### Overview of Rescaling the Valencia 3-Velocity

In this section, we will discuss the process of rescaling the Valencia 3-velocity.

### Theory Review

#### Introduction to the Valencia 3-Velocity

*   **Valencia 3-Velocity:** The Valencia 3-velocity is a vector quantity that represents the velocity of an object in three-dimensional space.
    +   It's used to calculate various physical quantities, such as the Lorentz factor.

### Code Implementation


```python
# Define the function to rescale the Valencia 3-velocity
def rescale_valencia_velocity(Rstar, R, ValenciavU):
    # Rescale the Valencia 3-velocity using the formula: v_{(n)}^i = sqrt{Rstar/R} v^i
    global rescaledValenciavU

    # Initialize the rescaled Valencia 3-velocity array
    rescaledValenciavU = ixp.zerorank1(DIM=3)

    # Rescale each component of the Valencia 3-velocity
    for i in range(3):
        rescaledValenciavU[i] = math.sqrt(Rstar/R) * ValenciavU[i]

    return rescaledValenciavU

# Define the function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define the function to calculate the Valencia 3-velocity components
def valencia_velocity(components):
    # Calculate the Valencia 3-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines a function `rescale_valencia_velocity` that takes the Lorentz factor (`Rstar`) and the original Valencia 3-velocity (`ValenciavU`) as**Rescaling the Valencia 3-Velocity**
=====================================

### Overview of Rescaling the Valencia 3-Velocity

In this section, we will discuss the process of rescaling the Valencia 3-velocity.

### Theory Review

#### Introduction to the Valencia 3-Velocity

*   **Valencia 3-Velocity:** The Valencia 3-velocity is a vector quantity that represents the velocity of an object in three-dimensional space.
    +   It's used to calculate various physical quantities, such as the Lorentz factor.

### Code Implementation


```python
# Define the function to rescale the Valencia 3-velocity
def rescale_valencia_velocity(Rstar, R, ValenciavU):
    # Check if R is zero
    if R == 0:
        # If R is zero, then Rstar must also be zero
        Rstar = 0
    else:
        # Rescale the Valencia 3-velocity using the formula: v_{(n)}^i = sqrt{Rstar/R} v^i
        rescaledValenciavU = ixp.zerorank1(DIM=3)
        for i in range(3):
            rescaledValenciavU[i] = math.sqrt(Rstar/(R+TINYDOUBLE)) * ValenciavU[i]

    return rescaledValenciavU
```

This code defines a function `rescale_valencia_velocity` that takes the Lorentz factor (`Rstar`) and the original Valencia 3-velocity (`ValenciavU`) as input. If `R` is zero, it sets `Rstar` to zero.

### Theory Review

#### Rescaling the Valencia 3-Velocity

*   **Rescaling:** Rescaling is used to limit the speed of an object to a maximum value, known as the speed limit.
    +   In this case, we're using the Lorentz factor (GAMMA_SPEED_LIMIT) and `Rstar` to rescale the Valencia 3-velocity.

### Mathematics


$$ v_{(n)}^i = \sqrt{\frac{Rstar}{R}} v^i $$

*   **Rescaled Velocity:** The rescaled velocity is calculated using the formula above.
    +   Where `v` is the original Valencia 3-velocity, and `Rstar` is the Lorentz factor.

### Code Implementation


```python
**Physically Valid Velocities**
=============================

### Overview of Physically Valid Velocities

In this section, we will discuss the concept of physically valid velocities.

### Theory Review

#### Introduction to Physically Valid Velocities

*   **Physically Valid Velocities:** Physically valid velocities are those that are consistent with the laws of physics and do not result in unphysical or singular behavior.
    +   They must be carefully chosen and calculated to ensure that they accurately represent the physical system being modeled.

### Code Implementation


```python
# Define a function to check if a velocity is physically valid
def is_physically_valid(v):
    # Check if the velocity is of order 1e-100 or smaller
    return math.isclose(v, 0.0)

# Define a function to rescale velocities that are not physically valid
def rescale_velocity(v):
    # Rescale the velocity by adding TINYDOUBLE
    return v + TINYDOUBLE

# Define a function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define a function to calculate the Valencia 3-velocity components
def valencia_velocity(components):
    # Calculate the Valencia 3-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines three functions: `is_physically_valid`, `rescale_velocity`, and `lorentz_factor`. The first function checks if a velocity is physically valid by comparing it to zero. If the velocity is not physically valid, the second function rescales it by adding TINYDOUBLE.

### Theory Review

#### Importance of Physically Valid Velocities

*   **Importance:** Physically valid velocities are essential for accurate and reliable simulations.
    +   They must be carefully chosen and calculated to ensure that they accurately represent the physical system being modeled.

### Mathematics**Rescaling the Valencia 3-Velocity**
=====================================

### Overview of Rescaling the Valencia 3-Velocity

In this section, we will discuss the process of rescaling the Valencia 3-velocity.

### Theory Review

#### Introduction to Unit Conversion

*   **Unit Conversion:** Unit conversion is an important aspect of physical calculations. It ensures that the units are consistent and correct.
    +   Inaccurate unit conversion can lead to incorrect results.

### Code Implementation


```python
# Define a function to rescale the Valencia 3-velocity
def rescale_valencia_velocity(Rstar, R, ValenciavU):
    # Check if R is zero
    if R == 0:
        # If R is zero, then Rstar must also be zero
        Rstar = 0

    # Rescale the Valencia 3-velocity using the formula: v_{(n)}^i = sqrt{Rstar/R} v^i
    rescaledValenciavU = [0.0, 0.0, 0.0]
    for i in range(3):
        rescaledValenciavU[i] = ValenciavU[i] * sp.sqrt(Rstar/(R + TINYDOUBLE))

    return rescaledValenciavU
```

This code defines a function `rescale_valencia_velocity` that takes the Lorentz factor (`Rstar`) and the original Valencia 3-velocity (`ValenciavU`) as input.

### Theory Review

#### Rescaling the Valencia 3-Velocity

*   **Rescaling:** Rescaling is used to limit the speed of an object to a maximum value, known as the speed limit.
    +   In this case, we're using the Lorentz factor (GAMMA_SPEED_LIMIT) and `Rstar` to rescale the Valencia 3-velocity.

### Mathematics


$$ v_{(n)}^i = \sqrt{\frac{Rstar}{R}} v^i $$

*   **Rescaled Velocity:** The rescaled velocity is calculated using the formula above.
    +   Where `v` is the original Valencia 3-velocity, and `Rstar` is the Lorentz factor.

### Code Implementation


```python
# Define a function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using**Computing the Four-Velocity**
================================

### Overview of Computing the Four-Velocity

In this section, we will discuss how to compute the four-velocity `u^mu` in terms of the Valencia three-velocity components `Valenciav^i`.

### Theory Review

#### Introduction to the Four-Velocity

*   **Four-Velocity:** The four-velocity is a fundamental concept in special relativity, representing the velocity of an object in spacetime.
    +   It's defined as the derivative of the position with respect to time.

### Code Implementation


```python
# Define a function to compute the four-velocity
def compute_four_velocity(ValenciavU):
    # Compute u^mu using the formula: u^mu = gamma * Valenciav^i
    u_mu = [0.0, 0.0, 0.0]
    for i in range(3):
        u_mu[i] = GAMMA * ValenciavU[i]

    return u_mu

# Define a function to calculate the Lorentz factor
def lorentz_factor(ValenciavU):
    # Calculate the Lorentz factor using the formula: R = gammaDD[i][j] * ValenciavU[i] * ValenciavU[j]
    R = 0.0
    for i in range(3):
        for j in range(3):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    return R

# Define a function to calculate the Valencia three-velocity components
def valencia_velocity(components):
    # Calculate the Valencia three-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines a function `compute_four_velocity` that takes the Valencia three-velocity components as input and returns the four-velocity.

### Theory Review

#### Computing the Four-Velocity

*   **Four-Velocity:** The four-velocity is computed using the formula: $$ u^mu = gamma * Valenciav^i $$
    +   Where `gamma` is the Lorentz factor, and `Valenciav^i` are the Valencia three-velocity components.

### Mathematics


$$ u^**Computing the Four-Velocity Components**
==========================================

### Overview of Computing the Four-Velocity Components

In this section, we will discuss how to compute the four-velocity components `u^mu`.

### Theory Review

#### Introduction to the Four-Velocity Components

*   **Four-Velocity Components:** The four-velocity is a fundamental concept in special relativity, representing the velocity of an object in spacetime.
    +   It's defined as the derivative of the position with respect to time.

### Code Implementation


```python
# Define a function to compute the four-velocity components
def compute_four_velocity_components(Rstar, alpha):
    # Compute u^0 using the formula: u^0 = 1/(alpha-sqrt(1-R^*))
    global u4U_ito_ValenciavU

    # Initialize the four-velocity array
    u4U_ito_ValenciavU = ixp.zerorank1(DIM=4)

    # Compute each component of the four-velocity
    u4U_ito_ValenciavU[0] = 1/(alpha*sp.sqrt(1-Rstar))
    for i in range(3):
        u4U_ito_ValenciavU[i+1] = ValenciavU[i]

    return u4U_ito_ValenciavU

# Define a function to calculate the Valencia three-velocity components
def valencia_velocity(components):
    # Calculate the Valencia three-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines a function `compute_four_velocity_components` that takes the Lorentz factor (`Rstar`) and the Valencia three-velocity components as input and returns the four-velocity components.

### Theory Review

#### Computing the Four-Velocity Components

*   **Four-Velocity Components:** The four-velocity is computed using the formula: $$ u^mu = \left( \frac{1}{\alpha - \sqrt{1-R^*}}, 0, 0, 0 \right) $$
    +   Where `u^0` is the time component of the four-velocity.

### Mathematics


$$ u^0 = \frac{1}{\alpha - \**Computing the Four-Velocity Components**
==========================================

### Overview of Computing the Four-Velocity Components

In this section, we will discuss how to compute the four-velocity components `u^mu`.

### Theory Review

#### Introduction to the Four-Velocity Components

*   **Four-Velocity Components:** The four-velocity is a fundamental concept in special relativity, representing the velocity of an object in spacetime.
    +   It's defined as the derivative of the position with respect to time.

### Code Implementation


```python
# Define a function to compute the four-velocity components
def compute_four_velocity_components(Rstar, alpha, betaU):
    # Compute u^0 using the formula: u^0 = 1/(alpha-sqrt(1-R^*))
    global u4U_ito_ValenciavU

    # Initialize the four-velocity array
    u4U_ito_ValenciavU = ixp.zerorank1(DIM=4)

    # Compute each component of the four-velocity
    u4U_ito_ValenciavU[0] = 1/(alpha*sp.sqrt(1-Rstar))

    # Compute each spatial component of the four-velocity
    for i in range(3):
        u4U_ito_ValenciavU[i+1] = u4U_ito_ValenciavU[0] * (alpha * rescaledValenciavU[i] - betaU[i])

    return u4U_ito_ValenciavU

# Define a function to calculate the Valencia three-velocity components
def valencia_velocity(components):
    # Calculate the Valencia three-velocity components using the formula: v = dx/dt
    v = [0.0, 0.0, 0.0]
    for i in range(3):
        v[i] = components[i]/dt

    return v
```

This code defines a function `compute_four_velocity_components` that takes the Lorentz factor (`Rstar`), the alpha parameter (`alpha`), and the beta parameters (`betaU`) as input and returns the four-velocity components.

### Theory Review

#### Computing the Four-Velocity Components

*   **Four-Velocity Components:** The four-velocity is computed using the formula: $$ u^mu = \left( \frac{1}{\alpha - \sqrt{1-R^*}}, 0**Step 6.b: Converting Valencia 3-Velocity to Four-Velocity**
===========================================================

### Overview of Converting Valencia 3-Velocity to Four-Velocity

In this section, we will discuss how to convert the Valencia 3-velocity `v^i` into the four-velocity `u^\mu`, and apply a speed limiter.

### Theory Review

#### Introduction to Valencia 3-Velocity and Four-Velocity

*   **Valencia 3-Velocity:** The Valencia 3-velocity is a vector quantity that represents the velocity of an object in three-dimensional space.
    +   It's used to calculate various physical quantities, such as the Lorentz factor.

*   **Four-Velocity:** The four-velocity is a fundamental concept in special relativity, representing the velocity of an object in spacetime.
    +   It's defined as the derivative of the position with respect to time.

### Code Implementation


```python
# Define a function to convert Valencia 3-velocity to four-velocity
def valencia_to_four_velocity(v_i):
    # Compute u^\mu using the formula: u^mu = gamma * v^i + beta^i
    u_mu = [0.0, 0.0, 0.0]
    for i in range(3):
        u_mu[i] = GAMMA * v_i[i]

    return u_mu

# Define a function to apply speed limiter
def apply_speed_limiter(u_mu):
    # Apply the speed limiter using the formula: u^\mu \rightarrow max(u^\mu, -1)
    for i in range(4):
        if u_mu[i] > 0:
            u_mu[i] = min(u_mu[i], 1)

    return u_mu
```

This code defines two functions: `valencia_to_four_velocity` and `apply_speed_limiter`. The first function converts the Valencia 3-velocity into the four-velocity, and the second function applies a speed limiter to ensure that the four-velocity does not exceed a certain value.

### Theory Review

#### Converting Valencia 3-Velocity to Four-Velocity

*   **Converting Valencia 3-Velocity:** The Valencia 3-velocity is converted into the four-velocity using the formula: $$ u^mu = gamma * v^i + beta^i $$
    +   Where `gamma` is the**Speed-Limiting the Valencia 3-Velocity**
==========================================

### Overview of Speed-Limiting the Valencia 3-Velocity

In this section, we will discuss how to speed-limit the Valencia 3-velocity `v^i` using the formula: $$ vU = \frac{1}{\alpha} (v + \beta) $$. We will also review the code implementation and provide a detailed explanation of each step.

### Theory Review

#### Introduction to Speed-Limiting the Valencia 3-Velocity

*   **Speed-Limiting:** The speed-limiter is used to ensure that the velocity of an object does not exceed a certain value. This is important in relativistic calculations, where high speeds can lead to singularities.
    +   The speed-limiter uses the formula: $$ vU = \frac{1}{\alpha} (v + \beta) $$

#### Applying the Speed-Limiter

*   **Applying the Speed-Limiter:** To apply the speed-limiter, we need to first compute the Valencia 3-velocity `ValenciavU` using the formula: $$ ValenciavU_i = \frac{1}{\alpha} (vU_i + \beta_i) $$
    +   We then use this value to rescale the Valencia 3-velocity by applying the speed-limiter.

### Code Implementation


```python
# Define a function to rescale the Valencia 3-velocity using the speed limiter
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, vU):
    # Initialize the Valencia 3-velocity array
    ValenciavU = ixp.zerorank1(DIM=3)

    # Compute each component of the Valencia 3-velocity
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i])/alpha

    # Rescale the Valencia 3-velocity using the speed limiter
    u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)

# Define a function to rescale the Valencia 3-velocity in terms of u^mu
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed**Rescaling the Valencia 3-Velocity**
=====================================

### Overview of Rescaling the Valencia 3-Velocity

In this section, we will discuss how to rescale the Valencia 3-velocity `ValenciavU` using the formula: $$ ValenciavU = \frac{1}{\alpha} (v + \beta) $$. We will also review the code implementation and provide a detailed explanation of each step.

### Theory Review

#### Introduction to Rescaling the Valencia 3-Velocity

*   **Rescaling:** The rescaling process is used to adjust the Valencia 3-velocity components in order to ensure that they are physically meaningful.
    +   This is achieved by applying a speed limiter to prevent velocities from exceeding the speed of light.

### Code Implementation


```python
# Define a function to rescale the Valencia 3-velocity
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU):
    # Rescale each component of the Valencia 3-velocity using the formula: ValenciavU = \frac{1}{\alpha} (v + \beta)
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i])/alpha

# Define a function to rescale the vU array
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, vU):
    # Initialize the Valencia 3-velocity array
    ValenciavU = ixp.zerorank1(DIM=3)

    # Compute each component of the Valencia 3-velocity
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i])/alpha

    # Rescale the vU array using the formula: vU = \frac{\alpha}{\gamma} (ValenciavU + \beta)
    for i in range(4):
        vU[i] = (alpha/gamma) * (ValenciavU[i] + betaU[i])
```

This code defines two functions: `u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit` and `u4U_in_terms_of_vU__rescale**Renaming the Four-Velocity Array**
=====================================

### Overview of Renaming the Four-Velocity Array

In this section, we will discuss how to rename the four-velocity array `u4U_ito_ValenciavU` to `u4U_ito_vU`.

### Theory Review

#### Introduction to Renaming the Four-Velocity Array

*   **Renaming Arrays:** In programming, it's common to use temporary arrays or variables to store intermediate results. However, once the calculation is complete, it's often helpful to rename these arrays to something more descriptive.
    +   This can make the code easier to read and understand.

### Code Implementation


```python
# Define a global variable for the four-velocity array in terms of vU
global u4U_ito_vU

# Initialize the four-velocity array in terms of vU
u4U_ito_vU = ixp.zerorank1(DIM=4)

# Copy the values from the original four-velocity array to the new one
for mu in range(4):
    u4U_ito_vU[mu] = u4U_ito_ValenciavU[mu]
```

This code defines a global variable `u4U_ito_vU` and initializes it as an array of rank 1 with dimension 4. It then copies the values from the original four-velocity array `u4U_ito_ValenciavU` to the new array.

### Theory Review

#### Mathematics behind Renaming the Four-Velocity Array

*   **Mathematics:** The mathematics behind renaming the four-velocity array is simple assignment.
    +   We're simply copying the values from one array to another, which doesn't change their mathematical meaning.
    +   However, by renaming the array, we make it easier to understand and work with in subsequent calculations.

### Mathematics


$$ u4U_ito_vU = \left\{ \begin{array}{c} u^0 \\ u^1 \\ u^2 \\ u^3 \end{array} \right\} $$

*   **Four-Velocity Array:** The four-velocity array `u4U_ito_vU` represents the velocity of an object in spacetime.
    +   It's a fundamental concept in special relativity and is used to describe the motion of objects in various physical systems.**Computing the Rescaled Velocity**
=====================================

### Overview of Computing the Rescaled Velocity

In this section, we will discuss how to compute the rescaled velocity `rescaledvU` using the speed-limited Valencia 3-velocity components.

### Theory Review

#### Introduction to Speed-Limiting the Valencia 3-Velocity

*   **Speed-Limiting:** The speed limiter is used to ensure that the velocity of an object does not exceed a certain value.
    +   This is important in relativistic calculations, where high speeds can lead to singularities.

### Code Implementation


```python
# Declare global variables for rescaled velocity and Valencia 3-velocity
global rescaledvU
global rescaledValenciavU

# Initialize the rescaled velocity array
rescaledvU = ixp.zerorank1(DIM=3)

# Compute each component of the rescaled velocity using the formula: vU_i = alpha * ValenciavU_i - beta_i
for i in range(3):
    rescaledvU[i] = alpha * rescaledValenciavU[i] - betaU[i]
```

This code defines two global variables, `rescaledvU` and `rescaledValenciavU`, and initializes the rescaled velocity array using the formula: $$ vU_i = \alpha \cdot ValenciavU_i - \beta_i $$

### Theory Review

#### Mathematics behind Computing the Rescaled Velocity

*   **Mathematics:** The mathematics behind computing the rescaled velocity is simple algebra.
    +   We're simply applying the speed limiter to each component of the Valencia 3-velocity.

### Mathematics


$$ vU = \left\{ \begin{array}{c} \alpha \cdot ValenciavU_0 - \beta_0 \\ \alpha \cdot ValenciavU_1 - \beta_1 \\ \alpha \cdot ValenciavU_2 - \beta_2 \end{array} \right\} $$

*   **Rescaled Velocity:** The rescaled velocity represents the speed-limited velocity of an object in spacetime.
    +   It's a fundamental concept in special relativity and is used to describe the motion of objects in various physical systems.**Step 7: Constructing the GRHD Equations**
==========================================

### Overview of Constructing the GRHD Equations

In this section, we will discuss how to declare the ADM (Arnowitt-Deser-Misner) and hydrodynamical input variables and construct the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the dynamics of fluid motion in a curved spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Declare ADM input variables
ADM_input_vars = {
    'grav_pot': 0.0,
    'density': 1.0,
    'pressure': 1.0,
    'velocity': [0.0, 0.0]
}

# Declare hydrodynamical input variables
hydro_input_vars = {
    'enthalpy': 1.0,
    'temperature': 1.0,
    'entropy': 1.0
}

# Construct GRHD equations using ADM and hydrodynamical input variables
GRHD_equations = {
    'Hamiltonian constraint': lambda x: ADM_input_vars['grav_pot']**2 - (ADM_input_vars['density'] + ADM_input_vars['pressure'])/x,
    'Momentum constraint': lambda x: ADM_input_vars['velocity'][0] * ADM_input_vars['density'],
    'Energy equation': lambda x: ADM_input_vars['enthalpy'] * ADM_input_vars['temperature']
}
```

This code defines three dictionaries, `ADM_input_vars`, `hydro_input_vars`, and `GRHD_equations`, which contain the input variables for the ADM equations, hydrodynamical variables, and GRHD equations, respectively.

### Theory Review

#### Mathematics behind Constructing the GRHD Equations

*   **Mathematics:** The mathematics behind constructing the GRHD equations involves using the ADM equations to describe the dynamics of fluid motion in a curved spacetime.
    +   The GRHD equations are then constructed by combining the ADM equations with hydrodynamical variables, such as enthalpy and temperature.

### Mathematics


$$ \frac{\partial g_{\mu\nu}}{\partial x^\sigma} = -g_{\alpha\beta}**Declaring Variables and Constructing GRHD Equations**
=====================================================

### Overview of Declaring Variables and Constructing GRHD Equations

In this section, we will discuss how to declare the variables used in the General Relativistic Hydrodynamics (GRHD) equations and construct these equations.

### Theory Review

#### Introduction to GRHD Equations

*   **GRHD Equations:** The GRHD equations are a set of partial differential equations that describe the dynamics of fluid motion in a curved spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Declare variables for ADM input parameters
ADM_input_params = {
    'grav_pot': 0.0,
    'density': 1.0,
    'pressure': 1.0,
    'velocity': [0.0, 0.0]
}

# Declare variables for hydrodynamical input parameters
hydro_input_params = {
    'enthalpy': 1.0,
    'temperature': 1.0,
    'entropy': 1.0
}

# Construct GRHD equations using ADM and hydrodynamical input parameters
GRHD_equations = {
    'Hamiltonian constraint': lambda x: ADM_input_params['grav_pot']**2 - (ADM_input_params['density'] + ADM_input_params['pressure'])/x,
    'Momentum constraint': lambda x: ADM_input_params['velocity'][0] * ADM_input_params['density'],
    'Energy equation': lambda x: ADM_input_params['enthalpy'] * ADM_input_params['temperature']
}
```

This code defines three dictionaries, `ADM_input_params`, `hydro_input_params`, and `GRHD_equations`, which contain the input parameters for the ADM equations, hydrodynamical variables, and GRHD equations, respectively.

### Theory Review

#### Mathematics behind Declaring Variables and Constructing GRHD Equations

*   **Mathematics:** The mathematics behind declaring variables and constructing GRHD equations involves using the ADM equations to describe the dynamics of fluid motion in a curved spacetime.
    +   The GRHD equations are then constructed by combining the ADM equations with hydrodynamical variables, such as enthalpy and temperature.

### Mathematics


$$ \frac{\partial g_{\mu\nu}}{\partial x^\sigma} = -g_{\alpha\beta}\left( \**Defining Hydrodynamical Quantities**
=====================================

### Overview of Defining Hydrodynamical Quantities

In this section, we will discuss how to define the hydrodynamical quantities used in the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to Hydrodynamical Quantities

*   **Hydrodynamical Quantities:** The hydrodynamical quantities are fundamental variables that describe the behavior of a fluid in a curved spacetime.
    +   They include the energy density (`rho_b`), pressure (`P`), and internal energy per unit mass (`epsilon`).

### Code Implementation


```python
# Define the four-velocity as an array with 4 components
u4U = ixp.declarerank1("u4U", DIM=4)

# Define hydrodynamical quantities as symbols
rho_b, P, epsilon = sp.symbols('rho_b P epsilon', real=True)
```

This code defines the four-velocity `u4U` as an array with 4 components and declares the hydrodynamical quantities `rho_b`, `P`, and `epsilon` as symbols.

### Theory Review

#### Mathematics behind Defining Hydrodynamical Quantities

*   **Mathematics:** The mathematics behind defining hydrodynamical quantities involves using symbolic manipulation to represent the variables.
    +   This allows for easy manipulation and calculation of the hydrodynamical quantities in subsequent equations.

### Mathematics


$$ u^{\mu} = \left( \begin{array}{c} u^0 \\ u^1 \\ u^2 \\ u^3 \end{array} \right) $$

*   **Four-Velocity:** The four-velocity `u^{\mu}` represents the velocity of an object in spacetime.
    +   It is a fundamental concept in special relativity and is used to describe the motion of objects in various physical systems.

$$ \rho_b = \frac{P}{\epsilon} $$

*   **Energy Density:** The energy density `rho_b` represents the total energy per unit volume of a fluid.
    +   It is an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.**Defining ADM Quantities**
==========================

### Overview of Defining ADM Quantities

In this section, we will discuss how to define the Arnowitt-Deser-Misner (ADM) quantities used in the General Relativistic Hydrodynamics (GRHD) equations.

### Theory Review

#### Introduction to ADM Quantities

*   **ADM Quantities:** The ADM quantities are a set of variables that describe the behavior of a fluid in a curved spacetime.
    +   They include the gamma tensor (`gammaDD`), Kappa tensor (`KDD`), and beta vector (`betaU`).

### Code Implementation


```python
# Define the gamma tensor as a 2D array with 3x3 components
gammaDD = ixp.declarerank2("gammaDD", "sym01", DIM=3)

# Define the Kappa tensor as a 2D array with 3x3 components
KDD     = ixp.declarerank2("KDD",    "sym01", DIM=3)

# Define the beta vector as a 1D array with 3 components
betaU   = ixp.declarerank1("betaU", DIM=3)

# Declare alpha as a real symbol
alpha   = sp.symbols('alpha', real=True)
```

This code defines the ADM quantities `gammaDD`, `KDD`, and `betaU` using the `ixp.declarerank2` and `ixp.declarerank1` functions, which create 2D and 1D arrays with specified dimensions. The variable `alpha` is declared as a real symbol.

### Theory Review

#### Mathematics behind Defining ADM Quantities

*   **Mathematics:** The mathematics behind defining ADM quantities involves using symbolic manipulation to represent the variables.
    +   This allows for easy manipulation and calculation of the ADM quantities in subsequent equations.

### Mathematics


$$ \gamma_{ij} = \frac{\partial u^k}{\partial x^i} \cdot \frac{\partial u^l}{\partial x^j} $$

*   **Gamma Tensor:** The gamma tensor `γ_{ij}` represents the metric tensor of a fluid in a curved spacetime.
    +   It is an important quantity in general relativity and is used to describe the behavior of fluids under various conditions.

$$ K_{ij} = \frac{1}{2}**Computing Stress-Energy Tensor Components**
=============================================

### Overview of Computing Stress-Energy Tensor Components

In this section, we will discuss how to compute the stress-energy tensor components `T4UU` and `T4UD`.

### Theory Review

#### Introduction to Stress-Energy Tensor

*   **Stress-Energy Tensor:** The stress-energy tensor is a fundamental concept in general relativity that describes the distribution of mass and energy in spacetime.
    +   It's used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Compute T4UU (stress-energy tensor components)
def compute_T4UU(gammaDD, betaU, alpha, rho_b, P, epsilon, u4U):
    # Compute each component of the stress-energy tensor
    for mu in range(4):
        T4UU[mu] = 0.0
        for nu in range(4):
            T4UU[mu] += gammaDD[nu, mu] * rho_b

# Compute T4UD (stress-energy tensor components)
def compute_T4UD(gammaDD, betaU, alpha, T4UU):
    # Compute each component of the stress-energy tensor
    for nu in range(4):
        T4UD[nu] = 0.0
        for mu in range(4):
            T4UD[nu] += gammaDD[mu, nu] * (T4UU[mu] + betaU[mu])
```

This code defines two functions: `compute_T4UU` and `compute_T4UD`. The first function computes the stress-energy tensor components `T4UU`, while the second function computes the stress-energy tensor components `T4UD`.

### Theory Review

#### Mathematics behind Computing Stress-Energy Tensor Components

*   **Mathematics:** The mathematics behind computing stress-energy tensor components involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the stress-energy tensor in subsequent equations.

### Mathematics


$$ T^{\mu\nu} = \rho u^{\mu} u^{\nu} - g^{\mu\nu}\frac{P}{\epsilon} $$

*   **Stress-Energy Tensor:** The stress-energy tensor `T^{\mu\nu}` represents the distribution of mass and energy in spacetime.
    +   It's**Computing the Square Root of the Gamma Determinant**
=====================================================

### Overview of Computing the Square Root of the Gamma Determinant

In this section, we will discuss how to compute the square root of the gamma determinant `sqrt(gamma)`.

### Theory Review

#### Introduction to the Gamma Determinant

*   **Gamma Determinant:** The gamma determinant is a fundamental quantity in general relativity that describes the curvature of spacetime.
    +   It's used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Compute sqrt(gamma)
def compute_sqrtgammaDET(gammaDD):
    # Calculate the determinant of gammaDD
    det_gamma = sp.det(gammaDD)

    # Compute the square root of the determinant
    sqrt_gamma = sp.sqrt(det_gamma)
```

This code defines a function `compute_sqrtgammaDET` that calculates the determinant of the gamma tensor and computes its square root.

### Theory Review

#### Mathematics behind Computing the Square Root of the Gamma Determinant

*   **Mathematics:** The mathematics behind computing the square root of the gamma determinant involves using matrix algebra and calculus.
    +   This allows for easy manipulation and calculation of the gamma determinant in subsequent equations.

### Mathematics


$$ \sqrt{\gamma} = \left( \frac{1}{2\pi} \int_{0}^{\infty} e^{-s^2/4} ds \right) ^3 $$

*   **Gamma Determinant:** The gamma determinant `√γ` represents the square root of the curvature of spacetime.
    +   It's used to study various phenomena, such as black hole formation and binary mergers.

### Mathematical Operations


```python
import sympy as sp

# Define the symbols
s = sp.symbols('s')

# Calculate the integral
integral = sp.integrate(sp.exp(-s**2/4), (s, 0, sp.oo))

# Raise the result to the power of 3
result = integral ** 3
```

This code demonstrates how to calculate the gamma determinant using mathematical operations in the `sympy` library.

### Note:

*   The `sp.det()` function is used to compute the determinant of a matrix.
    +   The `sp.sqrt()` function is used to compute the square root of a number.**Computing Conservative Variables**
=====================================

### Overview of Computing Conservative Variables

In this section, we will discuss how to compute the conservative variables in terms of primitive variables.

### Theory Review

#### Introduction to Conservative Variables

*   **Conservative Variables:** The conservative variables are a set of variables that describe the behavior of a fluid in a curved spacetime.
    +   They include the energy density (`rho_star`), pressure (`tau_tilde`), and stress tensor (`S_tildeD`).

### Code Implementation


```python
# Compute rho_star (energy density)
def compute_rho_star(alpha, sqrtgammaDET, rho_b, u4U):
    # Calculate the energy density using the formula: rho_star = alpha * rho_b / sqrt(gamma)
    rho_star = alpha * rho_b / sqrtgammaDET

# Compute tau_tilde (pressure)
def compute_tau_tilde(alpha, sqrtgammaDET, T4UU, rho_star):
    # Calculate the pressure using the formula: tau_tilde = alpha^2 * T^tt / rho_star
    tau_tilde = alpha**2 * T4UU[0] / rho_star

# Compute S_tildeD (stress tensor)
def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    # Calculate the stress tensor using the formula: S^i_j = alpha * T^i_j
    S_tildeD = [alpha * T4UD[i] for i in range(3)]
```

This code defines three functions that compute the conservative variables `rho_star`, `tau_tilde`, and `S_tildeD` in terms of primitive variables.

### Theory Review

#### Mathematics behind Computing Conservative Variables

*   **Mathematics:** The mathematics behind computing conservative variables involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the conservative variables in subsequent equations.

### Mathematics


$$ \rho^* = \alpha \rho $$

*   **Energy Density:** The energy density `ρ^*` represents the total energy per unit volume of a fluid.
    +   It's an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.

$$ \tau^\tilde{} = \frac{\alpha^2 T^{tt}}{\rho^*} $$

*   **Pressure:** The pressure `τ^~` represents the normal stress**Computing v^i from u^μ**
=========================

### Overview of Computing v^i from u^μ

In this section, we will discuss how to compute the components `v^i` of the velocity vector `v` from the components `u^μ` of the four-velocity vector `u`.

### Theory Review

#### Introduction to Velocity Vector

*   **Velocity Vector:** The velocity vector `v` represents the velocity of a fluid element in spacetime.
    +   It's an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.

### Code Implementation


```python
# Compute v^i from u^μ (no speed limit)
def compute_vU_from_u4U__no_speed_limit(u4U):
    # Calculate the velocity vector components using the formula: v^i = u^0 * u^i / sqrt(gamma)
    vU = [u4U[0] * u4U[i] / math.sqrt(math.det(u4U)) for i in range(3)]
```

This code defines a function that computes the velocity vector components `v^i` from the four-velocity vector components `u^μ`.

### Theory Review

#### Mathematics behind Computing v^i from u^μ

*   **Mathematics:** The mathematics behind computing `v^i` from `u^μ` involves using the formula: $$v^i = \frac{u^0 u^i}{\sqrt{\gamma}}$$
    +   This allows for easy manipulation and calculation of the velocity vector components in subsequent equations.

### Mathematics


$$ v^i = \frac{u^0 u^i}{\sqrt{\gamma}} $$

*   **Velocity Vector:** The velocity vector `v` represents the velocity of a fluid element in spacetime.
    +   It's an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.

### Note:

*   In this code, we assume that the four-velocity vector components `u^μ` are already computed using the formula: $$u^\mu = \frac{dx^\mu}{d\tau}$$
    +   The velocity vector components `v^i` are then computed using the above formula.**Computing Fluxes of Conservative Variables**
=============================================

### Overview of Computing Fluxes of Conservative Variables

In this section, we will discuss how to compute the fluxes of conservative variables.

### Theory Review

#### Introduction to Fluxes of Conservative Variables

*   **Fluxes:** The fluxes are a set of quantities that describe the flow of conservative variables through a surface.
    +   They are used to study various phenomena, such as shock waves and vortex dynamics.

### Code Implementation


```python
# Compute rho_star_fluxU (energy density flux)
def compute_rho_star_fluxU(vU, rho_star):
    # Calculate the energy density flux using the formula: rho_star_fluxU^i = v^i * rho_star
    rho_star_fluxU = [vU[i] * rho_star for i in range(3)]

# Compute tau_tilde_fluxU (pressure flux)
def compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU, T4UU, rho_star):
    # Calculate the pressure flux using the formula: tau_tilde_fluxU^i = alpha * (v^i * tau_tilde + T^ti)
    tau_tilde_fluxU = [alpha * (vU[i] * T4UU[0] / rho_star + T4UD[i][0]) for i in range(3)]

# Compute S_tilde_fluxUD (stress tensor flux)
def compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD):
    # Calculate the stress tensor flux using the formula: S_tilde_fluxUD^ij = alpha * v^i * S^j_0
    S_tilde_fluxUD = [[alpha * vU[i] * T4UD[j][0] for j in range(3)] for i in range(3)]
```

This code defines three functions that compute the fluxes of conservative variables.

### Theory Review

#### Mathematics behind Computing Fluxes of Conservative Variables

*   **Mathematics:** The mathematics behind computing fluxes involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the fluxes in subsequent equations.

### Mathematics


$$ \rho^*_f = v_i \rho^* $$

*   **Energy Density Flux:** The energy density flux `ρ^*_f` represents the flow of energy per unit volume through a surface.
    +   It's an**Declaring Derivatives and Computing Zero-Time Derivative**
==========================================================

### Overview of Declaring Derivatives and Computing Zero-Time Derivative

In this section, we will discuss how to declare derivatives of various variables and compute the zero-time derivative of the metric tensor `g4DD`.

### Theory Review

#### Introduction to Derivatives

*   **Derivatives:** Derivatives are used to describe the rate of change of a function with respect to one or more variables.
    +   They are essential in physics and engineering to model various phenomena, such as motion and energy transfer.

### Code Implementation


```python
# Declare derivatives
gammaDD_dD = ixp.declarerank3("gammaDD_dD", "sym01", DIM=3)
betaU_dD = ixp.declarerank2("betaU_dD", "nosym", DIM=3)
alpha_dD = ixp.declarerank1("alpha_dD", DIM=3)

# Compute zero-time derivative of metric tensor
def compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha, gammaDD_dD, betaU_dD, alpha_dD):
    # Calculate the zero-time derivative of the metric tensor using the formula: g4DD^{0i}_{,j} = 0
    for i in range(3):
        for j in range(3):
            g4DD_zerotimederiv_dD[i][j] = gammaDD_dD[0][i][j]
```

This code declares derivatives of the metric tensor `gammaDD`, beta vector `betaU`, and alpha scalar `alpha` using the `ixp.declarerank3`, `ixp.declarerank2`, and `ixp.declarerank1` functions. It then computes the zero-time derivative of the metric tensor `g4DD`.

### Theory Review

#### Mathematics behind Declaring Derivatives and Computing Zero-Time Derivative

*   **Mathematics:** The mathematics behind declaring derivatives involves using tensor notation to represent the rate of change of a function with respect to one or more variables.
    +   This allows for easy manipulation and calculation of derivatives in subsequent equations.

### Mathematics


$$ \gamma^{ij}_{,k} = 0 $$

*   **Zero-Time Derivative:** The zero-time derivative `γ^{ij}_{,k}` represents**Computing Source Terms**
==========================

### Overview of Computing Source Terms

In this section, we will discuss how to compute the source terms on the `tau_tilde` and `S_tilde` equations.

### Theory Review

#### Introduction to Source Terms

*   **Source Terms:** The source terms represent additional contributions to the evolution equations of the fluid variables.
    +   They are used to describe various physical processes, such as viscosity and heat conduction.

### Code Implementation


```python
# Compute source term on tau_tilde equation
def compute_s_source_term(KDD, betaU, alpha, sqrtgammaDET, alpha_dD, T4UU):
    # Calculate the source term using the formula: s = KDD \* betaU + alpha^2 * (T^tt / rho_star)
    s = KDD[0][1] * betaU[1] + alpha**2 * T4UU[0] / sqrtgammaDET

# Compute source term on S_tilde equation
def compute_S_tilde_source_termD(alpha, sqrtgammaDET, g4DD_zerotimederiv_dD, T4UU):
    # Calculate the source term using the formula: S^i_j = alpha * (g^{0i}_{,j} + betaU^i \* KDD_j)
    S_tilde_source_termD = [alpha * (g4DD_zerotimederiv_dD[0][i][1] + betaU[i] * KDD[1][1]) for i in range(3)]
```

This code defines two functions that compute the source terms on the `tau_tilde` and `S_tilde` equations.

### Theory Review

#### Mathematics behind Computing Source Terms

*   **Mathematics:** The mathematics behind computing source terms involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the source terms in subsequent equations.

### Mathematics


$$ s = KDD \cdot betaU + alpha^2 \frac{T^{tt}}{\rho^*} $$

*   **Source Term on Tau_tilde Equation:** The source term `s` represents additional contributions to the evolution equation of the pressure variable.
    +   It is used to describe various physical processes, such as viscosity and heat conduction.

$$ S^i_j = alpha \left( g^{0i}_{,j}**Computing 4-Velocities**
=========================

### Overview of Computing 4-Velocities

In this section, we will discuss how to compute the 4-velocities in terms of an input Valencia 3-velocity.

### Theory Review

#### Introduction to 4-Velocities

*   **4-Velocities:** The 4-velocities are a set of four components that describe the velocity of an object in spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Define input Valencia 3-velocity
testValenciavU = ixp.declarerank1("testValenciavU", DIM=3)

# Compute 4-velocities in terms of input Valencia 3-velocity
def u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, testValenciavU):
    # Apply speed limit to input Valencia 3-velocity
    rescaled_valenciavU = [min(max(u, -1), 1) for u in testValenciavU]

    # Compute 4-velocities using formula: u^i = valencia^i / sqrt(gamma)
    u4U = [u / math.sqrt(math.det(gammaDD)) for u in rescaled_valenciavU]
```

This code defines a function that computes the 4-velocities in terms of an input Valencia 3-velocity.

### Theory Review

#### Mathematics behind Computing 4-Velocities

*   **Mathematics:** The mathematics behind computing 4-velocities involves using the formula: $$u^i = \frac{valencia^i}{\sqrt{\gamma}}$$
    +   This allows for easy manipulation and calculation of the 4-velocities in subsequent equations.

### Mathematics


$$ u^i = \frac{valencia^i}{\sqrt{\gamma}} $$

*   **4-Velocities:** The 4-velocities $u^i$ represent the velocity of an object in spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Note:

*   In this code, we apply a speed limit to the input Valencia 3-velocity using the formula: $$\**Computing 4-Velocities**
=========================

### Overview of Computing 4-Velocities

In this section, we will discuss how to compute the 4-velocities in terms of an input 3-velocity.

### Theory Review

#### Introduction to 4-Velocities

*   **4-Velocities:** The 4-velocities are a set of four components that describe the velocity of an object in spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Define input 3-velocity
testvU = ixp.declarerank1("testvU", DIM=3)

# Compute 4-velocities in terms of input 3-velocity
def u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, testvU):
    # Apply speed limit to input 3-velocity
    rescaled_vU = [min(max(u, -1), 1) for u in testvU]

    # Compute 4-velocities using formula: u^i = v^i / beta^0
    u4U = [u / betaU[0] for u in rescaled_vU]
```

This code defines a function that computes the 4-velocities in terms of an input 3-velocity.

### Theory Review

#### Mathematics behind Computing 4-Velocities

*   **Mathematics:** The mathematics behind computing 4-velocities involves using the formula: $$u^i = \frac{v^i}{\beta^0}$$
    +   This allows for easy manipulation and calculation of the 4-velocities in subsequent equations.

### Mathematics


$$ u^i = \frac{v^i}{\beta^0} $$

*   **4-Velocities:** The 4-velocities $u^i$ represent the velocity of an object in spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Note:

*   In this code, we apply a speed limit to the input 3-velocity using the formula: $$\text{max}\left(\text{min}(v^i, -1), 1\right)$$
    +   This**Code Validation**
====================

### Overview of Code Validation

In this section, we will discuss how to validate the code written in previous steps against the `GRHD.equations` NRPy+ module.

### Theory Review

#### Introduction to NRPy+

*   **NRPy+:** The NRPy+ module is a set of Python routines for numerically solving partial differential equations (PDEs) and ordinary differential equations (ODEs).
    +   It provides a flexible and efficient framework for implementing numerical methods in computational physics.
*   **`GRHD.equations`:** The `GRHD.equations` module within NRPy+ contains the equations of general relativistic hydrodynamics (GRHD).

### Code Implementation


```python
# Import necessary modules
import NRPyplus as nr

# Define variables and parameters
gammaDD = nr.declarerank2("gammaDD", 3, "sym01")
betaU = nr.declarerank1("betaU", 3)
alpha = nr.declarerank0("alpha")
u4U = nr.declarerank1("u4U", 4)

# Compute equations of GRHD
equations = ["g00DD", "gijDD", "T00UU", "TijUU"]

# Compare computed results with `GRHD.equations` module
for equation in equations:
    # Get the result from NRPy+
    nrpy_result = getattr(nr, equation)

    # Get the result from our code
    local_result = eval(equation)

    # Check if results match
    if nrpy_result == local_result:
        print(f"Equation {equation} passed validation.")
    else:
        print(f"Equation {equation} failed validation. Expected {nrpy_result}, but got {local_result}.")
```

This code defines a function that compares the computed results with the `GRHD.equations` module.

### Theory Review

#### Mathematics behind Code Validation

*   **Mathematics:** The mathematics behind code validation involves comparing the computed results with the expected results from NRPy+.
    +   This ensures that our implementation is correct and consistent with established numerical methods in computational physics.

### Mathematics


$$ g^{00}_{DD} = \frac{1}{\alpha^2} $$

*   **Equations of GRHD:** The equations of general relativistic hydrodynamics (GRHD) describe the**Code Validation**
====================

### Overview of Code Validation

In this section, we will discuss how to validate the code written in previous steps against the `GRHD.equations` NRPy+ module.

### Theory Review

#### Introduction to Code Validation

*   **Code Validation:** Code validation is the process of verifying that our implementation matches the expected results.
    +   It ensures that our code is correct and consistent with established numerical methods in computational physics.

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Verify agreement between SymPy expressions for GRHD equations
def validate_GRHD_equations():
    # Get the SymPy expressions for GRHD equations from NRPy+ module
    nrpy_expressions = [Ge.g00DD, Ge.gijDD, Ge.T00UU, Ge.TijUU]

    # Get the SymPy expressions for GRHD equations from our code
    local_expressions = ["g00DD", "gijDD", "T00UU", "TijUU"]

    # Compare the two sets of expressions
    if nrpy_expressions == local_expressions:
        print("GRHD equations match.")
    else:
        print("GRHD equations do not match.")

# Run code validation check
validate_GRHD_equations()
```

This code defines a function that compares the SymPy expressions for the GRHD equations generated in this tutorial versus the NRPy+ module.

### Theory Review

#### Mathematics behind Code Validation

*   **Mathematics:** The mathematics behind code validation involves comparing the SymPy expressions for the GRHD equations.
    +   This ensures that our implementation is correct and consistent with established numerical methods in computational physics.

### Mathematics


$$ g^{00}_{DD} = \frac{1}{\alpha^2} $$

*   **GRHD Equations:** The general relativistic hydrodynamics (GRHD) equations describe the behavior of fluid dynamics in a curved spacetime.
    +   They are used to study various phenomena, such as black hole formation and binary mergers.

### Note:

*   In this code, we import the `GRHD.equations` module from NRPy+ to access the SymPy expressions for the GRHD equations.
    +   We then compare these expressions with those generated in this tutorial.**Computing Stress-Energy Tensor**
==================================

### Overview of Computing Stress-Energy Tensor

In this section, we will discuss how to compute the stress-energy tensor components `T4UU` and `T4UD`.

### Theory Review

#### Introduction to Stress-Energy Tensor

*   **Stress-Energy Tensor:** The stress-energy tensor is a fundamental quantity in general relativity that describes the energy-momentum distribution of a fluid.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Define variables and parameters
gammaDD = ...  # metric tensor
betaU = ...     # beta vector
alpha = ...     # alpha scalar
rho_b = ...     # rest-mass density
P = ...         # pressure
epsilon = ...   # specific internal energy
u4U = ...       # four-velocity

# Compute stress-energy tensor components
def compute_stress_energy_tensor():
    # Compute T4UU (stress-energy tensor components in the time direction)
    Ge.compute_T4UU(gammaDD, betaU, alpha, rho_b, P, epsilon, u4U)

    # Compute T4UD (stress-energy tensor components in the spatial directions)
    Ge.compute_T4UD(gammaDD, betaU, alpha, Ge.T4UU)

# Run computation
compute_stress_energy_tensor()
```

This code defines a function that computes the stress-energy tensor components `T4UU` and `T4UD`.

### Theory Review

#### Mathematics behind Computing Stress-Energy Tensor

*   **Mathematics:** The mathematics behind computing stress-energy tensor involves using Einstein's field equations and the fluid variables.
    +   This allows for easy manipulation and calculation of the stress-energy tensor in subsequent equations.

### Mathematics


$$ T^{\mu \nu} = (\rho + P) u^\mu u^\nu - g^{\mu \nu} (P + \epsilon) $$

*   **Stress-Energy Tensor:** The stress-energy tensor `T^{\mu \nu}` represents the energy-momentum distribution of a fluid.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.

### Note:

*   In this code, we import the `GRHD.equations` module from NRPy**Computing Square Root of Metric Tensor**
=========================================

### Overview of Computing Square Root of Metric Tensor

In this section, we will discuss how to compute the square root of the metric tensor `sqrt(gamma)`.

### Theory Review

#### Introduction to Metric Tensor

*   **Metric Tensor:** The metric tensor is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Define variables and parameters
gammaDD = ...  # metric tensor

# Compute square root of metric tensor
def compute_sqrt_gamma():
    # Compute sqrt(gamma) using Ge.compute_sqrtgammaDET function
    sqrt_gamma = Ge.compute_sqrtgammaDET(gammaDD)

# Run computation
compute_sqrt_gamma()
```

This code defines a function that computes the square root of the metric tensor `sqrt(gamma)`.

### Theory Review

#### Mathematics behind Computing Square Root of Metric Tensor

*   **Mathematics:** The mathematics behind computing square root of metric tensor involves using the determinant of the metric tensor.
    +   This allows for easy manipulation and calculation of the square root in subsequent equations.

### Mathematics


$$ \sqrt{\gamma} = \sqrt{\det(g_{\mu \nu})} $$

*   **Square Root of Metric Tensor:** The square root of the metric tensor `√γ` represents the square root of the determinant of the metric tensor.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.

### Note:

*   In this code, we import the `GRHD.equations` module from NRPy+ to access the function `Ge.compute_sqrtgammaDET`.
    +  This function computes the square root of the determinant of the metric tensor.**Computing Conservative Variables**
=====================================

### Overview of Computing Conservative Variables

In this section, we will discuss how to compute the conservative variables in terms of primitive variables.

### Theory Review

#### Introduction to Conservative Variables

*   **Conservative Variables:** The conservative variables are a set of quantities that describe the behavior of a fluid in spacetime.
    +   They include the energy density (`rho_star`), pressure (`tau_tilde`), and stress tensor (`S_tildeD`).

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Define variables and parameters
alpha = ...  # alpha scalar
Ge.sqrtgammaDET = ...  # square root of determinant of metric tensor
rho_b = ...  # rest-mass density
u4U = ...  # four-velocity
Ge.T4UU = ...  # stress-energy tensor components in the time direction
Ge.T4UD = ...  # stress-energy tensor components in the spatial directions

# Compute conservative variables
def compute_conservative_variables():
    # Compute rho_star (energy density) using Ge.compute_rho_star function
    Ge.rho_star = Ge.compute_rho_star(alpha, Ge.sqrtgammaDET, rho_b, u4U)

    # Compute tau_tilde (pressure) using Ge.compute_tau_tilde function
    Ge.tau_tilde = Ge.compute_tau_tilde(alpha, Ge.sqrtgammaDET, Ge.T4UU, Ge.rho_star)

    # Compute S_tildeD (stress tensor) using Ge.compute_S_tildeD function
    Ge.S_tildeD = Ge.compute_S_tildeD(alpha, Ge.sqrtgammaDET, Ge.T4UD)

# Run computation
compute_conservative_variables()
```

This code defines a function that computes the conservative variables in terms of primitive variables.

### Theory Review

#### Mathematics behind Computing Conservative Variables

*   **Mathematics:** The mathematics behind computing conservative variables involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the conservative variables in subsequent equations.

### Mathematics


$$ \rho^* = \alpha \rho $$

*   **Energy Density:** The energy density `ρ^*` represents the total energy per unit volume of a fluid.
    +   It's an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.

$$ \tau^\til**Computing Velocity Components**
==================================

### Overview of Computing Velocity Components

In this section, we will discuss how to compute the velocity components `v^i` from the four-velocity vector `u`.

### Theory Review

#### Introduction to Four-Velocity Vector

*   **Four-Velocity Vector:** The four-velocity vector `u` represents the velocity of a fluid element in spacetime.
    +   It's an important quantity in hydrodynamics and is used to study the behavior of fluids under various conditions.

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Define variables and parameters
u4U = ...  # four-velocity vector

# Compute velocity components from four-velocity vector
def compute_velocity_components():
    # Compute v^i (velocity components) using Ge.compute_vU_from_u4U__no_speed_limit function
    vU = Ge.compute_vU_from_u4U__no_speed_limit(u4U)

# Run computation
compute_velocity_components()
```

This code defines a function that computes the velocity components `v^i` from the four-velocity vector `u`.

### Theory Review

#### Mathematics behind Computing Velocity Components

*   **Mathematics:** The mathematics behind computing velocity components involves using the formula: $$v^i = \frac{u^i}{u^0}$$
    +   This allows for easy manipulation and calculation of the velocity components in subsequent equations.

### Mathematics


$$ v^i = \frac{u^i}{u^0} $$

*   **Velocity Components:** The velocity components `v^i` represent the velocity of a fluid element in spacetime.
    +   They're an important quantity in hydrodynamics and are used to study the behavior of fluids under various conditions.

### Note:

*   In this code, we import the `GRHD.equations` module from NRPy+ to access the function `Ge.compute_vU_from_u4U__no_speed_limit`.
    +   This function computes the velocity components `v^i` from the four-velocity vector `u`.**Computing Fluxes of Conservative Variables**
=============================================

### Overview of Computing Fluxes of Conservative Variables

In this section, we will discuss how to compute the fluxes of conservative variables.

### Theory Review

#### Introduction to Conservative Variables

*   **Conservative Variables:** The conservative variables are a set of quantities that describe the behavior of a fluid in spacetime.
    +   They include the energy density (`rho_star`), pressure (`tau_tilde`), and stress tensor (`S_tildeD`).

### Code Implementation


```python
# Import necessary modules
import GRHD.equations as Ge

# Define variables and parameters
alpha = ...  # alpha scalar
Ge.sqrtgammaDET = ...  # square root of determinant of metric tensor
vU = ...  # velocity vector
Ge.T4UU = ...  # stress-energy tensor components in the time direction
Ge.rho_star = ...  # energy density
Ge.T4UD = ...  # stress-energy tensor components in the spatial directions

# Compute fluxes of conservative variables
def compute_fluxes():
    # Compute flux of rho_star (energy density) using Ge.compute_rho_star_fluxU function
    Ge.rho_star_fluxU = Ge.compute_rho_star_fluxU(Ge.vU, Ge.rho_star)

    # Compute flux of tau_tilde (pressure) using Ge.compute_tau_tilde_fluxU function
    Ge.tau_tilde_fluxU = Ge.compute_tau_tilde_fluxU(alpha, Ge.sqrtgammaDET, Ge.vU, Ge.T4UU, Ge.rho_star)

    # Compute flux of S_tildeD (stress tensor) using Ge.compute_S_tilde_fluxUD function
    Ge.S_tilde_fluxUD = Ge.compute_S_tilde_fluxUD(alpha, Ge.sqrtgammaDET, Ge.T4UD)

# Run computation
compute_fluxes()
```

This code defines a function that computes the fluxes of conservative variables.

### Theory Review

#### Mathematics behind Computing Fluxes of Conservative Variables

*   **Mathematics:** The mathematics behind computing fluxes of conservative variables involves using Einstein's field equations and the stress-energy tensor.
    +   This allows for easy manipulation and calculation of the fluxes in subsequent equations.

### Mathematics


$$ \rho^*_F = v^\mu \rho^* $$

*   **Flux of Energy Density:** The flux of energy density `ρ^*_F` represents the flow of energy**Declaring Derivatives and Computing g4DD**
===========================================

### Overview of Declaring Derivatives and Computing g4DD

In this section, we will discuss how to declare the derivatives and compute `g4DD`, which represents the zero-time-derivative of the metric tensor.

### Theory Review

#### Introduction to Metric Tensor and Its Derivatives

*   **Metric Tensor:** The metric tensor is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Derivatives of Metric Tensor:** The derivatives of the metric tensor represent the changes in its components with respect to time or spatial coordinates.

### Code Implementation


```python
# Import necessary modules
import NRPyplus as nr

# Declare derivatives
def declare_derivatives():
    # Declare d/dt (time derivative) and d/dx^i (spatial derivatives)
    dDdtd = nr.declarerank2("dDdtd", 3, "sym01")
    dDdxDD = [nr.declarerank2(f"dDdx{i}D", 3, "sym01") for i in range(3)]

# Compute g4DD (zero-time-derivative of metric tensor)
def compute_g4DD():
    # Compute g4DD using NRPy+ function
    g4DD = nr.compute_zero_timederiv_dD(gammaDD, dDdtd, dDdxDD)

    return g4DD

# Run computations
declare_derivatives()
g4DD = compute_g4DD()
```

This code defines functions to declare the derivatives and compute `g4DD`.

### Theory Review

#### Mathematics behind Declaring Derivatives and Computing g4DD

*   **Mathematics:** The mathematics behind declaring derivatives involves using NRPy+ syntax to declare the time derivative (`dDdtd`) and spatial derivatives (`dDdxDD`).
    +   This allows for easy manipulation and calculation of the derivatives in subsequent equations.

### Mathematics


$$ \partial_t g_{\mu\nu} = \frac{\partial}{\partial t} g_{\mu\nu} $$

*   **Zero-Time-Derivative of Metric Tensor:** The zero-time-derivative of the metric tensor `g4DD` represents the change in its**Declaring Derivatives of Metric Tensor**
=========================================

### Overview of Declaring Derivatives of Metric Tensor

In this section, we will discuss how to declare the derivatives of the metric tensor.

### Theory Review

#### Introduction to Metric Tensor and Its Derivatives

*   **Metric Tensor:** The metric tensor is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Derivatives of Metric Tensor:** The derivatives of the metric tensor represent the changes in its components with respect to time or spatial coordinates.

### Code Implementation


```python
# Import necessary modules
import ixp

# Declare derivatives of metric tensor
def declare_gammaDD_dD():
    # Declare gammaDD_dD (derivative of metric tensor) using ixp.declarerank3 function
    gammaDD_dD = ixp.declarerank3("gammaDD_dD", "sym01", DIM=3)

# Run declaration
declare_gammaDD_dD()
```

This code defines a function to declare the derivatives of the metric tensor.

### Theory Review

#### Mathematics behind Declaring Derivatives of Metric Tensor

*   **Mathematics:** The mathematics behind declaring derivatives involves using ixp syntax to declare the derivative of the metric tensor (`gammaDD_dD`).
    +   This allows for easy manipulation and calculation of the derivatives in subsequent equations.

### Mathematics


$$ \partial_\mu g_{\nu\rho} = \frac{\partial}{\partial x^\mu} g_{\nu\rho} $$

*   **Derivative of Metric Tensor:** The derivative of the metric tensor (`gammaDD_dD`) represents the change in its components with respect to a spatial coordinate.

### Note:

*   In this code, we use the `ixp.declarerank3` function to declare the derivative of the metric tensor.
    +   This function takes three arguments: the name of the variable, the symmetry property (`sym01`), and the dimensionality (`DIM=3`).**Declaring Derivatives of Beta Vector**
=====================================

### Overview of Declaring Derivatives of Beta Vector

In this section, we will discuss how to declare the derivatives of the beta vector.

### Theory Review

#### Introduction to Beta Vector and Its Derivatives

*   **Beta Vector:** The beta vector is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Derivatives of Beta Vector:** The derivatives of the beta vector represent the changes in its components with respect to time or spatial coordinates.

### Code Implementation


```python
# Import necessary modules
import ixp

# Declare derivatives of beta vector
def declare_betaU_dD():
    # Declare betaU_dD (derivative of beta vector) using ixp.declarerank2 function
    betaU_dD = ixp.declarerank2("betaU_dD", "nosym", DIM=3)

# Run declaration
declare_betaU_dD()
```

This code defines a function to declare the derivatives of the beta vector.

### Theory Review

#### Mathematics behind Declaring Derivatives of Beta Vector

*   **Mathematics:** The mathematics behind declaring derivatives involves using ixp syntax to declare the derivative of the beta vector (`betaU_dD`).
    +   This allows for easy manipulation and calculation of the derivatives in subsequent equations.

### Mathematics


$$ \partial_\mu \beta^\nu = \frac{\partial}{\partial x^\mu} \beta^\nu $$

*   **Derivative of Beta Vector:** The derivative of the beta vector (`betaU_dD`) represents the change in its components with respect to a spatial coordinate.

### Note:

*   In this code, we use the `ixp.declarerank2` function to declare the derivative of the beta vector.
    +   This function takes three arguments: the name of the variable, the symmetry property (`nosym`), and the dimensionality (`DIM=3`).**Declaring Derivatives and Computing g4DD**
============================================

### Overview of Declaring Derivatives and Computing g4DD

In this section, we will discuss how to declare the derivatives of `alpha`, `gammaDD`, and `betaU` and compute `g4DD`.

### Theory Review

#### Introduction to Alpha, Metric Tensor, and Beta Vector

*   **Alpha:** The alpha scalar is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Metric Tensor:** The metric tensor is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Beta Vector:** The beta vector is a fundamental quantity in general relativity that describes the geometry of spacetime.
    +   It plays a crucial role in understanding various phenomena, such as black hole formation and binary mergers.

### Code Implementation


```python
# Import necessary modules
import ixp
import GRHD.equations as Ge

# Declare derivatives of alpha, gammaDD, and betaU
def declare_derivatives():
    # Declare alpha_dD (derivative of alpha) using ixp.declarerank1 function
    alpha_dD = ixp.declarerank1("alpha_dD", DIM=3)

    # Declare gammaDD_dD (derivative of gammaDD) using ixp.declarerank3 function
    gammaDD_dD = ixp.declarerank3("gammaDD_dD", "sym01", DIM=3)

    # Declare betaU_dD (derivative of betaU) using ixp.declarerank2 function
    betaU_dD = ixp.declarerank2("betaU_dD", "nosym", DIM=3)

# Compute g4DD
def compute_g4DD():
    # Compute g4DD_zerotimederiv_dD (zero-time-derivative of metric tensor) using Ge.compute_g4DD_zerotimederiv_dD function
    g4DD = Ge.compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha, gammaDD_dD, betaU_dD, alpha_dD)

# Run declarations and**Computing Source Terms and Testing Expressions**
=============================================

### Overview of Computing Source Terms and Testing Expressions

In this section, we will discuss how to compute the source terms on `tau_tilde` and `S_tilde` equations.

### Theory Review

#### Introduction to Source Terms and Testing Expressions

*   **Source Terms:** The source terms are used to represent the changes in the conservative variables due to various physical processes.
    +   They play a crucial role in understanding various phenomena, such as black hole formation and binary mergers.
*   **Testing Expressions:** The testing expressions are used to verify that the computed values of the conservative variables match with the expected values.

### Code Implementation


```python
# Import necessary modules
import ixp
import GRHD.equations as Ge

# Compute source terms on tau_tilde and S_tilde equations
def compute_source_terms():
    # Compute s_source_term (source term on tau_tilde equation)
    Ge.s_source_term = Ge.compute_s_source_term(KDD, betaU, alpha, Ge.sqrtgammaDET, alpha_dD, Ge.T4UU)

    # Compute S_tilde_source_termD (source term on S_tilde equation)
    Ge.S_tilde_source_termD = Ge.compute_S_tilde_source_termD(alpha, Ge.sqrtgammaDET, Ge.g4DD_zerotimederiv_dD, Ge.T4UU)

# Test expressions
def test_expressions():
    # Declare test variables
    GetestValenciavU = ixp.declarerank1("testValenciavU")
    Ge.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, GetestValenciavU)

    GetestvU = ixp.declarerank1("testvU")
    Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha, betaU, gammaDD, GetestvU)
```

This code defines functions to compute the source terms on `tau_tilde` and `S_tilde` equations.

### Theory Review

#### Mathematics behind Computing Source Terms and Testing Expressions

*   **Mathematics:** The mathematics behind computing source terms involves using GRHD.equations module to compute the source terms on `tau_tilde` and `S_tilde` equations.
    +   This allows**Outputting to LaTeX-formatted PDF File**
==========================================

### Overview of Outputting to LaTeX-formatted PDF File

In this final step, we will output the notebook to a LaTeX-formatted PDF file.

### Theory Review

#### Introduction to LaTeX and PDF Generation

*   **LaTeX:** LaTeX is a typesetting system that is commonly used for producing high-quality documents, especially those with mathematical content.
    +   It provides a way to create complex layouts and formatting using a markup language.
*   **PDF Generation:** PDF generation refers to the process of converting a document into a Portable Document Format (PDF) file.
    +   This allows the document to be preserved in its original layout and formatting, regardless of the device or software used to view it.

### Code Implementation


```python
# Import necessary modules
from IPython.core.display import display, HTML

# Output notebook to LaTeX-formatted PDF file
def output_to_pdf():
    # Use nbconvert to convert notebook to LaTeX-formatted PDF file
    display(HTML("<script>console.log('Outputting notebook to LaTeX-formatted PDF file...');</script>"))
    print("Outputting notebook to LaTeX-formatted PDF file...")
```

This code defines a function to output the notebook to a LaTeX-formatted PDF file.

### Theory Review

#### Mathematics behind Outputting to LaTeX-formatted PDF File

*   **Mathematics:** The mathematics behind outputting to a LaTeX-formatted PDF file involves using the `nbconvert` module to convert the notebook into a LaTeX document.
    +   This document can then be rendered into a PDF file using a LaTeX compiler.

### Mathematics


$$ \LaTeX \to \text{PDF} $$

*   **Outputting Notebook to LaTeX-formatted PDF File:** The code above uses the `nbconvert` module to convert the notebook into a LaTeX document.
    +   This document is then rendered into a PDF file using a LaTeX compiler.**Converting Jupyter Notebook to LaTeX-formatted PDF File**
===========================================================

### Overview of Converting Jupyter Notebook to LaTeX-formatted PDF File

In this section, we will discuss how to convert a Jupyter notebook into a proper, clickable LaTeX-formatted PDF file.

### Theory Review

#### Introduction to Jupyter Notebooks and LaTeX Generation

*   **Jupyter Notebooks:** Jupyter notebooks are an interactive environment for working with code, visualizations, and documentation.
    +   They provide a way to create reproducible research by combining code, output, and explanations in a single document.
*   **LaTeX Generation:** LaTeX generation refers to the process of converting a notebook into a LaTeX document that can be rendered into a PDF file.

### Code Implementation


```python
# Import necessary modules
import cmdline_helper as cmd

# Convert Jupyter notebook to LaTeX-formatted PDF file
def convert_to_pdf():
    # Use cmdline_helper module to convert notebook to LaTeX-formatted PDF file
    cmd.run("jupyter nbconvert --to pdf Tutorial-GRHD_Equations-Cartesian.ipynb")
```

This code defines a function to convert the Jupyter notebook into a LaTeX-formatted PDF file.

### Theory Review

#### Mathematics behind Converting Jupyter Notebook to LaTeX-formatted PDF File

*   **Mathematics:** The mathematics behind converting a Jupyter notebook to a LaTeX-formatted PDF file involves using the `cmdline_helper` module to run a command that converts the notebook into a PDF file.
    +   This process uses the `jupyter nbconvert` command to convert the notebook into a LaTeX document, which is then rendered into a PDF file.

### Mathematics


$$ \text{Jupyter Notebook} \to \LaTeX \to \text{PDF} $$

*   **Converting Jupyter Notebook to LaTeX-formatted PDF File:** The code above uses the `cmdline_helper` module to run a command that converts the notebook into a PDF file.
    +   This process uses the `jupyter nbconvert` command to convert the notebook into a LaTeX document, which is then rendered into a PDF file.

### Note:

*   After running this code cell, the generated PDF file can be found in the root NRPy+ tutorial directory with filename [Tutorial-GRHD_Equations-Cartesian.pdf](Tutorial-GRHD_Equations-Cartesian.pdf).
    +   The link to the PDF file**Outputting Jupyter Notebook to LaTeX-formatted PDF File using NRPy+**
=====================================================================

### Overview of Outputting Jupyter Notebook to LaTeX-formatted PDF File using NRPy+

In this section, we will discuss how to output a Jupyter notebook to a LaTeX-formatted PDF file using the NRPy+ command-line interface.

### Theory Review

#### Introduction to NRPy+ and LaTeX Generation

*   **NRPy+:** NRPy+ is a multi-platform Python command-line interface for numerical relativity.
    +   It provides a way to compute and manipulate expressions in general relativity, including those related to black holes and gravitational waves.
*   **LaTeX Generation:** LaTeX generation refers to the process of converting a notebook into a LaTeX document that can be rendered into a PDF file.

### Code Implementation


```python
# Import necessary modules
import cmdline_helper as cmd

# Output Jupyter notebook to LaTeX-formatted PDF file using NRPy+
def output_to_latex_pdf():
    # Use cmd.output_Jupyter_notebook_to_LaTeXed_PDF function to convert notebook to LaTeX-formatted PDF file
    cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-GRHD_Equations-Cartesian")
```

This code defines a function to output the Jupyter notebook to a LaTeX-formatted PDF file using NRPy+.

### Theory Review

#### Mathematics behind Outputting Jupyter Notebook to LaTeX-formatted PDF File using NRPy+

*   **Mathematics:** The mathematics behind outputting a Jupyter notebook to a LaTeX-formatted PDF file using NRPy+ involves using the `cmd.output_Jupyter_notebook_to_LaTeXed_PDF` function.
    +   This function takes as input the name of the Jupyter notebook and outputs a LaTeX document that can be rendered into a PDF file.

### Mathematics


$$ \text{Jupyter Notebook} \to \LaTeX \to \text{PDF} $$

*   **Outputting Jupyter Notebook to LaTeX-formatted PDF File using NRPy+:** The code above uses the `cmd.output_Jupyter_notebook_to_LaTeXed_PDF` function to convert the notebook into a LaTeX document.
    +   This document is then rendered into a PDF file using a LaTeX compiler.

### Note:

*   After running this code cell, the generated PDF file can be found in the current working directory with filename `Tutorial-GRHD_Equations