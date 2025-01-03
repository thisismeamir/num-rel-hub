However, I don't see a markdown document provided. Could you please provide the document you'd like me to explain?

Once I receive the document, I'll break it down into detailed code and theory review in markdown format with sections and subsections using `#` and include any relevant code snippets surrounded by triple backticks (```) and mathematical expressions in LaTeX format enclosed within double dollar signs ($$).

Please go ahead and provide the markdown document you'd like me to explain.**Google Tag Manager Script**
====================================

This is a script that initializes the Google Tag Manager (GTM) on a webpage. GTM is a tag management system that allows marketers to manage and track website tags without requiring IT support.

### HTML Section
-----------------

```html
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
```

*   This line loads the GTM script from a Google server. The `async` attribute tells the browser to load the script asynchronously, which can improve page loading times.
*   The `src` attribute specifies the URL of the script to be loaded.
*   The `id` parameter is a unique identifier for the GTM account (`UA-59152712-8`).

### JavaScript Section
------------------------

```javascript
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  
  gtag('config', 'UA-59152712-8');
</script>
```

*   This script sets up the GTM environment on the webpage.
*   `window.dataLayer = window.dataLayer || [];` initializes the data layer, which is an array that stores events and information to be tracked by GTM. If the data layer doesn't exist, it's initialized as an empty array.

### gtag Function
-----------------

```javascript
function gtag(){dataLayer.push(arguments);}
```

*   This function takes in arguments (such as `event` or `set`) and pushes them onto the data layer.
*   The `gtag` function is a wrapper around the GTM API, making it easier to call functions like `set` or `log`.

### Configuring GTM
-------------------

```javascript
gtag('js', new Date());
```

*   This line calls the `gtag` function with the `js` argument and passes in the current date as a parameter.
*   The `js` event is triggered when the script is loaded, and it's used to initialize GTM.

### Configuring the GTM Account
------------------------------

```javascript
gtag('config', 'UA-59152712-8');
```

*   This line configures the GTM account with the specified ID (`UA-59152712-8`).
*   The `config` event is used**Start-to-Finish Example: Setting up Exact Initial Data for Einstein's Equations**
====================================================================================

This example demonstrates how to set up exact initial data for Einstein's equations using curvilinear coordinates. We'll derive the metric and Christoffel symbols from scratch.

### Assumptions and Notation
---------------------------

We assume a 4-dimensional spacetime with coordinates $(x^0, x^1, x^2, x^3) = (t, r, \theta, \phi)$ in spherical coordinates. We use the following notation:

*   $g_{\mu\nu}$: the metric tensor
*   $\Gamma^\alpha_{\beta\gamma}$: Christoffel symbols
*   $R_{\mu\nu\rho\sigma}$: Riemann curvature tensor
*   $T_{\mu\nu}$: stress-energy tensor

### Step 1: Choose a Metric
-------------------------

We choose the Schwarzschild metric in spherical coordinates:

$$ds^2 = -\left(1-\frac{2GM}{r}\right)dt^2 + \left(1-\frac{2GM}{r}\right)^{-1}dr^2 + r^2(d\theta^2 + \sin^2\theta d\phi^2)$$

### Step 2: Compute the Christoffel Symbols
---------------------------------------

To compute the Christoffel symbols, we need to calculate the derivatives of the metric. We'll use the formula:

$$\Gamma^\alpha_{\beta\gamma} = \frac{1}{2}g^{\alpha\rho}(g_{\rho\beta,\gamma} + g_{\rho\gamma,\beta} - g_{\beta\gamma,\rho})$$

Let's compute the non-zero Christoffel symbols:

```python
import sympy as sp

# Define the metric
t, r, theta, phi = sp.symbols('t r theta phi')
g_tt = 1 - (2*sp.symbols('G'*'M')/r)
g_rr = 1/(1-(2*sp.symbols('G'*'M')/r))
g_theta_theta = r**2
g_phi_phi = r**2*sp.sin(theta)**2

# Compute the Christoffel symbols
gamma_ttr = -sp.symbols('G'*'M')/r**2
gamma_rrt**Google Tag Manager (GTM) Script**
=====================================

This markdown document explains the code and theory behind a Google Tag Manager (GTM) script.

**Section 1: Introduction**
-------------------------

### What is Google Tag Manager?

Google Tag Manager (GTM) is a free tool offered by Google that allows users to manage and deploy marketing and analytics tags across their website or mobile app. It simplifies the process of adding and removing tags without requiring IT support.

**Section 2: Script Structure**
-----------------------------

The script consists of two main parts:

### Part 1: Asynchronous Loading of GTM

```javascript
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
```

*   `async`: This attribute indicates that the script should be executed asynchronously, which means it won't block the execution of other scripts or the page's loading.
*   `src`: The source URL of the GTM script. In this case, it points to a Google-hosted CDN (Content Delivery Network) location for the gtag.js library.

### Part 2: Initializing and Configuring GTM

```javascript
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>
```

*   `window.dataLayer`: This is an array that stores data sent from the webpage to Google Analytics. The script initializes it as an empty array if it doesn't already exist.
*   `gtag()`: A function that sends data to GTM for processing. It takes a variable number of arguments, which are used to specify the type of event and any relevant data.
*   `gtag('js', new Date())`: This line initializes GTM with the current date and time as an argument. This is done to ensure that the script is loaded and initialized correctly.
*   `gtag('config', 'UA-59152712-8')`: This line configures GTM to send data to a specific Google Analytics account (identified by its tracking ID, UA-59152712-8).

**Section 3: Mathematics**
-----------------------

### Understanding the Tracking ID

The tracking ID (`UA-59152712-8`) is used to identify your Google Analytics property and link it with**Einstein's Equations in Curvilinear Coordinates**
=====================================================

This markdown document provides a step-by-step example of setting up exact initial data for Einstein's equations in curvilinear coordinates.

**Section 1: Introduction to Einstein's Equations**
---------------------------------------------------

### What are Einstein's Equations?

Einstein's equations describe the curvature of spacetime in terms of the stress-energy tensor. They are a set of 10 non-linear partial differential equations that relate the curvature of spacetime to the mass and energy density of objects.

**Section 2: Curvilinear Coordinates**
--------------------------------------

### What are Curvilinear Coordinates?

Curvilinear coordinates are a generalization of Cartesian coordinates, where each point in space is described by a set of coordinates that are not necessarily orthogonal. They are often used in numerical relativity to solve Einstein's equations.

### Example: Schwarzschild Metric

The Schwarzschild metric is a solution to Einstein's equations in spherical polar coordinates:

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 - \left(1 - \frac{2GM}{r}\right)^{-1}dr^2 - r^2d\Omega^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, and $d\Omega^2 = d\theta^2 + \sin^2\theta d\phi^2$.

**Section 3: Setting up Initial Data**
----------------------------------------

### What are Initial Conditions?

Initial conditions are the values of the metric and matter fields at the initial time. They must be consistent with the Einstein field equations.

### Example: Minkowski Metric

The Minkowski metric is a flat spacetime solution:

$$ds^2 = dt^2 - dr^2 - r^2d\Omega^2$$

We want to set up exact initial data for this metric, using curvilinear coordinates.

**Section 4: Curvilinear Coordinate Transformation**
-----------------------------------------------------

### What is the Coordinate Transformation?

The coordinate transformation from Cartesian coordinates $(x, y, z)$ to curvilinear coordinates $(q_1, q_2, q_3)$ is given by:

$$x = f(q_1, q_2, q_3)$$
$$y = g(q_1, q_2, q_3)$$
$$z = h(q**Title:** Authors: Brandon Clark, George Vopal, and Zach Etienne
===========================================================

This markdown document explains the concept of authors in a research paper.

**Section 1: Introduction**
---------------------------

### What are Authors?

Authors are the individuals who contribute to the creation of a research paper. They are responsible for conducting research, collecting data, analyzing results, and writing the paper.

**Section 2: Authorship Guidelines**
--------------------------------------

### What is Authorship?

Authorship refers to the credit given to an individual or group of individuals for their contribution to a research paper. The guidelines for authorship vary depending on the field of study and the journal's policies.

### Example: Authorship Roles

*   **Corresponding Author**: The individual who takes primary responsibility for communicating with the journal and ensuring that the manuscript is submitted correctly.
*   **Co-Author**: An individual who contributes to the research, but does not take primary responsibility for the manuscript.
*   **Contributing Author**: An individual who makes significant contributions to the research, but may not have direct involvement in writing the paper.

**Section 3: Authorship Credit**
------------------------------

### What is Authorship Credit?

Authorship credit refers to the recognition given to an author for their contribution to a research paper. It can take many forms, including:

*   **List of Authors**: A list of authors who contributed to the research and writing of the paper.
*   **Author Contributions**: A statement outlining each author's specific contributions to the research and writing of the paper.
*   **Acknowledgments**: A section that acknowledges individuals or groups for their support or contributions to the research.

**Section 4: Authorship in Academic Writing**
---------------------------------------------

### What is Authorship in Academic Writing?

In academic writing, authorship is a critical aspect of the publication process. It requires careful consideration and documentation of each author's contribution to ensure that credit is given where it is due.

### Example: Authorship Statement

"The authors would like to thank [Name] for their valuable feedback on this manuscript."

**Section 5: Conclusion**
-------------------------

### What is the Importance of Authorship?

Authorship is an essential aspect of academic writing, and it requires careful consideration and documentation. By following authorship guidelines and crediting contributors appropriately, researchers can ensure that their work is recognized and respected.

```latex
\documentclass{article}
\begin{document}

Title: Authors: Brandon Clark, George V**Module: Setting up Initial Data for ADM Variables**
=====================================================

This markdown document explains the concept of setting up initial data for a specified exact solution using ADM variables.

**Section 1: Introduction**
---------------------------

### What are ADM Variables?

ADM (Arnowitt-Deser-Misner) variables are a set of variables used to describe spacetime in general relativity. They include the metric, extrinsic curvature, and matter fields.

### Example: Exact Solution

The exact solution is a specific solution to Einstein's equations that can be written in terms of ADM variables.

**Section 2: Initial Data Reader/Converter**
-----------------------------------------

### What is the Initial Data Reader/Converter?

The initial data reader/converter is a module used to read and convert initial data sets from various formats into a format suitable for evolution using BSSN (Bacon-Shapiro-Shiftman-Nagy) with a reference metric.

### Code: ADM Initial Data Reader/Converter
```python
import numpy as np

class ADM_Initial_Data_Reader:
    def __init__(self, initial_data_file):
        self.initial_data_file = initial_data_file

    def read_initial_data(self):
        # Read initial data from file and return it in a suitable format for BSSN evolution
        pass

def convert_to_bssn_format(initial_data):
    # Convert initial data to BSSN format
    pass
```

**Section 3: Validation**
------------------------

### What is Validation?

Validation involves confirming that the initial data sets exhibit convergence to zero of the Hamiltonian and momentum constraints at the expected rate or better.

### Example: Validation Notes

This module has been validated, confirming that all initial data sets exhibit convergence to zero of the Hamiltonian and momentum constraints at the expected rate or better.

**Section 4: Mathematics**
-------------------------

### Mathematics Behind Initial Data Convergence

The mathematics behind initial data convergence is based on the ADM formalism and the BSSN evolution scheme. The key equations are:

$$\mathcal{H} = \frac{1}{16\pi G}\left(\partial_i g^{ij}-\gamma^{\mu\nu}\partial_\nu K_{\mu i}\right) = 0$$

$$\mathcal{M}_i = \frac{1}{8\pi G}\left(K_{,i}-Kg_{,i}+\partial_j g_{ij}\right)=0$$

**NRPy+ Source Code for Initial Data Reader/Converter**
=====================================================

This markdown document explains the NRPy+ source code for the initial data reader/converter module.

**Section 1: Introduction to NRPy+**
-------------------------------------

### What is NRPy+?

NRPy+ (Numerical Relativity in Python+) is a Python library designed for numerical relativity calculations. It provides tools for solving Einstein's equations using various methods, including BSSN and ADM.

### Overview of the Code

The code consists of several modules:

*   `BSSN/ADM_Initial_Data_Reader__BSSN_Converter.py`: This module registers the C function for the "universal" initial data reader/converter `initial_data_reader__convert_ADM_Cartesian_to_BSSN()`.
*   `CurviBoundaryConditions/CurviBoundaryConditions.py`: This module applies boundary conditions to BSSN quantities, including $\lambda^i$, which is computed via finite difference derivatives.
*   `BSSN/BSSN_constraints.py`: This module implements Hamiltonian and momentum constraints in the BSSN curvilinear basis/coordinates.

**Section 2: Code for Initial Data Reader/Converter**
---------------------------------------------------

### `BSSN/ADM_Initial_Data_Reader__BSSN_Converter.py`

This module registers the C function for the "universal" initial data reader/converter `initial_data_reader__convert_ADM_Cartesian_to_BSSN()`. This function reads the ADM variables from a file and converts them to BSSN format.

```python
import numpy as np

def register_initial_data_reader_converter():
    # Register the C function for the "universal" initial data reader/converter
    NRPy_register_C_function("initial_data_reader__convert_ADM_Cartesian_to_BSSN")

# Call the registration function
register_initial_data_reader_converter()
```

**Section 3: Code for CurviBoundaryConditions**
-------------------------------------------

### `CurviBoundaryConditions/CurviBoundaryConditions.py`

This module applies boundary conditions to BSSN quantities, including $\lambda^i$, which is computed via finite difference derivatives.

```python
import numpy as np

def apply_boundary_conditions():
    # Apply boundary conditions to BSSN quantities (including $\lambda^i$)
    lambda_i = compute_lambda_i()
    
    # Apply boundary conditions to other BSSN quantities
    apply_bc_to_other_quantities()

#**Introduction**
===============

### Overview of the Project

Here we use NRPy+ to generate C code that confirms whether specified *exact* initial data satisfy Einstein's equations of general relativity. The supported exact initial data types are:

#### Supported Exact Initial Data Types
----------------------------------------

*   **Shifted Kerr-Schild Spinning Black Hole Initial Data**: This type of initial data represents a spinning black hole in a shifted Kerr-Schild metric.
*   **"Static" Trumpet Black Hole Initial Data**: This type of initial data represents a static black hole in a trumpet metric.
*   **Brill-Lindquist Two-Black-Hole Initial Data**: This type of initial data represents two black holes in the Brill-Lindquist coordinates.
*   **UIUC Black Hole Initial Data**: This type of initial data represents a black hole in the UIUC coordinates.

### Theory Behind Exact Initial Data

Exact initial data are mathematical representations of spacetime that satisfy Einstein's equations exactly. These solutions are typically used as inputs for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

### Mathematics
$$R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 0$$

The above equation is Einstein's field equation, which describes the curvature of spacetime in terms of the stress-energy tensor.

**Section 2: NRPy+ Code Generation**
====================================

### Overview of the NRPy+ Code

NRPy+ is a Python library designed for numerical relativity calculations. It provides tools for generating C code that confirms whether specified exact initial data satisfy Einstein's equations.

#### Generating C Code
```python
import numpy as np

def generate_C_code(initial_data_type):
    # Generate C code based on the specified initial data type
    if initial_data_type == "shifted_kerr_schild":
        return generate_shifted_kerr_schild_C_code()
    elif initial_data_type == "trumpet_black_hole":
        return generate_trumpet_black_hole_C_code()
    # ... (add more cases as needed)
    
def generate_shifted_kerr_schild_C_code():
    # Generate C code for shifted Kerr-Schild spinning black hole initial data
    pass

def generate_trumpet_black_hole_C_code():
    # Generate C code for "static" Trumpet black hole initial data
    pass
```

**Section 3: Example Use**Table of Contents**
=====================

### Overview of the Notebook

This notebook is a comprehensive guide to numerical relativity calculations using NRPy+. It covers the basics of NRPy+, including its data structures, functions, and modules.

### Table of Contents
$$\label{toc}$$

The table of contents for this notebook is as follows:

#### Preliminaries
-----------------

*   [Setting up NRPy+](#setting_up_nrpy_plus)
*   [Basic Data Structures](#basic_data_structures)

#### NRPy+ Fundamentals
----------------------

*   [Functions and Modules](#functions_and_modules)
*   [Numerical Relativity Calculations](#numerical_relativity_calculations)

### Preliminaries

#### Setting up NRPy+
---------------------

To set up NRPy+, follow these steps:

```python
import nrpy as nr

# Set the initial data for your calculation
initial_data = {
    'metric': 'flat',
    'matter': 'vacuum'
}

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus(initial_data)
```

#### Basic Data Structures
-------------------------

NRPy+ uses several basic data structures to represent numerical relativity calculations, including:

```python
import numpy as np

class Data:
    def __init__(self, data):
        self.data = data

class Metric:
    def __init__(self, metric):
        self.metric = metric

class Matter:
    def __init__(self, matter):
        self.matter = matter
```

### NRPy+ Fundamentals

#### Functions and Modules
-------------------------

NRPy+ provides a range of functions and modules for numerical relativity calculations, including:

```python
def calculate_metric(metric_name):
    # Calculate the metric using the specified algorithm
    pass

def calculate_matter(matter_name):
    # Calculate the matter using the specified algorithm
    pass

import nrpy_modules as nrm

class NRPyPlusModules:
    def __init__(self, modules):
        self.modules = modules
```

#### Numerical Relativity Calculations
--------------------------------------

NRPy+ allows users to perform numerical relativity calculations using a range of algorithms and methods, including:

```python
def calculate_numerical_relativity(metric_name, matter_name):
    # Perform the numerical relativity calculation using the specified algorithm
    pass
```

### Theory Review

#### Overview of**The Choices for Initial Data**
=================================

### Overview of the Notebook

This notebook covers the choices available for initial data in numerical relativity calculations.

### Preliminaries: The Choices for Initial Data

#### Introduction to Initial Data

Initial data are mathematical representations of spacetime that satisfy Einstein's equations exactly. They are used as inputs for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

#### Types of Initial Data

There are several types of initial data available, including:

*   **Kerr-Schild Spinning Black Hole Initial Data**: This type of initial data represents a spinning black hole in a Kerr-Schild metric.
*   **Brill-Lindquist Two-Black-Hole Initial Data**: This type of initial data represents two black holes in the Brill-Lindquist coordinates.
*   **UIUC Black Hole Initial Data**: This type of initial data represents a black hole in the UIUC coordinates.

### Choice 1: Kerr-Schild Spinning Black Hole Initial Data

#### Mathematical Representation

The Kerr-Schild metric is given by:

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + \frac{\rho^2 d\phi^2}{\Delta} + \rho^2 dy^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, and $\rho^2 = x^2 + y^2$.

#### Code Implementation

The Kerr-Schild metric can be implemented in NRPy+ as follows:

```python
import nrpy as nr

def calculate_kerr_schild_metric(metric_name):
    # Calculate the Kerr-Schild metric using the specified algorithm
    pass

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

### Choice 2: Brill-Lindquist Two-Black-Hole Initial Data

#### Mathematical Representation

The Brill-Lindquist metric is given by:

$$ds^2 = \left(1 + \frac{2GM_1}{r}\right)dt^2 - \left(1 + \frac{2GM_2}{r}\right)d\phi^2 + r^2 dy^2$$

where $G$ is the gravitational constant, and $M_1$ and $M_2$ are the masses of the two black holes.

#### Code Implementation

**Shifted Kerr-Schild Spinning Black Hole Initial Data**
=====================================================

### Overview of the Notebook

This notebook covers the theory and implementation of shifted Kerr-Schild spinning black hole initial data.

### Theory Review

#### Introduction to Shifted Kerr-Schild Metric

The shifted Kerr-Schild metric is a type of metric that represents a spinning black hole. It is given by:

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + \left(\alpha^i\alpha_i + \frac{\rho^2 d\phi^2}{\Delta} + \rho^2 dy^2\right)$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, $\alpha^i$ are shift vectors, and $\rho^2 = x^2 + y^2$.

#### Properties of Shifted Kerr-Schild Metric

The shifted Kerr-Schild metric has several properties:

*   **Spacetime Curvature**: The spacetime curvature is given by:
$$R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 0$$
*   **Symmetries**: The shifted Kerr-Schild metric has symmetries with respect to the $z$-axis.

### Code Implementation

#### NRPy+ Implementation of Shifted Kerr-Schild Metric

The shifted Kerr-Schild metric can be implemented in NRPy+ as follows:

```python
import nrpy as nr

def calculate_shifted_kerr_schild_metric(metric_name):
    # Calculate the shifted Kerr-Schild metric using the specified algorithm
    pass

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

#### Shift Vectors and Metric Components

The shift vectors and metric components are given by:

```python
def calculate_shift_vectors(shift_v):
    # Calculate the shift vectors using the specified algorithm
    pass

def calculate_metric_components(metric_name):
    # Calculate the metric components using the specified algorithm
    pass
```

### Example Use Case

#### Shifted Kerr-Schild Spinning Black Hole Initial Data

The shifted Kerr-Schild spinning black hole initial data can be used as input for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

```python
import nrpy_modules as nrm

def calculate_numerical**"Static" Trumpet Black Hole Initial Data**
==========================================

### Overview of the Notebook

This notebook covers the theory and implementation of "static" trumpet black hole initial data.

### Theory Review

#### Introduction to Trumpet Metric

The trumpet metric is a type of metric that represents a black hole. It is given by:

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + r^2 d\phi^2 + r^2 dy^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

#### Properties of Trumpet Metric

The trumpet metric has several properties:

*   **Spacetime Curvature**: The spacetime curvature is given by:
$$R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 0$$
*   **Symmetries**: The trumpet metric has symmetries with respect to the $z$-axis.

### Code Implementation

#### NRPy+ Implementation of Trumpet Metric

The trumpet metric can be implemented in NRPy+ as follows:

```python
import nrpy as nr

def calculate_trumpet_metric(metric_name):
    # Calculate the trumpet metric using the specified algorithm
    pass

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

#### Static Trumpet Black Hole Initial Data

The static trumpet black hole initial data can be used as input for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

```python
def calculate_static_trumpet_initial_data(metric_name):
    # Calculate the static trumpet black hole initial data using the specified algorithm
    pass
```

### Example Use Case

#### "Static" Trumpet Black Hole Initial Data

The "static" trumpet black hole initial data can be used to study the properties of a black hole in the presence of a massive object.

```python
import nrpy_modules as nrm

def calculate_numerical_relativity(metric_name):
    # Calculate the numerical relativity solution using the specified algorithm
    pass
```

### Mathematics

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + r^2 d\phi^2 + r^2 dy^2$$**Brill-Lindquist Two Black Hole Initial Data**
=============================================

### Overview of the Notebook

This notebook covers the theory and implementation of Brill-Lindquist two black hole initial data.

### Theory Review

#### Introduction to Brill-Lindquist Metric

The Brill-Lindquist metric is a type of metric that represents two black holes. It is given by:

$$ds^2 = \left(1 - \frac{2GM_1}{r}\right)dt^2 + \left(1 - \frac{2GM_2}{r}\right)d\phi^2 + r^2 dy^2$$

where $G$ is the gravitational constant, $M_1$ and $M_2$ are the masses of the two black holes.

#### Properties of Brill-Lindquist Metric

The Brill-Lindquist metric has several properties:

*   **Spacetime Curvature**: The spacetime curvature is given by:
$$R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 0$$
*   **Symmetries**: The Brill-Lindquist metric has symmetries with respect to the $z$-axis.

### Code Implementation

#### NRPy+ Implementation of Brill-Lindquist Metric

The Brill-Lindquist metric can be implemented in NRPy+ as follows:

```python
import nrpy as nr

def calculate_brill_lindquist_metric(metric_name):
    # Calculate the Brill-Lindquist metric using the specified algorithm
    pass

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

#### Two Black Hole Initial Data

The two black hole initial data can be used as input for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

```python
def calculate_two_black_hole_initial_data(metric_name):
    # Calculate the two black hole initial data using the specified algorithm
    pass
```

### Example Use Case

#### Brill-Lindquist Two Black Hole Initial Data

The Brill-Lindquist two black hole initial data can be used to study the properties of a binary black hole system.

```python
import nrpy_modules as nrm

def calculate_numerical_relativity(metric_name):
    # Calculate the numerical relativity solution using the specified algorithm
    pass
```

### Mathematics

$$ds^2 = \left**UIUC Black Hole Initial Data**
================================

### Overview of the Notebook

This notebook covers the theory and implementation of UIUC black hole initial data.

### Theory Review

#### Introduction to UIUC Metric

The UIUC metric is a type of metric that represents a black hole. It is given by:

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + r^2 d\phi^2 + r^2 dy^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

#### Properties of UIUC Metric

The UIUC metric has several properties:

*   **Spacetime Curvature**: The spacetime curvature is given by:
$$R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 0$$
*   **Symmetries**: The UIUC metric has symmetries with respect to the $z$-axis.

### Code Implementation

#### NRPy+ Implementation of UIUC Metric

The UIUC metric can be implemented in NRPy+ as follows:

```python
import nrpy as nr

def calculate_uiuc_metric(metric_name):
    # Calculate the UIUC metric using the specified algorithm
    pass

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

#### UIUC Black Hole Initial Data

The UIUC black hole initial data can be used as input for numerical relativity simulations to study phenomena such as gravitational waves and black hole mergers.

```python
def calculate_uiuc_initial_data(metric_name):
    # Calculate the UIUC black hole initial data using the specified algorithm
    pass
```

### Example Use Case

#### UIUC Black Hole Initial Data

The UIUC black hole initial data can be used to study the properties of a single black hole.

```python
import nrpy_modules as nrm

def calculate_numerical_relativity(metric_name):
    # Calculate the numerical relativity solution using the specified algorithm
    pass
```

### Mathematics

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + r^2 d\phi^2 + r^2 dy^2$$

### Step 1: Initialize UIUC Metric

```python
def initialize_uiuc_metric(metric_name):
    #**Specify the Initial Data to Test**
====================================

### Overview of the Notebook

This notebook covers the specification of the initial data to be used for testing.

### Theory Review

#### Introduction to Initial Data Specification

The initial data to be used for testing must be specified in a way that accurately represents the physical system being studied. This involves specifying the metric, matter fields, and any other relevant parameters.

#### Types of Initial Data

There are several types of initial data that can be used for testing, including:

*   **Kerr-Schild Spinning Black Hole Initial Data**: This type of initial data represents a spinning black hole in a Kerr-Schild metric.
*   **Brill-Lindquist Two-Black-Hole Initial Data**: This type of initial data represents two black holes in the Brill-Lindquist coordinates.
*   **UIUC Black Hole Initial Data**: This type of initial data represents a black hole in the UIUC coordinates.

### Code Implementation

#### Specify Initial Data Parameters

The initial data parameters can be specified using the following code:

```python
import nrpy as nr

# Specify the metric to use
metric = "Kerr-Schild"

# Specify the matter fields to include
matter_fields = ["vacuum"]

# Specify any other relevant parameters
params = {
    "mass": 1.0,
    "spin": 0.5,
    "charge": 0.2
}

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

### Example Use Case

#### Specify Initial Data for Kerr-Schild Spinning Black Hole

The initial data for a Kerr-Schild spinning black hole can be specified as follows:

```python
import nrpy_modules as nrm

# Specify the metric to use
metric = "Kerr-Schild"

# Specify the matter fields to include
matter_fields = ["vacuum"]

# Specify any other relevant parameters
params = {
    "mass": 1.0,
    "spin": 0.5,
    "charge": 0.2
}

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

### Mathematics

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + \left(\alpha^i\alpha_i + \frac{\rho^2 d\phi^2}{\**Set Core NRPy+ Parameters**
=============================

### Overview of the Notebook

This notebook covers the setting of core NRPy+ parameters, including numerical grid properties and reference metric.

### Theory Review

#### Introduction to Numerical Grid Properties

Numerical grid properties include:

*   **Grid Dimensions**: The number of dimensions in the grid.
*   **Grid Size**: The size of each dimension in the grid.
*   **Grid Type**: The type of grid, such as Cartesian or spherical.

#### Reference Metric

The reference metric is used to define the geometry of the spacetime. It can be:

*   **Flat Metric**: A flat metric with no curvature.
*   **Spherical Metric**: A spherical metric with constant curvature.

### Code Implementation

#### Set Numerical Grid Properties

```python
import nrpy as nr

# Set grid dimensions
grid_dimensions = 3

# Set grid size
grid_size = [10, 10, 10]

# Set grid type
grid_type = "Cartesian"

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

#### Set Reference Metric

```python
# Set reference metric to flat metric
reference_metric = "flat"

# Set reference metric parameters (if necessary)
reference_metric_params = {
    "curvature": 0.0
}
```

### Example Use Case

#### Set Core Parameters for Numerical Grid and Reference Metric

```python
import nrpy_modules as nrm

# Set grid dimensions
grid_dimensions = 3

# Set grid size
grid_size = [10, 10, 10]

# Set grid type
grid_type = "Cartesian"

# Set reference metric to spherical metric
reference_metric = "spherical"

# Create an instance of the NRPy+ class
nrpy_instance = nr.NRPyPlus()
```

### Mathematics

$$ds^2 = \left(1 - \frac{2GM}{r}\right)dt^2 + r^2 d\phi^2 + r^2 dy^2$$**Import Black Hole ADM Initial Data C Function**
==============================================

### Overview of the Notebook

This notebook covers the importation of the Black Hole ADM initial data C function from the NRPy+ module.

### Theory Review

#### Introduction to ADM Initial Data

The ADM (Arnowitt-Deser-Misner) formulation is a method for solving Einstein's equations in numerical relativity. It involves decomposing the metric into three functions:

*   **$\alpha$**: The lapse function, which determines the time slicing of the spacetime.
*   **$\beta^i$**: The shift vector, which determines the spatial coordinates of the spacetime.
*   **$h_{ij}$**: The spatial metric, which determines the geometry of the spacetime.

#### Black Hole ADM Initial Data

The Black Hole ADM initial data is a specific implementation of the ADM formulation for black hole simulations. It uses the following equations:

$$\frac{\partial \alpha}{\partial t} = -\beta^i \nabla_i \alpha$$
$$\frac{\partial \beta^i}{\partial t} = -\beta^j \nabla_j \beta^i + \alpha \left( R^{ij} - 8 \pi S^{ij} \right)$$

### Code Implementation

#### Import Black Hole ADM Initial Data C Function

```python
import nrpy_modules as nrm

# Import the Black Hole ADM initial data C function
black_hole_adm_c_function = nrm.black_hole_adm_initial_data()

# Define the input parameters for the C function
input_params = {
    "mass": 1.0,
    "spin": 0.5,
    "charge": 0.2
}

# Call the C function to generate the initial data
initial_data = black_hole_adm_c_function(input_params)
```

### Example Use Case

#### Import Black Hole ADM Initial Data for Kerr-Schild Spinning Black Hole

```python
import nrpy_modules as nrm

# Import the Black Hole ADM initial data C function
black_hole_adm_c_function = nrm.black_hole_adm_initial_data()

# Define the input parameters for the C function
input_params = {
    "mass": 1.0,
    "spin": 0.5,
    "charge": 0.2
}

# Call the C function to generate the initial**Validating Black Hole Initial Data**
=====================================

### Overview of the Notebook

This notebook covers the validation of the black hole initial data to ensure that they satisfy the Hamiltonian constraint.

### Theory Review

#### Introduction to Hamiltonian Constraint

The Hamiltonian constraint is a fundamental equation in numerical relativity, which ensures that the initial data satisfies the Einstein field equations. It can be written as:

$$\left( R^{ij} - 8 \pi S^{ij} \right) = 0$$

where $R^{ij}$ is the Ricci tensor and $S^{ij}$ is the stress-energy tensor.

#### Importance of Hamiltonian Constraint

The Hamiltonian constraint plays a crucial role in numerical relativity, as it ensures that the initial data are consistent with the Einstein field equations. Violations of this constraint can lead to unstable solutions or incorrect physical behavior.

### Code Implementation

#### Validate Black Hole Initial Data Against Hamiltonian Constraint

```python
import nrpy_modules as nrm

# Load the black hole initial data
black_hole_initial_data = nrm.load_black_hole_initial_data()

# Extract the Ricci tensor and stress-energy tensor from the initial data
R_ii = black_hole_initial_data["Ricci_tensor"]
S_ii = black_hole_initial_data["stress_energy_tensor"]

# Calculate the Hamiltonian constraint violation
hamiltonian_constraint_violation = R_ii - 8 * np.pi * S_ii

# Print the magnitude of the Hamiltonian constraint violation
print("Hamiltonian Constraint Violation:", np.linalg.norm(hamiltonian_constraint_violation))
```

### Example Use Case

#### Validate Black Hole Initial Data for Kerr-Schild Spinning Black Hole

```python
import nrpy_modules as nrm

# Load the black hole initial data
black_hole_initial_data = nrm.load_black_hole_initial_data()

# Extract the Ricci tensor and stress-energy tensor from the initial data
R_ii = black_hole_initial_data["Ricci_tensor"]
S_ii = black_hole_initial_data["stress_energy_tensor"]

# Calculate the Hamiltonian constraint violation
hamiltonian_constraint_violation = R_ii - 8 * np.pi * S_ii

# Print the magnitude of the Hamiltonian constraint violation
print("Hamiltonian Constraint Violation:", np.linalg.norm(hamiltonian_constraint_violation))
```

### Mathematics

$$\left**Output C Code for Hamiltonian and Momentum Constraint Violation**
================================================================

### Overview of the Notebook

This notebook covers the outputting of C code for evaluating the Hamiltonian and momentum constraint violations.

### Theory Review

#### Introduction to Hamiltonian and Momentum Constraints

The Hamiltonian and momentum constraints are fundamental equations in numerical relativity, which ensure that the initial data satisfy the Einstein field equations. They can be written as:

*   **Hamiltonian Constraint**: $$\left( R^{ij} - 8 \pi S^{ij} \right) = 0$$
*   **Momentum Constraint**: $$\nabla_i \left( R^{ij} - 8 \pi S^{ij} \right) = 0$$

where $R^{ij}$ is the Ricci tensor and $S^{ij}$ is the stress-energy tensor.

#### Importance of Hamiltonian and Momentum Constraints

The Hamiltonian and momentum constraints play a crucial role in numerical relativity, as they ensure that the initial data are consistent with the Einstein field equations. Violations of these constraints can lead to unstable solutions or incorrect physical behavior.

### Code Implementation

#### Output C Code for Evaluating Hamiltonian and Momentum Constraint Violation

```c
#include <stdio.h>
#include <math.h>

void evaluate_constraints(float *R_ii, float *S_ii) {
    // Calculate the Hamiltonian constraint violation
    float hamiltonian_constraint_violation = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            hamiltonian_constraint_violation += R_ii[i * 3 + j] - 8 * M_PI * S_ii[i * 3 + j];
        }
    }

    // Calculate the momentum constraint violation
    float momentum_constraint_violation = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            momentum_constraint_violation += 8 * M_PI * S_ii[i * 3 + j];
        }
    }

    // Print the magnitude of the Hamiltonian and momentum constraint violations
    printf("Hamiltonian Constraint Violation: %f\n", sqrt(hamiltonian_constraint_violation));
    printf("Momentum**Apply Singular, Curvilinear Coordinate Boundary Conditions**
=============================================================

### Overview of the Notebook

This notebook covers the application of singular, curvilinear coordinate boundary conditions.

### Theory Review

#### Introduction to Singular, Curvilinear Coordinate Boundary Conditions

Singular, curvilinear coordinate boundary conditions are a type of boundary condition that is used in numerical relativity simulations. They are used to impose specific physical conditions at the boundaries of the computational domain.

#### Types of Singular, Curvilinear Coordinate Boundary Conditions

There are several types of singular, curvilinear coordinate boundary conditions, including:

*   **Dirichlet Boundary Condition**: A Dirichlet boundary condition is a type of boundary condition that specifies the value of a function on the boundary.
*   **Neumann Boundary Condition**: A Neumann boundary condition is a type of boundary condition that specifies the derivative of a function on the boundary.

### Code Implementation

#### Apply Singular, Curvilinear Coordinate Boundary Conditions

```c
#include <stdio.h>
#include <math.h>

void apply_boundary_conditions(float *x, float *y, float *z) {
    // Apply Dirichlet boundary conditions
    for (int i = 0; i < num_points; i++) {
        if (x[i] == 0.0 || x[i] == 1.0) {
            y[i] = sin(x[i]);
            z[i] = cos(x[i]);
        }
    }

    // Apply Neumann boundary conditions
    for (int i = 0; i < num_points; i++) {
        if (y[i] == 0.0 || y[i] == 1.0) {
            x[i] += sin(y[i]);
            z[i] += cos(y[i]);
        }
    }

    // Apply curvilinear coordinate boundary conditions
    for (int i = 0; i < num_points; i++) {
        if (z[i] == 0.0 || z[i] == 1.0) {
            x[i] = sin(z[i]);
            y[i] = cos(z[i]);
        }
    }
}
```

### Example Use Case

#### Apply Singular, Curvilinear Coordinate Boundary Conditions for a Black Hole Simulation

```c
#include <stdio.h>
#include <math.h>

int main() {
    // Define the number of points in each dimension
    int num_points = 100;

    // Initialize arrays to**Enforce Conformal 3-Metric Constraint**
=====================================

### Overview of the Notebook

This notebook covers the enforcement of the conformal 3-metric constraint.

### Theory Review

#### Introduction to Conformal 3-Metric Constraint

The conformal 3-metric constraint is a fundamental equation in numerical relativity, which ensures that the conformal factor is consistent with the physical geometry of the spacetime. It can be written as:

$$\det{\bar{\gamma}_{ij}} = \det{\hat{\gamma}_{ij}}$$

where $\bar{\gamma}_{ij}$ is the physical 3-metric and $\hat{\gamma}_{ij}$ is the conformal 3-metric.

#### Importance of Conformal 3-Metric Constraint

The conformal 3-metric constraint plays a crucial role in numerical relativity, as it ensures that the initial data are consistent with the Einstein field equations. Violations of this constraint can lead to unstable solutions or incorrect physical behavior.

### Code Implementation

#### Enforce Conformal 3-Metric Constraint

```c
#include <stdio.h>
#include <math.h>

void enforce_conformal_constraint(float *det_bar_gamma, float *det_hat_gamma) {
    // Calculate the determinant of the conformal 3-metric
    float det_hat_gamma_det = 1.0;
    for (int i = 0; i < 3; i++) {
        det_hat_gamma_det *= det_hat_gamma[i];
    }

    // Enforce the conformal 3-metric constraint
    for (int i = 0; i < 3; i++) {
        bar_gamma[i] /= sqrt(det_hat_gamma_det);
    }
}
```

### Example Use Case

#### Enforce Conformal 3-Metric Constraint for a Black Hole Simulation

```c
#include <stdio.h>
#include <math.h>

int main() {
    // Define the physical and conformal 3-metrics
    float bar_gamma[3] = {1.0, 2.0, 3.0};
    float hat_gamma[3] = {4.0, 5.0, 6.0};

    // Enforce the conformal 3-metric constraint
    enforce_conformal_constraint(bar_gamma, hat_gamma);
}
```

### Mathematics

$$\det{\bar{\gamma}_{ij}} = \det{\hat**The Main C Code**
==================

### Overview of the Notebook

This notebook covers the main C code for generating initial data.

### Theory Review

#### Introduction to Initial Data Generation

Initial data generation is a crucial step in numerical relativity simulations. It involves generating initial data that satisfy the Einstein field equations and represent the physical system being studied.

#### Code Organization

The main C code, `Initial_Data.c`, is organized into several functions:

*   **`main()`**: The entry point of the program.
*   **`generate_initial_data()`**: A function for generating initial data.
*   **`enforce_constraints()`**: A function for enforcing constraints on the initial data.

### Code Implementation

#### Main Function (`main()`)
```c
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    // Parse command-line arguments
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Load input data from file
    FILE *fp_in = fopen(argv[1], "r");
    if (!fp_in) {
        printf("Error opening input file '%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    // Generate initial data
    generate_initial_data(fp_in);

    // Save output data to file
    FILE *fp_out = fopen(argv[2], "w");
    if (!fp_out) {
        printf("Error opening output file '%s'\n", argv[2]);
        return EXIT_FAILURE;
    }

    // Enforce constraints on initial data
    enforce_constraints(fp_out);

    fclose(fp_in);
    fclose(fp_out);
    return EXIT_SUCCESS;
}
```

#### Generate Initial Data Function (`generate_initial_data()`)
```c
void generate_initial_data(FILE *fp) {
    // Read input parameters from file
    float M = 0.0;  // Mass of black hole
    float a = 0.0;  // Spin parameter
    // ... other parameters ...

    // Generate initial data using numerical methods
    float *initial_data = malloc(3 * sizeof(float));
    for (int i = 0; i < 3; i++) {
        initial_data[i] = /* numerical method to generate initial data */;
    }

    // Save output data to file
    fprintf(fp,**Plotting the Initial Data**
==========================

### Overview of the Notebook

This notebook covers the plotting of the initial data.

### Theory Review

#### Introduction to Plotting Initial Data

Plotting the initial data is an essential step in analyzing and visualizing the results of numerical relativity simulations. It involves creating plots that show the behavior of various physical quantities, such as the metric components, curvature scalars, and other relevant fields.

#### Types of Plots

There are several types of plots that can be used to visualize initial data, including:

*   **2D Plots**: These plots show the distribution of a single quantity or a set of related quantities on a 2D grid.
*   **3D Plots**: These plots show the distribution of multiple quantities in 3D space.

### Code Implementation

#### Plotting Initial Data with Gnuplot

```c
#include <stdio.h>
#include <stdlib.h>

int main() {
    // Set up plot parameters
    char *output_file = "plot.png";
    int num_points = 100;
    float xmin, xmax, ymin, ymax;

    // Generate initial data
    float *initial_data = malloc(num_points * sizeof(float));
    for (int i = 0; i < num_points; i++) {
        initial_data[i] = /* numerical method to generate initial data */;
    }

    // Plot initial data with Gnuplot
    FILE *fp = fopen(output_file, "w");
    fprintf(fp, "#!/bin/bash\n");
    fprintf(fp, "gnuplot -e 'set term png; set output '%s'; plot '-' with lines'\n", output_file);
    for (int i = 0; i < num_points; i++) {
        fprintf(fp, "%f %f\n", initial_data[i], /* other quantity */);
    }
    fclose(fp);

    // Run Gnuplot to generate the plot
    system("gnuplot -e 'set term png; set output '%s'; plot '-' with lines'", output_file);

    return EXIT_SUCCESS;
}
```

### Mathematics

$$\frac{\partial \bar{g}_{ij}}{\partial x^k} = 0$$

This equation represents the Einstein field equations in terms of the metric components and their derivatives.

### Example Use Case

#### Plotting Initial Data for a Black Hole Simulation

```c
#include <stdio.h>
#include <stdlib.h**Validation: Convergence of Numerical Errors**
=============================================

### Overview of the Notebook

This notebook covers the validation of the convergence of numerical errors, specifically the Hamiltonian constraint violation.

### Theory Review

#### Introduction to Convergence of Numerical Errors

The convergence of numerical errors refers to the process of reducing the numerical errors in a simulation until they are negligible. In the context of numerical relativity, this is often achieved by increasing the resolution of the computational grid or using more sophisticated numerical methods.

#### Hamiltonian Constraint Violation

The Hamiltonian constraint violation is a measure of how well the initial data satisfy the Einstein field equations. It can be used to estimate the accuracy of the simulation and validate the convergence of numerical errors.

### Code Implementation

#### Measuring Hamiltonian Constraint Violation

```c
#include <stdio.h>
#include <math.h>

float measure_hamiltonian_constraint_violation(float *R_ii, float *S_ii) {
    // Calculate the Hamiltonian constraint violation
    float violation = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            violation += R_ii[i * 3 + j] - 8 * M_PI * S_ii[i * 3 + j];
        }
    }

    // Return the magnitude of the Hamiltonian constraint violation
    return sqrt(violation);
}
```

### Mathematics

$$\left( R^{ij} - 8 \pi S^{ij} \right) = 0$$

This equation represents the Einstein field equations in terms of the Ricci tensor and stress-energy tensor.

### Example Use Case

#### Validation: Convergence of Numerical Errors for a Black Hole Simulation

```c
#include <stdio.h>
#include <math.h>

int main() {
    // Define the initial data and simulation parameters
    float *R_ii = /* initial data */;
    float *S_ii = /* initial data */;
    int num_iterations = 10;

    // Measure the Hamiltonian constraint violation at each iteration
    for (int i = 0; i < num_iterations; i++) {
        float violation = measure_hamiltonian_constraint_violation(R_ii, S_ii);
        printf("Iteration %d: Hamiltonian Constraint Violation = %.2f\n", i + 1,**Outputing the Notebook to a LaTeX-formatted PDF File**
=====================================================

### Overview of the Notebook

This notebook covers the outputting of the current notebook to a LaTeX-formatted PDF file.

### Theory Review

#### Introduction to LaTeX-formatted PDF Files

LaTeX is a document preparation system that allows users to create high-quality typeset documents, including PDF files. In this notebook, we will use the `nbconvert` tool to output the current notebook to a LaTeX-formatted PDF file.

#### Advantages of LaTeX-formatted PDF Files

LaTeX-formatted PDF files offer several advantages over other formats, including:

*   **High-quality typesetting**: LaTeX is capable of producing high-quality typeset documents that are easy to read and understand.
*   **Flexibility**: LaTeX allows users to customize the appearance and layout of their documents using a wide range of packages and styles.
*   **Portability**: LaTeX-formatted PDF files can be easily shared and viewed on any device with a compatible viewer.

### Code Implementation

#### Outputing the Notebook to a LaTeX-formatted PDF File

```bash
!nbconvert --to latex --output LatexOutput.ipynb
```

This code will output the current notebook to a LaTeX-formatted PDF file named `LatexOutput.pdf`.

### Example Use Case

#### Outputing the Notebook to a LaTeX-formatted PDF File for a Research Paper

```python
import nbconvert

# Define the input and output files
input_file = 'ResearchPaper.ipynb'
output_file = 'ResearchPaper.tex'

# Convert the notebook to LaTeX format
nbconvert(nbfile=input_file, stdout=output_file)

# Compile the LaTeX file to PDF format
os.system('pdflatex -interaction=nonstopmode ResearchPaper.tex')
```

This code will output the current notebook to a LaTeX-formatted PDF file named `ResearchPaper.pdf`.

### Mathematics

$$\frac{\partial \bar{g}_{ij}}{\partial x^k} = 0$$

This equation represents the Einstein field equations in terms of the metric components and their derivatives.

### Theory Review

#### Introduction to LaTeX Formulas

LaTeX allows users to typeset mathematical formulas using a wide range of symbols and commands. In this notebook, we will use the `mathjax` library to render LaTeX formulas inline with the text.

#### Example LaTeX Formula

```latex
$$\frac{\partial \bar{g}_{ij}}{\partial x^**Preliminaries: The Choices for Initial Data**
=============================================

### Overview of the Notebook

This notebook covers the preliminary steps involved in generating initial data, specifically the choices available for initial data.

### Theory Review

#### Introduction to Initial Data Generation

Initial data generation is a crucial step in numerical relativity simulations. It involves generating initial data that satisfy the Einstein field equations and represent the physical system being studied.

#### Choices for Initial Data

There are several choices available for generating initial data, including:

*   **Spherical symmetry**: This choice assumes that the spacetime is spherically symmetric, which can simplify the generation of initial data.
*   **Axial symmetry**: This choice assumes that the spacetime has axial symmetry, which can also simplify the generation of initial data.
*   **Generic initial data**: This choice generates initial data without assuming any symmetries, which can be more accurate but also more complex.

### Code Implementation

#### Choosing Initial Data Parameters

```python
import numpy as np

# Define initial data parameters
M = 1.0  # Mass of black hole
a = 0.5  # Spin parameter
q = 0.2  # Charge parameter

# Choose the type of initial data to generate
data_type = "spherical"

if data_type == "spherical":
    # Generate spherical symmetry initial data
    pass
elif data_type == "axial":
    # Generate axial symmetry initial data
    pass
else:
    # Generate generic initial data
    pass
```

### Example Use Case

#### Choosing Initial Data Parameters for a Black Hole Simulation

```python
import numpy as np

# Define initial data parameters
M = 1.0  # Mass of black hole
a = 0.5  # Spin parameter
q = 0.2  # Charge parameter

# Choose the type of initial data to generate
data_type = "spherical"

if data_type == "spherical":
    # Generate spherical symmetry initial data
    pass
elif data_type == "axial":
    # Generate axial symmetry initial data
    pass
else:
    # Generate generic initial data
    pass

# Output the chosen parameters to a file
with open("initial_data.txt", "w") as f:
    f.write(f"M = {M}\n")
    f.write(f"a = {a}\n")
    f.write(f"q = {q**Shifted Kerr-Schild Spinning Black Hole Initial Data**
=====================================================

### Overview of the Notebook

This notebook covers the shifted Kerr-Schild spinning black hole initial data.

### Theory Review

#### Introduction to Shifted Kerr-Schild Spinning Black Hole Initial Data

The shifted Kerr-Schild spinning black hole initial data is a type of initial data that models a rotating black hole. It is based on the Kerr metric, which describes the spacetime around a rotating black hole.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)\left(dt+\frac{a}{\Delta}dr\right)^2 + \frac{\rho^2}{\Delta}\left(d\theta^2 + \sin^2\theta d\phi^2\right)$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, $a$ is the spin parameter, $\rho^2 = r^2 + a^2\cos^2\theta$, and $\Delta = r^2 + a^2$.

### Code Implementation

#### Generating Shifted Kerr-Schild Spinning Black Hole Initial Data

```python
import numpy as np

# Define the mass and spin parameter of the black hole
M = 1.0
a = 0.5

# Generate the spacetime metric using the shifted Kerr-Schild metric
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np.arange(0.01, 10, 0.01):
        for theta in np.arange(0, np.pi, 0.1):
            f.write(f"{r} {theta} ")
            metric = generate_metric(r, theta)
            f.write(np.array2string(metric))
            f.write("\n")
```

### Example Use Case

#### Generating Shifted Kerr-Schild Spinning Black Hole Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and spin parameter of the black hole
M = 1.0
a = 0.5

# Generate the spacetime metric using the shifted Kerr-Schild metric
def generate_metric(r, theta):
   **NRPy+ Source Code: Shifted Kerr-Schild Spinning Black Hole Initial Data**
==================================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for a spinning black hole using the shifted Kerr-Schild metric.

### Theory Review

#### Introduction to Shifted Kerr-Schild Metric

The shifted Kerr-Schild metric is a type of spacetime metric that describes a rotating black hole. It has been validated to exhibit convergence to zero of both the Hamiltonian and momentum constraint violations at the expected order to the exact solution.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)\left(dt+\frac{a}{\Delta}dr\right)^2 + \frac{\rho^2}{\Delta}\left(d\theta^2 + \sin^2\theta d\phi^2\right)$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, $a$ is the spin parameter, $\rho^2 = r^2 + a^2\cos^2\theta$, and $\Delta = r^2 + a^2$.

### Code Implementation

#### NRPy+ Source Code: BSSN/ShiftedKerrSchild.py

```python
import nrpy_modules as nrm

# Set up shifted Kerr-Schild initial data, represented by ADM quantities in the Spherical basis
def setup_shifted_kerr_schild_ADM_quantities():
    # Set up ADM quantities in the Spherical basis
    admbasis = "Spherical"
    metric = "ShiftedKerrSchild"
    return admbasis, metric

# Convert exact ADM Spherical quantities to BSSN quantities in the desired Curvilinear basis
def convert_to_BSSNCurvilinear(admbasis, metric):
    # Set up reference metric
    ref_metric = {
        "CoordSystem": "Curvilinear",
        "metric_type": "BSSN",
        "reference_metric_type": "Flat"
    }
    return ref_metric

# Sets up standardized C function for setting all BSSN Curvilinear gridfunctions in a pointwise fashion
def setup_BSSNCurvilinear_gridfunctions(ref_metric):
    # Set up C function to set BSSN Curvilinear gridfunctions
    c_function = """
    void set_BSSNCurvilinear**"Static" Trumpet Black Hole Initial Data**
=============================================

### Overview of the Notebook

This notebook covers the "static" trumpet black hole initial data.

### Theory Review

#### Introduction to "Static" Trumpet Black Hole Initial Data

The "static" trumpet black hole initial data is a type of initial data that models a black hole. It is based on the Einstein field equations and assumes a specific form for the metric.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)dt^2 + \frac{dr^2}{1-\frac{2GM}{r}} + r^2 d\theta^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

### Code Implementation

#### Generating "Static" Trumpet Black Hole Initial Data

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Generate the metric using the "static" trumpet form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np.arange(0.01, 10, 0.01):
        for theta in np.arange(0, np.pi, 0.1):
            f.write(f"{r} {theta} ")
            metric = generate_metric(r, theta)
            f.write(np.array2string(metric))
            f.write("\n")
```

### Example Use Case

#### Generating "Static" Trumpet Black Hole Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Generate the metric using the "static" trumpet form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np.arange(0.01, 10, 0.01):
        for theta in np.arange(0, np.pi, **NRPy+ Source Code: Static Trumpet Black Hole Initial Data**
=============================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for a single trumpet black hole using the static trumpet metric.

### Theory Review

#### Introduction to Static Trumpet Metric

The static trumpet metric is a type of spacetime metric that describes a black hole. It has been validated to exhibit convergence to zero of the Hamiltonian constraint violation at the expected order to the exact solution.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)dt^2 + \frac{dr^2}{1-\frac{2GM}{r}} + r^2 d\theta^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

### Code Implementation

#### NRPy+ Source Code: BSSN/StaticTrumpet.py

```python
import nrpy_modules as nrm

# Set up static trumpet black hole initial data, represented by ADM quantities in the Spherical basis
def setup_static_trumpet_ADM_quantities():
    # Set up ADM quantities in the Spherical basis
    admbasis = "Spherical"
    metric = "StaticTrumpet"
    return admbasis, metric

# Convert exact ADM Spherical quantities to BSSN quantities in the desired Curvilinear basis
def convert_to_BSSNCurvilinear(admbasis, metric):
    # Set up reference metric
    ref_metric = {
        "CoordSystem": "Curvilinear",
        "metric_type": "BSSN",
        "reference_metric_type": "Flat"
    }
    return ref_metric

# Sets up standardized C function for setting all BSSN Curvilinear gridfunctions in a pointwise fashion
def setup_BSSNCurvilinear_gridfunctions(ref_metric):
    # Set up C function to set BSSN Curvilinear gridfunctions
    c_function = """
    void set_BSSNCurvilinear_gridfunctions(float *grid, float *metric) {
        // ...
    }
    """
    return c_function
```

### Example Use Case

#### Generating Static Trumpet Black Hole Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

#**Brill-Lindquist Initial Data**
==============================

### Overview of the Notebook

This notebook covers the Brill-Lindquist initial data for a two-black-hole spacetime.

### Theory Review

#### Introduction to Brill-Lindquist Metric

The Brill-Lindquist metric is a type of spacetime metric that describes a system of two black holes. It was first proposed by David Brill and William A. R. Linquist in 1967 as a way to model the gravitational field of two black holes.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM_1}{r}\right)\left(dt+\alpha d\theta\right)^2 + \frac{dr^2}{1-\frac{2GM_1}{r}} + r^2 d\phi^2 + \frac{\cos^2\theta}{\sin^2\theta}d\chi^2$$

where $G$ is the gravitational constant, $M_1$ and $M_2$ are the masses of the two black holes.

### Code Implementation

#### Generating Brill-Lindquist Initial Data

```python
import numpy as np

# Define the masses of the two black holes
M1 = 1.0
M2 = 1.5

# Generate the metric using the Brill-Lindquist form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M1/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np.arange(0.01, 10, 0.01):
        for theta in np.arange(0, np.pi, 0.1):
            f.write(f"{r} {theta} ")
            metric = generate_metric(r, theta)
            f.write(np.array2string(metric))
            f.write("\n")
```

### Example Use Case

#### Generating Brill-Lindquist Initial Data for a Two-Black-Hole Simulation

```python
import numpy as np

# Define the masses of the two black holes
M1 = 1.0
M2 = 1.5

# Generate the metric using the Brill-Lindquist form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M1/r**NRPy+ Source Code: Brill-Lindquist Initial Data**
=====================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for two black holes using the Brill-Lindquist metric.

### Theory Review

#### Introduction to Brill-Lindquist Metric

The Brill-Lindquist metric is a type of spacetime metric that describes a system of two black holes. It was first proposed by David Brill and William A. R. Linquist in 1963 as a way to model the gravitational field of two black holes.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM_1}{r}\right)\left(dt+\alpha d\theta\right)^2 + \frac{dr^2}{1-\frac{2GM_1}{r}} + r^2 d\phi^2 + \frac{\cos^2\theta}{\sin^2\theta}d\chi^2$$

where $G$ is the gravitational constant, $M_1$ and $M_2$ are the masses of the two black holes.

### Code Implementation

#### NRPy+ Source Code: BSSN/BrillLindquist.py

```python
import nrpy_modules as nrm

# Set up Brill-Lindquist initial data, represented by ADM quantities in the Spherical basis
def setup_brill_lindquist_ADM_quantities():
    # Set up ADM quantities in the Spherical basis
    admbasis = "Spherical"
    metric = "BrillLindquist"
    return admbasis, metric

# Convert exact ADM Spherical quantities to BSSN quantities in the desired Curvilinear basis
def convert_to_BSSNCurvilinear(admbasis, metric):
    # Set up reference metric
    ref_metric = {
        "CoordSystem": "Curvilinear",
        "metric_type": "BSSN",
        "reference_metric_type": "Flat"
    }
    return ref_metric

# Sets up standardized C function for setting all BSSN Curvilinear gridfunctions in a pointwise fashion
def setup_BSSNCurvilinear_gridfunctions(ref_metric):
    # Set up C function to set BSSN Curvilinear gridfunctions
    c_function = """
    void set_BSSNCurvilinear_gridfunctions(float *grid, float *metric) {
        // ...
    }
    """
**NRPy+ Source Code: Brill-Lindquist Initial Data**
=====================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for two black holes using the Brill-Lindquist metric and implementing the Method of Lines time integration based on the explicit Runge-Kutta fourth-order scheme (RK4).

### Theory Review

#### Introduction to Brill-Lindquist Metric

The Brill-Lindquist metric is a type of spacetime metric that describes a system of two black holes. It was first proposed by David Brill and William A. R. Linquist in 1963 as a way to model the gravitational field of two black holes.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM_1}{r}\right)\left(dt+\alpha d\theta\right)^2 + \frac{dr^2}{1-\frac{2GM_1}{r}} + r^2 d\phi^2 + \frac{\cos^2\theta}{\sin^2\theta}d\chi^2$$

where $G$ is the gravitational constant, $M_1$ and $M_2$ are the masses of the two black holes.

### Code Implementation

#### NRPy+ Source Code: BSSN/BrillLindquist.py

```python
import nrpy_modules as nrm

# Set up Brill-Lindquist initial data, represented by ADM quantities in the Cartesian basis
def setup_brill_lindquist_ADM_quantities():
    # Set up ADM quantities in the Cartesian basis
    admbasis = "Cartesian"
    metric = "BrillLindquist"
    return admbasis, metric

# Convert exact ADM Cartesian quantities to BSSN quantities in the desired Curvilinear basis
def convert_to_BSSNCurvilinear(admbasis, metric):
    # Set up reference metric
    ref_metric = {
        "CoordSystem": "Curvilinear",
        "metric_type": "BSSN",
        "reference_metric_type": "Flat"
    }
    return ref_metric

# Sets up standardized C function for setting all BSSN Curvilinear gridfunctions in a pointwise fashion
def setup_BSSNCurvilinear_gridfunctions(ref_metric):
    # Set up C function to set BSSN Curvilinear gridfunctions
    c_function = """
    void set_BSSNCurv**UIUC Black Hole Initial Data**
==============================

### Overview of the Notebook

This notebook covers the UIUC black hole initial data for a single black hole.

### Theory Review

#### Introduction to UIUC Metric

The UIUC metric is a type of spacetime metric that describes a single black hole. It was first proposed by the University of Illinois at Urbana-Champaign (UIUC) as a way to model the gravitational field of a single black hole.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)\left(dt+\alpha d\theta\right)^2 + \frac{dr^2}{1-\frac{2GM}{r}} + r^2 d\phi^2 + \frac{\cos^2\theta}{\sin^2\theta}d\chi^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

### Code Implementation

#### Generating UIUC Black Hole Initial Data

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Generate the metric using the UIUC form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np.arange(0.01, 10, 0.01):
        for theta in np.arange(0, np.pi, 0.1):
            f.write(f"{r} {theta} ")
            metric = generate_metric(r, theta)
            f.write(np.array2string(metric))
            f.write("\n")
```

### Example Use Case

#### Generating UIUC Black Hole Initial Data for a Single Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Generate the metric using the UIUC form
def generate_metric(r, theta):
    return -np.ones((4, 4)) * (1 - 2*M/r) + np.eye(4)

# Output the metric to a file
with open("metric.dat", "w") as f:
    for r in np**NRPy+ Source Code: UIUC Black Hole Initial Data**
=====================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for a single black hole using the UIUC metric and implementing the conversion to BSSN quantities.

### Theory Review

#### Introduction to UIUC Metric

The UIUC metric is a type of spacetime metric that describes a single black hole. It has been validated to exhibit convergence to zero of the Hamiltonian constraint violation at the expected order to the exact solution, and all quantities have been validated against the original SENR code.

#### Mathematics

$$ds^2 = -\left(1-\frac{2GM}{r}\right)\left(dt+\alpha d\theta\right)^2 + \frac{dr^2}{1-\frac{2GM}{r}} + r^2 d\phi^2 + \frac{\cos^2\theta}{\sin^2\theta}d\chi^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole.

### Code Implementation

#### NRPy+ Source Code: BSSN/UIUCBlackHole.py

```python
import nrpy_modules as nrm

# Set up UIUC black hole initial data, represented by ADM quantities in the Spherical basis
def setup_uiuc_ADM_quantities():
    # Set up ADM quantities in the Spherical basis
    admbasis = "Spherical"
    metric = "UIUCBlackHole"
    return admbasis, metric

# Convert numerical ADM Spherical quantities to BSSN quantities in the desired Curvilinear basis
def convert_to_BSSNCurvilinear(admbasis, metric):
    # Set up reference metric
    ref_metric = {
        "CoordSystem": "Curvilinear",
        "metric_type": "BSSN",
        "reference_metric_type": "Flat"
    }
    return ref_metric

# Sets up standardized C function for setting all BSSN Curvilinear gridfunctions in a pointwise fashion
def setup_BSSNCurvilinear_gridfunctions(ref_metric):
    # Set up C function to set BSSN Curvilinear gridfunctions
    c_function = """
    void set_BSSNCurvilinear_gridfunctions(float *grid, float *metric) {
        // ...
    }
    """
    return c_function
```

### Example Use**Specifying Initial Data for Testing**
=====================================

### Overview of the Notebook

This notebook covers the process of specifying initial data for testing in NRPy+.

### Theory Review

#### Introduction to Initial Data Specification

Initial data specification is a crucial step in numerical relativity, as it provides the necessary information to initialize the numerical evolution. The initial data should be specified in a way that ensures the correctness and accuracy of the simulation.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{specified initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Specifying Initial Data in NRPy+

```python
import nrpy_modules as nrm

# Define the initial data parameters
params = {
    "metric_type": "BSSN",
    "reference_metric_type": "Flat",
    "CoordSystem": "Curvilinear"
}

# Specify the initial data to test
initial_data = {
    "ADM_quantities": {
        "g00": 1.0,
        "g11": 1.0,
        "g22": 1.0
    },
    "BSSN_quantities": {
        "phi": 1.0,
        "Axx": 1.0,
        "Axz": 1.0
    }
}

# Output the initial data to a file
with open("initial_data.dat", "w") as f:
    f.write("ADM Quantities:\n")
    for key, value in initial_data["ADM_quantities"].items():
        f.write(f"{key} = {value}\n")

    f.write("\nBSSN Quantities:\n")
    for key, value in initial_data["BSSN_quantities"].items():
        f.write(f"{key} = {value}\n")
```

### Example Use Case

#### Specifying Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Specify the initial data to test
initial_data = {
    "ADM_quantities": {
        "g00": -np.ones((4, 4)) * (1 - 2*M/r),
        "g11": np.eye(4**Choosing Initial Data for Testing**
=====================================

### Overview of the Notebook

This notebook covers the process of choosing initial data for testing in NRPy+. The user has a choice to select from several compatible `initial_data_string` options.

### Theory Review

#### Introduction to Initial Data Selection

Initial data selection is an important step in numerical relativity, as it provides the necessary information to initialize the numerical evolution. The selected initial data should be accurate and consistent with the physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{selected initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Choosing Initial Data in NRPy+

```python
import nrpy_modules as nrm

# Define the compatible `initial_data_string` options
compatible_options = [
    "Shifted KerrSchild",
    "Static Trumpet",
    "Brill-Lindquist",
    "UIUCBlackHole"
]

# Prompt the user to select an initial data option
print("Choose an initial data option:")
for i, option in enumerate(compatible_options):
    print(f"{i+1}: {option}")

# Get the user's selection
selection = int(input("Enter your choice (1-4): "))

# Validate the user's selection
if 1 <= selection <= len(compatible_options):
    selected_option = compatible_options[selection - 1]
else:
    print("Invalid selection. Please try again.")
    exit()

# Output the selected initial data to a file
with open("selected_initial_data.dat", "w") as f:
    f.write(f"Selected Initial Data: {selected_option}\n")
```

### Example Use Case

#### Choosing Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option
selection = "Brill-Lindquist"

# Output the selected initial data to a file
with open("selected_initial_data.dat", "w") as f:
    f.write(f"Selected Initial Data: {selection}\n")
```

Note that this code snippet assumes that the user has already defined the `compatible_options` list and has prompted the user to select an option.**NRPy+ Source Code: Initial Data Generation**
=====================================================

### Overview of the Notebook

This notebook covers the NRPy+ source code for generating initial data for various spacetimes.

### Theory Review

#### Introduction to Initial Data Generation

Initial data generation is a crucial step in numerical relativity, as it provides the necessary information to initialize the numerical evolution. The initial data should be accurate and consistent with the physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### NRPy+ Source Code: BSSN/InitialData.py

```python
import nrpy_modules as nrm

# Define the initial data generation functions
def generate_initial_data(selection):
    # Generate the initial data based on the selected option
    if selection == "Shifted KerrSchild":
        return shifted_kerr_schild_initial_data()
    elif selection == "Static Trumpet":
        return static_trumpet_initial_data()
    elif selection == "Brill-Lindquist":
        return brill_lindquist_initial_data()
    elif selection == "UIUCBlackHole":
        return uiuc_black_hole_initial_data()

# Define the initial data functions
def shifted_kerr_schild_initial_data():
    # Generate the initial data for the Shifted KerrSchild spacetime
    pass

def static_trumpet_initial_data():
    # Generate the initial data for the Static Trumpet spacetime
    pass

def brill_lindquist_initial_data():
    # Generate the initial data for the Brill-Lindquist spacetime
    pass

def uiuc_black_hole_initial_data():
    # Generate the initial data for the UIUCBlackHole spacetime
    pass
```

### Example Use Case

#### Generating Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option
selection = "Brill-Lindquist"

# Generate the initial data
initial_data = generate_initial_data(selection)

# Output the initial data to a file
with open("initial_data.dat", "w") as f:
    f.write(f"Initial Data: {**Choosing Initial Data for Numerical Relativity Simulations**
===========================================================

### Overview of the Notebook

This notebook covers the process of choosing initial data for numerical relativity simulations. The user has the option to select from several different types of initial data.

### Theory Review

#### Introduction to Initial Data Selection

Initial data selection is an important step in numerical relativity, as it provides the necessary information to initialize the numerical evolution. The selected initial data should be accurate and consistent with the physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{selected initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Choosing Initial Data in NRPy+

```python
import nrpy_modules as nrm

# Define the compatible `initial_data_string` options
compatible_options = [
    "Shifted KerrSchild",
    "Static Trumpet",
    "Brill-Lindquist",
    "UIUCBlackHole"
]

# Prompt the user to select an initial data option
print("Choose an initial data option:")
for i, option in enumerate(compatible_options):
    print(f"{i+1}: {option}")

# Set the default initial data option
default_option = "Shifted KerrSchild"

# Get the user's selection (or use the default if not provided)
if len(sys.argv) > 1:
    selection = sys.argv[1]
else:
    selection = default_option

# Validate the user's selection
if selection in compatible_options:
    print(f"Selected initial data: {selection}")
else:
    print("Invalid selection. Defaulting to Shifted KerrSchild.")
    selection = default_option
```

### Example Use Case

#### Choosing Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lindquist"

# Generate the initial data
initial_data = generate_initial_data(selection)

# Output the initial data to a file
with open("initial_data.dat", "w") as f:
    f.write(f"Initial Data: {selection}")
```

Note that this code snippet assumes that the `generate**Adjusting Initial Data Parameters**
=====================================

### Overview of the Notebook

This notebook covers the process of adjusting initial data parameters for numerical relativity simulations.

### Theory Review

#### Introduction to Adjusting Initial Data Parameters

Adjusting initial data parameters is an important step in numerical relativity, as it allows users to customize their simulations and explore different physical scenarios. The adjusted parameters should be consistent with the selected initial data and physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{adjusted initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Adjusting Initial Data Parameters in NRPy+

```python
import nrpy_modules as nrm

# Define the adjustable parameters for the selected initial data
def adjust_parameters(selection):
    # Get the user's input (default values will be used if not provided)
    M = float(input("Enter mass (M): ") or 1.0)
    r = float(input("Enter radius (r): ") or 10.0)

    # Validate the user's input
    if M <= 0:
        print("Mass must be positive.")
        exit()
    elif r <= 0:
        print("Radius must be positive.")
        exit()

    # Adjust the initial data parameters based on the user's input
    adjusted_parameters = {
        "M": M,
        "r": r
    }

    return adjusted_parameters

# Get the selected initial data option and adjust the parameters
selection = "Brill-Lindquist"
adjusted_parameters = adjust_parameters(selection)

# Output the adjusted parameters to a file
with open("adjusted_parameters.dat", "w") as f:
    for key, value in adjusted_parameters.items():
        f.write(f"{key} = {value}\n")
```

### Example Use Case

#### Adjusting Initial Data Parameters for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lindquist"

# Adjust the initial data parameters
adjusted_parameters = adjust_parameters(selection)

# Generate the adjusted initial data
initial_data = generate_initial_data(selection, adjusted_parameters)

# Output**Adjusting Grid and Parameter Settings**
=====================================

### Overview of the Notebook

This notebook covers the process of adjusting grid and parameter settings for numerical relativity simulations.

### Theory Review

#### Introduction to Adjusting Grid and Parameter Settings

Adjusting grid and parameter settings is an important step in numerical relativity, as it allows users to customize their simulations and explore different physical scenarios. The adjusted parameters should be consistent with the selected initial data and physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{adjusted initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Adjusting Grid and Parameter Settings in NRPy+

```python
import nrpy_modules as nrm

# Define the adjustable parameters for the selected grid and parameter settings
def adjust_grid_and_parameters():
    # Get the user's input (default values will be used if not provided)
    DestGridCoordSystem = input("Enter destination grid coordinate system (Spherical/Curvilinear): ") or "Curvilinear"
    freeparams = float(input("Enter free parameters (e.g. gamma, beta): ") or 0.5)
    EnableMomentum = input("Enable momentum evolution? (y/n): ") or "n"

    # Validate the user's input
    if DestGridCoordSystem not in ["Spherical", "Curvilinear"]:
        print("Invalid destination grid coordinate system.")
        exit()
    elif freeparams < 0:
        print("Free parameters must be non-negative.")
        exit()

    # Adjust the grid and parameter settings based on the user's input
    adjusted_settings = {
        "DestGridCoordSystem": DestGridCoordSystem,
        "freeparams": freeparams,
        "EnableMomentum": EnableMomentum == "y"
    }

    return adjusted_settings

# Get the selected grid and parameter settings
adjusted_settings = adjust_grid_and_parameters()

# Output the adjusted settings to a file
with open("adjusted_settings.dat", "w") as f:
    for key, value in adjusted_settings.items():
        f.write(f"{key} = {value}\n")
```

### Example Use Case

#### Adjusting Grid and Parameter Settings for a Black Hole Simulation

```python
import numpy as np

# Define the mass**Adjusting Grid Coordinate System**
=====================================

### Overview of the Notebook

This notebook covers the process of adjusting the grid coordinate system for numerical relativity simulations.

### Theory Review

#### Introduction to Adjusting Grid Coordinate System

The grid coordinate system is an important parameter in numerical relativity, as it affects the accuracy and stability of the simulation. In this notebook, we will focus on adjusting the grid coordinate system to either Spherical or SinhSpherical coordinates.

#### Mathematics

$$\mathbf{x} = (x,y,z) \in D$$

where $D$ is the domain of the grid, and $\mathbf{x}$ represents the spatial coordinates.

### Code Implementation

#### Adjusting Grid Coordinate System in NRPy+

```python
import nrpy_modules as nrm

# Define the adjustable parameters for the selected grid coordinate system
def adjust_grid_coordinate_system():
    # Get the user's input (default values will be used if not provided)
    DestGridCoordSystem = input("Enter destination grid coordinate system (Spherical/SinhSpherical): ") or "Curvilinear"

    # Validate the user's input
    if DestGridCoordSystem not in ["Spherical", "SinhSpherical"]:
        print("Invalid destination grid coordinate system.")
        exit()

    # Adjust the grid coordinate system based on the user's input
    adjusted_settings = {
        "DestGridCoordSystem": DestGridCoordSystem
    }

    return adjusted_settings

# Get the selected grid coordinate system
adjusted_settings = adjust_grid_coordinate_system()

# Output the adjusted settings to a file
with open("adjusted_settings.dat", "w") as f:
    for key, value in adjusted_settings.items():
        f.write(f"{key} = {value}\n")
```

### Example Use Case

#### Adjusting Grid Coordinate System for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lindquist"

# Adjust the grid coordinate system
adjusted_settings = adjust_grid_coordinate_system()

# Generate the adjusted initial data
initial_data = generate_initial_data(selection, adjusted_settings)

# Output the adjusted initial data to a file
with open("adjusted_initial_data.dat", "w") as f:
    for key, value in initial_data.items**Customizing NRPy+ for Your Simulation**
=====================================

### Overview of the Notebook

This notebook covers the process of customizing NRPy+ for your simulation. NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity.

### Theory Review

#### Introduction to Numerical Relativity

Numerical relativity is a subfield of general relativity that involves solving Einstein's field equations using numerical methods. This allows us to simulate the behavior of black holes, neutron stars, and other astrophysical systems.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Customizing NRPy+ in Python

```python
import nrpy_modules as nrm

# Define a function to customize NRPy+
def customize_nrpy():
    # Get the user's input (default values will be used if not provided)
    customizations = {
        "grid": input("Enter grid type (Spherical/Curvilinear): ") or "Curvilinear",
        "equations_of_motion": input("Enter equations of motion (BSSN/1+log): ") or "BSSN"
    }

    # Validate the user's input
    if customizations["grid"] not in ["Spherical", "Curvilinear"]:
        print("Invalid grid type.")
        exit()
    elif customizations["equations_of_motion"] not in ["BSSN", "1+log"]:
        print("Invalid equations of motion.")
        exit()

    # Make the necessary modifications to NRPy+
    nrm.modify_nrpy(customizations)

# Call the function to customize NRPy+
customize_nrpy()
```

### Example Use Case

#### Customizing NRPy+ for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lind**Customizing NRPy+ for Different Coordinate Systems**
=====================================================

### Overview of the Notebook

This notebook covers the process of customizing NRPy+ for different coordinate systems. NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity.

### Theory Review

#### Introduction to Numerical Relativity and Coordinate Systems

Numerical relativity is a subfield of general relativity that involves solving Einstein's field equations using numerical methods. This allows us to simulate the behavior of black holes, neutron stars, and other astrophysical systems.

In numerical relativity, the choice of coordinate system is crucial for accurately simulating the behavior of spacetime. Different coordinate systems are suitable for different types of simulations.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Customizing NRPy+ in Python

```python
import nrpy_modules as nrm

# Define a function to customize NRPy+
def customize_nrpy():
    # Get the user's input (default values will be used if not provided)
    CoordSystem = input("Enter coordinate system (Spherical/Curvilinear): ") or "Curvilinear"

    # Validate the user's input
    if CoordSystem not in ["Spherical", "Curvilinear"]:
        print("Invalid coordinate system.")
        exit()

    # Make the necessary modifications to NRPy+
    nrm.modify_nrpy(CoordSystem)

# Call the function to customize NRPy+
customize_nrpy()
```

### Example Use Case

#### Customizing NRPy+ for a Black Hole Simulation in Spherical Coordinates

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lindquist"

# Customize NRPy+ for Spherical coordinates
customize_nrpy()

# Generate the initial data
initial_data = generate_initial_data(selection**Choosing Initial Data for Numerical Relativity Simulations**
===========================================================

### Overview of the Notebook

This notebook covers the process of choosing initial data for numerical relativity simulations.

### Theory Review

#### Introduction to Initial Data Selection

Initial data selection is an important step in numerical relativity, as it provides the necessary information to initialize the numerical evolution. The selected initial data should be accurate and consistent with the physical system being simulated.

#### Mathematics

$$\mathbf{u}(x,y,z,t=0) = \text{selected initial data}$$

where $\mathbf{u}$ represents the set of physical quantities to be evolved, such as the metric components or matter fields.

### Code Implementation

#### Choosing Initial Data in NRPy+

```python
import nrpy_modules as nrm

# Define the initial data string
initial_data_string = "Shifted KerrSchild"

# Print the selected initial data string
print("Selected Initial Data:", initial_data_string)
```

### Example Use Case

#### Choosing Initial Data for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Choose the initial data option (use the default if not provided)
selection = "Brill-Lindquist"

# Generate the initial data
initial_data = generate_initial_data(selection)

# Print the generated initial data
print("Generated Initial Data:", initial_data)
```

Note that this code snippet assumes that the `generate_initial_data` function is defined and implemented elsewhere in the notebook.

### Theory Review

#### Shifted KerrSchild Initial Data

The Shifted KerrSchild initial data is a type of initial data used for numerical relativity simulations. It is based on the Kerr metric, but with a shift term added to the coordinate system.

$$ds^2 = -\left(1-\frac{2GM}{r}\right)dt^2 + \frac{dr^2}{1-\frac{2GM}{r}} + r^2d\theta^2 + r^2\sin^2\theta d\phi^2$$

where $G$ is the gravitational constant, $M$ is the mass of the black hole, and $r$ is the radial coordinate.

### Code Implementation

#### Implementing Shifted KerrSchild Initial Data in NRPy+

```python
import nrpy_modules as**Initializing NRPy+ with Initial Data**
=====================================

### Overview of the Notebook

This notebook covers the process of initializing NRPy+ with initial data.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Initializing NRPy+ with Initial Data in Python

```python
import nrpy_modules as nrm
import collections

# Define a dictionary to store initial data information
dictID = {}

# Add initial data options to the dictionary
dictID['Shifted KerrSchild']  = IDmod_retfunc(
    modulename = "BSSN.ShiftedKerrSchild", functionname = "ShiftedKerrSchild",
    OrigCoordSystem = "Spherical", DestGridCoordSystem = "Spherical",
    freeparams = ["params.M   = 1.0;", "params.a   = 0.9;", "params.r0 = 1.0;"],
    EnableMomentum = True)

dictID['Static Trumpet'] = IDmod_retfunc(
    modulename = "BSSN.StaticTrumpet", functionname = "StaticTrumpet",
    OrigCoordSystem = "Spherical", DestGridCoordSystem = "Spherical",
    freeparams = ["params.M = 1.0;"],
    EnableMomentum = False)

dictID['Brill-Lindquist'] = IDmod_retfunc(
    modulename = "BSSN.BrillLindquist", functionname = "BrillLindquist",
    OrigCoordSystem = "Cartesian", DestGridCoordSystem = "SinhSpherical",
    freeparams = ["params.BH1_posn_x =+1.0; params.BH1_posn_y = 0.0; params.BH1_posn_z = 0**Setting Up NRPy+ Infrastructure**
=====================================

### Overview of the Notebook

This notebook covers the process of setting up the necessary NRPy+ infrastructure and declaring core grid functions.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Setting Up NRPy+ Infrastructure in Python

```python
import nrpy_modules as nrm

# Declare core grid functions
gridfunctions = [
    "metric",
    "christoffel",
    "Riemann",
    "Weyl"
]

# Set up the needed NRPy+ infrastructure
nrm.setup_nrpy_infrastructure(gridfunctions)
```

### Example Use Case

#### Setting Up NRPy+ Infrastructure for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Declare core grid functions
gridfunctions = [
    "metric",
    "christoffel",
    "Riemann",
    "Weyl"
]

# Set up the needed NRPy+ infrastructure
nrm.setup_nrpy_infrastructure(gridfunctions)

# Generate the initial data
initial_data = generate_initial_data("UIUCBlackHole")

# Output the generated initial data to a file
with open("initial_data.dat", "w") as f:
    for key, value in initial_data.items():
        f.write(f"{key} = {value}\n")
```

### Theory Review

#### Declaring Core Grid Functions

The core grid functions are the basic building blocks of any NRPy+ simulation. They are used to compute various quantities, such as the metric tensor, Christoffel symbols, Riemann tensor, and Weyl tensor.

```python
**Initializing NRPy+**
======================

### Overview of the Notebook

This notebook covers the process of initializing NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Importing Core NRPy+ Modules

```python
import nrpy_modules as nrm

# Import the core modules of NRPy that we will need
nrm.import_core_nrpy_modules()
```

#### Specifying Main Grid Functions

```python
# Specify the main grid functions we will need
main_gridfunctions = [
    "metric",
    "christoffel",
    "Riemann",
    "Weyl"
]
```

### Example Use Case

#### Initializing NRPy+ for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import the core modules of NRPy that we will need
nrm.import_core_nrpy_modules()

# Specify the main grid functions we will need
main_gridfunctions = [
    "metric",
    "christoffel",
    "Riemann",
    "Weyl"
]

# Generate the initial data
initial_data = generate_initial_data("UIUCBlackHole")

# Output the generated initial data to a file
with open("initial_data.dat", "w") as f:
    for key, value in initial_data.items():
        f.write(f"{key} = {value}\n")
```

### Theory Review

#### Importing Core NRPy+ Modules

The core modules of NRPy are imported using the `import_core_nrpy_modules()` function. This includes modules such as `nrpy_modules`, `gridfunctions`, and `numerical**Importing NRPy+ Core Modules**
====================================

### Overview of the Notebook

This notebook covers the process of importing necessary NRPy+ core modules.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Importing NRPy+ Core Modules in Python

```python
from outputC import add_to_Cfunction_dict
```

This line of code imports the `add_to_Cfunction_dict` function from the `outputC` module, which is part of the NRPy+ core modules. This function is used to add C functions to a dictionary for use in numerical simulations.

### Example Use Case

#### Importing NRPy+ Core Modules for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import necessary NRPy+ core modules
from outputC import add_to_Cfunction_dict

# Add C functions to dictionary
add_to_Cfunction_dict()
```

### Theory Review

#### Understanding NRPy+ Core Modules

NRPy+ core modules provide the foundation for numerical relativity simulations in Python. The `outputC` module, from which we imported the `add_to_Cfunction_dict` function, is responsible for outputting C functions to a dictionary.

```python
class outputC:
    def add_to_Cfunction_dict(self):
        # Add C functions to dictionary
        pass
```

This class provides a way to add C functions to a dictionary, which can then be used in numerical simulations. The `add_to_Cfunction_dict` method is where the actual addition of C functions takes place.

### Further Implementation

#### Implementing NRPy+ Core Modules

To implement NRPy+ core**NRPy+ Core C Code Output Module**
=====================================

### Overview of the Notebook

This notebook covers the process of importing the core C code output module of NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Mathematics

$$R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

where $R_{\mu\nu}$ is the Ricci tensor, $R$ is the Ricci scalar, $g_{\mu\nu}$ is the metric tensor, and $T_{\mu\nu}$ is the stress-energy tensor.

### Code Implementation

#### Importing NRPy+ Core C Code Output Module in Python

```python
import finite_difference as fin
```

This line of code imports the `finite_difference` module from the `fin` package, which provides a way to perform finite difference calculations. This module is part of the core C code output module of NRPy+, and is used to calculate derivatives and other quantities necessary for numerical simulations.

### Example Use Case

#### Importing NRPy+ Core C Code Output Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import finite difference module
from finite_difference import fin

# Perform finite difference calculations
fin.calculate_derivatives()
```

### Theory Review

#### Understanding Finite Difference Calculations

Finite difference calculations are used to approximate derivatives and other quantities necessary for numerical simulations. The `finite_difference` module provides a way to perform these calculations, and is an essential part of the core C code output module of NRPy+.

```python
class finite_difference:
    def calculate_derivatives(self):
        # Calculate derivatives using finite difference formula
        pass
```

This class provides a way to calculate derivatives using finite difference formulas. The `calculate_derivatives` method is where the actual calculation takes place.

### Further Implementation

#### Implementing Finite Difference Calculations in NRPy+

To implement finite difference calculations in NRPy+, we need to create a function that**NRPy+ Finite Difference C Code Generation Module**
=====================================================

### Overview of the Notebook

This notebook covers the process of importing the finite difference C code generation module of NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Finite Difference Calculations

Finite difference calculations are used to approximate derivatives and other quantities necessary for numerical simulations. The finite difference C code generation module of NRPy+ provides a way to generate C code that performs these calculations.

#### Mathematics

$$\frac{\partial u}{\partial x} \approx \frac{u(x+h) - u(x)}{h}$$

where $u$ is the function being approximated, and $x$ is the point at which the derivative is being calculated.

### Code Implementation

#### Importing NRPy+ Finite Difference C Code Generation Module in Python

```python
import NRPy_param_funcs as par
```

This line of code imports the `NRPy_param_funcs` module, which provides a way to generate parameters for finite difference calculations. This module is part of the finite difference C code generation module of NRPy+, and is used to generate C code that performs finite difference calculations.

### Example Use Case

#### Importing NRPy+ Finite Difference C Code Generation Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import NRPy_param_funcs module
from NRPy_param_funcs import par

# Generate parameters for finite difference calculations
par.generate_parameters()
```

### Theory Review

#### Understanding Finite Difference C Code Generation Module

The finite difference C code generation module of NRPy+ generates C code that performs finite difference calculations. This module uses the `NRPy_param_funcs` module to generate parameters for the calculations.

```python
class NRPy_param_funcs:
    def generate_parameters(self):
        # Generate parameters for finite difference calculations
        pass
```

This class provides a way to generate parameters for finite difference calculations using the `generate_parameters` method. The actual generation of parameters takes place in this method.

### Further Implementation

#### Implementing Finite Difference C Code Generation Module in NRPy+

To implement finite difference C code**NRPy+ Parameter Interface**
=============================

### Overview of the Notebook

This notebook covers the process of importing the parameter interface module of NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Mathematics

$$u(x) = f(x) \cdot v(x)$$

where $u$ is the solution function, $f$ and $v$ are input functions, and $x$ is the independent variable.

### Code Implementation

#### Importing NRPy+ Parameter Interface Module in Python

```python
import grid as gri                # Import grid module for parameter interface
```

This line of code imports the `grid` module from the `gri` package, which provides a way to interface with parameters used in numerical simulations.

### Example Use Case

#### Importing NRPy+ Parameter Interface Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import grid module for parameter interface
from grid import gri

# Set parameters for simulation
gri.set_parameters()
```

### Theory Review

#### Understanding NRPy+ Parameter Interface Module

The parameter interface module of NRPy+ provides a way to set and retrieve parameters used in numerical simulations. This module uses the `grid` module to store and manage these parameters.

```python
class grid:
    def set_parameters(self):
        # Set parameters for simulation
        pass
    
    def get_parameters(self):
        # Get parameters for simulation
        pass
```

This class provides a way to set and retrieve parameters using the `set_parameters` and `get_parameters` methods. The actual setting and retrieval of parameters takes place in these methods.

### Further Implementation

#### Implementing NRPy+ Parameter Interface Module

To implement the parameter interface module in NRPy+, we need to create functions that can set and retrieve parameters.

```python
def set_parameter(self, name, value):
    # Set a single parameter
    pass
    
def get_parameter(self, name):
    # Get a single parameter
    pass
```

These functions will be used to implement the `set_parameters` and `get_parameters` methods in the `grid` class**NRPy+ Numerical Grid Functions**
=====================================

### Overview of the Notebook

This notebook covers the process of importing the functions related to numerical grids in NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Numerical Grids

Numerical grids are used to discretize the continuous spacetime into a set of discrete points, which are then used to solve the Einstein field equations numerically.

#### Mathematics

$$\Gamma^{i}_{jk} = \frac{1}{2} g^{im} (\partial_j g_{mk} + \partial_k g_{mj})$$

where $\Gamma^{i}_{jk}$ is the Christoffel symbol, $g_{ij}$ is the metric tensor, and $\partial_i$ denotes partial differentiation with respect to the $i$-th coordinate.

### Code Implementation

#### Importing NRPy+ Numerical Grid Functions in Python

```python
import reference_metric as rfm    # Import reference metric module for numerical grids
```

This line of code imports the `reference_metric` module from the `rfm` package, which provides a way to define and manipulate numerical grids.

### Example Use Case

#### Importing NRPy+ Numerical Grid Functions for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import reference metric module for numerical grids
from reference_metric import rfm

# Define the numerical grid parameters
rfm.set_numerical_grid_parameters()
```

### Theory Review

#### Understanding NRPy+ Numerical Grid Functions

The numerical grid functions in NRPy+ provide a way to define and manipulate numerical grids, which are used to discretize the continuous spacetime into a set of discrete points.

```python
class reference_metric:
    def set_numerical_grid_parameters(self):
        # Set numerical grid parameters
        pass
    
    def get_numerical_grid_parameters(self):
        # Get numerical grid parameters
        pass
```

This class provides a way to set and retrieve numerical grid parameters using the `set_numerical_grid_parameters` and `get_numerical_grid_parameters` methods.

### Further Implementation

#### Implementing NRPy**NRPy+ Reference Metric Support**
=====================================

### Overview of the Notebook

This notebook covers the process of importing the reference metric support module in NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Reference Metric

The reference metric is the background metric used in numerical relativity simulations. It is used to define the spacetime geometry and to calculate various quantities such as the Christoffel symbols and the Riemann tensor.

#### Mathematics

$$g_{\mu\nu} = \frac{\partial x^{\alpha}}{\partial X^{\mu}} \frac{\partial x^{\beta}}{\partial X^{\nu}} g_{\alpha\beta}$$

where $g_{\mu\nu}$ is the background metric, $x^{\alpha}$ are the coordinates of the spacetime point, and $X^{\mu}$ are the coordinates of the reference point.

### Code Implementation

#### Importing NRPy+ Reference Metric Support Module in Python

```python
import cmdline_helper as cmd      # Import cmdline helper module for reference metric support
```

This line of code imports the `cmdline_helper` module from the `cmd` package, which provides a way to interact with the command line and to define the reference metric.

### Example Use Case

#### Importing NRPy+ Reference Metric Support Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import cmdline helper module for reference metric support
from cmdline_helper import cmd

# Define the reference metric parameters
cmd.define_reference_metric()
```

### Theory Review

#### Understanding NRPy+ Reference Metric Support

The reference metric support in NRPy+ provides a way to define and manipulate the background metric used in numerical relativity simulations.

```python
class cmdline_helper:
    def define_reference_metric(self):
        # Define reference metric parameters
        pass
    
    def get_reference_metric_parameters(self):
        # Get reference metric parameters
        pass
```

This class provides a way to set and retrieve reference metric parameters using the `define_reference_metric` and `get_reference_metric_parameters` methods.

### Further Implementation

**NRPy+ Multi-Platform Python Command-Line Interface**
=====================================================

### Overview of the Notebook

This notebook covers the process of importing the multi-platform Python command-line interface module in NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Command-Line Interface

The command-line interface (CLI) is a way of interacting with the NRPy+ code through the terminal or command prompt. This allows users to run simulations and analyze results without needing to modify the underlying code.

#### Mathematics

$$\mathcal{M} = \left( \begin{array}{ccc} g_{00} & g_{01} & g_{02} \\ g_{10} & g_{11} & g_{12} \\ g_{20} & g_{21} & g_{22} \end{array} \right)$$

where $\mathcal{M}$ is the metric tensor, and $g_{\mu\nu}$ are the components of the metric.

### Code Implementation

#### Importing NRPy+ Multi-Platform Python Command-Line Interface Module in Python

```python
import shutil, os, time      # Import necessary modules for multi-platform CLI
```

This line of code imports the `shutil`, `os`, and `time` modules, which are used to interact with the file system, manipulate files and directories, and handle timing and synchronization.

### Example Use Case

#### Importing NRPy+ Multi-Platform Python Command-Line Interface Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import necessary modules for multi-platform CLI
from shutil import copy2
from os import system
from time import sleep

# Run simulation using NRPy+
system("python run_simulation.py")
```

### Theory Review

#### Understanding NRPy+ Multi-Platform Python Command-Line Interface

The multi-platform Python command-line interface in NRPy+ provides a way to interact with the code through the terminal or command prompt. This allows users to run simulations and analyze results without needing to modify the underlying code.

```python
class shutil:
    def copy2(self, src, dst):
       **Standard Python Modules for Multi-Platform OS-Level Functions and Benchmarking**
=====================================================================================

### Overview of the Notebook

This notebook covers the process of importing standard Python modules for multi-platform OS-level functions and benchmarking.

### Theory Review

#### Introduction to Multi-Platform OS-Level Functions

Multi-platform OS-level functions are used to interact with the operating system, manipulate files and directories, and perform other low-level tasks. These functions are essential for creating cross-platform software that can run on multiple operating systems.

#### Benchmarking

Benchmarking is the process of measuring the performance of a program or algorithm by running it under controlled conditions and comparing its results against a baseline. Benchmarking is an important tool for optimizing code, identifying bottlenecks, and ensuring that programs meet their performance requirements.

#### Mathematics

$$T_{avg} = \frac{1}{N} \sum_{i=1}^{N} t_i$$

where $T_{avg}$ is the average time taken to execute a program or algorithm, and $t_i$ are the individual times taken for each execution.

### Code Implementation

#### Importing Standard Python Modules in Python

```python
import diagnostics_generic.process_2D_data as plot2D   # Import plotting module for 2D data
```

This line of code imports the `process_2D_data` module from the `diagnostics_generic` package, which is used to plot and analyze 2D data.

### Example Use Case

#### Using Standard Python Modules for Multi-Platform OS-Level Functions and Benchmarking

```python
import numpy as np
from time import perf_counter

# Define a function to be benchmarked
def my_function():
    x = np.arange(1000000)
    y = x**2
    return x, y

# Benchmark the function
start_time = perf_counter()
my_function()
end_time = perf_counter()

print("Time taken:", end_time - start_time)

# Plot 2D data using plot2D module
import matplotlib.pyplot as plt
x, y = my_function()
plt.plot(x, y)
plt.show()
```

### Theory Review

#### Understanding Standard Python Modules for Multi-Platform OS-Level Functions and Benchmarking

The standard Python modules for multi-platform OS-level functions and benchmarking provide a range of tools and functions for interacting with the operating system, manipulating files and directories, and measuring performance.

```python
class diagnostics_generic:
    def process_2**NRPy+: Analysis of Output Data**
=====================================

### Overview of the Notebook

This notebook covers the process of importing the module for analyzing output data in NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Output Data Analysis

Output data analysis is an essential part of any numerical simulation. In NRPy+, output data is analyzed using various diagnostic tools and plotting functions.

#### Mathematics

$$\delta g_{\mu\nu} = g_{\mu\nu} - \bar{g}_{\mu\nu}$$

where $\delta g_{\mu\nu}$ is the deviation of the metric tensor from its background value, $g_{\mu\nu}$ is the actual metric tensor, and $\bar{g}_{\mu\nu}$ is the background metric tensor.

### Code Implementation

#### Importing NRPy+ Output Data Analysis Module in Python

```python
import diagnostics_generic.output_yz_or_xy_plane as planar_diags   # Import plotting module for y-z or x-y plane
```

This line of code imports the `output_yz_or_xy_plane` module from the `diagnostics_generic` package, which is used to plot output data on the y-z or x-y plane.

### Example Use Case

#### Using NRPy+ Output Data Analysis Module for a Black Hole Simulation

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import plotting module for y-z or x-y plane
from diagnostics_generic import output_yz_or_xy_plane as planar_diags

# Plot output data on y-z or x-y plane
planar_diags.plot_planar_data()
```

### Theory Review

#### Understanding NRPy+ Output Data Analysis Module

The output data analysis module in NRPy+ provides a range of tools and functions for analyzing output data. These include plotting functions, diagnostic tools, and data manipulation functions.

```python
class diagnostics_generic:
    def plot_planar_data(self):
        # Plot planar data on y-z or x-y plane
        pass
    
    def calculate_diagnostic_quantities(self):
        # Calculate diagnostic quantities from output data
        pass**NRPy+: C Code Generation for Output Data**
=============================================

### Overview of the Notebook

This notebook covers the process of generating C code for output data in NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Output Data Generation

Output data generation is an essential part of any numerical simulation. In NRPy+, output data is generated using C code, which is then compiled and executed on the target machine.

#### Mathematics

$$\delta g_{\mu\nu} = \frac{\partial g_{\mu\nu}}{\partial t} + \mathcal{L}_{\xi}g_{\mu\nu}$$

where $\delta g_{\mu\nu}$ is the deviation of the metric tensor from its background value, $g_{\mu\nu}$ is the actual metric tensor, and $\mathcal{L}_{\xi}g_{\mu\nu}$ is the Lie derivative of the metric tensor.

### Code Implementation

#### Importing NRPy+ C Code Generation Module in Python

```python
import outputC as outc   # Import C code generation module for output data
```

This line of code imports the `outputC` module, which is used to generate C code for output data.

### Example Use Case

#### Generating C Code for Output Data using NRPy+

```python
import numpy as np

# Define the mass and radius of the black hole
M = 1.0
r = 10.0

# Import C code generation module for output data
from outputC import outc

# Generate C code for output data
outc.generate_output_data()
```

### Theory Review

#### Understanding NRPy+ C Code Generation Module

The C code generation module in NRPy+ provides a range of tools and functions for generating C code. These include functions for generating output data, calculating derivatives, and manipulating C code.

```python
class outputC:
    def generate_output_data(self):
        # Generate C code for output data
        pass
    
    def calculate_derivatives(self):
        # Calculate derivatives using C code
        pass
```

This class provides a way to generate C code and manipulate it using the `generate_output_data` and `**Step P2: Create C Code Output Directory**
=============================================

### Overview of the Notebook

This notebook covers the process of creating a directory for storing C code generated by NRPy+.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a Python module that provides a framework for solving Einstein's field equations in numerical relativity. It uses a modular approach, where each module can be easily integrated with others to create complex simulations.

#### Creating Output Directory

Creating an output directory is an essential step in the process of generating C code using NRPy+. This directory will store the generated C code and other related files.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import os   # Import os module for file system operations
```

This line of code imports the `os` module, which is used to interact with the file system.

#### Creating Output Directory using NRPy+

```python
Ccodesrootdir = os.path.join("Exact_InitialDataPlayground_Ccodes")   # Define output directory path

# Check if directory exists, create it otherwise
if not os.path.exists(Ccodesrootdir):
    os.makedirs(Ccodesrootdir)
```

This code checks if the `Ccodesrootdir` directory exists. If it does not exist, the `os.makedirs()` function is used to create it.

### Theory Review

#### Understanding NRPy+ Output Directory Creation

The output directory creation process in NRPy+ involves defining a path for the directory and checking if it exists. If the directory does not exist, it is created using the `os.makedirs()` function.

```python
class os:
    def path(self, *args):
        # Return a string representing a path
        
    def exists(self, path):
        # Check if a file or directory exists at a given path
        
    def makedirs(self, path):
        # Create a new directory and all its parents if they do not exist
```

This class provides a way to interact with the file system using functions such as `path()`, `exists()` and `makedirs()`.**Step P2.1: Remove C Code Output Directory if it Exists**
===========================================================

### Overview of the Notebook

This notebook covers the process of removing a directory for storing C code generated by NRPy+ if it already exists.

### Theory Review

#### Introduction to Directory Removal

Removing an existing directory is a common operation in file system management. In this case, we want to remove the C code output directory if it already exists before generating new C code.

#### Mathematics

$$\text{remove directory} \iff \text{directory does not exist or is empty}$$

where $\text{remove directory}$ means removing a directory from the file system and $\text{directory does not exist or is empty}$ means that the directory either does not exist at all or it is empty.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import shutil, os   # Import shutil and os modules for file system operations
```

This line of code imports the `shutil` and `os` modules, which are used to interact with the file system.

#### Removing C Code Output Directory using NRPy+

```python
Ccodesrootdir = os.path.join("Exact_InitialDataPlayground_Ccodes")   # Define output directory path

# Check if directory exists and remove it
if os.path.exists(Ccodesrootdir):
    shutil.rmtree(Ccodesrootdir)
```

This code checks if the `Ccodesrootdir` directory exists. If it does exist, the `shutil.rmtree()` function is used to remove it.

### Theory Review

#### Understanding NRPy+ Directory Removal

The directory removal process in NRPy+ involves checking if a directory exists and removing it using the `shutil.rmtree()` function if it does.

```python
class os:
    def path(self, *args):
        # Return a string representing a path
        
    def exists(self, path):
        # Check if a file or directory exists at a given path
        
    def makedirs(self, path):
        # Create a new directory and all its parents if they do not exist

class shutil:
    def rmtree(self, path):
        # Remove the entire directory tree rooted at path
```

This code uses the `os` and `shutil` classes to interact with the file system.**Removing a Non-Empty Directory in Python**
=============================================

### Overview of the Notebook

This notebook covers the process of removing a non-empty directory in Python using the `shutil.rmtree()` function.

### Theory Review

#### Introduction to Directory Removal

Removing a directory is a common operation in file system management. However, when the directory is not empty, it can be more complicated than simply deleting all its contents.

#### Mathematics

$$\text{remove directory} \iff \text{directory does not exist or is empty}$$

where $\text{remove directory}$ means removing a directory from the file system and $\text{directory does not exist or is empty}$ means that the directory either does not exist at all or it is empty.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import shutil   # Import shutil module for file system operations
```

This line of code imports the `shutil` module, which is used to interact with the file system.

#### Removing Non-Empty Directory using NRPy+

```python
Ccodesrootdir = os.path.join("Exact_InitialDataPlayground_Ccodes")   # Define output directory path

# Check if directory exists and remove it
if os.path.exists(Ccodesrootdir):
    shutil.rmtree(Ccodesrootdir, ignore_errors=True)
```

This code checks if the `Ccodesrootdir` directory exists. If it does exist, the `shutil.rmtree()` function is used to remove it, ignoring any errors that may occur.

### Theory Review

#### Understanding NRPy+ Directory Removal with Errors Ignored

The `shutil.rmtree()` function in NRPy+ provides a way to remove non-empty directories by ignoring any errors that may occur during removal. This makes it easier to handle cases where the directory contains files or subdirectories that cannot be deleted.

```python
class shutil:
    def rmtree(self, path, ignore_errors=False):
        # Remove the entire directory tree rooted at path
        
        if not ignore_errors:
            raise OSError("Directory not empty")
```

This code uses the `shutil` class to interact with the file system. The `rmtree()` function removes the specified directory and all its contents, ignoring any errors that may occur if `ignore_errors=True`.

### Further Implementation

#### Removing a Directory with Errors Ignored

To remove a directory with errors ignored, you can use the following code:

```python
import os
shutil.rmtree**Step P2.3: Create Fresh Directory for C Code Output**
=============================================================

### Overview of the Notebook

This notebook covers the process of creating a new directory for storing C code generated by NRPy+.

### Theory Review

#### Introduction to Directory Creation

Creating a new directory is an essential step in the process of generating C code using NRPy+. This allows us to store the generated C code and other related files.

#### Mathematics

$$\text{create directory} \iff \text{directory does not exist}$$

where $\text{create directory}$ means creating a new directory in the file system, and $\text{directory does not exist}$ means that the directory is not already present.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import cmdline_helper as cmd   # Import cmdline helper module for directory operations
```

This line of code imports the `cmdline_helper` module from the `cmd` package, which provides a way to interact with the command line and perform file system operations.

#### Creating Fresh Directory using NRPy+

```python
Ccodesrootdir = os.path.join("Exact_InitialDataPlayground_Ccodes")   # Define output directory path

# Create fresh directory for C code output
cmd.mkdir(Ccodesrootdir)
```

This code defines the `Ccodesrootdir` variable and uses the `mkdir()` function from the `cmdline_helper` module to create a new directory at the specified location.

### Theory Review

#### Understanding NRPy+ Directory Creation

The `cmdline_helper` module in NRPy+ provides a way to interact with the command line and perform file system operations, including creating directories. The `mkdir()` function is used to create a new directory at the specified location.

```python
class cmdline_helper:
    def mkdir(self, path):
        # Create a new directory at the specified location
        
        if not os.path.exists(path):
            os.makedirs(path)
```

This code uses the `cmdline_helper` class to interact with the command line and perform file system operations. The `mkdir()` function creates a new directory at the specified location.

### Further Implementation

#### Creating a Directory with Existing Subdirectories

To create a directory with existing subdirectories, you can use the following code:

```python
import os

os.makedirs(Ccodesrootdir, exist_ok=True)
```

This code uses the `makedirs()` function from the `os` module to create**Step P3: Create Executable Output Directory**
=============================================

### Overview of the Notebook

This notebook covers the process of creating an executable output directory for storing compiled C code generated by NRPy+.

### Theory Review

#### Introduction to Output Directory Creation

Creating an executable output directory is an essential step in the process of generating C code using NRPy+. This allows us to store the compiled C code and other related files that can be executed on a target machine.

#### Mathematics

$$\text{executable output directory} = \text{output directory} + \text{executables}$$

where $\text{executable output directory}$ is the final output directory, $\text{output directory}$ is the base directory for storing compiled C code and related files, and $\text{executables}$ are the compiled C code files that can be executed on a target machine.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import cmdline_helper as cmd   # Import cmdline helper module for directory operations
```

This line of code imports the `cmdline_helper` module from the `cmd` package, which provides a way to interact with the command line and perform file system operations.

#### Creating Executable Output Directory using NRPy+

```python
Ccodesrootdir = os.path.join("Exact_InitialDataPlayground_Ccodes")   # Define output directory path

# Create executable output directory for compiled C code
outdir = os.path.join(Ccodesrootdir, "output")
cmd.mkdir(outdir)
```

This code defines the `Ccodesrootdir` variable and uses the `mkdir()` function from the `cmdline_helper` module to create a new directory at the specified location.

### Theory Review

#### Understanding NRPy+ Output Directory Creation

The output directory creation process in NRPy+ involves creating a base directory for storing compiled C code and related files, and then creating an executable output directory within that base directory. The `mkdir()` function is used to create the new directories at the specified locations.

```python
class cmdline_helper:
    def mkdir(self, path):
        # Create a new directory at the specified location
        
        if not os.path.exists(path):
            os.makedirs(path)
```

This code uses the `cmdline_helper` class to interact with the command line and perform file system operations. The `mkdir()` function creates a new directory at the specified location.

### Further Implementation

#### Creating an Executable Output**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. This allows us to use finite-difference stencils for numerical differentiation and integration.

### Theory Review

#### Introduction to Finite-Difference Methods

Finite-difference methods are a type of numerical method used to approximate derivatives and integrals. They work by approximating the derivative or integral using the values of the function at nearby points, rather than trying to compute it exactly.

#### Mathematics

$$\frac{df}{dx} \approx \frac{\Delta x}{f(x + \Delta x) - f(x)}$$

where $\frac{df}{dx}$ is the derivative of the function $f$ with respect to $x$, and $\Delta x$ is a small change in $x$.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
```

This code uses the `NRPy_code_generator` class to interact with the finite-difference function implementation. The `enable_FD_functions()` function enables finite-difference functions in NRPy+ and generates C code for numerical differentiation and integration.

### Further Implementation**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to Inlined Static Functions

Inlined static functions are a way to implement small, self-contained functions that can be executed directly without the need for a function call. This can improve performance by reducing overhead associated with function calls.

#### Mathematics

$$f(x) \approx f_{\text{static}}(x) = \begin{cases} f(x + h), & x < 0 \\ f(x - h), & x > 0 \end{cases}$$

where $f(x)$ is the original function, and $f_{\text{static}}(x)$ is the inlined static function approximation.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which are then output as inlined static functions.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Output inlined static functions
        
        print("Inlined static functions enabled.")
```

This code uses the `NRPy_code_generator` class to interact with the finite-difference function implementation. The**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to GCC and Compilation

GCC (GNU Compiler Collection) is a compiler suite that includes front ends for C, C++, and other languages. It is widely used for compiling source code into machine-executable code.

#### Mathematics

$$\text{Compilation} = \text{Source Code} + \text{Compiler} + \text{Flags}$$

where $\text{Compilation}$ is the process of translating source code into executable code, $\text{Source Code}$ is the original code to be compiled, $\text{Compiler}$ is the tool used for compilation (in this case, GCC), and $\text{Flags}$ are options passed to the compiler.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and GCC

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which are then compiled using GCC.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Compile C code using GCC
        
        print("Compiling highly complex FD kernels with GCC.")
```

This code uses the `NRPy_code_generator` class to interact with**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to GCC and Compilation

GCC (GNU Compiler Collection) is a compiler suite that includes front ends for C, C++, and other languages. It is widely used for compiling source code into machine-executable code.

#### Mathematics

$$\text{Compilation} = \text{Source Code} + \text{Compiler} + \text{Flags}$$

where $\text{Compilation}$ is the process of translating source code into executable code, $\text{Source Code}$ is the original code to be compiled, $\text{Compiler}$ is the tool used for compilation (in this case, GCC), and $\text{Flags}$ are options passed to the compiler.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and GCC

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which are then compiled using GCC.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Compile C code using GCC
        
        print("Compiling highly complex FD kernels with GCC.")
```

This code uses the `NRPy_code_generator` class to interact with**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to Compilation Time

Compilation time is a critical aspect of code compilation. It refers to the amount of time taken by the compiler to translate source code into executable code.

#### Mathematics

$$\text{Compilation Time} = \text{Source Code Complexity} + \text{Compiler Efficiency} + \text{Flags}$$

where $\text{Compilation Time}$ is the time taken for compilation, $\text{Source Code Complexity}$ is the complexity of the source code being compiled, $\text{Compiler Efficiency}$ is the efficiency of the compiler in translating the source code into executable code, and $\text{Flags}$ are options passed to the compiler.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and GCC

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which are then compiled using GCC.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Compile C code using GCC
        
        print("Compiling highly complex FD kernels with GCC.")
```

This code uses the `NRPy_code_generator` class to interact with

### GCC**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to Code Performance

Code performance refers to the time taken by a program to execute a given task. In this case, we are concerned with the effect of enabling finite-difference functions on code performance.

#### Mathematics

$$\text{Code Performance} = \frac{\text{Execution Time}}{\text{Complexity of Task}}$$

where $\text{Code Performance}$ is a measure of how efficiently a program executes a task, $\text{Execution Time}$ is the time taken by the program to execute the task, and $\text{Complexity of Task}$ is a measure of the difficulty of the task being performed.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and Code Performance

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which do not slow down code performance but do add additional functionality.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Compile C code using GCC
        
        print("FD functions enabled.")
```

This code uses the `NRPy_code_generator` class to interact with

### FD Functions and**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to Header Files

Header files are used to declare function prototypes, macros, and other definitions that can be used by multiple source files. In this case, we need to add another header file to the C source tree.

#### Mathematics

$$\text{Header File} = \begin{cases} \text{Prototype declaration for functions} & \text{if function is declared in a separate file} \\ \text{No prototype declaration} & \text{otherwise} \end{cases}$$

where $\text{Header File}$ is the header file that contains the declarations and definitions, and $\text{Prototype declaration}$ is the declaration of a function without its definition.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
FD_functions_enabled = True  # Enable finite-difference functions

if FD_functions_enabled:
    NR.enable_FD_functions()   # Enable finite-difference functions in NRPy+
```

This code defines a variable `FD_functions_enabled` to indicate whether finite-difference functions are enabled. If they are, the `enable_FD_functions()` function is called to enable them in NRPy+.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and Header Files

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which add another header file to the C source tree.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        # Add another header file to the C source tree
        
        print("Adding another header file to the C source tree.")
```

This**Step 1.c: Enable Finite-Difference Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of enabling finite-difference functions in NRPy+. These functions are used to perform numerical differentiation and integration.

### Theory Review

#### Introduction to Performance Impact

Enabling finite-difference functions can have a significant impact on code performance. In this case, we will investigate the effect of enabling these functions with different versions of GCC (GNU Compiler Collection).

#### Mathematics

$$\text{Performance Impact} = \frac{\text{Execution Time}}{\text{Original Execution Time}} \times 100$$

where $\text{Performance Impact}$ is a measure of how much the code performance has changed, $\text{Execution Time}$ is the time taken by the program to execute with finite-difference functions enabled, and $\text{Original Execution Time}$ is the original execution time without finite-difference functions.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling Finite-Difference Functions using NRPy+

```python
enable_FD_functions = False  # Disable finite-difference functions

print("Original Execution Time:", original_execution_time)
```

This code defines a variable `enable_FD_functions` to indicate whether finite-difference functions are enabled. If they are, the `original_execution_time` will be used as a reference point for measuring performance impact.

### Theory Review

#### Understanding NRPy+ Finite-Difference Functions and Performance Impact

The finite-difference function implementation in NRPy+ involves using the `NRPy_code_generator` module to generate C code for numerical differentiation and integration. The `enable_FD_functions()` function is used to enable these functions in NRPy+, which can have a significant impact on code performance.

```python
class NRPy_code_generator:
    def enable_FD_functions(self):
        # Enable finite-difference functions in NRPy+
        
        self.FD_functions_enabled = True
        
        # Generate C code for finite-difference stencils
        
        print("Enabling finite-difference functions.")
```

This code uses the `NRPy_code_generator` class to interact with

### Performance Impact of Enabling Finite-Difference Functions

```python
enable_FD**Step 2: Set Core Parameters**
=====================================

### Overview of the Notebook

This notebook covers the process of setting core parameters for NRPy+ simulations. These parameters include the coordinate system, finite-difference (FD) order, and floating-point precision.

### Theory Review

#### Introduction to Coordinate Systems

Coordinate systems are used to describe the geometry of a problem in computational fluid dynamics (CFD). The most common coordinate systems are:

* Cartesian coordinates: `(x, y, z)`
* Cylindrical coordinates: `(r, \theta, z)`
* Spherical coordinates: `(r, \theta, \phi)`

#### Mathematics

$$\mathbf{x} = (x, y, z)$$
$$\mathbf{y} = (y, z, r)$$
$$\mathbf{z} = (z, r, \theta)$$

where $\mathbf{x}$, $\mathbf{y}$, and $\mathbf{z}$ are the coordinates in Cartesian, cylindrical, and spherical coordinate systems respectively.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Core Parameters using NRPy+

```python
# Set coordinate system
coord_system = "Cartesian"  # or "Cylindrical", "Spherical"

# Set FD order
fd_order = 4

# Set floating-point precision
fp_precision = 64  # in bits

print("Core parameters set.")
```

This code defines the core parameters, including the coordinate system, FD order, and floating-point precision.

### Theory Review

#### Understanding NRPy+ Core Parameters

The core parameters in NRPy+ are used to configure the simulation. The `coord_system` parameter determines the geometry of the problem, while the `fd_order` parameter determines the accuracy of the finite-difference scheme. Finally, the `fp_precision` parameter determines the precision of floating-point arithmetic.

```python
class NRPy_code_generator:
    def set_core_params(self):
        # Set core parameters
        
        self.coord_system = coord_system  # or "Cylindrical", "Spherical"
        
        self.fd_order = fd_order
        
        self.fp_precision =**Step 2: Set Core Parameters**
=====================================

### Overview of the Notebook

This notebook covers the process of setting core parameters for NRPy+ simulations. These parameters include the coordinate system, finite-difference (FD) order, and floating-point precision.

### Theory Review

#### Introduction to Coordinate Systems in NRPy+

NRPy+ supports several coordinate systems, including:

* **Spherical**: The most common coordinate system used in astrophysical simulations.
* **SinhSpherical**: A modified spherical coordinate system that can handle singularities at the origin.
* **SinhSphericalv2**: An improved version of SinhSpherical, with better performance and accuracy.
* **Cylindrical**: A coordinate system suitable for axisymmetric problems.
* **SinhCylindrical**: A modified cylindrical coordinate system that can handle singularities at the origin.

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$\mathbf{y} = (\sigma, \eta, \zeta)$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, and $\mathbf{y}$ are the coordinates in SinhSpherical and SinhCylindrical coordinate systems.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Coordinate System using NRPy+

```python
# Set coordinate system
coord_system = "Spherical"  # or "Cylindrical", "SinhSpherical", "SinhSphericalv2", "SinhCylindrical"

print("Coordinate system set to:", coord_system)
```

This code defines the coordinate system, which can be one of the above choices.

### Theory Review

#### Understanding NRPy+ Coordinate Systems

Each coordinate system has its own advantages and disadvantages. For example:

* **Spherical**: Suitable for most astrophysical simulations, but may not handle singularities well.
* **SinhSpherical**: Handles singularities at the origin better than Spherical, but may be slower in performance.
* **Cylindrical**: Suitable for axisymmetric problems, but may not handle non-axisymmetric effects.

```python
class NRPy_code_generator:
**Step 2: Set Core Parameters**
=====================================

### Overview of the Notebook

This notebook covers the process of setting core parameters for NRPy+ simulations. These parameters include the coordinate system, finite-difference (FD) order, and floating-point precision.

### Theory Review

#### Introduction to Coordinate Systems in NRPy+

NRPy+ supports several coordinate systems, including:

* **Spherical**: The most common coordinate system used in astrophysical simulations.
* **Cylindrical**: A coordinate system suitable for axisymmetric problems.
* **SymTP**: A symmetrized temporal parity (STP) coordinate system.
* **SinhSymTP**: A modified STP coordinate system that can handle singularities at the origin.

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$\mathbf{y} = (\sigma, \eta, \zeta)$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, and $\mathbf{y}$ are the coordinates in SymTP and SinhSymTP coordinate systems.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Coordinate System using NRPy+

```python
# Set coordinate system
CoordSystem = "Spherical"  # or "Cylindrical", "SymTP", "SinhSymTP"

# Set reference metric parameters
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)

# Generate reference metric code
rfm.reference_metric()
```

This code defines the coordinate system, which can be one of the above choices. The `par.set_parval_from_str()` function is used to set the parameter value from a string, and the `rfm.reference_metric()` function generates the reference metric code.

### Theory Review

#### Understanding NRPy+ Coordinate Systems

Each coordinate system has its own advantages and disadvantages. For example:

* **Spherical**: Suitable for most astrophysical simulations, but may not handle singularities well.
* **Cylindrical**: Suitable for axisymmetric problems, but may not handle non-axisymmetric effects.
* **SymTP**: Handles temporal parity symmetries better than**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Coordinate System Parameters

The coordinate system parameters in NRPy+ include:

* `CoordSystem`: The name of the coordinate system (e.g. "Spherical", "Cylindrical", etc.)
* `RGrid`: A flag indicating whether a radial grid is used
* `ThetaGrid`: A flag indicating whether a theta grid is used
* `PhiGrid`: A flag indicating whether a phi grid is used

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$R_{\text{min}} = 0.1M$$
$$R_{\text{max}} = 10M$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, $R_{\text{min}}$ is the minimum radial distance from the center of the grid, and $R_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Defaults for Coordinate System Parameters using NRPy+

```python
# Set defaults for coordinate system parameters
CoordSystem = "Spherical"  # or "Cylindrical", etc.
RGrid = True
ThetaGrid = True
PhiGrid = False

par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
par.set_parval_from_str("numerical_grid::RGrid", str(RGrid))
par.set_parval_from_str("numerical_grid::ThetaGrid", str(ThetaGrid))
par.set_parval_from_str("numerical_grid::PhiGrid", str(PhiGrid))

print("Coordinate system parameters set.")
```

This code sets the default values for coordinate system parameters, including `CoordSystem`, `RGrid`, `ThetaGrid`, and `PhiGrid`. The `par.set_parval_from_str()` function is used to set the parameter**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Commonly Adjusted Parameters

The following parameters are perhaps the most commonly adjusted:

* `CoordSystem`: The name of the coordinate system (e.g. "Spherical", "Cylindrical", etc.)
* `RGrid`: A flag indicating whether a radial grid is used
* `ThetaGrid`: A flag indicating whether a theta grid is used
* `PhiGrid`: A flag indicating whether a phi grid is used

These parameters are crucial in determining the geometry of the numerical grid and the reference metric.

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$R_{\text{min}} = 0.1M$$
$$R_{\text{max}} = 10M$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, $R_{\text{min}}$ is the minimum radial distance from the center of the grid, and $R_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Defaults for Coordinate System Parameters using NRPy+

```python
# Set defaults for coordinate system parameters
CoordSystem = "Spherical"  # or "Cylindrical", etc.
RGrid = True
ThetaGrid = True
PhiGrid = False

par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
par.set_parval_from_str("numerical_grid::RGrid", str(RGrid))
par.set_parval_from_str("numerical_grid::ThetaGrid", str(ThetaGrid))
par.set_parval_from_str("numerical_grid::PhiGrid", str(PhiGrid))

print("Coordinate system parameters set.")
```

This code sets the default values for coordinate system parameters, including `CoordSystem`, `RGrid`, `ThetaGrid`, and `**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to High-Level Modifications

By enabling modifications at a high level, we can easily modify various aspects of the simulation without having to delve into the underlying code. This allows us to focus on the physics of the problem rather than the implementation details.

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$R_{\text{min}} = 0.1M$$
$$R_{\text{max}} = 10M$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, $R_{\text{min}}$ is the minimum radial distance from the center of the grid, and $R_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Enabling High-Level Modifications using NRPy+

```python
# Enable high-level modifications at this level
high_level_modifications_enabled = True

par.set_parval_from_str("reference_metric::HighLevelModifications", str(high_level_modifications_enabled))

print("High-level modifications enabled.")
```

This code enables high-level modifications at the desired level, allowing us to modify various aspects of the simulation without having to delve into the underlying code.

### Theory Review

#### Understanding High-Level Modifications in NRPy+

By enabling high-level modifications, we can easily modify various aspects of the simulation, such as:

* Coordinate system parameters
* Numerical grid settings
* Reference metric parameters

This allows us to focus on the physics of the problem rather than the implementation details.

```python
class NRPy_code_generator:
    def enable_high_level_modifications(self):
        # Enable high-level modifications at this level
        
        self.high_level_modifications_enabled = True
        
        print("High-level modifications enabled.")
```

This code uses the `NRPy_code_generator` class to interact with**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Domain Size Parameter

The `domain_size` parameter sets the default value for various domain-related settings, including:

* Grid size
* Boundary conditions
* Initial conditions

By setting a default value for `domain_size`, we can easily modify these parameters without having to delve into the underlying code.

#### Mathematics

$$\mathbf{x} = (r, \theta, \phi)$$
$$R_{\text{min}} = 0.1M$$
$$R_{\text{max}} = 10M$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, $R_{\text{min}}$ is the minimum radial distance from the center of the grid, and $R_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Domain Size using NRPy+

```python
# Set default value for domain size
domain_size = 10.0  # units: M

par.set_parval_from_str("numerical_grid::DomainSize", str(domain_size))

print("Default value for domain size set.")
```

This code sets the default value for `domain_size` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Domain Size Parameter in NRPy+

The `domain_size` parameter determines the size of the numerical grid, which is crucial in setting up the simulation. By setting a default value for `domain_size`, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_domain_size(self):
        # Set default value for domain size
        
        self.domain_size = domain_size
        
        print("Default value for domain size set.")
```

This code uses the `NRPy_code_generator` class**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Spherical's RMAX Parameter

The `RMAX` parameter is a crucial component of the Spherical coordinate system. It determines the maximum radial distance from the center of the grid, which is essential for setting up the simulation.

$$\mathbf{x} = (r, \theta, \phi)$$
$$R_{\text{min}} = 0.1M$$
$$R_{\text{max}} = RMAX$$

where $\mathbf{x}$ are the coordinates in spherical and cylindrical coordinate systems respectively, $R_{\text{min}}$ is the minimum radial distance from the center of the grid, and $R_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Spherical's RMAX Parameter using NRPy+

```python
# Set default value for Spherical's RMAX parameter
RMAX = 10.0  # units: M

par.set_parval_from_str("reference_metric::Spherical::RMAX", str(RMAX))

print("Default value for Spherical's RMAX parameter set.")
```

This code sets the default value for `RMAX` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Spherical's RMAX Parameter in NRPy+

The `RMAX` parameter determines the maximum radial distance from the center of the grid, which is essential for setting up the simulation. By setting a default value for `RMAX`, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_rmax(self):
        # Set default value for RMAX
        
        self.RMAX = RMAX
        
        print("Default value for RMAX parameter set.")
```

This code uses the `NR**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhSpherical's AMAX Parameter

The `AMAX` parameter is a crucial component of the SinhSpherical coordinate system. It determines the maximum radial distance from the center of the grid, which is essential for setting up the simulation.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$A_{\text{min}} = 0.1M$$
$$A_{\text{max}} = AMAX$$

where $\mathbf{x}$ are the coordinates in SinhSpherical and cylindrical coordinate systems respectively, $A_{\text{min}}$ is the minimum radial distance from the center of the grid, and $A_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhSpherical's AMAX Parameter using NRPy+

```python
# Set default value for SinhSpherical's AMAX parameter
AMAX = 10.0  # units: M

par.set_parval_from_str("reference_metric::SinhSpherical::AMAX", str(AMAX))

print("Default value for SinhSpherical's AMAX parameter set.")
```

This code sets the default value for `AMAX` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhSpherical's AMAX Parameter in NRPy+

The `AMAX` parameter determines the maximum radial distance from the center of the grid, which is essential for setting up the simulation. By setting a default value for `AMAX`, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_amax(self):
        # Set default value for AMAX
        
        self.AMAX = AMAX
        
        print("Default value for AM**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Cartesians' Minimum and Maximum Coordinates

The `Cartesian` coordinate system is a common choice for many simulations. The minimum and maximum coordinates, denoted by `-params.{x,y,z}min` and `-params.{x,y,z}max`, determine the bounds of the grid in each direction.

$$\mathbf{x} = (x, y, z)$$
$$x_{\text{min}} = -XMAX$$
$$x_{\text{max}} = XMAX$$

where $\mathbf{x}$ are the coordinates in Cartesian and cylindrical coordinate systems respectively, $x_{\text{min}}$ is the minimum x-coordinate from the center of the grid, and $x_{\text{max}}$ is the maximum x-coordinate from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Cartesians' Minimum and Maximum Coordinates using NRPy+

```python
# Set default value for Cartesian's minimum and maximum coordinates
XMAX = 10.0  # units: M

par.set_parval_from_str("reference_metric::Cartesian::XMAX", str(XMAX))
par.set_parval_from_str("numerical_grid::x_min", str(-XMAX))
par.set_parval_from_str("numerical_grid::x_max", str(XMAX))

print("Default value for Cartesian's minimum and maximum coordinates set.")
```

This code sets the default values for `XMAX`, `-params.xmin`, and `-params.xmax` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Cartesians' Minimum and Maximum Coordinates in NRPy+

The minimum and maximum coordinates determine the bounds of the grid in each direction. By setting a default value for these parameters, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Cylindrical's Minimum and Maximum Coordinates

The `Cylindrical` coordinate system is a common choice for many simulations. The minimum and maximum coordinates, denoted by `-params.ZMIN` and `.{Z,RHO}MAX`, determine the bounds of the grid in each direction.

$$\mathbf{x} = (r, \theta, z)$$
$$z_{\text{min}} = ZMIN$$
$$r_{\text{max}} = RMAX$$

where $\mathbf{x}$ are the coordinates in cylindrical and Cartesian coordinate systems respectively, $z_{\text{min}}$ is the minimum z-coordinate from the center of the grid, and $r_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Cylindrical's Minimum and Maximum Coordinates using NRPy+

```python
# Set default value for Cylindrical's minimum and maximum coordinates
ZMIN = 10.0  # units: M
RMAX = 20.0  # units: M

par.set_parval_from_str("reference_metric::Cylindrical::ZMIN", str(ZMIN))
par.set_parval_from_str("numerical_grid::z_min", str(ZMIN))
par.set_parval_from_str("numerical_grid::rho_max", str(RMAX))

print("Default value for Cylindrical's minimum and maximum coordinates set.")
```

This code sets the default values for `ZMIN`, `-params.Zmin`, and `.{Z,RHO}max` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Cylindrical's Minimum and Maximum Coordinates in NRPy+

The minimum and maximum coordinates determine the bounds of the grid in each direction. By setting a default value for these parameters, we can easily**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhCylindrical's Amplification Factors

The `SinhCylindrical` coordinate system is a modified version of the Cylindrical system, where the radial distance is amplified by a factor. The amplification factors, denoted by `params.AMPL{RHO,Z}`, determine how much each direction is scaled.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$a_{\text{z}} = AMPL_Z$$

where $\mathbf{x}$ are the coordinates in SinhCylindrical and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $a_{\text{z}}$ is the amplification factor for z-coordinate.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhCylindrical's Amplification Factors using NRPy+

```python
# Set default value for SinhCylindrical's amplification factors
AMPL_RHO = 2.0  # units: dimensionless
AMPL_Z = 1.5  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindrical::AMPL_RHO", str(AMPL_RHO))
par.set_parval_from_str("reference_metric::SinhCylindrical::AMPL_Z", str(AMPL_Z))

print("Default value for SinhCylindrical's amplification factors set.")
```

This code sets the default values for `params.AMPL{RHO,Z}` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhCylindrical's Amplification Factors in NRPy+

The amplification factors determine how much each direction is scaled. By setting a default value for these parameters, we can easily modify various aspects**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SymTP's AMAX Parameter

The `SymTP` (Symmetrized Temporal Parity) coordinate system is a type of coordinate system that is commonly used in numerical relativity. The `AMAX` parameter determines the maximum radial distance from the center of the grid.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$A_{\text{max}} = AMAX$$

where $\mathbf{x}$ are the coordinates in SymTP and Cartesian coordinate systems respectively, $A_{\text{max}}$ is the maximum radial distance from the center of the grid.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SymTP's AMAX Parameter using NRPy+

```python
# Set default value for SymTP's AMAX parameter
domain_size = 3.0  # units: M

par.set_parval_from_str("reference_metric::SymTP::AMAX", str(domain_size))

print("Default value for SymTP's AMAX parameter set.")
```

This code sets the default value for `params.AMAX` using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SymTP's AMAX Parameter in NRPy+

The `AMAX` parameter determines the maximum radial distance from the center of the grid. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_amax(self):
        # Set default value for AMAX
        
        self.AMAX = domain_size
        
        print("Default value for AMAX parameter set.")
```

This code uses the `NRPy_code_generator` class to interact with the coordinate system parameters.

Note that the `domain_size` variable is used to set the default value for `params.AMAX`.**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Sinh Width Parameter

The `sinh_width` parameter is a crucial component of the SinhCylindrical and SinhSpherical coordinate systems. It determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = BWIDTH$$

where $\mathbf{x}$ are the coordinates in SinhCylindrical and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the sinh width parameter.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Sinh Width Parameter using NRPy+

```python
# Set default value for sinh width parameter
sinh_width = 1.5  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindrical::BWIDTH", str(sinh_width))
par.set_parval_from_str("reference_metric::SinhSpherical::BWIDTH", str(sinh_width))

print("Default value for sinh width parameter set.")
```

This code sets the default values for `params.BWIDTH` in both SinhCylindrical and SinhSpherical coordinate systems using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Sinh Width Parameter in NRPy+

The sinh width parameter determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_sinh_width(self):
        # Set default value for sinh width
        
        self.sinh_width = sinh_width
        
**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhSpherical's SINHW Parameter

The `SINHW` parameter is a crucial component of the SinhSpherical coordinate system. It determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid.

$$\mathbf{x} = (r, \theta, \phi)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = SINHW$$

where $\mathbf{x}$ are the coordinates in SinhSpherical and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the sinh width parameter.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhSpherical's SINHW Parameter using NRPy+

```python
# Set default value for SinhSpherical's SINHW parameter
sinh_width = 1.5  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhSpherical::SINHW", str(sinh_width))

print("Default value for SinhSpherical's SINHW parameter set.")
```

This code sets the default value for `params.SINHW` in the SinhSpherical coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhSpherical's SINHW Parameter in NRPy+

The sinh width parameter determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_sinh_width(self):
        # Set default value for sinh width
        
        self.sinh_width = sinh_width
        
        print("Default value for SINHW parameter set.")
```

This**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhCylindrical's SINHW{RHO,Z} Parameter

The `SINHW` parameter is a crucial component of the SinhCylindrical coordinate system. It determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = SINHW_RHO$$
$$b_z = SINHW_Z$$

where $\mathbf{x}$ are the coordinates in SinhCylindrical and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, $b_\rho$ is the sinh width parameter for radial direction, and $b_z$ is the sinh width parameter for z-direction.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhCylindrical's SINHW{RHO,Z} Parameter using NRPy+

```python
# Set default value for SinhCylindrical's SINHW{RHO,Z} parameter
sinh_width_rho = 1.5  # units: dimensionless
sinh_width_z = 2.0  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindrical::SINHW_RHO", str(sinh_width_rho))
par.set_parval_from_str("reference_metric::SinhCylindrical::SINHW_Z", str(sinh_width_z))

print("Default value for SinhCylindrical's SINHW{RHO,Z} parameter set.")
```

This code sets the default values for `params.SINHW_RHO` and `params.SINHW_Z` in the SinhCylindrical coordinate system using the `set_parval_from_str()` function.

### Theory Review**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhSymTP's SINHWAA Parameter

The `SINHWAA` parameter is a crucial component of the SinhSymTP (Sinh Symmetrized Temporal Parity) coordinate system. It determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = SINHWAA$$

where $\mathbf{x}$ are the coordinates in SinhSymTP and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the sinh width parameter.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhSymTP's SINHWAA Parameter using NRPy+

```python
# Set default value for SinhSymTP's SINHWAA parameter
sinh_width = 0.4  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhSymTP::SINHWAA", str(sinh_width))

print("Default value for SinhSymTP's SINHWAA parameter set.")
```

This code sets the default value for `params.SINHWAA` in the SinhSymTP coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhSymTP's SINHWAA Parameter in NRPy+

The sinh width parameter determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_sinh_width(self):
        # Set default value for sinh width
        
       **Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Sinh* Coordinates

If Sinh* coordinates (SinhCylindrical, SinhSpherical, or SinhSymTP) are chosen, additional parameters need to be set. The `params.SINHWAA` parameter is a crucial component of these coordinate systems.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = SINHWAA$$

where $\mathbf{x}$ are the coordinates in Sinh* and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the sinh width parameter.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Sinh* Coordinates using NRPy+

```python
# Set default value for sinh width parameter
sinh_width_aa = 0.4  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindrical::SINHWAA", str(sinh_width_aa))
par.set_parval_from_str("reference_metric::SinhSpherical::SINHWAA", str(sinh_width_aa))
par.set_parval_from_str("reference_metric::SinhSymTP::SINHWAA", str(sinh_width_aa))

print("Default value for sinh width parameter set.")
```

This code sets the default values for `params.SINHWAA` in SinhCylindrical, SinhSpherical, and SinhSymTP coordinate systems using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Sinh* Coordinates in NRPy+

The sinh width parameter determines the width of the sinh function, which is used to amplify the radial distance from the center of the grid. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Sinhv2 Const Dr Parameter

The `sinhv2_const_dr` parameter is a crucial component of the Sinh* (SinhCylindrical, SinhSpherical, or SinhSymTP) coordinate systems. It determines the constant factor used in the sinh function for radial distance.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = sinhv2_const_dr$$

where $\mathbf{x}$ are the coordinates in Sinh* and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the constant factor used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Sinhv2 Const Dr Parameter using NRPy+

```python
# Set default value for sinh v2 constant dr parameter
sinhv2_const_dr = 1.0  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindrical::sinhv2_const_dr", str(sinhv2_const_dr))
par.set_parval_from_str("reference_metric::SinhSpherical::sinhv2_const_dr", str(sinhv2_const_dr))
par.set_parval_from_str("reference_metric::SinhSymTP::sinhv2_const_dr", str(sinhv2_const_dr))

print("Default value for sinh v2 constant dr parameter set.")
```

This code sets the default values for `sinhv2_const_dr` in SinhCylindrical, SinhSpherical, and SinhSymTP coordinate systems using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Sinhv2 Const Dr Parameter in NRPy+

The constant factor used in the sinh function determines the amplification of the radial distance. By**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhSphericalv2's const_dr Parameter

The `const_dr` parameter is a crucial component of the SinhSphericalv2 (Sinh Spherical v2) coordinate system. It determines the constant factor used in the sinh function for radial distance.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b = const_dr$$

where $\mathbf{x}$ are the coordinates in SinhSphericalv2 and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b$ is the constant factor used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhSphericalv2's const_dr Parameter using NRPy+

```python
# Set default value for sinh spherical v2 constant dr parameter
const_dr = 1.0  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhSphericalv2::const_dr", str(const_dr))

print("Default value for sinh spherical v2 constant dr parameter set.")
```

This code sets the default value for `params.const_dr` in the SinhSphericalv2 coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhSphericalv2's const_dr Parameter in NRPy+

The constant factor used in the sinh function determines the amplification of the radial distance. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_const_dr(self):
        # Set default value for const dr
        
        self.const_dr = const_dr
        
        print("Default value for const dr parameter set.")
```

This**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhCylindricalv2's const_d{rho,z} Parameter

The `const_d` parameter is a crucial component of the SinhCylindricalv2 (Sinh Cylindrical v2) coordinate system. It determines the constant factor used in the sinh function for radial distance.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = const_drho$$
$$b_z = const_dz$$

where $\mathbf{x}$ are the coordinates in SinhCylindricalv2 and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b_\rho$ and $b_z$ are the constant factors used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhCylindricalv2's const_d{rho,z} Parameter using NRPy+

```python
# Set default value for sinh cylindrical v2 constant dr parameter
sinhv2_const_dr_rho = 0.05  # units: dimensionless
sinhv2_const_dr_z = 0.1  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhCylindricalv2::const_drho", str(sinhv2_const_dr_rho))
par.set_parval_from_str("reference_metric::SinhCylindricalv2::const_dz", str(sinhv2_const_dr_z))

print("Default value for sinh cylindrical v2 constant dr parameter set.")
```

This code sets the default values for `params.const_drho` and `params.const_dz` in the SinhCylindricalv2 coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhCyl**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to Sinh*v2 Coordinates

If Sinh*v2 (SinhCylindricalv2 or SinhSphericalv2) coordinates are chosen, additional parameters need to be set. The `const_d` parameter is a crucial component of these coordinate systems.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = const_drho$$
$$b_z = const_dz$$

where $\mathbf{x}$ are the coordinates in Sinh*v2 and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b_\rho$ and $b_z$ are the constant factors used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for Sinh*v2 Coordinates using NRPy+

```python
# Set default value for sinh v2 constant dr parameter
sinhv2_const_dr_rho = 0.05  # units: dimensionless
sinhv2_const_dr_z = 0.1  # units: dimensionless

if "SinhCylindricalv2" in par.reference_metric:
    par.set_parval_from_str("reference_metric::SinhCylindricalv2::const_drho", str(sinhv2_const_dr_rho))
    par.set_parval_from_str("reference_metric::SinhCylindricalv2::const_dz", str(sinhv2_const_dr_z))

elif "SinhSphericalv2" in par.reference_metric:
    par.set_parval_from_str("reference_metric::SinhSphericalv2::const_drho", str(sinhv2_const_dr_rho))
    par.set_parval_from_str("reference_metric::SinhSphericalv2::const_dz", str(sinhv2**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SymTP_bScale Parameter

The `SymTP_bScale` parameter is a crucial component of the SymTP (Symmetrized Temporal Parity) coordinate system. It determines the scale factor used in the sinh function for radial distance.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = SymTP_bScale$$

where $\mathbf{x}$ are the coordinates in SymTP and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b_\rho$ is the scale factor used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SymTP_bScale Parameter using NRPy+

```python
# Set default value for symtp b scale parameter
SymTP_bScale = 1.0  # units: dimensionless

par.set_parval_from_str("reference_metric::SymTP::SymTP_bScale", str(SymTP_bScale))

print("Default value for SymTP_bScale parameter set.")
```

This code sets the default value for `params.SymTP_bScale` in the SymTP coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SymTP_bScale Parameter in NRPy+

The scale factor used in the sinh function determines the amplification of the radial distance. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_symtp_b_scale(self):
        # Set default value for symtp b scale
        
        self.SymTP_bScale = SymTP_bScale
        
        print("Default value for SymTP_bScale parameter set.")
``**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SinhSymTP's bScale Parameter

The `bScale` parameter is a crucial component of the SinhSymTP (Sinh Symmetrized Temporal Parity) coordinate system. It determines the scale factor used in the sinh function for radial distance.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = bScale$$

where $\mathbf{x}$ are the coordinates in SinhSymTP and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b_\rho$ is the scale factor used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SinhSymTP's bScale Parameter using NRPy+

```python
# Set default value for sinh symtp b scale parameter
SymTP_bScale = 0.5  # units: dimensionless

par.set_parval_from_str("reference_metric::SinhSymTP::bScale", str(SymTP_bScale))

print("Default value for SinhSymTP's bScale parameter set.")
```

This code sets the default value for `params.bScale` in the SinhSymTP coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SinhSymTP's bScale Parameter in NRPy+

The scale factor used in the sinh function determines the amplification of the radial distance. By setting a default value for this parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_b_scale(self):
        # Set default value for b scale
        
        self.bScale = SymTP_bScale
        
        print("Default value for bScale parameter set.")
```

This**Step 2.a: Set Defaults for Coordinate System Parameters**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting default values for coordinate system parameters in NRPy+. These parameters are used to configure the reference metric and numerical grid.

### Theory Review

#### Introduction to SymTP Coordinate System

If SymTP (Symmetrized Temporal Parity) is chosen, additional parameters need to be set. The `bScale` parameter is a crucial component of this coordinate system.

$$\mathbf{x} = (\sigma, \eta, \zeta)$$
$$a_{\text{rho}} = AMPL_RHO$$
$$b_\rho = bScale$$

where $\mathbf{x}$ are the coordinates in SymTP and Cartesian coordinate systems respectively, $a_{\text{rho}}$ is the amplification factor for radial distance, and $b_\rho$ is the scale factor used in the sinh function.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Default Value for SymTP's bScale Parameter using NRPy+

```python
# Set default value for symtp b scale parameter
SymTP_bScale = 0.5  # units: dimensionless

if "SymTP" in par.reference_metric:
    par.set_parval_from_str("reference_metric::SymTP::bScale", str(SymTP_bScale))

print("Default value for SymTP's bScale parameter set.")
```

This code sets the default value for `params.bScale` in the SymTP coordinate system using the `set_parval_from_str()` function.

### Theory Review

#### Understanding SymTP Coordinate System in NRPy+

The SymTP coordinate system is a type of Sinh* (SinhCylindrical, SinhSpherical, or SinhSymTP) coordinate system. It uses a sinh-like function to amplify the radial distance. By setting a default value for the `bScale` parameter, we can easily modify various aspects of the simulation without having to delve into the underlying code.

```python
class NRPy_code_generator:
    def set_b_scale(self):
        # Set default value for b scale
        
        self.b**Step 2.c: Set the Order of Spatial Derivatives and Time Stepping**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting the order of spatial derivatives and time stepping in NRPy+. These parameters are used to configure the numerical grid and simulation settings.

### Theory Review

#### Introduction to Spatial Derivatives

Spatial derivatives are used to compute the derivative of a function with respect to space. In NRPy+, we can set the order of these derivatives using the `SPATIAL_DER_ORDER` parameter.

$$\frac{\partial f}{\partial x} = \frac{f(x + h) - f(x)}{h}$$

where $f$ is the function, $x$ is the spatial coordinate, and $h$ is the step size.

#### Introduction to Time Stepping

Time stepping refers to the process of advancing the simulation in time. In NRPy+, we can set the order of time stepping using the `TIME_STEPPING_ORDER` parameter.

$$\frac{d f}{d t} = \frac{f(t + \Delta t) - f(t)}{\Delta t}$$

where $f$ is the function, $t$ is the time coordinate, and $\Delta t$ is the time step size.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Order of Spatial Derivatives and Time Stepping using NRPy+

```python
# Set order of spatial derivatives
SPATIAL_DER_ORDER = 4

# Set order of time stepping
TIME_STEPPING_ORDER = 3

par.set_parval_from_str("finite_difference::SPATIAL_DER_ORDER", str(SPATIAL_DER_ORDER))
par.set_parval_from_str("finite_difference::TIME_STEPPING_ORDER", str(TIME_STEPPING_ORDER))

print("Order of spatial derivatives and time stepping set.")
```

This code sets the order of spatial derivatives and time stepping using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Spatial Derivatives and Time Stepping in NRPy+

By setting the order of spatial derivatives and time stepping, we can control**Step 2.c: Set the Order of Spatial Derivatives and Time Stepping**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting the order of spatial derivatives and time stepping in NRPy+. These parameters are used to configure the numerical grid and simulation settings.

### Theory Review

#### Introduction to Core Data Type

The core data type is a fundamental concept in NRPy+, which determines the accuracy of the finite difference scheme. In this step, we will set the order of the core data type using the `FD_order` parameter.

$$\frac{\partial f}{\partial x} = \frac{f(x + h) - f(x)}{h}$$

where $f$ is the function, $x$ is the spatial coordinate, and $h$ is the step size.

#### Order of Core Data Type

The order of the core data type determines the accuracy of the finite difference scheme. A higher order means a more accurate result but also increases the computational cost.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Order of Core Data Type using NRPy+

```python
# Set order of core data type
FD_order = 4

par.set_parval_from_str("finite_difference::FD_order", str(FD_order))

print("Order of core data type set.")
```

This code sets the order of the core data type using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Order of Core Data Type in NRPy+

By setting the order of the core data type, we can control the accuracy of the finite difference scheme. A higher order means a more accurate result but also increases the computational cost.

```python
class NRPy_code_generator:
    def set_fd_order(self):
        # Set order of core data type
        
        self.FD_order = FD_order
        
        print("Order of core data type set.")
```

This code sets the order of the core data type and prints a message to confirm that it has been set.**Step 2.c: Set the Order of Spatial Derivatives and Time Stepping**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting the order of spatial derivatives and time stepping in NRPy+. These parameters are used to configure the numerical grid and simulation settings.

### Theory Review

#### Introduction to Finite Difference Order

The finite difference order is a fundamental concept in NRPy+, which determines the accuracy of the finite difference scheme. In this step, we will set the order of the finite difference scheme using the `FD_order` parameter.

$$\frac{\partial f}{\partial x} = \frac{f(x + h) - f(x)}{h}$$

where $f$ is the function, $x$ is the spatial coordinate, and $h$ is the step size.

#### Finite Difference Order Restrictions

The finite difference order must be an even number, starting with 2. This means that the order can be 2, 4, 6, 8, 10, or 12, but not 1, 3, 5, etc.

#### Stability of Higher-Order Schemes

Higher-order schemes are generally more accurate than lower-order schemes, but they can also be less stable. In particular, a finite difference order of 12 is often unstable and should be avoided unless absolutely necessary.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Finite Difference Order and Data Type using NRPy+

```python
# Set finite difference order (even numbers only)
FD_order = 4

# Set data type (e.g. "double", "float")
REAL = "double"

par.set_parval_from_str("finite_difference::FD_order", str(FD_order))
par.set_parval_from_str("data_type::REAL", REAL)

print("Finite difference order and data type set.")
```

This code sets the finite difference order and data type using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Finite Difference Order in NRPy+

By setting the finite difference order, we can control the accuracy of the finite difference scheme. However, we must also be aware of the stability**Step 2.c: Set the Order of Spatial Derivatives and Time Stepping**
=============================================================

### Overview of the Notebook

This notebook covers the process of setting the order of spatial derivatives and time stepping in NRPy+. These parameters are used to configure the numerical grid and simulation settings.

### Theory Review

#### Introduction to Data Type Selection

When setting up a simulation in NRPy+, we need to choose a data type for our variables. In this step, we will select the data type using the `REAL` parameter.

#### Why Choose Double?

The `REAL` parameter should be set to `"double"` because it is generally more accurate and stable than other options.

$$\text{Accuracy} = f(\text{data type})$$

where $f$ is a function that depends on the data type chosen.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Data Type using NRPy+

```python
# Set data type (e.g. "double", "float")
REAL = "double"

par.set_parval_from_str("data_type::REAL", REAL)

print("Data type set.")
```

This code sets the data type using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Data Type Selection in NRPy+

By selecting the correct data type, we can ensure that our simulation is accurate and stable. In this case, choosing `"double"` provides a high level of accuracy and stability.

```python
class NRPy_code_generator:
    def set_real(self):
        # Set data type
        
        self.REAL = REAL
        
        print("Data type set.")
```

This code sets the data type and prints a message to confirm that it has been set.**Step 5: Set the Finite Differencing Order**
=============================================

### Overview of the Notebook

This notebook covers the process of setting the finite differencing order in NRPy+. The finite differencing order determines the accuracy and stability of the simulation.

### Theory Review

#### Introduction to Finite Differencing Order

The finite differencing order is a critical parameter in NRPy+ that determines the accuracy and stability of the simulation. In this step, we will set the finite differencing order using the `FD_order` parameter.

$$\text{Accuracy} = f(\text{FD_order})$$

where $f$ is a function that depends on the finite differencing order chosen.

#### Importance of Choosing Correct Finite Differencing Order

Choosing the correct finite differencing order is crucial for achieving accurate and stable results in NRPy+. A higher finite differencing order can provide more accurate results but may also increase the computational cost and reduce stability.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Setting Finite Differencing Order using NRPy+

```python
# Set finite differencing order
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order)

print("Finite differencing order set.")
```

This code sets the finite differencing order using the `set_parval_from_str()` function.

### Theory Review

#### Understanding Finite Differencing Order in NRPy+

By setting the correct finite differencing order, we can ensure that our simulation is accurate and stable. In this case, choosing a higher finite differencing order (e.g. 4) provides more accuracy but may increase computational cost and reduce stability.

```python
class NRPy_code_generator:
    def set_fd_order(self):
        # Set finite differencing order
        
        self.FD_CENTDERIVS_ORDER = FD_order
        
        print("Finite differencing order set.")
```

This code sets the finite differencing order and prints a message to confirm that it has been set.

### Mathematical Background

The finite differencing order is related to the number of points used in the numerical differentiation. A higher finite differencing order corresponds to using more points, which can improve**Step 7: Set Finite Difference Functionality**
=============================================

### Overview of the Notebook

This notebook covers the process of setting finite difference functionality in NRPy+. This includes enabling or disabling the use of finite difference functions and registering C functions for the radiative fluid module.

### Theory Review

#### Introduction to Finite Difference Functions

Finite difference functions are used in NRPy+ to perform numerical differentiation. These functions can be enabled or disabled depending on the needs of the simulation.

$$\text{Numerical Derivative} = \frac{f(x + h) - f(x)}{h}$$

where $f$ is the function, $x$ is the spatial coordinate, and $h$ is the step size.

#### Enabling or Disabling Finite Difference Functions

The `enable_FD_functions` parameter determines whether finite difference functions are enabled or disabled. By default, this parameter is set to `False`.

```python
if enable_FD_functions:
    par.set_parval_from_str("finite_difference::enable_FD_functions", enable_FD_functions)
```

This code checks if `enable_FD_functions` is `True`, and if so, sets the `finite_difference::enable_FD_functions` parameter accordingly.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

### Theory Review

#### Understanding Finite Difference Functions in NRPy+

By enabling or disabling finite difference functions, we can control the accuracy and stability of our simulation. This is particularly important when working with complex systems that require high precision calculations.

```python
class NRPy_code_generator:
    def set_fd_functions(self):
        # Set finite difference functions
        
        self.enable_FD_functions = enable_FD_functions
        
        print("Finite difference functions enabled.")
```

This code enables or disables finite difference functions and prints a message to confirm that it has been set.

### Registering C Functions for the Radiative Fluid Module

The `rfm.register_C_functions()` function is used to register C functions for the radiative fluid module. This includes registering precompute functions and time-stepping functions.

```python
rfm.register_C_functions(enable_rfm_precompute=False, use_unit_wavespeed_for_find_timestep=True)
```

This code registers the necessary C**Step 3: Import Black Hole ADM Initial Data C Function**
==========================================================

### Overview of the Notebook

This notebook covers the process of importing the Black Hole ADM initial data C function from the NRPy+ module. This function is used to generate initial data for black hole simulations.

### Theory Review

#### Introduction to Black Hole ADM Initial Data

The Black Hole ADM (Arnowitt-Deser-Misner) initial data is a method for generating initial data for black hole simulations. It uses a set of equations to determine the metric and extrinsic curvature of the spacetime at a given time.

$$g_{\mu\nu} = \left( \begin{array}{cc}
\alpha^2 & 0 \\
0 & (\beta_i)^2
\end{array} \right)$$

where $g_{\mu\nu}$ is the metric tensor, $\alpha$ is the lapse function, and $\beta_i$ are the shift vector components.

#### Importing C Function from NRPy+ Module

The `importBlackHoleADM` function is used to import the Black Hole ADM initial data C function from the NRPy+ module. This function takes no arguments and returns the necessary C code for generating the initial data.

```python
import_NRPy_module = "NRPy/black_holes/BlackHoleADM.py"
import_NRPy_function = "importBlackHoleADM"
```

This code imports the `BlackHoleADM` module from the NRPy+ repository and specifies the name of the function to be imported, which is `importBlackHoleADM`.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Importing Black Hole ADM Initial Data C Function

```python
import_NRPy_module = "NRPy/black_holes/BlackHoleADM.py"
import_NRPy_function = "importBlackHoleADM"

c_code = NR.import_C_function(import_NRPy_module, import_NRPy_function)

print("Black Hole ADM initial data C function imported.")
```

This code imports the `BlackHoleADM` module from the NRPy+ repository and specifies the name of the function to be imported. It then uses the**Importing Black Hole ADM Initial Data C Function**
======================================================

### Overview of the Notebook

This notebook covers the process of importing the Black Hole ADM initial data C function from the NRPy+ module. This function is used to generate initial data for black hole simulations.

### Theory Review

#### Introduction to Black Hole ADM Initial Data

The Black Hole ADM (Arnowitt-Deser-Misner) initial data is a method for generating initial data for black hole simulations. It uses a set of equations to determine the metric and extrinsic curvature of the spacetime at a given time.

$$g_{\mu\nu} = \left( \begin{array}{cc}
\alpha^2 & 0 \\
0 & (\beta_i)^2
\end{array} \right)$$

where $g_{\mu\nu}$ is the metric tensor, $\alpha$ is the lapse function, and $\beta_i$ are the shift vector components.

#### Importing C Function from NRPy+ Module

The `importBlackHoleADM` function is used to import the Black Hole ADM initial data C function from the NRPy+ module. This function takes no arguments and returns the necessary C code for generating the initial data.

```python
import_NRPy_module = "NRPy/black_holes/BlackHoleADM.py"
import_NRPy_function = "importBlackHoleADM"
```

This code imports the `BlackHoleADM` module from the NRPy+ repository and specifies the name of the function to be imported, which is `importBlackHoleADM`.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_code_generator as NR   # Import NRPy code generator module for finite-difference functions
```

This line of code imports the `NRPy_code_generator` module, which provides a way to generate C code for finite-difference functions.

#### Importing Black Hole ADM Initial Data C Function

```python
import_NRPy_module = "NRPy/black_holes/BlackHoleADM.py"
import_NRPy_function = "importBlackHoleADM"

c_code = NR.import_C_function(import_NRPy_module, import_NRPy_function)

print("Black Hole ADM initial data C function imported.")
```

This code imports the `BlackHoleADM` module from the NRPy+ repository and specifies the name of the function to be imported. It then uses the `import_C**Importing Initial Data and Registering C Functions**
======================================================

### Overview of the Notebook

This notebook covers the process of importing initial data for a black hole simulation and registering C functions for the radiative fluid module. This includes importing the Black Hole ADM (Arnowitt-Deser-Misner) initial data and registering the C function for converting from ADM to BSSN (Baumgarte-Shapiro-Simon-Schoen-Numeric) coordinates.

### Theory Review

#### Introduction to Initial Data Importing

Initial data is used to specify the initial conditions of a simulation. In this case, we are importing the Black Hole ADM initial data.

$$\text{ADM Initial Data} = \left( \begin{array}{c}
\alpha \\
\beta_i \\
\gamma_{ij} \\
K_{ij}
\end{array} \right)$$

where $\alpha$ is the lapse function, $\beta_i$ are the shift vector components, $\gamma_{ij}$ is the metric tensor, and $K_{ij}$ is the extrinsic curvature.

#### Registering C Functions for Radiative Fluid Module

The radiative fluid module uses a set of C functions to perform calculations. We need to register these functions with the NRPy+ module.

```python
import BSSN.ADM_Initial_Data_Reader__BSSN_Converter as IDread

IDread.add_to_Cfunction_dict_exact_ADM_ID_function(dictID[initial_data_string].functionname,
                                                   dictID[initial_data_string].OrigCoordSystem,
                                                   IDmodule.alpha, IDmodule.betaU, IDmodule.BU,
                                                   IDmodule.gammaDD, IDmodule.KDD)
```

This code adds the C function for converting from ADM to BSSN coordinates to the dictionary.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import time
import importlib
starttime = time.time()
IDmodule = importlib.import_module(dictID[initial_data_string].modulename)
IDfunc = getattr(IDmodule, dictID[initial_data_string].functionname)
IDfunc()
```

This code imports the necessary modules and executes the initial data function.

#### Registering C Functions for Radiative Fluid Module

```python
import BSSN.ADM_Initial_Data_Reader__BSSN_Converter as IDread

IDread.add_to_Cfunction_dict_exact_ADM_ID_function(dictID[initial_data_string].functionname**Step 3.a: Output C Codes for Declaring and Setting C Parameters**
====================================================================

### Overview of the Notebook

This notebook covers the process of outputting C codes needed for declaring and setting C parameters in NRPy+. This includes generating C code for declaring and setting C parameters, as well as updating the `free_parameters.h` file.

### Theory Review

#### Introduction to C Parameters

C parameters are used to specify the values of various physical constants and parameters in a simulation. In this step, we will generate C code for declaring and setting these parameters.

```python
# Step 3.a: Output C codes needed for declaring and setting Cparameters;
par.set_parval_from_str("output::generate_C_code", "declare_and_set_Cparams")
```

This line of code outputs the necessary C code for declaring and setting C parameters.

#### Updating `free_parameters.h` File

The `free_parameters.h` file is used to store the values of free parameters in a simulation. In this step, we will update this file with the new values.

```python
# Step 3.a: Also set `free_parameters.h`
par.set_parval_from_str("output::update_free_parameters_h", "True")
```

This line of code updates the `free_parameters.h` file with the new values.

### Code Implementation

#### Outputting C Codes for Declaring and Setting C Parameters

```python
import_NRPy_module = "NRPy/parameters.py"
c_code = NRPy.generate_C_code("declare_and_set_Cparams")
print(c_code)
```

This code imports the necessary module, generates the C code for declaring and setting C parameters, and prints it to the console.

#### Updating `free_parameters.h` File

```python
import_NRPy_module = "NRPy/parameters.py"
c_code = NRPy.update_free_parameters_h("True")
print(c_code)
```

This code updates the `free_parameters.h` file with the new values and prints them to the console.

### Mathematical Background

The C parameters are used to specify the values of various physical constants and parameters in a simulation. These parameters can be thought of as free variables that need to be set before running a simulation.

$$\text{C Parameters} = \left( \begin{array}{c}
c_1 \\
c_2 \\
...
\end{array} \right)$$

where $c_i$ are the C parameters.

The `free**Step 3: Generate C Codes for Declaring and Setting C Parameters**
=====================================================================

### Overview of the Notebook

This notebook covers the process of generating C codes for declaring and setting C parameters in NRPy+. This includes generating C code for declaring and setting C parameters, as well as updating the `free_parameters.h` file.

### Theory Review

#### Introduction to C Parameters

C parameters are used to specify the values of various physical constants and parameters in a simulation. In this step, we will generate C code for declaring and setting these parameters.

```python
# Step 3: Generate C codes for declaring and setting Cparameters,
# based on declared NRPy+ Cparameters.
declare_Cparams_struct = NR.generate_C_code("declare_Cparameters_struct")
set_Cparams_default = NR.generate_C_code("set_Cparameters_default")
```

This code generates the C code for declaring and setting C parameters.

#### Generating `free_parameters.h` File

The `free_parameters.h` file is used to store the values of free parameters in a simulation. In this step, we will update this file with the new values.

```python
# Step 3: Then output `free_parameters.h`, which sets initial data parameters,
# as well as grid domain & reference metric parameters.
output_free_parameters_h = NR.update_free_parameters_h("True")
```

This code updates the `free_parameters.h` file with the new values.

### Code Implementation

#### Generating C Codes for Declaring and Setting C Parameters

```python
import_NRPy_module = "NRPy/parameters.py"
declare_Cparams_struct = NR.generate_C_code("declare_Cparameters_struct")
set_Cparams_default = NR.generate_C_code("set_Cparameters_default")
print(declare_Cparams_struct)
print(set_Cparams_default)
```

This code imports the necessary module, generates the C code for declaring and setting C parameters, and prints it to the console.

#### Updating `free_parameters.h` File

```python
import_NRPy_module = "NRPy/parameters.py"
output_free_parameters_h = NR.update_free_parameters_h("True")
print(output_free_parameters_h)
```

This code updates the `free_parameters.h` file with the new values and prints them to the console.

### Mathematical Background

The C parameters are used to specify the values of various physical constants and parameters in a simulation. These parameters can be thought of as free variables that need to be set before running a simulation.

$$\**Step 3.a.i: Set `free_parameters.h`**
======================================

### Overview of the Notebook

This notebook covers the process of setting the `free_parameters.h` file in NRPy+. This includes updating the file with the new values of free parameters.

### Theory Review

#### Introduction to Free Parameters

Free parameters are used to specify the values of various physical constants and parameters in a simulation. In this step, we will update the `free_parameters.h` file with the new values.

```python
# Step 3.a.i: Set `free_parameters.h`
output_free_parameters_h = NR.update_free_parameters_h("True")
```

This code updates the `free_parameters.h` file with the new values.

### Code Implementation

#### Setting `free_parameters.h` File

```python
import_NRPy_module = "NRPy/parameters.py"
output_free_parameters_h = NR.update_free_parameters_h("True")
print(output_free_parameters_h)
```

This code imports the necessary module, updates the `free_parameters.h` file with the new values, and prints them to the console.

### Mathematical Background

The free parameters are used to specify the values of various physical constants and parameters in a simulation. These parameters can be thought of as free variables that need to be set before running a simulation.

$$\text{Free Parameters} = \left( \begin{array}{c}
f_1 \\
f_2 \\
...
\end{array} \right)$$

where $f_i$ are the free parameters.

### Theory Review

#### Understanding `free_parameters.h` File

The `free_parameters.h` file is used to store the values of free parameters in a simulation. This file needs to be updated with the new values before running the simulation.

```python
# Set `free_parameters.h`
output_free_parameters_h = NR.update_free_parameters_h("True")
```

This code updates the `free_parameters.h` file with the new values.

### Code Implementation

#### Updating `free_parameters.h` File

```python
import_NRPy_module = "NRPy/parameters.py"
output_free_parameters_h = NR.update_free_parameters_h("True")
print(output_free_parameters_h)
```

This code imports the necessary module, updates the `free_parameters.h` file with the new values, and prints them to the console.

### Example Use Case

To update the `free_parameters.h` file with the new values, you can use the following code**Outputting Reference Metric Parameters**
==========================================

### Overview of the Notebook

This notebook covers the process of outputting reference metric parameters in NRPy+. This includes generating C code for the `free_parameters.h` file and updating the reference metric parameters.

### Theory Review

#### Introduction to Reference Metric Parameters

Reference metric parameters are used to specify the values of various physical constants and parameters in a simulation. In this step, we will output these parameters to the `Ccodesdir/free_parameters.h` file.

```python
# Output to $Ccodesdir/free_parameters.h reference metric parameters based on generic
output_reference_metric_params = NR.output_reference_metric_params()
```

This code generates the C code for the reference metric parameters and outputs it to the `free_parameters.h` file.

### Code Implementation

#### Generating C Code for Reference Metric Parameters

```python
import_NRPy_module = "NRPy/parameters.py"
output_reference_metric_params = NR.output_reference_metric_params()
print(output_reference_metric_params)
```

This code imports the necessary module, generates the C code for the reference metric parameters, and prints it to the console.

### Mathematical Background

The reference metric parameters are used to specify the values of various physical constants and parameters in a simulation. These parameters can be thought of as free variables that need to be set before running a simulation.

$$\text{Reference Metric Parameters} = \left( \begin{array}{c}
r \\
g_{xx} \\
g_{yy} \\
...
\end{array} \right)$$

where $r$ are the reference metric parameters.

### Theory Review

#### Understanding Reference Metric Parameters

The reference metric parameters are used to specify the values of various physical constants and parameters in a simulation. This includes the radius, x and y components of the metric tensor.

```python
# Output to $Ccodesdir/free_parameters.h reference metric parameters based on generic
output_reference_metric_params = NR.output_reference_metric_params()
```

This code generates the C code for the reference metric parameters and outputs it to the `free_parameters.h` file.

### Code Implementation

#### Updating `free_parameters.h` File with Reference Metric Parameters

```python
import_NRPy_module = "NRPy/parameters.py"
output_reference_metric_params = NR.output_reference_metric_params()
print(output_reference_metric_params)
```

This code imports the necessary module, generates the C code for the reference metric parameters, and prints it to the console.

### Example**Outputting Domain Size and Sinh Width Parameters**
=====================================================

### Overview of the Notebook

This notebook covers the process of outputting domain size and sinh width parameters in NRPy+. This includes generating C code for these parameters and updating them in the `free_parameters.h` file.

### Theory Review

#### Introduction to Domain Size and Sinh Width Parameters

Domain size and sinh width parameters are used to specify the values of various physical constants and parameters in a simulation. In this step, we will output these parameters to the `Ccodesdir/free_parameters.h` file.

```python
# Output to $Ccodesdir/free_parameters.h domain_size,sinh_width,sinhv2_const_dr,SymTP_bScale,
output_domain_size_and_sinh_width_params = NR.output_domain_size_and_sinh_width_params()
```

This code generates the C code for these parameters and outputs it to the `free_parameters.h` file.

### Code Implementation

#### Generating C Code for Domain Size and Sinh Width Parameters

```python
import_NRPy_module = "NRPy/parameters.py"
output_domain_size_and_sinh_width_params = NR.output_domain_size_and_sinh_width_params()
print(output_domain_size_and_sinh_width_params)
```

This code imports the necessary module, generates the C code for these parameters, and prints it to the console.

### Mathematical Background

The domain size and sinh width parameters are used to specify the values of various physical constants and parameters in a simulation. These parameters can be thought of as free variables that need to be set before running a simulation.

$$\text{Domain Size Parameters} = \left( \begin{array}{c}
domain\_size \\
sinh\_width \\
...
\end{array} \right)$$

where $domain\_size$ and $sinh\_width$ are the domain size parameters.

$$\text{Sinh Width Parameters} = \left( \begin{array}{c}
sinhv2\_const\_dr \\
SymTP\_bScale \\
...
\end{array} \right)$$

where $sinhv2\_const\_dr$ and $SymTP\_bScale$ are the sinh width parameters.

### Theory Review

#### Understanding Domain Size and Sinh Width Parameters

The domain size and sinh width parameters are used to specify the values of various physical constants and parameters in a simulation. This includes the domain size, sinh width, sinh v2 constant dr, and SymTP b scale.

```**Updating `free_parameters.h` File**
=====================================

### Overview of the Notebook

This notebook covers the process of updating the `free_parameters.h` file in NRPy+. This includes generating C code for the parameters set above and writing it to the `free_parameters.h` file.

### Theory Review

#### Introduction to `free_parameters.h` File

The `free_parameters.h` file is used to store the values of various physical constants and parameters in a simulation. In this step, we will update this file with the new values of the domain size, sinh width, sinh v2 constant dr, and SymTP b scale.

```python
# Update $Ccodesrootdir/free_parameters.h based on parameters set above.
outstr = ""
for line in dictID[initial_data_string].freeparams:
    outstr += line + "\n"
```

This code generates the C code for the `free_parameters.h` file by reading from the `dictID[initial_data_string].freeparams` dictionary.

#### Adding Default Free Parameters for Radiative Fluid Module

```python
outstr += rfm.out_default_free_parameters_for_rfm("returnstring",
                                                  domain_size,sinh_width,sinhv2_const_dr,SymTP_bScale)
```

This code adds the default free parameters for the radiative fluid module to the C code.

### Code Implementation

#### Writing C Code to `free_parameters.h` File

```python
with open(os.path.join(Ccodesrootdir,"free_parameters.h"),"w") as file:
    file.write(outstr.replace("params.", "griddata.params."))
```

This code writes the C code to the `free_parameters.h` file.

### Mathematical Background

The domain size, sinh width, sinh v2 constant dr, and SymTP b scale parameters are used to specify the values of various physical constants and parameters in a simulation.

$$\text{Domain Size Parameters} = \left( \begin{array}{c}
domain\_size \\
sinh\_width \\
...
\end{array} \right)$$

where $domain\_size$ and $sinh\_width$ are the domain size parameters.

$$\text{Sinh Width Parameters} = \left( \begin{array}{c}
sinhv2\_const\_dr \\
SymTP\_bScale \\
...
\end{array} \right)$$

where $sinhv2\_const\_dr$ and**Step 4: Validating Black Hole Initial Data**
=============================================

### Overview of the Notebook

This notebook covers the process of validating the black hole initial data against the BSSN Hamiltonian and momentum constraints. This includes checking that the initial data satisfy the constraints using a numerical method.

### Theory Review

#### Introduction to BSSN Hamiltonian and Momentum Constraints

The BSSN (Baumgarte-Shapiro-Simon-Schoen-Numeric) formulation of general relativity is a widely used approach for solving the Einstein field equations numerically. The BSSN Hamiltonian and momentum constraints are a set of partial differential equations that must be satisfied by the initial data.

```python
# Step 4: Validating that the black hole initial data satisfy the BSSN Hamiltonian and momentum constraints.
import numpy as np

def validate_black_hole_initial_data(initial_data):
    # Check if the initial data satisfy the BSSN Hamiltonian constraint
    hamiltonian_constraint = np.linalg.norm(initial_data['A'] - 1) < 1e-6
    
    # Check if the initial data satisfy the BSSN momentum constraints
    momentum_constraints = np.all(np.abs(initial_data['K'] + np.dot(initial_data['B'], initial_data['B'])) < 1e-6)
    
    return hamiltonian_constraint and momentum_constraints

initial_data = dictID[initial_data_string]
is_valid = validate_black_hole_initial_data(initial_data)
```

This code checks if the initial data satisfy the BSSN Hamiltonian and momentum constraints using a numerical method.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import numpy as np
```

This line of code imports the `numpy` module, which is used for numerical computations.

#### Defining Function to Validate Initial Data

```python
def validate_black_hole_initial_data(initial_data):
    # Check if the initial data satisfy the BSSN Hamiltonian constraint
    hamiltonian_constraint = np.linalg.norm(initial_data['A'] - 1) < 1e-6
    
    # Check if the initial data satisfy the BSSN momentum constraints
    momentum_constraints = np.all(np.abs(initial_data['K'] + np.dot(initial_data['B'], initial_data['B'])) < 1e-6)
    
    return hamiltonian_constraint and momentum_constraints
```

This code defines a function `validate_black_hole_initial**Validating Black Hole Initial Data**
=====================================

### Overview of the Notebook

This notebook covers the process of validating the black hole initial data against the Hamiltonian constraint using numerical methods.

### Theory Review

#### Introduction to BSSN Hamiltonian Constraint

The BSSN (Baumgarte-Shapiro-Simon-Schoen-Numeric) formulation of general relativity is a widely used approach for solving the Einstein field equations numerically. The BSSN Hamiltonian constraint is a partial differential equation that must be satisfied by the initial data.

```python
# We will validate that the black hole initial data satisfy the Hamiltonian constraint, modulo numerical finite differencing error.
import BSSN.BSSN_Ccodegen_library as BCl
_ignore = BCl.add_Ricci_eval_to_Cfunction_dict(includes=["NRPy_basic_defines.h"], rel_path_to_Cparams=os.path.join("."),
                                               enable_rfm_precompute=False, enable_golden_kernels=False, enable_SIMD=False,
                                               enable_split_for_optimizations_doesnt_help=False, OMP_pragma_on="i2")

_ignore = BCl.add_BSSN_constraints_to_Cfunction_dict(includes=["NRPy_basic_defines.h"],
                                                     rel_path_to_Cparams=os.path.join("."), output_H_only=False,
                                                     enable_rfm_precompute=False, enable_SIMD=False,
                                                     leave_Ricci_symbolic=True)
```

This code generates the C code for the 3-Ricci tensor and BSSN constraints using numerical methods.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import BSSN.BSSN_Ccodegen_library as BCl
```

This line of code imports the `BSSN.BSSN_Ccodegen_library` module, which is used for generating C code for the 3-Ricci tensor and BSSN constraints.

#### Generating C Code for 3-Ricci Tensor

```python
_ignore = BCl.add_Ricci_eval_to_Cfunction_dict(includes=["NRPy_basic_defines.h"], rel_path_to_Cparams=os.path.join("."),
                                               enable_rfm_precompute=False, enable_golden_kernels=False, enable_SIMD=False,
                                               enable_split_for_optimizations_doesnt_help=False, OMP_pragma_on="i2")
```

This code generates the C code for the 3-Ricci tensor using numerical methods.

#### Generating C Code for BSSN**Step 5: `Initial_Data_Playground.c`**
=====================================

### Overview of the Notebook

This notebook covers the process of generating the main C code for the initial data playground, which is stored in the file `Initial_Data_Playground.c`. This code will be used to run the simulation and validate the initial data.

### Theory Review

#### Introduction to Initial Data Playground

The initial data playground is a C program that generates and validates the initial data for a black hole simulation. The program takes the input parameters, such as the mass of the black hole, the spin of the black hole, and the type of spacetime metric used, and uses them to generate the initial data.

```python
# Step 5: `Initial_Data_Playground.c`: The Main C Code

#include <stdio.h>
#include "NRPy_basic_defines.h"
#include "Ccodesdir/free_parameters.h"

int main() {
    // Read input parameters from command line
    char* initial_data_string = getenv("initial_data_string");

    // Initialize free parameters
    initialize_free_parameters();

    // Load initial data
    load_initial_data(initial_data_string);

    // Run simulation
    run_simulation();
}
```

This code includes the necessary header files, such as `NRPy_basic_defines.h` and `Ccodesdir/free_parameters.h`, which contain functions for initializing free parameters and loading initial data.

### Code Implementation

#### Including Header Files in C Code

```python
#include <stdio.h>
#include "NRPy_basic_defines.h"
#include "Ccodesdir/free_parameters.h"
```

This code includes the necessary header files, which contain functions for input/output operations, initializing free parameters, and loading initial data.

#### Initializing Free Parameters

```c
void initialize_free_parameters() {
    // Initialize free parameters using functions from NRPy+
    init_freedatavalues();
}
```

This function initializes the free parameters using functions from NRPy+.

#### Loading Initial Data

```c
void load_initial_data(char* initial_data_string) {
    // Load initial data using functions from NRPy+
    load_initialdata(initial_data_string);
}
```

This function loads the initial data using functions from NRPy+.

### Mathematical Background

The initial data playground uses numerical methods to generate and validate the initial data for a black hole simulation. The program takes the input parameters, such as the mass of the black hole, the spin of the black hole, and**Step 6: Adding Plane Diagnostics to `Initial_Data_Playground.c`**
=============================================================

### Overview of the Notebook

This notebook covers the process of adding plane diagnostics to the main C code for the initial data playground. This includes generating C code for diagnostic outputs at points closest to the xy plane.

### Theory Review

#### Introduction to Plane Diagnostics

Plane diagnostics are used to output data at specific points in a simulation. In this case, we will add plane diagnostics to output data on the xy plane.

```python
# Step 6: Adding plane diagnostics to `Initial_Data_Playground.c`

list_of_outputs = ["y_n_gfs[IDX4ptS(CFGF,idx)]",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(HGF,idx)]))",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU0GF,idx)])+1e-15)",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU1GF,idx)])+1e-15)",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU2GF,idx)])+1e-15)"]
planar_diags.add_to_Cfunction_dict__plane_diagnostics(plane="xy", include_ghosts=False,
                                                      list_of_outputs=list_of_outputs, num_sig_figs=4)
```

This code generates C code for diagnostic outputs on the xy plane.

### Code Implementation

#### Importing Necessary Modules in Python

```python
import NRPy_basic_defines as basic_defines
import planar_diags
```

This line of code imports the necessary modules, including `NRPy_basic_defines` and `planar_diags`.

#### Adding Plane Diagnostics to C Function Dictionary

```python
list_of_outputs = ["y_n_gfs[IDX4ptS(CFGF,idx)]",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(HGF,idx)]))",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU0GF,idx)])+1e-15)",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU1GF,idx)])+1e-15)",
                   "log10(fabs(diagnostic_output_gfs[IDX4ptS(MU2GF,idx)])+1e-15)"]
planar_diags.add_to_Cfunction**Step 7: Generating Main C Code**
================================

### Overview of the Notebook

This notebook covers the process of generating the main C code for the initial data playground. This includes reading command-line input, setting up grid structure, allocating memory for gridfunctions, and applying boundary conditions.

### Theory Review

#### Introduction to Command-Line Input Reading

The program reads three command-line arguments: `Nx0`, `Nx1`, and `Nx2`. These arguments specify the number of grid points in each direction.

```python
// Step 0.b: Read command-line input, error out if nonconformant
if((argc != 4) || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < 2 /* FIXME; allow for axisymmetric sims */) {
    fprintf(stderr,"Error: Expected three command-line arguments: ./BrillLindquist_Playground Nx0 Nx1 Nx2,\n");
    fprintf(stderr,"where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
    fprintf(stderr,"Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
    exit(1);
}
```

This code checks if the command-line input is valid.

#### Setting Up Grid Structure

The program sets up the grid structure by calling `set_Nxx_dxx_invdx_params__and__xx()` twice, once for the chosen eigen-coordinate system and once for the non-eigen coordinate system.

```c
// Step 0.d: Uniform coordinate grids are stored to *xx[3]
// Step 0.d.i: Set bcstruct
{
    int EigenCoord;
    EigenCoord = 1;
    // Step 0.d.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0.e: Find ghostzone mappings; set up bcstruct
    bcstruct_set_up(&griddata.params, griddata.xx, &griddata.bcstruct);
    // Step 0.e**Step 8: Registering C Functions and NRPy Basic Defines**
=====================================================

### Overview of the Notebook

This notebook covers the process of registering C functions and NRPy basic defines. This includes generating C code for outputting results to a file.

### Theory Review

#### Introduction to Outputting Results to a File

The program uses the `outputC` module to register C functions and NRPy basic defines. This allows the program to output results to a file.

```python
import outputC as outC
outC.outputC_register_C_functions_and_NRPy_basic_defines()
```

This code registers the C functions and NRPy basic defines.

### Code Implementation

#### Registering C Functions and NRPy Basic Defines

```python
import outputC as outC
outC.outputC_register_C_functions_and_NRPy_basic_defines()
```

This code registers the C functions and NRPy basic defines.

### Mathematical Background

The program uses numerical methods to generate and validate initial data for a black hole simulation. The program takes input parameters such as the mass of the black hole, the spin of the black hole, and the type of spacetime metric used, and uses them to generate the initial data.

$$\label{mathematical_background}$$

$$
\text{Initial Data Playground} = f(\text{mass}, \text{spin}, \text{spacetime metric})
$$

This equation represents the initial data playground as a function of input parameters.

### Example Use Cases

*   Running the program with different input parameters to generate and validate different sets of initial data.
*   Using the output file generated by the program to visualize and analyze the results.

```python
# Example usage:
initial_data_string = "BrillLindquist_Playground"
dictID[initial_data_string].functionname = "exact_brill_lindquist_solution"

mainC_function = """// main() function:
    // Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
    // Step 1: Set up initial data to an exact solution
    // Step 2: Output data on xy plane to file.
    // Step 3: Free all allocated memory

    ...
"""

add_to_Cfunction_dict_main__Exact_Initial_Data_Playground()
```

This code generates the C code for the main function using the `add_to_Cfunction_dict` function.**Defining Mathematical Constants**
=====================================

### Overview of the Notebook

This notebook covers the process of defining mathematical constants such as `M_PI`, which is a common constant used in mathematical formulas.

### Theory Review

#### Introduction to Mathematical Constants

Mathematical constants are values that do not change and are often used in mathematical formulas. In this case, we will define `M_PI` as 3.14159265358979323846.

```python
#define M_PI   3.14159265358979323846
```

This code defines `M_PI` as a constant value.

### Code Implementation

#### Defining Mathematical Constants in C

```c
#define M_PI   3.14159265358979323846
#define M_E    2.71828182845904523536
#define M_SQRT_2 1.41421356237309504880
```

This code defines several mathematical constants, including `M_PI`, `M_E`, and `M_SQRT_2`.

#### Using Mathematical Constants in C

```c
#include <stdio.h>

int main() {
    double x = M_PI / 4;
    printf("%f\n", x);
    return 0;
}
```

This code includes the mathematical constants defined earlier and uses them in a simple `main` function.

### Mathematical Background

The mathematical constants defined here are commonly used in mathematical formulas. For example, `M_PI` is often used to calculate the area of a circle:

$$\label{mathematical_background}$$

$$
A = \pi r^2
$$

This equation uses `M_PI` to calculate the area of a circle.

### Example Use Cases

*   Using mathematical constants in numerical computations.
*   Defining custom mathematical constants for specific applications.

```python
// Define custom mathematical constant
#define MY_PI 3.14159265358979323846

int main() {
    double x = MY_PI / 4;
    printf("%f\n", x);
    return 0;
}
```

This code defines a custom mathematical constant `MY_PI` and uses it in a simple `main` function.

### Theory Review

#### Introduction to Mathematical Constants in C

Mathematical constants are values that do not change and are often used in mathematical formulas. In C, these constants can be defined using the `#define` directive.

```c
#define M_PI**Defining Mathematical Constants**
=====================================

### Overview of the Notebook

This notebook covers the process of defining mathematical constants such as `M_PI`, which is a common constant used in mathematical formulas.

### Theory Review

#### Introduction to Mathematical Constants

Mathematical constants are values that do not change and are often used in mathematical formulas. In this case, we will define `M_PI` as 3.14159265358979323846.

```c
#define M_PI   3.14159265358979323846
```

This code defines `M_PI` as a constant value.

### Code Implementation

#### Defining Mathematical Constants in C

```c
#define M_PI   3.14159265358979323846
#define M_E    2.71828182845904523536
#define M_SQRT_2 1.41421356237309504880
```

This code defines several mathematical constants, including `M_PI`, `M_E`, and `M_SQRT_2`.

#### Using Mathematical Constants in C

```c
#include <stdio.h>

int main() {
    double x = M_PI / 4;
    printf("%f\n", x);
    return 0;
}
```

This code includes the mathematical constants defined earlier and uses them in a simple `main` function.

### Mathematical Background

The mathematical constants defined here are commonly used in mathematical formulas. For example, `M_PI` is often used to calculate the area of a circle:

$$\label{mathematical_background}$$

$$
A = \pi r^2
$$

This equation uses `M_PI` to calculate the area of a circle.

### Example Use Cases

*   Using mathematical constants in numerical computations.
*   Defining custom mathematical constants for specific applications.

```python
# Define custom mathematical constant
def define_custom_constant(name, value):
    return f"#define {name} {value}"

custom_constant = define_custom_constant("MY_PI", 3.14159265358979323846)
print(custom_constant)
```

This code defines a function to create custom mathematical constants and uses it to define `MY_PI`.

### Theory Review

#### Introduction to Mathematical Constants in C

Mathematical constants are values that do not change and are often used in mathematical formulas. In C, these constants can be defined using the `#define` directive.

```c
#define M_PI  **Setting up NRPy_basic_defines.h**
=====================================

### Overview of the Notebook

This notebook covers the process of setting up `NRPy_basic_defines.h`, which is a header file that defines various constants and data structures used in the numerical relativity code.

### Theory Review

#### Introduction to NRPy_basic_defines.h

`NRPy_basic_defines.h` is a header file that defines various constants and data structures used in the numerical relativity code. It is generated by the `outC` module, which is responsible for creating C code from Python modules.

```c
#include <stdio.h>

// Define NRPy_basic_defines.h
#define NRPY_BASIC_DEFINES_H

// Define gridfunction structure
typedef struct __gridfunctions_struct__ {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} gridfunctions_struct;

// Define MoL gridfunction structure
typedef struct __MoL_gridfunctions_struct__ {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;
```

This code defines the `gridfunctions_struct` and `MoL_gridfunctions_struct` structures, which are used to store gridfunction data.

### Code Implementation

#### Setting up NRPy_basic_defines.h using outC module

```python
import outC as outC

outC.outC_NRPy_basic_defines_h_dict["MoL"] = """
typedef struct __MoL_gridfunctions_struct__ {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;
"""
par.register_NRPy_basic_defines()
```

This code sets up the `NRPy_basic_defines.h` header file using the `outC` module. The `MoL` gridfunction structure is defined, and the `register_NRPy_basic_defines()` function is called to register the NRPy basic defines.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The `NRPy_basic_defines.h` header file provides a way to define these constants and structures in a platform-independent manner.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot \h**Modifying the Grid Data Structure**
=====================================

### Overview of the Notebook

This notebook covers the process of modifying the grid data structure by adding extra variables and structures.

### Theory Review

#### Introduction to Grid Data Structure

The grid data structure is used to store information about the grid, such as its dimensions, spacing, and coordinates. In this case, we will add an extra variable and structure to the grid data structure.

```python
struct griddata {
  // ... existing members ...
  MoL_gridfunctions_struct gridfuncs;
};
```

This code adds a `MoL_gridfunctions_struct` member to the `griddata` struct.

### Code Implementation

#### Modifying the Grid Data Structure in C

```c
#include <stdio.h>

// Define grid data structure
typedef struct {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;

struct griddata {
  // ... existing members ...
  MoL_gridfunctions_struct gridfuncs;
};
```

This code defines the `MoL_gridfunctions_struct` and modifies the `griddata` struct to include a member of this type.

#### Modifying the Grid Data Structure in Python

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Add extra variables and structures to the grid data structure.""",
    c_type="void",
    name="modify_grid_data_structure",
    params="paramstruct params",
    body="""
        // ... existing code ...
        MoL_gridfunctions_struct gridfuncs;
        griddata.gridfuncs = gridfuncs;
    """,
    rel_path_to_Cparams=os.path.join("."),
    enableCparameters=True
)
```

This code uses the `add_to_Cfunction_dict` function to modify the grid data structure by adding a member of type `MoL_gridfunctions_struct`.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The grid data structure provides a way to store information about the grid.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G**Modifying the `griddata` Struct**
=====================================

### Overview of the Notebook

This notebook covers the process of modifying the `griddata` struct to include extra variables and structures.

### Theory Review

#### Introduction to Modifying Data Structures

In C, data structures can be modified by adding new members or changing existing ones. This is done using the `typedef` keyword.

```c
typedef struct {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;

struct griddata {
  // ... existing members ...
  MoL_gridfunctions_struct gridfuncs;
};
```

This code defines a new data structure `MoL_gridfunctions_struct` and adds it to the `griddata` struct.

### Code Implementation

#### Modifying the `griddata` Struct in C

```c
#include <stdio.h>

typedef struct {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;

struct griddata {
  // ... existing members ...
  MoL_gridfunctions_struct gridfuncs;
};
```

This code defines the `MoL_gridfunctions_struct` and modifies the `griddata` struct.

#### Modifying the `griddata` Struct in Python

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Add extra variables and structures to the `griddata` struct.""",
    c_type="void",
    name="modify_grid_data_structure",
    params="paramstruct params",
    body="""
        // ... existing code ...
        MoL_gridfunctions_struct gridfuncs;
        griddata.gridfuncs = gridfuncs;
    """,
    rel_path_to_Cparams=os.path.join("."),
    enableCparameters=True
)
```

This code uses the `add_to_Cfunction_dict` function to modify the `griddata` struct.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The `griddata` struct provides a way to store information about the grid.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the**Defining `IDX3S()` and Other Functions**
=====================================

### Overview of the Notebook

This notebook covers the process of defining various functions used in the numerical relativity code, including `IDX3S()`, `GRIDINDEX()`, and others.

### Theory Review

#### Introduction to Indexing Functions

The indexing functions are used to calculate indices for grid points and other data structures. In this case, we will define `IDX3S()` and `GRIDINDEX()`.

```c
#define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])
```

This code defines the `IDX3S()` function, which calculates an index for a grid point.

### Code Implementation

#### Defining Indexing Functions in C

```c
#include <stdio.h>

#define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])

#define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
```

This code defines the `IDX3S()` and `GRIDINDEX()` functions.

#### Defining Indexing Functions in Python

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Define indexing functions used in the numerical relativity code.""",
    c_type="void",
    name="define_indexing_functions",
    params="paramstruct params",
    body="""
        // ... existing code ...
        #define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])
        #define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
    """,
    rel_path_to_Cparams=os.path.join("."),
    enableCparameters=True
)
```

This code uses the `add_to_Cfunction_dict` function to define the indexing functions.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The indexing functions provide a way to efficiently access grid points and other data structures.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can**Registering C Functions and NRPy Basic Defines**
=====================================================

### Overview of the Notebook

This notebook covers the process of registering C functions and NRPy basic defines.

### Theory Review

#### Introduction to Registering C Functions and NRPy Basic Defines

In this step, we will register the C functions and NRPy basic defines that have been defined in previous steps. This is done using the `register_C_functions_and_NRPy_basic_defines()` function.

```python
gri.register_C_functions_and_NRPy_basic_defines(list_of_extras_in_griddata_struct=list_of_extras_in_griddata_struct)
```

This code registers the C functions and NRPy basic defines.

### Code Implementation

#### Registering C Functions and NRPy Basic Defines in Python

```python
import gri as gri

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

gri.register_C_functions_and_NRPy_basic_defines(list_of_extras_in_griddata_struct=list_of_extras_in_griddata_struct)
```

This code registers the C functions and NRPy basic defines.

#### Registering C Functions and NRPy Basic Defines in C

```c
#include <stdio.h>

// ... existing code ...

void register_C_functions_and_NRPy_basic_defines(void) {
  // ... existing code ...
}
```

This code registers the C functions and NRPy basic defines.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The registration of C functions and NRPy basic defines is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Registering C functions and NRPy basic defines for a specific problem or simulation.
*   Modifying the registration process to accommodate different problem types or numerical methods.

```python
# ... existing code ...

gri.register_C_functions_and_NRPy_basic_defines(list_of_extras_in_griddata_struct=list_of_extras_in_griddata_struct)

# ... existing code ...
```

This code registers the C functions and NRPy basic defines.**Defining Indexing Functions**
=============================

### Overview of the Notebook

This notebook covers the process of defining various indexing functions used in the numerical relativity code.

### Theory Review

#### Introduction to Indexing Functions

Indexing functions are used to calculate indices for grid points and other data structures. In this case, we will define `IDX3S()`, `GRIDINDEX()`, and others.

```c
#define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])
```

This code defines the `IDX3S()` function, which calculates an index for a grid point.

### Code Implementation

#### Defining Indexing Functions in C

```c
#include <stdio.h>

#define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])

#define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
```

This code defines the `IDX3S()` and `GRIDINDEX()` functions.

#### Defining Indexing Functions in Python

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Define indexing functions used in the numerical relativity code.""",
    c_type="void",
    name="define_indexing_functions",
    params="paramstruct params",
    body="""
        // ... existing code ...
        #define IDX3S(zzz,i,j,k) ((zzz)[(i)*NXX[k] + (j)])
        #define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
    """,
    rel_path_to_Cparams=os.path.join("."),
    enableCparameters=True
)
```

This code uses the `add_to_Cfunction_dict` function to define the indexing functions.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The indexing functions provide a way to efficiently access grid points and other data structures.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P**Defining `GRIDINDEX()`**
=========================

### Overview of the Notebook

This notebook covers the process of defining the `GRIDINDEX()` function.

### Theory Review

#### Introduction to Grid Indexing

Grid indexing is a way to access grid points in numerical relativity simulations. The `GRIDINDEX()` function is used to calculate the index of a grid point.

```c
#define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
```

This code defines the `GRIDINDEX()` function, which calculates the index of a grid point.

### Code Implementation

#### Defining `GRIDINDEX()` in C

```c
#include <stdio.h>

#define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))

int main() {
    int i = 1;
    int j = 2;
    int k = 3;

    int index = GRIDINDEX(i, j, k);
    printf("%d\n", index);

    return 0;
}
```

This code defines the `GRIDINDEX()` function and uses it to calculate the index of a grid point.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The `GRIDINDEX()` function provides a way to efficiently access grid points.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Using `GRIDINDEX()` to access grid points in numerical relativity simulations.
*   Modifying the `GRIDINDEX()` function to accommodate different grid structures or indexing schemes.

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Define the `GRIDINDEX()` function used in numerical relativity simulations.""",
    c_type="void",
    name="define_GRIDINDEX",
    params="paramstruct params",
    body="""
        // ... existing code ...
        #define GRIDINDEX(i,j,k) (((i)*NXX[k] + (j)))
    """,
    rel_path_to_C**Declaring `paramstruct` and Registering C Parameters**
=====================================================

### Overview of the Notebook

This notebook covers the process of declaring a `paramstruct` and registering C parameters using the `set_Cparameters_to_default()` function.

### Theory Review

#### Introduction to Parametrized Code

The numerical relativity code is parametrized, meaning that certain variables are defined in terms of other variables. This allows for greater flexibility and customization of the code.

```c
typedef struct {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;

struct paramstruct {
  // ... existing members ...
};
```

This code defines a `paramstruct` that contains various members, including grid function structures and other variables.

### Code Implementation

#### Declaring `paramstruct` in C

```c
#include <stdio.h>

typedef struct {
  REAL *restrict y_n_gfs;
  REAL *restrict auxevol_gfs;
} MoL_gridfunctions_struct;

struct paramstruct {
  // ... existing members ...
};
```

This code defines the `paramstruct`.

#### Registering C Parameters using `set_Cparameters_to_default()` in Python

```python
import NRPy_basic_defines as basic_defines

list_of_extras_in_griddata_struct = ["MoL_gridfunctions_struct gridfuncs;"]

basic_defines.add_to_Cfunction_dict(
    includes=["NRPy_basic_defines.h"],
    desc="""Declare a `paramstruct` and register C parameters using the `set_Cparameters_to_default()` function.""",
    c_type="void",
    name="declare_paramstruct_and_register_C_parameters",
    params="paramstruct params",
    body="""
        // ... existing code ...
        struct paramstruct {
          // ... existing members ...
        };
        set_Cparameters_to_default();
    """,
    rel_path_to_Cparams=os.path.join("."),
    enableCparameters=True
)
```

This code uses the `add_to_Cfunction_dict` function to declare a `paramstruct` and register C parameters using the `set_Cparameters_to_default()` function.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The parametrized code allows for greater flexibility and customization of the code.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Plan**Registering C Functions and NRPy Basic Defines**
=====================================================

### Overview of the Notebook

This notebook covers the process of registering C functions and NRPy basic defines.

### Theory Review

#### Introduction to Registering C Functions and NRPy Basic Defines

In this step, we will register the C functions and NRPy basic defines using the `outC.NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines()` function.

```python
import outC as outC

outC.NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(os.path.join(Ccodesrootdir))
```

This code registers the C functions and NRPy basic defines.

### Code Implementation

#### Registering C Functions and NRPy Basic Defines in Python

```python
import os
import outC as outC

# Set the path to the NRPy root directory
NRPy_root_dir = os.path.join(Ccodesrootdir)

# Register C functions and NRPy basic defines
outC.NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(NRPy_root_dir)
```

This code registers the C functions and NRPy basic defines.

#### Registering C Functions and NRPy Basic Defines in C

```c
#include <stdio.h>

// ... existing code ...

void register_C_functions_and_NRPy_basic_defines(void) {
  // ... existing code ...
}
```

This code registers the C functions and NRPy basic defines.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The registration of C functions and NRPy basic defines is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Registering C functions and NRPy basic defines for a specific problem or simulation.
*   Modifying the registration process to accommodate different problem types or numerical methods.

```python
import os
import outC as outC

# Set the path to the NRPy root directory
NRPy_root_dir = os.path.join(Ccodesrootdir)

# Register C functions and NRPy basic defines
outC.NRPy_param_funcs_register_C_functions**Registering C Functions and NRPy Basic Defines**
=====================================================

### Overview of the Notebook

This notebook covers the process of registering C functions and NRPy basic defines.

### Theory Review

#### Introduction to Registering C Functions and NRPy Basic Defines

In this step, we will register the C functions and NRPy basic defines using the `fin.register_C_functions_and_NRPy_basic_defines()` function.

```python
import fin as fin

fin.register_C_functions_and_NRPy_basic_defines(NGHOSTS_account_for_onezone_upwind=False,
                                                enable_SIMD=False)
```

This code registers the C functions and NRPy basic defines.

### Code Implementation

#### Registering C Functions and NRPy Basic Defines in Python

```python
import os
import fin as fin

# Set the path to the NRPy root directory
NRPy_root_dir = os.path.join(Ccodesrootdir)

# Register C functions and NRPy basic defines
fin.register_C_functions_and_NRPy_basic_defines(NGHOSTS_account_for_onezone_upwind=False,
                                                 enable_SIMD=False)
```

This code registers the C functions and NRPy basic defines.

#### Registering C Functions and NRPy Basic Defines in C

```c
#include <stdio.h>

// ... existing code ...

void register_C_functions_and_NRPy_basic_defines(void) {
  // ... existing code ...
}
```

This code registers the C functions and NRPy basic defines.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The registration of C functions and NRPy basic defines is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Registering C functions and NRPy basic defines for a specific problem or simulation.
*   Modifying the registration process to accommodate different problem types or numerical methods.

```python
import os
import fin as fin

# Set the path to the NRPy root directory
NRPy_root_dir = os.path.join(Ccodesrootdir)

# Register C functions and NRPy basic defines
fin.register_C_functions_and_NRPy_basic_defines(NG**Defining `NGHOSTS` and the `UPWIND()` Macro**
=====================================================

### Overview of the Notebook

This notebook covers the process of defining `NGHOSTS` and the `UPWIND()` macro.

### Theory Review

#### Introduction to Grid Data Structures

In numerical relativity, grid data structures are used to store information about the grid. `NGHOSTS` is a parameter that defines the number of ghost zones in the grid.

```c
#define NGHOSTS 1
```

This code defines `NGHOSTS` as 1.

#### Introduction to the `UPWIND()` Macro

The `UPWIND()` macro is used to calculate the upwind direction of a quantity. If SIMD (Single Instruction, Multiple Data) is disabled, we need to define the `UPWIND()` macro.

```c
#define UPWIND(zzz,i,j,k,l) ((zzz)[(i)*NXX[k] + (j)])
```

This code defines the `UPWIND()` macro.

### Code Implementation

#### Defining `NGHOSTS` and the `UPWIND()` Macro in C

```c
#include <stdio.h>

#define NGHOSTS 1

#define UPWIND(zzz,i,j,k,l) ((zzz)[(i)*NXX[k] + (j)])

int main() {
    int i = 1;
    int j = 2;
    int k = 3;

    int upwind = UPWIND(griddata, i, j, k, 0);
    printf("%d\n", upwind);

    return 0;
}
```

This code defines `NGHOSTS` and the `UPWIND()` macro.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The definition of `NGHOSTS` and the `UPWIND()` macro is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Defining `NGHOSTS` for a specific problem or simulation.
*   Modifying the definition of**Outputting Finite-Difference Stencil Functions**
=====================================================

### Overview of the Notebook

This notebook covers the process of outputting functions for computing all finite-difference stencils.

### Theory Review

#### Introduction to Finite-Difference Methods

Finite-difference methods are used in numerical relativity to discretize spacetime and approximate derivatives. The stencil functions are a crucial part of this process, as they define how the grid values are combined to compute the derivatives.

```c
#include <stdio.h>

void output_stencil_functions(void) {
  // ... existing code ...
}
```

This code defines the `output_stencil_functions()` function, which outputs the finite-difference stencil functions.

### Code Implementation

#### Outputting Finite-Difference Stencil Functions in C

```c
#include <stdio.h>

#define NXX 3
#define NYY 3

void output_stencil_functions(void) {
  // Define the grid values and derivatives
  REAL *grid_data = (REAL *)malloc(NXX*NYY*sizeof(REAL));
  REAL *derivative_x = (REAL *)malloc(NXX*NYY*sizeof(REAL));
  REAL *derivative_y = (REAL *)malloc(NXX*NYY*sizeof(REAL));

  // Compute the grid values and derivatives using finite-difference methods
  for (int i = 0; i < NXX; i++) {
    for (int j = 0; j < NYY; j++) {
      grid_data[i*NYY + j] = /* compute grid value */;
      derivative_x[i*NYY + j] = /* compute x-derivative */;
      derivative_y[i*NYY + j] = /* compute y-derivative */;
    }
  }

  // Output the stencil functions
  printf("Stencil function for x-derivative: %f\n", derivative_x[0]);
  printf("Stencil function for y-derivative: %f\n", derivative_y[0]);

  free(grid_data);
  free(derivative_x);
  free(derivative_y);
}
```

This code outputs the finite-difference stencil functions by computing the grid values and derivatives using finite-difference methods.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The output of the finite-difference stencil functions is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant,**Outputting Finite-Difference Functions**
=========================================

### Overview of the Notebook

This notebook covers the process of outputting finite-difference functions.

### Theory Review

#### Introduction to Finite-Difference Methods

Finite-difference methods are used in numerical relativity to discretize spacetime and approximate derivatives. The stencil functions are a crucial part of this process, as they define how the grid values are combined to compute the derivatives.

```c
#include <stdio.h>

void output_finite_difference_functions_h(void) {
  // ... existing code ...
}
```

This code defines the `output_finite_difference_functions_h()` function, which outputs the finite-difference functions.

### Code Implementation

#### Outputting Finite-Difference Functions in C

```python
import fin as fin

def output_finite_difference_functions(path):
    # Define all functions depending on FD stencils
    # ...

    if enable_FD_functions:
        fin.output_finite_difference_functions_h(path=Ccodesrootdir)
```

This code outputs the finite-difference functions by calling the `output_finite_difference_functions_h()` function.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The output of the finite-difference functions is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Outputting finite-difference functions for a specific problem or simulation.
*   Modifying the output process to accommodate different problem types or numerical methods.

```python
import fin as fin

def output_finite_difference_functions(path):
    # Define all functions depending on FD stencils
    # ...

    if enable_FD_functions:
        fin.output_finite_difference_functions_h(path=Ccodesrootdir)

output_finite_difference_functions(Ccodesrootdir)
```

This code outputs the finite-difference functions by calling the `output_finite_difference_functions()` function.**Setting up NRPy Basic Defines and Function Prototypes**
=====================================================

### Overview of the Notebook

This notebook covers the process of setting up `NRPy_basic_defines.h` and `NRPy_function_prototypes.h`.

### Theory Review

#### Introduction to NRPy Basic Defines

`NRPy_basic_defines.h` is a header file that contains basic definitions for NRPy, such as mathematical constants and functions.

```c
#include <stdio.h>

void construct_NRPy_basic_defines_h(void) {
  // ... existing code ...
}
```

This code defines the `construct_NRPy_basic_defines_h()` function, which constructs `NRPy_basic_defines.h`.

#### Introduction to NRPy Function Prototypes

`NRPy_function_prototypes.h` is a header file that contains function prototypes for NRPy functions.

```c
#include <stdio.h>

void construct_NRPy_function_prototypes_h(void) {
  // ... existing code ...
}
```

This code defines the `construct_NRPy_function_prototypes_h()` function, which constructs `NRPy_function_prototypes.h`.

### Code Implementation

#### Constructing NRPy Basic Defines and Function Prototypes in C

```python
import outC as outC

def construct_NRPy_basic_defines_h(Ccodesrootdir, enable_SIMD=False):
    # ... existing code ...
    outC.construct_NRPy_basic_defines_h(Ccodesrootdir, enable_SIMD=enable_SIMD)

def construct_NRPy_function_prototypes_h(Ccodesrootdir):
    # ... existing code ...
    outC.construct_NRPy_function_prototypes_h(Ccodesrootdir)
```

This code constructs `NRPy_basic_defines.h` and `NRPy_function_prototypes.h`.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The construction of `NRPy_basic_defines.h` and `NRPy_function_prototypes.h` is an important step in preparing the code for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Constructing `NRPy_basic_defines.h` and `NRPy_function_prototypes.h` for a specific problem or**Compiling C Code with Command Line Helper**
=============================================

### Overview of the Notebook

This notebook covers the process of compiling C code using a command line helper.

### Theory Review

#### Introduction to Compilers and Compilation Options

A compiler is a program that translates source code into machine code. The compilation options determine how the compiler processes the source code.

```c
#include <stdio.h>

void new_C_compile(void) {
  // ... existing code ...
}
```

This code defines the `new_C_compile()` function, which compiles C code using a command line helper.

#### Introduction to Compilation Options

Compilation options are flags that control how the compiler processes the source code. Some common compilation options include:

*   `-O2` or `-O3`: Optimizes the code for performance
*   `-g`: Generates debugging information
*   `-Wall`: Enables all warnings

```c
#include <stdio.h>

void new_C_compile(void) {
  // ... existing code ...
  compiler_opt_option = "-O2";
}
```

This code sets the compilation option to optimize the code for performance.

### Code Implementation

#### Compiling C Code using Command Line Helper in C

```python
import cmdline_helper as cmd

def new_C_compile(Ccodesrootdir, uses_free_parameters_h=True, compiler_opt_option="fast"):
    # ... existing code ...
    cmd.new_C_compile(Ccodesrootdir, "Exact_Initial_Data_Playground",
                      uses_free_parameters_h=uses_free_parameters_h,
                      compiler_opt_option=compiler_opt_option)
```

This code compiles C code using a command line helper.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The compilation of the code is an important step in preparing it for execution.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Compiling C code with optimization for performance
*   Compiling C code with debugging information

```python
import cmdline_helper as cmd

def new_C_compile(Ccodesrootdir, uses_free_parameters_h=True, compiler_opt_option="fast"):
    # ... existing code ...
    cmd.new_C_compile(Ccodes**Executing and Plotting Results**
==================================

### Overview of the Notebook

This notebook covers the process of executing and plotting results.

### Theory Review

#### Introduction to Execution and Plotting

The execution and plotting steps are crucial parts of the numerical relativity code. They involve running the simulation and generating plots for visualization.

```c
#include <stdio.h>

void Execute(void) {
  // ... existing code ...
}
```

This code defines the `Execute()` function, which executes the simulation.

### Code Implementation

#### Changing to Output Directory in Python

```python
import os
import cmdline_helper as cmd
from glob import glob

# Change to output directory
os.chdir(Ccodesrootdir)

# Delete existing files
cmd.delete_existing_files("out*.txt")
cmd.delete_existing_files("out*.png")

# Define arguments and outputs
args_output_list = [["96 96 96", "out96.txt"], ["48 48 48", "out48.txt"]]

# Execute the code for each set of arguments
for args_output in args_output_list:
    cmd.Execute("Exact_Initial_Data_Playground", args_output[0], args_output[1])
```

This code changes to the output directory and executes the simulation.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The execution and plotting steps are important for analyzing the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Executing and plotting results for a specific simulation or problem.
*   Modifying the execution and plotting steps to accommodate different simulations or problems.

```python
import os
import cmdline_helper as cmd
from glob import glob

# Change to output directory
os.chdir(Ccodesrootdir)

# Delete existing files
cmd.delete_existing_files("out*.txt")
cmd.delete_existing_files("out*.png")

# Define arguments and outputs
args_output_list = [["96 96 96", "out96.txt"], ["48 48 48", "out48.txt"]]

# Execute the code for each set of arguments
for args_output in args_output_list:
    cmd.Execute("**Plotting Initial Data**
=========================

### Overview of the Notebook

This notebook covers the process of plotting the initial data.

### Theory Review

#### Introduction to Ploting Initial Data

The initial data is a crucial part of any numerical relativity simulation. It represents the initial conditions for the simulation, and it must be carefully chosen to ensure accurate results.

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_initial_data(data):
  # ... existing code ...
```

This code defines a function `plot_initial_data` that takes in the initial data and plots it.

### Code Implementation

#### Plotting Initial Data using Matplotlib

```python
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

# Get the initial data files
initial_data_files = glob("out*.txt")

# Loop over each file and plot the data
for file in initial_data_files:
  # Read in the data from the file
  data = np.loadtxt(file)
  
  # Plot the data
  plt.plot(data[:,0], data[:,1])
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Initial Data')
  plt.show()
```

This code plots the initial data using matplotlib.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The plotting of the initial data is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Plotting initial data for a specific simulation or problem.
*   Modifying the plotting step to accommodate different simulations or problems.

```python
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

# Get the initial data files
initial_data_files = glob("out*.txt")

# Loop over each file and plot the data
for file in initial_data_files:
  # Read in the data from the file
  data = np.loadtxt(file)
  
  # Plot the data
  plt.plot(data[:,0], data[:,1])
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Initial Data')
  plt**Plotting the Evolved Conformal Factor**
======================================

### Overview of the Notebook

This notebook covers the process of plotting the evolved conformal factor.

### Theory Review

#### Introduction to Conformal Factor and Gravitational Fields

The conformal factor is a crucial quantity in numerical relativity, as it encodes information about the gravitational field. By plotting the evolved conformal factor on a 2D grid, we can visualize the strength of the gravitational fields.

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_conformal_factor(data):
  # ... existing code ...
```

This code defines a function `plot_conformal_factor` that takes in the conformal factor data and plots it.

### Code Implementation

#### Plotting Conformal Factor using Matplotlib

```python
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

# Get the conformal factor files
conformal_factor_files = glob("out*.txt")

# Loop over each file and plot the data
for file in conformal_factor_files:
  # Read in the data from the file
  data = np.loadtxt(file)
  
  # Plot the data
  plt.imshow(data, cmap='plasma')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Conformal Factor')
  plt.show()
```

This code plots the conformal factor using matplotlib.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The plotting of the conformal factor is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Plotting conformal factor for a specific simulation or problem.
*   Modifying the plotting step to accommodate different simulations or problems.

```python
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

# Get the conformal factor files
conformal_factor_files = glob("out*.txt")

# Loop over each file and plot the data
for file in conformal_factor_files:
  # Read in the data from the file
 **Installing and Importing Required Libraries**
=============================================

### Overview of the Notebook

This notebook covers the process of installing and importing required libraries for numerical computations.

### Theory Review

#### Introduction to Scientific Computing with SciPy

SciPy is a powerful library for scientific computing in Python. It provides functions for tasks such as optimization, linear algebra, integration, and interpolation. In this notebook, we will use SciPy to perform interpolation on a grid of points.

```python
import scipy as sp
```

This code imports the SciPy library.

### Code Implementation

#### Installing SciPy if Not Already Installed

```bash
!pip install scipy
```

This code installs SciPy using pip if it is not already installed. If SciPy is already installed, this command will have no effect.

```python
Requirement already satisfied: scipy in /home/zetienne/jup311/lib/python3.11/site-packages (1.9.3)
Requirement already satisfied: numpy<1.26.0,>=1.18.5 in /home/zetienne/jup311/lib/python3.11/site-packages (from scipy) (1.24.0)
```

This output indicates that SciPy is already installed.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The interpolation of functions on a grid of points is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Installing SciPy for a specific simulation or problem.
*   Importing required libraries and modules for numerical computations.

```python
import numpy as np
from scipy.interpolate import griddata
from pylab import savefig
import matplotlib.pyplot as plt
from IPython.display import Image
```

This code imports the necessary libraries and modules for numerical computations.

### Interpolation on a Grid of Points

```python
xy_extent=domain_size
```

This code sets the extent of the x and y coordinates for interpolation.

The rest of the code is not shown in this snippet, but it would involve using the `griddata` function from SciPy to interpolate functions on a grid of**Plotting Initial Data on a 2D Grid**
=====================================

### Overview of the Notebook

This notebook covers the process of plotting initial data on a 2D grid.

### Theory Review

#### Introduction to Plotting Initial Data

The initial data is a crucial part of any numerical relativity simulation. It represents the initial conditions for the simulation, and it must be carefully chosen to ensure accurate results.

```python
import matplotlib.pyplot as plt
```

This code imports the matplotlib library.

### Code Implementation

#### Generating Uniform 2D Grids from Data

```python
data = np.loadtxt('out96.txt')
output_grid_data = []
for i in [3, 4, 5, 6, 7]:
    output_grid_x, output_grid_y, output_grid_data_i = \
        plot2D.generate_uniform_2D_grid(data[:,0], data[:,1], i, [-xy_extent,xy_extent], [-xy_extent,xy_extent])
    output_grid_data += [output_grid_data_i]
```

This code generates uniform 2D grids from the initial data.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The plotting of initial data on a 2D grid is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Plotting initial data on a 2D grid for a specific simulation or problem.
*   Modifying the plotting step to accommodate different simulations or problems.

```python
plt.clf()
plt.title(r"Initial Data, conformal factor $W$")
plt.xlabel(r"$x/M$")
plt.ylabel(r"$y/M$")

plt.imshow(output_grid_data[0], extent=(-xy_extent,xy_extent, -xy_extent,xy_extent))
savefig("ID.png")
plt.close()
Image("ID.png")
```

This code plots the initial data on a 2D grid using matplotlib.

### Output

The resulting plot will be saved as "ID.png".

![png](output_33_0.png)


<a id='convergence'></a>

**Convergence of Initial Data**
**Convergence of Numerical Errors**
====================================

### Overview of the Notebook

This notebook covers the process of validating the convergence of numerical errors.

### Theory Review

#### Introduction to Convergence and Error Analysis

The convergence of numerical errors is a crucial aspect of numerical relativity. It ensures that the simulation produces accurate results by minimizing the impact of numerical approximations.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Checking Hamiltonian and Momentum Constraint Violations

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract Hamiltonian and momentum constraint violations
hamiltonian_violation = data[:, 3]
momentum_constraint_violation = data[:, 4]

# Plot convergence of numerical errors
plt.plot(hamiltonian_violation, label='Hamiltonian Violation')
plt.plot(momentum_constraint_violation, label='Momentum Constraint Violation')
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Convergence of Numerical Errors')
plt.legend()
savefig('convergence.png')
```

This code checks the Hamiltonian and momentum constraint violations and plots their convergence.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The convergence of numerical errors is an important aspect of ensuring accurate results in the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Validating convergence of numerical errors for a specific simulation or problem.
*   Modifying the validation step to accommodate different simulations or problems.

```python
import numpy as np

# Load data from output file
data = np.loadtxt('output.txt')

# Extract Hamiltonian and momentum constraint violations
hamiltonian_violation = data[:, 3]
momentum_constraint_violation = data[:, 4]

# Plot convergence of numerical errors
plt.plot(hamiltonian_violation, label='Hamiltonian Violation')
plt.plot(momentum_constraint_violation, label='Momentum Constraint Violation')
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Convergence of Numerical Errors')
plt**Convergence of Numerical Errors**
====================================

### Overview of the Notebook

This notebook covers the process of validating the convergence of numerical errors.

### Theory Review

#### Introduction to Convergence and Error Analysis

The convergence of numerical errors is a crucial aspect of numerical relativity. It ensures that the simulation produces accurate results by minimizing the impact of numerical approximations.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Plotting Convergence of Numerical Errors

```python
plt.clf()

# Load data from output file
data = np.loadtxt('output.txt')

# Extract Hamiltonian constraint violation
hamiltonian_violation = data[:, 3]

# Plot convergence of numerical errors
plt.imshow(np.log10(np.abs(hamiltonian_violation)), cmap='plasma')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Convergence of Numerical Errors')
savefig('convergence.png')
```

This code plots the convergence of numerical errors using matplotlib.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The convergence of numerical errors is an important aspect of ensuring accurate results in the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Validating convergence of numerical errors for a specific simulation or problem.
*   Modifying the validation step to accommodate different simulations or problems.

```python
import numpy as np

# Load data from output file
data = np.loadtxt('output.txt')

# Extract Hamiltonian constraint violation
hamiltonian_violation = data[:, 3]

# Plot convergence of numerical errors
plt.imshow(np.log10(np.abs(hamiltonian_violation)), cmap='plasma')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Convergence of Numerical Errors')
savefig('convergence.png')
```

This code plots the convergence of numerical errors.

### Theory Review

#### Convergence Rate

The convergence rate is determined by the error term that drops with the uniform gridspacing to the fourth power: $\left(\Delta x^i\right)^4$.

**Plotting Hamiltonian and Momentum Constraints**
=============================================

### Overview of the Notebook

This notebook covers the process of plotting Hamiltonian and momentum constraints.

### Theory Review

#### Introduction to Hamiltonian and Momentum Constraints

Hamiltonian constraint is a fundamental concept in numerical relativity, ensuring that the simulation produces accurate results by minimizing the impact of numerical approximations. Similarly, momentum constraint is another crucial aspect of numerical relativity, which helps in understanding the behavior of the system under study.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Plotting Hamiltonian Constraint

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract Hamiltonian constraint violation
hamiltonian_violation = data[:, 3]

# Create plot for Hamiltonian constraint
plt.plot(hamiltonian_violation)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Hamiltonian Constraint Violation')
savefig('Hamiltonian_Constraint.png')
```

This code plots the Hamiltonian constraint violation.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The plotting of Hamiltonian and momentum constraints is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant. Then, we can define the following mathematical constants:

$$
M_G = G \cdot M_P
$$

where $M_P$ is the Planck mass.

### Example Use Cases

*   Plotting Hamiltonian constraint for a specific simulation or problem.
*   Modifying the plotting step to accommodate different simulations or problems.

#### Plotting Momentum Constraints

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract momentum constraint violation
momentum_violation = data[:, 4:7]

# Create plot for momentum constraints
plt.plot(momentum_violation)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Momentum Constraint Violation')
savefig('Momentum_Constraint.png')
```

This code plots the momentum constraint violation.

### Theory Review

#### Convergence Rate

The convergence rate is determined by the error term that drops with the uniform gridspacing to the fourth power: $\left(\Delta x^i\right)^4$.**Plotting Momentum Constraints**
=================================

### Overview of the Notebook

This notebook covers the process of plotting momentum constraints.

### Theory Review

#### Introduction to Momentum Constraints

Momentum constraint is a fundamental concept in numerical relativity, ensuring that the simulation produces accurate results by minimizing the impact of numerical approximations.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Plotting Momentum Constraint in Radial Direction (r)

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract momentum constraint violation in radial direction
momentum_violation_r = data[:, 4]

# Create plot for momentum constraint in radial direction
plt.plot(momentum_violation_r)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Momentum Constraint Violation (r)')
savefig('Momentum_Constraint_r.png')
```

This code plots the momentum constraint violation in the radial direction.

#### Plotting Momentum Constraint in Polar Direction (th)

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract momentum constraint violation in polar direction
momentum_violation_th = data[:, 5]

# Create plot for momentum constraint in polar direction
plt.plot(momentum_violation_th)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Momentum Constraint Violation (th)')
savefig('Momentum_Constraint_th.png')
```

This code plots the momentum constraint violation in the polar direction.

#### Plotting Momentum Constraint in Azimuthal Direction (ph)

```python
# Load data from output file
data = np.loadtxt('output.txt')

# Extract momentum constraint violation in azimuthal direction
momentum_violation_ph = data[:, 6]

# Create plot for momentum constraint in azimuthal direction
plt.plot(momentum_violation_ph)
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Momentum Constraint Violation (ph)')
savefig('Momentum_Constraint_ph.png')
```

This code plots the momentum constraint violation in the azimuthal direction.

### Mathematical Background

The numerical relativity code uses various mathematical constants and data structures to perform calculations. The plotting of momentum constraints is an important step in understanding the results of the simulation.

$$\label{mathematical_background}$$

Let $G$ be the gravitational constant, and let $\hbar$ be the reduced Planck constant**Defining Figure Size**
=======================

### Overview of the Notebook

This notebook covers the process of defining the size of a figure in matplotlib.

### Theory Review

#### Introduction to Matplotlib Figures

In matplotlib, figures are used to display plots and other visualizations. The `figsize` parameter is used to set the size of the figure.

```python
import matplotlib.pyplot as plt
```

This code imports the matplotlib library.

### Code Implementation

#### Defining Figure Size

```python
# Define the size of the overall figure
fig = plt.figure(figsize=(12, 12))
```

This code defines a figure with a width of 12 inches and a height of 12 inches.

### Mathematical Background

The size of the figure is typically measured in inches or centimeters. The `figsize` parameter takes a tuple of two values: the first value represents the width, and the second value represents the height.

$$\label{mathematical_background}$$

Let $W$ be the width of the figure and $H$ be the height of the figure. Then, we can define the following mathematical relationship:

$$
figsize = (W, H)
$$

where $W$ and $H$ are measured in inches or centimeters.

### Example Use Cases

*   Defining a figure size for a specific simulation or problem.
*   Modifying the figure size to accommodate different simulations or problems.

#### Setting Figure Size in Different Units

```python
# Set figure size in inches
fig = plt.figure(figsize=(12, 12))

# Set figure size in centimeters
fig = plt.figure(figsize=(30.48, 30.48))
```

This code sets the figure size in both inches and centimeters.

### Tips and Tricks

*   Use `plt.figure()` to create a new figure.
*   Use the `figsize` parameter to set the size of the figure.
*   Use inches or centimeters as units for figure size.**Configuring Plot Parameters**
==============================

### Overview of the Notebook

This notebook covers the process of configuring plot parameters for a set of plots.

### Theory Review

#### Introduction to Plot Configuration

Plot configuration is an essential step in creating multiple plots with different parameters. The `num_plots` variable determines the number of plots, and the `Labels` list stores the labels for each plot.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Defining Plot Parameters

```python
# Define the number of plots
num_plots = 4

# Define labels for each plot
Labels=[r"W",  r"\mathcal{H}",r"\mathcal{M}^r",r"\mathcal{M}^{\theta}",r"\mathcal{M}^{\phi}"]

# Check if the coordinate system is Cartesian
if "Cartesian" in CoordSystem:
    # Update labels for Cartesian coordinates
    Labels=[r"W",  r"\mathcal{H}",r"\mathcal{M}^x",r"\mathcal{M}^y",r"\mathcal{M}^z"]

# Define index for data
data_idx=[    0,                1,               2,               3,               4]

# Define plot list
plotlist = [1, 2, 3, 4]

# Check if momentum is enabled
if dictID[initial_data_string].EnableMomentum == False:
    # Update plot list to only include Hamiltonian
    plotlist = [1]
```

This code defines the number of plots, labels, data indices, and plot list.

### Mathematical Background

The `num_plots` variable determines the number of plots, which is used to configure the figure size and layout. The `Labels` list stores the labels for each plot, which are used as titles and axis labels.

$$\label{mathematical_background}$$

Let $N$ be the number of plots. Then, we can define the following mathematical relationship:

$$
num\_plots = N
$$

### Example Use Cases

*   Configuring plot parameters for a set of plots.
*   Modifying the `num_plots` variable to change the number of plots.

#### Updating Plot List Based on Momentum Enablement

```python
if dictID[initial_data_string].EnableMomentum ==**Initializing Axis/Plot Array**
================================

### Overview of the Notebook

This notebook covers the process of initializing an array to store axes/plots.

### Theory Review

#### Introduction to Plotting

In matplotlib, plots are created using a combination of figures and axes. The `plotlist` variable determines which plots to create.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Initializing Axis/Plot Array

```python
# Initialize axis/plot array
axN = []

# Loop through plot list
for i in plotlist:
    # Determine which plot to create (0-indexed)
    whichplot = i - 1
    
    # Create a new subplot for each plot
    axN.append(plt.subplot(num_plots, 1, whichplot + 1))
    
    # Set title and labels for each subplot
    axN[-1].set_title(r"Plot %d: %s" % (i, Labels[i - 1]))
    axN[-1].set_xlabel('x')
    axN[-1].set_ylabel(Labels[i-1])
```

This code initializes an array to store axes/plots and creates a new subplot for each plot in the `plotlist`.

### Mathematical Background

The `axN` array is used to store axes/plots, with each element corresponding to a separate plot. The `whichplot` variable determines which plot to create.

$$\label{mathematical_background}$$

Let $N$ be the number of plots and $i$ be the index of a specific plot. Then, we can define the following mathematical relationship:

$$
axN[i] = \text{axis/plot corresponding to } i^{th} \text{ plot}
$$

### Example Use Cases

*   Initializing an array to store axes/plots for multiple plots.
*   Modifying the `plotlist` variable to change which plots are created.

#### Creating Subplots with Different Titles and Labels

```python
for i in plotlist:
    # Determine which plot to create (0-indexed)
    whichplot = i - 1
    
    # Create a new subplot for each plot
    axN.append(plt.subplot(num_plots, 1, whichplot + 1))
    
    # Set title and labels for each subplot
    axN[-1].set_title(r"Plot %d: %s"**Generating Subplots for Constraints**
=====================================

### Overview of the Notebook

This notebook covers the process of generating subplots for constraints.

### Theory Review

#### Introduction to Subplots

Subplots are a way to display multiple plots in a single figure. In this case, we want to generate four subplots for each constraint (Hamiltonian, momentum radial, momentum polar, and momentum azimuthal).

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Generating Subplot for Each Constraint

```python
# Generate subplot for each constraint
for i in plotlist:
    # Determine which plot to create (0-indexed)
    whichplot = i - 1
    
    # Add new subplot to figure
    ax = fig.add_subplot(221+whichplot)
    
    # Append subplot to axis/plot array
    axN.append(ax)
```

This code generates a new subplot for each constraint using the `add_subplot()` method of the figure object.

### Mathematical Background

The subplots are arranged in a 2x2 grid, with each subplot corresponding to a different constraint. The `whichplot` variable determines which subplot to create.

$$\label{mathematical_background}$$

Let $N$ be the number of plots and $i$ be the index of a specific plot. Then, we can define the following mathematical relationship:

$$
ax = \text{subplot corresponding to } i^{th} \text{ constraint}
$$

### Example Use Cases

*   Generating subplots for multiple constraints.
*   Modifying the `plotlist` variable to change which constraints are plotted.

#### Adding Subplots to Figure

```python
# Generate subplot for each constraint
for i in plotlist:
    # Determine which plot to create (0-indexed)
    whichplot = i - 1
    
    # Add new subplot to figure
    ax = fig.add_subplot(221+whichplot)
    
    # Append subplot to axis/plot array
    axN.append(ax)
```

This code adds the subplots to the figure using the `add_subplot()` method.

### Theory Review

#### Subplot Layout

The subplots are arranged in a 2x2 grid, with each subplot corresponding to a different constraint. The `221+whichplot` argument determines which subplot to create.

$$
\begin{bmatrix}
1 & 2 \\
3 &**Configuring Subplot Layout and Plotting Data**
=============================================

### Overview of the Notebook

This notebook covers the process of configuring the subplot layout and plotting data.

### Theory Review

#### Introduction to Subplot Configuration

Subplots are a way to display multiple plots in a single figure. In this case, we want to configure each subplot to display the numerical error for each constraint.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Configuring Subplot Layout and Plotting Data

```python
# Configure subplot layout and plot data
for i in range(len(plotlist)):
    # Determine which plot to create (0-indexed)
    whichplot = i
    
    # Set x-axis label
    axN[whichplot].set_xlabel(r'$x/M$')
    
    # Set y-axis label
    axN[whichplot].set_ylabel(r'$y/M$')
    
    # Set title for subplot
    axN[whichplot].set_title(r"$96^3$ Numerical Err.: $log_{10}|"+Labels[i]+r"|$")
    
    # Plot data using imshow function
    figure = plt.imshow(output_grid_data[i], extent=(-xy_extent,xy_extent, -xy_extent,xy_extent))
    
    # Create colorbar for plot
    cb = plt.colorbar(figure)
```

This code configures each subplot to display the numerical error for each constraint.

### Mathematical Background

The subplots are arranged in a 2x2 grid, with each subplot corresponding to a different constraint. The `whichplot` variable determines which subplot to create.

$$\label{mathematical_background}$$

Let $N$ be the number of plots and $i$ be the index of a specific plot. Then, we can define the following mathematical relationship:

$$
ax_N[i] = \text{subplot corresponding to } i^{th} \text{ constraint}
$$

### Example Use Cases

*   Configuring subplot layout for multiple constraints.
*   Modifying the `plotlist` variable to change which constraints are plotted.

#### Plotting Data using Imshow Function

```python
# Configure subplot layout and plot data
for i in range(len(plotlist)):
    # Determine which plot to create (0-indexed)
    whichplot = i
    
    # Set x-axis label
    axN[whichplot].set_xlabel(r**Adjusting Plot Spacing**
=========================

### Overview of the Notebook

This notebook covers the process of adjusting the spacing between plots using the `tight_layout` function.

### Theory Review

#### Introduction to Plot Layout

Plot layout is an essential aspect of creating visually appealing and informative plots. The `tight_layout` function is used to adjust the spacing between plots.

```python
import matplotlib.pyplot as plt
```

This code imports the matplotlib library.

### Code Implementation

#### Adjusting Plot Spacing using Tight Layout Function

```python
# Adjust plot spacing using tight layout function
plt.tight_layout(pad=4)
```

This code adjusts the spacing between plots by setting a padding of 4 points.

### Mathematical Background

The `tight_layout` function is used to adjust the spacing between plots. The `pad` parameter determines the amount of padding between plots.

$$\label{mathematical_background}$$

Let $P$ be the number of plots and $p$ be the padding between plots. Then, we can define the following mathematical relationship:

$$
\text{padding} = p \cdot P
$$

### Example Use Cases

*   Adjusting plot spacing for multiple plots.
*   Modifying the `pad` parameter to change the amount of padding.

#### Plotting Convergence of Numerical Errors

```python
# Create a figure with four subplots
fig, ax = plt.subplots(2, 2)

# Set labels and titles for each subplot
ax[0, 0].set_title('Plot 1')
ax[0, 1].set_title('Plot 2')
ax[1, 0].set_title('Plot 3')
ax[1, 1].set_title('Plot 4')

# Adjust plot spacing using tight layout function
plt.tight_layout(pad=4)
```

This code creates a figure with four subplots and adjusts the spacing between plots.

### Theory Review

#### Convergence of Numerical Errors

The convergence of numerical errors is an essential aspect of numerical relativity. The `96^3` grid is used as a reference, and the `48^3` grid is used to demonstrate the expected order of convergence.

$$
\text{convergence rate} = \left(\frac{\Delta x^{i}_{96}}{\Delta x^{i}_{48}}\right)^4
$$**Plot Settings**
================

### Overview of the Notebook

This notebook covers the process of setting plot settings for a numerical relativity simulation.

### Theory Review

#### Introduction to Plot Settings

In a numerical relativity simulation, plot settings are used to configure various aspects of the plots generated by the code. The `x_extent` variable is used to set the extent of the x-axis in the plots.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Setting Plot Settings

```python
# Set plot settings
x_extent = domain_size
y_extent = domain_size
```

This code sets the extent of the x and y axes using the `domain_size` variable.

### Mathematical Background

The `x_extent` and `y_extent` variables are used to set the limits of the x and y axes in the plots. The `domain_size` variable is a parameter that controls the size of the domain.

$$\label{mathematical_background}$$

Let $D$ be the domain size, then we can define the following mathematical relationship:

$$
x\_extent = D
$$

### Example Use Cases

*   Setting plot settings for multiple plots.
*   Modifying the `domain_size` variable to change the size of the domain.

#### Plotting Data using Imshow Function

```python
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Generate sample data
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x, y)

# Set plot settings
x_extent = domain_size
y_extent = domain_size

# Plot data using imshow function
plt.imshow(X, extent=(x_extent, x_extent, y_extent, y_extent))
```

This code generates sample data and plots it using the `imshow` function.

### Theory Review

#### Extent of Axes

The `extent` parameter in the `imshow` function is used to set the limits of the axes. The `x_extent` and `y_extent` variables are used to set the extent of the x and y axes.

$$
\text{extent} = (x\_extent, y\_extent)
$$**Plot Settings**
================

### Overview of the Notebook

This notebook covers the process of setting plot settings for a numerical relativity simulation.

### Theory Review

#### Introduction to Plot Settings

In a numerical relativity simulation, plot settings are used to configure various aspects of the plots generated by the code. The `sample_numpts_x` variable is used to set the number of sample points in the x-direction for plotting data.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Setting Plot Settings

```python
# Set plot settings
x_extent = domain_size
y_extent = domain_size
sample_numpts_x = 100
```

This code sets the extent of the x and y axes using the `domain_size` variable, and sets the number of sample points in the x-direction.

### Mathematical Background

The `sample_numpts_x` variable is used to set the number of sample points in the x-direction. The number of sample points determines the resolution of the plot.

$$\label{mathematical_background}$$

Let $N$ be the number of sample points, then we can define the following mathematical relationship:

$$
sample\_numpts_x = N
$$

### Example Use Cases

*   Setting plot settings for multiple plots.
*   Modifying the `domain_size` variable to change the size of the domain.

#### Plotting Data using Imshow Function

```python
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Generate sample data
x = np.linspace(-x_extent, x_extent, sample_numpts_x)
y = np.linspace(-y_extent, y_extent, sample_numpts_x)
X, Y = np.meshgrid(x, y)

# Plot data using imshow function
plt.imshow(X, extent=(-x_extent, +x_extent, -y_extent, +y_extent))
```

This code generates sample data and plots it using the `imshow` function.

### Theory Review

#### Extent of Axes

The `extent` parameter in the `imshow` function is used to set the limits of the axes. The `-x_extent` to `+x_extent` extent sets the plot range from `-x_extent` to `+x_extent`.

$$
\text{extent} = (-x\_extent, +x\_extent)
$$

### Tips and Tricks

*   Use `np.linspace` function to generate**Interpolation Method**
=====================

### Overview of the Notebook

This notebook covers the process of setting the interpolation method for plotting data.

### Theory Review

#### Introduction to Interpolation Methods

In numerical relativity, interpolation methods are used to determine the number of points to plot. The `interp_method` variable is used to set the interpolation method.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Setting Interpolation Method

```python
# Set interpolation method
interp_method = "linear"
```

This code sets the interpolation method to linear.

### Mathematical Background

The `interp_method` variable is used to set the interpolation method. The available interpolation methods include:

*   Linear: This method uses linear interpolation between points.
*   Spline: This method uses spline interpolation between points.
*   Nearest: This method uses nearest neighbor interpolation between points.

$$\label{mathematical_background}$$

Let $M$ be the interpolation method, then we can define the following mathematical relationship:

$$
interp\_method = M
$$

### Example Use Cases

*   Setting interpolation method for multiple plots.
*   Modifying the `interp_method` variable to change the interpolation method.

#### Plotting Data using Imshow Function with Linear Interpolation

```python
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Generate sample data
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Set interpolation method to linear
interp_method = "linear"

# Plot data using imshow function with linear interpolation
plt.imshow(X, extent=(-1, 1, -1, 1), interpolation=interp_method)
```

This code generates sample data and plots it using the `imshow` function with linear interpolation.

### Theory Review

#### Interpolation Methods

The available interpolation methods include:

*   Linear: This method uses linear interpolation between points.
*   Spline: This method uses spline interpolation between points.
*   Nearest: This method uses nearest neighbor interpolation between points.

$$
\text{interpolation} = \left\{
\begin{array}{ll}
linear & M = "linear" \\
spline & M = "spline" \\
nearest & M = "nearest"
\end{array}
\right**Interpolation Method Options**
================================

### Overview of the Notebook

This notebook covers the process of setting the interpolation method for plotting data.

### Theory Review

#### Introduction to Interpolation Methods

In numerical relativity, interpolation methods are used to determine the number of points to plot. The `interp_method` variable is used to set the interpolation method.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Setting Interpolation Method Options

```python
# Set interpolation method options
interp_method = "linear"  # Recommended
# interp_method = "nearest"  # Not recommended; gridpoints are off-axis
# interp_method = "cubic"    # Optional
```

This code sets the interpolation method to linear, which is the recommended option.

### Mathematical Background

The `interp_method` variable is used to set the interpolation method. The available interpolation methods include:

*   Linear: This method uses linear interpolation between points.
*   Nearest: This method uses nearest neighbor interpolation between points (not recommended).
*   Cubic: This method uses cubic interpolation between points (optional).

$$\label{mathematical_background}$$

Let $M$ be the interpolation method, then we can define the following mathematical relationship:

$$
interp\_method = M
$$

### Example Use Cases

*   Setting interpolation method for multiple plots.
*   Modifying the `interp_method` variable to change the interpolation method.

#### Plotting Data using Imshow Function with Linear Interpolation

```python
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Generate sample data
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Set interpolation method to linear
interp_method = "linear"

# Plot data using imshow function with linear interpolation
plt.imshow(X, extent=(-1, 1, -1, 1), interpolation=interp_method)
```

This code generates sample data and plots it using the `imshow` function with linear interpolation.

### Theory Review

#### Interpolation Methods

The available interpolation methods include:

*   Linear: This method uses linear interpolation between points.
$$
\text{linear} = \mathcal{L}
$$

*   Nearest: This method uses nearest neighbor interpolation between points (not recommended**Extracting 1D Data from 2D Grids**
=====================================

### Overview of the Notebook

This notebook covers the process of extracting 1D data from 2D grids for plotting.

### Theory Review

#### Introduction to Extracting 1D Data

In numerical relativity, it's often necessary to extract 1D data from 2D grids for plotting. The `extract_1D_slice_from_2D_data` function is used to perform this task.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Extracting 1D Data from 2D Grids

```python
# Initialize lists to store 1D data
output_1D_grid_data48 = []
output_1D_grid_data96 = []

# Loop through different quantities (Hamiltonian, momenta)
for i in [4, 5, 6, 7]:
    # Extract 1D slice from 2D data for 48^3 grid
    output_grid_x96, output_1D_grid_data48_i = \
        plot2D.extract_1D_slice_from_2D_data('out48.txt', 0.0,
                                      0,1,i, [-x_extent, x_extent], sample_numpts_x=sample_numpts_x,
                                      interp_method=interp_method)
    
    # Extract 1D slice from 2D data for 96^3 grid
    output_grid_x48, output_1D_grid_data96_i = \
        plot2D.extract_1D_slice_from_2D_data('out96.txt', 0.0,
                                      0,1,i, [-x_extent, x_extent], sample_numpts_x=sample_numpts_x,
                                      interp_method=interp_method)
    
    # Append extracted data to lists
    output_1D_grid_data48 += [output_1D_grid_data48_i]
    output_1D_grid_data96 += [output_1D_grid_data96_i]

# Define plot titles based on coordinate system
PlotTitles=[r"\mathcal{H}",r"\mathcal{M}^r",r"\mathcal{M}^{\theta}",r"\mathcal{M}^{\phi}"]
if "Cartesian" in CoordSystem:
    PlotTitles=[r"\mathcal{H}",r"\mathcal{M}**Creating Plots**
================

### Overview of the Notebook

This notebook covers the process of creating plots for the Hamiltonian and momenta.

### Theory Review

#### Introduction to Plotting

In numerical relativity, it's often necessary to visualize data to understand its behavior. We will create four plots: one for the Hamiltonian and three for the momenta.

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Defining Plot Titles and Data

```python
# Define plot titles based on coordinate system
PlotTitles=[r"\mathcal{H}",r"\mathcal{M}^r",r"\mathcal{M}^{\theta}",r"\mathcal{M}^{\phi}"]
if "Cartesian" in CoordSystem:
    PlotTitles=[r"\mathcal{H}",r"\mathcal{M}^{x}",r"\mathcal{M}^{y}",r"\mathcal{M}^{z}"]

# Define data for plots
axN = []
```

This code defines the plot titles and initializes an array to store the axes.

### Mathematical Background

The Hamiltonian ($\mathcal{H}$) is a measure of the energy density in the system. The momenta ($\mathcal{M}$) are measures of the momentum density in the system.

$$\label{mathematical_background}$$

Let $\mathcal{H}$ be the Hamiltonian, and let $\mathcal{M}_i$ be the $i^{th}$ component of the momentum. Then we can define the following mathematical relationships:

$$
\text{Hamiltonian: } \mathcal{H}
$$

$$
\text{Momenta: } \mathcal{M}_i = \left\{
\begin{array}{ll}
r & i=1 \\
\theta & i=2 \\
\phi & i=3
\end{array}
\right.
$$

### Example Use Cases

*   Creating plots for the Hamiltonian and momenta in different coordinate systems.

#### Creating Plots with matplotlib

```python
# Create a figure with four subplots
fig, ax = plt.subplots(2, 2)

# Set plot titles
for i in range(len(axN)):
    ax[i // 2, i % 2].set_title(PlotTitles[i])

# Plot data using**Constraining Variables**
=======================

### Overview of the Notebook

This notebook covers the process of constraining variables in a numerical relativity simulation.

### Theory Review

#### Introduction to Constrained Variables

In numerical relativity, it's often necessary to constrain variables such as $r$, $\theta$, and $\phi$ to specific values or ranges. This is typically done using constraints such as:

*   Radial distance constraint: $r \in [0, r_{max}]$
*   Polar angle constraint: $\theta \in [0, \pi]$
*   Azimuthal angle constraint: $\phi \in [0, 2\pi)$

```python
import numpy as np
```

This code imports the NumPy library.

### Code Implementation

#### Defining Constraints

```python
# Define constraints for radial distance, polar angle, and azimuthal angle
r_constraint = [0, 10]
theta_constraint = [0, np.pi]
phi_constraint = [0, 2 * np.pi]

# Print constraints
print("Radial distance constraint:", r_constraint)
print("Polar angle constraint:", theta_constraint)
print("Azimuthal angle constraint:", phi_constraint)
```

This code defines the constraints for radial distance, polar angle, and azimuthal angle using lists.

### Mathematical Background

The constraints are used to define the range of values that the variables can take. For example, the radial distance constraint $r \in [0, r_{max}]$ means that the value of $r$ must be between 0 and $r_{max}$.

$$\label{mathematical_background}$$

Let $r$, $\theta$, and $\phi$ be the variables to be constrained. Then we can define the following mathematical relationships:

*   Radial distance constraint: $$ r \in [0, r_{max}]$$
*   Polar angle constraint: $$\theta \in [0, \pi]$$
*   Azimuthal angle constraint: $$\phi \in [0, 2\pi)$$

### Example Use Cases

*   Constraining variables in a numerical relativity simulation.
*   Modifying the constraints to change the range of values for the variables.

#### Using Constraints in Numerical Relativity Simulation

```python
# Define initial conditions
r_initial = 1
theta_initial = np.pi / 4
phi_initial = np.pi / **Defining Figure Size**
=====================

### Overview of the Notebook

This notebook covers the process of defining the size of the overall figure in a matplotlib plot.

### Theory Review

#### Introduction to Figure Size

In matplotlib, the figure size is an essential aspect of creating visually appealing and informative plots. The `figsize` parameter is used to set the size of the figure.

```python
import matplotlib.pyplot as plt
```

This code imports the `matplotlib.pyplot` module, which provides functions for creating static, animated, and interactive visualizations in python.

### Code Implementation

#### Defining Figure Size

```python
# Define the size of the overall figure
fig = plt.figure(figsize=(12, 12))
```

This code defines a new figure object using the `plt.figure()` function and sets its size to (12, 12) using the `figsize` parameter.

### Mathematical Background

The figure size is defined as a tuple (width, height) in inches. The width of the figure can be any positive real number, while the height must be non-negative.

$$\label{mathematical_background}$$

Let $w$ and $h$ be the width and height of the figure respectively. Then we can define the following mathematical relationship:

$$
\text{Figure size: } (w,h)
$$

### Example Use Cases

*   Defining the size of the overall figure in a matplotlib plot.
*   Modifying the `figsize` parameter to change the size of the figure.

#### Using Different Figure Sizes

```python
# Define different figure sizes
fig_small = plt.figure(figsize=(6, 6))
fig_medium = plt.figure(figsize=(8, 8))
fig_large = plt.figure(figsize=(16, 16))

print("Small figure size:", fig_small.get_size_inches())
print("Medium figure size:", fig_medium.get_size_inches())
print("Large figure size:", fig_large.get_size_inches())
```

This code defines three different figures with sizes (6, 6), (8, 8), and (16, 16) respectively and prints their sizes in inches.

### Tips and Tricks

*   Use the `figsize` parameter to set the size of the figure.
*   Define multiple figures with different sizes using the `plt.figure()` function.**Creating Multiple Subplots**
=============================

### Overview of the Notebook

This notebook covers the process of creating multiple subplots within a figure.

### Theory Review

#### Introduction to Subplots

In matplotlib, subplots are used to create multiple plots within a single figure. The `subplots()` function is used to create a figure with one or more subplots.

```python
import matplotlib.pyplot as plt
```

This code imports the `matplotlib.pyplot` module, which provides functions for creating static, animated, and interactive visualizations in python.

### Code Implementation

#### Creating Multiple Subplots

```python
# Define number of plots to create
num_plots = 4

# Create a figure with multiple subplots
fig, ax = plt.subplots(num_rows=2, num_cols=num_plots // 2, figsize=(8, 8))

# Iterate over the plots and set their titles
for p in range(num_plots):
    # Set title for each subplot
    ax[p].set_title(f"Plot {p+1}")
```

This code defines a figure with multiple subplots using the `subplots()` function. The number of rows is set to 2, and the number of columns is calculated by dividing the total number of plots by 2. The size of the figure is also specified.

### Mathematical Background

The number of subplots can be calculated as follows:

$$\label{mathematical_background}$$

Let $N$ be the total number of plots. Then we can define the following mathematical relationship:

$$
\text{Number of rows: } r = \left\lfloor \frac{N}{c} \right\rfloor
$$

$$
\text{Number of columns: } c = \left\lceil \frac{N}{r} \right\rceil
$$

where $\lfloor x \rfloor$ is the floor function and $\lceil x \rceil$ is the ceiling function.

### Example Use Cases

*   Creating multiple subplots within a single figure.
*   Modifying the `num_rows` and `num_cols` parameters to change the layout of the subplots.

#### Creating Different Subplot Layouts

```python
# Define different subplot layouts
fig, ax = plt.subplots(num_rows=3, num_cols=2, figsize=(8, 12))
fig, ax = plt.subplots(num_rows=num_plots // 2, num_cols=**Plotting Constraints**
=====================

### Overview of the Notebook

This notebook covers the process of plotting the constraints for radial distance, polar angle, and azimuthal angle.

### Theory Review

#### Introduction to Plotting Constraints

In numerical relativity, it's often necessary to plot the constraints for radial distance, polar angle, and azimuthal angle. This can be done using matplotlib.

```python
import numpy as np
import matplotlib.pyplot as plt
```

This code imports the NumPy library and the `matplotlib.pyplot` module.

### Code Implementation

#### Defining Constraints

```python
# Define constraints for radial distance, polar angle, and azimuthal angle
r_constraint = [0, 10]
theta_constraint = [0, np.pi]
phi_constraint = [0, 2 * np.pi]

# Print constraints
print("Radial distance constraint:", r_constraint)
print("Polar angle constraint:", theta_constraint)
print("Azimuthal angle constraint:", phi_constraint)
```

This code defines the constraints for radial distance, polar angle, and azimuthal angle using lists.

#### Plotting Constraints

```python
# Create a figure with multiple subplots
fig, ax = plt.subplots(3, 1, figsize=(8, 12))

# Loop to cycle through our constraints and plot the data
for i, (constraint_name, constraint) in enumerate(zip(['Radial Distance', 'Polar Angle', 'Azimuthal Angle'],
                                                     [r_constraint, theta_constraint, phi_constraint])):
    # Plot constraint
    ax[i].plot(constraint)
    
    # Set title for each subplot
    ax[i].set_title(constraint_name)

# Layout so plots do not overlap
plt.tight_layout()
```

This code creates a figure with multiple subplots and loops through the constraints to plot them. The `tight_layout()` function is used to ensure that the plots do not overlap.

### Mathematical Background

The constraints for radial distance, polar angle, and azimuthal angle can be represented as:

$$\label{mathematical_background}$$

Let $r$, $\theta$, and $\phi$ be the variables to be constrained. Then we can define the following mathematical relationships:

*   Radial distance constraint: $$ r \in [0, r_{max}]$$
*   Polar angle constraint: $$\theta \in [0, \pi]$$
*   Azimuthal angle constraint: $$\phi \**Generating Subplots for Constraints**
=====================================

### Overview of the Notebook

This notebook covers the process of generating subplots for each constraint in a numerical relativity simulation.

### Theory Review

#### Introduction to Generating Subplots

In matplotlib, the `add_subplot()` function is used to add a subplot to an existing figure. In this case, we are generating subplots for each constraint.

```python
import numpy as np
import matplotlib.pyplot as plt
```

This code imports the NumPy library and the `matplotlib.pyplot` module.

### Code Implementation

#### Defining Constraints

```python
# Define constraints for radial distance, polar angle, and azimuthal angle
r_constraint = [0, 10]
theta_constraint = [0, np.pi]
phi_constraint = [0, 2 * np.pi]

# Print constraints
print("Radial distance constraint:", r_constraint)
print("Polar angle constraint:", theta_constraint)
print("Azimuthal angle constraint:", phi_constraint)
```

This code defines the constraints for radial distance, polar angle, and azimuthal angle using lists.

#### Generating Subplots

```python
# Create a figure with multiple subplots
fig = plt.figure(figsize=(8, 12))

# Initialize list to store axes
axN = []

# Loop through each constraint and generate subplot
for p in range(3):
    # Generate subplot for current constraint
    ax = fig.add_subplot(221+p)
    
    # Append axis to list of axes
    axN.append(ax)

# Layout so plots do not overlap
plt.tight_layout()
```

This code creates a figure with multiple subplots and generates subplots for each constraint using the `add_subplot()` function. The `tight_layout()` function is used to ensure that the plots do not overlap.

### Mathematical Background

The constraints for radial distance, polar angle, and azimuthal angle can be represented as:

$$\label{mathematical_background}$$

Let $r$, $\theta$, and $\phi$ be the variables to be constrained. Then we can define the following mathematical relationships:

*   Radial distance constraint: $$ r \in [0, r_{max}]$$
*   Polar angle constraint: $$\theta \in [0, \pi]$$
*   Azimuthal angle constraint: $$\phi \in [0, 2\pi]$$

### Example Use Cases

*   Generating subplots for multiple**Customizing Subplot with Legend**
====================================

### Overview of the Notebook

This notebook covers the process of customizing a subplot with a legend.

### Theory Review

#### Introduction to Customizing Subplots

In matplotlib, subplots can be customized using various functions and methods. This includes setting titles, labels, and legends.

```python
import numpy as np
import matplotlib.pyplot as plt
```

This code imports the NumPy library and the `matplotlib.pyplot` module.

### Code Implementation

#### Customizing Subplot with Legend

```python
# Set title for current subplot
axN[p].set_title('Plot Demonstrating $4^{th}$-Order Convergence of $'+PlotTitles[p]+'$')

# Set x-axis label
axN[p].set_xlabel(r"$x/M$")

# Set y-axis label
axN[p].set_ylabel("$log_{10}$(Relative Error)")

# Plot data on current subplot
ax.plot(output_grid_x96, output_1D_grid_data96[p], 'k-', label='Nr=96')
ax.plot(output_grid_x48, output_1D_grid_data48[p] + 4*np.log10(48./96.), 'k--', label='Nr=48, mult by (48/96)^4')

# Set y-axis limits
ax.set_ylim([-15.2,4.5])

# Create legend for current subplot
legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large')

# Customize legend frame
legend.get_frame().set_facecolor('C1')
```

This code sets the title, labels, and plots data on a subplot with a legend.

### Mathematical Background

The legend is used to label each line in the plot. The `loc` parameter is set to `'lower right'` to place the legend at the lower right corner of the plot.

$$\label{mathematical_background}$$

Let $x$ be the x-axis, and let $y$ be the y-axis. Then we can define the following mathematical relationships:

*   Label for first line: $$ label_{1} = \text{Nr}=96 $$
*   Label for second line: $$ label_{2} = \text{Nr}=48, mult by (48/96)^4$$

### Example Use Cases

*   Customizing subplots with legends in matplotlib.
*   Modifying the `loc` parameter to**Adjusting Spacing between Plots**
==================================

### Overview of the Notebook

This notebook covers the process of adjusting the spacing between plots in a matplotlib figure.

### Theory Review

#### Introduction to Adjusting Plot Spacing

In matplotlib, the `tight_layout()` function is used to automatically adjust the layout so that the subplots fit well in the figure. This includes adjusting the spacing between plots.

```python
import matplotlib.pyplot as plt
```

This code imports the `matplotlib.pyplot` module.

### Code Implementation

#### Adjusting Plot Spacing

```python
# Set plot spacing using tight_layout function
plt.tight_layout(pad=4)
```

This code adjusts the spacing between plots in a figure using the `tight_layout()` function. The `pad` parameter is set to 4, which means that there will be a 4-point padding between the subplots.

### Mathematical Background

The spacing between plots can be represented as:

$$\label{mathematical_background}$$

Let $s$ be the spacing between plots, and let $p$ be the padding between plots. Then we can define the following mathematical relationship:

$$
s = p \times \text{font size}
$$

### Example Use Cases

*   Adjusting plot spacing in matplotlib.
*   Modifying the `pad` parameter to change the spacing between plots.

#### Plotting with Tight Layout

```python
# Create a figure with multiple subplots
fig, ax = plt.subplots(2, 2)

# Add data to subplots
ax[0, 0].plot([1, 2, 3])
ax[0, 1].plot([4, 5, 6])
ax[1, 0].plot([7, 8, 9])
ax[1, 1].plot([10, 11, 12])

# Adjust plot spacing using tight_layout function
plt.tight_layout(pad=4)

# Display the plot
plt.show()
```

This code creates a figure with multiple subplots and adjusts the spacing between plots using the `tight_layout()` function.

### Tips and Tricks

*   Use the `pad` parameter to adjust the spacing between plots.
*   Modify the `font.size` property to change the font size of the plot.**Outputting Notebook as LaTeX-formatted PDF**
=====================================================

### Overview of the Notebook

This notebook covers the process of outputting a Jupyter notebook as a LaTeX-formatted PDF file.

### Theory Review

#### Introduction to Outputting Notebooks

In Jupyter, notebooks can be outputted in various formats, including LaTeX-formatted PDF. This allows for easy conversion of notebooks into professional-looking documents.

```python
import nbformat
from IPython.display import display, Javascript
```

This code imports the necessary libraries for outputting notebooks as LaTeX-formatted PDF.

### Code Implementation

#### Outputting Notebook as LaTeX-formatted PDF

```python
# Read notebook contents
nb = nbformat.read('notebook.ipynb', as_version=4)

# Convert to LaTeX format
latex_nb = nbformat.writes(nb, version=4)

# Output as LaTeX-formatted PDF file
with open('output.pdf', 'w') as f:
    f.write(latex_nb)
```

This code reads the contents of the notebook, converts it into LaTeX format using the `nbformat` library, and outputs it as a LaTeX-formatted PDF file.

### Mathematical Background

The output of the notebook is represented as a string in LaTeX format. The conversion process involves using the `nbformat.writes()` function to convert the notebook into LaTeX format.

$$\label{mathematical_background}$$

Let $N$ be the notebook contents, and let $L$ be the LaTeX-formatted string. Then we can define the following mathematical relationship:

$$
L = \text{nbformat.writes}(N)
$$

### Example Use Cases

*   Outputting a Jupyter notebook as a LaTeX-formatted PDF file.
*   Modifying the conversion process to change the output format.

#### Converting Notebook to LaTeX-formatted PDF

```python
# Convert notebook to LaTeX-formatted PDF
nbformat.write(nb, 'output.pdf')

# Display the resulting PDF file
display(Javascript('
var pdf = new jsPDF();
pdf.text("Hello World!", 10, 20);
pdf.save("hello.pdf");
'))
```

This code converts the notebook into a LaTeX-formatted PDF file and displays the resulting PDF using the `jsPDF` library.

### Tips and Tricks

*   Use the `nbformat` library to convert notebooks into LaTeX format.
*   Modify the conversion process to change the output format.**Converting Jupyter Notebook to LaTeX-formatted PDF**
=====================================================

### Overview of the Notebook

This notebook covers the process of converting a Jupyter notebook into a LaTeX-formatted PDF file.

### Theory Review

#### Introduction to Converting Notebooks

In Jupyter, notebooks can be converted into various formats, including LaTeX-formatted PDF. This allows for easy conversion of notebooks into professional-looking documents.

```python
import cmdline_helper as cmd
```

This code imports the `cmdline_helper` library, which is used to convert the notebook into a LaTeX-formatted PDF file.

### Code Implementation

#### Converting Notebook to LaTeX-formatted PDF

```python
# Convert notebook to LaTeX-formatted PDF
cmd.run_command(cmd.Command("nbconvert --to latex --template=slides Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data.ipynb"))
```

This code uses the `run_command()` function from the `cmdline_helper` library to convert the notebook into a LaTeX-formatted PDF file.

### Mathematical Background

The conversion process involves using the `nbconvert` command with the `--to latex` and `--template=slides` options. The resulting PDF file is generated in the root NRPy+ tutorial directory.

$$\label{mathematical_background}$$

Let $N$ be the notebook contents, and let $P$ be the LaTeX-formatted PDF file. Then we can define the following mathematical relationship:

$$
P = \text{nbconvert}(N)
$$

### Example Use Cases

*   Converting a Jupyter notebook into a LaTeX-formatted PDF file.
*   Modifying the conversion process to change the output format.

#### Notes on Conversion Process

The `--to latex` option specifies that we want to convert the notebook into a LaTeX-formatted PDF file. The `--template=slides` option specifies the template to use for the LaTeX-formatted PDF file.

### Tips and Tricks

*   Use the `cmdline_helper` library to convert notebooks into LaTeX-format.
*   Modify the conversion process to change the output format.

#### Checking the Generated PDF File

The generated PDF file can be found in the root NRPy+ tutorial directory. To check if the conversion was successful, you can open the PDF file using a PDF viewer.

```python
# Check if PDF file exists
import os
pdf_file = "Tutorial-Start_to_Finish-BSSNCurvilinear**Converting Jupyter Notebook to LaTeX-formatted PDF using NRPy+**
==============================================================

### Overview of the Notebook

This notebook covers the process of converting a Jupyter notebook into a LaTeX-formatted PDF file using the NRPy+ library.

### Theory Review

#### Introduction to NRPy+

NRPy+ is a multi-platform Python command-line interface for numerical relativity. It provides tools for solving Einstein's field equations and simulating gravitational waves.

```python
import cmdline_helper as cmd
```

This code imports the `cmdline_helper` library from NRPy+, which is used to convert the Jupyter notebook into a LaTeX-formatted PDF file.

### Code Implementation

#### Converting Notebook to LaTeX-formatted PDF using NRPy+

```python
# Convert notebook to LaTeX-formatted PDF using NRPy+
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Start_to_Finish-BSSNCurvilinear-Exact_Initial_Data")
```

This code uses the `output_Jupyter_notebook_to_LaTeXed_PDF()` function from NRPy+ to convert the notebook into a LaTeX-formatted PDF file.

### Mathematical Background

The conversion process involves using the `output_Jupyter_notebook_to_LaTeXed_PDF()` function, which takes the name of the Jupyter notebook as input. The resulting LaTeX-formatted PDF file is generated in the current working directory.

$$\label{mathematical_background}$$

Let $N$ be the notebook contents, and let $P$ be the LaTeX-formatted PDF file. Then we can define the following mathematical relationship:

$$
P = \text{output\_Jupyter\_notebook\_to\_LaTeXed\_PDF}(N)
$$

### Example Use Cases

*   Converting a Jupyter notebook into a LaTeX-formatted PDF file using NRPy+.
*   Modifying the conversion process to change the output format.

#### Notes on Conversion Process

The `output_Jupyter_notebook_to_LaTeXed_PDF()` function generates a LaTeX-formatted PDF file by compiling the LaTeX code in the notebook. The resulting PDF file is then saved as a file with the same name as the input notebook, but with a `.pdf` extension.

### Tips and Tricks

*   Use NRPy+ to convert Jupyter notebooks into LaTeX-formatted PDF files.
*   Modify the conversion process to change the output format.

#### Checking the Generated PDF File

The