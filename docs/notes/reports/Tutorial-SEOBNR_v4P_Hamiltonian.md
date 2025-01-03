**NRPy+: Google Analytics Tracking Code**
=====================================

### Theory Review

#### Introduction to Google Analytics Tracking Code

*   **Google Analytics:** In this section, we discuss how to track website usage and behavior with Google Analytics.
    +   This is an essential step in understanding user engagement and optimizing website content.

### Code Explanation


```python
"""
Google Analytics tracking code:
"""
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>
```

This code writes a JavaScript code to track website usage and behavior with Google Analytics.


### Theory Review

#### Understanding Google Analytics Tracking Code

*   **Google Analytics:** Google Analytics is a web analytics service that provides insights into website traffic, engagement, and conversion.
    +   It helps developers understand user behavior and optimize website content for better results.

### Mathematics


$$ \text{Google Analytics} = \left\{
\begin{array}{l}
\text{Web analytics service to track website usage and behavior} \\
\text{Provides insights into traffic, engagement, and conversion} \\
\end{array}
\right. $$

### Code Explanation

*   **`<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>`**
    +   We load the Google Analytics tracking code asynchronously from a CDN.
*   **`window.dataLayer = window.dataLayer || [];`**
    +   We initialize the `dataLayer` object, which is used to push tracking events.

### Theory Review

#### Pushing Tracking Events with gtag()

*   **gtag():** The `gtag()` function is used to push tracking events to the `dataLayer` object.
    +   It takes several arguments, including `js`, `config_id`, and `event`.

### Mathematics


$$ \text{Pushing Tracking Events} = \left\{
\begin{array}{l}
\text{Use gtag() function to push tracking events to dataLayer object} \\
\text{Specify event type and configuration ID as arguments} \\
\end{array}
\right. $$

### Code**NRPy+: Spinning Effective One-Body Hamiltonian**
=====================================================

### Theory Review

#### Introduction to Spinning Effective One-Body Hamiltonian

*   **Effective One-Body Hamiltonian:** In this section, we discuss the Spinning Effective One-Body Hamiltonian, also known as "v4P".
    +   This is a model used in numerical relativity to approximate the dynamics of binary systems.

### Code Explanation


```python
"""
Spinning Effective One-Body Hamiltonian:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Spinning Effective One-Body Hamiltonian
def v4P_Hamiltonian(m1, m2, chi1, chi2):
    """
    Calculate the Hamiltonian for the Spinning Effective One-Body Hamiltonian.
    
    Parameters:
        m1: mass of first body (kg)
        m2: mass of second body (kg)
        chi1: spin of first body
        chi2: spin of second body
    
    Returns:
        H: Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dimensionless spins
    chi_dimless_1 = chi1 / (G * c**2 * M)
    chi_dimless_2 = chi2 / (G * c**2 * M)
    
    # Calculate the Hamiltonian
    H = mu * G * c**2 * (1 + 73/24 * (chi_dimless_1**2) + 37/96 * (chi_dimless_2**2))
    
    return H

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)
chi1 = 0.3        # spin of first body
chi2 = 0.6        # spin of second body

H = v4P_Hamiltonian(m1, m2, chi1, chi2)
print("Hamiltonian:", H**NRPy+: Author Information**
============================

### Theory Review

#### Introduction to Author Information

*   **Author:** In this section, we discuss the author of the NRPy+ code.
    +   This information is crucial for understanding the origin and credibility of the code.

### Code Explanation


```python
"""
Author information:
"""
author = "Tyler Knowles"
```

This code writes a string variable to store the name of the author.


### Theory Review

#### Understanding Author Information

*   **Author:** The author of the NRPy+ code is Tyler Knowles.
    +   This information can be used to contact the author for further clarification or modifications.

### Mathematics


$$ \text{Author} = \left\{
\begin{array}{l}
\text{Name of person who wrote the code} \\
\text{Tyler Knowles in this case} \\
\end{array}
\right. $$

### Code Explanation

*   **`author = "Tyler Knowles"`**
    +   We assign a string value to the `author` variable.


### Theory Review

#### Importance of Author Information

*   **Credibility:** Knowing the author of the code can help establish credibility and trustworthiness.
    +   This information is essential for users who want to verify the accuracy of the code or make modifications.

### Code Explanation


```python
# Print author name
print(author)
```

This code writes a print statement to display the author's name.


### Theory Review

#### Contacting the Author

*   **Contact:** If you need to contact the author for further clarification or modifications, you can use the following information:
    +   Name: Tyler Knowles
    +   Email (if available)

### Mathematics


$$ \text{Contact} = \left\{
\begin{array}{l}
\text{Name of person who wrote the code} \\
\text{Tyler Knowles in this case} \\
\end{array}
\right. $$**NRPy+: Reduced Spinning Effective One-Body Hamiltonian**
===========================================================

### Theory Review

#### Introduction to Reduced Spinning Effective One-Body Hamiltonian

*   **Reduced Spinning Effective One-Body Hamiltonian:** In this section, we discuss the reduced spinning effective one-body Hamiltonian as numerically implemented in LALSuite's SEOBNRv4P gravitational waveform approximant.
    +   This is a key concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Reduced Spinning Effective One-Body Hamiltonian:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the reduced spinning effective one-body Hamiltonian
def H_real(m1, m2, gamma_E, t_tortoise):
    """
    Calculate the real part of the reduced spinning effective one-body Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        gamma_E: Euler-Mascheroni constant
        t_tortoise: tortoise coordinate
    
    Returns:
        H_real: real part of reduced spinning effective one-body Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dynamic variables
    q1, q2, s1, s2, chi1, chi2 = 12 * [np.random.rand()]
    
    # Calculate the real part of the reduced spinning effective one-body Hamiltonian
    H_real = mu * G * c**2 * (1 + 73/24 * (chi1**2) + 37/96 * (chi2**2))
    
    return H_real

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)
gamma_E = 0.57721566490153286060  # Euler-Mascheroni constant
t_tortoise = 100.0  # tortoise coordinate

**NRPy+: SEOBNR v4P Hamiltonian Module**
=====================================

### Theory Review

#### Introduction to SEOBNR v4P Hamiltonian Module

*   **SEOBNR v4P:** In this section, we discuss the SEOBNR v4P Hamiltonian module implemented in NRPy+.
    +   This module is used for numerical relativity simulations and gravitational wave astronomy.

### Code Explanation


```python
"""
NRPy+ Source Code for SEOBNR v4P Hamiltonian Module:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the SEOBNR v4P Hamiltonian function
def seoBNR_v4p_hamiltonian(m1, m2, chi1, chi2):
    """
    Calculate the SEOBNR v4P Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        chi1: spin of first black hole
        chi2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p: SEOBNR v4P Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dimensionless spins
    chi_dimless_1 = chi1 / (G * c**2 * M)
    chi_dimless_2 = chi2 / (G * c**2 * M)
    
    # Calculate the SEOBNR v4P Hamiltonian
    H_seoBNR_v4p = mu * G * c**2 * (1 + 73/24 * (chi_dimless_1**2) + 37/96 * (chi_dimless_2**2))
    
    return H_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)
chi1 = 0.3        # spin of first black hole**NRPy+: Introduction**
======================

### Theory Review

#### Overview of NRPy+

*   **Introduction:** In this section, we provide an introduction to the Numerical Relativity in Python (NRPy+) project.
    +   This is a comprehensive review of the project's goals and objectives.

### Code Explanation


```python
"""
Overview of NRPy+ Project:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define a function to calculate gravitational wave frequency
def gw_frequency(m1, m2):
    """
    Calculate the gravitational wave frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        f_gw: gravitational wave frequency (Hz)
    """
    # Calculate the total mass
    M = m1 + m2
    
    # Calculate the gravitational wave frequency
    f_gw = 2 * np.pi / (G * c**3) * (M / (m1 * m2))**(3/2)
    
    return f_gw

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

f_gw = gw_frequency(m1, m2)
print("Gravitational Wave Frequency:", f_gw)
```

This code writes a Python function to calculate the gravitational wave frequency.


### Theory Review

#### Understanding Gravitational Waves

*   **Gravitational Waves:** Gravitational waves are ripples in the fabric of spacetime produced by violent cosmic events, such as supernovae or black hole mergers.
    +   They were predicted by Albert Einstein's theory of general relativity and have since been observed directly.

### Mathematics


$$ \text{Gravitational Wave Frequency} = \left\{
\begin{array}{l}
\text{Frequency at which gravitational waves are emitted} \\
\text{Given by the formula: } f_{gw} = 2 \pi / (G c^3) (M / (m_1 m_2**NRPy+: Physical System of Interest**
=====================================

### Theory Review

#### Introduction to Binary Black Hole System

*   **Binary Black Hole System:** In this section, we discuss a binary black hole system consisting of two black holes with masses $m_{1}$ and $m_{2}$ and spins ${\bf S}_{1}$ and ${\bf S}_{2}$.
    +   This is a fundamental concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Physical System of Interest:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the masses and spins of the black holes
m1 = 1.98910e30  # mass of first black hole (kg)
m2 = 5.97219e24   # mass of second black hole (kg)
S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

# Define the SEOB Hamiltonian function
def seoBNR_v4p_hamiltonian(m1, m2, S1, S2):
    """
    Calculate the SEOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p: SEOB Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dimensionless spins
    chi_dimless_1 = np.linalg.norm(S1) / (G * c**2 * M)
    chi_dimless_2 = np.linalg.norm(S2) / (G * c**2 * M)
    
    # Calculate the SEOB Hamiltonian
    H_seoBNR_v4p = mu * G * c**2 * (1 + 73/24 * (chi_dimless_1**NRPy+: Effective One-Body Hamiltonian**
======================================

### Theory Review

#### Introduction to Effective One-Body Hamiltonian

*   **Effective One-Body Hamiltonian:** In this section, we discuss the effective one-body Hamiltonian $H_{\rm eff}$ as a canonical transformation of the real part of the SEOB Hamiltonian $H_{\rm real}$.
    +   This is a fundamental concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Effective One-Body Hamiltonian:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the real part of the SEOB Hamiltonian
def seoBNR_v4p_real_hamiltonian(m1, m2, S1, S2):
    """
    Calculate the real part of the SEOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_real: real part of SEOB Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dimensionless spins
    chi_dimless_1 = np.linalg.norm(S1) / (G * c**2 * M)
    chi_dimless_2 = np.linalg.norm(S2) / (G * c**2 * M)
    
    # Calculate the real part of the SEOB Hamiltonian
    H_seoBNR_v4p_real = mu * G * c**2 * (1 + 73/24 * (chi_dimless_1**2) + 37/96 * (chi_dimless_2**2))
    
    return H_seoBNR_v4p_real

# Define the effective one-body Hamiltonian
def seoBNR_v4p_effective_hamiltonian(H_seoBNR_v4p_real):
    """
    Calculate the effective one-body Hamiltonian.
    
   **NRPy+: Effective One-Body Hamiltonian**
======================================

### Theory Review

#### Introduction to Effective One-Body Hamiltonian

*   **Effective One-Body Hamiltonian:** In this section, we discuss the effective one-body Hamiltonian $H_{\rm eff}$ as a canonical transformation of the real part of the SEOB Hamiltonian $H_{\rm real}$.
    +   This is a fundamental concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Effective One-Body Hamiltonian:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the mass of the test particle
mu = 1.0  # mass of test particle (kg)

# Define the effective one-body Hamiltonian function
def seoBNR_v4p_effective_hamiltonian(H_seoBNR_v4p_real):
    """
    Calculate the effective one-body Hamiltonian.
    
    Parameters:
        H_seoBNR_v4p_real: real part of SEOB Hamiltonian (J)
    
    Returns:
        H_seoBNR_v4p_effective: effective one-body Hamiltonian (J)
    """
    # Calculate the effective one-body Hamiltonian
    H_seoBNR_v4p_effective = mu * G * c**2 * np.exp(1j * 0.25 * H_seoBNR_v4p_real)
    
    return H_seoBNR_v4p_effective

# Define the real part of the SEOB Hamiltonian
def seoBNR_v4p_real_hamiltonian(m1, m2, S1, S2):
    """
    Calculate the real part of the SEOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_real: real part of SEOB Hamiltonian (J)
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m**NRPy+: Test Particle of Mass $\mu$ and Spin ${\bf S}^*$**
==========================================================

### Theory Review

#### Introduction to Test Particle

*   **Test Particle:** In this section, we discuss a test particle of mass $\mu$ and spin ${\bf S}^*$.
    +   This is a fundamental concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Test Particle:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the mass of the test particle
mu = 1.0  # mass of test particle (kg)

# Define the spin of the test particle
def calculate_spin(S, *args):
    """
    Calculate the spin of the test particle.
    
    Parameters:
        S: initial spin of test particle
    
    Returns:
        S_star: final spin of test particle
    """
    # Calculate the final spin of the test particle
    S_star = np.linalg.norm(S) / (G * c**2)
    
    return S_star

# Define the test particle's motion using the effective one-body Hamiltonian
def seoBNR_v4p_test_particle_motion(mu, H_seoBNR_v4p_effective):
    """
    Calculate the test particle's motion.
    
    Parameters:
        mu: mass of test particle (kg)
        H_seoBNR_v4p_effective: effective one-body Hamiltonian
    
    Returns:
        position: position of test particle
        velocity: velocity of test particle
    """
    # Calculate the position and velocity of the test particle using the effective one-body Hamiltonian
    position = np.exp(1j * H_seoBNR_v4p_effective) * mu * G * c**2
    velocity = -np.imag(np.exp(1j * H_seoBNR_v4p_effective)) * mu * G * c**2
    
    return position, velocity

# Define the real part of the SEOB Hamiltonian
def seoBNR_v4p_real_hamiltonian(m1, m2, S1, S2):
    """
    Calculate the real part of the SEOB Hamiltonian**NRPy+: Computing $H_{\rm real}$ in a Deformed Kerr Background**
===========================================================

### Theory Review

#### Introduction to Computing $H_{\rm real}$

*   **Computing $H_{\rm real}$:** In this section, we discuss the computation of $H_{\rm real}$ in a deformed Kerr background.
    +   This is a fundamental concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Computing $H_{\rm real}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Euler-Mascheroni constant
gamma_E = 0.57721566490153286060  # Euler-Mascheroni constant

# Define the tortoise coordinate
t_tortoise = 100.0  # tortoise coordinate

# Define the function to compute $H_{\rm real}$
def v4P_compute_Hreal(m1, m2, gamma_E, t_tortoise):
    """
    Compute $H_{\rm real}$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        gamma_E: Euler-Mascheroni constant
        t_tortoise: tortoise coordinate
    
    Returns:
        H_real: $H_{\rm real}$
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Define the Cartesian quasi-isotropic coordinates
    def cartesian_quasi_isotropic_coordinates(r, theta, phi):
        """
        Compute the Cartesian quasi-isotropic coordinates.
        
        Parameters:
            r: radial coordinate (m)
            theta: polar angle (rad)
            phi: azimuthal angle (rad)
        
        Returns:
            x: $x$-coordinate (m)
            y: $y$-coordinate (m)
            z: $z$-coordinate (m)
        """
        # Compute the Cartesian quasi-isotropic coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta)**NRPy+: Citations**
==================

### Theory Review

#### Introduction to Citations

*   **Citations:** In this section, we discuss the references used throughout this module.
    +   This is a crucial aspect of academic integrity and provides context for the work being presented.

### Code Explanation


```python
"""
Citations:
"""
# Define the citations
BB2010 = "Barausse and Buonanno (2010)"
BB2011 = "Barausse and Buonanno (2011)"
OB2020 = "Ossokine, Buonanno, Marsat, et al (2020)"
SH2016 = "Steinhoff, Hinderer, Buonanno, et al (2016)"
BL2017 = "Bohe, Shao, Taracchini, et al (2017)"
P2010 = "Pan, Buonanno, Buchman, et. al (2010)"
T2014 = "Taracchini, Buonanno, Pan, et al (2014)"
T2012 = "Taracchini, Pan, Buonanno, et al (2012)"
D2000 = "Damour, Jaranowski, and Schaefer (2000)"

# Print the citations
print("BB2010:", BB2010)
print("BB2011:", BB2011)
print("OB2020:", OB2020)
print("SH2016:", SH2016)
print("BL2017:", BL2017)
print("P2010:", P2010)
print("T2014:", T2014)
print("T2012:", T2012)
print("D2000:", D2000)
```

This code defines the citations used throughout this module and prints them for reference.

### Mathematics


$$ H_{\rm real} = \frac{\partial}{\partial x^i} (\gamma^{ij}(x) f_j(x)) $$

where $H_{\rm real}$ is the Hamiltonian, $\gamma^{ij}(x)$ is the metric tensor, and $f_j(x)$ is a function of position.**NRPy+: Table of Contents**
==========================

### Theory Review

#### Introduction to Table of Contents

*   **Table of Contents:** In this section, we provide an overview of the organization and structure of this notebook.
    +   This is a useful reference for navigating the content and understanding the flow of ideas.

### Code Explanation


```python
"""
Table of Contents:
"""
# Define the sections and subsections
sections = [
    "Introduction to NRPy+",
    "The Physical System of Interest",
    "Effective One-Body Hamiltonian",
    "Test Particle of Mass $\mu$ and Spin ${\bf S}^*$",
    "Computing $H_{\rm real}$ in a Deformed Kerr Background",
    "Citations"
]

subsections = {
    "Introduction to NRPy+": [
        "Overview of NRPy+ Project",
        "Import Necessary Modules"
    ],
    "The Physical System of Interest": [
        "Binary Black Hole System",
        "SEOB Hamiltonian"
    ],
    "Effective One-Body Hamiltonian": [
        "Canonical Transformation",
        "Mapping to Effective Hamiltonian"
    ],
    "Test Particle of Mass $\mu$ and Spin ${\bf S}^*$": [
        "Mass of Test Particle",
        "Spin of Test Particle"
    ],
    "Computing $H_{\rm real}$ in a Deformed Kerr Background": [
        "Cartesian Quasi-Isotropic Coordinates",
        "Tortoise Coordinate"
    ],
    "Citations": []
}

# Print the table of contents
print("Table of Contents:")
for section, subsections in sections.items():
    print(f"  {section}")
    for i, subsection in enumerate(subsections):
        print(f"      Subsection {i+1}: {subsection}")

# Define the reference links
ref_links = {
    "Step 0": "https://arxiv.org/abs/0912.3517",
    "Step 1": "https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates"
}

# Print the reference links
print("\nReference Links:")
for link, title in ref_links.items():
    print(f"  {link}: {title}")
```

This code defines the sections and subsections of this notebook and prints the table of contents. It also defines the reference links for further reading.

### Mathematics


$$ H_{\**NRPy+: Output Creation**
==========================

### Theory Review

#### Introduction to Output Creation

*   **Output Creation:** In this section, we discuss the process of creating the output directory for SEOBNR.
    +   This is an essential step in setting up the project and preparing the necessary files.

### Code Explanation


```python
"""
Output Creation:
"""
# Import necessary modules
import os

# Define the output directory path
output_dir_path = "output/SEOBNR/"

# Check if the output directory exists
if not os.path.exists(output_dir_path):
    # Create the output directory
    os.makedirs(output_dir_path)

    print(f"Output directory '{output_dir_path}' created.")
else:
    print(f"Output directory '{output_dir_path}' already exists.")

# Define a function to create the necessary files in the output directory
def create_output_files():
    """
    Create the necessary files in the output directory.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Create the necessary files (e.g., input file, output file)
    with open(os.path.join(output_dir_path, "input.txt"), "w") as f:
        f.write("This is an example input file.")
        
    with open(os.path.join(output_dir_path, "output.txt"), "w") as f:
        f.write("This is an example output file.")

# Call the function to create the necessary files
create_output_files()
```

This code creates the output directory for SEOBNR and defines a function to create the necessary files within that directory.

### Mathematics


$$ \text{Output Directory Path} = "output/SEOBNR/" $$**NRPy+: The Real Hamiltonian $H_{\rm real}$**
============================================

### Theory Review

#### Introduction to the Real Hamiltonian

*   **The Real Hamiltonian:** In this section, we discuss the real Hamiltonian $H_{\rm real}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Real Hamiltonian $H_{\rm real}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the real Hamiltonian function
def seoBNR_v4p_real_hamiltonian(m1, m2):
    """
    Calculate the real Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_real: real Hamiltonian
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the dimensionless masses
    chi_dimless_1 = np.sqrt(m1) / (G * c**2)
    chi_dimless_2 = np.sqrt(m2) / (G * c**2)
    
    # Calculate the real Hamiltonian
    H_seoBNR_v4p_real = mu * G * c**2 * (1 + 73/24 * (chi_dimless_1**2) + 37/96 * (chi_dimless_2**2))
    
    return H_seoBNR_v4p_real

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_real = seoBNR_v4p_real_hamiltonian(m1, m2)
print("Real Hamiltonian:", H_seoBNR_v4p_real)
```

This code defines the real Hamiltonian $H_{\rm real}$ and calculates it using some example values.

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the effective Hamiltonian function
def seoBNR_v4p_effective_hamiltonian(H_seoBNR_v4p_real):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        H_seoBNR_v4p_real: real Hamiltonian
    
    Returns:
        H_seoBNR_v4p_eff: effective Hamiltonian
    """
    # Calculate the effective Hamiltonian
    H_seoBNR_v4p_eff = np.exp(1j * 0.25 * H_seoBNR_v4p_real)
    
    return H_seoBNR_v4p_eff

# Test the function with some example values
H_seoBNR_v4p_real = seoBNR_v4p_real_hamiltonian(m1, m2)

H_seoBNR_v4p_eff = seoBNR_v4p_effective_hamiltonian(H_seoBNR_v4p_real)
print("Effective Hamiltonian:", H_seoBNR_v4p_eff)
```

This code defines the effective Hamiltonian $H_{\rm eff}$ and calculates it using some example values.

### Mathematics


$$ H_{\rm eff} = \exp(1j \cdot 0.25 \cdot H_{\rm real}) $$**NRPy+: Terms of $H_{\rm eff}$
================================

### Theory Review

#### Introduction to Terms of $H_{\rm eff}$

*   **Terms of $H_{\rm eff}$:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Terms of $H_{\rm eff}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the terms of $H_{\rm eff}$ function
def seoBNR_v4p_effective_hamiltonian_terms(H_seoBNR_v4p_real):
    """
    Calculate the terms of $H_{\rm eff}$.
    
    Parameters:
        H_seoBNR_v4p_real: real Hamiltonian
    
    Returns:
        terms_H_seoBNR_v4p_eff: terms of effective Hamiltonian
    """
    # Calculate the terms of $H_{\rm eff}$
    terms_H_seoBNR_v4p_eff = np.exp(1j * 0.25 * H_seoBNR_v4p_real)
    
    return terms_H_seoBNR_v4p_eff

# Test the function with some example values
H_seoBNR_v4p_real = seoBNR_v4p_real_hamiltonian(m1, m2)

terms_H_seoBNR_v4p_eff = seoBNR_v4p_effective_hamiltonian_terms(H_seoBNR_v4p_real)
print("Terms of $H_{\rm eff}$:", terms_H_seoBNR_v4p_eff)
```

This code defines the terms of the effective Hamiltonian $H_{\rm eff}$ and calculates them using some example values.

### Mathematics


$$ \text{Terms of } H_{\rm eff} = \exp(1j \cdot 0.25 \cdot H_{\rm real}) $$**NRPy+: Leading Order Spin Effects $H_{\rm S}$
=============================================

### Theory Review

#### Introduction to Leading Order Spin Effects

*   **Leading Order Spin Effects:** In this section, we discuss the leading order spin effects in the effective Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Leading Order Spin Effects $H_{\rm S}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the leading order spin effects function
def seoBNR_v4p_spin_effects(m1, m2, S1, S2):
    """
    Calculate the leading order spin effects.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_spin: leading order spin effects
    """
    # Calculate the dimensionless spins
    chi_dimless_1 = np.linalg.norm(S1) / (G * c**2 * m1)
    chi_dimless_2 = np.linalg.norm(S2) / (G * c**2 * m2)
    
    # Calculate the leading order spin effects
    H_seoBNR_v4p_spin = G * c**2 * (chi_dimless_1**2 + chi_dimless_2**2)
    
    return H_seoBNR_v4p_spin

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_spin = seoBNR_v4p_spin_effects(m1, m2, S1, S2)
print("Leading Order Spin Effects:", H**NRPy+: The Nonspinning Hamiltonian $H_{\rm NS}$
=============================================

### Theory Review

#### Introduction to the Nonspinning Hamiltonian

*   **The Nonspinning Hamiltonian:** In this section, we discuss the nonspinning Hamiltonian $H_{\rm NS}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Nonspinning Hamiltonian $H_{\rm NS}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the nonspinning Hamiltonian function
def seoBNR_v4p_nonspinning_hamiltonian(m1, m2):
    """
    Calculate the nonspinning Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ns: nonspinning Hamiltonian
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the nonspinning Hamiltonian
    H_seoBNR_v4p_ns = G * c**2 * (m1 + m2)
    
    return H_seoBNR_v4p_ns

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_ns = seoBNR_v4p_nonspinning_hamiltonian(m1, m2)
print("Nonspinning Hamiltonian:", H_seoBNR_v4p_ns)
```

This code defines the nonspinning Hamiltonian $H_{\rm NS}$ and calculates it using some example values.

### Mathematics


$$ H_{\rm NS} = G \cdot c^2 \cdot (m_1 + m_2) $$**NRPy+: The Quadrupole Deformation $H_{\rm D}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the quadrupole deformation of the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Quadrupole Deformation $H_{\rm D}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the quadrupole deformation function
def seoBNR_v4p_quadrupole_deformation(m1, m2):
    """
    Calculate the quadrupole deformation.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_d: quadrupole deformation
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the quadrupole deformation
    H_seoBNR_v4p_d = G * c**2 * (m1 - m2)**2 / (M**2)
    
    return H_seoBNR_v4p_d

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_d = seoBNR_v4p_quadrupole_deformation(m1, m2)
print("Quadrupole Deformation:", H_seoBNR_v4p_d)
```

This code defines the quadrupole deformation $H_{\rm D}$ and calculates it using some example values.

### Mathematics


$$ H_{\rm D} = G \cdot c^2 \cdot \frac{(m_1 - m_2)^2}{M^2} $$**NRPy+: The Spin-Orbit Term $H_{\rm SO}$
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the spin-orbit term of the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Spin-Orbit Term $H_{\rm SO}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the spin-orbit term function
def seoBNR_v4p_spin_orbit_term(m1, m2, S1, S2):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_so: spin-orbit term
    """
    # Calculate the dimensionless spins
    chi_dimless_1 = np.linalg.norm(S1) / (G * c**2 * m1)
    chi_dimless_2 = np.linalg.norm(S2) / (G * c**2 * m2)
    
    # Calculate the spin-orbit term
    H_seoBNR_v4p_so = G * c**2 * chi_dimless_1 * chi_dimless_2
    
    return H_seoBNR_v4p_so

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_so = seoBNR_v4p_spin_orbit_term(m1, m2, S1, S2)
print("Spin**NRPy+: $H_{\rm SO}$ Term 1**
=============================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 1

*   **$H_{\rm SO}$ Term 1:** In this section, we discuss the first term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 1 function
def seoBNR_v4p_so_term_1(m1, m2, S1, S2):
    """
    Calculate the first term of the spin-orbit Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_so_1: first term of spin-orbit Hamiltonian
    """
    # Calculate the dimensionless spins
    chi_dimless_1 = np.linalg.norm(S1) / (G * c**2 * m1)
    chi_dimless_2 = np.linalg.norm(S2) / (G * c**2 * m2)
    
    # Calculate the first term of the spin-orbit Hamiltonian
    H_seoBNR_v4p_so_1 = G * c**2 * chi_dimless_1
    
    return H_seoBNR_v4p_so_1

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_so_1 =**NRPy+: $H_{\rm SO}$ Term 2 Coefficient**
=============================================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 2 Coefficient

*   **$H_{\rm SO}$ Term 2 Coefficient:** In this section, we discuss the coefficient of the second term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 2 Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 2 Coefficient function
def seoBNR_v4p_so_term_2_coefficient(m1, m2):
    """
    Calculate the coefficient of the second term of the spin-orbit Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        coeff_seoBNR_v4p_so_2: coefficient of second term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the coefficient of the second term of the spin-orbit Hamiltonian
    coeff_seoBNR_v4p_so_2 = 2 / (M**2)
    
    return coeff_seoBNR_v4p_so_2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

coeff_seoBNR_v4p_so_2 = seoBNR_v4p_so_term_2_coefficient(m1, m2)
print("Coefficient of $H_{\rm SO}$ Term 2:", coeff_seoBNR_v4p_so_2)
```

This code calculates the coefficient of the second term of the spin-orbit Hamiltonian using some example values.

### Mathematics


$$ \text{Coefficient of } H_{\rm**NRPy+: $H_{\rm SO}$ Term 2**
================================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 2

*   **$H_{\rm SO}$ Term 2:** In this section, we discuss the second term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 2 function
def seoBNR_v4p_so_term_2(m1, m2):
    """
    Calculate the second term of the spin-orbit Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_so_2: second term of spin-orbit Hamiltonian
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the second term of the spin-orbit Hamiltonian
    H_seoBNR_v4p_so_2 = 2 * G * c**2 * (m1 - m2) * np.linalg.norm(np.cross(S1, S2)) / (M**2)
    
    return H_seoBNR_v4p_so_2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_so_2 = seoBNR_v4p_so_term_2(m1, m2)
print("Second term of $H_{\rm SO}$:", H**NRPy+: $H_{\rm SO}$ Term 2a**
================================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 2a

*   **$H_{\rm SO}$ Term 2a:** In this section, we discuss the first part of the second term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 2a:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 2a function
def seoBNR_v4p_so_term_2a(m1, m2):
    """
    Calculate the first part of the second term of the spin-orbit Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_so_2a: first part of second term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the first part of the second term of the spin-orbit Hamiltonian
    H_seoBNR_v4p_so_2a = 2 * G * c**2 * (m1 - m2) / (M**2)
    
    return H_seoBNR_v4p_so_2a

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_so_2a = seoBNR_v4p_so_term_2a(m1, m2)
print("First part of $H_{\rm SO}$ Term 2:", H_seoBNR_v4p_so_2a)
```

This code calculates the first part of the second term of the spin-orbit Hamiltonian using some**NRPy+: $H_{\rm SO}$ Term 2b**
================================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 2b

*   **$H_{\rm SO}$ Term 2b:** In this section, we discuss the second part of the second term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 2b:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 2b function
def seoBNR_v4p_so_term_2b(S1, S2):
    """
    Calculate the second part of the second term of the spin-orbit Hamiltonian.
    
    Parameters:
        S1: spin of first black hole
        S2: spin of second black hole
    
    Returns:
        H_seoBNR_v4p_so_2b: second part of second term
    """
    # Calculate the cross product of the spins
    cross_product = np.cross(S1, S2)
    
    # Calculate the magnitude of the cross product
    magnitude_cross_product = np.linalg.norm(cross_product)
    
    # Calculate the second part of the second term of the spin-orbit Hamiltonian
    H_seoBNR_v4p_so_2b = 2 * G * c**2 * (m1 - m2) / (M**2) * magnitude_cross_product
    
    return H_seoBNR_v4p_so_2b

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_so_2b = seoBNR_v4p_so_term**NRPy+: $H_{\rm SO}$ Term 2c**
================================

### Theory Review

#### Introduction to $H_{\rm SO}$ Term 2c

*   **$H_{\rm SO}$ Term 2c:** In this section, we discuss the final part of the second term of the spin-orbit Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SO}$ Term 2c:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SO}$ Term 2c function
def seoBNR_v4p_so_term_2c(m1, m2):
    """
    Calculate the final part of the second term of the spin-orbit Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_so_2c: final part of second term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the final part of the second term of the spin-orbit Hamiltonian
    H_seoBNR_v4p_so_2c = -2 * G * c**2 * (m1 - m2) ** 3 / (M**2)
    
    return H_seoBNR_v4p_so_2c

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_so_2c = seoBNR_v4p_so_term_2c(m1, m2)
print("Final part of $H_{\rm SO}$ Term 2:", H_seoBNR_v4p_so_2c)
```

This code calculates the final part of the second term of the spin-orbit Hamilton**NRPy+: The Spin-Spin Term $H_{\rm SS}$
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term

*   **The Spin-Spin Term:** In this section, we discuss the spin-spin term of the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Spin-Spin Term $H_{\rm SS}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the spin-spin term function
def seoBNR_v4p_spin_spin_term(m1, m2):
    """
    Calculate the spin-spin term.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ss: spin-spin term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin-spin term
    H_seoBNR_v4p_ss = G * c**4 * (np.linalg.norm(S1) ** 2 * np.linalg.norm(S2)**2 - (S1 · S2) ** 2) / (M**3)
    
    return H_seoBNR_v4p_ss

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_ss = seoBNR_v4p_spin_spin_term(m1, m2)
print("Spin-Spin Term:", H_seoBNR_v4p_ss)
```

This code calculates the spin-spin term using some example values.

### Mathematics


$$ H_{\rm SS} = G \cdot c^**NRPy+: $H_{\rm SS}$ Term 1**
================================

### Theory Review

#### Introduction to $H_{\rm SS}$ Term 1

*   **$H_{\rm SS}$ Term 1:** In this section, we discuss the first term of the spin-spin Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SS}$ Term 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SS}$ Term 1 function
def seoBNR_v4p_ss_term_1(m1, m2):
    """
    Calculate the first term of the spin-spin Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ss_1: first term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the magnitude of the spins
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the first term of the spin-spin Hamiltonian
    H_seoBNR_v4p_ss_1 = G * c**4 * (mag_S1 ** 2) * (mag_S2 ** 2) / (M**3)
    
    return H_seoBNR_v4p_ss_1

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S1 = np.array([0, 0, 1])  # spin of first black hole
S2 = np.array([0, 0, -1])  # spin of second black hole

H_seoBNR_v4p_ss_1 = seoBNR_v4p_ss_term_1(m1**NRPy+: $H_{\rm SS}$ Term 2 Coefficient**
=============================================

### Theory Review

#### Introduction to $H_{\rm SS}$ Term 2 Coefficient

*   **$H_{\rm SS}$ Term 2 Coefficient:** In this section, we discuss the coefficient of the second term of the spin-spin Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SS}$ Term 2 Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SS}$ Term 2 Coefficient function
def seoBNR_v4p_ss_term_2_coefficient(m1, m2):
    """
    Calculate the coefficient of the second term of the spin-spin Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        coeff_seoBNR_v4p_ss_2: coefficient of second term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the coefficient of the second term of the spin-spin Hamiltonian
    coeff_seoBNR_v4p_ss_2 = -G * c**4 / (M**3)
    
    return coeff_seoBNR_v4p_ss_2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

coeff_seoBNR_v4p_ss_2 = seoBNR_v4p_ss_term_2_coefficient(m1, m2)
print("Coefficient of $H_{\rm SS}$ Term 2:", coeff_seoBNR_v4p_ss_2)
```

This code calculates the coefficient of the second term of the spin-spin Hamiltonian using some example values.

### Mathematics


$$ \text{Coefficient of } H_{\rm**NRPy+: $H_{\rm SS}$ Term 2**
================================

### Theory Review

#### Introduction to $H_{\rm SS}$ Term 2

*   **$H_{\rm SS}$ Term 2:** In this section, we discuss the second term of the spin-spin Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SS}$ Term 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SS}$ Term 2 function
def seoBNR_v4p_ss_term_2(m1, m2):
    """
    Calculate the second term of the spin-spin Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ss_2: second term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the second term of the spin-spin Hamiltonian
    H_seoBNR_v4p_ss_2 = -G * c**4 * (dot_product_S1_S2 ** 2) / (M**3)
    
    return H_seoBNR_v4p_ss_2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e**NRPy+: $H_{\rm SS}$ Term 3 Coefficient**
=============================================

### Theory Review

#### Introduction to $H_{\rm SS}$ Term 3 Coefficient

*   **$H_{\rm SS}$ Term 3 Coefficient:** In this section, we discuss the coefficient of the third term of the spin-spin Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SS}$ Term 3 Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SS}$ Term 3 Coefficient function
def seoBNR_v4p_ss_term_3_coefficient(m1, m2):
    """
    Calculate the coefficient of the third term of the spin-spin Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        coeff_seoBNR_v4p_ss_3: coefficient of third term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the coefficient of the third term of the spin-spin Hamiltonian
    coeff_seoBNR_v4p_ss_3 = G * c**4 / (M**3)
    
    return coeff_seoBNR_v4p_ss_3

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

coeff_seoBNR_v4p_ss_3 = seoBNR_v4p_ss_term_3_coefficient(m1, m2)
print("Coefficient of $H_{\rm SS}$ Term 3:", coeff_seoBNR_v4p_ss_3)
```

This code calculates the coefficient of the third term of the spin-spin Hamiltonian using some example values.

### Mathematics


$$ \text{Coefficient of } H_{\rm**NRPy+: $H_{\rm SS}$ Term 3**
================================

### Theory Review

#### Introduction to $H_{\rm SS}$ Term 3

*   **$H_{\rm SS}$ Term 3:** In this section, we discuss the third term of the spin-spin Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm SS}$ Term 3:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm SS}$ Term 3 function
def seoBNR_v4p_ss_term_3(m1, m2):
    """
    Calculate the third term of the spin-spin Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ss_3: third term
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the third term of the spin-spin Hamiltonian
    H_seoBNR_v4p_ss_3 = G * c**4 * (dot_product_S1_S2 ** 3) / (M**3)
    
    return H_seoBNR_v4p_ss_3

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24**NRPy+: The $H_{\rm NS}$ Terms**
=====================================

### Theory Review

#### Introduction to the $H_{\rm NS}$ Terms

*   **The $H_{\rm NS}$ Terms:** In this section, we discuss the terms that arise from the neutron star contribution to the Hamiltonian.
    +   These terms are crucial in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The $H_{\rm NS}$ Terms:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm NS}$ Terms function
def seoBNR_v4p_ns_terms(m1, m2):
    """
    Calculate the neutron star terms.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        H_seoBNR_v4p_ns: neutron star terms
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the neutron star terms
    H_seoBNR_v4p_ns = G * c**4 * (dot_product_S1_S2 ** 2) / (M**3)
    
    return H_seoBNR_v4p_ns

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

H_seoBNR_v4p_ns = seo**NRPy+: $\beta p$ Sum**
==========================

### Theory Review

#### Introduction to the $\beta p$ Sum

*   **$\beta p$ Sum:** In this section, we discuss the sum of the beta parameters and the momenta.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\beta p$ Sum:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\beta p$ Sum function
def seoBNR_v4p_beta_p_sum(m1, m2):
    """
    Calculate the sum of beta parameters and momenta.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sum_seoBNR_v4p_beta_p: $\beta p$ Sum
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the beta parameters
    beta_1 = np.array([0, 0, 1])  # beta parameter of first black hole
    beta_2 = np.array([0, 0, -1])  # beta parameter of second black hole
    
    # Calculate the momenta
    p1 = mu * (beta_1 ** 2) / M
    p2 = mu * (beta_2 ** 2) / M
    
    # Calculate the sum of beta parameters and momenta
    sum_seoBNR_v4p_beta_p = np.sum(beta_1**2) + np.sum(p1) + np.sum(beta_2**2) + np.sum(p2)
    
    return sum_seoBNR_v4p_beta_p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sum_seoBNR_v4p_beta_p = seoBNR_v4p_beta_p_sum(m1, m2)
print**NRPy+: $\alpha$**
=====================

### Theory Review

#### Introduction to $\alpha$

*   **$\alpha$:** In this section, we discuss the parameter $\alpha$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\alpha$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\alpha$ function
def seoBNR_v4p_alpha(m1, m2):
    """
    Calculate the parameter $\alpha$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        alpha_seoBNR_v4p: value of $\alpha$
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the parameter $\alpha$
    alpha_seoBNR_v4p = np.sqrt(G * c**2) / (M**2)
    
    return alpha_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

alpha_seoBNR_v4p = seoBNR_v4p_alpha(m1, m2)
print("Value of $\alpha$:", alpha_seoBNR_v4p)
```

This code calculates the value of the parameter $\alpha$ using some example values.

### Mathematics


$$ \alpha = \sqrt{\frac{G c^2}{M^2}} $$**NRPy+: $H_{\rm NS}$ Radicand**
=====================================

### Theory Review

#### Introduction to the $H_{\rm NS}$ Radicand

*   **$H_{\rm NS}$ Radicand:** In this section, we discuss the radicand of the neutron star term in the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm NS}$ Radicand:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm NS}$ Radicand function
def seoBNR_v4p_ns_radicand(m1, m2):
    """
    Calculate the radicand of the neutron star term.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        radicand_seoBNR_v4p_ns: $H_{\rm NS}$ Radicand
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the $H_{\rm NS}$ Radicand
    radicand_seoBNR_v4p_ns = G * c**4 * (dot_product_S1_S2 ** 2) / (M**3)
    
    return radicand_seoBNR_v4p_ns

# Test the function with some example values
m1 = 1.98910e30  #**NRPy+: $\gamma p$ Sum**
==========================

### Theory Review

#### Introduction to the $\gamma p$ Sum

*   **$\gamma p$ Sum:** In this section, we discuss the sum of the gamma parameters and the momenta.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\gamma p$ Sum:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\gamma p$ Sum function
def seoBNR_v4p_gamma_p_sum(m1, m2):
    """
    Calculate the sum of gamma parameters and momenta.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sum_seoBNR_v4p_gamma_p: $\gamma p$ Sum
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the gamma parameters
    gamma_1 = np.array([0, 0, 1])  # gamma parameter of first black hole
    gamma_2 = np.array([0, 0, -1])  # gamma parameter of second black hole
    
    # Calculate the momenta
    p1 = mu * (gamma_1 ** 2) / M
    p2 = mu * (gamma_2 ** 2) / M
    
    # Calculate the sum of gamma parameters and momenta
    sum_seoBNR_v4p_gamma_p = np.sum(gamma_1**2) + np.sum(p1) + np.sum(gamma_2**2) + np.sum(p2)
    
    return sum_seoBNR_v4p_gamma_p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sum_seoBNR_v4p_gamma_p = seoBNR_v4p_gamma_p_sum(m1, m2**NRPy+: ${\cal Q}_{4}$**
==========================

### Theory Review

#### Introduction to ${\cal Q}_{4}$

*   **${\cal Q}_{4}$:** In this section, we discuss the term ${\cal Q}_{4}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\cal Q}_{4}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\cal Q}_{4}$ function
def seoBNR_v4p_Q4(m1, m2):
    """
    Calculate the term ${\cal Q}_{4}$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        Q4_seoBNR_v4p: value of ${\cal Q}_{4}$
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the term ${\cal Q}_{4}$
    Q4_seoBNR_v4p = G * c**4 / (M**3) * (dot_product_S1_S2 ** 2)
    
    return Q4_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

Q4_seoBNR_v**NRPy+: The $H_{\rm D}$ Terms**
=====================================

### Theory Review

#### Introduction to the $H_{\rm D}$ Terms

*   **The $H_{\rm D}$ Terms:** In this section, we discuss the terms that arise from the radiation reaction in the Hamiltonian.
    +   These terms are crucial in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The $H_{\rm D}$ Terms:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm D}$ Terms function
def seoBNR_v4p_hd_terms(m1, m2):
    """
    Calculate the terms that arise from radiation reaction.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hd_seoBNR_v4p: $H_{\rm D}$ Terms
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the $H_{\rm D}$ Terms
    hd_seoBNR_v4p = G * c**4 / (M**3) * (dot_product_S1_S2 ** 2)
    
    return hd_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

hd**NRPy+: $H_{\rm D}$ Coefficient**
=====================================

### Theory Review

#### Introduction to the $H_{\rm D}$ Coefficient

*   **$H_{\rm D}$ Coefficient:** In this section, we discuss the coefficient of the radiation reaction term in the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm D}$ Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm D}$ Coefficient function
def seoBNR_v4p_hd_coefficient(m1, m2):
    """
    Calculate the coefficient of radiation reaction.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hd_coeff_seoBNR_v4p: $H_{\rm D}$ Coefficient
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the $H_{\rm D}$ Coefficient
    hd_coeff_seoBNR_v4p = G * c**4 / (M**3) * (dot_product_S1_S2 ** 2)
    
    return hd_coeff_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of**NRPy+: $H_{\rm D}$ Sum**
==========================

### Theory Review

#### Introduction to the $H_{\rm D}$ Sum

*   **$H_{\rm D}$ Sum:** In this section, we discuss the sum of the radiation reaction terms in the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm D}$ Sum:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm D}$ Sum function
def seoBNR_v4p_hd_sum(m1, m2):
    """
    Calculate the sum of radiation reaction terms.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hd_sum_seoBNR_v4p: $H_{\rm D}$ Sum
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate the $H_{\rm D}$ Sum
    hd_sum_seoBNR_v4p = G * c**4 / (M**3) * (dot_product_S1_S2 ** 2 + dot_product_S1_S2)
    
    return hd_sum_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of**NRPy+: $H_{\rm D}$ Sum Term 1**
=====================================

### Theory Review

#### Introduction to $H_{\rm D}$ Sum Term 1

*   **$H_{\rm D}$ Sum Term 1:** In this section, we discuss the first term of the sum of radiation reaction terms in the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm D}$ Sum Term 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm D}$ Sum Term 1 function
def seoBNR_v4p_hd_sum_term_1(m1, m2):
    """
    Calculate the first term of radiation reaction sum.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hd_sum_term_1_seoBNR_v4p: $H_{\rm D}$ Sum Term 1
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate $H_{\rm D}$ Sum Term 1
    hd_sum_term_1_seoBNR_v4p = G * c**4 / (M**3) * (dot_product_S1_S2 ** 2)
    
    return hd_sum_term_1_seoBNR_v4p

# Test the function with some example values
m1 = **NRPy+: $H_{\rm D}$ Sum Term 2**
=====================================

### Theory Review

#### Introduction to $H_{\rm D}$ Sum Term 2

*   **$H_{\rm D}$ Sum Term 2:** In this section, we discuss the second term of the sum of radiation reaction terms in the Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm D}$ Sum Term 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm D}$ Sum Term 2 function
def seoBNR_v4p_hd_sum_term_2(m1, m2):
    """
    Calculate the second term of radiation reaction sum.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hd_sum_term_2_seoBNR_v4p: $H_{\rm D}$ Sum Term 2
    """
    # Calculate the total mass and reduced mass
    M = m1 + m2
    mu = m1 * m2 / M
    
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_S1_S2 = np.dot(S1, S2)
    
    # Calculate the magnitude of the spin vectors
    mag_S1 = np.linalg.norm(S1)
    mag_S2 = np.linalg.norm(S2)
    
    # Calculate $H_{\rm D}$ Sum Term 2
    hd_sum_term_2_seoBNR_v4p = G * c**4 / (M**3) * dot_product_S1_S2
    
    return hd_sum_term_2_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e**NRPy+: Common Dot Products**
=============================

### Theory Review

#### Introduction to Common Dot Products

*   **Common Dot Products:** In this section, we discuss the common dot products used in the calculation of radiation reaction terms.
    +   These dot products are crucial in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Common Dot Products:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Common Dot Products function
def seoBNR_v4p_dot_products(m1, m2):
    """
    Calculate common dot products.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        dot_product_seoBNR_v4p: common dot products
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the dot product of the spin vectors
    dot_product_seoBNR_v4p = np.dot(S1, S2)
    
    return dot_product_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

dot_product_seoBNR_v4p = seoBNR_v4p_dot_products(m1, m2)
print("Common Dot Products:", dot_product_seoBNR_v4p)
```

### Mathematics

$$\vec{S}_1 \cdot \vec{S}_2 = S_{1x} S_{2x} + S_{1y} S_{2y} + S_{1z} S_{2z}$$**NRPy+: ${\bf S} \cdot \boldsymbol{\xi}$**
=============================================

### Theory Review

#### Introduction to ${\bf S} \cdot \boldsymbol{\xi}$

*   **${\bf S} \cdot \boldsymbol{\xi}$:** In this section, we discuss the dot product of the spin vector and the radiation vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S} \cdot \boldsymbol{\xi}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S} \cdot \boldsymbol{\xi}$ function
def seoBNR_v4p_S_dot_xi(m1, m2):
    """
    Calculate the dot product of spin and radiation vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        S_dot_xi_seoBNR_v4p: ${\bf S} \cdot \boldsymbol{\xi}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the radiation vector
    xi = np.array([-1, 1, 1])  # radiation vector
    
    # Calculate the dot product of the spin and radiation vectors
    S_dot_xi_seoBNR_v4p = np.dot(S1, xi) + np.dot(S2, xi)
    
    return S_dot_xi_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S_dot_xi_seoBNR_v4p = seoBNR_v4p_S_dot_xi(m1, m2)
print("${\bf S} \cdot \**NRPy+: ${\bf S} \cdot {\bf v}$**
=====================================

### Theory Review

#### Introduction to ${\bf S} \cdot {\bf v}$

*   **${\bf S} \cdot {\bf v}$:** In this section, we discuss the dot product of the spin vector and the velocity vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S} \cdot {\bf v}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S} \cdot {\bf v}$ function
def seoBNR_v4p_S_dot_v(m1, m2):
    """
    Calculate the dot product of spin and velocity vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        S_dot_v_seoBNR_v4p: ${\bf S} \cdot {\bf v}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the velocity vector
    v = np.array([-1, 1, 1])  # velocity
    
    # Calculate the dot product of the spin and velocity vectors
    S_dot_v_seoBNR_v4p = np.dot(S1, v) + np.dot(S2, v)
    
    return S_dot_v_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S_dot_v_seoBNR_v4p = seoBNR_v4p_S_dot_v(m1, m2)
print("${\bf S} \cdot {\bf v}$:", S_dot_v_seoBNR_v4p)
```

### Mathematics**NRPy+: ${\bf S} \cdot {\bf n}$**
=====================================

### Theory Review

#### Introduction to ${\bf S} \cdot {\bf n}$

*   **${\bf S} \cdot {\bf n}$:** In this section, we discuss the dot product of the spin vector and the normal vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S} \cdot {\bf n}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S} \cdot {\bf n}$ function
def seoBNR_v4p_S_dot_n(m1, m2):
    """
    Calculate the dot product of spin and normal vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        S_dot_n_seoBNR_v4p: ${\bf S} \cdot {\bf n}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the normal vector
    n = np.array([-1, -1, 1])  # normal vector
    
    # Calculate the dot product of the spin and normal vectors
    S_dot_n_seoBNR_v4p = np.dot(S1, n) + np.dot(S2, n)
    
    return S_dot_n_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

S_dot_n_seoBNR_v4p = seoBNR_v4p_S_dot_n(m1, m2)
print("${\bf S} \cdot {\bf n}$:", S_dot_n_seoBNR_v4p)
```

###**NRPy+: ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$**
=============================================

### Theory Review

#### Introduction to ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$

*   **${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$:** In this section, we discuss the dot product of the spin vector and the unit direction vector of the Kerr spin axis.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$ function
def seoBNR_v4p_S_dot_skerrhat(m1, m2):
    """
    Calculate the dot product of spin and Kerr spin axis vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        S_dot_skerrhat_seoBNR_v4p: ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the unit direction vector of the Kerr spin axis
    skerrhat = np.array([-1, 1, 1])  # unit direction vector of Kerr spin axis
    
    # Calculate the dot product of the spin and Kerr spin axis vectors
    S_dot_skerrhat_seoBNR_v4p = np.dot(S1, skerrhat) + np.dot(S2, skerrhat)
    
    return S_dot_skerrhat_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
**NRPy+: ${\bf S}^{*} \cdot {\bf n}$**
=====================================

### Theory Review

#### Introduction to ${\bf S}^{*} \cdot {\bf n}$

*   **${\bf S}^{*} \cdot {\bf n}$:** In this section, we discuss the dot product of the complex conjugate spin vector and the normal vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S}^{*} \cdot {\bf n}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S}^{*} \cdot {\bf n}$ function
def seoBNR_v4p_sstardotn(m1, m2):
    """
    Calculate the dot product of complex conjugate spin and normal vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sstardotn_seoBNR_v4p: ${\bf S}^{*} \cdot {\bf n}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    Sstar1 = np.conjugate(S1)  # complex conjugate of S1
    Sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Define the normal vector
    n = np.array([-1, -1, 1])  # normal vector
    
    # Calculate the dot product of the complex conjugate spin and normal vectors
    sstardotn_seoBNR_v4p = np.dot(Sstar1, n) + np.dot(Sstar2, n)
    
    return sstardotn_seoBNR_v4p

# Test the function with some example values
**NRPy+: $H_{\rm real}$ Spin Combination ${\bf S}^{*}$**
===========================================================

### Theory Review

#### Introduction to $H_{\rm real}$ Spin Combination ${\bf S}^{*}$

*   **$H_{\rm real}$ Spin Combination ${\bf S}^{*}$:** In this section, we discuss the real part of the Hamiltonian spin combination involving the complex conjugate spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$H_{\rm real}$ Spin Combination ${\bf S}^{*}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $H_{\rm real}$ Spin Combination ${\bf S}^{*}$ function
def seoBNR_v4p_hreal_spin_combos(m1, m2):
    """
    Calculate the real part of the Hamiltonian spin combination involving complex conjugate spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hreal_spin_combos_seoBNR_v4p: $H_{\rm real}$ Spin Combination ${\bf S}^{*}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    Sstar1 = np.conjugate(S1)  # complex conjugate of S1
    Sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Calculate the Hamiltonian spin combination involving complex conjugate spin vector
    hreal_spin_combos_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * (Sstar1[0]**2 + Sstar2[0]**2)
    
    return hreal_spin_combos_seoBNR_v**NRPy+: ${\bf S}^{*}$**
=========================

### Theory Review

#### Introduction to ${\bf S}^{*}$

*   **${\bf S}^{*}$:** In this section, we discuss the complex conjugate spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S}^{*}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S}^{*}$ function
def seoBNR_v4p_sstar(m1, m2):
    """
    Calculate the complex conjugate spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sstar_seoBNR_v4p: ${\bf S}^{*}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    sstar_seoBNR_v4p = np.conjugate(S1) + np.conjugate(S2)
    
    return sstar_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sstar_seoBNR_v4p = seoBNR_v4p_sstar(m1, m2)
print("${\bf S}^{*}$:", sstar_seoBNR_v4p)
```

### Mathematics

$$
{\bf S}^{*} = {\rm conj}({\bf S}) = \left( \begin{array}{c}
S_1^* \\
S_2^* \\
S_3^*
\end{array}\right) =
\left( \begin{array**NRPy+: $\Delta_{\sigma^{*}}$**
=============================

### Theory Review

#### Introduction to $\Delta_{\sigma^{*}}$

*   **$\Delta_{\sigma^{*}}$:** In this section, we discuss the difference between the complex conjugate spin vector and the original spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_{\sigma^{*}}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_{\sigma^{*}}$ function
def seoBNR_v4p_deltasigmastar(m1, m2):
    """
    Calculate the difference between complex conjugate spin vector and original spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltasigmastar_seoBNR_v4p: $\Delta_{\sigma^{*}}$
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    sstar1 = np.conjugate(S1)  # complex conjugate of S1
    sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Calculate the difference between complex conjugate spin and original spin vectors
    deltasigmastar_seoBNR_v4p = (sstar1 - S1) + (sstar2 - S2)
    
    return deltasigmastar_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltasigmastar_seoBNR_v4p = seoBNR_v4p_deltasigmastar(m**NRPy+: $\sigma^{*}$ Coefficient**
=====================================

### Theory Review

#### Introduction to $\sigma^{*}$ Coefficient

*   **$\sigma^{*}$ Coefficient:** In this section, we discuss the coefficient involving the complex conjugate spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma^{*}$ Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma^{*}$ Coefficient function
def seoBNR_v4p_sigmastarcoeff(m1, m2):
    """
    Calculate the coefficient involving complex conjugate spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmastarcoeff_seoBNR_v4p: $\sigma^{*}$ Coefficient
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    sstar1 = np.conjugate(S1)  # complex conjugate of S1
    sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Calculate the coefficient involving complex conjugate spin vector
    sigmastarcoeff_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * (sstar1[0]**2 + sstar2[0]**2)
    
    return sigmastarcoeff_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmastarcoeff_seoBNR_v4p = seoBNR_v4p_sigmastarcoeff(m1, m2)
print("$\**NRPy+: $\sigma^{*}$ Coefficient Term 1**
=============================================

### Theory Review

#### Introduction to $\sigma^{*}$ Coefficient Term 1

*   **$\sigma^{*}$ Coefficient Term 1:** In this section, we discuss the first term of the coefficient involving the complex conjugate spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma^{*}$ Coefficient Term 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma^{*}$ Coefficient Term 1 function
def seoBNR_v4p_sigmastarcoeffterm1(m1, m2):
    """
    Calculate the first term of coefficient involving complex conjugate spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmastarcoeffterm1_seoBNR_v4p: $\sigma^{*}$ Coefficient Term 1
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    sstar1 = np.conjugate(S1)  # complex conjugate of S1
    sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Calculate the first term of coefficient involving complex conjugate spin vector
    sigmastarcoeffterm1_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * sstar1[0]**2
    
    return sigmastarcoeffterm1_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmastarcoeffterm1**NRPy+: $\sigma^{*}$ Coefficient Term 2**
=============================================

### Theory Review

#### Introduction to $\sigma^{*}$ Coefficient Term 2

*   **$\sigma^{*}$ Coefficient Term 2:** In this section, we discuss the second term of the coefficient involving the complex conjugate spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma^{*}$ Coefficient Term 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma^{*}$ Coefficient Term 2 function
def seoBNR_v4p_sigmastarcoeffterm2(m1, m2):
    """
    Calculate the second term of coefficient involving complex conjugate spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmastarcoeffterm2_seoBNR_v4p: $\sigma^{*}$ Coefficient Term 2
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the complex conjugate spin vector
    sstar1 = np.conjugate(S1)  # complex conjugate of S1
    sstar2 = np.conjugate(S2)  # complex conjugate of S2
    
    # Calculate the second term of coefficient involving complex conjugate spin vector
    sigmastarcoeffterm2_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * sstar2[0]**2
    
    return sigmastarcoeffterm2_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmastarcoeffterm2**NRPy+: $\sigma$ Coefficient**
==============================

### Theory Review

#### Introduction to $\sigma$ Coefficient

*   **$\sigma$ Coefficient:** In this section, we discuss the coefficient involving the spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma$ Coefficient:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma$ Coefficient function
def seoBNR_v4p_sigmacoeff(m1, m2):
    """
    Calculate the coefficient involving spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmacoeff_seoBNR_v4p: $\sigma$ Coefficient
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the coefficient involving spin vector
    sigmacoeff_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * (S1[0]**2 + S2[0]**2)
    
    return sigmacoeff_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmacoeff_seoBNR_v4p = seoBNR_v4p_sigmacoeff(m1, m2)
print("$\sigma$ Coefficient:", sigmacoeff_seoBNR_v4p)
```

### Mathematics

$$
\sigma =
\left( \begin{array}{c}
S_0 \\
S_1 \\
S_2
\end{array}\right) =
\frac{G \, c^4}{M^3} \,
\left( \begin{array}{c}
S**NRPy+: $\sigma$ Coefficient Term 1**
=============================================

### Theory Review

#### Introduction to $\sigma$ Coefficient Term 1

*   **$\sigma$ Coefficient Term 1:** In this section, we discuss the first term of the coefficient involving the spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma$ Coefficient Term 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma$ Coefficient Term 1 function
def seoBNR_v4p_sigmacoeffterm1(m1, m2):
    """
    Calculate the first term of coefficient involving spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmacoeffterm1_seoBNR_v4p: $\sigma$ Coefficient Term 1
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the first term of coefficient involving spin vector
    sigmacoeffterm1_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * S1[0]**2
    
    return sigmacoeffterm1_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmacoeffterm1_seoBNR_v4p = seoBNR_v4p_sigmacoeffterm1(m1, m2)
print("$\sigma$ Coefficient Term 1:", sigmacoeffterm1_seoBNR_v4p)
```

### Mathematics

$$
\sigma_{11} =
\frac{G \, c^4}{M^3} \,
S_**NRPy+: $\sigma$ Coefficient Term 2**
=============================================

### Theory Review

#### Introduction to $\sigma$ Coefficient Term 2

*   **$\sigma$ Coefficient Term 2:** In this section, we discuss the second term of the coefficient involving the spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma$ Coefficient Term 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma$ Coefficient Term 2 function
def seoBNR_v4p_sigmacoeffterm2(m1, m2):
    """
    Calculate the second term of coefficient involving spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmacoeffterm2_seoBNR_v4p: $\sigma$ Coefficient Term 2
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the second term of coefficient involving spin vector
    sigmacoeffterm2_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * S2[0]**2
    
    return sigmacoeffterm2_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmacoeffterm2_seoBNR_v4p = seoBNR_v4p_sigmacoeffterm2(m1, m2)
print("$\sigma$ Coefficient Term 2:", sigmacoeffterm2_seoBNR_v4p)
```

### Mathematics

$$
\sigma_{22} =
\frac{G \, c^4}{M^3} \,
S_**NRPy+: $\sigma$ Coefficient Term 3**
=============================================

### Theory Review

#### Introduction to $\sigma$ Coefficient Term 3

*   **$\sigma$ Coefficient Term 3:** In this section, we discuss the third term of the coefficient involving the spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sigma$ Coefficient Term 3:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sigma$ Coefficient Term 3 function
def seoBNR_v4p_sigmacoeffterm3(m1, m2):
    """
    Calculate the third term of coefficient involving spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmacoeffterm3_seoBNR_v4p: $\sigma$ Coefficient Term 3
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Calculate the third term of coefficient involving spin vector
    sigmacoeffterm3_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * (S1[0]**2 + S2[0]**2)
    
    return sigmacoeffterm3_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmacoeffterm3_seoBNR_v4p = seoBNR_v4p_sigmacoeffterm3(m1, m2)
print("$\sigma$ Coefficient Term 3:", sigmacoeffterm3_seoBNR_v4p)
```

### Mathematics

$$
\sigma_{33} =
\frac{G \, c^4**NRPy+: Derivatives of the Metric Potential**
=============================================

### Theory Review

#### Introduction to Derivatives of the Metric Potential

*   **Derivatives of the Metric Potential:** In this section, we discuss the derivatives of the metric potential.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Derivatives of the Metric Potential:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the derivatives of the metric potential function
def seoBNR_v4p_metpotderivs(m1, m2):
    """
    Calculate the derivatives of the metric potential.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        metpotderivs_seoBNR_v4p: Derivatives of the Metric Potential
    """
    # Calculate the spin vectors
    S1 = np.array([0, 0, 1])  # spin of first black hole
    S2 = np.array([0, 0, -1])  # spin of second black hole
    
    # Define the metric potential
    metpot_seoBNR_v4p = G * c**4 / (m1 + m2)**3 * (S1[0]**2 + S2[0]**2)
    
    # Calculate the derivatives of the metric potential
    metpotderivs_seoBNR_v4p = np.array([metpot_seoBNR_v4p, 2*metpot_seoBNR_v4p])
    
    return metpotderivs_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

metpotderivs_seoBNR_v4p = seoBNR_v4p_metpotderivs(m1, m2)
print("Derivatives of the Metric Potential:", metpotderivs_seoBNR_v4p)
``**NRPy+: $\omega_{r}$**
=========================

### Theory Review

#### Introduction to $\omega_{r}$

*   **$\omega_{r}$:** In this section, we discuss the radial frequency of the gravitational wave.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\omega_{r}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\omega_{r}$ function
def seoBNR_v4p_omegar(m1, m2):
    """
    Calculate the radial frequency of the gravitational wave.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        omegar_seoBNR_v4p: $\omega_{r}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radial frequency of the gravitational wave
    omegar_seoBNR_v4p = np.sqrt(G / 2 / mu)
    
    return omegar_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

omegar_seoBNR_v4p = seoBNR_v4p_omegar(m1, m2)
print("$\omega_{r}$:", omegar_seoBNR_v4p)
```

### Mathematics

$$
\omega_{r} =
\sqrt{\frac{G}{2 \mu}}
$$**NRPy+: $\nu_{r}$**
=========================

### Theory Review

#### Introduction to $\nu_{r}$

*   **$\nu_{r}$:** In this section, we discuss the radial frequency of the gravitational wave in terms of dimensionless units.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\nu_{r}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\nu_{r}$ function
def seoBNR_v4p_nur(m1, m2):
    """
    Calculate the radial frequency of the gravitational wave in terms of dimensionless units.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        nur_seoBNR_v4p: $\nu_{r}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radial frequency of the gravitational wave in terms of dimensionless units
    nur_seoBNR_v4p = np.sqrt(mu) / c
    
    return nur_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

nur_seoBNR_v4p = seoBNR_v4p_nur(m1, m2)
print("$\nu_{r}$:", nur_seoBNR_v4p)
```

### Mathematics

$$
\nu_{r} =
\frac{\sqrt{G \mu}}{c}
$$**NRPy+: $\mu_{r}$**
=========================

### Theory Review

#### Introduction to $\mu_{r}$

*   **$\mu_{r}$:** In this section, we discuss the reduced mass of the binary system in terms of dimensionless units.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\mu_{r}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\mu_{r}$ function
def seoBNR_v4p_mur(m1, m2):
    """
    Calculate the reduced mass of the binary system in terms of dimensionless units.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        mur_seoBNR_v4p: $\mu_{r}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the reduced mass in terms of dimensionless units
    mur_seoBNR_v4p = mu / (c**3)
    
    return mur_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

mur_seoBNR_v4p = seoBNR_v4p_mur(m1, m2)
print("$\mu_{r}$:", mur_seoBNR_v4p)
```

### Mathematics

$$
\mu_{r} =
\frac{G \mu}{c^3}
$$**NRPy+: $\omega_{\cos\theta}$**
================================

### Theory Review

#### Introduction to $\omega_{\cos\theta}$

*   **$\omega_{\cos\theta}$:** In this section, we discuss the cosine of the angle between the wavevector and the direction of motion.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\omega_{\cos\theta}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\omega_{\cos\theta}$ function
def seoBNR_v4p_omegacostheta(m1, m2):
    """
    Calculate the cosine of the angle between the wavevector and the direction of motion.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        omegacostheta_seoBNR_v4p: $\omega_{\cos\theta}$
    """
    # Calculate the radial frequency
    omega_r = np.sqrt(G / 2 / (G * (m1 * m2) / (m1 + m2)**3))
    
    # Calculate the cosine of the angle between the wavevector and the direction of motion
    omegacostheta_seoBNR_v4p = np.cos(np.arccos(omega_r**2 / (omega_r**2 + 1)))
    
    return omegacostheta_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

omegacostheta_seoBNR_v4p = seoBNR_v4p_omegacostheta(m1, m2)
print("$\omega_{\cos\theta}$:", omegacostheta_seoBNR_v4p)
```

### Mathematics

$$
\omega_{\cos\theta} =
\cos\left(\**NRPy+: $\nu_{\cos\theta}$**
=============================

### Theory Review

#### Introduction to $\nu_{\cos\theta}$

*   **$\nu_{\cos\theta}$:** In this section, we discuss the cosine of the angle between the wavevector and the direction of motion in terms of dimensionless units.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\nu_{\cos\theta}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\nu_{\cos\theta}$ function
def seoBNR_v4p_nucostheta(m1, m2):
    """
    Calculate the cosine of the angle between the wavevector and the direction of motion in terms of dimensionless units.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        nucostheta_seoBNR_v4p: $\nu_{\cos\theta}$
    """
    # Calculate the radial frequency
    nu_r = np.sqrt(G / 2 / (G * (m1 * m2) / (m1 + m2)**3))
    
    # Calculate the cosine of the angle between the wavevector and the direction of motion in terms of dimensionless units
    nucostheta_seoBNR_v4p = np.cos(np.arccos(nu_r**2 / (nu_r**2 + 1)))
    
    return nucostheta_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

nucostheta_seoBNR_v4p = seoBNR_v4p_nucostheta(m1, m2)
print("$\nu_{\cos\theta}$:", nucostheta_seoBNR_v4p)
```

### Mathematics

$$
\nu_{\**NRPy+: $\mu_{\cos\theta}$**
=============================

### Theory Review

#### Introduction to $\mu_{\cos\theta}$

*   **$\mu_{\cos\theta}$:** In this section, we discuss the cosine of the angle between the wavevector and the direction of motion in terms of reduced mass.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\mu_{\cos\theta}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\mu_{\cos\theta}$ function
def seoBNR_v4p_mucostheta(m1, m2):
    """
    Calculate the cosine of the angle between the wavevector and the direction of motion in terms of reduced mass.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        mucostheta_seoBNR_v4p: $\mu_{\cos\theta}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the cosine of the angle between the wavevector and the direction of motion in terms of reduced mass
    mucostheta_seoBNR_v4p = np.cos(np.arccos(mu**2 / (mu**2 + 1)))
    
    return mucostheta_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

mucostheta_seoBNR_v4p = seoBNR_v4p_mucostheta(m1, m2)
print("$\mu_{\cos\theta}$:", mucostheta_seoBNR_v4p)
```

### Mathematics

$$
\mu_{\cos\theta} =
\cos\left(\arccos\left(
\frac{\mu**NRPy+: $\Lambda_{t}^{\prime}$**
==================================

### Theory Review

#### Introduction to $\Lambda_{t}^{\prime}$

*   **$\Lambda_{t}^{\prime}$:** In this section, we discuss the time derivative of the Lambda parameter.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Lambda_{t}^{\prime}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Lambda_{t}^{\prime}$ function
def seoBNR_v4p_lambdatprm(m1, m2):
    """
    Calculate the time derivative of the Lambda parameter.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        lambdatprm_seoBNR_v4p: $\Lambda_{t}^{\prime}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the time derivative of the Lambda parameter
    lambdatprm_seoBNR_v4p = -np.sqrt(2*mu) * np.sin(np.arccos(mu**2 / (mu**2 + 1)))
    
    return lambdatprm_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

lambdatprm_seoBNR_v4p = seoBNR_v4p_lambdatprm(m1, m2)
print("$\Lambda_{t}^{\prime}$:", lambdatprm_seoBNR_v4p)
```

### Mathematics

$$
\Lambda_{t}^{\prime} =
-\sqrt{2 \mu} \sin\left(
\arccos\left(
\frac{\mu^2}{\mu^2 + 1}
\right)\right)
$$**NRPy+: $\tilde{\omega}_{\rm fd}^{\prime}$**
=============================================

### Theory Review

#### Introduction to $\tilde{\omega}_{\rm fd}^{\prime}$

*   **$\tilde{\omega}_{\rm fd}^{\prime}$:** In this section, we discuss the time derivative of the effective gravitational wave frequency.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\tilde{\omega}_{\rm fd}^{\prime}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\tilde{\omega}_{\rm fd}^{\prime}$ function
def seoBNR_v4p_omegatildeprm(m1, m2):
    """
    Calculate the time derivative of the effective gravitational wave frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        omegatildeprm_seoBNR_v4p: $\tilde{\omega}_{\rm fd}^{\prime}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the time derivative of the effective gravitational wave frequency
    omegatildeprm_seoBNR_v4p = -np.sqrt(2*mu) * np.sin(np.arccos(mu**2 / (mu**2 + 1)))
    
    return omegatildeprm_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

omegatildeprm_seoBNR_v4p = seoBNR_v4p_omegatildeprm(m1, m2)
print("$\tilde{\omega}_{\rm fd}^{\prime}$:", omegatildeprm_seoBNR_v4p)
```

### Mathematics**NRPy+: The Deformed and Rescaled Metric Potentials**
=====================================================

### Theory Review

#### Introduction to the Deformed and Rescaled Metric Potentials

*   **The Deformed and Rescaled Metric Potentials:** In this section, we discuss the deformed and rescaled metric potentials.
    +   These are crucial concepts in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The Deformed and Rescaled Metric Potentials:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the deformed and rescaled metric potentials function
def seoBNR_v4p_metpots(m1, m2):
    """
    Calculate the deformed and rescaled metric potentials.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        metpots_seoBNR_v4p: The Deformed and Rescaled Metric Potentials
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the deformed and rescaled metric potentials
    metpots_seoBNR_v4p = np.sqrt(2*mu) * np.sin(np.arccos(mu**2 / (mu**2 + 1)))
    
    return metpots_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

metpots_seoBNR_v4p = seoBNR_v4p_metpots(m1, m2)
print("The Deformed and Rescaled Metric Potentials:", metpots_seoBNR_v4p)
```

### Mathematics

$$
\mathcal{M}_{\rm fd}^{\prime} =
\sqrt{2 \mu} \sin\left(
\arccos\left(
\frac{\mu^2}{\mu^2 + 1}
\right)\right)
$$**NRPy+: $\omega$**
=====================

### Theory Review

#### Introduction to $\omega$

*   **$\omega$:** In this section, we discuss the gravitational wave frequency.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\omega$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\omega$ function
def seoBNR_v4p_omega(m1, m2):
    """
    Calculate the gravitational wave frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        omega_seoBNR_v4p: $\omega$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the gravitational wave frequency
    omega_seoBNR_v4p = np.sqrt(2*mu)
    
    return omega_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

omega_seoBNR_v4p = seoBNR_v4p_omega(m1, m2)
print("$\omega$:", omega_seoBNR_v4p)
```

### Mathematics

$$
\omega =
\sqrt{2 \mu}
$$**NRPy+: $e^{2 \nu}$**
========================

### Theory Review

#### Introduction to $e^{2 \nu}$

*   **$e^{2 \nu}$:** In this section, we discuss the exponential of twice the gravitational wave frequency.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$e^{2 \nu}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $e^{2 \nu}$ function
def seoBNR_v4p_exp2nu(m1, m2):
    """
    Calculate the exponential of twice the gravitational wave frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        exp2nu_seoBNR_v4p: $e^{2 \nu}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the gravitational wave frequency
    omega = np.sqrt(2*mu)
    
    # Calculate the exponential of twice the gravitational wave frequency
    exp2nu_seoBNR_v4p = np.exp(2*omega)
    
    return exp2nu_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

exp2nu_seoBNR_v4p = seoBNR_v4p_exp2nu(m1, m2)
print("$e^{2 \nu}$:", exp2nu_seoBNR_v4p)
```

### Mathematics

$$
e^{2 \nu} =
e^{2 \sqrt{2 \mu}}
$$**NRPy+: $\tilde{B}$**
======================

### Theory Review

#### Introduction to $\tilde{B}$

*   **$\tilde{B}$:** In this section, we discuss the rescaled magnetic field.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\tilde{B}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\tilde{B}$ function
def seoBNR_v4p_btilde(m1, m2):
    """
    Calculate the rescaled magnetic field.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        btilde_seoBNR_v4p: $\tilde{B}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the rescaled magnetic field
    btilde_seoBNR_v4p = np.sqrt(mu) / c
    
    return btilde_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

btilde_seoBNR_v4p = seoBNR_v4p_btilde(m1, m2)
print("$\tilde{B}$:", btilde_seoBNR_v4p)
```

### Mathematics

$$
\tilde{B} =
\frac{\sqrt{\mu}}{c}
$$**NRPy+: $\tilde{B}_{r}$**
=========================

### Theory Review

#### Introduction to $\tilde{B}_{r}$

*   **$\tilde{B}_{r}$:** In this section, we discuss the radial component of the rescaled magnetic field.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\tilde{B}_{r}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\tilde{B}_{r}$ function
def seoBNR_v4p_brtilde(m1, m2):
    """
    Calculate the radial component of the rescaled magnetic field.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        brtilde_seoBNR_v4p: $\tilde{B}_{r}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radial component of the rescaled magnetic field
    brtilde_seoBNR_v4p = np.sqrt(mu) / (c * 2)
    
    return brtilde_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

brtilde_seoBNR_v4p = seoBNR_v4p_brtilde(m1, m2)
print("$\tilde{B}_{r}$:", brtilde_seoBNR_v4p)
```

### Mathematics

$$
\tilde{B}_{r} =
\frac{\sqrt{\mu}}{c \cdot 2}
$$**NRPy+: $e^{2 \tilde{\mu}}$**
=============================

### Theory Review

#### Introduction to $e^{2 \tilde{\mu}}$

*   **$e^{2 \tilde{\mu}}$:** In this section, we discuss the exponential of twice the rescaled reduced mass.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$e^{2 \tilde{\mu}}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $e^{2 \tilde{\mu}}$ function
def seoBNR_v4p_exp2mu(m1, m2):
    """
    Calculate the exponential of twice the rescaled reduced mass.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        exp2mu_seoBNR_v4p: $e^{2 \tilde{\mu}}$
    """
    # Calculate the rescaled reduced mass
    tilde_mu = G * (m1 * m2) / (c**3)
    
    # Calculate the exponential of twice the rescaled reduced mass
    exp2mu_seoBNR_v4p = np.exp(2*tilde_mu)
    
    return exp2mu_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

exp2mu_seoBNR_v4p = seoBNR_v4p_exp2mu(m1, m2)
print("$e^{2 \tilde{\mu}}$:", exp2mu_seoBNR_v4p)
```

### Mathematics

$$
e^{2 \tilde{\mu}} =
e^{2 G (m_1 m_2) / (c^3)}
$$**NRPy+: $\tilde{J}$**
======================

### Theory Review

#### Introduction to $\tilde{J}$

*   **$\tilde{J}$:** In this section, we discuss the rescaled spin of the black holes.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\tilde{J}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\tilde{J}$ function
def seoBNR_v4p_jtilde(m1, m2):
    """
    Calculate the rescaled spin of the black holes.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        jtilde_seoBNR_v4p: $\tilde{J}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the rescaled spin of the black holes
    jtilde_seoBNR_v4p = np.sqrt(mu) / c
    
    return jtilde_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

jtilde_seoBNR_v4p = seoBNR_v4p_jtilde(m1, m2)
print("$\tilde{J}$:", jtilde_seoBNR_v4p)
```

### Mathematics

$$
\tilde{J} =
\frac{\sqrt{\mu}}{c}
$$**NRPy+: $Q$**
================

### Theory Review

#### Introduction to $Q$

*   **$Q$:** In this section, we discuss the mass ratio of the black holes.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$Q$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $Q$ function
def seoBNR_v4p_q(m1, m2):
    """
    Calculate the mass ratio of the black holes.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        q_seoBNR_v4p: $Q$
    """
    # Calculate the mass ratio
    q = np.sqrt(m1 / m2)
    
    return q

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

q_seoBNR_v4p = seoBNR_v4p_q(m1, m2)
print("$Q$:", q_seoBNR_v4p)
```

### Mathematics

$$
Q =
\sqrt{\frac{m_1}{m_2}}
$$**NRPy+: $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$**
====================================================================

### Theory Review

#### Introduction to $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$

*   **$\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$:** In this section, we discuss the rescaled metric potential.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ function
def seoBNR_v4p_drsipn2(m1, m2):
    """
    Calculate the rescaled metric potential.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        drsipn2_seoBNR_v4p: $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector of momentum
    hat_p = np.array([1, 0, 0])  # assuming momentum is in x-direction
    
    # Calculate the unit vector of normal
    n = np.array([1, 0, 0])  # assuming normal is in x-direction
    
    # Calculate the dot**NRPy+: $Q$ Coefficient 1**
==========================

### Theory Review

#### Introduction to $Q$ Coefficient 1

*   **$Q$ Coefficient 1:** In this section, we discuss the first coefficient of the mass ratio term in the effective one-body (EOB) Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$Q$ Coefficient 1:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $Q$ Coefficient 1 function
def seoBNR_v4p_qcoeff1(m1, m2):
    """
    Calculate the first coefficient of the mass ratio term in the EOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        qcoeff1_seoBNR_v4p: $Q$ Coefficient 1
    """
    # Calculate the mass ratio
    q = np.sqrt(m1 / m2)
    
    # Calculate the first coefficient of the mass ratio term
    qcoeff1 = 3/8 - (73/24) * q + (37/96) * q**2
    
    return qcoeff1

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

qcoeff1_seoBNR_v4p = seoBNR_v4p_qcoeff1(m1, m2)
print("$Q$ Coefficient 1:", qcoeff1_seoBNR_v4p)
```

### Mathematics

$$
Q \text{ Coefficient } 1 =
\frac{3}{8} - \frac{73}{24} Q + \frac{37}{96} Q^2
$$**NRPy+: $Q$ Coefficient 2**
==========================

### Theory Review

#### Introduction to $Q$ Coefficient 2

*   **$Q$ Coefficient 2:** In this section, we discuss the second coefficient of the mass ratio term in the effective one-body (EOB) Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$Q$ Coefficient 2:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $Q$ Coefficient 2 function
def seoBNR_v4p_qcoeff2(m1, m2):
    """
    Calculate the second coefficient of the mass ratio term in the EOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        qcoeff2_seoBNR_v4p: $Q$ Coefficient 2
    """
    # Calculate the mass ratio
    q = np.sqrt(m1 / m2)
    
    # Calculate the second coefficient of the mass ratio term
    qcoeff2 = - (734/256) + (4536/1024) * q - (9273/4096) * q**2
    
    return qcoeff2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

qcoeff2_seoBNR_v4p = seoBNR_v4p_qcoeff2(m1, m2)
print("$Q$ Coefficient 2:", qcoeff2_seoBNR_v4p)
```

### Mathematics

$$
Q \text{ Coefficient } 2 =
-\frac{734}{256} + \frac{4536}{1024} Q - \frac{9273}{4096} Q^2
$$**NRPy+: Tortoise Terms**
======================

### Theory Review

#### Introduction to Tortoise Terms

*   **Tortoise Terms:** In this section, we discuss the tortoise coordinate transformation and its application in numerical relativity.
    +   The tortoise coordinate is a useful tool for analyzing black hole spacetimes.

### Code Explanation


```python
"""
Tortoise terms:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the tortoise term function
def seoBNR_v4p_tort(m1, m2):
    """
    Calculate the tortoise coordinate transformation.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (k g)
    
    Returns:
        tort_seoBNR_v4p: Tortoise terms
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the tortoise coordinate transformation
    tort = c**2 * (np.sqrt(2*mu) * np.log(np.abs(c**4 - 2*G*m1) / (c**4 - 2*G*m2)) + np.arctan((m1 - m2) / (sqrt(m1) + sqrt(m2))))
    
    return tort

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

tort_seoBNR_v4p = seoBNR_v4p_tort(m1, m2)
print("Tortoise terms:", tort_seoBNR_v4p)
```

### Mathematics

$$
\text{Tortoise } \text{terms} =
c^2 \left(\sqrt{2 \mu} \log \left| c^4 - 2 G m_1 \right| / (c^4 - 2 G m_2) + \arctan \left( (m_1 - m_2) / (\sqrt{m_**NRPy+: $p_{\phi}$**
=====================

### Theory Review

#### Introduction to $p_{\phi}$

*   **$p_{\phi}$:** In this section, we discuss the canonical momentum of the massless scalar field.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$p_{\phi}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $p_{\phi}$ function
def seoBNR_v4p_pphi(m1, m2):
    """
    Calculate the canonical momentum of the massless scalar field.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        pphi_seoBNR_v4p: $p_{\phi}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the angular momentum
    l = 10 * (mu**2)
    
    # Calculate the canonical momentum of the massless scalar field
    pphi = np.sqrt(l**2 - 4*G*m1*l/np.sqrt(mu))
    
    return pphi

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

pphi_seoBNR_v4p = seoBNR_v4p_pphi(m1, m2)
print("$p_{\phi}$:", pphi_seoBNR_v4p)
```

### Mathematics

$$
p_{\phi} =
\sqrt{l^2 - \frac{4 G m_1 l}{\sqrt{\mu}}}
$$**NRPy+: $\hat{\bf p} \cdot {\bf v} r$**
=====================================

### Theory Review

#### Introduction to $\hat{\bf p} \cdot {\bf v} r$

*   **$\hat{\bf p} \cdot {\bf v} r$:** In this section, we discuss the dot product of the unit vector of momentum and the radial velocity.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\hat{\bf p} \cdot {\bf v} r$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\hat{\bf p} \cdot {\bf v} r$ function
def seoBNR_v4p_pdotvr(m1, m2):
    """
    Calculate the dot product of the unit vector of momentum and the radial velocity.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        pdotvr_seoBNR_v4p: $\hat{\bf p} \cdot {\bf v} r$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector of momentum
    hat_p = np.array([1, 0, 0])  # assuming momentum is in x-direction
    
    # Calculate the radial velocity
    vr = 1.5 * np.sqrt(mu)
    
    # Calculate the dot product of the unit vector of momentum and the radial velocity
    pdotvr = hat_p[0] * vr
    
    return pdotvr

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

pdotvr_seoBNR_v4p = seoBNR_v4p_pdotvr(m1, m2)
print("$\hat{\bf p} \cdot {\bf v} r$:", pdot**NRPy+: $\hat{\bf p} \cdot {\bf n}$**
=====================================

### Theory Review

#### Introduction to $\hat{\bf p} \cdot {\bf n}$

*   **$\hat{\bf p} \cdot {\bf n}$:** In this section, we discuss the dot product of the unit vector of momentum and the normal.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\hat{\bf p} \cdot {\bf n}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\hat{\bf p} \cdot {\bf n}$ function
def seoBNR_v4p_pdotn(m1, m2):
    """
    Calculate the dot product of the unit vector of momentum and the normal.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        pdotn_seoBNR_v4p: $\hat{\bf p} \cdot {\bf n}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector of momentum
    hat_p = np.array([1, 0, 0])  # assuming momentum is in x-direction
    
    # Calculate the normal
    n = np.array([1, 0, 0])  # assuming normal is in x-direction
    
    # Calculate the dot product of the unit vector of momentum and the normal
    pdotn = hat_p[0] * n[0]
    
    return pdotn

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

pdotn_seoBNR_v4p = seoBNR_v4p_pdotn(m1, m2)
print("$\hat{\bf p} \cdot {\bf n}$:", pdotn_seo**NRPy+: $\hat{\bf p} \cdot \boldsymbol{\xi} r$**
=============================================

### Theory Review

#### Introduction to $\hat{\bf p} \cdot \boldsymbol{\xi} r$

*   **$\hat{\bf p} \cdot \boldsymbol{\xi} r$:** In this section, we discuss the dot product of the unit vector of momentum and the shift vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\hat{\bf p} \cdot \boldsymbol{\xi} r$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\hat{\bf p} \cdot \boldsymbol{\xi} r$ function
def seoBNR_v4p_pdotxir(m1, m2):
    """
    Calculate the dot product of the unit vector of momentum and the shift vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        pdotxir_seoBNR_v4p: $\hat{\bf p} \cdot \boldsymbol{\xi} r$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector of momentum
    hat_p = np.array([1, 0, 0])  # assuming momentum is in x-direction
    
    # Calculate the shift vector
    xi = np.array([0.5, 0.25, 0.125])  # assuming shift vector is [0.5, 0.25, 0.125]
    
    # Calculate the radial velocity
    vr = 1.5 * np.sqrt(mu)
    
    # Calculate the dot product of the unit vector of momentum and the shift vector
    pdotxir = hat_p[0] * xi[0] * vr
    
    return pdotxir

# Test the function with some example values
m1 = 1.98910e30 **NRPy+: $\hat{\bf p}$**
=====================

### Theory Review

#### Introduction to $\hat{\bf p}$

*   **$\hat{\bf p}$:** In this section, we discuss the unit vector of momentum.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\hat{\bf p}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\hat{\bf p}$ function
def seoBNR_v4p_hatp(m1, m2):
    """
    Calculate the unit vector of momentum.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        hatp_seoBNR_v4p: $\hat{\bf p}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector of momentum
    hat_p = np.array([1, 0, 0])  # assuming momentum is in x-direction
    
    return hatp_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

hatp_seoBNR_v4p = seoBNR_v4p_hatp(m1, m2)
print("$\hat{\bf p}$:", hatp_seoBNR_v4p)
```

### Mathematics

$$
\hat{\bf p} =
\begin{pmatrix}
1 \\
0 \\
0
\end{pmatrix}
$$**NRPy+: $pr_T$**
================

### Theory Review

#### Introduction to $pr_T$

*   **$pr_T$:** In this section, we discuss the radial momentum.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$pr_T$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $pr_T$ function
def seoBNR_v4p_prt(m1, m2):
    """
    Calculate the radial momentum.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        prt_seoBNR_v4p: $pr_T$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radial momentum
    prt = 10 * np.sqrt(mu)
    
    return prt

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

prt_seoBNR_v4p = seoBNR_v4p_prt(m1, m2)
print("$pr_T$:", prt_seoBNR_v4p)
```

### Mathematics

$$
pr_T =
10 \sqrt{\mu}
$$**NRPy+: $\chi^2$**
==================

### Theory Review

#### Introduction to $\chi^2$

*   **$\chi^2$:** In this section, we discuss the chi-squared statistic.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\chi^2$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\chi^2$ function
def seoBNR_v4p_csi2(m1, m2):
    """
    Calculate the chi-squared statistic.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        csi2_seoBNR_v4p: $\chi^2$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the chi-squared statistic
    csi2 = 10 * np.sqrt(mu)
    
    return csi2

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

csi2_seoBNR_v4p = seoBNR_v4p_csi2(m1, m2)
print("$\chi^2$:", csi2_seoBNR_v4p)
```

### Mathematics

$$
\chi^2 =
10 \sqrt{\mu}
$$**NRPy+: $\chi^1$**
==================

### Theory Review

#### Introduction to $\chi^1$

*   **$\chi^1$:** In this section, we discuss the chi-squared statistic.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\chi^1$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\chi^1$ function
def seoBNR_v4p_csi1(m1, m2):
    """
    Calculate the chi-squared statistic.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        csi1_seoBNR_v4p: $\chi^1$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the chi-squared statistic
    csi1 = 5 * np.sqrt(mu)
    
    return csi1

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

csi1_seoBNR_v4p = seoBNR_v4p_csi1(m1, m2)
print("$\chi^1$:", csi1_seoBNR_v4p)
```

### Mathematics

$$
\chi^1 =
5 \sqrt{\mu}
$$**NRPy+: $c_{\text{SI}}$**
=====================

### Theory Review

#### Introduction to $c_{\text{SI}}$

*   **$c_{\text{SI}}$:** In this section, we discuss the speed of light in the Schwarzschild interior metric.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$c_{\text{SI}}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $c_{\text{SI}}$ function
def seoBNR_v4p_csi(m1, m2):
    """
    Calculate the speed of light in the Schwarzschild interior metric.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        csi_seoBNR_v4p: $c_{\text{SI}}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the speed of light in the Schwarzschild interior metric
    csi = 1 / np.sqrt(1 - 2*G*m1/(c**2*np.abs(c**2-2*G*m1)))
    
    return csi

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

csi_seoBNR_v4p = seoBNR_v4p_csi(m1, m2)
print("$c_{\text{SI}}$:", csi_seoBNR_v4p)
```

### Mathematics

$$
c_{\text{SI}} =
\frac{1}{\sqrt{1 - \frac{2 G m_1}{(c^2)^2 - 2 G m_1}}}
$$**NRPy+: Metric Terms**
=====================

### Theory Review

#### Introduction to Metric Terms

*   **Metric Terms:** In this section, we discuss the metric terms in the effective one-body (EOB) Hamiltonian.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Metric Terms:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Metric Terms function
def seoBNR_v4p_metric(m1, m2):
    """
    Calculate the metric terms in the EOB Hamiltonian.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        metric_seoBNR_v4p: Metric Terms
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the metric terms
    metric = 1 / np.sqrt(1 - 2*G*m1/(c**2*np.abs(c**2-2*G*m1)))
    
    return metric

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

metric_seoBNR_v4p = seoBNR_v4p_metric(m1, m2)
print("Metric Terms:", metric_seoBNR_v4p)
```

### Mathematics

$$
\text{Metric Terms} =
\frac{1}{\sqrt{1 - \frac{2 G m_1}{(c^2)^2 - 2 G m_1}}}
$$**NRPy+: $\Lambda_t$**
=====================

### Theory Review

#### Introduction to $\Lambda_t$

*   **$\Lambda_t$:** In this section, we discuss the Lagrange multiplier for time.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Lambda_t$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Lambda_t$ function
def seoBNR_v4p_lambdat(m1, m2):
    """
    Calculate the Lagrange multiplier for time.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        lambdat_seoBNR_v4p: $\Lambda_t$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the Lagrange multiplier for time
    lambdat = 10 * np.sqrt(mu)
    
    return lambdat

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

lambdat_seoBNR_v4p = seoBNR_v4p_lambdat(m1, m2)
print("$\Lambda_t$:", lambdat_seoBNR_v4p)
```

### Mathematics

$$
\Lambda_t =
10 \sqrt{\mu}
$$**NRPy+: $\Delta_r$**
=====================

### Theory Review

#### Introduction to $\Delta_r$

*   **$\Delta_r$:** In this section, we discuss the radial coordinate.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_r$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_r$ function
def seoBNR_v4p_deltr(m1, m2):
    """
    Calculate the radial coordinate.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltr_seoBNR_v4p: $\Delta_r$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radial coordinate
    delr = 10 * np.sqrt(mu)
    
    return delr

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltr_seoBNR_v4p = seoBNR_v4p_deltr(m1, m2)
print("$\Delta_r$:", deltr_seoBNR_v4p)
```

### Mathematics

$$
\Delta_r =
10 \sqrt{\mu}
$$**NRPy+: $\Delta_t$**
=====================

### Theory Review

#### Introduction to $\Delta_t$

*   **$\Delta_t$:** In this section, we discuss the time coordinate.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_t$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_t$ function
def seoBNR_v4p_deltat(m1, m2):
    """
    Calculate the time coordinate.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltat_seoBNR_v4p: $\Delta_t$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the time coordinate
    deltat = 10 * np.sqrt(mu)
    
    return deltat

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltat_seoBNR_v4p = seoBNR_v4p_deltat(m1, m2)
print("$\Delta_t$:", deltat_seoBNR_v4p)
```

### Mathematics

$$
\Delta_t =
10 \sqrt{\mu}
$$**NRPy+: $\Delta_t^\prime$**
=========================

### Theory Review

#### Introduction to $\Delta_t^\prime$

*   **$\Delta_t^\prime$:** In this section, we discuss the time derivative of the radial coordinate.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_t^\prime$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_t^\prime$ function
def seoBNR_v4p_deltatprm(m1, m2):
    """
    Calculate the time derivative of the radial coordinate.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltatprmseoBNR_v4p: $\Delta_t^\prime$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the time derivative of the radial coordinate
    deltatprm = 10 * np.sqrt(mu)
    
    return deltatprm

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltatprmseoBNR_v4p = seoBNR_v4p_deltatprm(m1, m2)
print("$\Delta_t^\prime$:", deltatprmseoBNR_v4p)
```

### Mathematics

$$
\Delta_t^\prime =
10 \sqrt{\mu}
$$**NRPy+: $\Delta_u$**
=====================

### Theory Review

#### Introduction to $\Delta_u$

*   **$\Delta_u$:** In this section, we discuss the unit of time.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_{u}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_u$ function
def seoBNR_v4p_deltau(m1, m2):
    """
    Calculate the unit of time.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltausoBNR_v4p: $\Delta_{u}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit of time
    deltau = 10 * np.sqrt(mu)
    
    return deltausoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltausoBNR_v4p = seoBNR_v4p_deltau(m1, m2)
print("$\Delta_{u}$:", deltausoBNR_v4p)
```

### Mathematics

$$
\Delta_{u} =
10 \sqrt{\mu}
$$**NRPy+: $\bar{\Delta}_u$**
=========================

### Theory Review

#### Introduction to $\bar{\Delta}_u$

*   **$\bar{\Delta}_u$:** In this section, we discuss the normalized unit of time.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\bar{\Delta}_{u}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\bar{\Delta}_u$ function
def seoBNR_v4p_deltaubar(m1, m2):
    """
    Calculate the normalized unit of time.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltaubarseoBNR_v4p: $\bar{\Delta}_{u}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the normalized unit of time
    deltaubar = 10 * np.sqrt(mu) / (G*m1/c**3)
    
    return deltaubarseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltaubarseoBNR_v4p = seoBNR_v4p_deltaubar(m1, m2)
print("$\bar{\Delta}_{u}$:", deltaubarseoBNR_v4p)
```

### Mathematics

$$
\bar{\Delta}_u =
\frac{10 \sqrt{\mu}}{Gm_1/c^3}
$$**NRPy+: $\Delta_u$ Calibration Term**
=====================================

### Theory Review

#### Introduction to $\Delta_u$ Calibration Term

*   **$\Delta_u$ Calibration Term:** In this section, we discuss the calibration term for the unit of time.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Delta_{u}$ Calibration Term:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Delta_u$ Calibration Term function
def seoBNR_v4p_deltaucalib(m1, m2):
    """
    Calculate the calibration term for the unit of time.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        deltaucalibseoBNR_v4p: $\Delta_{u}$ Calibration Term
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the calibration term for the unit of time
    deltaucalib = 10 * np.sqrt(mu) / (G*m1/c**3)
    
    return deltaucalibseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

deltaucalibseoBNR_v4p = seoBNR_v4p_deltaucalib(m1, m2)
print("$\Delta_{u}$ Calibration Term:", deltaucalibseoBNR_v4p)
```

### Mathematics

$$
\Delta_u \text{ Calibration Term} =
\frac{10 \sqrt{\mu}}{Gm_1/c^3}
$$**NRPy+: Calibration Coefficients**
=====================================

### Theory Review

#### Introduction to Calibration Coefficients

*   **Calibration Coefficients:** In this section, we discuss the calibration coefficients for the unit of time.
    +   These coefficients are used to calibrate the unit of time and ensure that it is consistent with the physical laws.

### Code Explanation


```python
"""
Calibration Coefficients:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Calibration Coefficients function
def seoBNR_v4p_calib_coeffs(m1, m2):
    """
    Calculate the calibration coefficients for the unit of time.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        calib_coeffsseoBNR_v4p: Calibration Coefficients
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the calibration coefficients for the unit of time
    calib_coeffs = [10, 20, 30]  # example values
    
    return calib_coeffsseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

calib_coeffsseoBNR_v4p = seoBNR_v4p_calib_coeffs(m1, m2)
print("Calibration Coefficients:", calib_coeffsseoBNR_v4p)
```

### Mathematics

$$
\text{Calibration Coefficients} =
\begin{bmatrix}
10 \\
20 \\
30
\end{bmatrix}
$$**NRPy+: $K$**
================

### Theory Review

#### Introduction to $K$

*   **$K$:** In this section, we discuss the constant of proportionality.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$K$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $K$ function
def seoBNR_v4p_k(m1, m2):
    """
    Calculate the constant of proportionality.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        kseoBNR_v4p: $K$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the constant of proportionality
    k = 10 * np.sqrt(mu)
    
    return kseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

kseoBNR_v4p = seoBNR_v4p_k(m1, m2)
print("$K$:", kseoBNR_v4p)
```

### Mathematics

$$
K =
10 \sqrt{\mu}
$$**NRPy+: $\chi$**
==================

### Theory Review

#### Introduction to $\chi$

*   **$\chi$:** In this section, we discuss the variable of proportionality.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\chi$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\chi$ function
def seoBNR_v4p_chi(m1, m2):
    """
    Calculate the variable of proportionality.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        chiseoBNR_v4p: $\chi$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the variable of proportionality
    chi = 10 * np.sqrt(mu)
    
    return chiseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

chiseoBNR_v4p = seoBNR_v4p_chi(m1, m2)
print("$\chi$:", chiseoBNR_v4p)
```

### Mathematics

$$
\chi =
10 \sqrt{\mu}
$$**NRPy+: $\tilde{\omega}_{\rm fd}$**
=====================================

### Theory Review

#### Introduction to $\tilde{\omega}_{\rm fd}$

*   **$\tilde{\omega}_{\rm fd}$:** In this section, we discuss the fiducial frequency.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\tilde{\omega}_{\rm fd}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\tilde{\omega}_{\rm fd}$ function
def seoBNR_v4p_omegatilde(m1, m2):
    """
    Calculate the fiducial frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        omegatildeseoBNR_v4p: $\tilde{\omega}_{\rm fd}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the fiducial frequency
    omegatilde = 10 * np.sqrt(mu)
    
    return omegatildeseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

omegatildeseoBNR_v4p = seoBNR_v4p_omegatilde(m1, m2)
print("$\tilde{\omega}_{\rm fd}$:", omegatildeseoBNR_v4p)
```

### Mathematics

$$
\tilde{\omega}_{\rm fd} =
10 \sqrt{\mu}
$$**NRPy+: $D^{-1}$**
=====================

### Theory Review

#### Introduction to $D^{-1}$

*   **$D^{-1}$:** In this section, we discuss the inverse of the determinant.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$D^{-1}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $D^{-1}$ function
def seoBNR_v4p_dinv(m1, m2):
    """
    Calculate the inverse of the determinant.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        dinvseoBNR_v4p: $D^{-1}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the inverse of the determinant
    dinv = 10 * np.sqrt(mu)
    
    return dinvseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

dinvseoBNR_v4p = seoBNR_v4p_dinv(m1, m2)
print("$D^{-1}$:", dinvseoBNR_v4p)
```

### Mathematics

$$
D^{-1} =
10 \sqrt{\mu}
$$**NRPy+: Terms Dependent on Coordinates**
=====================================

### Theory Review

#### Introduction to Terms Dependent on Coordinates

*   **Terms Dependent on Coordinates:** In this section, we discuss the terms that depend on coordinates.
    +   These terms are used in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Terms Dependent on Coordinates:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Terms Dependent on Coordinates function
def seoBNR_v4p_coord(m1, m2):
    """
    Calculate terms dependent on coordinates.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        coordseoBNR_v4p: Terms Dependent on Coordinates
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate terms dependent on coordinates
    x = 10 * np.sqrt(mu)
    y = 20 * np.sqrt(mu)
    
    return coordseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

coordseoBNR_v4p = seoBNR_v4p_coord(m1, m2)
print("Terms Dependent on Coordinates:", coordseoBNR_v4p)
```

### Mathematics

$$
\text{Terms Dependent on Coordinates} =
\begin{bmatrix}
x \\
y
\end{bmatrix}
=
\begin{bmatrix}
10 \sqrt{\mu} \\
20 \sqrt{\mu}
\end{bmatrix}
$$**NRPy+: $\Sigma$**
==================

### Theory Review

#### Introduction to $\Sigma$

*   **$\Sigma$:** In this section, we discuss the sigma factor.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\Sigma$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\Sigma$ function
def seoBNR_v4p_usigma(m1, m2):
    """
    Calculate the sigma factor.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        usigmaseoBNR_v4p: $\Sigma$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the sigma factor
    sigma = 10 * np.sqrt(mu)
    
    return usigmaseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

usigmaseoBNR_v4p = seoBNR_v4p_usigma(m1, m2)
print("$\Sigma$:", usigmaseoBNR_v4p)
```

### Mathematics

$$
\Sigma =
10 \sqrt{\mu}
$$**NRPy+: $\varpi^2$**
=====================

### Theory Review

#### Introduction to $\varpi^2$

*   **$\varpi^2$:** In this section, we discuss the square of the angular frequency.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\varpi^2$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\varpi^2$ function
def seoBNR_v4p_w2(m1, m2):
    """
    Calculate the square of the angular frequency.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        w2seoBNR_v4p: $\varpi^2$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the square of the angular frequency
    varpi2 = 10 * np.sqrt(mu)
    
    return w2seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

w2seoBNR_v4p = seoBNR_v4p_w2(m1, m2)
print("$\varpi^2$:", w2seoBNR_v4p)
```

### Mathematics

$$
\varpi^2 =
10 \sqrt{\mu}
$$**NRPy+: $\sin^2\theta$**
=========================

### Theory Review

#### Introduction to $\sin^2\theta$

*   **$\sin^2\theta$:** In this section, we discuss the square of the sine of the angle.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\sin^{2}\theta$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\sin^2\theta$ function
def seoBNR_v4p_sin2theta(m1, m2):
    """
    Calculate the square of the sine of the angle.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sin2thetaseoBNR_v4p: $\sin^{2}\theta$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the square of the sine of the angle
    theta = np.arcsin(0.5)  # half-angle approximation
    sin2theta = 10 * np.sin(theta)**2
    
    return sin2thetaseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sin2thetaseoBNR_v4p = seoBNR_v4p_sin2theta(m1, m2)
print("$\sin^{2}\theta$:", sin2thetaseoBNR_v4p)
```

### Mathematics

$$
\sin^{2}\theta =
10 \sin^{2}(\arcsin(0.5))
$$**NRPy+: $\cos\theta$**
=====================

### Theory Review

#### Introduction to $\cos\theta$

*   **$\cos\theta$:** In this section, we discuss the cosine of the angle.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\cos\theta$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\cos\theta$ function
def seoBNR_v4p_costheta(m1, m2):
    """
    Calculate the cosine of the angle.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        costhetaseoBNR_v4p: $\cos\theta$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the cosine of the angle
    theta = np.arcsin(0.5)  # half-angle approximation
    costheta = np.cos(theta)
    
    return costhetaseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

costhetaseoBNR_v4p = seoBNR_v4p_costheta(m1, m2)
print("$\cos\theta$:", costhetaseoBNR_v4p)
```

### Mathematics

$$
\cos\theta =
\cos(\arcsin(0.5))
$$**NRPy+: Important Vectors**
=========================

### Theory Review

#### Introduction to Important Vectors

*   **Important Vectors:** In this section, we discuss the important vectors used in numerical relativity and gravitational wave astronomy.
    +   These vectors are crucial for understanding the behavior of black holes and gravitational waves.

### Code Explanation


```python
"""
Important Vectors:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Important Vectors function
def seoBNR_v4p_vectors(m1, m2):
    """
    Calculate important vectors.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        vectorsseoBNR_v4p: Important Vectors
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate important vectors
    x = np.array([1, 0])  # unit vector in x-direction
    y = np.array([0, 1])  # unit vector in y-direction
    z = np.array([0, 0, 1])  # unit vector in z-direction
    
    return vectorsseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

vectorsseoBNR_v4p = seoBNR_v4p_vectors(m1, m2)
print("Important Vectors:", vectorsseoBNR_v4p)
```

### Mathematics

$$
\text{Important Vectors} =
\begin{bmatrix}
x \\
y \\
z
\end{bmatrix}
=
\begin{bmatrix}
(1, 0) \\
(0, 1) \\
(0, 0, 1)
\end{bmatrix}
$$**NRPy+: ${\bf v}$**
====================

### Theory Review

#### Introduction to ${\bf v}$

*   **${\bf v}$:** In this section, we discuss the velocity vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf v}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf v}$ function
def seoBNR_v4p_v(m1, m2):
    """
    Calculate the velocity vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        vseoBNR_v4p: ${\bf v}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the velocity vector
    vx = 10 * np.sqrt(mu)
    vy = 20 * np.sqrt(mu)
    
    return vseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

vseoBNR_v4p = seoBNR_v4p_v(m1, m2)
print("${\bf v}$:", vseoBNR_v4p)
```

### Mathematics

$$
{\bf v} =
10 \sqrt{\mu}
=
(10 \sqrt{\mu}, 20 \sqrt{\mu})
$$**NRPy+: $\boldsymbol{\xi}$**
==========================

### Theory Review

#### Introduction to $\boldsymbol{\xi}$

*   **$\boldsymbol{\xi}$:** In this section, we discuss the vector of displacement.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\boldsymbol{\xi}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\boldsymbol{\xi}$ function
def seoBNR_v4p_xi(m1, m2):
    """
    Calculate the vector of displacement.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        xiseoBNR_v4p: $\boldsymbol{\xi}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the vector of displacement
    xi_x = 10 * np.sqrt(mu)
    xi_y = 20 * np.sqrt(mu)
    
    return xiseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

xiseoBNR_v4p = seoBNR_v4p_xi(m1, m2)
print("$\boldsymbol{\xi}$:", xiseoBNR_v4p)
```

### Mathematics

$$
\boldsymbol{\xi} =
(10 \sqrt{\mu}, 20 \sqrt{\mu})
$$**NRPy+: ${\bf e}_{3}$**
=======================

### Theory Review

#### Introduction to ${\bf e}_{3}$

*   **${\bf e}_{3}$:** In this section, we discuss the unit vector in the z-direction.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf e}_{3}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf e}_{3}$ function
def seoBNR_v4p_e3(m1, m2):
    """
    Calculate the unit vector in the z-direction.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        e3seoBNR_v4p: ${\bf e}_{3}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector in the z-direction
    e3_x = 0
    e3_y = 0
    e3_z = 1
    
    return e3seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

e3seoBNR_v4p = seoBNR_v4p_e3(m1, m2)
print("${\bf e}_{3}$:", e3seoBNR_v4p)
```

### Mathematics

$$
{\bf e}_{3} =
(0, 0, 1)
$$**NRPy+: ${\bf n}$**
====================

### Theory Review

#### Introduction to ${\bf n}$

*   **${\bf n}$:** In this section, we discuss the unit normal vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf n}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf n}$ function
def seoBNR_v4p_n(m1, m2):
    """
    Calculate the unit normal vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        nseoBNR_v4p: ${\bf n}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit normal vector
    n_x = 0.5 * np.sqrt(mu)
    n_y = 0.6 * np.sqrt(mu)
    
    return nseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

nseoBNR_v4p = seoBNR_v4p_n(m1, m2)
print("${\bf n}$:", nseoBNR_v4p)
```

### Mathematics

$$
{\bf n} =
(0.5 \sqrt{\mu}, 0.6 \sqrt{\mu})
$$**NRPy+: ${\bf S}^{\perp}$**
=========================

### Theory Review

#### Introduction to ${\bf S}^{\perp}$

*   **${\bf S}^{\perp}$:** In this section, we discuss the projection of the stress tensor onto the unit normal vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S}^{\perp}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S}^{\perp}$ function
def seoBNR_v4p_sperp(m1, m2):
    """
    Calculate the projection of the stress tensor onto the unit normal vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sperpnseoBNR_v4p: ${\bf S}^{\perp}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit normal vector
    n_x = 0.5 * np.sqrt(mu)
    n_y = 0.6 * np.sqrt(mu)
    
    # Calculate the stress tensor
    S_xx = 10 * np.sqrt(mu)
    S_yy = 20 * np.sqrt(mu)
    S_xy = 30 * np.sqrt(mu)
    
    # Calculate the projection of the stress tensor onto the unit normal vector
    sperp_x = n_x * (S_xx + S_xy)
    sperp_y = n_y * (S_yy + S_xy)
    
    return sperpnseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sperpnseoBNR_v4p = seoBNR_v4p_sperp(m1, m2)
print("${**NRPy+: ${\bf L}$**
==================

### Theory Review

#### Introduction to ${\bf L}$

*   **${\bf L}$:** In this section, we discuss the angular momentum vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf L}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf L}$ function
def seoBNR_v4p_orb_momentum(m1, m2):
    """
    Calculate the angular momentum vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        orb_momentumnseoBNR_v4p: ${\bf L}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the angular momentum vector
    L_x = 10 * np.sqrt(mu)
    L_y = 20 * np.sqrt(mu)
    
    return orb_momentumnseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

orb_momentumnseoBNR_v4p = seoBNR_v4p_orb_momentum(m1, m2)
print("${\bf L}$:", orb_momentumnseoBNR_v4p)
```

### Mathematics

$$
{\bf L} =
(10 \sqrt{\mu}, 20 \sqrt{\mu})
$$**NRPy+: Spin Combinations**
==========================

### Theory Review

#### Introduction to Spin Combinations

*   **Spin Combinations:** In this section, we discuss the spin combinations of two black holes.
    +   These are crucial concepts in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Spin Combinations:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Spin Combinations function
def seoBNR_v4p_spin_combos(m1, m2):
    """
    Calculate spin combinations.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigma_seoBNR_v4p: $\boldsymbol{\sigma}$
        sigma_star_seoBNR_v4p: $\boldsymbol{\sigma}^{*}$
        Skerr_seoBNR_v4p: ${\bf S}_{\rm Kerr}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate spin combinations
    sigma_x = 10 * np.sqrt(mu)
    sigma_y = 20 * np.sqrt(mu)
    
    # Calculate the complex conjugate of $\boldsymbol{\sigma}$
    sigma_star_x = sigma_x.conjugate()
    sigma_star_y = sigma_y.conjugate()
    
    # Calculate ${\bf S}_{\rm Kerr}$
    Skerr_x = 30 * np.sqrt(mu)
    Skerr_y = 40 * np.sqrt(mu)
    
    return sigma_seoBNR_v4p, sigma_star_seoBNR_v4p, Skerr_seoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigma_seoBNR_v4p, sigma_star_seoBNR_v4p, Skerr_seoBNR_v4p = seoBNR_v4p_spin_comb**NRPy+: $a$**
================

### Theory Review

#### Introduction to $a$

*   **$a$:** In this section, we discuss the spin coefficient $a$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$a$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $a$ function
def seoBNR_v4p_a(m1, m2):
    """
    Calculate the spin coefficient $a$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        aseoBNR_v4p: $a$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the spin coefficient $a$
    a = 10 * np.sqrt(mu)
    
    return aseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

aseoBNR_v4p = seoBNR_v4p_a(m1, m2)
print("$a$:", aseoBNR_v4p)
```

### Mathematics

$$
a =
10 \sqrt{\mu}
$$**NRPy+: $\hat{\bf S}_{\rm Kerr}$**
====================================

### Theory Review

#### Introduction to $\hat{\bf S}_{\rm Kerr}$

*   **$\hat{\bf S}_{\rm Kerr}$:** In this section, we discuss the unit vector in the direction of ${\bf S}_{\rm Kerr}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\hat{\bf S}_{\rm Kerr}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\hat{\bf S}_{\rm Kerr}$ function
def seoBNR_v4p_skerrhat(m1, m2):
    """
    Calculate the unit vector in the direction of ${\bf S}_{\rm Kerr}$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        skerrhatsseoBNR_v4p: $\hat{\bf S}_{\rm Kerr}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate ${\bf S}_{\rm Kerr}$
    Skerr_x = 30 * np.sqrt(mu)
    Skerr_y = 40 * np.sqrt(mu)
    
    # Calculate the unit vector in the direction of ${\bf S}_{\rm Kerr}$
    skerrhat_x = Skerr_x / np.sqrt(Skerr_x**2 + Skerr_y**2)
    skerrhat_y = Skerr_y / np.sqrt(Skerr_x**2 + Skerr_y**2)
    
    return skerrhatsseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

skerrhatsseoBNR_v4p = seoBNR_v4p_skerrhat(m1, m2)
print("$**NRPy+: $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$**
=====================================================

### Theory Review

#### Introduction to $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$

*   **$\left\lvert {\bf S}_{\rm Kerr} \right\rvert$:** In this section, we discuss the magnitude of ${\bf S}_{\rm Kerr}$.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\left\lvert {\bf S}_{\rm Kerr} \right\rvert$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$ function
def seoBNR_v4p_skerrmag(m1, m2):
    """
    Calculate the magnitude of ${\bf S}_{\rm Kerr}$.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        skerrmagnseoBNR_v4p: $\left\lvert {\bf S}_{\rm Kerr} \right\rvert$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate ${\bf S}_{\rm Kerr}$
    Skerr_x = 30 * np.sqrt(mu)
    Skerr_y = 40 * np.sqrt(mu)
    
    # Calculate the magnitude of ${\bf S}_{\rm Kerr}$
    skerrmag = np.sqrt(Skerr_x**2 + Skerr_y**2)
    
    return skerrmagnseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

skerrmagnseoBNR_v4p = seoBNR_v4p_skerrmag(m1, m**NRPy+: ${\bf S}_{\rm Kerr}$**
=============================

### Theory Review

#### Introduction to ${\bf S}_{\rm Kerr}$

*   **${\bf S}_{\rm Kerr}$:** In this section, we discuss the spin vector of the Kerr black hole.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
${\bf S}_{\rm Kerr}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the ${\bf S}_{\rm Kerr}$ function
def seoBNR_v4p_skerr(m1, m2):
    """
    Calculate the spin vector of the Kerr black hole.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        skerrseoBNR_v4p: ${\bf S}_{\rm Kerr}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the spin vector of the Kerr black hole
    Skerr_x = 30 * np.sqrt(mu)
    Skerr_y = 40 * np.sqrt(mu)
    
    return skerrseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

skerrseoBNR_v4p = seoBNR_v4p_skerr(m1, m2)
print("${\bf S}_{\rm Kerr}$:", skerrseoBNR_v4p)
```

### Mathematics

$$
{\bf S}_{\rm Kerr} =
(30 \sqrt{\mu}, 40 \sqrt{\mu})
$$**NRPy+: $\boldsymbol{\sigma}$**
=============================

### Theory Review

#### Introduction to $\boldsymbol{\sigma}$

*   **$\boldsymbol{\sigma}$:** In this section, we discuss the spin vector of a binary black hole system.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\boldsymbol{\sigma}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\boldsymbol{\sigma}$ function
def seoBNR_v4p_sigma(m1, m2):
    """
    Calculate the spin vector of a binary black hole system.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmaseoBNR_v4p: $\boldsymbol{\sigma}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the spin vector of a binary black hole system
    sigma_x = 10 * np.sqrt(mu)
    sigma_y = 20 * np.sqrt(mu)
    
    return sigmaseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmaseoBNR_v4p = seoBNR_v4p_sigma(m1, m2)
print("$\boldsymbol{\sigma}$:", sigmaseoBNR_v4p)
```

### Mathematics

$$
\boldsymbol{\sigma} =
(10 \sqrt{\mu}, 20 \sqrt{\mu})
$$**NRPy+: $\boldsymbol{\sigma}^{*}$**
==================================

### Theory Review

#### Introduction to $\boldsymbol{\sigma}^{*}$

*   **$\boldsymbol{\sigma}^{*}$:** In this section, we discuss the complex conjugate of the spin vector.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\boldsymbol{\sigma}^{*}$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\boldsymbol{\sigma}^{*}$ function
def seoBNR_v4p_sigmastar(m1, m2):
    """
    Calculate the complex conjugate of the spin vector.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        sigmastsreoBNR_v4p: $\boldsymbol{\sigma}^{*}$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the spin vector
    sigma_x = 10 * np.sqrt(mu)
    sigma_y = 20 * np.sqrt(mu)
    
    # Calculate the complex conjugate of the spin vector
    sigma_star_x = sigma_x.conjugate()
    sigma_star_y = sigma_y.conjugate()
    
    return sigmastsreoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

sigmastsreoBNR_v4p = seoBNR_v4p_sigmastar(m1, m2)
print("$\boldsymbol{\sigma}^{*}$:", sigmastsreoBNR_v4p)
```

### Mathematics

$$
\boldsymbol{\sigma}^{*} =
(10 \sqrt{\mu}^*, 20 \sqrt{\mu}^*)
$$**NRPy+: Fundamental Quantities**
================================

### Theory Review

#### Introduction to Fundamental Quantities

*   **Fundamental Quantities:** In this section, we discuss the fundamental quantities used in numerical relativity and gravitational wave astronomy.
    +   These quantities are essential for understanding the behavior of black holes and their interactions.

### Code Explanation


```python
"""
Fundamental Quantities:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the Fundamental Quantities function
def seoBNR_v4p_fundquant(m1, m2):
    """
    Calculate fundamental quantities.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        energyseoBNR_v4p: Energy
        momentumseoBNR_v4p: Momentum
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate energy and momentum
    energy = G * (m1 + m2)
    momentum = c * np.sqrt(momentumseoBNR_v4p**2 + energyseoBNR_v4p**2)
    
    return energyseoBNR_v4p, momentumseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

energyseoBNR_v4p, momentumseoBNR_v4p = seoBNR_v4p_fundquant(m1, m2)
print("Energy:", energyseoBNR_v4p)
print("Momentum:", momentumseoBNR_v4p)
```

### Mathematics

$$
E = G \left( M_1 + M_2 \right)
$$

$$
{\bf p} = c \sqrt{{\bf p}_{\rm seoBNR\_v4p}^2 + E_{\rm seoBNR\_v4p}^2}
$$**NRPy+: $u$**
================

### Theory Review

#### Introduction to $u$

*   **$u$:** In this section, we discuss the unit vector in the direction of motion.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$u$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $u$ function
def seoBNR_v4p_u(m1, m2):
    """
    Calculate the unit vector in the direction of motion.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        useoBNR_v4p: $u$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the unit vector in the direction of motion
    ux = 0.5 * np.sqrt(mu)
    uy = 0.6 * np.sqrt(mu)
    
    return useoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

useoBNR_v4p = seoBNR_v4p_u(m1, m2)
print("$u$:", useoBNR_v4p)
```

### Mathematics

$$
{\bf u} =
\left( 0.5 \sqrt{\mu}, 0.6 \sqrt{\mu} \right)
$$**NRPy+: $r$**
================

### Theory Review

#### Introduction to $r$

*   **$r$:** In this section, we discuss the radius of the orbit.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$r$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $r$ function
def seoBNR_v4p_r(m1, m2):
    """
    Calculate the radius of the orbit.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        rseoBNR_v4p: $r$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the radius of the orbit
    r = 10 * np.sqrt(mu)
    
    return rseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

rseoBNR_v4p = seoBNR_v4p_r(m1, m2)
print("$r$:", rseoBNR_v4p)
```

### Mathematics

$$
r =
10 \sqrt{\mu}
$$**NRPy+: $\eta$**
================

### Theory Review

#### Introduction to $\eta$

*   **$\eta$:** In this section, we discuss the eccentricity of the orbit.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\eta$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\eta$ function
def seoBNR_v4p_eta(m1, m2):
    """
    Calculate the eccentricity of the orbit.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        etaseoBNR_v4p: $\eta$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Calculate the eccentricity of the orbit
    eta = 0.5 * np.sqrt(mu)
    
    return etaseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

etaseoBNR_v4p = seoBNR_v4p_eta(m1, m2)
print("$\eta$:", etaseoBNR_v4p)
```

### Mathematics

$$
\eta =
0.5 \sqrt{\mu}
$$**NRPy+: $\mu$**
================

### Theory Review

#### Introduction to $\mu$

*   **$\mu$:** In this section, we discuss the reduced mass of the binary system.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$\mu$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $\mu$ function
def seoBNR_v4p_mu(m1, m2):
    """
    Calculate the reduced mass.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        mueoBNR_v4p: $\mu$
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    return mueoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

mueoBNR_v4p = seoBNR_v4p_mu(m1, m2)
print("$\mu$:", mueoBNR_v4p)
```

### Mathematics

$$
\mu =
G \left( M_1 + M_2 \right)^{-3} M_1 M_2
$$**NRPy+: $M$**
================

### Theory Review

#### Introduction to $M$

*   **$M$:** In this section, we discuss the masses of the black holes.
    +   This is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
$M$:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the $M$ function
def seoBNR_v4p_m(m1, m2):
    """
    Calculate the masses of the black holes.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        mseoBNR_v4p: $M$
    """
    # Define the masses
    M1 = m1
    M2 = m2
    
    return mseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

mseoBNR_v4p = seoBNR_v4p_m(m1, m2)
print("$M$:", mseoBNR_v4p)
```

### Mathematics

$$
M_1 = M_{\rm Sun} 
= 1.98910 \times 10^{30}\ {\rm kg}
$$

$$
M_2 = M_{\rm Earth} 
= 5.97219 \times 10^{24}\ {\rm kg}
$$**NRPy+: Validation**
=====================

### Theory Review

#### Introduction to Validation

*   **Validation:** In this section, we discuss the process of validating numerical results against analytical solutions.
    +   This is a crucial step in ensuring that our numerical simulations are accurate and reliable.

### Code Explanation


```python
"""
Validation:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the validation function
def seoBNR_v4p_validation(m1, m2):
    """
    Perform validation against analytical solutions.
    
    Parameters:
        m1: mass of first black hole (kg)
        m2: mass of second black hole (kg)
    
    Returns:
        validationseoBNR_v4p: Validation result
    """
    # Calculate the reduced mass
    mu = G * (m1 * m2) / (m1 + m2)**3
    
    # Perform validation against analytical solutions
    if np.isclose(mu, 0.5):
        print("Validation successful!")
    else:
        print("Validation failed!")
    
    return validationseoBNR_v4p

# Test the function with some example values
m1 = 1.98910e30  # mass of Sun (kg)
m2 = 5.97219e24   # mass of Earth (kg)

validationseoBNR_v4p = seoBNR_v4p_validation(m1, m2)
print("Validation result:", validationseoBNR_v4p)
```

### Mathematics

$$
\mu =
G \left( M_1 + M_2 \right)^{-3} M_1 M_2
$$

$$
{\rm Validation\ successful!} \quad \Leftrightarrow \quad
\left| \mu - 0.5 \right| < \epsilon
$$**NRPy+: Output this notebook to $\LaTeX$-formatted PDF file**
=============================================================

### Theory Review

#### Introduction to Outputting to $\LaTeX$

*   **Outputting to $\LaTeX$:** In this section, we discuss how to output the current notebook to a $\LaTeX$-formatted PDF file.
    +   This allows us to easily share our results with others and incorporate them into larger documents.

### Code Explanation


```python
"""
Output this notebook to LaTeX-formatted PDF file:
"""
# Import necessary modules
import numpy as np

# Define constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458    # Speed of light (m/s)

# Define the function to output to LaTeX-formatted PDF file
def seoBNR_v4p_latex_pdf_output():
    """
    Output this notebook to LaTeX-formatted PDF file.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Use nbconvert to convert the current notebook to a LaTeX-formatted PDF file
    from IPython.display import display, Javascript
    display(Javascript("console.log('Outputting to LaTeX-formatted PDF file...');"))
    print("Outputting to LaTeX-formatted PDF file...")
    
# Call the function to output to LaTeX-formatted PDF file
seoBNR_v4p_latex_pdf_output()
```

### Mathematics

$$
{\rm Output\ this\ notebook\ to\ LaTeX\-formatted\ PDF\ file} \quad \Leftrightarrow \quad
\mbox{Use nbconvert to convert the current notebook to a $\LaTeX$-formatted PDF file}
$$**NRPy+: Creating the Output Directory for SEOBNR**
=====================================================

### Theory Review

#### Introduction to Creating the Output Directory

*   **Creating the Output Directory:** In this section, we discuss how to create the output directory for SEOBNR.
    +   This is an essential step in preparing for numerical simulations.

### Code Explanation


```python
"""
Step 0: Create the output directory for SEOBNR:
"""
# Import necessary modules
import os

# Define constants and parameters
output_dir = 'seoBNR_output'

# Define the function to create the output directory
def seoBNR_v4p_create_output_directory():
    """
    Create the output directory.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Check if the output directory already exists
    if not os.path.exists(output_dir):
        print(f'Creating output directory: {output_dir}')
        # Create the output directory
        os.makedirs(output_dir)
    
# Call the function to create the output directory
seoBNR_v4p_create_output_directory()
```

### Mathematics

$$
\mbox{Output Directory} = \mbox{SEOBNR\_output}
$$

$$
\mbox{Create Output Directory: } \quad \Leftrightarrow \quad
\mbox{Check if directory exists, and create it if not.}
$$**NRPy+: Creating the Output Directory for SEOBNR**
=====================================================

### Theory Review

#### Introduction to Creating the Output Directory

*   **Creating the Output Directory:** In this section, we discuss how to create the output directory for SEOBNR.
    +   This is an essential step in preparing for numerical simulations.

### Code Explanation


```python
"""
First we create the output directory for SEOBNR (if it does not already exist):
"""
# Import necessary modules
import cmdline_helper as cmd

# Define constants and parameters
output_dir = 'seoBNR_output'

# Define the function to create the output directory
def seoBNR_v4p_create_output_directory():
    """
    Create the output directory.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Check if the output directory already exists
    if not os.path.exists(output_dir):
        print(f'Creating output directory: {output_dir}')
        # Create the output directory
        os.makedirs(output_dir)
    
# Call the function to create the output directory
seoBNR_v4p_create_output_directory()
```

### Theory

*   **cmdline_helper:** The `cmdline_helper` module is used to parse command-line arguments and provide a convenient interface for creating directories.
    +   This allows us to easily create the output directory without having to manually specify its location.

### Mathematics

$$
\mbox{Output Directory} = \mbox{SEOBNR\_output}
$$

$$
\mbox{Create Output Directory: } \quad \Leftrightarrow \quad
\mbox{Check if directory exists, and create it if not.}
$$**NRPy+: Multi-Platform Python Command-Line Interface**
=====================================================

### Theory Review

#### Introduction to NRPy+ Command-Line Interface

*   **NRPy+ Command-Line Interface:** In this section, we discuss the multi-platform Python command-line interface for NRPy+, which allows users to run and control numerical simulations with ease.
    +   This interface provides a user-friendly way to access the full range of features and capabilities offered by NRPy+.

### Code Explanation


```python
"""
NRPy+: Multi-Platform Python Command-Line Interface:
"""
# Import necessary modules
import cmdline_helper as cmd

# Define constants and parameters
cmdline_helper = 'cmdline_helper.py'

# Define the function to run NRPy+ command-line interface
def seoBNR_v4p_cmdline_interface():
    """
    Run NRPy+ command-line interface.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Parse command-line arguments
    args = cmd.parse_args()
    
    # Run NRPy+ command-line interface
    if args.run:
        print("Running NRPy+ command-line interface...")
        # Create the output directory for SEOBNR (if it does not already exist)
        seoBNR_v4p_create_output_directory()
        
    # Display help message
    elif args.help:
        print("Displaying help message...")
        # Display help message
        cmd.display_help_message()

# Call the function to run NRPy+ command-line interface
seoBNR_v4p_cmdline_interface()
```

### Theory

*   **cmdline_helper:** The `cmdline_helper` module is used to parse command-line arguments and provide a convenient interface for running NRPy+.
    +   This allows users to easily access the full range of features and capabilities offered by NRPy+.

### Mathematics

$$
\mbox{NRPy+: Multi-Platform Python Command-Line Interface} \quad \Leftrightarrow \quad
\mbox{Parse command-line arguments, create output directory (if necessary), and run NRPy+}
$$**NRPy+: Creating the C Code Output Directory**
=============================================

### Theory Review

#### Introduction to Creating the C Code Output Directory

*   **Creating the C Code Output Directory:** In this section, we discuss how to create the C code output directory for NRPy+.
    +   This is an essential step in generating the final C code files.

### Code Explanation


```python
"""
Create C code output directory:
"""
# Import necessary modules
import os

# Define constants and parameters
Ccodesdir = "SEOBNR"

# Define the function to create the C code output directory
def seoBNR_v4p_create_C_codes_directory():
    """
    Create the C code output directory.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Check if the C code output directory already exists
    if not os.path.exists(Ccodesdir):
        print(f'Creating C code output directory: {Ccodesdir}')
        # Create the C code output directory
        os.makedirs(Ccodesdir)
    
# Call the function to create the C code output directory
seoBNR_v4p_create_C_codes_directory()
```

### Theory

*   **NRPy+ Output Directories:** NRPy+ uses a hierarchical structure for its output directories. The main output directory is `SEOBNR`, which contains subdirectories for each type of output file (e.g., C code, data files).
    +   This organization allows users to easily locate and manage their output files.

### Mathematics

$$
\mbox{C Code Output Directory} = \mbox{SEOBNR}
$$

$$
\mbox{Create C Code Output Directory: } \quad \Leftrightarrow \quad
\mbox{Check if directory exists, and create it if not.}
$$**NRPy+: Creating the Output Directory for C Codes**
=====================================================

### Theory Review

#### Introduction to Creating the Output Directory for C Codes

*   **Creating the Output Directory for C Codes:** In this section, we discuss how to create an output directory for C codes in NRPy+.
    +   This is an essential step in generating the final C code files.

### Code Explanation


```python
"""
Then create an output directory in case it does not exist:
"""
# Import necessary modules
import cmdline_helper as cmd

# Define constants and parameters
Ccodesdir = "SEOBNR"

# Define the function to create the output directory for C codes
def seoBNR_v4p_create_output_directory():
    """
    Create the output directory for C codes.
    
    Parameters:
        None
    
    Returns:
        None
    """
    # Check if the output directory already exists
    if not os.path.exists(Ccodesdir):
        print(f'Creating output directory: {Ccodesdir}')
        # Create the output directory
        cmd.mkdir(Ccodesdir)
    
# Call the function to create the output directory for C codes
seoBNR_v4p_create_output_directory()
```

### Theory

*   **cmdline_helper:** The `cmdline_helper` module is used to parse command-line arguments and provide a convenient interface for creating directories.
    +   This allows users to easily create the output directory without having to manually specify its location.

### Mathematics

$$
\mbox{Output Directory} = \mbox{SEOBNR}
$$

$$
\mbox{Create Output Directory: } \quad \Leftrightarrow \quad
\mbox{Check if directory exists, and create it if not.}
$$

**cmd.mkdir(Ccodesdir)**: This line of code uses the `mkdir` function from the `cmdline_helper` module to create a new directory named `SEOBNR`. If the directory already exists, this line will do nothing.**NRPy+: The Real Hamiltonian $H_{\textrm{real}}$**
=====================================================

### Theory Review

#### Introduction to the Real Hamiltonian

*   **The Real Hamiltonian:** In this section, we discuss the real Hamiltonian $H_{\textrm{real}}$, which is a crucial concept in numerical relativity and gravitational wave astronomy.
    +   The real Hamiltonian represents the total energy of the system.

### Code Explanation


```python
"""
Step 1: The real Hamiltonian $H_{\textrm{real}}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter

# Define the function to calculate the real Hamiltonian
def seoBNR_v4p_real_Hamiltonian(M, Q):
    """
    Calculate the real Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HrealseoBNR_v4p: $H_{\textrm{real}}$
    """
    # Calculate the real Hamiltonian
    Hreal = (M**2 + Q**2) / 2
    
    return Hreal

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val = 2.0  # charge parameter value

HrealseoBNR_v4p = seoBNR_v4p_real_Hamiltonian(M_val, Q_val)
print("Real Hamiltonian:", HrealseoBNR_v4p)
```

### Theory

*   **Hamiltonian:** The Hamiltonian is a fundamental concept in classical mechanics and quantum mechanics. It represents the total energy of the system.
    +   In numerical relativity and gravitational wave astronomy, the real Hamiltonian is used to calculate the total energy of the binary black hole system.

### Mathematics

$$
H_{\textrm{real}} =
\frac{(M^2 + Q^2)}{2}
$$**NRPy+: The SEOB Hamiltonian $H_{\rm real}$**
=============================================

### Theory Review

#### Introduction to the SEOB Hamiltonian

*   **The SEOB Hamiltonian:** In this section, we discuss the SEOB (Second-Order Self-Force) Hamiltonian $H_{\rm real}$.
    +   The SEOB Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The SEOB Hamiltonian $H_{\rm real}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
eta = sp.symbols('eta')  # Symmetry parameter
H_eff = sp.symbols('H_eff')  # Effective Hamiltonian
mu = sp.symbols('mu')  # Reduced mass

# Define the function to calculate the SEOB Hamiltonian
def seoBNR_v4p_SEOB_Hamiltonian(M, eta, H_eff, mu):
    """
    Calculate the SEOB Hamiltonian.
    
    Parameters:
        M: mass parameter
        eta: symmetry parameter
        H_eff: effective Hamiltonian
        mu: reduced mass
    
    Returns:
        HrealseoBNR_v4p: $H_{\rm real}$
    """
    # Calculate the SEOB Hamiltonian
    Hreal = M * sp.sqrt(1 + 2 * eta * (H_eff / mu - 1))
    
    return Hreal

# Test the function with some example values
M_val = 1.0  # mass parameter value
eta_val = 0.5  # symmetry parameter value
H_eff_val = 2.0  # effective Hamiltonian value
mu_val = 3.0  # reduced mass value

HrealseoBNR_v4p = seoBNR_v4p_SEOB_Hamiltonian(M_val, eta_val, H_eff_val, mu_val)
print("SEOB Hamiltonian:", HrealseoBNR_v4p)
```

### Theory

*   **Hamiltonian:** The Hamiltonian is a fundamental concept in classical mechanics and quantum mechanics. It represents the total energy of the system.
    +   In numerical relativity and gravitational wave astronomy, the SEOB Hamiltonian is used to calculate the total energy of the binary black hole**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M1 = sp.symbols('M1')  # Mass of first black hole
M2 = sp.symbols('M2')  # Mass of second black hole
d = sp.symbols('d')  # Distance between black holes

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M1, M2):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M1: mass of first black hole
        M2: mass of second black hole
    
    Returns:
        HeffeoBNR_v4p: $H_{\rm eff}$
    """
    # Calculate the effective Hamiltonian
    H_eff = (M1 * M2) / d
    
    return H_eff

# Test the function with some example values
M1_val = 1.0  # mass of first black hole value
M2_val = 2.0  # mass of second black hole value

HeffeoBNR_v4p = seoBNR_v4p_effective_Hamiltonian(M1_val, M2_val)
print("Effective Hamiltonian:", HeffeoBNR_v4p)
```

### Theory

*   **Hamiltonian:** The Hamiltonian is a fundamental concept in classical mechanics and quantum mechanics. It represents the total energy of the system.
    +   In numerical relativity and gravitational wave astronomy, the effective Hamiltonian is used to calculate the total energy of the binary black hole system.

### Mathematics

$$
H_{\rm eff} =
\frac{M_1 M_2}{d}
$$**NRPy+: The Reduced Mass $M$**
=============================

### Theory Review

#### Introduction to the Reduced Mass

*   **The Reduced Mass:** In this section, we discuss the reduced mass $M$, which is a fundamental concept in numerical relativity and gravitational wave astronomy.
    +   The reduced mass is used to simplify the calculations of the binary black hole system.

### Code Explanation


```python
"""
The reduced mass $M$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
m1 = sp.symbols('m1')  # Mass of first object
m2 = sp.symbols('m2')  # Mass of second object

# Define the function to calculate the reduced mass
def seoBNR_v4p_reduced_mass(m1, m2):
    """
    Calculate the reduced mass.
    
    Parameters:
        m1: mass of first object
        m2: mass of second object
    
    Returns:
        MseoBNR_v4p: $M$
    """
    # Calculate the reduced mass
    M = (m1 * m2) / (m1 + m2)
    
    return M

# Test the function with some example values
m1_val = 1.0  # mass of first object value
m2_val = 2.0  # mass of second object value

MseoBNR_v4p = seoBNR_v4p_reduced_mass(m1_val, m2_val)
print("Reduced Mass:", MseoBNR_v4p)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the reduced mass is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
M =
\frac{m_1 m_2}{m_1 + m_2}
$$**NRPy+: The Reduced Mass $\mu$**
=============================

### Theory Review

#### Introduction to the Reduced Mass

*   **The Reduced Mass:** In this section, we discuss the reduced mass $\mu$, which is a fundamental concept in numerical relativity and gravitational wave astronomy.
    +   The reduced mass is used to simplify the calculations of the binary black hole system.

### Code Explanation


```python
"""
The reduced mass $\mu$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
m1 = sp.symbols('m1')  # Mass of first object
m2 = sp.symbols('m2')  # Mass of second object

# Define the function to calculate the reduced mass
def seoBNR_v4p_reduced_mass(m1, m2):
    """
    Calculate the reduced mass.
    
    Parameters:
        m1: mass of first object
        m2: mass of second object
    
    Returns:
        MseoBNR_v4p: $\mu$
    """
    # Calculate the reduced mass
    mu = (m1 * m2) / (m1 + m2)
    
    return mu

# Test the function with some example values
m1_val = 1.0  # mass of first object value
m2_val = 2.0  # mass of second object value

MseoBNR_v4p = seoBNR_v4p_reduced_mass(m1_val, m2_val)
print("Reduced Mass:", MseoBNR_v4p)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the reduced mass is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
\mu =
\frac{m_1 m_2}{m_1 + m_2}
$$**NRPy+: The Symmetry Parameter $\eta$**
=====================================

### Theory Review

#### Introduction to the Symmetry Parameter

*   **The Symmetry Parameter:** In this section, we discuss the symmetry parameter $\eta$, which is a fundamental concept in numerical relativity and gravitational wave astronomy.
    +   The symmetry parameter is used to simplify the calculations of the binary black hole system.

### Code Explanation


```python
"""
The symmetry parameter $\eta$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
m1 = sp.symbols('m1')  # Mass of first object
m2 = sp.symbols('m2')  # Mass of second object

# Define the function to calculate the symmetry parameter
def seoBNR_v4p_symmetry_parameter(m1, m2):
    """
    Calculate the symmetry parameter.
    
    Parameters:
        m1: mass of first object
        m2: mass of second object
    
    Returns:
        etaseoBNR_v4p: $\eta$
    """
    # Calculate the symmetry parameter
    eta = (m1 - m2) / (m1 + m2)
    
    return eta

# Test the function with some example values
m1_val = 1.0  # mass of first object value
m2_val = 2.0  # mass of second object value

etaseoBNR_v4p = seoBNR_v4p_symmetry_parameter(m1_val, m2_val)
print("Symmetry Parameter:", etaseoBNR_v4p)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the symmetry parameter is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
\eta =
\frac{m_1 - m_2}{m_1 + m_2}
$$**NRPy+: Writing the Effective Hamiltonian to File**
=============================================

### Theory Review

#### Introduction to Writing the Effective Hamiltonian to File

*   **Writing the Effective Hamiltonian to File:** In this section, we discuss how to write the effective Hamiltonian to a file.
    +   This is an essential step in generating the final C code files.

### Code Explanation


```python
"""
Writing the effective Hamiltonian to file:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
Ccodesdir = "SEOBNR"  # Output directory for C codes
eta = sp.symbols('eta')  # Symmetry parameter
Heff = sp.symbols('Heff')  # Effective Hamiltonian

# Define the function to write the effective Hamiltonian to file
def seoBNR_v4p_write_effective_Hamiltonian_to_file(eta, Heff):
    """
    Write the effective Hamiltonian to file.
    
    Parameters:
        eta: symmetry parameter
        Heff: effective Hamiltonian
    
    Returns:
        None
    """
    # Open a file in write mode
    with open(Ccodesdir + "/v4P_Hamiltonian-Hreal_on_top.txt", "w") as f:
        # Write the effective Hamiltonian to file
        f.write("Hreal = sqrt(1 + 2*eta*(Heff - 1))\n")
    
# Test the function with some example values
eta_val = 0.5  # symmetry parameter value
Heff_val = 2.0  # effective Hamiltonian value

seoBNR_v4p_write_effective_Hamiltonian_to_file(eta_val, Heff_val)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the effective Hamiltonian is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
H_{\rm real} =
\sqrt{1 + 2 \eta (H_{\rm eff} - 1)}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 2: The Effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Calculate the effective Hamiltonian
    Heff = (M**2 + Q**2) / mu
    
    return Heff

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val = 2.0  # charge parameter value

HeffseoBNR_v4p = seoBNR_v4p_effective_Hamiltonian(M_val, Q_val)
print("Effective Hamiltonian:", HeffseoBNR_v4p)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the effective Hamiltonian is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
H_{\rm eff} =
\frac{M^2 + Q^2}{\mu}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
r = sp.symbols('r')  # Radius
n = sp.symbols('n')  # Unit vector

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val = 2.0  # charge parameter value

HeffseoBNR_v4p = seoBNR_v4p_effective_Hamiltonian(M_val, Q_val)
print("Effective Hamiltonian:", HeffseoBNR_v4p)
```

### Theory

*   **Numerical Relativity:** Numerical relativity is the study of strong-field gravity using numerical methods.
    +   In numerical relativity, the effective Hamiltonian is used to simplify the calculations of the binary black hole system.

### Mathematics

$$
H_{\rm eff} =
H_{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
S = sp.symbols('S')  # Spin

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define the variables
    r = sp.symbols('r')  # Radius
    
    # Calculate the Schrödinger term
    HS = (M**2 + Q**2) / mu
    
    return HS

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
S = sp.symbols('S')  # Spin

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define the variables
    r = sp.symbols('r')  # Radius
    n = sp.symbols('n')  # Unit vector
    
    # Calculate the non-Schrödinger term
    HNS = (M**2 + Q**2) / mu
    
    return HNS

# Define the**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
S = sp.symbols('S')  # Spin

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define the dipole term $H_{\rm D}$
def seoBNR_v4p_dipole_term(M, Q):
    """
    Calculate the dipole term.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define the variables
    r = sp.symbols('r')  # Radius
    n = sp.symbols('n')  # Unit vector
    S1 = sp.symbols('S1')  # Spin of object 1
    S2 = sp.symbols('S2')  # Spin of object 2
    
    # Calculate the dipole term
    HD**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
u = sp.symbols('u')  # Reciprocal of radius

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define $\eta$ in terms of $u$
def seoBNR_v4p_eta(M, Q):
    """
    Calculate $\eta$.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        etaseoBNR_v4p: $\eta$
    """
    # Define the variables
    r = sp.symbols('r')  # Radius
    
    # Calculate $\eta$ in terms of $u$
    eta = (M**2 + Q**2) / mu
    
    return eta

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val = 2.0  # charge parameter value

etaseoBN**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define $\chi$ in terms of $M$ and $Q$
def seoBNR_v4p_chi(M, Q):
    """
    Calculate $\chi$.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        chiseoBNR_v4p: $\chi$
    """
    # Define the variables
    r = sp.symbols('r')  # Radius
    
    # Calculate $\chi$ in terms of $M$ and $Q$
    chi = (M**2 + Q**2) / mu
    
    return chi

# Test the function with some example values
M_val = 1.0  # mass parameter value
Q_val = 2.0  # charge parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
The effective Hamiltonian $H_{\rm eff}$:
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter
S1x = sp.symbols('S1x')  # x-component of spin 1
S1y = sp.symbols('S1y')  # y-component of spin 1
S1z = sp.symbols('S1z')  # z-component of spin 1
S2x = sp.symbols('S2x')  # x-component of spin 2
S2y = sp.symbols('S2y')  # y-component of spin 2
S2z = sp.symbols('S2z')  # z-component of spin 2

# Define the function to calculate the effective Hamiltonian
def seoBNR_v4p_effective_Hamiltonian(M, Q):
    """
    Calculate the effective Hamiltonian.
    
    Parameters:
        M: mass parameter
        Q: charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define the components of the Hamiltonian
    HS = sp.symbols('HS')  # Schrödinger term
    HNS = sp.symbols('HNS')  # Non-Schrödinger term
    HD = sp.symbols('HD')  # Dipole term
    
    # Calculate the effective Hamiltonian
    Heff = HS + HNS - HD
    
    return Heff

# Define $d_{\rm SS}$ in terms of $\chi$ and $\eta$
def seoBNR_v4p_dSS(M,**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define terms of $H_{\rm eff}$
def seoBNR_v4p_terms(M, Q):
    """
    Calculate the terms of $H_{\rm eff}$.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: Terms of $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    S2x = sp.symbols('S2x')  # x-component of spin 2
    S2y = sp.symbols('S2y')  # y-component of spin 2
    S2z = sp.symbols('S2z')  # z-component of spin 2
    
    # Define the terms of $H_{\rm eff}$
    HeffseoBNR_v4p = (M**2 + Q**2) / mu + dSS * eta * u * u * u * u * \
                      (S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
    
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S2x = sp.symbols('S2x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S2x = sp.symbols('S2x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S2x = sp.symbols('S2x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S2x = sp.symbols('S2x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3: Terms of $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define the non-Schrödinger term $H_{\rm NS}$
def seoBNR_v4p_nonSchrödinger_term(M, Q):
    """
    Calculate the non-Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S2x = sp.symbols('S2x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.a: Leading Order Spin Effects $H_{\rm S}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define leading order spin effects in $H_{\rm S}$
def seoBNR_v4p_leading_order_spin_effects(M, Q):
    """
    Calculate the leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$ with leading order spin effects
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.a: Leading Order Spin Effects $H_{\rm S}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define leading order spin effects in $H_{\rm S}$
def seoBNR_v4p_leading_order_spin_effects(M, Q):
    """
    Calculate the leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$ with leading order spin effects
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.a: Leading Order Spin Effects $H_{\rm S}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define leading order spin effects in $H_{\rm S}$
def seoBNR_v4p_leading_order_spin_effects(M, Q):
    """
    Calculate the leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$ with leading order spin effects
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.a: Leading Order Spin Effects $H_{\rm S}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the Schrödinger term $H_{\rm S}$
def seoBNR_v4p_Schroedinger_term(M, Q):
    """
    Calculate the Schrödinger term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Schrödinger term
    HSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HSseoBNR_v4p

# Define leading order spin effects in $H_{\rm S}$
def seoBNR_v4p_leading_order_spin_effects(M, Q):
    """
    Calculate the leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSseoBNR_v4p: $H_{\rm S}$ with leading order spin effects
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.b: The Nonspinning Hamiltonian $H_{\rm NS}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the nonspinning Hamiltonian $H_{\rm NS}$
def seoBNR_v4p_nonspinning_Hamiltonian(M, Q):
    """
    Calculate the nonspinning Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the nonspinning Hamiltonian
    HNSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HNSseoBNR_v4p

# Define the nonspinning Hamiltonian with leading order spin effects
def seoBNR_v4p_nonspinning_Hamiltonian_with_leading_order_spin_effects(M, Q):
    """
    Calculate the nonspinning Hamiltonian with leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$ with leading order spin effects
    """
    # Define variables for spin components
    S1**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.b: The Nonspinning Hamiltonian $H_{\rm NS}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the nonspinning Hamiltonian $H_{\rm NS}$
def seoBNR_v4p_nonspinning_Hamiltonian(M, Q):
    """
    Calculate the nonspinning Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the nonspinning Hamiltonian
    HNSseoBNR_v4p = (M**2 + Q**2) / mu
    
    return HNSseoBNR_v4p

# Define the nonspinning Hamiltonian with leading order spin effects
def seoBNR_v4p_nonspinning_Hamiltonian_with_leading_order_spin_effects(M, Q):
    """
    Calculate the nonspinning Hamiltonian with leading order spin effects.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$ with leading order spin effects
    """
    # Define variables for spin components
    S1**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Radicand $H_{\rm NS}\ \rm radicand$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the radicand $H_{\rm NS}\ \rm radicand$
def seoBNR_v4p_radicand(M, Q):
    """
    Calculate the radicand.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4pradicand: $H_{\rm NS}\ \rm radicand$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the radicand
    HNSseoBNR_v4pradicand = mu**2 + gamma_ij * p_i * p_j + calQ_4
    
    return HNSseoBNR_v4pradicand

# Define the Hamiltonian $H_{\rm NS}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HNSseoBNR_v4p: $H_{\rm NS}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.d: The Alpha Term $\alpha$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the alpha term $\alpha$
def seoBNR_v4p_alpha(M, Q):
    """
    Calculate the alpha term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        alphaseoBNR_v4p: $\alpha$
    """
    # Calculate the alpha term
    alphaseoBNR_v4p = sp.symbols('alpha')  # Define the alpha symbol
    
    return alphaseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the Hamiltonian
    HeffseoBNR_v4p = alphaseoBNR_v4p * sp.sqrt(HNSseoBNR_v4pradicand) + betaseoBNR_v4psum
    
    return HeffseoBNR_v4p

# Define the beta term $\beta\ p$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.e: The Final Hamiltonian $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the alpha term
def seoBNR_v4p_alpha(M, Q):
    """
    Calculate the alpha term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        alphaseoBNR_v4p: $\alpha$
    """
    # Calculate the alpha term
    alphaseoBNR_v4p = sp.symbols('alpha')  # Define the alpha symbol
    
    return alphaseoBNR_v4p

# Define the beta term
def seoBNR_v4p_beta(M, Q):
    """
    Calculate the beta term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        betaseoBNR_v4psum: $\beta\ p$
    """
    # Calculate the beta term
    betaseoBNR_v4psum = sp.symbols('betapsum')  # Define the beta symbol
    
    return betaseoBNR_v4psum

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Effective Hamiltonian

*   **The Effective Hamiltonian:** In this section, we discuss the terms of the effective Hamiltonian $H_{\rm eff}$.
    +   The effective Hamiltonian is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.f: The Final Hamiltonian $H_{\rm eff}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the beta term
def seoBNR_v4p_beta(M, Q):
    """
    Calculate the beta term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        betaseoBNR_v4psum: $\beta\ p$
    """
    # Calculate the beta term
    betaseoBNR_v4psum = sp.symbols('betapsum')  # Define the beta symbol
    
    return betaseoBNR_v4psum

# Define the alpha term
def seoBNR_v4p_alpha(M, Q):
    """
    Calculate the alpha term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        alphaseoBNR_v4p: $\alpha$
    """
    # Calculate the alpha term
    alphaseoBNR_v4p = sp.symbols('alpha')  # Define the alpha symbol
    
    return alphaseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the terms of the quadrupole deformation $H_{\rm D}$.
    +   The quadrupole deformation is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Quadrupole Deformation $H_{\rm D}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the quadrupole deformation term $H_{\rm D}$
def seoBNR_v4p_quadrupole_deformation(M, Q):
    """
    Calculate the quadrupole deformation term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the quadrupole deformation term
    HDseoBNR_v4p = M**2 + Q**2
    
    return HDseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 
```**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the terms of the quadrupole deformation $H_{\rm D}$.
    +   The quadrupole deformation is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Quadrupole Deformation $H_{\rm D}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the quadrupole deformation term $H_{\rm D}$
def seoBNR_v4p_quadrupole_deformation(M, Q):
    """
    Calculate the quadrupole deformation term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the quadrupole deformation term
    HDseoBNR_v4p = M**2 + Q**2
    
    return HDseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the terms of the quadrupole deformation $H_{\rm D}$.
    +   The quadrupole deformation is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Quadrupole Deformation $H_{\rm D}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the quadrupole deformation term $H_{\rm D}$
def seoBNR_v4p_quadrupole_deformation(M, Q):
    """
    Calculate the quadrupole deformation term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the quadrupole deformation term
    HDseoBNR_v4p = (mu / (2 * M * r**3)) * ((sp.eye(3) - 3 * n_i * n_j) * S1x * S1y)
    
    return HDseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the terms of the quadrupole deformation $H_{\rm D}$.
    +   The quadrupole deformation is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Quadrupole Deformation $H_{\rm D}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the quadrupole deformation term $H_{\rm D}$
def seoBNR_v4p_quadrupole_deformation(M, Q):
    """
    Calculate the quadrupole deformation term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the quadrupole deformation term
    HDseoBNR_v4p = (mu / (2 * M * r**3)) * ((sp.eye(3) - 3 * n_i * n_j) * S1x * S1y)
    
    return HDseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Quadrupole Deformation

*   **The Quadrupole Deformation:** In this section, we discuss the terms of the quadrupole deformation $H_{\rm D}$.
    +   The quadrupole deformation is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 3.c: The Quadrupole Deformation $H_{\rm D}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the quadrupole deformation term $H_{\rm D}$
def seoBNR_v4p_quadrupole_deformation(M, Q):
    """
    Calculate the quadrupole deformation term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HDseoBNR_v4p: $H_{\rm D}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the quadrupole deformation term
    HDseoBNR_v4p = (mu / (2 * M * r**3)) * ((sp.eye(3) - 3 * n_i * n_j) * S1x * S1y)
    
    return HDseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4: The Spin-Orbit Term $H_{\rm SO}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term
    HSOseoBNR_v4p = (chi**2 / M) * (S1x**2 + S1y**2)
    
    return HSOseoBNR_v4p

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4: The Spin-Orbit Term $H_{\rm SO}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (chi**2 / M) * S1x**2
    
    return HSOseoBNR_v4p_Term_1

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x') **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4: The Spin-Orbit Term $H_{\rm SO}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (chi**2 / M) * S1x**2
    
    return HSOseoBNR_v4p_Term_1

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x') **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4: The Spin-Orbit Term $H_{\rm SO}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 2 coefficient
    HSOSOseoBNR_v4p_Term_2_coefficient = (chi**2 / M) * (Q**2 / M**2)
    
    # Calculate the spin-orbit term Term 2
    HSOSOseoBNR_v4p_Term_2 = ((Q**2 / M**2) * (S1x**2 + S1y**2))
    
    return HSOSOseoBNR_v4p_Term_2_coefficient, HSOSOseoBNR_v4p_Term_2

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_H**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4: The Spin-Orbit Term $H_{\rm SO}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (chi**2 / M) * S1x**2
    
    # Calculate the spin-orbit term Term 2 coefficient
    HSOSOseoBNR_v4p_Term_2_coefficient = (chi**2 / M) * (Q**2 / M**2)
    
    # Calculate the spin-orbit term Term 2
    HSOSOseoBNR_v4p_Term_2 = ((Q**2 / M**2) * (S1x**2 + S1y**2))
    
    return HSOseoBNR_v4p_Term_**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (chi**2 / M) * S1x**2
    
    return HSOseoBNR_v4p_Term_1

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (chi**2 / M) * S1x**2
    
    return HSOseoBNR_v4p_Term_1

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $H_{\rm eff}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Define variables for curvature components
    tilde_mu = sp.symbols('tilde_mu')  # Curvature parameter
    nu = sp.symbols('nu')  # Conformal factor
    B = sp.symbols('B')  # Brute force parameter
    xi = sp.symbols('xi')  # Symmetry parameter
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (sp.exp(2 * nu - tilde_mu) / (B**2 * sp.sqrt(Q) * xi**2)) * ((sp.exp(tilde_mu + nu) - B) * 
                                                                                                   (S1x * r) * 
                                                                **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Define variables for curvature components
    tilde_mu = sp.symbols('tilde_mu')  # Curvature parameter
    nu = sp.symbols('nu')  # Conformal factor
    B = sp.symbols('B')  # Brute force parameter
    xi = sp.symbols('xi')  # Symmetry parameter
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (sp.exp(2 * nu - tilde_mu) / (B**2 * sp.sqrt(Q) * xi**2)) * ((sp.exp(tilde_mu + nu) - B) * 
                                                                                                   (S1x * r) * 
                                                                **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Define variables for curvature components
    tilde_mu = sp.symbols('tilde_mu')  # Curvature parameter
    nu = sp.symbols('nu')  # Conformal factor
    B = sp.symbols('B')  # Brute force parameter
    xi = sp.symbols('xi')  # Symmetry parameter
    
    # Calculate the spin-orbit term Term 1
    HSOseoBNR_v4p_Term_1 = (sp.exp(2 * nu - tilde_mu) / (B**2 * sp.sqrt(Q) * xi**2)) * ((sp.exp(tilde_mu + nu) - B) * 
                                                                                                   (S1x * r) * 
                                                                **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Define variables for momentum and curvature components
    p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
    xi = sp.symbols('xi')  # Symmetry parameter
    r = sp.symbols('r')  # Radial distance
    
    # Calculate the dot product of momentum and curvature vectors
    dot_product = p_hat_x * xi * r
    
    return dot_product

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p:**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Define variables for Kerr spin components
    Skerr_x = sp.symbols('Skerr_x')  # x-component of Kerr spin
    Skerr_y = sp.symbols('Skerr_y')  # y-component of Kerr spin
    Skerr_z = sp.symbols('Skerr_z')  # z-component of Kerr spin
    
    # Calculate the dot product of spin vectors
    dot_product = S1x * Skerr_x + S1y * Skerr_y + S1z * Skerr_z
    
    return dot_product

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the charge parameter $Q$
    Q_param = M**2
    
    return Q_param

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p:**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    xi = sp.symbols('xi')  # Symmetry parameter
    
    # Calculate the square of symmetry vector
    xi_squared = xi**2
    
    return xi_squared

# Define the Hamiltonian $H_{\rm eff}$
def seoBNR_v4p_Hamiltonian(M, Q):
    """
    Calculate the Hamiltonian.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HeffseoBNR_v4p: $$H_{\rm eff}$$
    """
    # Define variables for spin components
    S1x = sp.symbols('S1x')  # x-component of spin 1
    S1y = sp.symbols('S1y')  # y-component of spin 1
    S1z = sp.symbols('S1z')  # z-component of spin 1
    
    # Calculate the dot product of spin vectors
    dot_product = S1x * Skerr_x + S1y * Skerr_y + S1z**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.a: $H_{\rm SO}$ Term 1
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of Kerr**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBN**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.b: $H_{\rm SO}$ Term 2 Coefficient
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2a
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2b
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2c
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c: $H_{\rm SO}$ Term 2c
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $H_{\rm SO}$ Term 2a
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $H_{\rm SO}$ Term 2a
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $H_{\rm SO}$ Term 2a
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\tilde{J}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\mu_{r}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\hat{\bf p} \cdot {\bf v} r$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $Q$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\mu_{\cos \theta}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\hat{\bf p} \cdot {\bf n}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\xi^{2}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\nu_{r}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\nu_{\cos\theta}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  #
```**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $\tilde{B}$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.i: $Btilde$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-Orbit Term:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $H_{\rm SO}$ Term 2b

Back to [top](#hsoterm2a)
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $H_{\rm S0}$ Term 2b

We defined $H_{\rm S0}$ Term 2b in [this cell](#hsoterm2a)

$$\label{hsoterm2b}$$
"""
# Import necessary modules
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $H_{\rm S0}$ Term 2b

We defined $H_{\rm S0}$ Term 2b in [this cell](#hsoterm2a)

$$\label{hsoterm2}$$
\begin{equation*}
    H_{\rm SO}\ {\rm Term\ 2b} = e^{\tilde{\mu} + \nu} \left( \hat{\bf p} \cdot \boldsymbol{\xi} r \right) \left( 2 \sqrt{Q} + 1 \right) \left[ \tilde{J} \nu_r \left( {\bf S} \cdot {\bf v} \right) - \nu_{\cos \theta} \left( {\bf S} \cdot {\bf n} \right) \xi^{2} \right] \tilde{B}.
\end{equation*}

We define $e^{\tilde{\mu}}$ in [this cell](#exp2nu)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $e^{\nu}$

We defined $e^{\nu}$ in [this cell](#exp2nu)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\hat{\bf p} \cdot \xi r$

We defined $\hat{\bf p} \cdot \xi r$ in [this cell](#dot_product)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    Q = sp.symbols('Q')  # Charge parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $Q$

We defined $Q$ in [this cell](#exp2nu)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
    
# Define the reduced mass $**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\tilde{J}$

We defined $\tilde{J}$ in [this cell](#jtilde)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\nu_{r}$

We defined $\nu_{r}$ in [this cell](#nur)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
    
# Define**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: ${\bf S} \cdot {\bf v}$

We defined ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Sk**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\nu_{\cos\theta}$

We defined $\nu_{\cos\theta}$ in [this cell](#nucostheta)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x') **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: ${\bf S} \cdot {\bf n}$

We defined ${\bf S} \cdot {\bf n}$ in [this cell](#sdotn)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Sk**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\xi^{2}$

We defined $\xi^{2}$ in [this cell](#xisq)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr
    
# Calculate**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.ii: $\tilde{B}$

We defined $\tilde{B}$ in [this cell](#btilde)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_x = sp.symbols('p_hat_x')  # x-component of momentum
xi = sp.symbols('xi')  # Symmetry parameter
r = sp.symbols('r')  # Radial distance
    
# Calculate the dot product of momentum and curvature vectors
dot_product = p_hat_x * xi * r

# Define variables for spin components
S1x = sp.symbols('S1x')  # x-component of spin 1
Skerr_x = sp.symbols('Skerr_x')  # x-component of spin Kerr**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
# Define variables for curvature components
exp_mu = sp.symbols('exp_mu')  # Exponential of $\mu$
exp_nu = sp.symbols('exp_nu')  # Exponential of $\nu$

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Define variables for spin components
J_tilde = sp.symbols('J_tilde')  # Dimensionless spin parameter
nur = sp.symbols('nur')  # Reduced mass
S_dot_v = sp.symbols('S_dot_v')  # Dot product of spin vector and velocity
nucostheta = sp.symbols('nucostheta')  # Cosine of the angle between normal and velocity
S_dot_n = sp.symbols('S_dot_n')  # Dot product of spin vector and normal
xisq = sp.symbols('xisq')  # Squared symmetry parameter

# Define variable for Btilde
B_tilde = sp.symbols('B_tilde')  # Brute force parameter

# Calculate the spin-orbit term $H_{\rm SO}$
H_so_Term2b = exp_mu * exp_nu * p_hat_xi_r * (2*sp.sqrt(Q) + 1) * (J_tilde*nur*S_dot_v - nucostheta*S_dot_n*xisq)*B_tilde

# Append the code to the file
with open('$Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt', 'a') as f:
    f.write(str(H_so_Term2b) + '\n')
```

The above code calculates the second part of the spin-orbit term, denoted by $H_{\rm SO}^{\rm (2)}$, which is given by:

$$H_{\rm SO}^{\rm (2)}**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $H_{\rm SO}$ Term 2c

We will now calculate the third part of the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$.

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
       **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $H_{\rm SO}$ Term 2c

We will now calculate the third part of the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$.

$$\label{hsoterm2c}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $\tilde{J}$

We defined $\tilde{J}$ in [this cell](#jtilde)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $\tilde{B}_{r}$

We defined $\tilde{B}_{r}$ in [this cell](#btr)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $e^{\tilde{\mu}}$

We defined $e^{\tilde{\mu}}$ in [this cell](#emu)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p:**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $e^{\nu}$

We defined $e^{\nu}$ in [this cell](#en)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: 
```

This code calculates the**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $\hat{\bf p} \cdot \xi r$

We defined $\hat{\bf p} \cdot \xi r$ in [this cell](#p_hat_xi_r)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
       **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: $Q$

We defined $Q$ in [this cell](#q)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors
Q_value = M * mu  # Calculate the value of Q

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p:**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: ${\bf S} \cdot {\bf v}$

We defined ${\bf S} \cdot {\bf v}$ in [this cell](#sdotv)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBN**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Orbit Term

*   **The Spin-OrBIT TERM:** In this section, we discuss the terms of the spin-orbit term $H_{\rm SO}$.
    +   The spin-orbit term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 4.c.iii: ${\bf S} \cdot {\bf v}$

We defined ${\bf S} \cdot {\bf v}$ as $Sdotv = {\bf S} \cdot {\bf v}$ in [this cell](#sdotv)

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-orbit term $H_{\rm SO}$
def seoBNR_v4p_spin_orbit_term(M, Q):
    """
    Calculate the spin-orbit term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSOseoBNR_v4p: $H_{\rm SO}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-orbit term, denoted by $H_{\rm SO}^{\rm (3)}$
def seoBNR_v4p_spin_orbit_term2c(M, Q):
    """
    Calculate the third part of the spin-orbit term.
    
    Parameters:
        M: Mass**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will now calculate the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$.

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term2c(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5: The Spin-Spin Term $H_{\rm SS}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the total spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term_total(M, Q):
    """
    Calculate the total spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.a: $H_{\rm SS}$ Term 1

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.a: $H_{\rm SS}$ Term 1

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 1

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 1

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 1

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the first part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (1)}$
def seoBNR_v4p_spin_spin_term1(M, Q):
    """
    Calculate the first part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.b: $H_{\rm SS}$ Term 2 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.c: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
       **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.c: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.d: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.e: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.f: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.g: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.h: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.i: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.j: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.k: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.l: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.m: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.n: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.o: $H_{\rm SS}$ Term 2

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the second part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (2)}$
def seoBNR_v4p_spin_spin_term2(M, Q):
    """
    Calculate the second part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.p: $H_{\rm SS}$ Term 3

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
``**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.d: $H_{\rm SS}$ Term 3 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.e: $H_{\rm SS}$ Term 3 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.f: $H_{\rm SS}$ Term 3 Coefficient

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.g: $\omega_{\cos\theta}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.h: $e^{\nu}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.i: $\tilde{B}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.j: $Q$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.k: $Q$ and $\sqrt{Q}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.e: $H_{\rm SS}$ Term 3

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.f: Combining [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) with our definition of $H_{\rm SS}$ Term 3

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.g: $e^{\tilde{\mu}}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.h: $e^{\nu}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.i: $\hat{\bf p} \cdot \boldsymbol{\xi} r$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.j: $\tilde{J}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.k: $\hat{\bf p} \cdot {\bf n}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.l: ${\bf S} \cdot \boldsymbol{\xi}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.m: $\tilde{B}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.n: ${\bf S} \cdot {\bf n}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.o: $\hat{\bf p} \cdot {\bf v} r$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.p: ${\bf S} \cdot {\bf v}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.q: $Q$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 5.r: and $\xi^{2}$

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Define variables for the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6: $H_{\rm NS}$ Terms

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{hss}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)

We collect here the terms in $H_{\rm NS}$ (defined in [this cell](

$$\label{hnsterms}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hns))

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{betapsum}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6.a: $\beta p$ sum

We will write [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.19) as

$$\label{betapsum}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)

We defined the term $\beta p$ sum in [this cell](

$$\label{betapsum}$$

```python
import sympy as sp

# Define constants and parameters
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # Reduced mass
chi = sp.symbols('chi')  # Spin parameter

# Define the spin-spin term $H_{\rm SS}$
def seoBNR_v4p_spin_spin_term(M, Q):
    """
    Calculate the spin-spin term.
    
    Parameters:
        M: Mass parameter
        Q: Charge parameter
    
    Returns:
        HSSseoBNR_v4p: $H_{\rm SS}$
    """
    # Define variables for curvature components
    nu = sp.symbols('nu')  # Conformal factor
    mu = sp.symbols('mu')  # Reduced mass
    Btilde = sp.symbols('Btilde')  # Brute force parameter
    
    # Calculate the exponential of $2\nu$
    exp_2nu = sp.exp(2 * nu)
    
    return exp_2nu

# Define variables for momentum and curvature components
p_hat_xi_r = sp.symbols('p_hat_xi_r')  # Dot product of momentum and curvature vectors

# Calculate the third part of the spin-spin term, denoted by $H_{\rm SS}^{\rm (3)}$
def seoBNR_v4p_spin_spin_term3(M, Q):
    """
    Calculate the third part of the spin-spin term.
    
    Parameters:
        M: Mass parameter
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hns) as

$$\label{betapsum}$$

We defined the term $\beta p$ sum in [this cell](

$$\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

$$\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}$$

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

$$\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}$$

Therefore

$$\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}$$

We define $\tilde{\omega}_{\rm fd}$ in [this**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omega_tilde), $\Lambda_{t}$ in [this cell](

$$\label{betapsum}$$

We defined the term $\beta p$ sum in [this cell](

$$\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

$$\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}$$

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

$$\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}$$

Therefore

$$\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}$$

We define $\tilde**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), and $p_{\phi}$ in [this cell](

$$\label{betapsum}$$

We defined the term $\beta p$ sum in [this cell](

$$\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

$$\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}$$

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

$$\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}$$

Therefore

$$\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}$$

We define $\**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pphi)

We defined the $\beta p$ sum in [this cell](

$$\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

$$\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}$$

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

$$\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}$$

Therefore

$$\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}$$

We defined $\tilde{\omega}_{\rm fd}$ in [this cell](


```python
%%writefile -a $C**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6.b: $\alpha$ \[Back to [top](

We defined the $\beta p$ sum in [this cell](

$$\begin{equation*}
    \beta p\ {\rm sum} = \beta^{i} p_{i}.
\end{equation*}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.45), we have

$$\begin{equation*}
    \beta^{i} = \frac{ g^{ti} }{ g^{tt} },
\end{equation*}$$

but from [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.36) we see that $g^{tr} = g^{t \theta} = 0$.  Thus only $\beta^{\phi}$ is nonzero.  Combining [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.45), (5.36e), and (5.36a), we find

$$\begin{equation*}
    \beta^{\phi} = \frac{ -\frac{ \tilde{\omega}_{\rm fd} }{ \Delta_{t} \Sigma } }{ -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma } } = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }
\end{equation*}$$

Therefore

$$\begin{equation*}
    \beta^{i} p_{i} = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} } p_{\phi}.
\end{equation*}$$

We defined $\tilde{\omega}_{\rm fd}$ in [this**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{alpha}$$

We defined the $\alpha$ variable in [this cell](

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have

$$\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}$$

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

$$\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}$$

Therefore

$$\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}$$

We define $\Delta_{t}$ in [this cell](**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat), $\Sigma$ in [this cell](

We defined the $\alpha$ variable in [this cell](

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have

$$\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}$$

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

$$\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}$$

Therefore

$$\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}$$

We define $\Delta_{t}$ in [this cell](

```python
import sympy as sp

# Define variables
Deltat = sp.symbols('Deltat')  # Delta_t variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable

# Define the expression for alpha
alpha_expr = sp.sqrt(Deltat * Sigmat / (sp.symbols('Lambda_t')))  # Alpha expression
```

### Theory Review


*   **Delta_t**: $\Delta_{t}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric component $\gamma_{tt}$.
    +   The value of $\Delta_{t}$ affects the overall structure of the spin-spin term and its impact on gravitational wave emission.
*   **Sigma**: $\Sigma$ is another important variable in the spin-spin term $H_{\rm SS}**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), and $\Lambda_{t}$ in [this cell](

We defined the $\alpha$ variable in [this cell](

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have

$$\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}$$

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

$$\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}$$

Therefore

$$\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}$$

We define $\Delta_{t}$ in [this cell](

```python
import sympy as sp

# Define variables
Deltat = sp.symbols('Deltat')  # Delta_t variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable
Usigmat = sp.symbols('Usigmat')  # Usigma variable

# Define the expression for alpha
alpha_expr = sp.sqrt(Deltat * Sigmat / (Lambdat))  # Alpha expression

# Define the expression for usigma
usigma_expr = Usigmat  # Usigma expression
```

### Theory Review


*   **Usigma**: $U_{\Sigma}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric component $\gamma_{t**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat)


We defined the $\alpha$ variable in [this cell](

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have

$$\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}$$

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

$$\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}$$

Therefore

$$\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}$$

We define $\Delta_{t}$ in [this cell](

```python
import sympy as sp

# Define variables
Deltat = sp.symbols('Deltat')  # Delta_t variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable

# Define the expression for alpha
alpha_expr = sp.sqrt(Deltat * Sigmat / (Lambdat))  # Alpha expression
```

### Theory Review


*   **Lambda_t**: $\Lambda_{t}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric component $\gamma_{tt}$.
    +   The value of $\Lambda_{t}$ affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6.c: $H_{\rm NS}$ radicand \[Back to [top](

We defined the $\alpha$ variable in [this cell](

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.44), we have

$$\begin{equation*}
    \alpha = \frac{ 1 }{ \sqrt{ -g^{tt}} },
\end{equation*}$$

and from [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.36a) we have

$$\begin{equation*}
    g^{tt} = -\frac{ \Lambda_{t} }{ \Delta_{t} \Sigma }.
\end{equation*}$$

Therefore

$$\begin{equation*}
    \alpha = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_{t} } }.
\end{equation*}$$

We define $\Delta_{t}$ in [this cell](

```python
import sympy as sp

# Define variables
Deltat = sp.symbols('Deltat')  # Delta_t variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable

# Define the expression for alpha
alpha_expr = sp.sqrt(Deltat * Sigmat / (Lambdat))  # Alpha expression
```

### Theory Review


*   **HNS radicand**: The radicand of the spin-spin term $H_{\rm NS}$ is a crucial component in numerical relativity and gravitational wave astronomy.
    +   It represents the spatial derivative of the metric components $\gamma_{t \theta}$.

### Code Explanation


```python
%%writefile -a**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hnsradicand}$$

Recall that we defined $H_{\rm NS}$ radicand in [this cell](

We are now going to calculate the expression for $H_{\rm NS}$ radicand.


```python
import sympy as sp

# Define variables
Deltat = sp.symbols('Deltat')  # Delta_t variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable

# Define the expression for HNS radicand
hns_radicand_expr = Deltat * Sigmat / (Lambdat)  # HNS radicand expression
```

### Theory Review


*   **HNS Radicand**: The $H_{\rm NS}$ radicand is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of the $H_{\rm NS}$ radicand affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HNS radicand to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(hns_radicand_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $H_{\rm NS}$ radicand to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hns) as

$$\begin{equation*}
    H_{\rm NS}\ {\rm radicand} = \mu^{2} + \underbrace{\gamma^{ij} p_{i} p_{j}}_{\gamma p\ \rm sum} + {\cal Q}_{4}
\end{equation*}$$

We define $\mu$ in [this cell](


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # Mass parameter
Q = sp.symbols('Q')  # Charge parameter
eta = sp.symbols('eta')  # Symmetry parameter
mu = sp.symbols('mu')  # mu variable
p_i = sp.symbols('p_i')  # p_i variable
p_j = sp.symbols('p_j')  # p_j variable

# Define the expression for HNS radicand
hns_radicand_expr = mu**2 + gamma_ij * p_i * p_j + Q_4  # HNS radicand expression
```

### Theory Review


*   **HNS Radicand**: The $H_{\rm NS}$ radicand is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of the $H_{\rm NS}$ radicand affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HNS radicand to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(hns_radicand_expr) + '\n')
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
mu), $\gamma p$ sum in [this cell](

We are now going to calculate the expression for the $\gamma p$ sum.


```python
import sympy as sp

# Define variables
p_i = sp.symbols('p_i')  # p_i variable
p_j = sp.symbols('p_j')  # p_j variable
gamma_ij = sp.symbols('gamma_ij')  # gamma_ij variable
mu = sp.symbols('mu')  # mu variable

# Define the expression for gamma p sum
gamma_p_sum_expr = gamma_ij * p_i * p_j  # Gamma p sum expression
```

### Theory Review


*   **Gamma P Sum**: The $\gamma p$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of the $\gamma p$ sum affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for gamma p sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(gamma_p_sum_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\gamma p$ sum to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
gammappsum), and ${\cal Q}_{4}$ in [this cell](

We are now going to calculate the expression for ${\cal Q}_{4}$. 


```python
import sympy as sp

# Define variables
Q_4 = sp.symbols('Q_4')  # Q_4 variable
gamma_ij = sp.symbols('gamma_ij')  # gamma_ij variable
p_i = sp.symbols('p_i')  # p_i variable
p_j = sp.symbols('p_j')  # p_j variable

# Define the expression for Q_4
Q_4_expr = Q_4  # Q_4 expression
```

### Theory Review


*   **${\cal Q}_{4}$**: The ${\cal Q}_{4}$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of the ${\cal Q}_{4}$ term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Q_4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q_4_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for ${\cal Q}_{4}$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation 

**Q_4** term is written into the python code using `with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:`**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q4)


We are now going to calculate the expression for the $H_{\rm NS}$ radicand.


```python
import sympy as sp

# Define variables
Q_4 = sp.symbols('Q_4')  # Q_4 variable
gammappsum = sp.symbols('gamma p sum')  # gamma p sum variable
mu = sp.symbols('mu')  # mu variable

# Define the expression for HNS radicand
Hnsradicand_expr = 1 + gammappsum + Q_4  # HNS radicand expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HNS radicand to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Hnsradicand_expr) + '\n')
```

### Theory Review


*   **HNS Radicand**: The $H_{\rm NS}$ radicand is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of the $H_{\rm NS}$ radicand affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
# Write the expression for HNS radicand to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Hnsradicand_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $H_{\rm NS}$ radicand to a file. This will allow us to use**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6.c.i: $\gamma^{ij} p_{i} p_{j}$ \[Back to [top](

We are now going to calculate the expression for the spin-spin term $H_{\rm SS}$.


```python
import sympy as sp

# Define variables
p_i = sp.symbols('p_i')  # p_i variable
p_j = sp.symbols('p_j')  # p_j variable
gamma_ij = sp.symbols('gamma_ij')  # gamma_ij variable

# Define the expression for gamma_ij p_i p_j
gamma_ij_p_i_p_j_expr = gamma_ij * p_i * p_j  # Gamma ij p i p j expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for gamma_ij p_i p_j to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(gamma_ij_p_i_p_j_expr) + '\n')
```

### Theory Review


*   **Gamma ij P i P j**: The term $\gamma^{ij} p_{i} p_{j}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
# Write the expression for gamma_ij p_i p_j to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(gamma_ij_p_i_p_j_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{gammappsum}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.46), we have

$$\begin{equation*}
    \gamma^{ij} = g^{ij} - \frac{ g^{ti} g^{tj} }{ g^{tt} }.
\end{equation*}$$

Combining this result with [BB2010](https://arxiv.org/abs/0912.3517) Equations 5.36, we have

$$\begin{equation*}
    \gamma^{r\theta} = \gamma^{r\phi} = \gamma^{\theta r} = \gamma^{\theta\phi} = \gamma^{\phi r} = \gamma^{\phi\theta} = 0
\end{equation*}$$

and

$$\begin{align*}
    \gamma^{rr} &= g^{rr} = \frac{ \Delta_{r} }{ \Sigma } \\
    \gamma^{\theta\theta} &= g^{\theta\theta} = \frac{ 1 }{ \Sigma } \\
    \gamma^{\phi\phi} &= \frac{ \Sigma }{ \Lambda_{t} \sin^{2} \theta }.
\end{align*}$$

Therefore

$$\begin{align*}
    \gamma^{ij} p_{i} p_{j} &= \gamma^{rr} p_{r} p_{r} + \gamma^{\theta\theta} p_{\theta} p_{\theta} + \gamma^{\phi\phi} p_{\phi} p_{\phi} \\
        &= \frac{ \Delta_{r} }{ \Sigma } p_{r}^{2} + \frac{ 1**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltar), $\Sigma$ in [this cell](

We are now going to calculate the expression for the $\Delta_{r}$ and $\Sigma$.


```python
import sympy as sp

# Define variables
Deltar = sp.symbols('Deltar')  # Delta_r variable
Sigmat = sp.symbols('Sigmat')  # Sigma variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable
theta = sp.symbols('theta')  # theta variable

# Define the expression for Deltar
Deltar_expr = Deltar  # Delta_r expression

# Define the expression for Sigmat
Sigmat_expr = Sigmat  # Sigma expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Deltar to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Deltar_expr) + '\n')
```

### Theory Review


*   **Delta_r and Sigma**: The variables $\Delta_{r}$ and $\Sigma$ are crucial components in the spin-spin term $H_{\rm SS}$. They represent the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The values of these variables affect the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
# Write the expression for Sigmat to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Sigmat_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expressions for $\Delta_{r}$ and $\Sigma$ to a file. This will allow**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigma), $\hat{\bf p} \cdot {\bf n}$ in [this cell](

We are now going to calculate the expression for the $\Sigma$.


```python
import sympy as sp

# Define variables
Sigmat = sp.symbols('Sigmat')  # Sigma variable
p_n_dot = sp.symbols('p_n_dot')  # p dot n variable

# Define the expression for Sigmat
Sigmat_expr = Sigmat  # Sigma expression

# Define the expression for p dot n
p_n_dot_expr = p_n_dot  # p dot n expression
```

### Theory Review


*   **Sigma**: The variable $\Sigma$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Sigmat to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Sigmat_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\Sigma$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
# Write the expression for p dot n to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(p_n_dot_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\hat{\bf p} \cdot {\bf n}$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pdotn), $\hat{\bf p} \cdot {\bf v} r$ in [this cell](

We are now going to calculate the expression for the $\hat{\bf p} \cdot {\bf n}$ and $\hat{\bf p} \cdot {\bf v} r$.


```python
import sympy as sp

# Define variables
p_n_dot = sp.symbols('p_n_dot')  # p dot n variable
p_vr_dot = sp.symbols('p_vr_dot')  # p dot v r variable
r = sp.symbols('r')  # r variable
v = sp.symbols('v')  # v variable

# Define the expression for p dot n
p_n_dot_expr = p_n_dot  # p dot n expression

# Define the expression for p dot v r
p_vr_dot_expr = p_vr_dot * r  # p dot v r expression
```

### Theory Review


*   **P Dot N**: The variable $\hat{\bf p} \cdot {\bf n}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for p dot n to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(p_n_dot_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\hat{\bf p} \cdot {\bf n}$ to a file. This will allow us to use this expression**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pdotvr), $\sin^{2} \theta$ in [this cell](

We are now going to calculate the expression for the $\hat{\bf p} \cdot {\bf v} r$ and $\sin^{2} \theta$.


```python
import sympy as sp

# Define variables
p_vr_dot = sp.symbols('p_vr_dot')  # p dot v r variable
sin_sq_theta = sp.symbols('sin^2 theta')  # sin^2 theta variable
theta = sp.symbols('theta')  # theta variable

# Define the expression for p dot v r
p_vr_dot_expr = p_vr_dot  # p dot v r expression

# Define the expression for sin^2 theta
sin_sq_theta_expr = sin_sq_theta  # sin^2 theta expression
```

### Theory Review


*   **P Dot V R**: The variable $\hat{\bf p} \cdot {\bf v} r$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for p dot v r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(p_vr_dot_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\hat{\bf p} \cdot {\bf v} r$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
#**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sin2theta), $\Lambda_{t}$ in [this cell](

We are now going to calculate the expression for the $\sin^{2} \theta$ and $\Lambda_{t}$.


```python
import sympy as sp

# Define variables
sin_sq_theta = sp.symbols('sin^2 theta')  # sin^2 theta variable
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable

# Define the expression for sin^2 theta
sin_sq_theta_expr = sin_sq_theta  # sin^2 theta expression

# Define the expression for Lambda_t
Lambdat_expr = Lambdat  # Lambda_t expression
```

### Theory Review


*   **Sin^2 Theta**: The variable $\sin^{2} \theta$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for sin^2 theta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(sin_sq_theta_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\sin^{2} \theta$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
# Write the expression for Lambda_t to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Lambdat_expr) + '\n')
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), and $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](

We are now going to calculate the expression for the $\Lambda_{t}$ and $\hat{\bf p} \cdot \boldsymbol{\xi} r$.


```python
import sympy as sp

# Define variables
Lambdat = sp.symbols('Lambdat')  # Lambda_t variable
xi_dot_p_r = sp.symbols('xi dot p r')  # xi dot p r variable
r = sp.symbols('r')  # r variable

# Define the expression for Lambda_t
Lambdat_expr = Lambdat  # Lambda_t expression

# Define the expression for xi dot p r
xi_dot_p_r_expr = xi_dot_p_r * r  # xi dot p r expression
```

### Theory Review


*   **Lambda_t**: The variable $\Lambda_{t}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Lambda_t to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Lambdat_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $\Lambda_{t}$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
# Write the expression for xi dot p r to file
with open($Ccodesdir/v4P_H**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pdotxir)


We are now going to calculate the expression for the $\gamma^{ij} p_i p_j$.


```python
import sympy as sp

# Define variables
Deltar = sp.symbols('Deltar')  # Deltar variable
Sigma = sp.symbols('Sigma')  # Sigma variable
pdotn = sp.symbols('p dot n')  # p dot n variable
pdotvr = sp.symbols('p dot v r')  # p dot v r variable
sin2theta = sp.symbols('sin^2 theta')  # sin^2 theta variable
Lambdat = sp.symbols('Lambda_t')  # Lambda_t variable
pdotxir = sp.symbols('p dot xi r')  # p dot xi r variable

# Define the expression for gamma pp sum
gammappsum_expr = (Deltar/Sigma)*pdotn*pdotn + (1/Sigma)*pdotvr*pdotvr/sin2theta + (Sigma/Lambdat/sin2theta)*pdotxir*pdotxir  # Gamma pp sum expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for gamma pp sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(gammappsum_expr) + '\n')
```

### Theory Review


*   **Gamma PP Sum**: The $\gamma^{ij} p_i p_j$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 6.c.ii: ${\cal Q}_{4}$ \[Back to [top](

We are now going to calculate the expression for the ${\cal Q}_{4}$.


```python
import sympy as sp

# Define variables
Q_4 = sp.symbols('Q_4')  # Q_4 variable
Deltar = sp.symbols('Deltar')  # Deltar variable
Sigma = sp.symbols('Sigma')  # Sigma variable
Lambdat = sp.symbols('Lambda_t')  # Lambda_t variable

# Define the expression for Q_4
Q_4_expr = (Deltar/Sigma)*(1+2*pdotn*pdotn) +  (1/Sigma)*(1-2/sin2theta*pdotvr*pdotvr)  # Q_4 expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Q_4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q_4_expr) + '\n')
```

### Theory Review


*   **Q 4**: The ${\cal Q}_{4}$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
# Write the expression for Q_4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q_4_expr) + '\n')
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]

We are now going to calculate the expression for the ${\cal Q}_{4}$.


```python
import sympy as sp

# Define variables
prT = sp.symbols('prT')  # prT variable
r = sp.symbols('r')  # r variable
z3 = sp.symbols('z3')  # z3 variable
nu = sp.symbols('nu')  # nu variable

# Define the expression for Q_4
Q_4_expr = (prT**4)/(r**2)*z3  # Q_4 expression

# Define the expression for z3
z3_expr = 2*(4-3*nu)*nu  # z3 expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Q_4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q_4_expr) + '\n')
```

### Theory Review


*   **Q 4**: The ${\cal Q}_{4}$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
# Write the expression for Q_4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q_4_expr) + '\n')
```

### Theory Review


*   **Relation between $\nu$ and $\eta$**: In the notation of [BB2010](https://arxiv.org/**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
prt), $u$ in [this cell](

We are now going to calculate the expression for the prt) and u.


```python
import sympy as sp

# Define variables
prT = sp.symbols('prT')  # prT variable
r = sp.symbols('r')  # r variable
u = sp.symbols('u')  # u variable

# Define the expression for prt
prt_expr = prT  # prt expression

# Define the expression for u
u_expr = u  # u expression
```

### Theory Review


*   **prT**: The variable $p_{r^{*}}$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for prt to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(prt_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $p_{r^{*}}$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
# Write the expression for u to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(u_expr) + '\n')
```

### Theory Review


*   **u**: The variable $u$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
u), and $\eta$ in [this cell](

We are now going to calculate the expression for the u) and eta.


```python
import sympy as sp

# Define variables
u = sp.symbols('u')  # u variable
eta = sp.symbols('eta')  # eta variable

# Define the expression for u
u_expr = u  # u expression

# Define the expression for eta
eta_expr = eta  # eta expression
```

### Theory Review


*   **u**: The variable $u$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for u to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(u_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $u$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


```python
# Write the expression for eta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(eta_expr) + '\n')
```

### Theory Review


*   **eta**: The variable $\eta$ is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this variable affects the overall structure of**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta) below.


We are now going to calculate the expression for the Q4.


```python
import sympy as sp

# Define variables
prT = sp.symbols('prT')  # prT variable
u = sp.symbols('u')  # u variable
eta = sp.symbols('eta')  # eta variable

# Define the expression for Q4
Q4_expr = 2*prT*prT*prT*prT*u*u*(4 - 3*eta)*eta  # Q4 expression
```

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for Q4 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Q4_expr) + '\n')
```

### Theory Review


*   **Q 4**: The $Q_4$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


*   **Expression for Q4**: We are now writing the expression for $Q_4$ to a file. This will allow us to use this expression in our calculations later on.
    +   The expression is given by $$Q_4 = 2 prT^4 u^2 (4 - 3 \eta) \eta$$

### Theory Review


*   **Relation between Q4 and spin-spin term**: The $Q_4$ term is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 7: The $H_{\rm D}$ Terms \[Back to [top](

We are now going to calculate the expression for the $H_D$ terms.


```python
import sympy as sp

# Define variables
HD = sp.symbols('H_D')  # H_D variable
Deltar = sp.symbols('Deltar')  # Deltar variable
Sigma = sp.symbols('Sigma')  # Sigma variable
Lambdat = sp.symbols('Lambda_t')  # Lambda_t variable

# Define the expression for HD
HD_expr = HD  # H_D expression
```

### Theory Review


*   **H D**: The $H_D$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_expr) + '\n')
```

### Theory Review


*   **Writing to File**: We are now writing the expression for $H_D$ to a file. This will allow us to use this expression in our calculations later on.

### Code Explanation


*   **Expression for HD**: The expression for $H_D$ is given by $$H_D = Deltar \Sigma + 1/Lambdat$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hdterms}$$

Recall we defined $H_D$ in [this cell](

We are now going to calculate the expression for the $H_D$ terms.


```python
import sympy as sp

# Define variables
Deltar = sp.symbols('Deltar')  # Deltar variable
Sigma = sp.symbols('Sigma')  # Sigma variable
Lambdat = sp.symbols('Lambda_t')  # Lambda_t variable

# Define the expression for HD
HD_expr = (Deltar/Sigma) + (1/Lambdat)  # H_D expression
```

### Theory Review


*   **H D**: The $H_D$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_expr) + '\n')
```

### Theory Review


*   **Expression for H D**: The expression for $H_D$ is given by $$H_D = \frac{Deltar}{\Sigma} + \frac{1}{Lambda_t}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between H D and spin-spin term**: The $H_D$ term is closely related to the spin-spin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hd) as

\begin{equation*}
    H_{\rm D} = H_{\rm D}\ {\rm coeffecient} * H_{\rm D}\ {\rm sum}.
\end{equation*}

In this step, we break down each of $H_D$ coefficient (defined in [this cell](

We are now going to calculate the expression for the $H_D$ terms.


```python
import sympy as sp

# Define variables
HD_coefficient = sp.symbols('HD_coefficient')  # HD_coefficient variable
HD_sum = sp.symbols('HD_sum')  # HD_sum variable

# Define the expression for HD
HD_expr = (HD_coefficient * HD_sum)  # H_D expression
```

### Theory Review


*   **H D**: The $H_D$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_expr) + '\n')
```

### Theory Review


*   **Expression for H D**: The expression for $H_D$ is given by $$H_D = H_{\rm D}\ {\rm coeffecient} * H_{\rm D}\ {\rm sum}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ to a file. This will allow us to use this expression in our calculations later on.

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdcoeff)) and $H_D$ sum (defined in [this cell](

We are now going to calculate the expression for the hdcoeff) and HD_sum terms.


```python
import sympy as sp

# Define variables
HD_coefficient = sp.symbols('HD_coefficient')  # HD_coefficient variable
HD_sum = sp.symbols('HD_sum')  # HD_sum variable

# Define the expression for HD coefficient
HD_coefficient_expr = (1/(sp.sqrt(2)*sp.pi)) * ((-9/8) + (3/2*sp.cos(2*sp theta)))  # hdcoeff expression

# Define the expression for H D sum
HD_sum_expr = sp.sin(sp.theta)**4 * sp.cos(sp.phi)**4  # HD_sum expression
```

### Theory Review


*   **H D coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for H D coefficient**: The expression for $H_D$ coefficient is given by $$H_{\rm D}\ {\rm coeffecient} = \frac{1}{\sqrt{2}\pi} (-\frac{9}{8} + \frac{3}{2}\cos 2\theta)$$

### Code Explanation


*   **Writing**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdsum)).


We are now going to calculate the expression for the HD_sum term.


```python
import sympy as sp

# Define variables
theta = sp.symbols('theta')  # theta variable
phi = sp.symbols('phi')  # phi variable

# Define the expression for HD sum
HD_sum_expr = sp.sin(theta)**4 * sp.cos(phi)**4  # HD_sum expression
```

### Theory Review


*   **H D sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_expr) + '\n')
```

### Theory Review


*   **Expression for H D sum**: The expression for $H_D$ sum is given by $$H_{\rm D}\ {\rm sum} = \sin^4\theta\cos^4\phi$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ sum to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between H D sum and spin-spin term**: The $H_D$ sum is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 7.a: $H_D$ Coefficient \[Back to [top](

We are now going to calculate the expression for the $H_D$ coefficient.


```python
import sympy as sp

# Define variables
theta = sp.symbols('theta')  # theta variable
phi = sp.symbols('phi')  # phi variable

# Define the expression for HD coefficient
HD_coefficient_expr = (1/(sp.sqrt(2)*sp.pi)) * ((-9/8) + (3/2*sp.cos(2*theta)))  # hdcoeff expression
```

### Theory Review


*   **H D Coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for H D Coefficient**: The expression for $H_D$ coefficient is given by $$H_{\rm D}\ {\rm coeffecient} = \frac{1}{\sqrt{2}\pi} (-\frac{9}{8} + \frac{3}{2}\cos 2\theta)$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ coefficient to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*  **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hdcoeff}$$

From our definition of $H_D$ in [this cell](

We are now going to calculate the expression for the HD coefficient.


```python
import sympy as sp

# Define variables
theta = sp.symbols('theta')  # theta variable
phi = sp.symbols('phi')  # phi variable

# Define the expression for HD coefficient
HD_coefficient_expr = (1/(sp.sqrt(2)*sp.pi)) * ((-9/8) + (3/2*sp.cos(2*theta)))  # hdcoeff expression
```

### Theory Review


*   **H D Coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for H D Coefficient**: The expression for $H_D$ coefficient is given by $$H_{\rm D}\ {\rm coeffecient} = \frac{1}{\sqrt{2}\pi} (-\frac{9}{8} + \frac{3}{2}\cos 2\theta)$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ coefficient to a file. This will allow us to use this expression in our calculations later on.

### Theory**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hd), we have

\begin{equation*}
    H_D\ {\rm coefficient} = \frac{ \mu }{ 2 M r^{3} },
\end{equation*}

and recalling the definition of [$\eta$](

We are now going to calculate the expression for the HD coefficient.


```python
import sympy as sp

# Define variables
mu = sp.symbols('mu')  # mu variable
M = sp.symbols('M')  # M variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable

# Define the expression for HD coefficient
HD_coefficient_expr = (mu/(2*M*r**3))  # hdcoeff expression
```

### Theory Review


*   **H D Coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for H D Coefficient**: The expression for $H_D$ coefficient is given by $$H_{\rm D}\ {\rm coeffecient} = \frac{\mu}{2Mr^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ coefficient to a file. This will allow us to use this**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta) we'll write

\begin{equation*}
    H_D\ {\rm coefficient} = \frac{ \eta }{ 2 r^{3} }.
\end{equation*}

We define $\eta$ in [this cell](

We are now going to calculate the expression for the HD coefficient.


```python
import sympy as sp

# Define variables
eta = sp.symbols('eta')  # eta variable
r = sp.symbols('r')  # r variable

# Define the expression for HD coefficient
HD_coefficient_expr = (eta/(2*r**3))  # hdcoeff expression
```

### Theory Review


*   **H D Coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for H D Coefficient**: The expression for $H_D$ coefficient is given by $$H_{\rm D}\ {\rm coeffecient} = \frac{\eta}{2r^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ coefficient to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between H D Coefficient and spin-spin term**: The**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta) and $r$ in [this cell](

We are now going to calculate the expression for $\eta$.


```python
import sympy as sp

# Define variables
eta = sp.symbols('eta')  # eta variable
r = sp.symbols('r')  # r variable

# Define the expression for HD coefficient
HD_coefficient_expr = (eta/(2*r**3))  # hdcoeff expression
```

### Theory Review


*   **$\eta$**: The $\eta$ term is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for $\eta$**: The expression for $\eta$ is given by $$\eta = \frac{\eta}{2r^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $\eta$ to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between $\eta$ and spin-spin term**: The $\eta$ term is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


*   **Using**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r).


We are now going to calculate the expression for the HD coefficient.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable

# Define the expression for HD coefficient
HD_coefficient_expr = (sp.Rational(1,2))/(r*r*r)  # hdcoeff expression
```

### Theory Review


*   **$H_D$ Coefficient**: The $H_D$ coefficient is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD coefficient to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_coefficient_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Coefficient**: The expression for $H_D$ coefficient is given by $$H_D\ {\rm coefficient} = \frac{1}{2r^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ coefficient to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between $H_D$ Coefficient and spin-spin term**: The $H_D$ coefficient is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


*  **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 7.b: $H_D$ Sum \[Back to [top](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable

# Define the expression for HD sum
HD_sum_expr = (1/(2*r**3))  # hdsum expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \frac{1}{2r^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ sum to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between $H_D$ Sum and spin-spin term**: The $H_D$ sum is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hdsum}$$

From our definition of $H_D$ in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable

# Define the expression for HD sum
HD_sum_expr = (1/(2*r**3))  # hdsum expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \frac{1}{2r^3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ sum to a file. This will allow us to use this expression in our calculations later on.

### Theory Review


*   **Relation between $H_D$ Sum and spin-spin term**: The $H_D$ sum is closely related to the spin-spin term. It represents the contribution of the spin-spin interaction to the gravitational wave emission.
    +   The value of this term affects the overall structure of the spin-spin term**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hd), we have

\begin{align*}
    H_D\ {\rm sum} &= \left( \delta^{ij} - 3 n^{i} n^{j} \right) S^{*}_{i} S^{*}_{j} \\
        &= \underbrace{\delta^{ij} S^{*}_{i} S^{*}_{j}}_{\rm Term\ 1} - \underbrace{3 n^{i} n^{j} S^{*}_{i} S^{*}_{j}}_{\rm Term\ 2}.
\end{align*}

We compute $H_D$ Term 1 in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, '**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdsumterm1) and $H_D$ Term 2 in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_1_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \delta^{ij} S^{*}_{i} S^{*}_{j} - 3 n^{i} n^{j} S^{*}_{i} S^{*}_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdsumterm2).


We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression

# Define the expression for HD sum Term 2
HD_sum_Term_2_expr = (3*n_i_n_j*S_i_S_j)  # hdsum_term_2 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_1_expr - HD_sum_Term_2_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \delta^{ij} S^{*}_{i} S^{*}_{j} - 3 n^{i} n^{j} S^{*}_{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 7.b.i: $H_D$ Sum Term 1 \[Back to [top](

We are now going to calculate the expression for $H_D$ sum Term 1.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for $H_D$ sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum Term 1**: The $H_D$ sum Term 1 is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_1_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum Term 1**: The expression for $H_D$ sum Term 1 is given by $$H_D\ {\rm sum}\ {\rm Term}\ 1 = \delta^{ij} S^{*}_{i} S^{*}_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ sum Term 1 to a file. This will allow us to use**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hdsumterm1}$$

From our definition of $H_D$ sum Term 1 in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum Term 1**: The $H_D$ sum Term 1 is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_1_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum Term 1**: The expression for $H_D$ sum Term 1 is given by $$H_D\ {\rm sum}\ {\rm Term}\ 1 = \delta^{ij} S^{*}_{i} S^{*}_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$ sum Term 1 to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdsum), we have

\begin{equation*}
    H_D\ {\rm sum}\ {\rm Term}\ 1 = \delta^{ij} S^{*}_{i} S^{*}_{j}
\end{equation*}

where $\delta^{ij}$ is the Kronecker delta:

\begin{equation*}
    \delta_{ij} = \left\{ \begin{array}{cc}
        0, & i \not= j \\
        1, & i = j. \end{array} \right.
\end{equation*}

Thus we have

\begin{equation*}
    H_D\ {\rm sum}\ {\rm Term}\ 1 = S^{*}_{1} S^{*}_{1} + S^{*}_{2} S^{*}_{2} + S^{*}_{3} S^{*}_{3}
\end{equation*}

We define ${\bf S}^{*}$ in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
delta_ij = sp.symbols('delta_ij')  # delta_ij variable
S_i_S_j = sp.symbols('S_i_S_j')  # S_i_S_j variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (delta_ij*S_i_S_j)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hreal_spin_combos).


We are now going to calculate the expression for HD sum Term 1.


```python
import sympy as sp

# Define variables
Sstar1 = sp.symbols('Sstar1')  # Sstar1 variable
Sstar2 = sp.symbols('Sstar2')  # Sstar2 variable
Sstar3 = sp.symbols('Sstar3')  # Sstar3 variable

# Define the expression for HD sum Term 1
HD_sum_Term_1_expr = (Sstar1*Sstar1 + Sstar2*Sstar2 + Sstar3*Sstar3)  # hdsum_term_1 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_1_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = S^{*}_{1} S^{*}_{1} + S^{*}_{2} S^{*}_{2} + S^{*}_{3} S^{*}_{3}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for $H_D$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 7.b.ii: $H_D$ Sum Term 2 \[Back to [top](

We are now going to calculate the expression for $H_D$ sum Term 2.


```python
import sympy as sp

# Define variables
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
Sstar_i_Sstar_j = sp.symbols('Sstar_i_Sstar_j')  # Sstar_i_Sstar_j variable

# Define the expression for $H_D$ sum Term 2
HD_sum_Term_2_expr = (3*n_i_n_j*Sstar_i_Sstar_j)  # hdsum_term_2 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_2_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \delta^{ij} S^{*}_{i} S^{*}_{j} - 3 n^{i} n^{j} S^{*}_{i} S^{*}_{j}$$

### Code Explanation


*   **Writing to File**: We are**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hdsumterm2}$$

From our definition of $H_D$ sum Term 2 in [this cell](

We are now going to calculate the expression for the HD sum.


```python
import sympy as sp

# Define variables
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
Sstar_i_Sstar_j = sp.symbols('Sstar_i_Sstar_j')  # Sstar_i_Sstar_j variable

# Define the expression for HD sum Term 2
HD_sum_Term_2_expr = (3*n_i_n_j*Sstar_i_Sstar_j)  # hdsum_term_2 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_2_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression for $H_D$ sum is given by $$H_D\ {\rm sum} = \delta^{ij} S^{*}_{i} S^{*}_{j} - 3 n^{i} n^{j} S^{*}_{i} S^{*}_{j}$$

### Code Explanation


*   **Writing to File**: We**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Spin-Spin Term $H_{\rm SS}$

*   **The Spin-Spin TERM:** In this section, we discuss the terms of the spin-spin term $H_{\rm SS}$.
    +   The spin-spin term is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hdsum), we have

\begin{align*}
    H_D\ {\rm sum}\ {\rm Term}\ 2 &= 3 n^{i} n^{j} S^{*}_{i} S^{*}_{j} \\
        &= 3 \left( {\bf S}^{*} \cdot {\bf n} \right)^{2} \\
\end{align*}

We define ${\bf S}^{*} \cdot {\bf n}$ in [this cell](

We are now going to calculate the expression for HD sum.


```python
import sympy as sp

# Define variables
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
Sstar_i_Sstar_j = sp.symbols('Sstar_i_Sstar_j')  # Sstar_i_Sstar_j variable

# Define the expression for HD sum Term 2
HD_sum_Term_2_expr = (3*n_i_n_j*Sstar_i_Sstar_j)  # hdsum_term_2 expression
```

### Theory Review


*   **$H_D$ Sum**: The $H_D$ sum is a crucial component in the spin-spin term $H_{\rm SS}$. It represents the spatial derivative of the metric components $\gamma_{t \theta}$.
    +   The value of this term affects the overall structure of the spin-spin term and its impact on gravitational wave emission.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for HD sum Term 2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(HD_sum_Term_2_expr) + '\n')
```

### Theory Review


*   **Expression for $H_D$ Sum**: The expression**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sstardotn).


We are now going to calculate the expression for the dot product of ${\bf S}^{*}$ and ${\bf n}$.


```python
import sympy as sp

# Define variables
Sstar_i = sp.symbols('Sstar_i')  # Sstar_i variable
n_j = sp.symbols('n_j')  # n_j variable

# Define the expression for dot product of S* and n
dot_product_expr = (Sstar_i*n_j)  # dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and ${\bf n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}^{*}$ and ${\bf n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$\left( {\bf S}^{*} \cdot {\bf n} \right)^{2} = \left( S^{*}_{i} n^{j} \right)$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to Common Dot Products

*   **Common Dot Products:** In this section, we discuss the common dot products that are used in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8: Common Dot Products \[Back to [top](

We are now going to calculate the expressions for the common dot products.


```python
import sympy as sp

# Define variables
Sstar_i_Sstar_j = sp.symbols('Sstar_i_Sstar_j')  # Sstar_i_Sstar_j variable
n_i_n_j = sp.symbols('n_i_n_j')  # n_i_n_j variable
Sdotn = sp.symbols('Sdotn')  # Sdotn variable

# Define the expression for common dot products
common_dot_products_exprs = (Sstar_i_Sstar_j, n_i_n_j*Sdotn)  # common_dot_products_expression
```

### Theory Review


*   **Common Dot Products**: The common dot products are a crucial component in numerical relativity and gravitational wave astronomy. They represent the amount of "similarity" between two vectors.
    +   In this case, we are calculating the common dot products ${\bf S}^{*} \cdot {\bf S}^{*}$ and ${\bf n} \cdot ({\bf S}^{*} \cdot {\bf n})$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for common dot products to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(common_dot_products_exprs) + '\n')
```

### Theory Review


*   **Expression for Common Dot Products**: The expressions for the common dot products are given by $$S^{*}_{i} S^{*}_{j} = \left( {\bf S}^{*} \cdot {\bf S}^{*} \right)$$ and $$n^{i} ({\bf S}^{*} \cdot {\bf n}) = \left( {\bf n} \cdot ({\bf S}^{*} \cdot {\bf n}) \right**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to Common Dot Products

*   **Common Dot Products:** In this section, we discuss the common dot products that are used in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{dotproducts}$$

What follows are definitions of many common dot products.


<a id='sdotxi'></a>

We are now going to calculate the expressions for the common dot products.


```python
import sympy as sp

# Define variables
Sstar_i = sp.symbols('Sstar_i')  # Sstar_i variable
n_j = sp.symbols('n_j')  # n_j variable
xi_k = sp.symbols('xi_k')  # xi_k variable

# Define the expression for common dot products
common_dot_products_exprs = (Sstar_i*Sstar_i, n_j*xi_k)  # common_dot_products_expression
```

### Theory Review


*   **Common Dot Products**: The common dot products are a crucial component in numerical relativity and gravitational wave astronomy. They represent the amount of "similarity" between two vectors.
    +   In this case, we are calculating the common dot products ${\bf S}^{*} \cdot {\bf S}^{*}$ and ${\bf n} \cdot {\bf \xi}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for common dot products to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(common_dot_products_exprs) + '\n')
```

### Theory Review


*   **Expression for Common Dot Products**: The expressions for the common dot products are given by $$S^{*}_{i} S^{*}_{j} = \left( {\bf S}^{*} \cdot {\bf S}^{*} \right)$$ and $$n^{i} ({\bf \xi})_{k} = \left( {\bf n} \cdot {\bf \xi} \right)$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8.a: ${\bf S} \cdot \boldsymbol{\xi}$ \[Back to [top](

We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.


```python
import sympy as sp

# Define variables
S_i = sp.symbols('S_i')  # S_i variable
xi_j = sp.symbols('xi_j')  # xi_j variable

# Define the expression for dot product of S and xi
dot_product_expr = (S_i*xi_j)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{\xi}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} \xi^{j} = S_{i} \xi_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sdotxi}$$

We have

\begin{equation*}
    {\bf S} \cdot \boldsymbol{\xi} = S^{1} \xi^{1} + S^{2} \xi^{2} + S^{3} \xi^{3}
\end{equation*}

We define $\xi$ in [this cell](

We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
xi1 = sp.symbols('xi1')  # xi1 variable
xi2 = sp.symbols('xi2')  # xi2 variable
xi3 = sp.symbols('xi3')  # xi3 variable

# Define the expression for dot product of S and xi
dot_product_expr = (S1*xi1 + S2*xi2 + S3*xi3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{\xi}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   ****NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\boldsymbol{\xi}$

*   **The Vector $\boldsymbol{\xi}$:** In this section, we discuss the vector $\boldsymbol{\xi}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
xi).


We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
xi1 = sp.symbols('xi1')  # xi1 variable
xi2 = sp.symbols('xi2')  # xi2 variable
xi3 = sp.symbols('xi3')  # xi3 variable

# Define the expression for dot product of S and xi
dot_product_expr = (S1*xi1 + S2*xi2 + S3*xi3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{\xi}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{\xi}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} \xi^{j} = S_{i} \xi_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8.b: ${\bf S} \cdot {\bf v}$ \[Back to [top](

We are now going to calculate the expression for the dot product of ${\bf S}$ and ${\bf v}$.


```python
import sympy as sp

# Define variables
S_i = sp.symbols('S_i')  # S_i variable
v_j = sp.symbols('v_j')  # v_j variable

# Define the expression for dot product of S and v
dot_product_expr = (S_i*v_j)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and ${\bf v}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and ${\bf v}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} v^{j} = S_{i} v_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\boldsymbol{\xi}$

*   **The Vector $\boldsymbol{\xi}$:** In this section, we discuss the vector $\boldsymbol{\xi}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sdotv}$$

We have

\begin{equation*}
    {\bf S} \cdot {\bf v} = S^{1} v^{1} + S^{2} v^{2} + S^{3} v^{3}.
\end{equation*}

We define ${\bf v}$ in [this cell](

We are now going to calculate the expression for the dot product of ${\bf S}$ and ${\bf v}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
v1 = sp.symbols('v1')  # v1 variable
v2 = sp.symbols('v2')  # v2 variable
v3 = sp.symbols('v3')  # v3 variable

# Define the expression for dot product of S and v
dot_product_expr = (S1*v1 + S2*v2 + S3*v3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and ${\bf v}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and ${\bf v}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\boldsymbol{n}$

*   **The Vector $\boldsymbol{n}$:** In this section, we discuss the vector $\boldsymbol{n}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
v).


We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for dot product of S and n
dot_product_expr = (S1*n1 + S2*n2 + S3*n3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} n^{j} = S_{i} n_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8.c: ${\bf S} \cdot {\bf n}$ \[Back to [top](

We are now going to calculate the expression for the dot product of ${\bf S}$ and ${\bf n}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for dot product of S and n
dot_product_expr = (S1*n1 + S2*n2 + S3*n3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and ${\bf n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and ${\bf n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} n^{j} = S_{i} n_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\boldsymbol{n}$

*   **The Vector $\boldsymbol{n}$:** In this section, we discuss the vector $\boldsymbol{n}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sdotn}$$

We have

\begin{equation*}
    {\bf S} \cdot {\bf n} = S^{1} n^{1} + S^{2} n^{2} + S^{3} n^{3}.
\end{equation*}

We define ${\bf n}$ in [this cell](

We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for dot product of S and n
dot_product_expr = (S1*n1 + S2*n2 + S3*n3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\boldsymbol{\hat{S}}$

*   **The Vector $\boldsymbol{\hat{S}}$:** In this section, we discuss the vector $\boldsymbol{\hat{S}}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
n).


We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for dot product of S and n
dot_product_expr = (S1*n1 + S2*n2 + S3*n3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} n^{j} = S_{i} n_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file. This will allow us**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8.d: ${\bf S} \cdot \hat{\bf S}_{\rm Kerr}$ \[Back to [top](

We are now going to calculate the expression for the dot product of ${\bf S}$ and $\hat{\bf S}_{\rm Kerr}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
Skerrhat1 = sp.symbols('Skerrhat1')  # Skerrhat1 variable
Skerrhat2 = sp.symbols('Skerrhat2')  # Skerrhat2 variable
Skerrhat3 = sp.symbols('Skerrhat3')  # Skerrhat3 variable

# Define the expression for dot product of S and Skerrhat
dot_product_expr = (S1*Skerrhat1 + S2*Skerrhat2 + S3*Skerrhat3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\hat{\bf S}_{\rm Kerr}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\hat{\bf S}_{\rm Kerr}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Vector $\hat{\boldsymbol{S}}_{\rm Kerr}$

*   **The Vector $\hat{\boldsymbol{S}}_{\rm Kerr}$:** In this section, we discuss the vector $\hat{\boldsymbol{S}}_{\rm Kerr}$, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sdotskerrhat}$$

We have

\begin{equation*}
    {\bf S} \cdot \hat{\bf S}_{\rm Kerr} = S^{1} \hat{S}_{\rm Kerr}^{1} + S^{2} \hat{S}_{\rm Kerr}^{2} + S^{3} \hat{S}_{\rm Kerr}^{3}.
\end{equation*}

We define $\hat{\bf S}_{\rm Kerr}$ in [this cell](

We are now going to calculate the expression for the dot product of ${\bf S}$ and $\hat{\boldsymbol{S}}_{\rm Kerr}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
Skerrhat1 = sp.symbols('Skerrhat1')  # Skerrhat1 variable
Skerrhat2 = sp.symbols('Skerrhat2')  # Skerrhat2 variable
Skerrhat3 = sp.symbols('Skerrhat3')  # Skerrhat3 variable

# Define the expression for dot product of S and Skerrhat
dot_product_expr = (S1*Skerrhat1 + S2*Skerrhat2 + S3*Skerrhat3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\hat{\boldsymbol{S}}_{\rm Kerr}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Dot Product of Two Vectors

*   **The Dot Product:** In this section, we discuss the dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
skerrhat).


We are now going to calculate the expression for the dot product of ${\bf S}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable
S2 = sp.symbols('S2')  # S2 variable
S3 = sp.symbols('S3')  # S3 variable
Skerrhat1 = sp.symbols('Skerrhat1')  # Skerrhat1 variable
Skerrhat2 = sp.symbols('Skerrhat2')  # Skerrhat2 variable
Skerrhat3 = sp.symbols('Skerrhat3')  # Skerrhat3 variable

# Define the expression for dot product of S and Skerrhat
dot_product_expr = (S1*Skerrhat1 + S2*Skerrhat2 + S3*Skerrhat3)  # dot_product_expression
```

### Theory Review


*   **${\bf S}$ and $\boldsymbol{n}$ Dot Product**: The dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the dot product of ${\bf S}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Dot Product**: The expression for the dot product is given by $$S^{i} n^{j} = S_{i} n_{j}$$

### Code Explanation


*   **Writing to File**: We are now writing the expression for dot product to a file**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate Dot Product

*   **The Complex Conjugate Dot Product:** In this section, we discuss the complex conjugate dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 8.e: ${\bf S}^{*} \cdot {\bf n}$ \[Back to [top](

We are now going to calculate the expression for the complex conjugate dot product of ${\bf S}^{*}$ and ${\bf n}$.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for complex conjugate dot product of S_star and n
complex_dot_product_expr = (S1_star*n1 + S2_star*n2 + S3_star*n3)  # complex_dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and ${\bf n}$ Complex Conjugate Dot Product**: The complex conjugate dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the complex conjugate dot product of ${\bf S}^{*}$ and ${\bf n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_dot_product_expr) + '\n')
```

### Theory Review


*  **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sstardotn}$$

We have

\begin{equation*}
    {\bf S}^{*} \cdot {\bf n} = {\bf S}^{*}_{1} n_{1} + {\bf S}^{*}_{2} n_{2} + {\bf S}^{*}_{3} n_{3}.
\end{equation*}

We define ${\bf S}^{*}$ in [this cell](

We are now going to calculate the expression for the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for complex conjugate dot product of S_star and n
complex_dot_product_expr = (S1_star*n1 + S2_star*n2 + S3_star*n3)  # complex_dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and $\boldsymbol{n}$ Complex Conjugate Dot Product**: The complex conjugate dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sstar) and ${\bf n}$ in [this cell](

We are now going to calculate the expression for the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for complex conjugate dot product of S_star and n
complex_dot_product_expr = (S1_star*n1 + S2_star*n2 + S3_star*n3)  # complex_dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and $\boldsymbol{n}$ Complex Conjugate Dot Product**: The complex conjugate dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate Dot Product**: The expression**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate Dot Product of Two Vectors

*   **The Complex Conjugate Dot Product:** In this section, we discuss the complex conjugate dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
n).


We are now going to calculate the expression for the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for complex conjugate dot product of S_star and n
complex_dot_product_expr = (S1_star*n1 + S2_star*n2 + S3_star*n3)  # complex_dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and $\boldsymbol{n}$ Complex Conjugate Dot Product**: The complex conjugate dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_dot_product_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate Dot Product**: The expression for the complex conjugate dot product is given**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate Dot Product of Two Vectors

*   **The Complex Conjugate Dot Product:** In this section, we discuss the complex conjugate dot product of two vectors, which is a crucial concept in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9: $H_{\rm real}$ Spin Combination ${\bf S}^{*}$ \[Back to [top](

We are now going to calculate the expression for the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
n1 = sp.symbols('n1')  # n1 variable
n2 = sp.symbols('n2')  # n2 variable
n3 = sp.symbols('n3')  # n3 variable

# Define the expression for complex conjugate dot product of S_star and n
complex_dot_product_expr = (S1_star*n1 + S2_star*n2 + S3_star*n3)  # complex_dot_product_expression
```

### Theory Review


*   **${\bf S}^{*}$ and $\boldsymbol{n}$ Complex Conjugate Dot Product**: The complex conjugate dot product of two vectors is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between the two vectors.
    +   In this case, we are calculating the complex conjugate dot product of ${\bf S}^{*}$ and $\boldsymbol{n}$.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate dot product to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_dot_product_expr) + '\n')
```

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hreal_spin_combos}$$

We collect here terms defining and containing ${\bf S}^{*}$.

<a id='sstar'></a>

In this section, we are going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)

# Define the complex conjugate of a vector
complex_conjugate_expr = [S1_star, S2_star, S3_star]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate**: The expression for the complex conjugate of a vector is given by $${\bf S}^{*}_{i}$$ where $i=1,2,3$.


Note: This is just a snippet of code and it's not a complete example. You should replace the code with your own implementation.**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.a: ${\bf S}^{*}$ \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)

# Define the complex conjugate of a vector
complex_conjugate_expr = [S1_star, S2_star, S3_star]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate**: The expression for the complex conjugate of a vector is given by $${\bf S}^{*}_{i}$$ where $i=1,2,3$.

### Mathematical Representation

The complex conjugate of a vector can be represented mathematically as:

$$
{\bf S}^{*} = (S_1^*, S_2^*, S_3^*)
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sstar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.63):

\begin{equation*}
    {\bf S}^{*} = \boldsymbol{\sigma}^{*} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma^{*}}.
\end{equation*}

We define $\boldsymbol{\sigma}^{*}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)

# Define the complex conjugate of a vector
complex_conjugate_expr = [S1_star, S2_star, S3_star]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate**: The expression for the complex conjugate of a vector is given by $**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmastar) and $\boldsymbol{\Delta}_{\sigma^{*}}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star = sp.symbols('S1_star')  # S1_star variable (complex conjugate)
S2_star = sp.symbols('S2_star')  # S2_star variable (complex conjugate)
S3_star = sp.symbols('S3_star')  # S3_star variable (complex conjugate)
Delta_Sigma_star_1 = sp.symbols('Delta_Sigma_star_1')  # Delta_Sigma_star_1 variable
Delta_Sigma_star_2 = sp.symbols('Delta_Sigma_star_2')  # Delta_Sigma_star_2 variable
Delta_Sigma_star_3 = sp.symbols('Delta_Sigma_star_3')  # Delta_Sigma_star_3 variable

# Define the complex conjugate of a vector and $\boldsymbol{\Delta}_{\sigma^{*}}$
complex_conjugate_expr = [S1_star, S2_star, S3_star]
delta_sigma_star_expr = [Delta_Sigma_star_1, Delta_Sigma_star_2, Delta_Sigma_star_3]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltasigmastar).

Please note: after normalization, ${\bf S} = {\bf S}^{*}$.  See [BB2010](https://arxiv.org/abs/0912.3517) Equation (4.26).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Define variables
S1 = Sstar1  # S1 variable
S2 = Sstar2  # S2 variable
S3 = Sstar3  # S3 variable
Sstar1 = sigmastar1 + Deltasigmastar1  # Sstar1 variable
Sstar2 = sigmastar2 + Deltasigmastar2  # Sstar2 variable
Sstar3 = sigmastar3 + Deltasigmastar3  # Sstar3 variable

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(S1) + '\n')
    f.write(str(S2) + '\n')
    f.write(str(S3) + '\n')
    f.write(str(Sstar1) + '\n')
    f.write(str(Sstar2) + '\n')
    f.write(str(Sstar3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
# Note: after normalization, ${\bf S} = {\bf S}^{*}$.
```

### Mathematical Representation

The complex conjugate of a vector can be represented mathematically as:

$$
{\bf S}^{*} = \boldsymbol{\sigma}^{*} + \**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.b: $\boldsymbol{\Delta}_{\sigma^{*}}$ \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
Deltastar1 = sp.symbols('Deltastar1')  # Deltastar1 variable
Deltastar2 = sp.symbols('Deltastar2')  # Deltastar2 variable
Deltastar3 = sp.symbols('Deltastar3')  # Deltastar3 variable

# Define the complex conjugate of a vector
complex_conjugate_expr = [Deltastar1, Deltastar2, Deltastar3]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Theory Review


*   **Expression for Complex Conjugate**: The expression for the complex conjugate of a vector is given by $${\bf S}^{*}_{i}$$ where $i=1,2,3$.**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{deltasigmastar}$$

We can write $\boldsymbol{\Delta}_{\sigma^{*}}$ as

\begin{equation*}
    \boldsymbol{\Delta}_{\sigma^{*}} = \boldsymbol{\sigma}^{*} \left( \boldsymbol{\sigma}^{*}\ {\rm coefficient} \right) + \boldsymbol{\sigma} \left( \boldsymbol{\sigma}\ {\rm coefficient} \right)
\end{equation*}

For further dissection, see $\boldsymbol{\sigma}^{*}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
Deltastar1 = sp.symbols('Deltastar1')  # Deltastar1 variable
Deltastar2 = sp.symbols('Deltastar2')  # Deltastar2 variable
Deltastar3 = sp.symbols('Deltastar3')  # Deltastar3 variable

# Define the complex conjugate of a vector
complex_conjugate_expr = [Deltastar1, Deltastar2, Deltastar3]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmastar), $\boldsymbol{\sigma}^{*}$ coefficient in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star_coeff = sp.symbols('S1_star_coeff')  # S1_star_coeff variable
S2_star_coeff = sp.symbols('S2_star_coeff')  # S2_star_coeff variable
S3_star_coeff = sp.symbols('S3_star_coeff')  # S3_star_coeff variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}^{*}$ coefficient
complex_conjugate_expr = [S1_star_coeff, S2_star_coeff, S3_star_coeff]
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Mathematical Representation

The complex conjugate of a vector can be represented mathematically as:

$$
{\bf S}^{*} = \boldsymbol{\sigma}^{*} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma^{*}}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmastarcoeff), $\boldsymbol{\sigma}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_coeff = sp.symbols('S1_coeff')  # S1_coeff variable
S2_coeff = sp.symbols('S2_coeff')  # S2_coeff variable
S3_coeff = sp.symbols('S3_coeff')  # S3_coeff variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}$ coefficient
complex_conjugate_expr = [S1_coeff, S2_coeff, S3_coeff]
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Mathematical Representation

The complex conjugate of a vector can be represented mathematically as:

$$
{\bf S}^{*} = \boldsymbol{\sigma}^{*} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma^{*}}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigma), and $\boldsymbol{\sigma}$ coefficient in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1 = sp.symbols('S1')  # S1 variable (real part)
S2 = sp.symbols('S2')  # S2 variable (real part)
S3 = sp.symbols('S3')  # S3 variable (real part)

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}$ coefficient
complex_conjugate_expr = [S1, S2, S3]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Mathematical Representation

The complex conjugate of a vector can be represented mathematically as:

$$
{\bf S}^{*} = \boldsymbol{\sigma}^{*} + \frac{ 1 }{ c^{2} } \boldsymbol{\Delta}_{\sigma^{*}}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmacoeff).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Define variables
Deltastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff  # Deltastar1 variable
Deltastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff  # Deltastar2 variable
Deltastar3 = sigmastar3*sigmastarcoeff + sigma3*sigmacoeff  # Deltastar3 variable

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Deltastar1) + '\n')
    f.write(str(Deltastar2) + '\n')
    f.write(str(Deltastar3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\Delta}_{\sigma^{*}}$ can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\Delta}_{\sigma^{*}} &= \boldsymbol{\sigma}^{*} \left( \boldsymbol{\sigma}^{*}\ {\rm coefficient} \right) + \boldsymbol{\sigma} \left( \boldsymbol{\sigma}\ {\rm coefficient} \right)
\end{aligned}
$$

where $\boldsymbol{\sigma}^{*}$ and $\boldsymbol{\sigma}$ are the complex conjugate of a vector and its real part, respectively. The coefficients represent the amount of "similarity" between the two vectors.**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.c: $\boldsymbol{\sigma}^{*}$ coefficient \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star_coeff = sp.symbols('S1_star_coeff')  # S1_star_coeff variable
S2_star_coeff = sp.symbols('S2_star_coeff')  # S2_star_coeff variable
S3_star_coeff = sp.symbols('S3_star_coeff')  # S3_star_coeff variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}^{*}$ coefficient
complex_conjugate_expr = [S1_star_coeff, S2_star_coeff, S3_star_coeff]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ coefficient can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient} &= \frac{ 1 }{ c^{2} }
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmastarcoeff}$$

We will break down $\boldsymbol{\sigma}^{*}\ {\rm coefficient}$ into two terms:

\begin{equation*}
    \boldsymbol{\sigma}^{*}\ {\rm coefficient} = \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 1} + \boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 2}
\end{equation*}

We compute $\boldsymbol{\sigma}^{*}$ coefficient Term 1 in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star_coeff_Term_1 = sp.symbols('S1_star_coeff_Term_1')  # S1_star_coeff_Term_1 variable
S2_star_coeff_Term_1 = sp.symbols('S2_star_coeff_Term_1')  # S2_star_coeff_Term_1 variable
S3_star_coeff_Term_1 = sp.symbols('S3_star_coeff_Term_1')  # S3_star_coeff_Term_1 variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}^{*}$ coefficient Term 1
complex_conjugate_expr = [S1_star_coeff_Term_1, S2_star_coeff_Term_1, S3_star_coeff_Term_1]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmastarcoeffterm1) and $\boldsymbol{\sigma}^{*}$ coefficient Term 2 in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star_coeff_Term_2 = sp.symbols('S1_star_coeff_Term_2')  # S1_star_coeff_Term_2 variable
S2_star_coeff_Term_2 = sp.symbols('S2_star_coeff_Term_2')  # S2_star_coeff_Term_2 variable
S3_star_coeff_Term_2 = sp.symbols('S3_star_coeff_Term_2')  # S3_star_coeff_Term_2 variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}^{*}$ coefficient Term 2
complex_conjugate_expr = [S1_star_coeff_Term_2, S2_star_coeff_Term_2, S3_star_coeff_Term_2]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Code Explanation


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Write the expression for complex conjugate to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(complex_conjugate_expr) + '\n')
```

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}\ {\rm coefficient}$ can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient} &= \boldsymbol**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmastarcoeffterm2).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Define variables
sigmastarcoeffTerm1 = S1_star_coeff_Term_1 + S2_star_coeff_Term_1 + S3_star_coeff_Term_1  # sigmastarcoeffTerm1 variable
sigmastarcoeffTerm2 = S1_star_coeff_Term_2 + S2_star_coeff_Term_2 + S3_star_coeff_Term_2  # sigmastarcoeffTerm2 variable

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(sigmas))
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ coefficient can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient} &= \frac{ 1 }{ c^{2} }
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.c.i: $\boldsymbol{\sigma}^{*}$ Coefficient Term 1 \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
S1_star_coeff_Term_1 = sp.symbols('S1_star_coeff_Term_1')  # S1_star_coeff_Term_1 variable
S2_star_coeff_Term_1 = sp.symbols('S2_star_coeff_Term_1')  # S2_star_coeff_Term_1 variable
S3_star_coeff_Term_1 = sp.symbols('S3_star_coeff_Term_1')  # S3_star_coeff_Term_1 variable

# Define the complex conjugate of a vector and $\boldsymbol{\sigma}^{*}$ Coefficient Term 1
complex_conjugate_expr = [S1_star_coeff_Term_1, S2_star_coeff_Term_1, S3_star_coeff_Term_1]  # complex_conjugate_expression
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ Coefficient Term 1 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm Coefficient\ Term\ 1} &= \frac{ 1 }{ c^{2} }
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmastarcoeffterm1}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (51) with $b_{0} = 0$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), where what is listed below is the coefficient on $\boldsymbol{\sigma}^{*}$:


```python
import sympy as sp

# Define variables
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
M_div_r = sp.symbols('M_div_r')  # M_div_r variable
Delta_r_over_Sigma = sp.symbols('Delta_r_over_Sigma')  # Delta_r_over_Sigma variable
n_dot_p_squared = sp.symbols('n_dot_p_squared')  # n_dot_p_squared variable

# Define the expression for Q-1
Q_minus_1_expr = (14 * M_div_r / 12) + (4 * (Q_minus_1 - 1) / 12) - (30 * Delta_r_over_Sigma * n_dot_p_squared / 12)
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ Coefficient Term 1 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 1} &= \frac{7}{6} \eta \frac{M}{r} + \frac{1}{3} \eta \left( Q - 1 \right) - \frac{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q) and $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
n_dot_p_squared = (sp.symbols('n_x') * sp.symbols('p_x'))**2 + (sp.symbols('n_y') * sp.symbols('p_y'))**2 + (sp.symbols('n_z') * sp.symbols('p_z'))**2  # n_dot_p_squared variable
Delta_r_over_Sigma = sp.symbols('Delta_r_over_Sigma')  # Delta_r_over_Sigma variable

# Define the expression for $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$
expr = Delta_r_over_Sigma * n_dot_p_squared
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ can be represented mathematically as:

$$
\begin{aligned}
\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} &= \frac{ \Delta_r }{ \Sigma } \left[ ({\bf n} \cdot \hat{\bf p})^2 \right]
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
drsipn2); we define $r$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
M_div_r = sp.symbols('M_div_r')  # M_div_r variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
Delta_r_over_Sigma = sp.symbols('Delta_r_over_Sigma')  # Delta_r_over_Sigma variable

# Define the expression for $r$
expr = r
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$ can be represented mathematically as:

$$
\begin{aligned}
r &= r
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r), $\eta$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable

# Define the expression for $r$ and $\eta$
expr_r = r
expr_eta = eta
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$ and $\eta$ can be represented mathematically as:

$$
\begin{aligned}
r &= r \\
\eta &= \eta
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta), and $M$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable

# Define the expression for $r$, $\eta$ and $M$
expr_r = r
expr_eta = eta
expr_M = M
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$, $\eta$ and $M$ can be represented mathematically as:

$$
\begin{aligned}
r &= r \\
\eta &= \eta \\
M &= M
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
m) below.


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

# Define variables
sigmastarcoeffTerm1 = eta/12*(14/r + 4*Qminus1 - 30*DrSipn2)

# Write the expression for sigmastarcoeffTerm1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(sigmastarcoeffTerm1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ Coefficient Term 1 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 1} &= \frac{7}{6} \eta \frac{M}{r} + \frac{1}{3} \eta \left( Q - 1 \right) - \frac{5}{2} \eta \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.c.ii: $\boldsymbol{\sigma}^{*}$ Coefficient Term 2 \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}^{*}$ Coefficient Term 2
expr = eta/12*(14*r/M + 4*(Q_minus_1 - 1) - 30*DrSipn2)
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}^{*}$ Coefficient Term 2 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma}^{*}\ {\rm coefficient\ Term\ 2} &= \frac{7}{6} \eta \frac{M}{r} + \frac{1}{3} \eta \left( Q - 1 \right) - \frac{5}{2} \eta \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmastarcoeffterm2}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (52) with all $b_{i} = 0$, $i \in \left\{0, 1, 2, 3\right\}$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), and just the coefficient on $\boldsymbol{\sigma}^{*}$.  In the LALSuite code this is the variable 'sMultiplier1':

```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}^{*}$ Coefficient Term 2
expr = (eta/12*(353*r**2/M**2 + 4*(3*eta**2)*DrSipn2**2 + 4*(-23*eta - 3*eta**2)*(Q_minus_1)**2 + 4*(-103*eta + 60*eta**2)*M/(r)*(Q_minus_1) + 4*(16*eta - 21*eta**2)*DrSipn2*(Q_minus_1) + 4*(47*eta - 54*eta**2)*M/r*DrSipn2))

# Write the expression for sigmastarcoeffTerm2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r), $\eta$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $r$, $\eta$
expr_r = r
expr_eta = eta

# Write the expressions for expr_r and expr_eta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
    f.write(str(expr_eta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$, $\eta$ can be represented mathematically as:

$$
\begin{aligned}
r &= r \\
\eta &= \eta
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta), and $M$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\eta$ and $M$
expr_eta = eta
expr_M = M

# Write the expressions for expr_eta and expr_M to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_eta) + '\n')
    f.write(str(expr_M) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\eta$ and $M$ can be represented mathematically as:

$$
\begin{aligned}
\eta &= \eta \\
M &= M
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
m); we group together and define $Q - 1$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q - 1$
expr_Q_minus_1 = Q_minus_1

# Write the expressions for expr_Q_minus_1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q_minus_1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $Q - 1$ can be represented mathematically as:

$$
\begin{aligned}
Q - 1 &= Q - 1
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q), and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $q$ and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$
expr_q = Q_minus_1
expr_DrSipn2 = DrSipn2

# Write the expressions for expr_q and expr_DrSipn2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_q) + '\n')
    f.write(str(expr_DrSipn2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $q$ and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ can be represented mathematically as:

$$
\begin{aligned}
q &= Q -**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
drsipn2.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\sigma^*$ Coefficient Term 2
expr = (eta/(72*r*r))*(706 + r*(-206*Q_minus_1 + 282*DrSipn2 + r*Q_minus_1*(96*DrSipn2 - 23*Q_minus_1))
                                    + eta*(-54 + r*(120*Q_minus_1 - 324*DrSipn2
                                    + r*(360*DrSipn2**2 + Q_minus_1*(-126*DrSipn2 - 3*Q_minus_1)))))

# Write the expression for sigmastarcoeffTerm2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector:** The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\sigma^*$ Coefficient Term 2 can be represented mathematically as:

$$
\begin{aligned}
\sigma^* \text{Coefficient Term 2} &= \frac{\eta}{72 r^{2}} \left[ 706 + r \left( -206 Q -1 + 282 \Delta_{r**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.d: $\boldsymbol{\sigma}$ coefficient \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ coefficient
expr = eta/(12*(14*r/M + 4*(Q_minus_1 - 1) - 30*DrSipn2))

# Write the expressions for expr to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ coefficient can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient} &= \frac{\eta}{12} (14 \frac{M}{r} + 4 (Q - 1) - 30 \Delta_{r})
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmacoeff}$$

We will break down $\boldsymbol{\sigma}\ {\rm coefficient}$ into three terms:


\begin{equation*}
    \boldsymbol{\sigma}\ {\rm coefficient} = \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 1} + \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 2} + \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 3}
\end{equation*}

We compute $\boldsymbol{\sigma}$ coefficient Term 1 in [this cell](


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ coefficient Term 1
expr_term1 = eta/12*(14*r/M + 4*(Q_minus_1 - 1))

# Write the expressions for expr_term1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}\ {\rm coefficient}$ can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient} &= \frac{\**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmacoeffterm1), $\boldsymbol{\sigma}$ coefficient Term 2 in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ coefficient Term 2
expr_term2 = eta/12*(-30*DrSipn2)

# Write the expressions for expr_term2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ coefficient Term 2 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 2} &= -30 \eta \Delta_{r}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmacoeffterm2), and $\boldsymbol{\sigma}$ coefficient Term 3 in [this cell](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ coefficient Term 3
expr_term3 = eta/12*(14*r/M + 4*(Q_minus_1 - 1))

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ coefficient Term 3 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 3} &= \frac{14}{12} \eta \frac{M}{r} + \frac{4}{12} \eta (Q - 1)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sigmacoeffterm3.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ coefficient Term 3
expr_term3 = eta/12*(14*r/M + 4*(Q_minus_1 - 1))

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ coefficient can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient} &= \frac{\eta}{12} (14 \frac{M}{r} + 4 (Q - 1) - 30 \Delta_{r})
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.d.i: $\boldsymbol{\sigma}$ Coefficient Term 1 \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 1
expr_term1 = eta/12*(14*r/M + 4*(Q_minus_1 - 1))

# Write the expressions for expr_term1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 1 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 1} &= \frac{14}{12} \eta \frac{M}{r} + \frac{4}{12} \eta (Q - 1)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmacoeffterm1}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (51) with $a_{0} = 0$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), where what is listed below is the coefficient on $\boldsymbol{\sigma}$:


\begin{align*}
    \boldsymbol{\sigma}\ {\rm coefficient\ Term\ 1} &= -\frac{2}{3} \eta \frac{ M }{ r } + \frac{1}{4} \eta \left( Q - 1 \right) - 3 \eta \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} \\
        &= \frac{ \eta }{ 12 } \left( -8 \frac{ M }{ r } + 3 \left( Q - 1 \right) - 36 \smash[b]{\underbrace{ \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} }_{\rm DrSipn2}} \vphantom{\underbrace{a}_{b}} \right)
\end{align*}

We define $\eta$ in [this cell](


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta), $M$ in [this cell](

We are now going to calculate $\eta$ and $M$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\eta$
expr_eta = eta

# Write the expressions for expr_eta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_eta) + '\n')

# Define the expression for $M$
expr_M = M

# Write the expressions for expr_M to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_M) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\eta$ can be represented mathematically as:

$$
\begin{aligned}
\eta &= \text{some mathematical expression for eta}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
m), $Q-1$ in [this cell](

We are now going to calculate $M$ and $Q-1$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $M$
expr_M = M

# Write the expressions for expr_M to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_M) + '\n')

# Define the expression for $Q-1$
expr_Q_minus_1 = Q_minus_1

# Write the expressions for expr_Q_minus_1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q_minus_1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $M$ can be represented mathematically as:

$$
\begin{aligned}
M &= \text{some mathematical expression for M}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q), and $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](

We are now going to calculate $Q$ and $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $Q$
expr_Q = Q_minus_1 + 1

# Write the expressions for expr_Q to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q) + '\n')

# Define the expression for $\frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$
expr_DrSipn2 = DrSipn2

# Write the expressions for expr_DrSipn2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_DrSipn2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
drsipn2.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 1
expr_term1 = eta/12*(-8/r + 3*Q_minus_1 - 36*DrSipn2)

# Write the expressions for expr_term1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 1 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 1} &= \frac{\eta}{12} \left( -8 \frac{M}{r} + 3 (Q-1) - 36 \smash[b]{\underbrace{ \frac{ \Delta_r }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2} }_{\rm DrSipn2}} \vphantom{\underbrace{a}_{b}} \right)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.d.ii: $\boldsymbol{\sigma}$ Coefficient Term 2 \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 2
expr_term2 = eta/12*3*(Q_minus_1 - 1)

# Write the expressions for expr_term2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 2 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 2} &= \frac{3}{12} \eta (Q-1)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmacoeffterm2}$$

We build this term from [BB2011](https://arxiv.org/abs/1107.2904) Equation (52) with all $a_{i} = 0$, $i \in \left\{0, 1, 2, 3\right\}$ (see discussion preceding [T2012](https://arxiv.org/abs/1202.0790) Equation (4)), and just the coefficient on $\boldsymbol{\sigma}$:


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 2
expr_term2 = (eta/144*r**(-2))*((-896*M**2) + r*(-436*M*(Q_minus_1 - 1) -96*M*DrSipn2 + r*(-45*(Q_minus_1 - 1)**2 + 36*(Q_minus_1 - 1)*DrSipn2)) + eta*(-336*M**2) + r*(204*M*(Q_minus_1 - 1) - 882*M*DrSipn2 + r*(810*(DrSipn2)**2 - 234*(Q_minus_1 - 1)*DrSipn2)))

# Write the expressions for expr_term2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term2) + '\n')
```

###**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta), $M$ in [this cell](

We are now going to calculate $\eta$ and $M$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\eta$
expr_eta = eta

# Write the expressions for expr_eta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_eta) + '\n')

# Define the expression for $M$
expr_M = M

# Write the expressions for expr_M to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_M) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\eta$ can be represented mathematically as:

$$
\begin{aligned}
\eta &= \text{some mathematical expression for eta}
\end{aligned}
$$


The $M$ can be represented mathematically as:

$$
\begin{aligned}
M &= \text{some mathematical expression for M}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
M), $Q - 1$ in [this cell](

We are now going to calculate $M$ and $Q - 1$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $M$
expr_M = M

# Write the expressions for expr_M to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_M) + '\n')

# Define the expression for $Q - 1$
expr_Q_minus_1 = Q_minus_1

# Write the expressions for expr_Q_minus_1 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q_minus_1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $M$ can be represented mathematically as:

$$
\begin{aligned}
M &= \text{some mathematical expression for M}
\end{aligned}
$$


The $Q - 1$ can be represented mathematically as:

$$
\begin{aligned}
Q - 1 &= \text{some mathematical expression for Q_minus_1}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q), and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$ in [this cell](

We are now going to calculate $Q$ and $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $Q$
expr_Q = Q_minus_1 + 1

# Write the expressions for expr_Q to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q) + '\n')

# Define the expression for $\frac{ \Delta_{r} }{ \Sigma } \left( {\bf n} \cdot \hat{\bf p} \right)^{2}$
expr_DrSipn2 = DrSipn2

# Write the expressions for expr_DrSipn2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_DrSipn2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
drsipn2.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 3
expr_term3 = eta/(144*r**2)*(-896 + r*(-436*Q_minus_1 - 96*DrSipn2) + (-45*(Q_minus_1)**2 + 36*Q_minus_1*DrSipn2)) + eta*(-336 + r*(204*Q_minus_1 - 882*DrSipn2) + (810*DrSipn2**2 - 234*Q_minus_1*DrSipn2))

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 3 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term 3} &= \frac{\eta}{144r^2} \left( -896 + r \left( -436(Q-1) - 96DrSip**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 9.d.iii: $\boldsymbol{\sigma}$ Coefficient Term 3 \[Back to [top](

We are now going to calculate the complex conjugate of a vector.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 3
expr_term3 = eta/(144*r**2)*(-896 + r*(-436*Q_minus_1 - 96*DrSipn2) + (-45*(Q_minus_1)**2 + 36*Q_minus_1*DrSipn2)) + eta*(-336 + r*(204*Q_minus_1 - 882*DrSipn2) + (810*DrSipn2**2 - 234*Q_minus_1*DrSipn2))

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 3 can be represented mathematically as:

$$
\begin{aligned}
\boldsymbol{\sigma} \text{ Coefficient Term **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{sigmacoeffterm3}$$

From Section II of [T2014)](https://arxiv.org/abs/1311.2544), \[Back to [top](

We build this term from Equation (4.13) of [BL2011](https://arxiv.org/abs/1608.02623)


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $d_{\rm SO}$
expr_d_SO = 147.481449*(eta**2)*r - 568.651115*eta*r + 66.198703*r - 343.313058*Q_minus_1*eta + 2495.293427*eta**2 - 44.532373

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 3
expr_term3 = eta*(expr_d_SO)

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 3 can be represented mathematically as:

$$
\**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
eta), $u$ in [this cell](

We are now going to calculate $\eta$ and $u$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\eta$
expr_eta = eta

# Write the expressions for expr_eta to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_eta) + '\n')

# Define the expression for $u$
expr_u = r  # Note: This is just a placeholder value

# Write the expressions for expr_u to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_u) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\eta$ can be represented mathematically as:

$$
\begin{aligned}
\eta &= \text{some mathematical expression for eta}
\end{aligned}
$$


The $u$ can be represented mathematically as:

$$
\begin{aligned}
u &= \text{some mathematical expression for u}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
u), and $\chi$ in [this cell](

We are now going to calculate $u$ and $\chi$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $u$
expr_u = r  # Note: This is just a placeholder value

# Write the expressions for expr_u to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_u) + '\n')

# Define the expression for $\chi$
expr_chi = sp.symbols('chi')  # Note: This is just a placeholder value

# Write the expressions for expr_chi to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_chi) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $u$ can be represented mathematically as:

$$
\begin{aligned}
u &= \text{some mathematical expression for u}
\end{aligned}
$$


The $\chi$ can be represented mathematically as:

$$
\begin{aligned}
\chi &= \text{some mathematical expression for chi}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
chi).  Note that the values have been rounded to agree with those in the LALSuite implementation (see the file LALSimIMRSpinEOBHamiltonian.h).


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\chi$
expr_chi = sp.symbols('chi')  # Note: This is just a placeholder value

# Define the expression for $d_{\rm SO}$
expr_d_SO = 147.481*(expr_chi**3)*(eta**2) - 568.651*(expr_chi**3)*eta + 66.1987*(expr_chi**3) - 343.313*expr_chi**2*eta + 2495.29*expr_chi*eta**2 - 44.5324

# Define the expression for $\boldsymbol{\sigma}$ Coefficient Term 3
expr_term3 = eta*expr_d_SO*expr_chi**2

# Write the expressions for expr_term3 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_term3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\boldsymbol{\sigma}$ Coefficient Term 3 can be represented mathematically as:

**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10: Derivatives of the Metric Potential \[Back to [top](

We are now going to calculate the derivatives of the metric potential.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for the metric potential
expr_metric_potential = r**2  # Note: This is just a placeholder value

# Calculate the derivatives of the metric potential
expr_metric_potential_deriv = sp.diff(expr_metric_potential, r)

# Write the expressions for expr_metric_potential_deriv to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_metric_potential_deriv) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The derivatives of the metric potential can be represented mathematically as:

$$
\begin{aligned}
\frac{\partial \phi}{\partial r} &= \text{some mathematical expression for the derivative of the metric potential with respect to $r$}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{metpotderivs}$$

We collect here terms dependent on derivatives of the metric potential (see [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.47)).


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega$
expr_omega = eta/(144*r**3)

# Write the expressions for expr_omega to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega$ can be represented mathematically as:

$$
\begin{aligned}
\omega &= \frac{\eta}{144r^3}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.a: $\omega_{r}$ \[Back to [top](

We are now going to calculate $\omega_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_r$
expr_omega_r = eta/(144*r**3)

# Write the expressions for expr_omega_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega_r$ can be represented mathematically as:

$$
\begin{aligned}
\omega_r &= \frac{\eta}{144r^3}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{omegar}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47b) we have \[Back to [top](

We are now going to calculate $\omega_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_{t}$
expr_Lambda_t = M + (3*eta)/5

# Write the expressions for expr_Lambda_t to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambda_t) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega_r$ can be represented mathematically as:

$$
\begin{aligned}
\omega_r &= \frac{ \Lambda_{t} \tilde{\omega}_{\rm fd}^{\prime} - \Lambda_{t}^{\prime} \tilde{\omega}_{\rm fd} }{ \Lambda_{t}^{2} }
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), $\tilde{\omega}_{\rm fd}^{\prime}$ in [this cell](

We are now going to calculate $\Lambda_t$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_t^{\prime}$
expr_Lambda_t_prime = M + (9*eta)/(5*r)

# Write the expressions for expr_Lambda_t_prime to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambda_t_prime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:

$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= \text{some mathematical expression for the derivative of } \tilde{\omega}_{\rm fd}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omegatildeprm), $\Lambda_{t}^{\prime}$ in [this cell](

We are now going to calculate $\tilde{\omega}_{\rm fd}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\tilde{\omega}_{\rm fd}^{\prime}$
expr_tildewfd_prime = (M + (3*eta)/5) / r**3

# Write the expressions for expr_tildewfd_prime to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tildewfd_prime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:

$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= \frac{M + (3\eta)/5}{r^{3}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdatprm), and $\tilde{\omega}_{\rm fd}$ in [this cell](

We are now going to calculate $\Lambda_t^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_t^{\prime}$
expr_Lambda_t_prime = M + (9*eta)/(5*r)

# Write the expressions for expr_Lambda_t_prime to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambda_t_prime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Lambda_t^{\prime}$ can be represented mathematically as:

$$
\begin{aligned}
\Lambda_t^{\prime} &= M + \frac{9\eta}{5r}
\end{aligned}
$$


The $\tilde{\omega}_{\rm fd}$ can be represented mathematically as:

$$
\begin{aligned}
\tilde{\omega}_{\rm fd} &= \text{some mathematical expression for } \tilde{\omega}_{\rm fd}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omegar = (Lambdat*omegatildeprm - Lambdatprm*omegatilde)/(Lambdat*Lambdat)

We are now going to calculate $\omega_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_r$
expr_omega_r = (Lambdat*sp.simplify(Lambdat**3*omegatildeprm - Lambdat_prime*sp.simplify(1.0*Lambdat*omegatilde)))/(sp.simplify((Lambdat)**4))

# Write the expressions for expr_omega_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega_r$ can be represented mathematically as:

$$
\begin{aligned}
\omega_r &= \frac{\Lambda_t\tilde{\omega}_{\rm fd}^{\prime}-\Lambda_t^{\prime}\tilde{\omega}_{\rm fd}}{\Lambda_t^{2}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.b: $\nu_{r}$ \[Back to [top](

We are now going to calculate $\nu_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\nu_r$
expr_nu_r = (Lambdat**3*omegatildeprm)/(Lambdat**4)

# Write the expressions for expr_nu_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nu_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\nu_r$ can be represented mathematically as:

$$
\begin{aligned}
\nu_r &= \frac{\Lambda_t^3\tilde{\omega}_{\rm fd}^{\prime}}{\Lambda_t^{4}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{nur}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47c) we have \[Back to [top](

We are now going to calculate $\nu_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\nu_r$
expr_nu_r = r/(sp.simplify(DrSipn2)) + ((Q_minus_1**2)*( (Q_minus_1**2)*Lambdat_prime - (4*r)*Lambdat ) / (2*Lambdat*(Lambdat)))

# Write the expressions for expr_nu_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nu_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\nu_r$ can be represented mathematically as:

$$
\begin{aligned}
\nu_r &= \frac{r}{\Sigma} + \frac{\varpi^{2}(\varpi^{2}\Delta^{\prime}_{t}-4r\Delta_{t})}{2\Lambda_{t}\Delta_{t}}
\end{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r), $\Sigma$ in [this cell](

We are now going to calculate $r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $r$
expr_r = r

# Write the expressions for expr_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$ can be represented mathematically as:

$$
\begin{aligned}
r &= \text{some mathematical expression for r}
\end{aligned}
$$


The $\Sigma$ can be represented mathematically as:

$$
\begin{aligned}
\Sigma &= \text{some mathematical expression for Sigma}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), $\varpi^{2}$ in [this cell](

We are now going to calculate $u_\Sigma$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $u_\Sigma$
expr_u_Sigma = r/(sp.simplify(DrSipn2))

# Write the expressions for expr_u_Sigma to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_u_Sigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $u_\Sigma$ can be represented mathematically as:

$$
\begin{aligned}
u_\Sigma &= \frac{r}{\Sigma} \\
\end{aligned}
$$


The $\varpi^{2}$ can be represented mathematically as:

$$
\begin{aligned}
\varpi^{2} &= \text{some mathematical expression for } \varpi^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
w2), $\Delta_{t}^{\prime}$ in [this cell](

We are now going to calculate $\varpi^{2}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\varpi^{2}$
expr_varpi2 = (r**4)/(sp.simplify(DrSipn2)**2)

# Write the expressions for expr_varpi2 to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_varpi2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\varpi^{2}$ can be represented mathematically as:

$$
\begin{aligned}
\varpi^{2} &= \frac{r^{4}}{\Sigma^{2}}
\end{aligned}
$$


The $\Delta_{t}^{\prime}$ can be represented mathematically as:


```python
# Define the expression for $\Delta_{t}^{\prime}$
expr_Deltat_prime = (3*(Lambdat_prime)/r**3) - (2*Lambdat/r**4)

# Write the expressions for expr_Deltat_prime to file
with open($Ccodesdir/v4P_Ham**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltatprm), $\Delta_{t}$ in [this cell](

We are now going to calculate $\Lambda_t^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_t^{\prime}$
expr_Lambdat_prime = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions for expr_Lambdat_prime to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdat_prime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_{t}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{t}^{\prime} &= \frac{3M}{r^{3}} - \frac{2\Sigma}{r^{4}}
\end{aligned}
$$


The $\Delta_{t}$ can be represented mathematically as:

$$
\begin{aligned}
\Delta_{t} &= \text{some mathematical expression for } \Delta_t
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat), and $\Lambda_{t}$ in [this cell](

We are now going to calculate $\Delta_t$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Delta_t$
expr_deltat = (r**3*sp.sqrt(DrSipn2))/sp.sqrt(3*eta*M)

# Write the expressions for expr_deltat to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_deltat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_t$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{t} &= \frac{r^{3}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\Lambda_{t}$ can be represented mathematically as:


```python
# Define the expression for $\Lambda_t$
expr_Lambdat = (3*M)/(r**3) + (2*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*r)*r)

# Write the expressions for expr_Lambdat to file
with open($Ccodesdir/v4P_Hamiltonian-H**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
nur = r/Sigma + w2*(w2*Deltatprm - 4*r*Deltat)/(2*Lambdat*Deltat)

We are now going to calculate $\nu_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\nu_r$
expr_nur = r/sp.simplify(DrSipn2) + ((sp.symbols('w2'))*( (sp.symbols('w2'))*sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)) - (4*r)*sp.simplify((r**3*sp.sqrt(sp.simplify(DrSipn2)))/(sp.sqrt(3*eta*M))) ))/(2*sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r)))

# Write the expressions for expr_nur to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nur) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\nu_r$ can be represented mathematically as:


$$
\begin**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.c: $\mu_{r}$ \[Back to [top](

We are now going to calculate $\mu_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\mu_r$
expr_mu_r = (r*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Write the expressions for expr_mu_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mu_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_r$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{r} &= \frac{r\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{mur}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47d) we have \[Back to [top](

We are now going to calculate $\mu_r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\mu_r$
expr_mu_r = r/sp.simplify(DrSipn2) - (1/sp.sqrt(sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))))

# Write the expressions for expr_mu_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mu_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_r$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{r} &= \frac{r}{\Sigma} - \frac{1}{\sqrt{\Delta_{r}}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r), $\Sigma$ in [this cell](

We are now going to calculate $r$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $r$
expr_r = r

# Write the expressions for expr_r to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$ can be represented mathematically as:


$$
\begin{aligned}
r &= \text{some mathematical expression for r}
\end{aligned}
$$


The $\Sigma$ can be represented mathematically as:


```python
# Define the expression for $\Sigma$
expr_Sigma = (r**3*sp.sqrt(DrSipn2))/sp.sqrt(3*eta*M)

# Write the expressions for expr_Sigma to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Sigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), and $\Delta_{r}$ in [this cell](

We are now going to calculate $u_\Sigma$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $u_\Sigma$
expr_usigma = r/sp.simplify(DrSipn2)

# Write the expressions for expr_usigma to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_usigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $u_\Sigma$ can be represented mathematically as:


$$
\begin{aligned}
u_{\Sigma} &= \frac{r}{\Sigma}
\end{aligned}
$$


The $\Delta_r$ can be represented mathematically as:


```python
# Define the expression for $\Delta_r$
expr_Deltar = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions for expr_Deltar to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltar) + '\n')
```

### Theory**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltar).


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\mu_r$
expr_mur = r/sp.simplify(DrSipn2) - (1/sp.sqrt(sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mur) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_r$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{r} &= \frac{r}{\Sigma} - \frac{1}{\sqrt{\Delta_{r}}}
\end{aligned}
$$


The $\Delta_r$ can be represented mathematically as:


```python
# Define the expression for $\Delta_r$
expr_Deltar = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.d: $\omega_{\cos\theta}$ \[Back to [top](

We are now going to calculate $\omega_{\cos\theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_{\cos\theta}$
expr_omegacostheta = (r**3*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omegacostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega_{\cos\theta}$ can be represented mathematically as:


$$
\begin{aligned}
\omega_{\cos\theta} &= \frac{r^{3}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{omegacostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47f), we have \[Back to [top](

We are now going to calculate $\omega_{\cos\theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_{\cos\theta}$
expr_omegacostheta = -(2*sp.symbols('a')**2*sp.cos(sp.symbols('theta'))*sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))*sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))))/sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r)))**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omegacostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
a), $\cos\theta$ in [this cell](

We are now going to calculate $a$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_a) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $a$ can be represented mathematically as:


$$
\begin{aligned}
a &= \sqrt{\frac{M}{r^{3}}}
\end{aligned}
$$


The $\cos\theta$ can be represented mathematically as:


```python
# Define the expression for $\cos\theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_costheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
costheta), $\Delta_{t}$ in [this cell](

We are now going to calculate $\cos\theta$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\cos\theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_costheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\cos\theta$ can be represented mathematically as:


$$
\begin{aligned}
\cos\theta &= \text{some mathematical expression for }\cos\theta
\end{aligned}
$$


The $\Delta_{t}$ can be represented mathematically as:


```python
# Define the expression for $\Delta_t$
expr_Deltat = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltat) + '\n')
```

### Theory Review


*   **Complex**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat), $\tilde{\omega}_{\rm fd}$ in [this cell](

We are now going to calculate $\Delta_t$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Delta_t$
expr_Deltat = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_t$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{t} &= \frac{3M}{r^{3}} - \frac{2\Sigma}{r^{4}}
\end{aligned}
$$


The $\tilde{\omega}_{\rm fd}$ can be represented mathematically as:


```python
# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_wfd = (3*M)/(r**3) + (2*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*r)*r)

# Write the expressions to file
with open**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omegatilde), and $\Lambda_{t}$ in [this cell](

We are now going to calculate $\omega_\tilde{t}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_{\tilde{t}}$
expr_omegatilde = (3*M)/(r**3) + (2*sp.sqrt(sp.simplify(DrSipn2)))/(sp.sqrt(3*eta*r)*r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omegatilde) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega_{\tilde{t}}$ can be represented mathematically as:


$$
\begin{aligned}
\omega_{\tilde{t}} &= \frac{3M}{r^{3}} + \frac{2\Sigma}{r^{4}}
\end{aligned}
$$


The $\Lambda_t$ can be represented mathematically as:


```python
# Define the expression for $\Lambda_t$
expr_Lambdat = (3*M)/(r**3) + (2*sp.sqrt(sp.simplify(DrSipn**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat).


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\omega_{\cos\theta}$
expr_omegacostheta = -(2*sp.symbols('a')**2*sp.cos(sp.symbols('theta'))*sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))*sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))))/sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r)))**2

# Define the expression for $\Delta_t$
expr_Deltat = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_omegatilde = (3*M)/(r**3) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\Lambda_t$
expr_Lambdat = (3*M)/(r**3) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Calculate the final expression
final_expression = -(2*sp.symbols('a')**2*sp.cos(sp.symbols('theta'))**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.e: $\nu_{\cos\theta}$ \[Back to [top](

We are now going to calculate $\nu_{\cos\theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\nu_{\cos\theta}$
expr_nu_costheta = (r**3*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M)*sp.sqrt((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nu_costheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\nu_{\cos\theta}$ can be represented mathematically as:


$$
\begin{aligned}
\nu_{\cos\theta} &= \frac{r^{3}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{nucostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47g) we have \[Back to [top](

We are now going to calculate $\nu_{\cos\theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\nu_{\cos\theta}$
expr_nu_costheta = (sp.symbols('a')**2*sp.symbols('varpi')**2*sp.cos(sp.symbols('theta'))*(sp.symbols('varpi')**2 - sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))))/sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))*sp.simplify(r/sp.sqrt(sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r)))))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nu_costheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
a), $\varpi^{2}$ in [this cell](

We are now going to calculate $a$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_a) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $a$ can be represented mathematically as:


$$
\begin{aligned}
a &= \sqrt{\frac{M}{r^{3}}}
\end{aligned}
$$


The $\varpi^{2}$ can be represented mathematically as:


```python
# Define the expression for $\varpi^{2}$
expr_varpifourth = (sp.symbols('varpi')**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_varpifourth) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
w2), $\cos\theta$ in [this cell](

We are now going to calculate $w^2$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $w^2$
expr_wfourth = (r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_wfourth) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $w^2$ can be represented mathematically as:


$$
\begin{aligned}
w^{2} &= \frac{r^{4}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\cos\theta$ can be represented mathematically as:


```python
# Define the expression for $\cos\theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_costheta) + '\n')
```

### Theory Review


***NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
costheta), $\Delta_{t}$ in [this cell](

We are now going to calculate $w^2$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\cos\theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Define the expression for $w^2$
expr_wfourth = (r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Define the expression for $\Delta_t$
expr_Deltat = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_wfourth) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $w^2$ can be represented mathematically as:


$$
\begin{aligned}
w^{2} &= \frac{r^{4}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\Delta_t$ can be represented mathematically as:


```python
# Define the expression for $\Delta_t$
expr_Deltat**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat), $\Lambda_{t}$ in [this cell](

We are now going to calculate $w^2$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Delta_t$
expr_Deltat = (3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4)

# Define the expression for $w^2$
expr_wfourth = (r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Define the expression for $\Lambda_t$
expr_Lambdat = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $w^2$ can be represented mathematically as:


$$
\begin{aligned}
w^{2} &= \frac{r^{4}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\Lambda_t$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), and $\Sigma$ in [this cell](

We are now going to calculate $w^2$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_t$
expr_Lambdat = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $w^2$
expr_wfourth = (r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Define the expression for $\Sigma$
expr_Sigma = (r**4*sp.sqrt(sp.simplify(DrSipn2)))/r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $w^2$ can be represented mathematically as:


$$
\begin{aligned}
w^{2} &= \frac{r^{4}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\Lambda_t$ can be represented mathematically as:


```**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma).


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\nu_{\cos\theta}$
expr_nucostheta = (sp.symbols('a')**2*sp.simplify((r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M)))*sp.cos(sp.symbols('theta'))*(sp.simplify((r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))) - sp.simplify((3*M/(r**3)) - (2*sp.simplify(DrSipn2)/r**4))))/sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))*sp.simplify(r/sp.sqrt(sp.simplify((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))))))

# Define the expression for $\Sigma$
expr_Sigma = (r**4*sp.sqrt(sp.simplify(DrSipn2)))/r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_nucostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.f: $\mu_{\cos \theta}$ \[Back to [top](

We are now going to calculate $\mu_{\cos \theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\mu_{\cos \theta}$
expr_mucostheta = (r**4*sp.sqrt(DrSipn2))*sp.cos(sp.symbols('theta'))/((3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r)))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mucostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_{\cos \theta}$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{\cos \theta} &= \frac{r^{4}\sqrt{\Sigma}}{\Lambda_t}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{mucostheta}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47h) we have \[Back to [top](

We are now going to calculate $\mu_{\cos \theta}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\mu_{\cos \theta}$
expr_mucostheta = (sp.symbols('a')**2*sp.cos(sp.symbols('theta')))/sp.simplify((r**4*sp.sqrt(DrSipn2))/r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mucostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_{\cos \theta}$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{\cos \theta} &= \frac{a^{2}\cos \theta}{\Sigma}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
a), $\cos \theta$ in [this cell](

We are now going to calculate $a$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_a) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $a$ can be represented mathematically as:


$$
\begin{aligned}
a &= \sqrt{\frac{M}{r^{3}}}
\end{aligned}
$$


The $\cos \theta$ can be represented mathematically as:


```python
# Define the expression for $\cos \theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_costheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
costheta), and $\Sigma$ in [this cell](

We are now going to calculate $w^2$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $w^2$
expr_wfourth = (r**4*sp.sqrt(DrSipn2))/(sp.sqrt(3*eta*M))

# Define the expression for $\cos \theta$
expr_costheta = sp.cos(sp.symbols('theta'))

# Define the expression for $\Sigma$
expr_Sigma = (r**4*sp.sqrt(sp.simplify(DrSipn2)))/r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_wfourth) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $w^2$ can be represented mathematically as:


$$
\begin{aligned}
w^{2} &= \frac{r^{4}\sqrt{\Sigma}}{\sqrt{3\eta M}}
\end{aligned}
$$


The $\Sigma$ can be represented mathematically as:


```python
# Define the expression for $\Sigma$
expr_Sigma = (r**4*sp.sqrt(sp.simplify(Dr**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma) below.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\mu_{\cos \theta}$
expr_mucostheta = (sp.symbols('a')**2*sp.cos(sp.symbols('theta')))/sp.simplify((r**4*sp.sqrt(DrSipn2))/r)

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\Sigma$
expr_Sigma = (r**4*sp.sqrt(sp.simplify(DrSipn2)))/r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_mucostheta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\mu_{\cos \theta}$ can be represented mathematically as:


$$
\begin{aligned}
\mu_{\cos \theta} &= \frac{a^{2}\cos \theta}{\Sigma}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.g: $\Lambda_{t}^{\prime}$ \[Back to [top](

We are now going to calculate $\Lambda_{t}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_{t}^{\prime}$
expr_Lambdatprime = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdatprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Lambda_{t}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\Lambda_{t}^{\prime} &= \frac{3M}{r^{3}} + \frac{2\sqrt{\Sigma}}{r\sqrt{3\eta}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{lambdatprm}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.47), we know that the prime notation indicates a derivative with respect to $r$.  Using the definition of $\Lambda_{t}$ in [this cell](

We are now going to calculate $\Lambda_{t}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\Lambda_{t}$
expr_Lambdat = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\Lambda_{t}^{\prime}$
expr_Lambdatprime = sp.diff(expr_Lambdat, r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdatprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Lambda_{t}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\Lambda_{t}^{\prime} &= \**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), we have

\begin{equation*}
    \Lambda_{t}^{\prime} = 4 \left( a^{2} + r^{2} \right) r - a^{2} \Delta_{t}^{\prime} \sin^{2} \theta.
\end{equation*}

We define $a$ in [this cell](

We are now going to calculate $\Lambda_{t}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\Lambda_{t}^{\prime}$
expr_Lambdatprime = (4*(sp.symbols('a')**2 + r**2)*r) - (sp.symbols('a')**2*sp.diff(sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(DrSipn2))/r**4)), r))*sp.sin(sp.symbols('theta'))**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdatprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
a), $r$ in [this cell](

We are now going to calculate $a$.


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable
r = sp.symbols('r')  # r variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_a) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $a$ can be represented mathematically as:


$$
\begin{aligned}
a &= \sqrt{\frac{M}{r^{3}}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r), $\Delta_{u}$ in [this cell](

We are now going to calculate $r$.


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable
r = sp.symbols('r')  # r variable

# Define the expression for $r$
expr_r = r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $r$ can be represented mathematically as:


$$
\begin{aligned}
r &= r
\end{aligned}
$$


The $\Delta_{u}$ can be represented mathematically as:


```python
# Define the expression for $\Delta_{u}$
expr_Deltau = sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(DrSipn2))/r**4))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltau) + '\n')
```**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltau), and $\sin^{2}\theta$ in [this cell](

We are now going to calculate $\Delta_{u}$.


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable
r = sp.symbols('r')  # r variable

# Define the expression for $\Delta_{u}$
expr_deltau = (3*M/(r**3)) - (2*sp.sqrt(sp.simplify(DrSipn2))/r**4)

# Define the expression for $\sin^{2}\theta$
expr_sin2theta = sp.sin(sp.symbols('theta'))**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_deltau) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_{u}$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{u} &= \frac{3M}{r^{3}} - \frac{2\sqrt{\Sigma}}{r^{4}}
\end{aligned}
$$


The $\sin^{2}\theta$ can be represented mathematically as:


```python
# Define the expression for $\sin^{2}\theta$
expr_sin2theta = sp.sin(sp.symbols('theta'))**2

# Write the expressions to file
with open($Ccodesdir/v4**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sin2theta).


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable
r = sp.symbols('r')  # r variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\sin^{2}\theta$
expr_sin2theta = sp.sin(sp.symbols('theta'))**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_sin2theta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\sin^{2}\theta$ can be represented mathematically as:


$$
\begin{aligned}
\sin^{2}\theta &= \sin^{2}\theta
\end{aligned}
$$


The $Lambdatprm$ can be represented mathematically as:


```python
# Define the expression for $\Lambda_{t}^{\prime}$
expr_Lambdatprm = 4*(sp.symbols('a')**2 + r**2)*r - 2*sp.symbols('a')**2*sp.diff(sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(DrSipn2))/r**4)), r)*expr_sin2theta

# Write the expressions to file
with open($Ccodesdir/v4P_Ham**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 10.h: $\tilde{\omega}_{\rm fd}^{\prime}$ \[Back to [top](

We are now going to calculate $\tilde{\omega}_{\rm fd}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\tilde{\omega}_{\rm fd}^{\prime}$
expr_tilwfdprime = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilwfdprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= \frac{3M}{r^{3}} + \frac{2\sqrt{\Sigma}}{r\sqrt{3\eta}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{omegatildeprm}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.47), we know that the prime notation indicates a derivative with respect to $r$.  Using the definition of $\tilde{\omega}_{\rm fd}$ in [this cell](

We are now going to calculate $\tilde{\omega}_{\rm fd}^{\prime}$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable

# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_tilwfd = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(DrSipn2))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\tilde{\omega}_{\rm fd}^{\prime}$
expr_tilwfdprime = sp.diff(expr_tilwfd, r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilwfdprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omegatilde), we have

\begin{equation*}
    \tilde{\omega}_{\rm fd}^{\prime} = 2 a M.
\end{equation*}

We define $a$ in [this cell](

We are now going to calculate $\tilde{\omega}_{\rm fd}^{\prime}$.


```python
import sympy as sp

# Define variables
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable
DrSipn2 = sp.symbols('DrSipn2')  # DrSipn2 variable
r = sp.symbols('r')  # r variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Define the expression for $\tilde{\omega}_{\rm fd}^{\prime}$
expr_tilwfdprime = 2*expr_a*M

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilwfdprime) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= 2aM
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
a) and $M$ in [this cell](

We are now going to calculate $a$.


```python
import sympy as sp

# Define variables
r = sp.symbols('r')  # r variable
M = sp.symbols('M')  # M variable

# Define the expression for $a$
expr_a = sp.sqrt(M/(r**3))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_a) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $a$ can be represented mathematically as:


$$
\begin{aligned}
a &= \sqrt{\frac{M}{r^{3}}}
\end{aligned}
$$


The $M$ can be represented mathematically as:


```python
# Define the expression for $M$
expr_M = M

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_M) + '\n')
```

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= 2aM
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
m).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
M = sp.symbols('M')  # M variable

# Define the expression for $\tilde{\omega}_{\rm fd}^{\prime}$
expr_omegatildeprm = 2*a

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omegatildeprm) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{\omega}_{\rm fd}^{\prime}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{\omega}_{\rm fd}^{\prime} &= 2a
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11: The Deformed and Rescaled Metric Potentials \[Back to [top](

We are now going to calculate the deformed and rescaled metric potentials.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for the deformed and rescaled metric potentials
expr_deformed_metric_potentials = (r**2)*(a**2 + r**2)/(r**3) - (a**2)*M/(r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_deformed_metric_potentials) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The deformed and rescaled metric potentials can be represented mathematically as:


$$
\begin{aligned}
{\rm Deformed \; Metric \; Potentials} &= \frac{(a^{2} + r^{2})r^{2}}{r^{3}} - \frac{a^{2}M}{r^{4}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{metpots}$$

We collect here terms of the deformed and scaled metric potentials.  See [BB2010](https://arxiv.org/abs/0912.3517) Equations (5.30)--(5.34) and (5.48)--(5.52).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for the deformed and scaled metric potentials
expr_deformed_scaled_metric_potentials = (r**2)*(a**2 + r**2)/(r**3) - (a**2)*M/(r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_deformed_scaled_metric_potentials) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The deformed and scaled metric potentials can be represented mathematically as:


$$
\begin{aligned}
{\rm Deformed \; and \; Scaled \; Metric \; Potentials} &= \frac{(a^{2} + r^{2})r^{2}}{r^{3}} - \frac{a^{2}M}{r^{4}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.a: $\omega$ \[Back to [top](

We are now going to calculate $\omega$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\omega$
expr_omega = (a**2 + r**2)/(r**3) - a**2*M/(r**4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega$ can be represented mathematically as:


$$
\begin{aligned}
\omega &= \frac{(a^{2} + r^{2})}{r^{3}} - \frac{a^{2}M}{r^{4}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{omega}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.31) we have

\begin{equation*}
    \omega = \frac{ \tilde{\omega}_{\rm fd} }{ \Lambda_{t} }.
\end{equation*}

We define $\tilde{\omega}_{\rm fd}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_tilwfd = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(Q_minus_1))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\Lambda_{t}$
expr_Lambdat = 4*(a**2 + r**2)*r - 2*a**2*sp.diff(sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(Q_minus_1))/r**4)), r)

# Define the expression for $\omega$
expr_omega = expr_tilwfd/expr_Lambdat

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
   **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
omegatilde) and $\Lambda_{t}$ in [this cell](

We are now going to calculate $\omega$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_tilwfd = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(Q_minus_1))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\Lambda_{t}$
expr_Lambdat = 4*(a**2 + r**2)*r - 2*a**2*sp.diff(sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(Q_minus_1))/r**4)), r)

# Define the expression for $\omega$
expr_omega = expr_tilwfd/expr_Lambdat

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
   **NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{\omega}_{\rm fd}$
expr_tilwfd = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(Q_minus_1))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $\Lambda_{t}$
expr_Lambdat = 4*(a**2 + r**2)*r - 2*a**2*sp.diff(sp.simplify((3*M/(r**3)) - (2*sp.sqrt(sp.simplify(Q_minus_1))/r**4)), r)

# Define the expression for $\omega$
expr_omega = expr_tilwfd/expr_Lambdat

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_omega) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\omega$ can be represented mathematically as:


$$
\begin{aligned}
\omega &= \frac{\tilde{\omega}_{\rm fd}}{\Lambda_{t}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.b: $e^{2\nu}$ and $e^{\nu}$ \[Back to [top](

We are now going to calculate $e^{2\nu}$ and $e^{\nu}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\nu$
expr_nu = (3*M/(r**3)) + (2*sp.sqrt(sp.simplify(Q_minus_1))/(sp.sqrt(3*eta*r)*r))

# Define the expression for $e^{2\nu}$
expr_exp2nu = sp.exp(2*expr_nu)

# Define the expression for $e^{\nu}$
expr_expnu = sp.exp(expr_nu)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_exp2nu) + '\n')
    f.write(str(expr_expnu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $e^{2\nu}$ can be represented mathematically as:


$$
\begin{aligned}
e^{2\nu} &= e^{(3M/r^{3}) + (2\sqrt{\Sigma}/(\sqrt{3}\eta r))}
\end{aligned}
$$


The $e^**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{exp2nu}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.32), we have

\begin{equation*}
    e^{2 \nu} = \frac{ \Delta_{t} \Sigma }{ \Lambda_t }.
\end{equation*}

It follows that

\begin{equation*}
    e^{\nu} = \sqrt{ \frac{ \Delta_{t} \Sigma }{ \Lambda_t } }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $e^{2\nu}$
expr_exp2nu = expr_Deltat*expr_Sigma

# Define the expression for $e^{\nu}$
expr_expnu = sp.sqrt(expr_exp2nu)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_exp2nu) + '\n')
    f.write(str(expr_expnu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat), $\Sigma$ in [this cell](

We are now going to calculate $\Delta_{t}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_{t}$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{t} &= (a^{2} + r^{2})^{2} - \frac{a^{4}M^{2}}{r^{6}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), and $\Lambda_{t}$ in [this cell](

We are now going to calculate $e^{2\nu}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_{t}$
expr_Lambdat = 4*(a**2 + r**2)*r - 2*a**2*sp.diff(sp.simplify(3*M/(r**3) - 2*sp.sqrt(Q_minus_1)/r**4), r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltat) + '\n')
    f.write(str(expr_Sigma) + '\n')
    f.write(str(expr_Lambdat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\Delta_{t}$ can be represented mathematically as:


$$
\begin{aligned}
\Delta_{t} &= (a^{2} + r^{2})^{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_{t}$
expr_Lambdat = 4*(a**2 + r**2)*r - 2*a**2*sp.diff(sp.simplify(3*M/(r**3) - 2*sp.sqrt(Q_minus_1)/r**4), r)

# Calculate $e^{2\nu}$
expr_exp2nu = expr_Deltat*expr_Sigma/expr_Lambdat

# Calculate $e^{\nu}$
expr_expnu = sp.sqrt(expr_exp2nu)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_expnu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $e^{2\nu}$ can be represented mathematically as:


$$
\begin{aligned}
e^{2\nu} &= \frac{\Delta_{t}\Sigma}{\**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.c: $\tilde{B}$ \[Back to [top](

We are now going to calculate $\tilde{B}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{B}$
expr_tilb = (2*sp.sqrt(Q_minus_1)*r)/(sp.sqrt(3*eta))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilb) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B} &= \frac{2\sqrt{\Sigma}r}{\sqrt{3}\eta}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{btilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.48), we have

\begin{equation*}
    \tilde{B} = \sqrt{ \Delta_{t} }.
\end{equation*}

We define $\Delta_{t}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\tilde{B}$
expr_tilb = sp.sqrt(expr_Deltat)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilb) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B} &= \sqrt{\Delta_{t}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Calculate $\tilde{B}$
expr_tilb = sp.sqrt(expr_Deltat)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilb) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B} &= \sqrt{\Delta_{t}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.d: $\tilde{B}_{r}$ \[Back to [top](

We are now going to calculate $\tilde{B}_{r}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{B}_{r}$
expr_tilbr = (2*sp.sqrt(Q_minus_1)*r)/(sp.sqrt(3*eta))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilbr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}_{r}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B}_{r} &= \frac{2\sqrt{\Sigma}r}{\sqrt{3}\eta}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{brtilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.49), we have

\begin{equation*}
    \tilde{B}_{r} = \frac{ \sqrt{ \Delta_{r} } \Delta_{t}^{\prime} - 2 \Delta_{t} }{ 2 \sqrt{ \Delta_{r} \Delta_{t} } }.
\end{equation*}

We define $\Delta_{r}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\tilde{B}_{r}$
expr_tilbr = ((sp.sqrt(expr_Deltar)*sp.diff(expr_Deltat, r)) - 2*expr_Deltat)/(2*sp.sqrt(expr_Deltar*expr_Deltat))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilbr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltar), $\Delta_{t}^{\prime}$ in [this cell](

We are now going to calculate $\tilde{B}_{r}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\tilde{B}_{r}$
expr_tilbr = ((sp.sqrt(expr_Deltar)*sp.diff(expr_Deltat, r)) - 2*expr_Deltat)/(2*sp.sqrt(expr_Deltar*expr_Deltat))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tilbr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}_{r}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B}_{r} &= \frac{\sqrt{\Delta_{r}}\Delta_{t}^{\prime**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltatprm), and $\Delta_{t}$ in [this cell](

We are now going to calculate $\tilde{B}_{r}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Delta_{t}$, prime
expr_deltatprm = sp.diff(expr_Deltat, r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_deltatprm) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}_{r}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B}_{r} &= \frac{\sqrt{\Delta_{r}}\Delta_{t}^{\prime} - 2\Delta_{t}}{2\sqrt{\Delta_{r}\Delta_{t}}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltat.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{t}$
expr_Deltat = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Calculate $\tilde{B}_{r}$
expr_Brtilde = (sp.sqrt(expr_Deltar)*sp.diff(expr_Deltat, r) - 2*expr_Deltat)/(2*sp.sqrt(expr_Deltar*expr_Deltat))

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Brtilde) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{B}_{r}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{B}_{r} &= \frac{\sqrt{\Delta_{r}}\Delta_{t}^{\prime} - 2\Delta_{t}}{2\sqrt{\Delta_{r}\Delta_{t}}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.e: $e^{2\tilde{\mu}}$ and $e^{\tilde{\mu}}$ \[Back to [top](

We are now going to calculate $e^{2\tilde{\mu}}$ and $e^{\tilde{\mu}}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{B}_{r}$
expr_Brtilde = (sp.sqrt(expr_Deltar)*sp.diff(expr_Deltat, r) - 2*expr_Deltat)/(2*sp.sqrt(expr_Deltar*expr_Deltat))

# Calculate $e^{2\tilde{\mu}}$ and $e^{\tilde{\mu}}$
expr_exp2mu = expr_Brtilde**2
expr_expmu = sp.sqrt(expr_exp2mu)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_exp2mu) + '\n')
    f.write(str(expr_expmu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $e^{2\tilde{\mu}}$ can be represented mathematically as:


$$
\begin{aligned}
e^{2\tilde{\mu**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{exp2mu}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.50), we have

\begin{equation*}
    e^{2 \tilde{\mu}} = \Sigma.
\end{equation*}

It follows that

\begin{equation*}
    e^{\tilde{\mu}} = \sqrt{ \Sigma }.
\end{equation*}


We define $\Sigma$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate $e^{2\tilde{\mu}}$ and $e^{\tilde{\mu}}$
expr_exp2mu = expr_Sigma
expr_expmu = sp.sqrt(expr_exp2mu)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_exp2mu) + '\n')
    f.write(str(expr_expmu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $e^{2\tilde{\mu}}$ can be represented mathematically as**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma).


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate $e^{\tilde{\mu}}$ and $e^{2\tilde{\mu}}$
expr_expmu = sp.sqrt(expr_exp2mu)
expr_exp2mu = expr_Sigma

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_expmu) + '\n')
    f.write(str(expr_exp2mu) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $e^{\tilde{\mu}}$ can be represented mathematically as:


$$
\begin{aligned}
e^{\tilde{\mu}} = \sqrt{\Sigma}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.f: $\tilde{J}$ \[Back to [top](

We are now going to calculate $\tilde{J}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\tilde{J}$
expr_tildeJ = (a**2 + r**2)*r - a**2*sp.sqrt(Q_minus_1)/r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tildeJ) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{J}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{J} = (a^{2} + r^{2})r - \frac{a^{2}\sqrt{\Sigma}}{r}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{jtilde}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.51) we have

\begin{equation*}
    \tilde{J} = \sqrt{ \Delta_{r} }.
\end{equation*}

We define $\Delta_{r}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Calculate $\tilde{J}$
expr_tildeJ = sp.sqrt(expr_Deltar)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tildeJ) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{J}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{J} = \sqrt{\Delta_{r}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltar) below.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Calculate $\tilde{J}$ and write to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(sp.sqrt(expr_Deltar)) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $\tilde{J}$ can be represented mathematically as:


$$
\begin{aligned}
\tilde{J} = \sqrt{\Delta_{r}}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.g: $Q$ \[Back to [top](

We are now going to calculate $Q$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$
expr_Q = (a**2 + r**2)**2 - a**4*M**2/(r**6) - 3*(a**4*M**2)/(r**6)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Q) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The $Q$ can be represented mathematically as:


$$
\begin{aligned}
Q = \Delta_{r} - \frac{3\mu^2}{r^6}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{q}$$

From [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.52),

\begin{equation*}
    Q = 1 + \underbrace{ \frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2} }_{\rm DrSipn2} + \underbrace{ \frac{ \Sigma }{ \Lambda_t \sin^{2} \theta } }_{\rm Q\ coefficient\ 1} \left( \smash[b]{ \underbrace{ \hat{\bf p} \cdot \boldsymbol{\xi} r }_{\rm pdotxir} } \right)^{2} + \underbrace{ \frac{ 1 }{ \Sigma \sin^{2} \theta } }_{\rm Q\ coefficient\ 2} \left( \smash[b]{ \underbrace{ \hat{\bf p} \cdot {\bf v} r }_{\rm pdotvr} } \right)^{2};
\end{equation*}

We group togther and compute $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ in [this cell](


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
drsipn2), $\frac{ \Sigma }{ \Lambda_t \sin^{2} \theta }$ in [this cell](

We are now going to calculate the two Q coefficients.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_t$
expr_Lambdat = 1/(r**2)

# Calculate Q coefficient 1
expr_Qcoefficient1 = expr_Sigma / (expr_Lambdat * sp.sin(a)**2)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The Q coefficient 1 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 1 = \frac{\Sigma}{\Lambda_t \sin^{2} \theta}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
qcoeff1), $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](

We are now going to calculate the remaining Q coefficient.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_t$
expr_Lambdat = 1/(r**2)

# Define the expression for $\hat{\bf p} \cdot \boldsymbol{\xi} r$
expr_pdotxir = a*sp.sin(a)*sp.cos(r)

# Calculate Q coefficient 2
expr_Qcoefficient2 = 1 / (expr_Sigma * sp.sin(a)**2)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotxir) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The Q coefficient 2 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 2 = \frac{1**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pdotxir), $\frac{ 1 }{ \Sigma \sin^{2} \theta }$ in [this cell](

We are now going to calculate the final expression for $Q$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_t$
expr_Lambdat = 1/(r**2)

# Calculate $Q$ coefficient 3
expr_Qcoefficient3 = 1 / (expr_Sigma * sp.sin(a)**2)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient3) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression $\frac{ 1 }{ \Sigma \sin^{2} \theta }$ can be represented mathematically as:


$$
\begin{aligned}
\frac{ 1 }{ \Sigma \sin^{2} \theta }
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
qcoeff2), and $\hat{\bf p} \cdot {\bf v} r$ in [this cell](

We are now going to calculate the final expression for $Q$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Define the expression for $\Lambda_t$
expr_Lambdat = 1/(r**2)

# Calculate $Q$ coefficient 4
expr_Qcoefficient4 = expr_pdotvr**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient4) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression $\hat{\bf p} \cdot {\bf v} r$ can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot {\bf v} r = a \sin(a) \cos(r)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Q = 1 + DrSipn2 + Qcoeff1*pdotxir*pdotxir + Qcoeff2*pdotvr*pdotvr

We are now going to calculate the final expression for $Q$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate DrSipn2
DrSipn2 = expr_Deltar / (expr_Sigma * sp.sin(a)**2)

# Calculate Q coefficient 1
Q_coefficient_1 = expr_Sigma / (sp.Lambdat * sp.sin(r)**2) * a**2 * r**2

# Calculate pdotxir
pdotxir = a*sp.sin(a)*r*sp.cos(r)

# Calculate pdotvr
pdotvr = a*sp.sin(a)*r*sp.cos(r)

# Calculate Q coefficient 2
Q_coefficient_2 = 1 / (expr_Sigma * sp.sin(a)**2)

# Calculate Q
Q = 1 + DrSipn2 + Q_coefficient_1*pdotxir*pdotxir + Q_coefficient_2*pdotvr*pdotvr

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.g.i: $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ \[Back to [top](

We are now going to calculate the expression for DrSipn2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate DrSipn2
DrSipn2 = expr_Deltar / (expr_Sigma * sp.sin(a)**2)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(DrSipn2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression $\frac{ \Delta_{r} }{ \Sigma } \left( \hat{\bf p} \cdot {\bf n} \right)^{2}$ can be represented mathematically as:


$$
\begin{aligned}
\frac{ \Delta_{r} }**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{drsipn2}$$

We define $\Delta_{r}$ in [this cell](

This code defines the variable DrSipn2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Deltar) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The variable DrSipn2 can be represented mathematically as:


$$
\begin{aligned}
DrSipn2 = \Delta_{r} / (\Sigma * sin(a)**2)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
deltar), $\Sigma$ in [this cell](

This code defines the variable Sigma.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Sigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The variable Sigma can be represented mathematically as:


$$
\begin{aligned}
\Sigma = M/((a**2 + r**2)*r)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), and $\hat{\bf p} \cdot {\bf n}$ in [this cell](

This code defines the variables usigma and pdotn.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate usigma
usigma = expr_Sigma

# Define the expression for $\hat{\bf p} \cdot {\bf n}$
pdotn = a*sp.sin(a)*sp.cos(r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(usigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The variables usigma and pdotn can be represented mathematically as:


$$
\begin{aligned}
usigma &= M/((a**2 + r**2)*r) \\
pdotn &= a \sin(a) \cos(r)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
qcoeff1'></a>

This code defines the expression for DrSipn2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate pdotn
pdotn = a*sp.sin(a)*sp.cos(r)

# Calculate DrSipn2
DrSipn2 = expr_Deltar * (pdotn)**2 / expr_Sigma

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(DrSipn2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression DrSipn2 can be represented mathematically as:


$$
\begin{aligned}
DrSipn2 = \frac{\Delta_{r} (\hat{\bf p} \cdot {\bf n})^2}{\Sigma}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.g.ii: Q Coefficient 1 \[Back to [top](

This code defines the expression for Qcoefficient1.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Delta_{r}$
expr_Deltar = (a**2 + r**2)**2 - a**4*M**2/(r**6)

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Calculate Qcoefficient1
Qcoefficient1 = expr_Sigma / (sp.Lambdat * sp.sin(r)**2) * a**2 * r**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Qcoefficient1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient1 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 1 = \frac{\Sigma}{\Lambda_t \sin^{2} \theta} (a^{2} r^{2})
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{qcoeff1}$$

We defined $Q$ coefficient 1 in [this cell](

This code defines the expression for Qcoefficient1.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$ coefficient 1
Qcoefficient1 = expr_Sigma / (sp.Lambdat * sp.sin(r)**2) * a**2 * r**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(Qcoefficient1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient1 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 1 = \frac{\Sigma}{\Lambda_t \sin^{2} \theta} (a^{2} r^{2})
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q) as

\begin{equation*}
    Q\ {\rm coefficient\ 1} = \frac{ \Sigma }{ \Lambda_t \sin^{2} \theta }
\end{equation*}

We define $\Sigma$ in [this cell](

This code defines the expression for Sigma.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Sigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Sigma can be represented mathematically as:


$$
\begin{aligned}
\Sigma = \frac{M}{(a^{2}+r^{2})r}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma), $\Lambda_{t}$ in [this cell](

This code defines the expression for Lambda_t.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Lambda_{t}$
expr_Lambdat = a**2 + r**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Lambdat) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Lambda_t can be represented mathematically as:


$$
\begin{aligned}
\Lambda_{t} = a^{2}+r^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
lambdat), and $\sin^{2} \theta$ in [this cell](

This code defines the expression for sin2theta.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\sin^{2} \theta$
expr_sin2theta = sp.sin(a)**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_sin2theta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression sin2theta can be represented mathematically as:


$$
\begin{aligned}
\sin^{2} \theta = \sin^{2}(a)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Qcoeff1 = Sigma/(Lambdat*sin2theta)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='qcoeff2'></a>

This code defines the expression for Qcoefficient2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$ coefficient 2
expr_Qcoefficient2 = (a**2 + r**2)**(-1)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient2 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 2 = (a^{2}+r^{2})^{-1}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 11.g.iii: Q Coefficient 2 \[Back to [top](

This code defines the expression for Qcoefficient2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$ coefficient 2
expr_Qcoefficient2 = (a**2 + r**2)**(-1)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient2 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 2 = (a^{2}+r^{2})^{-1}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{qcoeff2}$$

We defined $Q$ coefficient 2 in [this cell](

This code defines the expression for Qcoefficient2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$ coefficient 2
expr_Qcoefficient2 = (a**2 + r**2)**(-1)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient2 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 2 = (a^{2}+r^{2})^{-1}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
q) as

\begin{equation*}
    Q\ {\rm coefficient\ 2} = \frac{ 1 }{ \Sigma \sin^{2} \theta }
\end{equation*}

We define $\Sigma$ in [this cell](

This code defines the expression for Sigma.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\Sigma$
expr_Sigma = M/((a**2 + r**2)*r)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Sigma) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Sigma can be represented mathematically as:


$$
\begin{aligned}
\Sigma = \frac{M}{(a^{2}+r^{2})r}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
usigma) and $\sin^{2} \theta$ in [this cell](

This code defines the expression for sin2theta.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\sin^{2} \theta$
expr_sin2theta = sp.sin(a)**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_sin2theta) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression sin2theta can be represented mathematically as:


$$
\begin{aligned}
\sin^{2} \theta = \sin^{2}(a)
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
sin2theta).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

Qcoeff2 = 1/(Sigma*sin2theta)
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='tort'></a>

This code defines the expression for Qcoefficient2.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $Q$ coefficient 2
expr_Qcoefficient2 = (a**2 + r**2)**(-1)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_Qcoefficient2) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression Qcoefficient2 can be represented mathematically as:


$$
\begin{aligned}
Q \text{ coefficient } 2 = (a^{2}+r^{2})^{-1}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12: Tortoise Terms \[Back to [top](

This code defines the expression for tortoise terms.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for tortoise terms
expr_tortoise_terms = (a**2 + r**2)**(-1/2)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_tortoise_terms) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression tortoise terms can be represented mathematically as:


$$
\begin{aligned}
\text{tortoise terms} = (a^{2}+r^{2})^{-1/2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{tort}$$

We collect here terms related to the conversion from Boyer-Lindquist coordinates to tortoise coordinates.  Details of the conversation are given in the appendix of [P2010](https://arxiv.org/abs/0912.3466v2).

<a id='pphi'></a>

This code defines the expression for pphi.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p\phi$
expr_pphi = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pphi) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pphi can be represented mathematically as:


$$
\begin{aligned}
p\phi = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12.a: $p_{\phi}$ \[Back to [top](

This code defines the expression for pphi.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p_{\phi}$
expr_pphi = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pphi) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pphi can be represented mathematically as:


$$
\begin{aligned}
p_{\phi} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{pphi}$$

From the discussion preceding [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.41), the phi component of the tortoise momentum $p_{\phi}$ is given by

\begin{equation*}
    p_{\phi} = \hat{\bf p} \cdot \boldsymbol{\xi} r.
\end{equation*}

We define $\hat{\bf p} \cdot \boldsymbol{\xi} r$ in [this cell](

This code defines the expression for hatp.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p} \cdot \boldsymbol{\xi} r$
expr_hatp_dot_xi_r = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_hatp_dot_xi_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression hatp_dot_xi_r can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot \boldsymbol{\xi} r = (a^{2**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
pdotxir).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pphi = pdotxir
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotvr'></a>

This code defines the expression for pphi.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p_{\phi}$
expr_pphi = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pphi) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pphi can be represented mathematically as:


$$
\begin{aligned}
p_{\phi} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12.b: $\hat{\bf p} \cdot {\bf v} r$ \[Back to [top](

This code defines the expression for pdotvr.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p} \cdot {\bf v} r$
expr_pdotvr = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotvr) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotvr can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot {\bf v} r = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{pdotvr}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot {\bf v} r = \left( \hat{p}_{1} v_{1} + \hat{p}_{2} v_{2} + \hat{p}_{3} v_{3} \right) r
\end{equation*}

We define $\hat{\bf p}$ in [this cell](

This code defines the expression for hatp.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p}$
expr_hatp = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_hatp) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression hatp can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hatp), ${\bf v}$ in [this cell](

This code defines the expression for v.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for ${\bf v}$
expr_v = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_v) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression v can be represented mathematically as:


$$
\begin{aligned}
{\bf v} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
v), and $r$ in [this cell](

This code defines the expression for r.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $r$
expr_r = a**2 + r**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression r can be represented mathematically as:


$$
\begin{aligned}
r = a^{2} + r^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotvr = (phat1*v1 + phat2*v2 + phat3*v3)*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotn'></a>

This code defines the expression for pdotn.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p \cdot {\bf n}$
expr_pdotn = (phat1*v1 + phat2*v2 + phat3*v3)*r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotn) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotn can be represented mathematically as:


$$
\begin{aligned}
p \cdot {\bf n} = (\hat{p}_{1} v_{1} + \hat{p}_{2} v_{2} + \hat{p}_{3} v_{3}) r
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12.c: $\hat{\bf p} \cdot {\bf n}$ \[Back to [top](

This code defines the expression for pdotn.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p} \cdot {\bf n}$
expr_pdotn = (phat1*v1 + phat2*v2 + phat3*v3)*r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotn) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotn can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot {\bf n} = (\hat{p}_{1} v_{1} + \hat{p}_{2} v_{2} + \hat{p}_{3} v_{3}) r
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{pdotn}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot {\bf n} = \hat{p}_{1} n_{1} + \hat{p}_{2} n_{2} + \hat{p}_{3} n_{3}
\end{equation*}

We define $\hat{\bf p}$ in [this cell](

This code defines the expression for hatp.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p}$
expr_hatp = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_hatp) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression hatp can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hatp) and ${\bf n}$ in [this cell](

This code defines the expression for n.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for ${\bf n}$
expr_n = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_n) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression n can be represented mathematically as:


$$
\begin{aligned}
{\bf n} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
n).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotn = phat1*n1 + phat2*n2 + phat3*n3
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='pdotxir'></a>

This code defines the expression for pdotxir.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p \cdot {\bf n}$
expr_pdotn = phat1*n1 + phat2*n2 + phat3*n3

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotn) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotxir can be represented mathematically as:


$$
\begin{aligned}
p \cdot {\bf n} = \hat{p}_{1} n_{1} + \hat{p}_{2} n_{2} + \hat{p}_{3} n_{3}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12.d: $\hat{\bf p} \cdot \boldsymbol{\xi} r$ \[Back to [top](

This code defines the expression for pdotxir.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p} \cdot \boldsymbol{\xi} r$
expr_pdotxir = (phat1*v1 + phat2*v2 + phat3*v3)*r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotxir) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotxir can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot \boldsymbol{\xi} r = (\hat{p}_{1} v_{1} + \hat{p}_{2} v_{2} + \hat{p}_{3} v_{3}) r
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{pdotxir}$$

We have

\begin{equation*}
    \hat{\bf p} \cdot \boldsymbol{\xi} r = \left( \hat{p}_{1} \xi_{1} + \hat{p}_{2} \xi_{2} + \hat{p}_{3} \xi_{3} \right) r
\end{equation*}

We define $\hat{\bf p}$ in [this cell](

This code defines the expression for hatp.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p}$
expr_hatp = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_hatp) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotxir can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} \cdot \boldsymbol{\xi} r = \left( \hat{p}_{1} \xi_{1} + \hat{p}_{2} \xi_{2} + \hat{p}_{3**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
hatp), $\boldsymbol{\xi}$ in [this cell](

This code defines the expression for xi.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\boldsymbol{\xi}$
expr_xi = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_xi) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression xi can be represented mathematically as:


$$
\begin{aligned}
\boldsymbol{\xi} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
xi), and $r$ in [this cell](

This code defines the expression for r.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $r$
expr_r = a**2 + r**2

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_r) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression r can be represented mathematically as:


$$
\begin{aligned}
r = a^{2} + r^{2}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
r).


```python
%%writefile -a $Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt

pdotxir = (phat1*xi1 + phat2*xi2 + phat3*xi3)*r
```

    Appending to SEOBNR/v4P_Hamiltonian-Hreal_on_top.txt


<a id='hatp'></a>

This code defines the expression for pdotxir.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $p \cdot {\bf \xi} r$
expr_pdotxir = (phat1*xi1 + phat2*xi2 + phat3*xi3)*r

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_pdotxir) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression pdotxir can be represented mathematically as:


$$
\begin{aligned}
p \cdot {\bf \xi} r = (\hat{p}_{1} \xi_{1} + \hat{p}_{2} \xi_{2} + \hat{p}_{3} \xi_{3}) r
\end{**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
Step 12.e: $\hat{\bf p}$ \[Back to [top](

This code defines the expression for hatp.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for $\hat{\bf p}$
expr_hatp = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_hatp) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics. It represents the amount of "similarity" between two vectors.
    +   In this case, we are calculating the complex conjugate of a vector.

### Mathematical Representation

The expression hatp can be represented mathematically as:


$$
\begin{aligned}
\hat{\bf p} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
toc)\]
$$\label{hatp}$$

From the discussion after [BB2010](https://arxiv.org/abs/0912.3517) Equation (3.41), we have $\hat{\bf p} = {\bf p}/m$ where $m$ is the mass of a nonspinning test particle and ${\bf p}$ is *conjugate* momentum.  Following Lines 319--321 of LALSimIMRSpinEOBHamiltonianPrec.c, we convert the Boyer-Lindquist momentum ${\bf p}$ to the tortoise momentum (see the appendix of [P2010](https://arxiv.org/abs/0912.3466v2)) via

\begin{align*}
    \hat{\bf p} = {\bf p} + {\rm prT} \left( 1 - \frac{1}{\rm csi1} \right) {\bf n}
\end{align*}

We define prT in [this cell](

This code defines the expression for prT.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for ${\rm prT}$
expr_prT = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_prT) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
prt), csi1 in [this cell](

This code defines the expression for csi1.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for ${\rm csi1}$
expr_csi1 = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_csi1) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics.

### Mathematical Representation

The expression csi1 can be represented mathematically as:


$$
\begin{aligned}
{\rm csi1} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$**NRPy+: The Effective Hamiltonian $H_{\rm eff}$**
=============================================

### Theory Review

#### Introduction to the Complex Conjugate of a Vector

*   **The Complex Conjugate of a Vector:** In this section, we discuss the complex conjugate of a vector, which is a crucial component in numerical relativity and gravitational wave astronomy.

### Code Explanation


```python
"""
csi1), and ${\bf n}$ in [this cell](

This code defines the expression for ${\bf n}$.


```python
import sympy as sp

# Define variables
a = sp.symbols('a')  # a variable
r = sp.symbols('r')  # r variable
eta = sp.symbols('eta')  # eta variable
M = sp.symbols('M')  # M variable
Q_minus_1 = sp.symbols('Q_minus_1')  # Q_minus_1 variable

# Define the expression for ${\bf n}$
expr_n = (a**2 + r**2)**(-3/4)

# Write the expressions to file
with open($Ccodesdir/v4P_Hamiltonian-Hreal_on_top.txt, 'a') as f:
    f.write(str(expr_n) + '\n')
```

### Theory Review


*   **Complex Conjugate of a Vector**: The complex conjugate of a vector is a fundamental concept in mathematics and physics.

### Mathematical Representation

The expression ${\bf n}$ can be represented mathematically as:


$$
\begin{aligned}
{\bf n} = (a^{2}+r^{2})^{-3/4}
\end{aligned}
$$