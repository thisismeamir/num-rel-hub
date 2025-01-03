<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# NRPy+'s Reference Metric Interface

## Author: Zach Etienne
### Formatting improvements courtesy Brandon Clark

### NRPy+ Source Code for this module: [reference_metric.py](../edit/reference_metric.py)

## This notebook expounds on NRPy+'s Reference Metric Interface, underlining its efficacy in reducing computational complexity when modeling geometrically specific problems. By utilizing appropriate coordinate systems—Spherical, Cylindrical, Cartesian, or Prolate Spheroidal—geometric entities are optimally defined, thereby mitigating the 'Curse of Dimensionality'.

## Introduction:
### Why use a reference metric? Benefits of choosing the best coordinate system for the problem

When solving a partial differential equation on the computer, it is useful to first pick a coordinate system well-suited to the geometry of the problem. For example, if we are modeling a spherically-symmetric star, it would be hugely wasteful to model the star in 3-dimensional Cartesian coordinates ($x$,$y$,$z$). This is because, in Cartesian coordinates, we would need to choose high sampling in all three Cartesian directions. If instead, we chose to model the star in spherical coordinates  ($r$,$\theta$,$\phi$), so long as the star is centered at $r=0$, we would not need to model the star with more than one point in the $\theta$ and $\phi$ directions!

A similar argument holds for stars that are *nearly* spherically symmetric. Such stars may exhibit density distributions that vary slowly in $\theta$ and $\phi$ directions (e.g., isolated neutron stars or black holes). In these cases, the number of points needed to sample the angular directions will still be much smaller than in the radial direction.

Thus the choice of an appropriate reference metric may directly mitigate the [Curse of Dimensionality](https://en.wikipedia.org/wiki/Curse_of_dimensionality).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follow

1. [Step 1](#define_ref_metric): Defining a reference metric, [`reference_metric.py`](../edit/reference_metric.py)
1. [Step 2](#define_geometric): Defining geometric quantities, **`ref_metric__hatted_quantities()`**
1. [Step 3](#prescribed_ref_metric): Prescribed reference metrics in [`reference_metric.py`](../edit/reference_metric.py)
    1. [Step 3.a](#sphericallike): Spherical-like coordinate systems
        1. [Step 3.a.i](#spherical): **`reference_metric::CoordSystem = "Spherical"`**
        1. [Step 3.a.ii](#sinhspherical): **`reference_metric::CoordSystem = "SinhSpherical"`**
        1. [Step 3.a.iii](#sinhsphericalv2): **`reference_metric::CoordSystem = "SinhSphericalv2"`**
    1. [Step 3.b](#cylindricallike): Cylindrical-like coordinate systems
        1. [Step 3.b.i](#cylindrical): **`reference_metric::CoordSystem = "Cylindrical"`**
        1. [Step 3.b.ii](#sinhcylindrical): **`reference_metric::CoordSystem = "SinhCylindrical"`**
        1. [Step 3.b.iii](#sinhcylindricalv2): **`reference_metric::CoordSystem = "SinhCylindricalv2"`**
    1. [Step 3.c](#cartesianlike): Cartesian-like coordinate systems
        1. [Step 3.c.i](#cartesian): **`reference_metric::CoordSystem = "Cartesian"`**
        1. [Step 3.c.ii](#sinhcartesian): **`reference_metric::CoordSystem = "SinhCartesian"`**
    1. [Step 3.d](#prolatespheroidal): Prolate spheroidal coordinates
        1. [Step 3.d.i](#symtp): **`reference_metric::CoordSystem = "SymTP"`**
        1. [Step 3.d.ii](#sinhsymtp): **`reference_metric::CoordSystem = "SinhSymTP"`**
1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='define_ref_metric'></a>

# Step 1: Defining a reference metric, [`reference_metric.py`](../edit/reference_metric.py) \[Back to [top](#toc)\]
$$\label{define_ref_metric}$$

***Note that currently only orthogonal reference metrics of dimension 3 or fewer are supported. This can be extended if desired.***

NRPy+ assumes all curvilinear coordinate systems map directly from a uniform, Cartesian numerical grid with coordinates $(x,y,z)$=(`xx[0]`,`xx[1]`,`xx[2]`). Thus, when defining reference metrics, all defined coordinate quantities must be in terms of the `xx[]` array. As we will see, this adds a great deal of flexibility

For example,  [**reference_metric.py**](../edit/reference_metric.py) requires that the *orthogonal coordinate scale factors* be defined. As described [here](https://en.wikipedia.org/wiki/Curvilinear_coordinates), the $i$th scale factor is the positive root of the metric element $g_{ii}$. In ordinary spherical coordinates $(r,\theta,\phi)$, with line element $ds^2 = g_{ij} dx^i dx^j = dr^2+ r^2 d \theta^2 + r^2 \sin^2\theta \ d\phi^2$, we would first define
* $r = xx_0$
* $\theta = xx_1$
* $\phi = xx_2$,

so that the scale factors are defined as
* `scalefactor_orthog[0]` = $1$
* `scalefactor_orthog[1]` = $r$
* `scalefactor_orthog[2]` = $r \sin \theta$.

Here is the corresponding code:


```python
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support

r = rfm.xx[0]
th = rfm.xx[1]
ph = rfm.xx[2]

rfm.scalefactor_orthog[0] = 1
rfm.scalefactor_orthog[1] = r
rfm.scalefactor_orthog[2] = r*sp.sin(th)

# Notice that the scale factor will be given
#    in terms of the fundamental Cartesian
#    grid variables, and not {r,th,ph}:
print("r*sin(th) = "+str(rfm.scalefactor_orthog[2]))
```

    r*sin(th) = xx0*sin(xx1)


Next suppose we wish to modify our radial coordinate $r(xx_0)$ to be an exponentially increasing function so that our numerical grid $(xx_0,xx_1,xx_2)$ will map to a spherical grid with radial grid spacing ($\Delta r$) that *increases* with $r$. Generally, we will find it useful to define $r(xx_0)$ to be an odd function, so let's choose

$$r(xx_0) = a \sinh(xx_0/s),$$

where $a$ is an overall radial scaling factor, and $s$ denotes the scale (in units of $xx_0$) over which exponential growth will take place. In our implementation below, note that we use the relation

$$\sinh(x) = \frac{e^x - e^{-x}}{2},$$

as SymPy finds it easier to evaluate exponentials than hyperbolic trigonometric functions.


```python
a,s = sp.symbols('a s',positive=True)
xx0_rescaled = rfm.xx[0] / s
r = a*(sp.exp(xx0_rescaled) - sp.exp(-xx0_rescaled))/2

# Must redefine the scalefactors since 'r' has been updated!
rfm.scalefactor_orthog[0] = 1
rfm.scalefactor_orthog[1] = r
rfm.scalefactor_orthog[2] = r*sp.sin(th)

print(rfm.scalefactor_orthog[2])
```

    a*(exp(xx0/s) - exp(-xx0/s))*sin(xx1)/2


Often we will find it useful to also define the appropriate mappings from (`xx[0]`,`xx[1]`,`xx[2]`) to Cartesian coordinates (for plotting purposes) and ordinary spherical coordinates (e.g., in case of initial data when solving a PDE are naturally written in spherical coordinates). For this purpose, `reference_metric.py` also declares lists **`xx_to_Cart[]`** and **`xxSph[]`**, which in this case are defined as


```python
rfm.xxSph[0] = r
rfm.xxSph[1] = th
rfm.xxSph[2] = ph

rfm.xx_to_Cart[0] = r*sp.sin(th)*sp.cos(ph)
rfm.xx_to_Cart[1] = r*sp.sin(th)*sp.sin(ph)
rfm.xx_to_Cart[2] = r*sp.cos(th)

# Here we show off SymPy's pretty_print()
#   and simplify() functions. Nice, no?
sp.pretty_print(sp.simplify(rfm.xx_to_Cart[0]))
```

                            ⎛xx₀⎞
    a⋅sin(xx₁)⋅cos(xx₂)⋅sinh⎜───⎟
                            ⎝ s ⎠


<a id='define_geometric'></a>

# Step 2: Define geometric quantities, `ref_metric__hatted_quantities()` \[Back to [top](#toc)\]
$$\label{define_geometric}$$

Once `scalefactor_orthog[]` has been defined, the function **`ref_metric__hatted_quantities()`** within [reference_metric.py](../edit/reference_metric.py) can be called to define a number of geometric quantities useful for solving PDEs in curvilinear coordinate systems. 

Adopting the notation of [Baumgarte, Montero, Cordero-Carrión, and Müller, PRD 87, 044026 (2012)](https://arxiv.org/abs/1211.6632), geometric quantities related to the reference metric are named "hatted" quantities. For example, the reference metric is defined as $\hat{g}_{ij}$=`ghatDD[i][j]`:


```python
rfm.ref_metric__hatted_quantities()

sp.pretty_print(sp.Matrix(rfm.ghatDD))
```

    ⎡1           0                         0              ⎤
    ⎢                                                     ⎥
    ⎢                     2                               ⎥
    ⎢      ⎛ xx₀    -xx₀ ⎞                                ⎥
    ⎢      ⎜ ───    ─────⎟                                ⎥
    ⎢    2 ⎜  s       s  ⎟                                ⎥
    ⎢   a ⋅⎝ℯ    - ℯ     ⎠                                ⎥
    ⎢0  ───────────────────                0              ⎥
    ⎢            4                                        ⎥
    ⎢                                                     ⎥
    ⎢                                          2          ⎥
    ⎢                           ⎛ xx₀    -xx₀ ⎞           ⎥
    ⎢                           ⎜ ───    ─────⎟           ⎥
    ⎢                         2 ⎜  s       s  ⎟     2     ⎥
    ⎢                        a ⋅⎝ℯ    - ℯ     ⎠ ⋅sin (xx₁)⎥
    ⎢0           0           ─────────────────────────────⎥
    ⎣                                      4              ⎦


In addition to $\hat{g}_{ij}$, **`ref_metric__hatted_quantities()`** also provides the following. 
* The rescaling "matrix" `ReDD[i][j]` is used for separating singular (due to chosen coordinate system) pieces of smooth rank-2 tensor components from the smooth parts, so that the smooth parts can be used within temporal and spatial differential operators.
* The inverse reference metric: $\hat{g}^{ij}$=`ghatUU[i][j]`.
* The reference metric determinant: $\det\left(\hat{g}_{ij}\right)$=`detgammahat`.
* The first and second derivatives of the reference metric: $\hat{g}_{ij,k}$=`ghatDD_dD[i][j][k]`; $\hat{g}_{ij,kl}$=`ghatDD_dDD[i][j][k][l]`.
* The Christoffel symbols associated with the reference metric, $\hat{\Gamma}^i_{jk}$ = `GammahatUDD[i][j][k]` and their first derivatives $\hat{\Gamma}^i_{jk,l}$ = `GammahatUDD_dD[i][j][k][l]`.

For example, the Christoffel symbol $\hat{\Gamma}^{xx_1}_{xx_2 xx_2}=\hat{\Gamma}^1_{22}$ is given by `GammahatUDD[1][2][2]`:


```python
sp.pretty_print(sp.simplify(rfm.GammahatUDD[1][2][2]))
```

    -sin(2⋅xx₁) 
    ────────────
         2      


Given the trigonometric identity $2\sin(x)\cos(x) = \sin(2x)$, notice that the above expression is equivalent to Eq. 18 of [Baumgarte, Montero, Cordero-Carrión, and Müller, PRD 87, 044026 (2012)](https://arxiv.org/abs/1211.6632). This is expected since the sinh-radial spherical coordinate system is equivalent to ordinary spherical coordinates in the angular components.

<a id='prescribed_ref_metric'></a>

# Step 3: Prescribed reference metrics in [`reference_metric.py`](../edit/reference_metric.py) \[Back to [top](#toc)\]
$$\label{prescribed_ref_metric}$$

One need not manually define scale factors or other quantities for reference metrics, as a number of prescribed reference metrics are already defined in [reference_metric.py](../edit/reference_metric.py). These can be accessed by first setting the parameter **reference_metric::CoordSystem** to one of the following, and then calling the function **`rfm.reference_metric()`**.


```python
import indexedexp as ixp    # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri          # NRPy+: Functions having to do with numerical grids

# Step 0a: Initialize parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "Spherical"))

# Step 0b: Declare global variables
xx = gri.xx
xx_to_Cart = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cart_to_xx = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cartx,Carty,Cartz = sp.symbols("Cartx Carty Cartz", real=True)
Cart = [Cartx,Carty,Cartz]
xxSph  = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
scalefactor_orthog = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
have_already_called_reference_metric_function = False



CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
M_PI,M_SQRT1_2 = par.Cparameters("#define",thismodule,["M_PI","M_SQRT1_2"],"")

UnitVectors = ixp.zerorank2(DIM=3)
```

We will find the following plotting function useful for analyzing coordinate systems in which the radial coordinate is rescaled.


```python
def create_r_of_xx0_plots(CoordSystem, r_of_xx0,rprime_of_xx0):
    import matplotlib.pyplot as plt  # matplotlib: Python module specializing in plotting capabilities
    plt.clf()
    Nr = 20
    dxx0 = 1.0 / float(Nr)
    xx0s    = []
    rs      = []
    deltars = []
    rprimes = []
    for i in range(Nr):
        xx0 = (float(i) + 0.5)*dxx0
        xx0s.append(xx0)
        rs.append(     sp.sympify(str(r_of_xx0     ).replace("xx0",str(xx0))))
        rprimes.append(sp.sympify(str(rprime_of_xx0).replace("xx0",str(xx0))))
        if i>0:
            deltars.append(sp.log(rs[i]-rs[i-1],10))
        else:
            deltars.append(sp.log(2*rs[0],10))

    # fig, ax = plt.subplots()
    fig = plt.figure(figsize=(12,12)) # 8 in x 8 in

    ax = fig.add_subplot(221)
    ax.set_title(r"$r(xx_0)$ for "+CoordSystem,fontsize='x-large')
    ax.set_xlabel(r"$xx_0$",fontsize='x-large')
    ax.set_ylabel(r"$r(xx_0)$",fontsize='x-large')
    ax.plot(xx0s, rs, 'k.', label='Spacing between\nadjacent gridpoints')
    # legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large')
    # legend.get_frame().set_facecolor('C1')

    ax = fig.add_subplot(222)
    ax.set_title('Grid spacing for '+CoordSystem,fontsize='x-large')
    ax.set_xlabel(r"$xx_0$",fontsize='x-large')
    ax.set_ylabel(r"$\log_{10}(\Delta r)$",fontsize='x-large')
    ax.plot(xx0s, deltars, 'k.', label='Spacing between\nadjacent gridpoints\nin $r(xx_0)$ plot')
    legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('C1')

    ax = fig.add_subplot(223)
    ax.set_title(r"$r'(xx_0)$ for "+CoordSystem,fontsize='x-large')
    ax.set_xlabel(r"$xx_0$",fontsize='x-large')
    ax.set_ylabel(r"$r'(xx_0)$",fontsize='x-large')
    ax.plot(xx0s, rprimes, 'k.', label='Nr=96')
    # legend = ax.legend(loc='upper left', shadow=True, fontsize='x-large')
    # legend.get_frame().set_facecolor('C1')

    plt.tight_layout(pad=2)
    plt.show()
```

<a id='sphericallike'></a>

## Step 3.a: Spherical-like coordinate systems \[Back to [top](#toc)\]
$$\label{sphericallike}$$

<a id='spherical'></a>

### Step 3.a.i: **`reference_metric::CoordSystem = "Spherical"`** \[Back to [top](#toc)\]
$$\label{spherical}$$

Standard spherical coordinates, with $(r,\theta,\phi)=(xx_0,xx_1,xx_2)$


```python
if CoordSystem == "Spherical":
    # Adding assumption real=True can help simplify expressions involving xx[0] & xx[1] below.
    xx[0] = sp.symbols("xx0", real=True)
    xx[1] = sp.symbols("xx1", real=True)

    RMAX = par.Cparameters("REAL", thismodule, ["RMAX"],10.0)
    xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
    xxmax = [         RMAX,          M_PI,  M_PI]

    r  = xx[0]
    th = xx[1]
    ph = xx[2]

    Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
    Cart_to_xx[1] = sp.acos(Cartz / Cart_to_xx[0])
    Cart_to_xx[2] = sp.atan2(Carty, Cartx)

    xxSph[0] = r
    xxSph[1] = th
    xxSph[2] = ph

    # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
    # Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
    xx_to_Cart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
    xx_to_Cart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
    xx_to_Cart[2] = xxSph[0]*sp.cos(xxSph[1])

    scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
    scalefactor_orthog[1] = xxSph[0]
    scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

    # Set the unit vectors
    UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                   [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                   [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]
```

Now let's analyze $r(xx_0)$ for **"Spherical"** coordinates.


```python
%matplotlib inline

CoordSystem = "Spherical"
par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric()

RMAX     = 10.0
r_of_xx0      = sp.sympify(str(rfm.xxSph[0]                   ).replace("RMAX",str(RMAX)))
rprime_of_xx0 = sp.sympify(str(sp.diff(rfm.xxSph[0],rfm.xx[0])).replace("RMAX",str(RMAX)))

create_r_of_xx0_plots(CoordSystem, r_of_xx0,rprime_of_xx0)
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_21_1.png)
    


<a id='sinhspherical'></a>

### Step 3.a.ii: **`reference_metric::CoordSystem = "SinhSpherical"`** \[Back to [top](#toc)\]
$$\label{sinhspherical}$$

Spherical coordinates, but with $$r(xx_0) = \text{AMPL} \frac{\sinh\left(\frac{xx_0}{\text{SINHW}}\right)}{\sinh\left(\frac{1}{\text{SINHW}}\right)}.$$

SinhSpherical uses two parameters: `AMPL` and `SINHW`. `AMPL` sets the outer boundary distance; and `SINHW` sets the focusing of the coordinate points near $r=0$, where a small `SINHW` ($\sim 0.125$) will greatly focus the points near $r=0$ and a large `SINHW` will look more like an ordinary spherical polar coordinate system.


```python
if CoordSystem == "SinhSpherical":
    xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
    xxmax = [sp.sympify(1),          M_PI,  M_PI]

    AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"],[10.0,0.2])
    # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
    r = AMPL * (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) / \
               (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
    th = xx[1]
    ph = xx[2]

    Cart_to_xx[0] = SINHW*sp.asinh(sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)*sp.sinh(1/SINHW)/AMPL)
    Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
    Cart_to_xx[2] = sp.atan2(Carty, Cartx)

    xxSph[0] = r
    xxSph[1] = th
    xxSph[2] = ph

    # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
    # Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
    xx_to_Cart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
    xx_to_Cart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
    xx_to_Cart[2] = xxSph[0]*sp.cos(xxSph[1])

    scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
    scalefactor_orthog[1] = xxSph[0]
    scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

    # Set the unit vectors
    UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                   [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                   [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]
```

Now we explore $r(xx_0)$ for `SinhSpherical` assuming `AMPL=10.0` and `SINHW=0.2`:


```python
%matplotlib inline

CoordSystem = "SinhSpherical"
par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric()

AMPL     = 10.0
SINHW    = 0.2
r_of_xx0      = sp.sympify(str(rfm.xxSph[0]                   ).replace("AMPL",str(AMPL)).replace("SINHW",str(SINHW)))
rprime_of_xx0 = sp.sympify(str(sp.diff(rfm.xxSph[0],rfm.xx[0])).replace("AMPL",str(AMPL)).replace("SINHW",str(SINHW)))

create_r_of_xx0_plots(CoordSystem, r_of_xx0,rprime_of_xx0)
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_25_1.png)
    


<a id='sinhsphericalv2'></a>

### Step 3.a.iii: **`reference_metric::CoordSystem = "SinhSphericalv2"`** \[Back to [top](#toc)\]
$$\label{sinhsphericalv2}$$

The same as SinhSpherical coordinates, but with an additional `AMPL*const_dr*xx_0` term:
$$r(xx_0) = \text{AMPL} \left[\text{const_dr}\ xx_0 + \frac{\sinh\left(\frac{xx_0}{\text{SINHW}}\right)}{\sinh\left(\frac{1}{\text{SINHW}}\right)}\right].$$


```python
if CoordSystem == "SinhSphericalv2":
    # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near xx[0]=0 to have
    # constant radial resolution of const_dr, provided the sinh() term does not dominate near xx[0]=0.
    xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
    xxmax = [sp.sympify(1),          M_PI,  M_PI]

    AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"],[10.0,0.2])
    const_dr = par.Cparameters("REAL",thismodule,["const_dr"],0.0625)

    r = AMPL*( const_dr*xx[0] + (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) /
               (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW)) )
    th = xx[1]
    ph = xx[2]

    # NO CLOSED-FORM EXPRESSION FOR RADIAL INVERSION.
    # Cart_to_xx[0] = "NewtonRaphson"
    # Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
    # Cart_to_xx[2] = sp.atan2(Carty, Cartx)

    xxSph[0] = r
    xxSph[1] = th
    xxSph[2] = ph

    # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
    # Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
    xx_to_Cart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
    xx_to_Cart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
    xx_to_Cart[2] = xxSph[0]*sp.cos(xxSph[1])

    scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
    scalefactor_orthog[1] = xxSph[0]
    scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

    # Set the unit vectors
    UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                   [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                   [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]
```

Now we explore $r(xx_0)$ for `SinhSphericalv2` assuming `AMPL=10.0`, `SINHW=0.2`, and `const_dr=0.05`. Notice that the `const_dr` term significantly increases the grid spacing near $xx_0=0$ relative to `SinhSpherical` coordinates.


```python
%matplotlib inline

CoordSystem = "SinhSphericalv2"
par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric()

AMPL     = 10.0
SINHW    = 0.2
const_dr = 0.05
r_of_xx0      = sp.sympify(str(rfm.xxSph[0]                   ).replace("AMPL",str(AMPL)).replace("SINHW",str(SINHW)).replace("const_dr",str(const_dr)))
rprime_of_xx0 = sp.sympify(str(sp.diff(rfm.xxSph[0],rfm.xx[0])).replace("AMPL",str(AMPL)).replace("SINHW",str(SINHW)).replace("const_dr",str(const_dr)))

create_r_of_xx0_plots(CoordSystem, r_of_xx0,rprime_of_xx0)
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_29_1.png)
    


<a id='cylindricallike'></a>

## Step 3.b: Cylindrical-like coordinate systems \[Back to [top](#toc)\]
$$\label{cylindricallike}$$

<a id='cylindrical'></a>

### Step 3.b.i: **`reference_metric::CoordSystem = "Cylindrical"`** \[Back to [top](#toc)\]
$$\label{cylindrical}$$

Standard cylindrical coordinates, with $(\rho,\phi,z)=(xx_0,xx_1,xx_2)$


```python
if CoordSystem == "Cylindrical":
    # Assuming the cylindrical radial coordinate
    #   is positive makes nice simplifications of
    #   unit vectors possible.
    xx[0] = sp.symbols("xx0", real=True)

    RHOMAX,ZMIN,ZMAX = par.Cparameters("REAL",thismodule,["RHOMAX","ZMIN","ZMAX"],[10.0,-10.0,10.0])
    xxmin = [sp.sympify(0), -M_PI, ZMIN]
    xxmax = [       RHOMAX,  M_PI, ZMAX]

    RHOCYL = xx[0]
    PHICYL = xx[1]
    ZCYL   = xx[2]

    Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2)
    Cart_to_xx[1] = sp.atan2(Carty, Cartx)
    Cart_to_xx[2] = Cartz

    xx_to_Cart[0] = RHOCYL*sp.cos(PHICYL)
    xx_to_Cart[1] = RHOCYL*sp.sin(PHICYL)
    xx_to_Cart[2] = ZCYL

    xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
    xxSph[1] = sp.acos(ZCYL / xxSph[0])
    xxSph[2] = PHICYL

    scalefactor_orthog[0] = sp.diff(RHOCYL,xx[0])
    scalefactor_orthog[1] = RHOCYL
    scalefactor_orthog[2] = sp.diff(ZCYL,xx[2])

    # Set the unit vectors
    UnitVectors = [[ sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
                   [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
                   [ sp.sympify(0),  sp.sympify(0),  sp.sympify(1)]]
```

Next let's plot **"Cylindrical"** coordinates.


```python
%matplotlib inline

import numpy as np                  # NumPy: A numerical methods module for Python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for projection=3d below.

R = np.linspace(0, 2, 24)
h = 2
u = np.linspace(0,  2*np.pi, 24)

x = np.outer(R, np.cos(u))
y = np.outer(R, np.sin(u))
z = h * np.outer(np.ones(np.size(u)), np.ones(np.size(u)))

r = np.arange(0,2,0.25)
theta = 2*np.pi*r*0

fig = plt.figure(figsize=(12,12)) # 8 in x 8 in
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1 = plt.axes(projection='polar')

ax1.set_rmax(2)

ax1.set_rgrids(r,labels=[])

thetas = np.linspace(0,360,24, endpoint=True)
ax1.set_thetagrids(thetas,labels=[])

# ax.grid(True)
ax1.grid(True,linewidth='1.0')
ax1.set_title("Top Down View")
plt.show()

ax2 = plt.axes(projection='3d', xticklabels=[], yticklabels=[], zticklabels=[])
#ax2.plot_surface(x,y,z, alpha=.75, cmap = 'viridis') # z in case of disk which is parallel to XY plane is constant and you can directly use h

x=np.linspace(-2, 2, 100)
z=np.linspace(-2, 2, 100)
Xc, Zc=np.meshgrid(x, z)
Yc = np.sqrt(4-Xc**2)

rstride = 10
cstride = 10
ax2.plot_surface(Xc, Yc, Zc, alpha=1.0, rstride=rstride, cstride=cstride, cmap = 'viridis')
ax2.plot_surface(Xc, -Yc, Zc, alpha=1.0, rstride=rstride, cstride=cstride, cmap = 'viridis')
ax2.set_title("Standard Cylindrical Grid in 3D")
ax2.grid(False)
plt.axis('off')

plt.show()
```


    <Figure size 864x864 with 0 Axes>



    
![png](output_34_1.png)
    



    
![png](output_34_2.png)
    


<a id='sinhcylindrical'></a>

### Step 3.b.ii" **`reference_metric::CoordSystem = "SinhCylindrical"`** \[Back to [top](#toc)\]
$$\label{sinhcylindrical}$$

Cylindrical coordinates, but with
$$\rho(xx_0) = \text{AMPLRHO} \frac{\sinh\left(\frac{xx_0}{\text{SINHWRHO}}\right)}{\sinh\left(\frac{1}{\text{SINHWRHO}}\right)}$$
and 
$$z(xx_2) = \text{AMPLZ} \frac{\sinh\left(\frac{xx_2}{\text{SINHWZ}}\right)}{\sinh\left(\frac{1}{\text{SINHWZ}}\right)}.$$


```python
if CoordSystem == "SinhCylindrical":
    # Assuming the cylindrical radial coordinate
    #   is positive makes nice simplifications of
    #   unit vectors possible.
    xx[0] = sp.symbols("xx0", real=True)

    xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
    xxmax = [sp.sympify(1),  M_PI, sp.sympify(+1)]

    AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = par.Cparameters("REAL",thismodule,
                                                       ["AMPLRHO","SINHWRHO","AMPLZ","SINHWZ"],
                                                       [     10.0,       0.2,   10.0,    0.2])

    # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
    RHOCYL = AMPLRHO * (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO))
    # phi coordinate remains unchanged.
    PHICYL = xx[1]
    ZCYL   = AMPLZ   * (sp.exp(xx[2] / SINHWZ)   - sp.exp(-xx[2] / SINHWZ))   / (sp.exp(1 / SINHWZ)   - sp.exp(-1 / SINHWZ))
    Cart_to_xx[0] = SINHWRHO*sp.asinh(sp.sqrt(Cartx ** 2 + Carty ** 2)*sp.sinh(1/SINHWRHO)/AMPLRHO)
    Cart_to_xx[1] = sp.atan2(Carty, Cartx)
    Cart_to_xx[2] = SINHWZ*sp.asinh(Cartz*sp.sinh(1/SINHWZ)/AMPLZ)

    xx_to_Cart[0] = RHOCYL*sp.cos(PHICYL)
    xx_to_Cart[1] = RHOCYL*sp.sin(PHICYL)
    xx_to_Cart[2] = ZCYL

    xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
    xxSph[1] = sp.acos(ZCYL / xxSph[0])
    xxSph[2] = PHICYL

    scalefactor_orthog[0] = sp.diff(RHOCYL,xx[0])
    scalefactor_orthog[1] = RHOCYL
    scalefactor_orthog[2] = sp.diff(ZCYL,xx[2])

    # Set the unit vectors
    UnitVectors = [[ sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
                   [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
                   [ sp.sympify(0),  sp.sympify(0),  sp.sympify(1)]]
```

Next let's plot **"SinhCylindrical"** coordinates.


```python
fig=plt.figure()


plt.clf()

fig = plt.figure()

ax = plt.subplot(1,1,1, projection='polar')

ax.set_rmax(2)

Nr = 20
xx0s = np.linspace(0,2,Nr, endpoint=True) + 1.0/(2.0*Nr)

rs = []
AMPLRHO = 1.0
SINHW = 0.4
for i in range(Nr):
    rs.append(AMPLRHO * (np.exp(xx0s[i] / SINHW) - np.exp(-xx0s[i] / SINHW)) / \
                        (np.exp(1.0 / SINHW) - np.exp(-1.0 / SINHW)))

ax.set_rgrids(rs,labels=[])

thetas = np.linspace(0,360,25, endpoint=True)
ax.set_thetagrids(thetas,labels=[])

# ax.grid(True)
ax.grid(True,linewidth='1.0')

plt.show()


```


    <Figure size 432x288 with 0 Axes>



    
![png](output_38_1.png)
    


<a id='sinhcylindricalv2'></a>

### Step 3.b.iii: **`reference_metric::CoordSystem = "SinhCylindricalv2"`** \[Back to [top](#toc)\]
$$\label{sinhcylindricalv2}$$

Cylindrical coordinates, but with
$$\rho(xx_0) = \text{AMPLRHO} \left[\text{const_drho}\ xx_0 + \frac{\sinh\left(\frac{xx_0}{\text{SINHWRHO}}\right)}{\sinh\left(\frac{1}{\text{SINHWRHO}}\right)}\right]$$
and 
$$z(xx_2) = \text{AMPLZ} \left[\text{const_dz}\ xx_2 + \frac{\sinh\left(\frac{xx_2}{\text{SINHWZ}}\right)}{\sinh\left(\frac{1}{\text{SINHWZ}}\right)}\right].$$


```python
if CoordSystem == "SinhCylindricalv2":
    # Assuming the cylindrical radial coordinate
    #   is positive makes nice simplifications of
    #   unit vectors possible.
    xx[0] = sp.symbols("xx0", real=True)

    # SinhCylindricalv2 adds the parameters "const_drho", "const_dz", which allows for regions near xx[0]=0
    # and xx[2]=0 to have constant rho and z resolution of const_drho and const_dz, provided the sinh() terms
    # do not dominate near xx[0]=0 and xx[2]=0.
    xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
    xxmax = [sp.sympify(1),  M_PI, sp.sympify(+1)]
    AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = par.Cparameters("REAL",thismodule,
                                                       ["AMPLRHO","SINHWRHO","AMPLZ","SINHWZ"],
                                                       [     10.0,       0.2,   10.0,    0.2])
    const_drho, const_dz = par.Cparameters("REAL",thismodule,["const_drho","const_dz"],[0.0625,0.0625])

    RHOCYL = AMPLRHO * ( const_drho*xx[0] + (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO)) )
    PHICYL = xx[1]
    ZCYL   = AMPLZ   * ( const_dz  *xx[2] + (sp.exp(xx[2] / SINHWZ  ) - sp.exp(-xx[2] / SINHWZ  )) / (sp.exp(1 / SINHWZ  ) - sp.exp(-1 / SINHWZ  )) )

    # NO CLOSED-FORM EXPRESSION FOR RADIAL OR Z INVERSION.
    # Cart_to_xx[0] = "NewtonRaphson"
    # Cart_to_xx[1] = sp.atan2(Carty, Cartx)
    # Cart_to_xx[2] = "NewtonRaphson"

    xx_to_Cart[0] = RHOCYL*sp.cos(PHICYL)
    xx_to_Cart[1] = RHOCYL*sp.sin(PHICYL)
    xx_to_Cart[2] = ZCYL

    xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
    xxSph[1] = sp.acos(ZCYL / xxSph[0])
    xxSph[2] = PHICYL

    scalefactor_orthog[0] = sp.diff(RHOCYL,xx[0])
    scalefactor_orthog[1] = RHOCYL
    scalefactor_orthog[2] = sp.diff(ZCYL,xx[2])

    # Set the unit vectors
    UnitVectors = [[ sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
                   [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
                   [ sp.sympify(0),  sp.sympify(0),  sp.sympify(1)]]
```

For example, let's set up **`SinhCylindricalv2`** coordinates and output the Christoffel symbol $\hat{\Gamma}^{xx_2}_{xx_2 xx_2}$, or more simply $\hat{\Gamma}^2_{22}$.


```python
par.set_parval_from_str("reference_metric::CoordSystem","SinhCylindricalv2")

rfm.reference_metric()

sp.pretty_print(sp.simplify(rfm.GammahatUDD[2][2][2]))
```

                             ⎛ 2⋅xx₂     ⎞    1                             
                             ⎜ ──────    ⎟  ──────                          
                             ⎜ SINHWZ    ⎟  SINHWZ                          
                            -⎝ℯ       - 1⎠⋅ℯ                                
    ────────────────────────────────────────────────────────────────────────
           ⎛                  ⎛   2       ⎞   xx₂     ⎛ 2⋅xx₂     ⎞    1   ⎞
           ⎜                  ⎜ ──────    ⎟  ──────   ⎜ ──────    ⎟  ──────⎟
           ⎜                  ⎜ SINHWZ    ⎟  SINHWZ   ⎜ SINHWZ    ⎟  SINHWZ⎟
    SINHWZ⋅⎝- SINHWZ⋅const_dz⋅⎝ℯ       - 1⎠⋅ℯ       - ⎝ℯ       + 1⎠⋅ℯ      ⎠


As we will soon see, defining these "hatted" quantities will be quite useful when expressing hyperbolic ([wave-equation](https://en.wikipedia.org/wiki/Wave_equation)-like) PDEs in non-Cartesian coordinate systems.

<a id='cartesianlike'></a>

## Step 3.c: Cartesian-like coordinate systems \[Back to [top](#toc)\]
$$\label{cartesianlike}$$

<a id='cartesian'></a>

### Step 3.c.i: **`reference_metric::CoordSystem = "Cartesian"`** \[Back to [top](#toc)\]
$$\label{cartesian}$$

Standard Cartesian coordinates, with $(x,y,z)=$ `(xx0,xx1,xx2)`.


```python
if CoordSystem == "Cartesian":
    xmin, xmax, ymin, ymax, zmin, zmax = par.Cparameters("REAL",thismodule,
                                                         ["xmin","xmax","ymin","ymax","zmin","zmax"],
                                                         [ -10.0,  10.0, -10.0,  10.0, -10.0,  10.0])
    xxmin = ["xmin", "ymin", "zmin"]
    xxmax = ["xmax", "ymax", "zmax"]

    xx_to_Cart[0] = xx[0]
    xx_to_Cart[1] = xx[1]
    xx_to_Cart[2] = xx[2]

    xxSph[0] = sp.sqrt(xx[0] ** 2 + xx[1] ** 2 + xx[2] ** 2)
    xxSph[1] = sp.acos(xx[2] / xxSph[0])
    xxSph[2] = sp.atan2(xx[1], xx[0])

    Cart_to_xx[0] = Cartx
    Cart_to_xx[1] = Carty
    Cart_to_xx[2] = Cartz

    scalefactor_orthog[0] = sp.sympify(1)
    scalefactor_orthog[1] = sp.sympify(1)
    scalefactor_orthog[2] = sp.sympify(1)

    # Set the transpose of the matrix of unit vectors
    UnitVectors = [[sp.sympify(1), sp.sympify(0), sp.sympify(0)],
                   [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                   [sp.sympify(0), sp.sympify(0), sp.sympify(1)]]
```


```python
%matplotlib inline

import numpy as np               # NumPy: A numerical methods module for Python
import matplotlib.pyplot as plt  # matplotlib: Python module specializing in plotting capabilities

plt.clf()

fig = plt.figure()
ax = fig.gca()
Nx = 16
ax.set_xticks(np.arange(0, 1., 1./Nx))
ax.set_yticks(np.arange(0, 1., 1./Nx))

for tick in ax.get_xticklabels():
    tick.set_rotation(60)

# plt.scatter(x, y)
ax.set_aspect('equal')

plt.grid()

# plt.savefig("Cartgrid.png",dpi=300)
plt.show()
# plt.close(fig)
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_47_1.png)
    


<a id='sinhcartesian'></a>

### Step 3.c.ii: **`reference_metric::CoordSystem = "SinhCartesian"`** \[Back to [top](#toc)\]
$$\label{sinhcartesian}$$

In this coordinate system, all three coordinates behave like the $z$-coordinate in SinhCylindrical coordinates, i.e.

$$
\begin{align}
x(xx_0) &= \text{AMPLXYZ} \left[\frac{\sinh\left(\frac{xx_0}{\text{SINHWXYZ}}\right)}{\sinh\left(\frac{1}{\text{SINHWXYZ}}\right)}\right]\ ,\\
y(xx_1) &= \text{AMPLXYZ} \left[\frac{\sinh\left(\frac{xx_1}{\text{SINHWXYZ}}\right)}{\sinh\left(\frac{1}{\text{SINHWXYZ}}\right)}\right]\ ,\\
z(xx_2) &= \text{AMPLXYZ} \left[\frac{\sinh\left(\frac{xx_2}{\text{SINHWXYZ}}\right)}{\sinh\left(\frac{1}{\text{SINHWXYZ}}\right)}\right].
\end{align}
$$


```python
if CoordSystem == "SinhCartesian":
    # SinhCartesian coordinates allows us to push the outer boundary of the
    # computational domain a lot further away, while keeping reasonably high
    # resolution towards the center of the computational grid.

    # Set default values for min and max (x,y,z)
    xxmin = [sp.sympify(-1), sp.sympify(-1), sp.sympify(-1)]
    xxmax = [sp.sympify(+1), sp.sympify(+1), sp.sympify(+1)]

    # Declare basic parameters of the coordinate system and their default values
    AMPLXYZ, SINHWXYZ = par.Cparameters("REAL", thismodule,
                                        ["AMPLXYZ", "SINHWXYZ"],
                                        [     10.0,        0.2])

    # Compute (xx_to_Cart0,xx_to_Cart1,xx_to_Cart2) from (xx0,xx1,xx2)
    for ii in [0, 1, 2]:
        xx_to_Cart[ii] = AMPLXYZ*(sp.exp(xx[ii]/SINHWXYZ) - sp.exp(-xx[ii]/SINHWXYZ))/(sp.exp(1/SINHWXYZ) - sp.exp(-1/SINHWXYZ))

    # Compute (r,th,ph) from (xx_to_Cart2,xx_to_Cart1,xx_to_Cart2)
    xxSph[0] = sp.sqrt(xx_to_Cart[0] ** 2 + xx_to_Cart[1] ** 2 + xx_to_Cart[2] ** 2)
    xxSph[1] = sp.acos(xx_to_Cart[2] / xxSph[0])
    xxSph[2] = sp.atan2(xx_to_Cart[1], xx_to_Cart[0])

    # Compute (xx0,xx1,xx2) from (Cartx,Carty,Cartz)
    Cart_to_xx[0] = SINHWXYZ*sp.asinh(Cartx*sp.sinh(1/SINHWXYZ)/AMPLXYZ)
    Cart_to_xx[1] = SINHWXYZ*sp.asinh(Carty*sp.sinh(1/SINHWXYZ)/AMPLXYZ)
    Cart_to_xx[2] = SINHWXYZ*sp.asinh(Cartz*sp.sinh(1/SINHWXYZ)/AMPLXYZ)

    # Compute scale factors
    scalefactor_orthog[0] = sp.diff(xx_to_Cart[0], xx[0])
    scalefactor_orthog[1] = sp.diff(xx_to_Cart[1], xx[1])
    scalefactor_orthog[2] = sp.diff(xx_to_Cart[2], xx[2])

    # Set the transpose of the matrix of unit vectors
    UnitVectors = [[sp.sympify(1), sp.sympify(0), sp.sympify(0)],
                   [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                   [sp.sympify(0), sp.sympify(0), sp.sympify(1)]]

```


```python
%matplotlib inline

import numpy as np               # NumPy: A numerical methods module for Python
import matplotlib.pyplot as plt  # matplotlib: Python module specializing in plotting capabilities

plt.clf()

fig = plt.figure(dpi=160)
ax = fig.gca()

# Set plot title
ax.set_title(r"$z=0$ slice of the 3D grid")

# Set SINH parameters. Here we assume:
#
# AMPLX  = AMPLY  = SINHA
# SINHWX = SINHWY = SINHW
SINHA = 10.0
SINHW = 0.45

# Set number of points. We assume the same point
# distribution along the (x,y)-directions
Nxxs = 24
xxis = np.linspace(-1,1,Nxxs, endpoint=True)

# Compute axis ticks by evaluating x and y using SinhCartesian coordinates
axis_ticks = []
for i in range(Nxxs):
    axis_ticks.append(SINHA * (np.exp(xxis[i] / SINHW) - np.exp(-xxis[i] / SINHW)) / \
                        (np.exp(1.0 / SINHW) - np.exp(-1.0 / SINHW)))

# Set the axis ticks
ax.set_xticks(axis_ticks)
ax.set_yticks(axis_ticks)

# Set x and y labels. Initialize array with empty strings
labelsx = ["" for i in range(Nxxs)]
labelsy = ["" for i in range(Nxxs)]

# Set x_min and x_max tick label
labelsx[0] = r"-AMPLX"
labelsx[-1] = r"AMPLX"

# Set y_min and y_max tick label
labelsy[0] = r"-AMPLY"
labelsy[-1] = r"AMPLY"

# Set tick labels
ax.set_xticklabels(labelsx)
ax.set_yticklabels(labelsy)

# Rotate x labels by 60 degrees
for tick in ax.get_xticklabels():
    tick.set_rotation(60)

# Draw the x=0 and y=0 ticklabel
ax.text(0,-11,"0",ha="center",va="center")
ax.text(-11,0,"0",ha="center",va="center")

# plt.scatter(x, y)
ax.set_aspect('equal')

plt.grid(color='black',linewidth=0.3)

plt.show()
# plt.savefig("Cartgrid.png",dpi=400)
# plt.close(fig)
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_50_1.png)
    


<a id='prolatespheroidal'></a>

## Step 3.d: [Prolate spheroidal](https://en.wikipedia.org/wiki/Prolate_spheroidal_coordinates)-like coordinate systems \[Back to [top](#toc)\]
$$\label{prolatespheroidal}$$

<a id='symtp'></a>

### Step 3.d.i: **`reference_metric::CoordSystem = "SymTP"`** \[Back to [top](#toc)\]
$$\label{symtp}$$

The Symmetric TwoPuncture (SymTP) coordinate system is obtained by slightly modifying [prolate spheroidal coordinates](https://en.wikipedia.org/wiki/Prolate_spheroidal_coordinates) (PSC). Standard PSC are related to Cartesian coordinates $(x,y,z)$ via

$$
\begin{aligned}
x &= a\sinh\mu\sin\nu\cos\varphi,\\
y &= a\sinh\mu\sin\nu\sin\varphi,\\
z &= a\cosh\mu\cos\nu = \left(a^{2}\sinh^{2}\mu + a^{2}\right)^{1/2}\cos\nu,
\end{aligned}
$$

where $\mu\in[0,\infty)$, $\nu\in[0,\pi]$, and $\varphi\in[0,2\pi]$, and we have used the [identity](https://en.wikipedia.org/wiki/Hyperbolic_functions)

$$
\cosh^{2}\mu - \sinh^{2}\mu = 1.
$$

PSC have two foci located at $z=\pm a$, where the grid lines get more dense and the grid resolution increases. However, note that the parameter $a$ controls both the foci position *and* the grid scaling, which is not a particularly interesting property.

In order to remedy this, we introduce new coordinates $(xx_{0},xx_{1},xx_{2})$, such that

$$
\begin{aligned}
xx_{0} &= \frac{1}{a}\sinh\mu,\\
xx_{1} &= \nu,\\
xx_{2} &= \varphi,
\end{aligned}
$$

and change $a\to \text{bScale}$, so that we obtain

$$
\begin{aligned}
x &= xx_{0}\sin(xx_{1})\cos(xx_{2}),\\
y &= xx_{0}\sin(xx_{1})\sin(xx_{2}),\\
z &= \left(xx_{0}^{2}+\text{bScale}^{2}\right)\cos(xx_{1}),
\end{aligned}
$$

Note that, for numerical convenience, we change the range of $xx_{2}$ to $[-\pi,\pi]$. Comparing SymTP coordinates with Cylindrical coordinates, we find that $(\rho,\phi,z)=(xx_{0}\sin(xx_{1}), xx_{2}, \sqrt{xx_{0}^2 + \text{bScale}^2}\cos(xx_{1}))$.


```python
if CoordSystem == "SymTP":

    var1, var2= sp.symbols('var1 var2',real=True)
    bScale, AW, AMAX, RHOMAX, ZMIN, ZMAX = par.Cparameters("REAL",thismodule,
                                                           ["bScale","AW","AMAX","RHOMAX","ZMIN","ZMAX"],
                                                           [0.5,     0.2,   10.0,    10.0, -10.0,  10.0])

    # Assuming xx0, xx1, and bScale
    #   are positive makes nice simplifications of
    #   unit vectors possible.
    xx[0],xx[1] = sp.symbols("xx0 xx1", real=True)

    xxmin = [sp.sympify(0), sp.sympify(0),-M_PI]
    xxmax = [         AMAX,          M_PI, M_PI]

    AA = xx[0]

    var1 = sp.sqrt(AA**2 + (bScale * sp.sin(xx[1]))**2)
    var2 = sp.sqrt(AA**2 + bScale**2)

    RHOSYMTP = AA*sp.sin(xx[1])
    PHSYMTP = xx[2]
    ZSYMTP = var2*sp.cos(xx[1])

    xx_to_Cart[0] = AA  *sp.sin(xx[1])*sp.cos(xx[2])
    xx_to_Cart[1] = AA  *sp.sin(xx[1])*sp.sin(xx[2])
    xx_to_Cart[2] = ZSYMTP

    xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
    xxSph[1] = sp.acos(ZSYMTP / xxSph[0])
    xxSph[2] = PHSYMTP

    rSph  = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
    thSph = sp.acos(Cartz / rSph)
    phSph = sp.atan2(Carty, Cartx)

    # Mathematica script to compute Cart_to_xx[]
    #             AA = x1;
    #             var2 = Sqrt[AA^2 + bScale^2];
    #             RHOSYMTP = AA*Sin[x2];
    #             ZSYMTP = var2*Cos[x2];
    #             Solve[{rSph == Sqrt[RHOSYMTP^2 + ZSYMTP^2],
    #                    thSph == ArcCos[ZSYMTP/Sqrt[RHOSYMTP^2 + ZSYMTP^2]],
    #                    phSph == x3},
    #                   {x1, x2, x3}]
    Cart_to_xx[0] = sp.sqrt(-bScale**2 + rSph**2 +
                            sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 -
                                    4*bScale**2*rSph**2*sp.cos(thSph)**2))*M_SQRT1_2 # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

    # The sign() function in the following expression ensures the correct root is taken.
    Cart_to_xx[1] = sp.acos(sp.sign(Cartz)*(
                              sp.sqrt(1 + rSph**2/bScale**2 -
                                      sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 -
                                              4*bScale**2*rSph**2*sp.cos(thSph)**2)/bScale**2)*M_SQRT1_2)) # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

    Cart_to_xx[2] = phSph
```


```python
%matplotlib inline

Nxxs   = 24
xx0    = np.linspace(-2,2,Nxxs,endpoint=True)
xx1    = np.linspace(0,np.pi,Nxxs,endpoint=True)
xx2    = 0
bScale = 1

def x(xx0,xx1):
    return xx0 * np.sin(xx1)
def z(xx0,xx1,bScale):
    return np.sqrt(xx0**2 + bScale**2) * np.cos(xx1)

import numpy as np               # NumPy: A numerical methods module for Python
import matplotlib.pyplot as plt  # matplotlib: Python module specializing in plotting capabilities

plt.clf()

fig = plt.figure(dpi=160)
ax = fig.gca()

# Set plot title
ax.set_title(r"""SymTP Coordinates: zx-plane ($xx_{2}$=0 and $xx_{2}=\pi$)
Blue (red) lines have constant $xx_{0}$ ($xx_{1}$)""")

ax.set_xlim(-2.5,2.5)
ax.set_ylim(-2.5,2.5)
ax.set_aspect('equal')
ax.axis('off')

for i0 in range(Nxxs):
    plt.plot(x(xx0[i0],xx1),z(xx0[i0],xx1,bScale),'b',lw=0.75)

for i1 in range(Nxxs):
    plt.plot(x(xx0,xx1[i1]),z(xx0,xx1[i1],bScale),'r',lw=0.75)

plt.show()
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_54_1.png)
    


<a id='sinhsymtp'></a>

### Step 3.d.ii: **`reference_metric::CoordSystem = "SinhSymTP"`** \[Back to [top](#toc)\]
$$\label{sinhsymtp}$$

SinhSymTP coordinates are obtained from SymTP coordinates by making the substitution

$$
xx0 \to \mathcal{A}\frac{\sinh(xx_{0}/w)}{\sinh(1/w)},
$$

with the further modification that $xx_{0}\in[0,1]$. The parameter $\mathcal{A}$ controls the scale of the grid, while the parameter $w$ controls how densily sampled the region around the foci are.


```python
if CoordSystem == "SinhSymTP":

    var1, var2= sp.symbols('var1 var2',real=True)
    bScale, AW, AMAX, RHOMAX, ZMIN, ZMAX = par.Cparameters("REAL",thismodule,
                                                           ["bScale","AW","AMAX","RHOMAX","ZMIN","ZMAX"],
                                                           [0.5,     0.2,   10.0,    10.0, -10.0,  10.0])

    # Assuming xx0, xx1, and bScale
    #   are positive makes nice simplifications of
    #   unit vectors possible.
    xx[0],xx[1] = sp.symbols("xx0 xx1", real=True)

    xxmin = [sp.sympify(0), sp.sympify(0),-M_PI]
    xxmax = [sp.sympify(1),          M_PI, M_PI]

    AA = AMAX * (sp.exp(xx[0]/SINHWAA) - sp.exp(-xx[0]/SINHWAA)) / (sp.exp(1/SINHWAA) - sp.exp(-1/SINHWAA))

    var1 = sp.sqrt(AA**2 + (bScale * sp.sin(xx[1]))**2)
    var2 = sp.sqrt(AA**2 + bScale**2)

    RHOSYMTP = AA*sp.sin(xx[1])
    PHSYMTP = xx[2]
    ZSYMTP = var2*sp.cos(xx[1])

    xx_to_Cart[0] = AA*sp.sin(xx[1])*sp.cos(xx[2])
    xx_to_Cart[1] = AA*sp.sin(xx[1])*sp.sin(xx[2])
    xx_to_Cart[2] = ZSYMTP

    xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
    xxSph[1] = sp.acos(ZSYMTP / xxSph[0])
    xxSph[2] = PHSYMTP

    rSph  = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
    thSph = sp.acos(Cartz / rSph)
    phSph = sp.atan2(Carty, Cartx)

    # Mathematica script to compute Cart_to_xx[]
    #             AA = x1;
    #             var2 = Sqrt[AA^2 + bScale^2];
    #             RHOSYMTP = AA*Sin[x2];
    #             ZSYMTP = var2*Cos[x2];
    #             Solve[{rSph == Sqrt[RHOSYMTP^2 + ZSYMTP^2],
    #                    thSph == ArcCos[ZSYMTP/Sqrt[RHOSYMTP^2 + ZSYMTP^2]],
    #                    phSph == x3},
    #                   {x1, x2, x3}]
    Cart_to_xx[0] = sp.sqrt(-bScale**2 + rSph**2 +
                            sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 -
                                    4*bScale**2*rSph**2*sp.cos(thSph)**2))*M_SQRT1_2 # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

    # The sign() function in the following expression ensures the correct root is taken.
    Cart_to_xx[1] = sp.acos(sp.sign(Cartz)*(
                              sp.sqrt(1 + rSph**2/bScale**2 -
                                      sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 -
                                              4*bScale**2*rSph**2*sp.cos(thSph)**2)/bScale**2)*M_SQRT1_2)) # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

    Cart_to_xx[2] = phSph

    scalefactor_orthog[0] = sp.diff(AA,xx[0]) * var1 / var2
    scalefactor_orthog[1] = var1
    scalefactor_orthog[2] = AA * sp.sin(xx[1])

    # Set the transpose of the matrix of unit vectors
    UnitVectors = [[sp.sin(xx[1]) * sp.cos(xx[2]) * var2 / var1,
                    sp.sin(xx[1]) * sp.sin(xx[2]) * var2 / var1,
                    AA * sp.cos(xx[1]) / var1],
                   [AA * sp.cos(xx[1]) * sp.cos(xx[2]) / var1,
                    AA * sp.cos(xx[1]) * sp.sin(xx[2]) / var1,
                        -sp.sin(xx[1]) * var2 / var1],
                   [-sp.sin(xx[2]), sp.cos(xx[2]), sp.sympify(0)]]
```


```python
%matplotlib inline

import numpy as np               # NumPy: A numerical methods module for Python
import matplotlib.pyplot as plt  # matplotlib: Python module specializing in plotting capabilities

Nxx0   = 12
Nxx1   = 24
xx0    = np.linspace(0,1,Nxx0,endpoint=True)
xx1    = np.linspace(0,np.pi,Nxx1,endpoint=True)
xx2    = 0
bScale = 1
sinhA  = 2.2
sinhW  = 0.3

def rtilde(xx0,sinhA,sinhW):
    return sinhA * np.sinh(xx0/sinhW) / np.sinh(1.0/sinhW)
def x(xx0,xx1,sinhA,sinhW):
    return rtilde(xx0,sinhA,sinhW) * np.sin(xx1)
def z(xx0,xx1,bScale,sinhA,sinhW):
    return np.sqrt(rtilde(xx0,sinhA,sinhW)**2 + bScale**2) * np.cos(xx1)

plt.clf()

fig = plt.figure(dpi=160)
ax = fig.gca()

# Set plot title
ax.set_title(r"""SinhSymTP Coordinates: zx-plane ($xx_{2}$=0 and $xx_{2}=\pi$)
Blue (red) lines have constant $xx_{0}$ ($xx_{1}$)""")

ax.set_xlim(-2.5,2.5)
ax.set_ylim(-2.5,2.5)
ax.set_aspect('equal')
ax.axis('off')

for i0 in range(Nxx0):
    plt.plot(x(xx0[i0],xx1,sinhA,sinhW),z(xx0[i0],xx1,bScale,sinhA,sinhW),'b',lw=0.75)
    plt.plot(-x(xx0[i0],xx1,sinhA,sinhW),z(xx0[i0],xx1,bScale,sinhA,sinhW),'b',lw=0.75)

for i1 in range(Nxx1):
    plt.plot(x(xx0,xx1[i1],sinhA,sinhW),z(xx0,xx1[i1],bScale,sinhA,sinhW),'r',lw=0.75)
    plt.plot(-x(xx0,xx1[i1],sinhA,sinhW),z(xx0,xx1[i1],bScale,sinhA,sinhW),'r',lw=0.75)

plt.show()
```


    <Figure size 432x288 with 0 Axes>



    
![png](output_57_1.png)
    


<a id='latex_pdf_output'></a>

# Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Reference_Metric.pdf](Tutorial-Reference_Metric.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Reference_Metric")
```

    Created Tutorial-Reference_Metric.tex, and compiled LaTeX file to PDF file
        Tutorial-Reference_Metric.pdf

