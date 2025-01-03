<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# C code generation of GRHD Equations

## Authors: Phil Chang & Zach Etienne
### Formatting improvements courtesy Brandon Clark

## This notebook demonstrates the C code generation of General Relativistic HydroDynamics (GRHD) equations in conservative form using a specific state vector ${\boldsymbol{\mathcal{U}}}=(\rho_*,\tilde{S},\tilde{\tau})$. The process of transitioning between primitive and conservative variables, computing flux and source terms, as well as performing Lorentz boosts and other transformations is detailed.
[comment]: <> (omit white space somehow?)
$\newcommand{\be}{\begin{equation}}$
$\newcommand{\ee}{\end{equation}}$
$\newcommand{\grad}{{\boldsymbol{\nabla}}}$
$\newcommand{\vel}{{\boldsymbol{v}}}$
$\newcommand{\mom}{{\boldsymbol{p}}}$
$\newcommand{\ddt}[1]{{\frac{\partial #1}{\partial t}}}$
$\newcommand{\ddx}[1]{{\frac{\partial #1}{\partial x}}}$
$\newcommand{\state}{{\boldsymbol{\mathcal{U}}}}$
$\newcommand{\charge}{{\boldsymbol{U}}}$
$\newcommand{\psicharge}{{\boldsymbol{\psi}}}$
$\newcommand{\lapse}{\alpha}$
$\newcommand{\shift}{\boldsymbol{\beta}}$
$\newcommand{\rhostar}{{\rho_*}}$
$\newcommand{\tautilde}{{\tilde{\tau}}}$
$\newcommand{\Svectilde}{{\tilde{\boldsymbol{S}}}}$
$\newcommand{\rtgamma}{{\sqrt{\gamma}}}$
$\newcommand{\T}[2]{{T^{#1 #2}}}$
$\newcommand{\uvec}{{\boldsymbol{u}}}$
$\newcommand{\Vvec}{{\boldsymbol{\mathcal{V}}}}$
$\newcommand{\vfluid}{{\boldsymbol{v}_{\rm f}}}$
$\newcommand{\vVal}{{\tilde{\boldsymbol{v}}}}$

$\newcommand{\flux}{{\boldsymbol{\mathcal{F}}}}$
$\newcommand{\fluxV}{{\boldsymbol{F}}}$
$\newcommand{\source}{{\boldsymbol{\mathcal{S}}}}$
$\newcommand{\sourceV}{{\boldsymbol{S}}}$

$\newcommand{\area}{{\boldsymbol{A}}}$
$\newcommand{\normal}{{\hat{\boldsymbol{n}}}}$
$\newcommand{\pt}{{\boldsymbol{p}}}$
$\newcommand{\nb}{{\boldsymbol{n}}}$
$\newcommand{\meshv}{{\boldsymbol{w}}}$
$\newcommand{\facev}{{\boldsymbol{\tilde{w}}_{ij}}}$
$\newcommand{\facer}{{\boldsymbol{\tilde{r}}_{ij}}}$
$\newcommand{\meshr}{{\boldsymbol{r}}}$
$\newcommand{\cmr}{{\boldsymbol{c}}}$

## Introduction: 
We start out with the ** GRHD ** equations in conservative form with the state vector $\state=(\rhostar, \Svectilde, \tautilde)$:
\begin{equation}
\ddt{\state} + \grad\cdot\flux = \source,
\end{equation}
where $\rhostar = \lapse\rho\rtgamma u^0$, $\Svectilde = \rhostar h \uvec$, $\tautilde = \lapse^2\rtgamma \T00 - \rhostar$. The associated set of primitive variables are $(\rho, \vel, \epsilon)$, which are the rest mass density, fluid 3-velocity, and internal energy (measured in the rest frame).  

The flux, $\flux$ is given by
\begin{equation}
 \flux=(\rhostar \vel, \lapse\rtgamma\T{j}{\beta}g_{\beta i}, \lapse^2\rtgamma\T0j - \rhostar\vel
\end{equation}
where $\vel$ is the 3-velocity, and $\source = (0, \frac 1 2 \lapse\rtgamma \T{\lapse}{\beta}g_{\lapse\beta,i}, s)$ is the source function, and
\begin{equation}
s = \lapse\rtgamma\left[\left(\T00\beta^i\beta^j + 2\T0i\beta^j\right)K_{ij} - \left(\T00\beta^i + \T0i\right)\partial_i\lapse\right]
\end{equation}
The stress-energy tensor for a perfect fluid is written as 
\begin{equation}
\T{\mu}{\nu} = \rho h u^{\mu} u^{\nu} + P g^{\mu\nu},
\end{equation}
where $h = 1 + \epsilon + P/\rho$ is the specific enthalpy and $u^{\mu}$ are the respective components of the four velocity.  

Noting that the mass $\flux$ is defined in terms of $\rhostar$ and $\vel$, we need to first find a mapping between $\vel$ and $u$.  

### Alternative formulation

The Athena++ folks have an alternative formulation that might be superior.  
Begin with the continuity equation
\begin{equation}
\grad_{\mu}\rho u^{\mu} = 0,
\end{equation}
where $\grad$ is the covariant derivative.  This can be mapped directly to 
\begin{equation}
\partial_{0} \sqrt{-g}\rho u^0 + \partial_i\sqrt{-g} \rho u^0 v^i = 0 
\end{equation}
which we can identify with $\rhostar = \alpha\rtgamma \rho u^0$ because $\sqrt{-g} = \alpha\rtgamma$.

Now the second equation is the conservation of energy-momentum which we write as
\begin{equation}
\grad_{\nu}T^{\nu}_{\mu} = 0 
\end{equation}
writing this out we have 
\begin{equation}
\partial_0 g_{\mu\alpha}T^{\alpha 0} + \partial_i g_{\mu\alpha}T^{\alpha i} - \Gamma_{\mu\alpha}^{\gamma} g_{\gamma\beta}T^{\alpha\beta} = 0 
\end{equation}
Noting that
\begin{equation}
\Gamma^{\alpha}_{\beta\gamma} = \frac 1 2 g^{\alpha\delta}\left(\partial_{\gamma}g_{\beta\delta} + \partial_{\beta}g_{\gamma\delta} - \partial_{\delta}g_{\beta\gamma}\right)
\end{equation}
Writing this all out, we note the last term is
\begin{equation}
\Gamma_{\mu\alpha}^{\gamma} g_{\gamma\beta}T^{\alpha\beta} =
\frac 1 2 g^{\gamma\delta}\left(\partial_{\alpha}g_{\mu\delta} + \partial_{\mu}g_{\alpha \delta} - \partial_{\delta}g_{\mu \alpha}\right) T_{\gamma}^{\alpha} = 
\frac 1 2 \left(\partial_{\alpha}g_{\mu\delta} + \partial_{\mu}g_{\alpha \delta} - \partial_{\delta}g_{\mu \alpha}\right)
T^{\alpha\delta}
\end{equation}
We sum over $\alpha$ and $\delta$, but note that we are antisymmetric in first and last terms in $\alpha$ and $\delta$ in the () but symmetric in $T_{\alpha\delta}$ so we have
\begin{equation}
\Gamma_{\mu\alpha}^{\gamma} g_{\gamma\beta}T^{\alpha\beta} = \frac 1 2 \partial_{\mu}g_{\alpha \delta} T^{\alpha\delta}
\end{equation}

Thus we have 
\begin{equation}
\partial_0 T^{0}_{\mu} + \partial_i T^{i}_{\mu} = \frac 1 2 \partial_{\mu}g_{\alpha \delta} T^{\alpha\delta}
\end{equation}
The $\mu = (1,2,3)$, we almost get back the equations in the standard formulation
\begin{equation}
\partial_0 \rho h u^0 u_i + \partial_j T^j_i = \frac 1 2 \partial_{i}g_{\alpha \delta} T^{\alpha\delta},
\end{equation}
which modulo factors of $\lapse\rtgamma$ in front is the same as the "standard" equations.

The $T^0_0$ term is more interesting. Here we have
\begin{equation}
\partial_0 (\rho h u^0 u_0 + + \partial_j T^j_i = \frac 1 2 \partial_{0}g_{\alpha \delta} T^{\alpha\delta},
\end{equation}

However, the disadvantage is that we need the time derivative of the metric.

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This notebook is organized as follows

1. [Step 1](#mapping): Primitive to Conservative Mapping
1. [Step 2](#zach): Compute $u^0$ from the Valencia 3-velocity (Zach step)
1. [Step 3](#flux): Compute the flux
1. [Step 4](#source): Source Terms
1. [Step 5](#rotation): Rotation
1. [Step 6](#solver): Conservative to Primitive Solver
1. [Step 7](#lorentz): Lorentz Boosts
1. [Step 8](#TmunuSph): Compute $T^{\mu\nu}$ in Cartesian Coordinates
1. [Step 9](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='mapping'></a>

# Step 1: Primitive to Conservative Mapping
$$\label{mapping}$$

We want to make a mapping from the primitives to conserved variables following Zach notebook:
\begin{equation}
(\rho, \vel, \epsilon) \rightarrow (\rhostar = \lapse\rho\rtgamma u^0, \Svectilde = \rhostar h \uvec, \tautilde = \lapse^2\rtgamma \T00 - \rhostar).
\end{equation}



```python
import GRHD.equations as Ge  # NRPy: Implementation of GRHD equations in Cartesian coordinates
import indexedexp as ixp     # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sympy as sp           # SymPy: The Python computer algebra package upon which NRPy+ depends
from outputC import outputC  # NRPy+: Basic C code output functionality

# declare gammaDD
gammaDD = ixp.zerorank2()
components = ["xx", "xy", "xz", "yy", "yz", "zz"]
names = ""
for comp in components :
    names = names + "mi.gamDD{0} ".format(comp)

gxx, gxy, gxz, gyy, gyz, gzz = sp.symbols( names)

gammaDD[0][0] = gxx
gammaDD[0][1] = gxy
gammaDD[0][2] = gxz
gammaDD[1][0] = gxy
gammaDD[1][1] = gyy
gammaDD[1][2] = gyz
gammaDD[2][0] = gxz
gammaDD[2][1] = gyz
gammaDD[2][2] = gzz

#declare alpha
alpha = sp.symbols( "mi.alpha")

#declare beta
betaU = ixp.zerorank1()
for i, comp in enumerate(["X", "Y", "Z"]) :
    betaU[i] = sp.symbols( "mi.beta{0}".format(comp), real=True)

#now get the primitives
rho_b, epsilon, P = sp.symbols("rho ie p")

#get the 3-velocities
vU = ixp.zerorank1()
for i, comp in enumerate( ["vx", "vy", "vz"]) :
    vU[i] = sp.symbols("{0}".format(comp))

Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha,betaU,gammaDD, vU)

u4U = Ge.u4U_ito_vU
# Zach says: Probably want to adopt speed-limited vU[i], Ge.rescaledvU[i], here, a la:
# for i in range(3):
# ... vU[i] = Ge.rescaledvU[i]

# First compute stress-energy tensor T4UU and T4UD:
Ge.compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)
Ge.compute_T4UD(gammaDD,betaU,alpha, Ge.T4UU)

# Next sqrt(gamma)
Ge.compute_sqrtgammaDET(gammaDD)

# Compute conservative variables in terms of primitive variables
Ge.compute_rho_star( alpha, Ge.sqrtgammaDET, rho_b, u4U)
Ge.compute_tau_tilde(alpha, Ge.sqrtgammaDET, Ge.T4UU,Ge.rho_star)
Ge.compute_S_tildeD( alpha, Ge.sqrtgammaDET, Ge.T4UD)

# Zach says: Why only output u^x? Debugging reasons?
outputC([u4U[1],     Ge.rho_star, Ge.S_tildeD[0], Ge.S_tildeD[1], Ge.S_tildeD[2], Ge.tau_tilde],
        ["u4U1", "con[iRhoStar]",     "con[iSx]",     "con[iSy]",     "con[iSz]",  "con[iTau]"],
        filename="NRPY+prim2Con.h", params="outCverbose=False")
!cat NRPY+prim2Con.h

outputC([Ge.sqrtgammaDET*alpha], ["detg"], filename="NRPY+detg.h")
!cat NRPY+detg.h


gammaUU, gammabarDet = ixp.symm_matrix_inverter3x3(gammaDD)
outputC([gammaUU[0][0],gammaUU[0][1],gammaUU[0][2],gammaUU[1][1],gammaUU[1][2],gammaUU[2][2]],
        [    "gamUUxx",    "gamUUxy",    "gamUUxz",    "gamUUyy",    "gamUUyz",    "gamUUzz"],
        filename="NRPY+gamUU.h")
!cat NRPY+gamUU.h
```

    Wrote to file "NRPY+prim2Con.h"
    {
      const double tmp_0 = (1.0/((GAMMA_SPEED_LIMIT)*(GAMMA_SPEED_LIMIT)));
      const double tmp_3 = (1.0/((mi.alpha)*(mi.alpha)));
      const double tmp_4 = mi.betaX + vx;
      const double tmp_5 = mi.gamDDxx*tmp_3*((tmp_4)*(tmp_4));
      const double tmp_7 = mi.betaY + vy;
      const double tmp_8 = mi.gamDDyy*tmp_3*((tmp_7)*(tmp_7));
      const double tmp_10 = mi.betaZ + vz;
      const double tmp_11 = mi.gamDDzz*((tmp_10)*(tmp_10))*tmp_3;
      const double tmp_14 = mi.gamDDxy*tmp_3*tmp_4*tmp_7;
      const double tmp_15 = mi.gamDDxz*tmp_10*tmp_3*tmp_4;
      const double tmp_16 = mi.gamDDyz*tmp_10*tmp_3*tmp_7;
      const double tmp_20 = (1.0/2.0)*fabs(-tmp_0 - tmp_11 - 2*tmp_14 - 2*tmp_15 - 2*tmp_16 - tmp_5 - tmp_8 + 1);
      const double tmp_21 = (1.0/2.0)*tmp_0 - 1.0/2.0*tmp_11 - tmp_14 - tmp_15 - tmp_16 + tmp_20 - 1.0/2.0*tmp_5 - 1.0/2.0*tmp_8 + 1.0/2.0;
      const double tmp_22 = (1.0/sqrt(tmp_21));
      const double tmp_23 = sqrt((-1.0/2.0*tmp_0 + (1.0/2.0)*tmp_11 + tmp_14 + tmp_15 + tmp_16 - tmp_20 + (1.0/2.0)*tmp_5 + (1.0/2.0)*tmp_8 + 1.0/2.0)/(TINYDOUBLE + tmp_11 + 2*tmp_14 + 2*tmp_15 + 2*tmp_16 + tmp_5 + tmp_8));
      const double tmp_24 = -mi.betaX + tmp_23*tmp_4;
      const double tmp_25 = sqrt(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy);
      const double tmp_26 = rho*tmp_22*tmp_25;
      const double tmp_27 = p*tmp_3;
      const double tmp_28 = rho*tmp_3*(ie + p/rho + 1)/tmp_21;
      const double tmp_29 = -tmp_27 + tmp_28;
      const double tmp_30 = mi.betaX*tmp_27 + tmp_24*tmp_28;
      const double tmp_31 = mi.betaY*tmp_27 + tmp_28*(-mi.betaY + tmp_23*tmp_7);
      const double tmp_32 = mi.betaZ*tmp_27 + tmp_28*(-mi.betaZ + tmp_10*tmp_23);
      const double tmp_33 = mi.alpha*tmp_25;
      u4U1 = tmp_22*tmp_24/mi.alpha;
      con[iRhoStar] = tmp_26;
      con[iSx] = tmp_33*(mi.gamDDxx*tmp_30 + mi.gamDDxy*tmp_31 + mi.gamDDxz*tmp_32 + tmp_29*(mi.betaX*mi.gamDDxx + mi.betaY*mi.gamDDxy + mi.betaZ*mi.gamDDxz));
      con[iSy] = tmp_33*(mi.gamDDxy*tmp_30 + mi.gamDDyy*tmp_31 + mi.gamDDyz*tmp_32 + tmp_29*(mi.betaX*mi.gamDDxy + mi.betaY*mi.gamDDyy + mi.betaZ*mi.gamDDyz));
      con[iSz] = tmp_33*(mi.gamDDxz*tmp_30 + mi.gamDDyz*tmp_31 + mi.gamDDzz*tmp_32 + tmp_29*(mi.betaX*mi.gamDDxz + mi.betaY*mi.gamDDyz + mi.betaZ*mi.gamDDzz));
      con[iTau] = ((mi.alpha)*(mi.alpha))*tmp_25*tmp_29 - tmp_26;
    }
    Wrote to file "NRPY+detg.h"
    /*
     *  Original SymPy expression:
     *  "detg = mi.alpha*sqrt(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy)"
     */
    {
      detg = mi.alpha*sqrt(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy);
    }
    Wrote to file "NRPY+gamUU.h"
    /*
     *  Original SymPy expressions:
     *  "[gamUUxx = (mi.gamDDyy*mi.gamDDzz - mi.gamDDyz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy),
     *    gamUUxy = (-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy),
     *    gamUUxz = (mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy),
     *    gamUUyy = (mi.gamDDxx*mi.gamDDzz - mi.gamDDxz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy),
     *    gamUUyz = (-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy),
     *    gamUUzz = (mi.gamDDxx*mi.gamDDyy - mi.gamDDxy**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy)]"
     */
    {
      const double tmp_5 = (1.0/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy));
      gamUUxx = tmp_5*(mi.gamDDyy*mi.gamDDzz - ((mi.gamDDyz)*(mi.gamDDyz)));
      gamUUxy = tmp_5*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz);
      gamUUxz = tmp_5*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy);
      gamUUyy = tmp_5*(mi.gamDDxx*mi.gamDDzz - ((mi.gamDDxz)*(mi.gamDDxz)));
      gamUUyz = tmp_5*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz);
      gamUUzz = tmp_5*(mi.gamDDxx*mi.gamDDyy - ((mi.gamDDxy)*(mi.gamDDxy)));
    }


<a id='flux'></a>

# Step 3: Compute the flux
$$\label{flux}$$

The fluxes are as follows
\begin{equation}
\frac{\partial}{\partial t} 
\begin{pmatrix}
\rhostar\\
\Svectilde\\
\tautilde
\end{pmatrix} + \frac{\partial}{\partial x^j}\begin{pmatrix} \rhostar v^j\\
\lapse\rtgamma T^j_i\\ \lapse^2\rtgamma T^{0j} - \rhostar v^j
\end{pmatrix}  = \begin{pmatrix} 0 \\ \frac 1 2 \lapse\rtgamma T^{\alpha\beta}g_{\alpha\beta,i} \\ s \end{pmatrix}
\end{equation}
so the flux is 
\begin{equation}
\mathcal{F} = \begin{pmatrix} \rhostar v^i \\ \lapse\rtgamma T^i_k \\ \lapse^2\rtgamma T^{0i} - \rhostar v^i
\end{pmatrix}
\end{equation}
In the moving-mesh formalism, the flux is just taken along the x directions so we have
\begin{equation}
\mathcal{F} = \begin{pmatrix} \rhostar v^1 \\ \lapse\rtgamma T^1_k \\ \lapse^2\rtgamma T^{01} - \rhostar v^1
\end{pmatrix}
\end{equation}
Note that we will need to rotate $T^{\mu\nu}$ and $g_{\mu\nu}$ to get the right orientation.
In order to do this, we must first compute the stress energy tensor:
\begin{equation}
T^{\mu\nu} = \rho h u^{\mu}u^{\nu} + Pg^{\mu\nu} = \rho h (u^0)^2v^iv^j + P g^{\mu\nu}
\end{equation}



```python
# Next compute fluxes of conservative variables
Ge.compute_rho_star_fluxU(                            vU,        Ge.rho_star)
Ge.compute_tau_tilde_fluxU(alpha, Ge.sqrtgammaDET,    vU,Ge.T4UU,Ge.rho_star)
Ge.compute_S_tilde_fluxUD( alpha, Ge.sqrtgammaDET,       Ge.T4UD)


normD = ixp.zerorank1()
normD[0], normD[1], normD[2] = sp.symbols("norm[0] norm[1] norm[2]", real=True)

faceVelU = ixp.zerorank1()
faceVelU[0], faceVelU[1], faceVelU[2] = sp.symbols("faceVelocity[0] faceVelocity[1] faceVelocity[2]", real=True)
# Zach says: don't forget to limit the velocities after they are computed!

faceVelNorm = sp.sympify(0)
for i in range(3) :
    faceVelNorm += normD[i]*faceVelU[i]

exprArray = []
nameArray = []
exprArray.append( Ge.rho_star)
nameArray.append( "temp_rho_star")
exprArray.append( Ge.T4UU[0][1])
nameArray.append( "temp_T4UU01")

rho_star_flux = sp.sympify(0)
for i in range(3) :
    rho_star_flux += Ge.rho_star_fluxU[i]*normD[i]
rho_star_flux -= Ge.rho_star*faceVelNorm

exprArray.append( rho_star_flux)
nameArray.append( "flux[iRhoStar]")


tau_tilde_flux = sp.sympify(0)
for i in range(3) :
    tau_tilde_flux += Ge.tau_tilde_fluxU[i]*normD[i]

tau_tilde_flux -= Ge.tau_tilde*faceVelNorm

S_tilde_fluxD = ixp.zerorank1()

for i in range(3) :
    S_tilde_fluxD[i] -= Ge.S_tildeD[i]*faceVelNorm
    for j in range(3) :
        S_tilde_fluxD[i] += Ge.S_tilde_fluxUD[j][i]*normD[j]

for j, comp in enumerate(["x","y", "z"]) :
    exprArray.append( S_tilde_fluxD[j])
    nameArray.append( "flux[iS{0}]".format(comp))

exprArray.append( tau_tilde_flux)
nameArray.append( "flux[iTau]")

#for expr, name in zip( exprArray, nameArray) :
#    print( name)
outputC(exprArray, nameArray, filename="NRPY+calFlux.h", params="outCverbose=False")

!cat NRPY+calFlux.h
```

    Wrote to file "NRPY+calFlux.h"
    {
      const double tmp_6 = mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy;
      const double tmp_7 = sqrt(tmp_6);
      const double tmp_8 = (1.0/((GAMMA_SPEED_LIMIT)*(GAMMA_SPEED_LIMIT)));
      const double tmp_11 = (1.0/((mi.alpha)*(mi.alpha)));
      const double tmp_12 = mi.betaX + vx;
      const double tmp_13 = mi.gamDDxx*tmp_11*((tmp_12)*(tmp_12));
      const double tmp_15 = mi.betaY + vy;
      const double tmp_16 = mi.gamDDyy*tmp_11*((tmp_15)*(tmp_15));
      const double tmp_18 = mi.betaZ + vz;
      const double tmp_19 = mi.gamDDzz*tmp_11*((tmp_18)*(tmp_18));
      const double tmp_22 = tmp_11*tmp_12*tmp_15;
      const double tmp_24 = mi.gamDDxz*tmp_11*tmp_12*tmp_18;
      const double tmp_25 = mi.gamDDyz*tmp_11*tmp_15*tmp_18;
      const double tmp_26 = 2*mi.gamDDxy*tmp_22;
      const double tmp_29 = (1.0/2.0)*fabs(-tmp_13 - tmp_16 - tmp_19 - 2*tmp_24 - 2*tmp_25 - tmp_26 - tmp_8 + 1);
      const double tmp_30 = -mi.gamDDxy*tmp_22 - 1.0/2.0*tmp_13 - 1.0/2.0*tmp_16 - 1.0/2.0*tmp_19 - tmp_24 - tmp_25 + tmp_29 + (1.0/2.0)*tmp_8 + 1.0/2.0;
      const double tmp_31 = rho*tmp_7/sqrt(tmp_30);
      const double tmp_32 = p*tmp_11;
      const double tmp_33 = sqrt((mi.gamDDxy*tmp_22 + (1.0/2.0)*tmp_13 + (1.0/2.0)*tmp_16 + (1.0/2.0)*tmp_19 + tmp_24 + tmp_25 - tmp_29 - 1.0/2.0*tmp_8 + 1.0/2.0)/(TINYDOUBLE + tmp_13 + tmp_16 + tmp_19 + 2*tmp_24 + 2*tmp_25 + tmp_26));
      const double tmp_34 = -mi.betaX + tmp_12*tmp_33;
      const double tmp_35 = rho*tmp_11*(ie + p/rho + 1)/tmp_30;
      const double tmp_36 = tmp_34*tmp_35;
      const double tmp_37 = mi.betaX*tmp_32 + tmp_36;
      const double tmp_41 = faceVelocity[0]*norm[0] + faceVelocity[1]*norm[1] + faceVelocity[2]*norm[2];
      const double tmp_42 = mi.betaX*mi.gamDDxx + mi.betaY*mi.gamDDxy + mi.betaZ*mi.gamDDxz;
      const double tmp_43 = -tmp_32 + tmp_35;
      const double tmp_44 = -mi.betaY + tmp_15*tmp_33;
      const double tmp_46 = mi.betaY*tmp_32 + tmp_35*tmp_44;
      const double tmp_47 = -mi.betaZ + tmp_18*tmp_33;
      const double tmp_48 = mi.betaZ*tmp_32 + tmp_35*tmp_47;
      const double tmp_49 = mi.alpha*tmp_7;
      const double tmp_50 = tmp_41*tmp_49;
      const double tmp_51 = (1.0/(tmp_6));
      const double tmp_52 = p*(-((mi.betaX)*(mi.betaX))*tmp_11 + tmp_51*(mi.gamDDyy*mi.gamDDzz - ((mi.gamDDyz)*(mi.gamDDyz)))) + ((tmp_34)*(tmp_34))*tmp_35;
      const double tmp_54 = p*(-mi.betaX*mi.betaY*tmp_11 + tmp_51*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz)) + tmp_36*tmp_44;
      const double tmp_56 = p*(-mi.betaX*mi.betaZ*tmp_11 + tmp_51*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy)) + tmp_36*tmp_47;
      const double tmp_58 = norm[0]*tmp_49;
      const double tmp_59 = p*(-((mi.betaY)*(mi.betaY))*tmp_11 + tmp_51*(mi.gamDDxx*mi.gamDDzz - ((mi.gamDDxz)*(mi.gamDDxz)))) + tmp_35*((tmp_44)*(tmp_44));
      const double tmp_60 = p*(-mi.betaY*mi.betaZ*tmp_11 + tmp_51*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz)) + tmp_35*tmp_44*tmp_47;
      const double tmp_61 = norm[1]*tmp_49;
      const double tmp_62 = p*(-((mi.betaZ)*(mi.betaZ))*tmp_11 + tmp_51*(mi.gamDDxx*mi.gamDDyy - ((mi.gamDDxy)*(mi.gamDDxy)))) + tmp_35*((tmp_47)*(tmp_47));
      const double tmp_63 = norm[2]*tmp_49;
      const double tmp_64 = mi.betaX*mi.gamDDxy + mi.betaY*mi.gamDDyy + mi.betaZ*mi.gamDDyz;
      const double tmp_66 = mi.betaX*mi.gamDDxz + mi.betaY*mi.gamDDyz + mi.betaZ*mi.gamDDzz;
      const double tmp_67 = ((mi.alpha)*(mi.alpha))*tmp_7;
      temp_rho_star = tmp_31;
      temp_T4UU01 = tmp_37;
      flux[iRhoStar] = norm[0]*tmp_31*vx + norm[1]*tmp_31*vy + norm[2]*tmp_31*vz - tmp_31*tmp_41;
      flux[iSx] = -tmp_50*(mi.gamDDxx*tmp_37 + mi.gamDDxy*tmp_46 + mi.gamDDxz*tmp_48 + tmp_42*tmp_43) + tmp_58*(mi.gamDDxx*tmp_52 + mi.gamDDxy*tmp_54 + mi.gamDDxz*tmp_56 + tmp_37*tmp_42) + tmp_61*(mi.gamDDxx*tmp_54 + mi.gamDDxy*tmp_59 + mi.gamDDxz*tmp_60 + tmp_42*tmp_46) + tmp_63*(mi.gamDDxx*tmp_56 + mi.gamDDxy*tmp_60 + mi.gamDDxz*tmp_62 + tmp_42*tmp_48);
      flux[iSy] = -tmp_50*(mi.gamDDxy*tmp_37 + mi.gamDDyy*tmp_46 + mi.gamDDyz*tmp_48 + tmp_43*tmp_64) + tmp_58*(mi.gamDDxy*tmp_52 + mi.gamDDyy*tmp_54 + mi.gamDDyz*tmp_56 + tmp_37*tmp_64) + tmp_61*(mi.gamDDxy*tmp_54 + mi.gamDDyy*tmp_59 + mi.gamDDyz*tmp_60 + tmp_46*tmp_64) + tmp_63*(mi.gamDDxy*tmp_56 + mi.gamDDyy*tmp_60 + mi.gamDDyz*tmp_62 + tmp_48*tmp_64);
      flux[iSz] = -tmp_50*(mi.gamDDxz*tmp_37 + mi.gamDDyz*tmp_46 + mi.gamDDzz*tmp_48 + tmp_43*tmp_66) + tmp_58*(mi.gamDDxz*tmp_52 + mi.gamDDyz*tmp_54 + mi.gamDDzz*tmp_56 + tmp_37*tmp_66) + tmp_61*(mi.gamDDxz*tmp_54 + mi.gamDDyz*tmp_59 + mi.gamDDzz*tmp_60 + tmp_46*tmp_66) + tmp_63*(mi.gamDDxz*tmp_56 + mi.gamDDyz*tmp_60 + mi.gamDDzz*tmp_62 + tmp_48*tmp_66);
      flux[iTau] = norm[0]*(-tmp_31*vx + tmp_37*tmp_67) + norm[1]*(-tmp_31*vy + tmp_46*tmp_67) + norm[2]*(-tmp_31*vz + tmp_48*tmp_67) - tmp_41*(-tmp_31 + tmp_43*tmp_67);
    }


<a id='source'></a>

# Step 4: Source Terms
$$\label{source}$$

The sources terms are for mass, momentum and energy are: 
\begin{equation}
\source = (0, \frac 1 2 \lapse\rtgamma \T{\alpha}{\beta}g_{\alpha\beta,i}, s),
\end{equation}
For a time stationary metric $s\neq 0$, so we will ignore this until the next section.  As for the rest, we need to define derivatives of the metric.  Suppose I have done this already.  Then the code for the source terms is:



```python
# FIXME: Assume static spacetime with KDD = betaU = betaU_dD = 0
KDD = ixp.zerorank2()
betaU = ixp.zerorank1()
betaU_dD = ixp.zerorank2()

# Set+evaluate derivatives of alpha, performing 2nd-order finite difference
alpha_dD = ixp.zerorank1()
h = sp.symbols("h")
for i in range(3) :
    alpha_plus, alpha_minus = sp.symbols("mi_plus[{0}].alpha mi_minus[{0}].alpha".format(i))
    alpha_dD[i] = (alpha_plus - alpha_minus)/(2*h)

# Set+evaluate derivatives of gamma_{ij}, performing 2nd-order finite difference
gammaDD_dD = ixp.zerorank3()
components = ["xx", "xy", "xz", "yy", "yz", "zz"]
for i in range(3) :
    names_plus = ""
    names_minus = ""
    for comp in components :
        names_plus  = names_plus  +  "mi_plus[{0}].gamDD{1} ".format(i, comp)
        names_minus = names_minus + "mi_minus[{0}].gamDD{1} ".format(i, comp)

    gxx_plus, gxy_plus, gxz_plus, gyy_plus, gyz_plus, gzz_plus = sp.symbols( names_plus)
    gxx_minus, gxy_minus, gxz_minus, gyy_minus, gyz_minus, gzz_minus = sp.symbols( names_minus)

    gammaDD_dD[0][0][i] = (gxx_plus - gxx_minus)/(2*h)
    gammaDD_dD[0][1][i] = (gxy_plus - gxy_minus)/(2*h)
    gammaDD_dD[0][2][i] = (gxz_plus - gxz_minus)/(2*h)
    gammaDD_dD[1][0][i] = (gxy_plus - gxy_minus)/(2*h)
    gammaDD_dD[1][1][i] = (gyy_plus - gyy_minus)/(2*h)
    gammaDD_dD[1][2][i] = (gyz_plus - gyz_minus)/(2*h)
    gammaDD_dD[2][0][i] = (gxz_plus - gxz_minus)/(2*h)
    gammaDD_dD[2][1][i] = (gyz_plus - gyz_minus)/(2*h)
    gammaDD_dD[2][2][i] = (gzz_plus - gzz_minus)/(2*h)

# Compute g_{mu nu, i} based on ADM quantities & derivatives defined above
Ge.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

# Compute source terms for tau tilde & S tilde:
Ge.compute_s_source_term(KDD,betaU,alpha, Ge.sqrtgammaDET,alpha_dD,  Ge.T4UU)
Ge.compute_S_tilde_source_termD(   alpha, Ge.sqrtgammaDET,Ge.g4DD_zerotimederiv_dD, Ge.T4UU)

exprArray = []
nameArray = []

#momentum terms
for i in range(3) :
    exprArray.append( Ge.S_tilde_source_termD[i])
    nameArray.append( "vSource[{0}]".format(i))

#tau term
exprArray.append( Ge.s_source_term)
nameArray.append( "eSource")
outputC( exprArray, nameArray, filename="NRPY+calSources.h", params="outCverbose=False")

!cat NRPY+calSources.h
```

    Wrote to file "NRPY+calSources.h"
    {
      const double tmp_0 = (1.0/(h));
      const double tmp_1 = (1.0/2.0)*tmp_0;
      const double tmp_2 = tmp_1*(-mi_minus[0].alpha + mi_plus[0].alpha);
      const double tmp_9 = mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy;
      const double tmp_10 = sqrt(tmp_9);
      const double tmp_12 = (1.0/((mi.alpha)*(mi.alpha)));
      const double tmp_13 = p*tmp_12;
      const double tmp_14 = (1.0/((GAMMA_SPEED_LIMIT)*(GAMMA_SPEED_LIMIT)));
      const double tmp_16 = mi.betaX + vx;
      const double tmp_17 = mi.gamDDxx*tmp_12*((tmp_16)*(tmp_16));
      const double tmp_19 = mi.betaY + vy;
      const double tmp_20 = mi.gamDDyy*tmp_12*((tmp_19)*(tmp_19));
      const double tmp_22 = mi.betaZ + vz;
      const double tmp_23 = mi.gamDDzz*tmp_12*((tmp_22)*(tmp_22));
      const double tmp_26 = tmp_12*tmp_16*tmp_19;
      const double tmp_28 = mi.gamDDxz*tmp_12*tmp_16*tmp_22;
      const double tmp_29 = mi.gamDDyz*tmp_12*tmp_19*tmp_22;
      const double tmp_30 = 2*mi.gamDDxy*tmp_26;
      const double tmp_33 = (1.0/2.0)*fabs(-tmp_14 - tmp_17 - tmp_20 - tmp_23 - 2*tmp_28 - 2*tmp_29 - tmp_30 + 1);
      const double tmp_34 = rho*tmp_12*(ie + p/rho + 1)/(-mi.gamDDxy*tmp_26 + (1.0/2.0)*tmp_14 - 1.0/2.0*tmp_17 - 1.0/2.0*tmp_20 - 1.0/2.0*tmp_23 - tmp_28 - tmp_29 + tmp_33 + 1.0/2.0);
      const double tmp_35 = ((mi.alpha)*(mi.alpha))*tmp_10*(-tmp_13 + tmp_34);
      const double tmp_36 = (1.0/(tmp_9));
      const double tmp_37 = sqrt((mi.gamDDxy*tmp_26 - 1.0/2.0*tmp_14 + (1.0/2.0)*tmp_17 + (1.0/2.0)*tmp_20 + (1.0/2.0)*tmp_23 + tmp_28 + tmp_29 - tmp_33 + 1.0/2.0)/(TINYDOUBLE + tmp_17 + tmp_20 + tmp_23 + 2*tmp_28 + 2*tmp_29 + tmp_30));
      const double tmp_38 = -mi.betaX + tmp_16*tmp_37;
      const double tmp_39 = mi.alpha*tmp_10;
      const double tmp_40 = (1.0/4.0)*tmp_0*tmp_39;
      const double tmp_41 = tmp_40*(p*(-((mi.betaX)*(mi.betaX))*tmp_12 + tmp_36*(mi.gamDDyy*mi.gamDDzz - ((mi.gamDDyz)*(mi.gamDDyz)))) + tmp_34*((tmp_38)*(tmp_38)));
      const double tmp_42 = -mi.betaY + tmp_19*tmp_37;
      const double tmp_43 = tmp_40*(p*(-((mi.betaY)*(mi.betaY))*tmp_12 + tmp_36*(mi.gamDDxx*mi.gamDDzz - ((mi.gamDDxz)*(mi.gamDDxz)))) + tmp_34*((tmp_42)*(tmp_42)));
      const double tmp_44 = -mi.betaZ + tmp_22*tmp_37;
      const double tmp_45 = tmp_40*(p*(-((mi.betaZ)*(mi.betaZ))*tmp_12 + tmp_36*(mi.gamDDxx*mi.gamDDyy - ((mi.gamDDxy)*(mi.gamDDxy)))) + tmp_34*((tmp_44)*(tmp_44)));
      const double tmp_47 = tmp_34*tmp_38;
      const double tmp_48 = tmp_1*tmp_39;
      const double tmp_49 = tmp_48*(p*(-mi.betaX*mi.betaY*tmp_12 + tmp_36*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz)) + tmp_42*tmp_47);
      const double tmp_50 = tmp_48*(p*(-mi.betaX*mi.betaZ*tmp_12 + tmp_36*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy)) + tmp_44*tmp_47);
      const double tmp_52 = tmp_48*(p*(-mi.betaY*mi.betaZ*tmp_12 + tmp_36*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz)) + tmp_34*tmp_42*tmp_44);
      const double tmp_53 = tmp_1*(-mi_minus[1].alpha + mi_plus[1].alpha);
      const double tmp_54 = tmp_1*(-mi_minus[2].alpha + mi_plus[2].alpha);
      vSource[0] = -tmp_2*tmp_35 + tmp_41*(-mi_minus[0].gamDDxx + mi_plus[0].gamDDxx) + tmp_43*(-mi_minus[0].gamDDyy + mi_plus[0].gamDDyy) + tmp_45*(-mi_minus[0].gamDDzz + mi_plus[0].gamDDzz) + tmp_49*(-mi_minus[0].gamDDxy + mi_plus[0].gamDDxy) + tmp_50*(-mi_minus[0].gamDDxz + mi_plus[0].gamDDxz) + tmp_52*(-mi_minus[0].gamDDyz + mi_plus[0].gamDDyz);
      vSource[1] = -tmp_35*tmp_53 + tmp_41*(-mi_minus[1].gamDDxx + mi_plus[1].gamDDxx) + tmp_43*(-mi_minus[1].gamDDyy + mi_plus[1].gamDDyy) + tmp_45*(-mi_minus[1].gamDDzz + mi_plus[1].gamDDzz) + tmp_49*(-mi_minus[1].gamDDxy + mi_plus[1].gamDDxy) + tmp_50*(-mi_minus[1].gamDDxz + mi_plus[1].gamDDxz) + tmp_52*(-mi_minus[1].gamDDyz + mi_plus[1].gamDDyz);
      vSource[2] = -tmp_35*tmp_54 + tmp_41*(-mi_minus[2].gamDDxx + mi_plus[2].gamDDxx) + tmp_43*(-mi_minus[2].gamDDyy + mi_plus[2].gamDDyy) + tmp_45*(-mi_minus[2].gamDDzz + mi_plus[2].gamDDzz) + tmp_49*(-mi_minus[2].gamDDxy + mi_plus[2].gamDDxy) + tmp_50*(-mi_minus[2].gamDDxz + mi_plus[2].gamDDxz) + tmp_52*(-mi_minus[2].gamDDyz + mi_plus[2].gamDDyz);
      eSource = tmp_39*(tmp_2*(-mi.betaX*tmp_13 - tmp_47) + tmp_53*(-mi.betaY*tmp_13 - tmp_34*tmp_42) + tmp_54*(-mi.betaZ*tmp_13 - tmp_34*tmp_44));
    }


<a id='solver'></a>

# Step 6: Conservative to Primitive Solver
$$\label{solver}$$

We now discuss reverse mapping from conservative to primitive variables.
Given the lapse, shift vector, and $\rtgamma$, the mapping between primitive and conserved variables is straightforward.  However, the reverse is not as simple.  In GRMHD, the conservative to primitive solver is amplified by the inclusion of the magnetic field, leading to rather sophisticated root-finding strategies.  The failure rates of these algorithms are low (??), but since this algorithm may be executed several times per timestep for every gridpoint, even a low failure can give unacceptable collective failure rates.  However, for purely polytropic equations of state, e.g., $P\propto\rho^{\Gamma_1}$, the conservative to primitive variable solver is greatly simplified.  

To construct the conservative-to-primitive variable solver, we restrict ourselves to polytropic equations of states
\begin{equation}
P = P_0\left(\frac{\rho}{\rho_0}\right)^{\Gamma_1} \quad\textrm{and}\quad \epsilon = \epsilon_0\left(\frac{\rho}{\rho_0}\right)^{\Gamma_1-1},
\end{equation}
where $P_0$, $\rho_0$, and $\epsilon_0$ are the fiducial pressure, density, and internal energy, and we have used the relation $P = (\Gamma_1 - 1)\rho\epsilon$.  

For such a polytropic equation of state, the energy equation is redundant and effectively we are only concerned with the continuity and momentum equations. The conservative variables of concern are $\rhostar$ and $\Svectilde$.  Noting that the shift, $\alpha$, and $\rtgamma$ are provided by the Einsteins field equation solver, we can write
\begin{equation}
u^0 = \frac{\rhostar}{\alpha\rtgamma\rho} = u^0(\rho)  \quad\textrm{and}\quad \uvec = \frac{\Svectilde}{\alpha\rtgamma\rho h} = \uvec(\rho).
\end{equation}
Noting that the four-velocity $u^2 = g_{\mu\nu}u^{\mu}u^{\nu} = g^{00}u^0u^0 + 2g^{0i}u^0\uvec^i + g_{ij}\uvec^i\uvec^j = -1$, we have
\begin{equation}
 0 = f(\rho)\equiv \alpha^2\gamma\rho^2h^2 + \left(-\lapse^2 + \shift\cdot\shift\right)\rhostar^2h^2 + 2h\rhostar\shift\cdot\Svectilde + \Svectilde\cdot\Svectilde,
\end{equation}
which is an implicit equation of either $\rho$ or $u^0$, where $h(\rho = \rhostar/(\alpha\rtgamma u^0)) = 1 + \gamma_1 \epsilon$ which can be inverted by standard nonlinear root-finding algorithms, e.g., Newton-Raphson. 

We put this all together to define a function, $f(\rho)$, whose root is zero that we will find via Newton-Raphson.  

Several checks must be performed:

1. $\rhostar > 0$ : This check is performed at the very beginning

2. $\rho > \rho_{\rm min}$ : This check is performed after the fact

3. $u_0 < \alpha^{-1}\Gamma_{\rm max}$ : This check is performed after the fact as well


```python
DIM = 3
# Declare rank-1 contravariant ("v") vector
vU = ixp.declarerank1("vU")
shiftU = ixp.zerorank1()
rho, gamma1 = sp.symbols("rho gamma")
Sx, Sy, Sz = sp.symbols("con[iSx] con[iSy] con[iSz]")
p0, rho0, rhostar = sp.symbols("p_0 rho_0 rhostar")
# Declare rank-2 covariant gmunu
#gammaDD = ixp.declarerank2("gammaDD","sym01")
StildeD = ixp.declarerank1("StildeD")
lapse, beta_x, beta_y, beta_z = sp.symbols( "mi.alpha mi.betaX mi.betaY mi.betaZ")
rtgamma = Ge.sqrtgammaDET

shiftU[0] = beta_x
shiftU[1] = beta_y
shiftU[2] = beta_z

StildeD[0] = Sx
StildeD[1] = Sy
StildeD[2] = Sz

# gamma = rtgamma*rtgamma <- unused
lapse2 = lapse*lapse
uU0 = rhostar/(lapse*rtgamma*rho)
epsilon = p0/rho0*(rho/rho0)**(gamma1 - 1)/(gamma1 - 1)
h = 1 + gamma1*epsilon

beta2 = 0.

for i in range(DIM) :
    for j in range(DIM) :
        beta2 += gammaDD[i][j] * shiftU[i]*shiftU[j]

betaDotStilde = 0
for i in range(DIM) :
    betaDotStilde += shiftU[i]*StildeD[i]

# Note that this is |Stilde|^2, where the absolute value denotes
#    that this is not a proper contraction of Stilde_i, as
#    Stilde^i is NOT equal to gamma^{ij} Stilde_j (to understand
#    why this is, notice that Stilde_i is proportional to the
#    *4D* stress-energy tensor.)
Stilde2 = 0
for i in range(DIM) :
    for j in range(DIM) :
        Stilde2 += gammaUU[i][j] * StildeD[i]*StildeD[j]

f = rhostar**2*h**2 + (-lapse2 + beta2)*rhostar**2.*h**2.*uU0**2 + 2.*h*rhostar*betaDotStilde*uU0 + Stilde2

outputC(f,"rootRho",filename="NRPY+rhoRoot.h")
outputC(Stilde2, "Stilde2", filename="NRPY+Stilde2.h")
!cat NRPY+rhoRoot.h
!cat NRPY+Stilde2.h
```

    Wrote to file "NRPY+rhoRoot.h"
    Wrote to file "NRPY+Stilde2.h"
    /*
     *  Original SymPy expression:
     *  "rootRho = con[iSx]**2*(mi.gamDDyy*mi.gamDDzz - mi.gamDDyz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSx]*con[iSy]*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSx]*con[iSz]*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + con[iSy]**2*(mi.gamDDxx*mi.gamDDzz - mi.gamDDxz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSy]*con[iSz]*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + con[iSz]**2*(mi.gamDDxx*mi.gamDDyy - mi.gamDDxy**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + rhostar**2*(gamma*p_0*(rho/rho_0)**(gamma - 1)/(rho_0*(gamma - 1)) + 1)**2 + rhostar**2*(2.0*gamma*p_0*(rho/rho_0)**(gamma - 1)/(rho_0*(gamma - 1)) + 2.0)*(con[iSx]*mi.betaX + con[iSy]*mi.betaY + con[iSz]*mi.betaZ)/(mi.alpha*rho*sqrt(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy)) + rhostar**4.0*(gamma*p_0*(rho/rho_0)**(gamma - 1)/(rho_0*(gamma - 1)) + 1)**2.0*(-mi.alpha**2 + mi.betaX**2*mi.gamDDxx + 2*mi.betaX*mi.betaY*mi.gamDDxy + 2*mi.betaX*mi.betaZ*mi.gamDDxz + mi.betaY**2*mi.gamDDyy + 2*mi.betaY*mi.betaZ*mi.gamDDyz + mi.betaZ**2*mi.gamDDzz)/(mi.alpha**2*rho**2*(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy))"
     */
    {
      const double tmp_1 = (1.0/(rho_0));
      const double tmp_3 = gamma*p_0*tmp_1*pow(rho*tmp_1, gamma - 1)/(gamma - 1);
      const double tmp_11 = mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy;
      const double tmp_12 = (1.0/(tmp_11));
      const double tmp_13 = 2*con[iSx]*tmp_12;
      rootRho = ((con[iSx])*(con[iSx]))*tmp_12*(mi.gamDDyy*mi.gamDDzz - ((mi.gamDDyz)*(mi.gamDDyz))) + ((con[iSy])*(con[iSy]))*tmp_12*(mi.gamDDxx*mi.gamDDzz - ((mi.gamDDxz)*(mi.gamDDxz))) + 2*con[iSy]*con[iSz]*tmp_12*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz) + con[iSy]*tmp_13*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz) + ((con[iSz])*(con[iSz]))*tmp_12*(mi.gamDDxx*mi.gamDDyy - ((mi.gamDDxy)*(mi.gamDDxy))) + con[iSz]*tmp_13*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy) + ((rhostar)*(rhostar))*((tmp_3 + 1)*(tmp_3 + 1)) + ((rhostar)*(rhostar))*(2.0*tmp_3 + 2.0)*(con[iSx]*mi.betaX + con[iSy]*mi.betaY + con[iSz]*mi.betaZ)/(mi.alpha*rho*sqrt(tmp_11)) + ((rhostar)*(rhostar)*(rhostar)*(rhostar))*tmp_12*((tmp_3 + 1)*(tmp_3 + 1))*(-((mi.alpha)*(mi.alpha)) + ((mi.betaX)*(mi.betaX))*mi.gamDDxx + 2*mi.betaX*mi.betaY*mi.gamDDxy + 2*mi.betaX*mi.betaZ*mi.gamDDxz + ((mi.betaY)*(mi.betaY))*mi.gamDDyy + 2*mi.betaY*mi.betaZ*mi.gamDDyz + ((mi.betaZ)*(mi.betaZ))*mi.gamDDzz)/(((mi.alpha)*(mi.alpha))*((rho)*(rho)));
    }
    /*
     *  Original SymPy expression:
     *  "Stilde2 = con[iSx]**2*(mi.gamDDyy*mi.gamDDzz - mi.gamDDyz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSx]*con[iSy]*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSx]*con[iSz]*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + con[iSy]**2*(mi.gamDDxx*mi.gamDDzz - mi.gamDDxz**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + 2*con[iSy]*con[iSz]*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy) + con[iSz]**2*(mi.gamDDxx*mi.gamDDyy - mi.gamDDxy**2)/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*mi.gamDDyz**2 - mi.gamDDxy**2*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - mi.gamDDxz**2*mi.gamDDyy)"
     */
    {
      const double tmp_5 = (1.0/(mi.gamDDxx*mi.gamDDyy*mi.gamDDzz - mi.gamDDxx*((mi.gamDDyz)*(mi.gamDDyz)) - ((mi.gamDDxy)*(mi.gamDDxy))*mi.gamDDzz + 2*mi.gamDDxy*mi.gamDDxz*mi.gamDDyz - ((mi.gamDDxz)*(mi.gamDDxz))*mi.gamDDyy));
      const double tmp_6 = 2*con[iSx]*tmp_5;
      Stilde2 = ((con[iSx])*(con[iSx]))*tmp_5*(mi.gamDDyy*mi.gamDDzz - ((mi.gamDDyz)*(mi.gamDDyz))) + ((con[iSy])*(con[iSy]))*tmp_5*(mi.gamDDxx*mi.gamDDzz - ((mi.gamDDxz)*(mi.gamDDxz))) + 2*con[iSy]*con[iSz]*tmp_5*(-mi.gamDDxx*mi.gamDDyz + mi.gamDDxy*mi.gamDDxz) + con[iSy]*tmp_6*(-mi.gamDDxy*mi.gamDDzz + mi.gamDDxz*mi.gamDDyz) + ((con[iSz])*(con[iSz]))*tmp_5*(mi.gamDDxx*mi.gamDDyy - ((mi.gamDDxy)*(mi.gamDDxy))) + con[iSz]*tmp_6*(mi.gamDDxy*mi.gamDDyz - mi.gamDDxz*mi.gamDDyy);
    }


The root solve above finds $\rho$, which then allows us to get 
\begin{equation}
u^0 = \frac{\rhostar}{\alpha\rtgamma\rho}\quad\textrm{and}\quad \vel = \frac{\uvec}{u^0} = \frac{\Svectilde}{\rhostar h(\rho)}.
\end{equation}
and thus we can find the rest of the primitives.


```python
#rhostar = sp.symbols("rhostar")
#StildeU = ixp.declarerank1("StildeU")
velU = ixp.zerorank1()

#lapse, rtgamma, rho, gamma1, c = sp.symbols("lapse rtgamma rho gamma1 c")
rho, rhostar = sp.symbols("testPrim[iRho] con[iRhoStar]")

u0 = rhostar/(lapse*rtgamma*rho)
epsilon = p0/rho0*(rho/rho0)**(gamma1 - 1)/(gamma1 - 1)

h = 1. + gamma1*epsilon

for i in range(DIM) :
    for j in range(DIM) :
        velU[i] += gammaUU[i][j]*StildeD[j]/(rhostar * h)/u0

outputC([h,u0,velU[0],velU[1],velU[2]], ["h", "u0","testPrim[ivx]", "testPrim[ivy]", "testPrim[ivz]"],filename="NRPY+getv.h")
```

    Wrote to file "NRPY+getv.h"


<a id='lorentz'></a>

# Step 7: Lorentz Boosts
$$\label{lorentz}$$

We need to boost to the frame of the moving face.  The boost is
\begin{equation}
B(\beta) =\begin{pmatrix}
\gamma & -\beta\gamma n_x & -\beta\gamma n_y & -\beta\gamma n_z  \\
-\beta\gamma n_x  & 1 + (\gamma-1)n_x^2 & (\gamma-1)n_x n_y & (\gamma-1)n_x n_z\\
-\beta\gamma n_x  & (\gamma-1)n_y n_x & 1 + (\gamma-1)n_y^2 & (\gamma-1)n_y n_z\\
-\beta\gamma n_x & (\gamma-1) n_z n_x & (\gamma-1)n_z n_x & 1 + (\gamma-1)n_z^2 
\end{pmatrix} 
\end{equation}
And the boost is $X' = B(\beta) X$, where $X'$ and $X$ are four vectors.

So the rest of this is straightforward.  

<a id='TmunuSph'></a>

# Step 8: Compute $T^{\mu\nu}$ in Cartesian Coordinates
$$\label{TmunuSph}$$




```python
# declare gammaDD
gammaDD = ixp.zerorank2()
components = ["xx", "xy", "xz", "yy", "yz", "zz"]
names = ""
for comp in components :
    names = names + "mi.gamDD{0} ".format(comp)

g11, g12, g13, g22, g23, g33 = sp.symbols( names)

gammaDD[0][0] = g11
gammaDD[0][1] = g12
gammaDD[0][2] = g13
gammaDD[1][0] = g12
gammaDD[1][1] = g22
gammaDD[1][2] = g23
gammaDD[2][0] = g13
gammaDD[2][1] = g23
gammaDD[2][2] = g33

#declare alpha
alpha = sp.symbols( "mi.alpha")

#declare beta
betaU = ixp.zerorank1()
for i, comp in enumerate(["X", "Y", "Z"]) :
    betaU[i] = 0 #NEED A BETTER WAY sp.symbols( "mi.beta{0}".format(comp), real=True)

#now get the primitives
rho_b, epsilon, P = sp.symbols("rho ie press")

#get the 3-velocities
vU = ixp.zerorank1()
for i, comp in enumerate( ["vx", "vy", "vz"]) :
    vU[i] = sp.symbols("{0}".format(comp))

Ge.u4U_in_terms_of_vU__rescale_vU_by_applying_speed_limit(alpha,betaU,gammaDD, vU)

u4U = Ge.u4U_ito_vU

# First compute stress-energy tensor T4UU and T4UD in Spherical Coordinates:
Ge.compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U)
Ge.compute_T4UD(gammaDD,betaU,alpha, Ge.T4UU)

outputC([Ge.T4UU[0][0],Ge.T4UU[1][1],Ge.T4UU[2][2],Ge.T4UU[3][3] ], ["T4UU_diag[0]", "T4UU_diag[1]","T4UU_diag[2]", "T4UU_diag[3]"],filename="NRPY+getT4UU.h")
```

    Wrote to file "NRPY+getT4UU.h"


<a id='latex_pdf_output'></a>

# Step 9: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-GRHD_Equations-Cartesian-c-code.pdf](Tutorial-GRHD_Equations-Cartesian-c-code.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-GRHD_Equations-Cartesian-c-code")
```

    Created Tutorial-GRHD_Equations-Cartesian-c-code.tex, and compiled LaTeX
        file to PDF file Tutorial-GRHD_Equations-Cartesian-c-code.pdf

