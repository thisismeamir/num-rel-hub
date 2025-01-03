<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Derivation of the GRMHD Evolution Equations

## Author: Samuel Cupp

## This notebook provides an in-depth derivation of the evolution equations for GRMHD variables, based on the work of [Duez et al. (2005)](https://arxiv.org/abs/astro-ph/0503420). The notebook explains the two sets of variables often used in a general relativistic setting: the primitive variables ($\mathbf{P} = \left\{\rho_b,P,v^i,B^i\right\}$) and the conservative variables ($\mathbf{C} = \left\{\tilde{\rho},\tilde{\tau},\tilde{S}_i,\tilde{B}^i\right\}$). Detailed steps are presented for deriving the evolution equations for these conservative variables, along with insights into the electromagnetic contribution to the stress-energy tensor, and the decision not to directly evolve the energy variable due to its potential inaccuracy.

## Introduction:

When considering the evolution of magnetohydrodynamics in a general relativistic setting, there are generally two sets of variables, the primitive variables

$$
\mathbf{P} = \left\{\rho_b,P,v^i,B^i\right\}
$$

and the conservative variables

$$
\mathbf{C} = \left\{\tilde{\rho},\tilde{\tau},\tilde{S}_i,\tilde{B}^i\right\}.
$$

The details of the individual variables in each set are discussed more in [Step 2](#suitable_equations). While the primitive variables are the quantities we associate with the usual description of physical properties, the conservative variables have far better numerical stability. As such, IllinoisGRMHD evolves $\mathbf{C}$, and so the following evolution equations are in terms of these variables. The derivation in this document is based on the work by Duez _et al_ (citation below), and I have rederived their results. They label their primitive density as $\rho_0$ (0 for rest mass density), while we label ours as $\rho_b$ ($b$ for baryonic density). These are the same quantity despite the different labels.

### Original Source
* M. D. Duez, Y. T. Liu, S. L. Shapiro, and B. C. Stephens. Relativistic magnetohydrodynamics in dynamical spacetimes: Numerical methods and tests. Phys. Rev. D 72, 024028 (2005). ([arxiv:astro-ph/0503420](https://arxiv.org/abs/astro-ph/0503420)).

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

This module is organized as follows

1. [Step 1](#EM_tensor): **Deriving $T_{EM}^{\mu\nu}$**
    1. [Step 1.1](#F_scalar): Solving for $F_{\alpha\beta}F^{\alpha\beta}$
    1. [Step 1.2](#F_tensor): Solving for $F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda}$
    1. [Step 1.3](#T_solve): Finishing Touches
1. [Step 2](#suitable_equations): **Evolution Equations for the Conservative Variables**
    1. [Step 2.1](#density_eqn): Evolution Equation for Fluid Density $\tilde{\rho}$
    1. [Step 2.2](#momentum_eqn): Evolution Equation for Momentum Density Variable $\tilde{S}_i$
    1. [Step 2.3](#energy_eqn): Evolution Equation for Energy Density Variable $\tilde{\tau}$
    1. [Step 2.4](#induction_eqn): Evolution Equation for the Magnetic Field $\tilde{B}^i$
    1. [Step 2.5](#summary): Summary of the Conservative Variable Evolution Equations
1. [Appendix A](#levi_civita): **Levi-Civita Contractions**
1. [Appendix B](#curv_tensor): **Extrinsic Curvature Identities**
1. [Step 3](#latex_pdf_output): **Output this notebook to $\LaTeX$-formatted PDF file**

<a id='EM_tensor'></a>

# Step 1: Deriving $T_{EM}^{\mu\nu}$ \[Back to [top](#toc)\]
$$\label{EM_tensor}$$

In order to find evolution equations for the conserved quantities for the system, we start by finding a convenient form of the electromagnetic energy-momentum tensor $T_{EM}^{\mu\nu}$, which is defined as

$$
T_{EM}^{\mu\nu} = \frac{1}{4\pi}\left( F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda} - \frac{1}{4}g^{\mu\nu}F_{\alpha\beta}F^{\alpha\beta} \right)
$$

where $F^{\mu\nu}$ is the Faraday tensor. We want to examine this with respect to an arbitrary observer with normalized four-velocity $\xi^\mu$. If we consider the electromagnetic fields in the rest frame of this observer,

$$
\xi_\mu E^\mu = 0 = \xi_\mu B^\mu
$$

since $E$ and $B$ are purely spatial. The Faraday tensor can then be decomposed in terms of $E$ and $B$:

$$
F^{\mu\nu} = \xi^\mu E^\nu - \xi^\nu E^\mu + \xi_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta}.
$$

##TODO: All the EM derivation stuff has a lot in common with Terrence's section. Further review and combination is probably necessary.

<a id='F_scalar'></a>

## Step 1.1: Solving for $F_{\alpha\beta}F^{\alpha\beta}$ \[Back to [top](#toc)\]
$$\label{F_scalar}$$

Now, we simply want to find an expression for of $T_{EM}^{\mu\nu}$ in terms of $\xi$, $E$, and $B$. Let's first consider the relatively simpler term

\begin{align}
F_{\alpha\beta}F^{\alpha\beta} &= \left(\xi_\alpha E_\beta - \xi_\beta E_\alpha 
  + \xi_\gamma {\epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta}} B_{\delta}\right)
  \left(\xi^\alpha E^\beta - \xi^\beta E^\alpha
  + \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu}\right) \\
%
&= \xi_\alpha \xi^\alpha E_\beta E^\beta - \xi_\alpha E^\alpha E_\beta \xi^\beta - E_\alpha \xi^\alpha \xi_\beta E^\beta + \xi_\beta \xi^\beta E_\alpha E^\alpha 
  + \xi_\alpha E_\beta \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu} - \xi_\beta E_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} B_{\mu} \\
&\ \ \ \ \ + \xi^\alpha E^\beta \xi_\gamma {\epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta}} B_{\delta} - \xi^\beta E^\alpha \xi_\gamma {\epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta}} B_{\delta} +  \xi_\sigma \xi_\gamma {\epsilon^{\gamma\hphantom{\alpha\beta}\delta}_{\hphantom{\gamma}\alpha\beta}} \epsilon^{\sigma\alpha\beta\mu} B_{\delta} B_{\mu} \\
%
&= -2E_\beta E^\beta + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
  + \xi_\sigma \xi^\gamma \epsilon_{\gamma\alpha\beta\delta} \epsilon^{\sigma\alpha\beta\mu} B^{\delta} B_{\mu}
\end{align}

where we have used $\xi^\mu \xi_\mu=-1$ and the permutation rules of the Levi-Civita tensor. Further simplification requires abusing the properties of the Levi-Civita tensor, which is given in [Appendix A](#levi_civita):

\begin{align}
F_{\alpha\beta}F^{\alpha\beta} &= -2E_\beta E^\beta + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
  + 2\xi_\sigma \xi^\gamma (\delta^\sigma_\delta \delta^\mu_\gamma - \delta^\sigma_\gamma \delta^\mu_\delta) B^{\delta} B_{\mu} \\
%
&= -2E_\beta E^\beta + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
  + 2\xi_\sigma \xi^\mu B^{\sigma} B_{\mu} - 2\xi_\sigma \xi^\sigma B^{\mu} B_{\mu} \\
%
&= -2E_\beta E^\beta + 4\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} E_\beta B_{\mu}
  + 2B^{\mu} B_{\mu}
\end{align}

Now, since we are in the rest frame of the observer, the spatial components of the observer's velocity are zero. Then, 
\begin{align}
\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} &= \xi_0 \xi_\sigma \epsilon^{\sigma 0\beta\mu} \\
&= \xi_0^2 \epsilon^{0 0\beta\mu} \\
&= 0
\end{align}

since any components of the Levi-Civita tensor with repeated indices are zero. Therefore,

$$
F_{\alpha\beta}F^{\alpha\beta} = -2E^\beta E_\beta + 2B^{\mu} B_{\mu}
$$

<a id='F_tensor'></a>

## Step 1.2: Solving for $F^{\mu\lambda}F^{\nu}_{\ \ \lambda}$ \[Back to [top](#toc)\]
$$\label{F_tensor}$$

In this section, we tackle the term $F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda}$ in the energy-momentum tensor.

\begin{align}
F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda} &= \left(\xi^\mu E^\lambda - \xi^\lambda E^\mu + \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta}\right)
  \left(\xi^\nu E_\lambda - \xi_\lambda E^\nu + \xi_\alpha \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\beta}\right) \\
%
&= \xi^\mu E^\lambda \xi^\nu E_\lambda - \xi^\mu E^\lambda \xi_\lambda E^\nu -\xi^\lambda E^\mu \xi^\nu E_\lambda + \xi^\lambda E^\mu \xi_\lambda E^\nu \\
&\ \ \ \ \ + \xi^\mu E_\lambda \xi_\alpha \epsilon^{\alpha\nu\lambda\beta} B_{\beta} - \xi_\lambda E^\mu \xi_\alpha \epsilon^{\alpha\nu\lambda\beta} B_{\beta}
+ \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi^\nu E_\lambda \\
&\ \ \ \ \ - \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi_\lambda E^\nu
+ \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} B_{\delta} \xi_\alpha \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\beta} \\
%
&= \xi^\mu \xi^\nu E^\lambda E_\lambda - E^\mu E^\nu +  \xi_\alpha E_\lambda B_{\beta} (\xi^\mu \epsilon^{\alpha\nu\lambda\beta} + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}) \\
&\ \ \ \ \ - \xi_\lambda \xi_\alpha B_{\beta} (\epsilon^{\alpha\nu\lambda\beta} E^\mu + \epsilon^{\alpha\mu\lambda\beta} E^\nu) + \xi_\alpha \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\delta} B_{\beta}
\end{align}

where we have again used the relations $\xi^\mu \xi_\mu=-1$ and $\xi_\mu E^\mu = 0 = \xi_\mu B^\mu$. We've also taken the liberty to relabel dummy indices to simplify the resulting expression. To simplify the final term, we must (unfortunately) once again rely on the generalized Kronecker delta. This time we use the second (far more complicated) expression in [Appendix A](#levi_civita):

\begin{align}
\xi_\alpha \xi_\gamma \epsilon^{\gamma\mu\lambda\delta} \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} B_{\delta} B_{\beta} &= -g^{\nu\phi} \xi^\alpha \xi_\gamma B_{\delta} B^{\beta} \left( \delta^\gamma_\alpha (\delta^\mu_\phi \delta^\delta_\beta - \delta^\mu_\beta \delta^\delta_\phi)
  - \delta^\gamma_\phi (\delta^\mu_\alpha \delta^\delta_\beta - \delta^\mu_\beta \delta^\delta_\alpha)
  + \delta^\gamma_\beta (\delta^\mu_\alpha  \delta^\delta_\phi - \delta^\mu_\phi \delta^\delta_\alpha) \right) \\
%
&= -g^{\nu\phi} \xi^\alpha \xi_\alpha B_{\delta} B^{\beta} (\delta^\mu_\phi \delta^\delta_\beta - \delta^\mu_\beta \delta^\delta_\phi)
  + g^{\nu\gamma} \xi^\alpha \xi_\gamma B_{\delta} B^{\beta} (\delta^\mu_\alpha \delta^\delta_\beta - \delta^\mu_\beta \delta^\delta_\alpha)
  - g^{\nu\phi} \xi^\alpha \xi_\beta B_{\delta} B^{\beta} (\delta^\mu_\alpha  \delta^\delta_\phi - \delta^\mu_\phi \delta^\delta_\alpha) \\
%
&= -g^{\mu\nu} \xi^\alpha \xi_\alpha B_{\beta} B^{\beta}
  + g^{\nu\delta} \xi^\alpha \xi_\alpha B_{\delta} B^{\mu}
  + g^{\nu\gamma} \xi^\mu \xi_\gamma B_{\beta} B^{\beta}
  - g^{\nu\gamma} \xi^\alpha \xi_\gamma B_{\alpha} B^{\mu}
  - g^{\nu\delta} \xi^\mu \xi_\beta B_{\delta} B^{\beta}
  + g^{\mu\nu} \xi^\alpha \xi_\beta B_{\alpha} B^{\beta} \\
%
&= g^{\mu\nu} B_{\beta} B^{\beta} - B^{\mu} B^{\nu} + \xi^\mu \xi^\nu B_{\beta} B^{\beta} \\
\end{align}

which leaves us with

\begin{align}
F^{\mu\lambda}F^{\nu}_{\hphantom{\nu}\lambda} &= \xi^\mu \xi^\nu E^\lambda E_\lambda - E^\mu E^\nu +  \xi_\alpha E_\lambda B_{\beta} (\xi^\mu \epsilon^{\alpha\nu\lambda\beta} + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}) \\
&\ \ \ \ \ - \xi_\lambda \xi_\alpha B_{\beta} (\epsilon^{\alpha\nu\lambda\beta} E^\mu + \epsilon^{\alpha\mu\lambda\beta} E^\nu) + g^{\mu\nu} B_{\beta} B^{\beta} - B^{\mu} B^{\nu} + \xi^\mu \xi^\nu B_{\beta} B^{\beta} \\
%
&= \xi^\mu \xi^\nu E^\lambda E_\lambda + \left(g^{\mu\nu} + \xi^\mu \xi^\nu\right) B_{\beta} B^{\beta} - E^\mu E^\nu - B^{\mu} B^{\nu} +  \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta} + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}\right)
\end{align}

where we have again used the fact that $\xi_\alpha \xi_\sigma \epsilon^{\sigma\alpha\beta\mu} = 0$.

<a id='T_solve'></a>

## Step 1.3: Finishing Touches \[Back to [top](#toc)\]
$$\label{T_solve}$$

Now we just need to put all these pieces together. The final expression for the electromagnetic stress-energy tensor for an arbitrary observer is

\begin{align}
T_{EM}^{\mu\nu}|_{\xi\rightarrow n} &= \frac{1}{4\pi}\left( \xi^\mu \xi^\nu E^\lambda E_\lambda + \left(g^{\mu\nu} + \xi^\mu \xi^\nu\right) B_{\beta} B^{\beta} - E^\mu E^\nu - B^{\mu} B^{\nu} +  \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta} + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}\right) \\
  - \frac{1}{4}g^{\mu\nu}(-2E^\beta E_\beta + 2B^{\mu} B_{\mu}) \right) \\
&= \frac{1}{8\pi}\left(g^{\mu\nu} + 2\xi^\mu \xi^\nu\right)\left(E^\alpha E_\alpha + B^{\alpha} B_{\alpha}\right)
  - \frac{1}{4\pi}\left(E^\mu E^\nu + B^{\mu} B^{\nu}\right)
  + \frac{1}{4\pi} + \xi_\alpha E_\lambda B_{\beta} \left(\xi^\mu \epsilon^{\alpha\nu\lambda\beta} + \xi^\nu \epsilon^{\alpha\mu\lambda\beta}\right)
\end{align}

Now, we can consider two special cases. First, let's consider the normal observer, in which case the normalized 4-velocity becomes the unit normal vector $\xi\rightarrow n$

\begin{align}
n^\nu &= \left( \frac{1}{\alpha}, -\frac{\beta^i}{\alpha} \right) \\
n_\nu &= \left( -\alpha, 0 \right).
\end{align}

If we also define

$$
\xi_\alpha \epsilon^{\alpha\nu\lambda\beta} \equiv \epsilon^{\nu\lambda\beta}
$$

then the stress-energy tensor becomes

\begin{align}
T_{EM}^{\mu\nu}|_{\xi\rightarrow n} &= \frac{1}{8\pi}\left(g^{\mu\nu} + 2n^\mu n^\nu\right)\left(E^\alpha E_\alpha + B_{\alpha} B^{\alpha}\right)
  - \frac{1}{4\pi}\left(E^\mu E^\nu + B^{\mu} B^{\nu}\right)
  + \frac{1}{4\pi}E_\alpha B_{\beta} n^{(\mu} \epsilon^{\nu)\alpha\beta}
\end{align}

We can also consider the observer co-moving with the fluid. We can call this co-moving velocity $u_\mu$. In the case of ideal conditions (that is to say, a fluid with perfect conductivity), we have (from Ohm's law) that

$$
u_\mu F^{\mu\nu}=0
$$

which implies that

$$
E^\mu = u_\nu F^{\mu\nu} = 0
$$

Since $E^\mu = 0$, the stress-energy tensor simplifies significantly in this frame:

$$
T_{EM}^{\mu\nu}|_{\xi\rightarrow u} = \frac{1}{8\pi}\left(g^{\mu\nu} + 2u^\mu u^\nu\right)B_{\beta} B^{\beta}
  - \frac{1}{4\pi}B^{\mu} B^{\nu}
$$

However, variable $b$ is often used instead of $B^\mu_{(u)}$, given by

$$
b^\mu = \frac{B^\mu_{(u)}}{\sqrt{4\pi}}
$$

In terms of $b$, the electromagnetic energy-momentum tensor is given by

$$
T_{EM}^{\mu\nu}|_{\xi\rightarrow u} = b^2\left(\frac{1}{2}g^{\mu\nu} + u^\mu u^\nu\right)
  - b^{\mu} b^{\nu}
$$

<a id='suitable_equations'></a>

# Step 2: Evolution Equations for the Conservative Variables \[Back to [top](#toc)\]
$$\label{suitable_equations}$$

Before discussing the evolution variables, recall the 3+1 decomposition which is used for GRMHD. The metric is given by
$$/
ds^2 = -\alpha^2 dt^2 + \gamma_{i j}\left(dx^i + \beta^j dt\right)\left(dx^j + \beta^i dt\right)
$$
where $\alpha$, $\beta^i$, and $\gamma_{i j}$ are the lapse, shift, and spatial metric, respectively. These quantities are used in the definitions of the conservative variables which we will soon introduce.

In GRMHD simulations, there are two sets of variables which can be used for evolution. The first set, called primitive variables, are

$$
\mathbf{P} = \left\{\rho_b,P,v^i,B^i\right\}
$$

where $\rho_b$ is the fluid rest-mass density, $P$ is the pressure, $v^i$ is the rescaled fluid three-velocity $u^i/u^0$, and $B^i$ is the magnetic field. An important sidenote here is that different GRMHD codes define different $v^i$. As an example, the widely used GRHydro thorn in the Einstein Toolkit uses the Valencia formalism

$$
v^i = \frac{u^i}{W} + \frac{\beta^i}{\alpha}
$$

where $W=\sqrt{1-v^i v_i}$ is the Lorentz factor. The quantity $v^i$ used by IllinoisGRMHD is **not** the Valencia three-velocity, so those interested in seeing more about the Valencia formulation of the evolution equations can look to the [GRHydro paper](https://arxiv.org/abs/1304.5544) and its references.

The second set, the conservative variables, are

$$
\mathbf{C} = \left\{\tilde{\rho},\tilde{\tau},\tilde{S}_i,\tilde{B}^i\right\}
$$

where

\begin{align}
\tilde{\rho} &= \sqrt{-g}\rho_b u^0 = \alpha\sqrt{\gamma}\rho_b u^0 \\
\tilde{\tau} &= \sqrt{\gamma} n_\mu n_\nu T^{\mu\nu} - \tilde{\rho} = \alpha^2 \sqrt{\gamma} T^{00} - \tilde{\rho} \\
\tilde{S}_i &= \sqrt{\gamma}S_i = \alpha \sqrt{\gamma}T^0_i \\
\tilde{B}^i &= \sqrt{\gamma}B^i
\end{align}

One can design numerical codes to evolve either set of variables, and we choose to evolve $\mathbf{C}$. Conservative variables are used to evolve these systems because the evolution of primitive variables leads to highly unstable simulations. This change of variables for the sake of stability has some similarity to the change from ADM variables to BSSN variables when evolving space-time, where one set of variables has substantially better numerical stability.

Another question that may arise is why we do not simply evolve the energy variable $\alpha^2 \sqrt{\gamma} T^{00}$ and instead subtract $\tilde{\rho}$. It turns out that evolving $T^{00}$ directly is inaccurate because the magnetic and internal energy density can be orders of magnitude smaller than the rest mass density (this is mentioned in e.g. the [`HARM` software paper](https://arxiv.org/pdf/astro-ph/0301509.pdf)). Because of this, we subtract the baryonic density evolution equation from the energy evolution equation to find a more accurate and stable evolution equation.

All of our previous derivations, of course, only considered the electromagnetic contribution. We do not derive the hydrodynamic part of the stress-energy tensor and instead merely state that in the ideal MHD limit the full stress-energy tensor takes the form

$$
T^{\mu\nu} = \rho_b h u^\mu u^\nu + P g^{\mu\nu} + T_{EM}^{\mu\nu}
$$

Our next step is to derive the actual evolution equations for our designated evolution variables $\mathbf{C}$. To this end, we will use the identity

$$
\Gamma^\nu_{\gamma\nu} = \frac{1}{\sqrt{-g}} \partial_\gamma\sqrt{-g}
$$

to simplify the equations. We will be frequently using this to condense two terms into one in the following way:

\begin{align}
\partial_\nu T_\mu^\nu + \Gamma^\nu_{\gamma\nu} T^\gamma_\mu &= \partial_\nu T_\mu^\nu + \frac{1}{\sqrt{-g}} T^\gamma_\mu \partial_\gamma\sqrt{-g} \\
%
&= \partial_\nu T_\mu^\nu + \frac{1}{\sqrt{-g}} T^\gamma_\mu \partial_\gamma\sqrt{-g} \\
%
&= \frac{1}{\sqrt{-g}} \left(\sqrt{-g} \partial_\nu T_\mu^\nu + T^\gamma_\mu \partial_\gamma\sqrt{-g}\right) \\
%
&= \frac{1}{\sqrt{-g}} \partial_\gamma \left(\sqrt{-g}T^\gamma_\mu\right)
\end{align}

<a id='density_eqn'></a>

## Step 2.1: Evolution Equation for the Fluid Density $\tilde{\rho}$ \[Back to [top](#toc)\]
$$\label{density_eqn}$$

The first evolution equation comes from the baryon number conservation equation

\begin{align}
0 &= \nabla_\nu\left( \rho_b u^\nu \right) \\
&= \partial_\nu\left( \rho_b u^\nu \right) + \Gamma^\nu_{\gamma\nu} \rho_b u^\gamma \\
&= \partial_\nu\left( \rho_b u^\nu \right) + \rho_b u^\gamma \frac{1}{\sqrt{-g}} \partial_\gamma\sqrt{-g} \\
&= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g}\rho_b u^\nu \right) \\
&= \partial_t \left(\alpha\sqrt{\gamma}\rho_b u^0 \right) + \partial_i \left(\alpha\sqrt{\gamma}\rho_b u^i \right) \\
&= \partial_t \tilde{\rho} + \partial_i \left(\frac{\tilde{\rho}}{u^0} u^i \right) \\
&= \partial_t \tilde{\rho} + \partial_i \left(\tilde{\rho} v^i \right)
\end{align}

<a id='momentum_eqn'></a>

## Step 2.2: Evolution Equation for Momentum Density Variable $\tilde{S}_i$ \[Back to [top](#toc)\]
$$\label{momentum_eqn}$$

The next two sections derive their equations from the energy-momentum conservation equation

\begin{align}
\nabla_\nu T_\mu^\nu &= 0 \\
&= \partial_\nu T_\mu^\nu + \Gamma^\nu_{\gamma\nu} T^\gamma_\mu - \Gamma^\gamma_{\mu\nu} T^\nu_\gamma
\end{align}

Rearranging and using the identity for $\Gamma^\nu_{\gamma\nu}$, we get

$$
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T_\mu^\nu \right) = \Gamma^\gamma_{\mu\nu} T^\nu_\gamma
$$

For the momentum density variable $\tilde{S}_i$, we consider the indices $\mu\in\{1,2,3\}$. Replacing $\mu$ with $i$ to find the evolution equation,

\begin{align}
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T_i^\nu \right) &= \Gamma^\gamma_{i\nu} T^\nu_\gamma \\
%
\frac{1}{\sqrt{-g}}\partial_t \left(\sqrt{-g} T_i^0 \right) + \frac{1}{\sqrt{-g}}\partial_j \left(\sqrt{-g} T_i^j \right) &= \frac{1}{2}T^\nu_\gamma g^{\gamma\beta}\left( g_{\beta\nu,i} + g_{i\beta,\nu} - g_{i\nu,\beta} \right) \\
%
\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right) &= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu,i} + \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}\left( g_{i\beta,\nu} - g_{i\nu,\beta} \right) \\
%
\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right) &= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu,i}
\end{align}

<a id='energy_eqn'></a>

## Step 2.3: Evolution Equation for Energy Density Variable $\tilde{\tau}$ \[Back to [top](#toc)\]
$$\label{energy_eqn}$$

The energy evolution equation comes from the same conservation equation as the previous section, but with $\mu=0$. However, we choose to start with $\mu$ raised instead of lowered. Then,

\begin{align}
0 &= \nabla_\nu T^{\mu\nu} \\
&= \partial_\nu T^{\mu\nu} + \Gamma^\mu_{\sigma\nu} T^{\sigma\nu} + \Gamma^\nu_{\sigma\nu} T^{\mu\sigma} \\
&= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T^{\mu\nu} \right) + \Gamma^\mu_{\sigma\nu} T^{\sigma\nu}
\end{align}

Setting $\mu=0$,

\begin{align}
\frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} T^{0\nu} \right) &= -\Gamma^0_{\sigma\nu} T^{\sigma\nu} \\
\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right) + \partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right) &= -\alpha\sqrt{\gamma}\Gamma^0_{\sigma\nu} T^{\sigma\nu} \\
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{\sigma\nu} g^{0\beta}\left( g_{\beta\sigma,\nu} + g_{\beta\nu,\sigma} - g_{\sigma\nu,\beta} \right) \\
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{\sigma\nu} g^{0\beta}\left(2 g_{\beta\sigma,\nu} - g_{\sigma\nu,\beta} \right)
\end{align}

where the final step follows from the fact that $T^{\sigma\nu}$ is symmetric. Now, we will consider the right-hand side for various components of $T^{\sigma\nu}$. For these derivations, we will use derivations relating to the extrinsic curvature in [Appendix B](#curv_tensor). For $T^{00}$,

\begin{align}
-\frac{\alpha\sqrt{\gamma}}{2} T^{00} g^{0\beta}\left(2 g_{\beta 0,0} - g_{00,\beta} \right) &= -\frac{\alpha\sqrt{\gamma}}{2} T^{00} \left[2 \left(g^{00} g_{00,0} + g^{0i}g_{0i,0}\right) - g^{00}g_{00,0} - g^{0i}g_{00,i} \right] \\
%
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{00} \left( g^{00} g_{00,0} + 2g^{0i}g_{0i,0}  - g^{0i}g_{00,i} \right) \\
%
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{00} \left[ -\alpha^{-2} \partial_t \left(\beta^2 - \alpha^2\right) + \frac{2\beta^i}{\alpha^2} \partial_t \beta_i - \frac{\beta^i}{\alpha^2} \partial_i \left(\beta^2 - \alpha^2\right) \right] \\
%
&= -\frac{\sqrt{\gamma}}{2\alpha} T^{00} \left[ -\partial_t \left(\beta^2 - \alpha^2\right) + 2\beta^i \partial_t \beta_i - \beta^i \partial_i \left(\beta^2 - \alpha^2\right) \right] \\
%
&= -\frac{\sqrt{\gamma}}{\alpha} T^{00} \left( \alpha\partial_t \alpha + \beta^i \alpha\partial_i \alpha - \beta^i \beta^j \partial_i \beta_j \right) \\
%
&= \sqrt{\gamma} T^{00} \left( \beta^i \beta^j K_{ij} - \partial_t \alpha - \beta^i \partial_i \alpha \right)
\end{align}

where we have used the identity in [Appendix B](#curv_tensor). Next, we look at the mixed term $T^{0i} + T^{i0}$:

\begin{align}
-\frac{\alpha\sqrt{\gamma}}{2} T^{0i} g^{0\beta}\left(2 g_{\beta 0,i} - g_{0i,\beta} + 2 g_{\beta i,0} - g_{i 0,\beta} \right) &= -\alpha\sqrt{\gamma} T^{0i} \left[ g^{00}\left( g_{00,i} + g_{0i,0} \right) + g^{0j}\left( g_{0j,i} + g_{ij,0} \right) - \left(g^{00} g_{0i,0} + g^{0j} g_{0i,j}\right) \right] \\
%
&= -\alpha\sqrt{\gamma} T^{0i} \left[ g^{00} g_{00,i} + g^{0j}\left( g_{0j,i} + g_{ij,0} - g_{0i,j} \right) \right] \\
%
&= -\alpha\sqrt{\gamma} T^{0i} \left[ -\alpha^{-2} \partial_i \left( \beta^2 - \alpha^2 \right) + \frac{\beta^j}{\alpha^2} \left( \partial_i \beta_j + \partial_t \gamma_{ij} - \partial_j \beta_i \right) \right] \\
%
&= -\frac{\sqrt{\gamma}}{\alpha} T^{0i} \left[ 2\alpha\partial_i \alpha + \beta^j \left(  + \partial_t \gamma_{ij} - \partial_i \beta_j - \partial_j \beta_i \right) \right] \\
%
&= -\frac{\sqrt{\gamma}}{\alpha} T^{0i} \left[ 2\alpha\partial_i \alpha + \beta^j \left( - 2\alpha K_{ij} - 2\Gamma^k_{ij}\beta_k \right) \right] \\
%
&= 2\sqrt{\gamma} T^{0i} \left( \beta^j K_{ij} - \partial_i \alpha \right) + \frac{2\sqrt{\gamma}}{\alpha} T^{0i} \beta^j \Gamma^k_{ij}\beta_k \\
%
&= 2\sqrt{\gamma} T^{0i} \left( \beta^j K_{ij} - \partial_i \alpha \right)
\end{align}

Finally, for $T^{ij}$ we have

\begin{align}
-\frac{\alpha\sqrt{\gamma}}{2} T^{ij} g^{0\beta}\left(2 g_{\beta i,j} - g_{ij,\beta} \right) &= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[ 2\left(g^{00} g_{0i,j} + g^{0k} g_{ki,j}\right) - g^{0\beta} g_{ij,\beta} \right] \\
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[ 2\left(g^{00} \partial_j \beta_i + g^{0k} \partial_j \gamma_{ki}\right) - g^{0\beta} \partial_\beta \gamma_{ij} \right] \\
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[ g^{00}\left( 2\partial_j \beta_i - \partial_t \gamma_{ij} \right) + g^{0k}\left( 2\partial_j \gamma_{ki} - \partial_k \gamma_{ij} \right) \right] \\
&= -\frac{\alpha\sqrt{\gamma}}{2} T^{ij} \left[ -\alpha^{-2}\left( 2\partial_j \beta_i - \partial_t \gamma_{ij} \right) + \frac{\beta^k}{\alpha^2}\left( 2\partial_j \gamma_{ki} - \partial_k \gamma_{ij} \right) \right] \\
\end{align}

We can see that, thanks to the symmetry of $T^{ij}$,

\begin{align}
T^{ij}\beta^k\left(2\partial_j \gamma_{ki} - \partial_k \gamma_{ij}\right) &= T^{ij}\beta_n g^{kn}\left(\partial_j \gamma_{ki} + \partial_i \gamma_{kj} - \partial_k \gamma_{ij}\right) \\
&= 2T^{ij}\Gamma^k_{ij}\beta_k
\end{align}

Using this and using the definition of $K_{ij}$ to replace $\partial_t \gamma_{ij}$,

\begin{align}
-\frac{\alpha\sqrt{\gamma}}{2} T^{ij} g^{0\beta}\left(2 g_{\beta i,j} - g_{ij,\beta} \right) &= \frac{\sqrt{\gamma}}{2\alpha} T^{ij} \left( 2\partial_j \beta_i + 2\alpha K_{ij} - \partial_i\beta_j - \partial_j\beta_i + 2\Gamma^k_{ij}\beta_k - 2\Gamma^k_{ij}\beta_k \right) \\
&= \frac{\sqrt{\gamma}}{2\alpha} T^{ij} \left( 2\alpha K_{ij} + \partial_j \beta_i - \partial_i\beta_j \right) \\
&= \sqrt{\gamma} T^{ij} K_{ij} \\
\end{align}

where the final step again follows from the symmetry of $T^{ij}$. Putting all this together and multiplying through by $\alpha$,

\begin{align}
\alpha\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right) + \alpha\partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right) &= \alpha\sqrt{\gamma} T^{00} \left( \beta^i \beta^j K_{ij} - \partial_t \alpha - \beta^i \partial_i \alpha \right) \\
&\ \ \ \ \ \ \ + 2\alpha\sqrt{\gamma} T^{0i} \left( \beta^j K_{ij} - \partial_i \alpha \right) + \alpha\sqrt{\gamma} T^{ij} K_{ij} \\
%
\alpha\sqrt{\gamma} T^{00}\partial_t \alpha + 2\alpha\sqrt{\gamma} T^{0i}\partial_i \alpha + \alpha\partial_t \left(\alpha\sqrt{\gamma} T^{00} \right) + \alpha\partial_i \left(\alpha\sqrt{\gamma} T^{0i} \right) &= \alpha\sqrt{\gamma} T^{00} \left( \beta^i \beta^j K_{ij} - \beta^i \partial_i \alpha \right) + 2\alpha\sqrt{\gamma} T^{0i}\beta^j K_{ij} \\
&\ \ \ \ \ \ \ - \alpha\sqrt{\gamma} T^{0i} \partial_i \alpha + \alpha\sqrt{\gamma} T^{ij} K_{ij} \\
%
\partial_t \left(\alpha^2\sqrt{\gamma} T^{00} \right) + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} \right) &= \alpha\sqrt{\gamma} \left[ \left( T^{00} \beta^i \beta^j + 2T^{0i}\beta^j + T^{ij} \right) K_{ij} - \left( T^{00} \beta^i + T^{0i} \right)\partial_i \alpha \right]
\end{align}

We have finally arrived at the evolution equation. However, we still need to get it in terms of the evolution variable $\tilde{\tau} = \alpha^2 \sqrt{\gamma} T^{00} - \tilde{\rho}$. To do so, we can add the fluid density evolution equation to the left-hand side. We also define the source term

$$
s = \alpha\sqrt{\gamma} \left[ \left( T^{00} \beta^i \beta^j + 2T^{0i}\beta^j + T^{ij} \right) K_{ij} - \left( T^{00} \beta^i + T^{0i} \right)\partial_i \alpha \right]
$$

Then, the energy evolution equation becomes

\begin{align}
\partial_t \left(\alpha^2\sqrt{\gamma} T^{00} \right) + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} \right) - \partial_t \tilde{\rho} - \partial_i \left(\tilde{\rho} v^i \right) &= s \\
\partial_t \tilde{\tau} + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} - \tilde{\rho} v^i \right) &= s
\end{align}

<a id='induction_eqn'></a>

## Step 2.4: Evolution Equation for the Magnetic Field $\tilde{B}^i$ \[Back to [top](#toc)\]
$$\label{induction_eqn}$$


To find a convenient start for the derivation, I find the dual of Maxwell's equation. Recall that the dual of the Faraday tensor is

$$
F^{*\mu\nu} = \frac{1}{2} \epsilon^{\mu\nu\alpha\beta}F_{\alpha\beta}
$$

Also, since the covariant derivative of the metric is 0,

$$
\epsilon^{\alpha\lambda\mu\nu} F_{\mu\nu;\lambda} = \nabla_\lambda \left(\epsilon^{\alpha\lambda\mu\nu} F_{\mu\nu}\right) = 2\nabla_\lambda F^{*\alpha\lambda}
$$

Therefore, Maxwell's equation becomes

\begin{align}
0 &= \epsilon^{\alpha\lambda\mu\nu} F_{[\mu\nu;\lambda]} \\
&= \epsilon^{\alpha\lambda\mu\nu} \left( F_{\mu\nu;\lambda} - F_{\lambda\nu;\mu} + F_{\nu\lambda;\mu} - F_{\mu\lambda;\nu} + F_{\lambda\mu;\nu} - F_{\nu\mu;\lambda} \right) \\
%
&= \epsilon^{\alpha\lambda\mu\nu}F_{\mu\nu;\lambda} + \epsilon^{\alpha\lambda\mu\nu} \nabla_\mu \left( F_{\nu\lambda} - F_{\lambda\nu} \right) + \epsilon^{\alpha\lambda\mu\nu} \nabla_\nu \left( F_{\lambda\mu} - F_{\mu\lambda}\right) - \epsilon^{\alpha\lambda\mu\nu}\nabla_\lambda F_{\nu\mu} \\
%
&= 2\nabla_\lambda F^{*\alpha\lambda} + 2\epsilon^{\alpha\lambda\mu\nu} \nabla_\mu F_{\nu\lambda} + 2\epsilon^{\alpha\lambda\mu\nu} \nabla_\nu F_{\lambda\mu} - \epsilon^{\alpha\lambda\mu\nu}\nabla_\lambda F_{\nu\mu}
\end{align}

By exploiting the symmetries of the Levi-Civita tensor,

\begin{align}
0 &= 2\nabla_\lambda F^{*\alpha\lambda} + 2\epsilon^{\alpha\mu\nu\lambda} \nabla_\mu F_{\nu\lambda} + 2\epsilon^{\alpha\nu\lambda\mu} \nabla_\nu F_{\lambda\mu} + \epsilon^{\alpha\lambda\nu\mu}\nabla_\lambda F_{\nu\mu} \\
%
&= 2\nabla_\lambda F^{*\alpha\lambda} + 4\nabla_\mu F^{*\alpha\mu} + 4\nabla_\nu F^{*\alpha\nu} + 2\nabla_\lambda F^{*\alpha\lambda} \\
%
&= 12\nabla_\lambda F^{*\alpha\lambda} \\
0 &= \nabla_\lambda F^{*\alpha\lambda} = \partial_\lambda F^{*\alpha\lambda} + \Gamma^\alpha_{\beta\lambda} F^{*\beta\lambda} + \Gamma^\lambda_{\beta\lambda} F^{*\alpha\beta} \\
0 &= \frac{1}{\sqrt{-g}}\partial_\lambda \left(\sqrt{-g} F^{*\alpha\lambda}\right)
\end{align}

where the first $\Gamma$ term dissapears because it is summing over the multiplication of symmetric and anti-symmetric objects.

We can now find the magnetic equations from this. Taking the time component simply gives us the no-monopole constraint. To see this, we need to first consider the individual components of the dual. By the definition of $B$,

$$
B^\mu = \frac{1}{2} \epsilon^{\mu\nu\beta\alpha}n_\nu F_{\alpha\beta} = n_\nu F^{*\nu\mu}
$$

Since $B$ for the normal observer is purely spatial ($B^\mu n_\mu=0$), this implies that $F^{*00}=0$. Then, the time component of the magnetic equations is

\begin{align}
0 &= \frac{1}{\sqrt{-g}}\partial_\lambda \left(\sqrt{-g} F^{*0\lambda}\right) \\
&= \partial_i \left(\alpha\sqrt{\gamma} F^{*0i}\right) \\
&= \partial_i \left(\alpha\sqrt{\gamma} \frac{B^i}{\alpha}\right) \\
0 &= \partial_i \left(\tilde{B}^i \right)
\end{align}


Before examining the spatial components of Maxwell's equations, we will do some preliminary work to make our lives easier later. In the co-moving frame,

\begin{align}
F^{\mu\nu} &= u^\mu E^\nu - u^\nu E^\mu + u_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta} \\
&= u_\gamma \epsilon^{\gamma\mu\nu\delta} B_{\delta}
\end{align}

Then, the dual is

\begin{align}
F^{*\mu\nu} &= \frac{1}{2} \epsilon^{\mu\nu\alpha\beta} F_{\alpha\beta} \\
&= \frac{1}{2} \epsilon^{\mu\nu\alpha\beta} u^\gamma \epsilon_{\gamma\alpha\beta\delta} B^{\delta} \\
&= \frac{1}{2} \epsilon^{\alpha\beta\mu\nu}\epsilon_{\alpha\beta\gamma\delta} u^\gamma  B^{\delta} \\
&= \left( \delta^\mu_\delta \delta^\nu_\gamma - \delta^\mu_\gamma \delta^\nu_\delta \right) u^\gamma  B^{\delta} \\
&= u^\nu  B^{\mu} - u^\mu  B^{\nu}
\end{align}

where we have again used the Levi-Civita identity from [Appendix A](#levi_civita). Next, we need to find a relationship between the magnetic field of the normal observer and the co-moving observer. For clarity, let the co-moving magnetic field be $B^\mu_{(u)}$ and the normal magnetic field be $B^\mu$. Then, we can define a projection operator

$$
P^{\mu\nu} = g_{\mu\nu} + u_\mu u_\nu
$$

Naturally, the projection of the co-moving magnetic field should simply project back into the same field:

\begin{align}
P^\mu_\nu B^\nu_{(u)} &= \left( \delta^\mu_\nu + u^\mu u_\nu \right)B^\nu_{(u)} \\
&= B^\mu_{(u)}
\end{align}

where we have used the orthogonality relation $u_\nu B^\nu_{(u)}=0$. Projecting the normal observer's magnetic field,

\begin{align}
P^\mu_\nu B^\nu &= P^\mu_\nu n_\alpha F^{*\alpha\nu} \\
&= P^\mu_\nu n_\alpha \left( u^\nu  B^{\alpha}_{(u)} - u^\alpha  B^{\nu}_{(u)} \right) \\
&= n_\alpha \left( \delta^\mu_\nu + u^\mu u_\nu \right) \left( u^\nu  B^{\alpha}_{(u)} - u^\alpha  B^{\nu}_{(u)} \right) \\
%
&= n_\alpha \left( u^\mu  B^{\alpha}_{(u)} +  u^\mu u_\nu u^\nu  B^{\alpha}_{(u)} - u^\alpha  B^{\mu}_{(u)} - u^\mu u_\nu u^\alpha  B^{\nu}_{(u)} \right) \\
%
&= n_\alpha u^\mu B^{\alpha}_{(u)} \left( 1 + u_\nu u^\nu \right) - n_\alpha u^\alpha \left( B^{\mu}_{(u)} + u^\mu u_\nu B^{\nu}_{(u)} \right) \\
%
&= - \alpha u^0 B^{\mu}_{(u)} \\
\Rightarrow B^{\mu}_{(u)} &= \frac{B^\mu + u^\mu u_\nu B^\nu}{\alpha u^0}
\end{align}

Since we will only need the spatial component for our purposes,

\begin{align}
B^{i}_{(u)} &= \frac{B^i + u^i u_\nu B^\nu}{\alpha u^0} \\
&= \frac{B^i + u^i u_j B^j}{\alpha u^0} \\
&= \frac{B^i}{\alpha u^0} + \frac{v^i u_j B^j}{\alpha}
\end{align}

where $v^i$ is the conservative variable $u^i/u^0$. Finally, the spatial components of the dual of Maxwell's equations gives the induction equation for the magnetic field:

\begin{align}
0 &= \frac{1}{\sqrt{-g}}\partial_\nu \left(\sqrt{-g} F^{*i \nu}\right) \\
0 &= \partial_t \left(\alpha\sqrt{\gamma} F^{*i 0}\right) + \partial_j \left(\alpha\sqrt{\gamma} F^{*i j}\right) \\
%
0 &= \partial_t \left(\sqrt{\gamma} B^i \right) + \partial_j \left(\alpha\sqrt{\gamma} \left[ u^j B^i_{(u)} - u^i B^j_{(u)} \right] \right) \\
%
0 &= \partial_t \tilde{B}^i + \partial_j \left(\alpha\sqrt{\gamma} \left[ u^j \left( \frac{B^i}{\alpha u^0} + \frac{v^i u_k B^k}{\alpha} \right) - u^i \left( \frac{B^j}{\alpha u^0} + \frac{v^j u_k B^k}{\alpha} \right) \right] \right) \\
%
0 &= \partial_t \tilde{B}^i  + \partial_j \left(v^j \tilde{B}^i - v^i \tilde{B}^j + \frac{\sqrt{\gamma}}{u^0}\left[ u^j u^i u_k B^k - u^i u^j u_k B^k \right] \right) \\
%
0 &= \partial_t \tilde{B}^i  + \partial_j \left(v^j \tilde{B}^i - v^i \tilde{B}^j \right)
\end{align}



<a id='summary'></a>

## Step 2.5: Summary of the Conservative Variable Evolution Equations \[Back to [top](#toc)\]
$$\label{summary}$$

In the previous sections, we have derived the evolution equations for the conservative variables \mathbf{C}. To summarize, these are

\begin{align}
\partial_t \tilde{\rho} + \partial_i \left(\tilde{\rho} v^i \right) &= 0 \\
\partial_t \tilde{\tau} + \partial_i \left(\alpha^2\sqrt{\gamma} T^{0i} - \tilde{\rho} v^i \right) &= s \\
\partial_t \tilde{S}_i + \partial_j \left(\alpha\sqrt{\gamma} T_i^j \right) &= \frac{\alpha\sqrt{\gamma}}{2} T^{\beta\nu}g_{\beta\nu,i} \\
\partial_t \tilde{B}^i  + \partial_j \left(v^j \tilde{B}^i - v^i \tilde{B}^j \right) &= 0
\end{align}

where

$$
s = \alpha\sqrt{\gamma} \left[ \left( T^{00} \beta^i \beta^j + 2T^{0i}\beta^j + T^{ij} \right) K_{ij} - \left( T^{00} \beta^i + T^{0i} \right)\partial_i \alpha \right]
$$

<a id='levi_civita'></a>

# Appendix A: Levi-Civita Contractions \[Back to [top](#toc)\]
$$\label{levi_civita}$$

In order to simplify the various combinations of $F^{\mu\nu}$, contractions of the Levi-Civita tensor are required. First, summing over the first indices of two such tensors yields

\begin{align}
\epsilon^{\alpha\beta\sigma\mu} \epsilon_{\alpha\beta\gamma\delta} &= -\delta^{\beta\sigma\mu}_{\hphantom{\beta\sigma\nu}\beta\gamma\delta} \\
%
&= -\delta^\beta_\beta \delta^{\sigma\mu}_{\gamma\delta} 
  + \delta^\beta_\gamma \delta^{\sigma\mu}_{\beta\delta}
  - \delta^\beta_\delta \delta^{\sigma\mu}_{\beta\gamma} \\
%
&= -4 \left(\delta^\sigma_\gamma \delta^\mu_\delta - \delta^\sigma_\delta \delta^\mu_\gamma\right)
  + \delta^\beta_\gamma \left(\delta^\sigma_\beta \delta^\mu_\delta - \delta^\sigma_\delta \delta^\mu_\beta\right)
  - \delta^\beta_\delta \left(\delta^\sigma_\beta \delta^\mu_\gamma - \delta^\sigma_\gamma \delta^\mu_\beta\right) \\
%
&= -4 \left(\delta^\sigma_\gamma \delta^\mu_\delta - \delta^\sigma_\delta \delta^\mu_\gamma\right)
  + 2\left(\delta^\sigma_\gamma \delta^\mu_\delta - \delta^\sigma_\delta \delta^\mu_\gamma\right) \\
%
&= 2\left(\delta^\sigma_\delta \delta^\mu_\gamma - \delta^\sigma_\gamma \delta^\mu_\delta\right)
\end{align}

Second, summing over the third index (with all other components raised) yields

\begin{align}
\epsilon^{\gamma\mu\lambda\delta} \epsilon^{\alpha\nu\hphantom{\lambda}\beta}_{\hphantom{\alpha\nu}\lambda} &= g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta} \epsilon^{\gamma\mu\lambda\delta} \epsilon_{\tau\phi\lambda\theta} \\
%
&= -g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta} \delta^{\gamma\mu\delta}_{\tau\phi\theta} \\
%
&= -g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta}\left( \delta^\gamma_\tau\delta^{\mu\delta}_{\phi\theta}
  - \delta^\gamma_\phi \delta^{\mu\delta}_{\tau\theta}
  + \delta^\gamma_\theta \delta^{\mu\delta}_{\tau\phi} \right) \\
%
&= -g^{\alpha\tau}g^{\nu\phi}g^{\beta\theta}
  \left( \delta^\gamma_\tau (\delta^\mu_\phi \delta^\delta_\theta - \delta^\mu_\theta \delta^\delta_\phi)
  - \delta^\gamma_\phi (\delta^\mu_\tau \delta^\delta_\theta - \delta^\mu_\theta \delta^\delta_\tau)
  + \delta^\gamma_\theta (\delta^\mu_\tau  \delta^\delta_\phi - \delta^\mu_\phi \delta^\delta_\tau) \right) \\
\end{align}

<a id='curv_tensor'></a>

# Appendix B: Extrinsic Curvature Identities \[Back to [top](#toc)\]
$$\label{curv_tensor}$$

The source term of the energy equation involves the extrinsic curvature. Defined in terms of the 3+1 formalism quantities, the extrinsic curvature is
\begin{align}
K_{ij} &= -\frac{1}{2\alpha}\left( \partial_t \gamma_{ij} - \nabla_i\beta_j - \nabla_j\beta_i \right) \\
&= -\frac{1}{2\alpha}\left( \partial_t \gamma_{ij} - \partial_i\beta_j - \partial_j\beta_i + 2\Gamma^k_{ij}\beta_k \right)
\end{align}

As for how the extrinsic curvature enters into the energy equations, it involves several relations between it and the other quantities which we must show. First, consider the contraction of the Christoffel symbol with $\beta^i \beta^j$:

\begin{align}
\beta^j \Gamma^k_{ij}\beta_k &= \beta^j \beta^k \left( \gamma_{k i,j} + \gamma_{k j,i} - \gamma_{i j,k} \right) \\
&= \beta^j \beta^k \gamma_{k j,i} \\
&= \alpha^2 \left( \gamma^{jk} - g^{jk} \right) \gamma_{k j,i} \\
&\propto \gamma^{jk}\partial_i \gamma_{k j} - g_{jk} \partial_i \gamma^{k j} \\
&\propto \gamma^{jk}\partial_i \gamma_{k j} - \gamma_{jk} \partial_i \gamma^{k j} \\
&\propto \gamma^{jk} \partial_i \gamma_{k j} - \gamma^{jk} \partial_i \gamma_{k j} \\
&= 0
\end{align}

The change $g_{jk}\rightarrow\gamma_{jk}$ is allowed because the 3-metric is identical to the 4-metric with only spatial indices. We can use this to find a relationship between the partial derivative of $\beta_j$ and the extrinsic curvature:

\begin{align}
\beta^i \beta^j \partial_i \beta_j &= \beta^i \beta^j \left( 2\alpha K_{ij} + \partial_t \gamma_{ij} - \partial_j\beta_i + 2\Gamma^k_{ij}\beta_k \right) \\
&= \beta^i \beta^j \left( 2\alpha K_{ij} + \partial_t \gamma_{ij} - \partial_i\beta_j + 2\Gamma^k_{ij}\beta_k \right) \\
&= \beta^i \beta^j \left( 2\alpha K_{ij} + \partial_t \gamma_{ij} - \partial_i\beta_j \right) \\
\Rightarrow \beta^i \beta^j \partial_i \beta_j &= \beta^i \beta^j \left( \alpha K_{ij} + \frac{1}{2}\partial_t \gamma_{ij} \right) \\
\Rightarrow \beta^i \beta^j \partial_i \beta_j &= \beta^i \beta^j \alpha K_{ij}
\end{align}

where in the second step we use the symmetry of $\beta^i \beta^j$. The final step uses the same trick as with $\beta^j \Gamma^k_{ij}\beta_k$, just with the derivative being $\partial_t$ instead of $\partial_i$. 

<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Derivation_of_GRMHD_Evolution_Equations")
```

    Created Tutorial-Derivation_of_GRMHD_Evolution_Equations.tex, and compiled
        LaTeX file to PDF file Tutorial-
        Derivation_of_GRMHD_Evolution_Equations.pdf

