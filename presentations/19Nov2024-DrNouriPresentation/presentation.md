# OVER\nVIEW
3+1 decomposition is an alternative view on the equations of general relativity. This view is the basics of numerical relativity and hence it's the most important basis of understanding the algorithms, methods, and simulations that we do.
---
# INTRO-\nDUCTION
The intrinsically "covariant view" of Einstein's theory of general relativity is based on the concept that all coordinates are equivalent and, hence, the distinction between spatial and time coordinates is more an organizational matter than a strict requirement of the theory. Yet, our experience, and the laws of physics on sufficiently large scales, do suggest that a distinction between the time coordinate and from the spatial ones is the most natural one in describing physical processes.
---
# DECOM-\nPOSITION
Given such constant-time hypersurface, $\Sigma_t$, belonging to the foliation $\Sigma$, we can introduce a timelike four-vector $n$ normal to the hypersurface at each event in spacetime and such that its dual one-form $\Omega := \nabla_t$ is parallel to the gradient of the coordinate $t$.

$$
n_\mu = A \Omega_\mu = \nabla_\mu t
$$

with $n_\mu = \{A,0,0,0\}$ and $A$ a constant to be determined. 
---
If we now require that the four-vector $n$ defines an observer and thus that it measures the corresponding four-velocity, then from the normalization condition on timelike four-vectors, $n^\mu n_\mu = -1$.
$$
\begin{align}
n^\mu n_\mu &= g^{\mu\nu} n_\mu n_\nu \\
&= g^{00} A^2 \\
&= -\frac1{\alpha^2}A^2,
\end{align}
$$
where we defined $g^{00} = -\frac{1}{\alpha^2}$. Now that:
$$
\begin{align}
-\frac{1}{\alpha^2}A^2 &= -1\\
A &= \pm\alpha\\
\end{align}
$$
$$
\left\{
\begin{matrix}+ & \text{ Past Direction} \\
- & \text{ Future Direction}
\end{matrix}
\right.
$$
and thus:
$$
A = -\alpha
$$
to showcase the future. 
---
# SPATIAL\nMETRIC
Using this normal vector we can define the spatial metric for the hypersurface:
$$
\begin{matrix}
\gamma_{\mu\nu} = g_{\mu\nu} + n_\mu n_\nu , & \gamma^{\mu\nu} = g^{\mu\nu} + n^\mu n^\nu
\end{matrix}
$$
here $\gamma^{0\mu} = 0, \gamma_{ij} = g_{ij}$ but in general  $\gamma^{ij} \not = g^{ij}$. It's important to note that:
$$
\gamma^{ij}\gamma_{jk} = \delta^i_k
$$
Which means that these are the inverse of each other, so that the spacial metric can be used for raising and lowering indices of purely spacial vectors and tensors.
---
# DECOMPOSING\nTOOLKIT
The spatial metric and the hypersurface vector provide us a good toolkit to decompose any tensors into a purely spatial part and a purely timelike part.
---
# SPATIAL\nPART
To obtain the spatial part we contract with the spatial projection operator:
$$
\gamma^{\mu}_{\cdot\nu} := g^{\mu\alpha}\gamma_{\alpha\nu}=g^\mu_{\cdot\nu} + n^\mu n_\nu = \delta^{\mu}_\nu + n^\mu n_\nu
$$
---
# TIMELIKE\nPART
To obtain the timelike part we contract with the time projection operator:
$$
N^\mu_{\cdot \nu} := -n^\mu n_\nu
$$
---
# GENERAL\nCASE
Therefore, a general four-vector can be decomposed as:
$$
U^\mu = \gamma^\mu_{\cdot \nu}U^\nu + N^{\mu}_{\cdot \nu}U^\nu
$$
Note that the spatial part is still a four-vector but now $V^t = 0$ whereas it has the covariant time component, $V_t = g_{\mu t} V^{\mu}$, which is non-zero in general. Analogous considerations can be done about tensors of any rank.
---
