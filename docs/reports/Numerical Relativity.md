# Numerical Relativity
## A Report on Foundations

Numerical relativity is a branch of general relativity that concerns theoretical and computational aspects of the theory. In this report I will explain the main concepts and general theoretical approaches.

Numerical Relativity [NR], is based on 3+1 decomposition, an alternative formulation of general relativity, which in principle "foliates" spacetime, and thus, separates time from space. This view in general is more pleasing in computational physics since we are more used to simulate or solve physical systems in a time evolution.

## Table of Contents
- Decomposing Spacetime
- Understanding The Construction of Metric and Normal Vector
- Construction of the Time Vector
- The Generic Line Element in 3+1 Decomposition
- Intrinsic and Extrinsic Curvature of Hyper-surfaces
- Decomposing Field Equations
- Einstein Tensor
- ADM-York Formulation
- Bianchi Identities
- How to Find Solutions
- Initial Data: The Elliptic Problem
- Main Approaches
- Conformal Decomposition
- Conclusion
- Further Options

## Decomposing Spacetime
### Understanding The Construction of Metric and Normal Vector

Einstein's theory of general relativity holds all the dimensions, be it time, or space dimensions, with the same regard, making any distinction between space and time coordinates more a matter of convention than a fundamental aspect of the theory. However, both our experiences and the laws of physics—especially on large scales—indicate that differentiating time from spatial coordinates is the most natural approach for describing physical processes[1].

By following this approach, we foliate spacetime with a series of non-overlapping spacelike hypersurfaces $\Sigma(t)$, each defined by a constant coordinate value $t$. This effectively separates the three spatial coordinates from the single temporal coordinate, a method known as the 3+1-decomposition of spacetime[1].

Assuming a hypersurface of such, we can introduce a four-vector normal to the hypersurface at each event (which is timelike because of orthogonality). We also would want it's dual one-form parallel to the gradient of the coordinate $t$.

### Normal Vector Normalization

The normal vector $n$ is defined with the normalization condition:

$$n^a n_a = g_{00} = -\frac{1}{2\alpha^2}$$

Where $\alpha$ is the lapse function.

## Construction of the Time Vector

The time-vector is introduced to carry out time evolutions and is defined as:

$$t^a = e^a_t = n^a + \beta^a$$

Where $\beta^a$ is the shift vector.

## The Generic Line Element in 3+1 Decomposition

The line element is defined as:

$$ds^2 = -(\alpha^2 - \beta_i\beta^i) dt^2 + 2\beta_i dx^i dt + \gamma_{ij} dx^i dx^j$$

## Intrinsic and Extrinsic Curvature of Hyper-surfaces

### Intrinsic Curvature

The intrinsic curvature is the three-dimensional Riemann tensor:

$${}^{(3)}R_{ilkm} = \gamma_{il}\gamma_{km} - \gamma_{ik}\gamma_{lm} + \gamma_{in}\gamma_{ln}\gamma_{km} - \gamma_{im}\gamma_{nn}\gamma_{lk}$$

### Extrinsic Curvature

The extrinsic curvature is given by:

$$K_{ij} = \frac{1}{2}(\partial_i\beta_j + \partial_j\beta_i - \alpha \partial_i\partial_j)$$

## Decomposing Field Equations

### Constraint Equations

**Hamiltonian Constraint:**
$$n^n G_{\mu\nu} = \frac{1}{2}H = R - K_{ij}K^{ij} + K^2 = 0$$

**Momentum Constraint:**
$$-i^n G_{\mu\nu} = M_i = D_j K^{ij} - D_i K = 0$$

**Dynamic Evolution Equations:**
$${}^{ij}G_{\mu\nu} = {}^{ij} = L_n K_{ij} + D_i D_j - (R_{ij} + K K_{ij} - 2K_{ik}K^{k}_{j}) = 0$$

## Conformal Decomposition

The main idea is to decompose the spatial metric into a conformal factor and metric:

$$\gamma_{ij} = \Psi^4 \bar{\gamma}_{ij}$$

The extrinsic curvature is split into:

$$K_{ij} = A_{ij} + \frac{1}{3}\gamma_{ij}K$$

Where:
- $A_{ij}$ is the trace-free part of the extrinsic curvature
- $K = \gamma^{ij}K_{ij}$ is the trace of the extrinsic curvature

## Conclusion

Numerical relativity has both deep theoretical and numerical problems. The main approach today is to decompose the Einstein tensor, and given other equations (such as the laws of evolution for matter), constraints and ansatz we form a basic process to the analysis of general relativity, numerically.

## References

[1] Y. Mizuno and L. Rezzolla, "General-Relativistic Magnetohydrodynamic Equations: the bare essential," Apr. 22, 2024, arXiv: arXiv:2404.13824