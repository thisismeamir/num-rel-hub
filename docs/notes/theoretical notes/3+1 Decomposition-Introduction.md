# Overview
3+1 decomposition is an alternative view on the equations of general relativity. This view is the basics of numerical relativity and hence it's the most important basis of understanding the algorithms, methods, and simulations that we do.

In this small note we'll investigate and learn it with detailed mathematics. This can be an easy starting point for and undergraduate that has a basic understanding of [[General Relativity Basics]].

# Hamiltonian Formulation
The normal approaches to general relativity, as a field theory is based on Lagrangian formulation, for some reasons this approach is not necessarily good for numerical relativity. Therefore, a Hamiltonian formulation is preferred. 

The Hamiltonian formulation of general relativity requires a separation of time and space coordinates, known as a $3+1$ decomposition. 

## Formulation
We start off by decomposing the metric as below:
$$
\begin{align}
g_{00} &= -\alpha^2 + \gamma^{ij}\beta_i\beta_j\\
g_{0i} &= \beta_i\\
g_{ij} &= \gamma_{ij} 
\end{align}
$$
This decomposition of the metric replaces the 10 independent metric components by the **lapse function** $\alpha(x)$, the shift vector $\beta_i(x)$, and the symmetric spatial metric $\gamma_{ij}(x)$. The inverse spacetime metric components are:
$$
\begin{align}
g^{00} &= -\frac{1}{\alpha^2}\\
g^{0i} &= \frac{\beta^i}{\alpha^2}\\
g^{ij} &= \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
\end{align}$$
The $3+1$ decomposition separates the treatment of time and space coordinates. In place of four dimensional gradients, we use time derivatives and three-dimensional gradients. Thus, as an example:
$$
\nabla_j A^i = \partial_j A ^i + \gamma^i_{jk}A^{k}
$$
where:
$$
\gamma^i_{jk} \equiv \frac12 \gamma^{il} \left(\partial_j\gamma_{kl} + \partial_k \gamma_{jl} + \partial_l \gamma_{jk}\right)
$$

## Intrinsic and Extrinsic Curvature pf Hyper-surfaces
In $3+1$ decomposition spacetime is described as a set of three-dimensional hyper-surfaces of constant time. These hyper-surfaces have intrinsic curvature, given by the three dimensional Riemann tensor:
$$
\!^{(3)}R^{i}_{lkm} = \partial_k \gamma^i_{kl} -\partial_m\gamma^i_{kl} +\gamma^i_{kn}\gamma^n_{lm} - \gamma^i_{mn}\gamma^n_{kl}
$$
In addition to the intrinsic curvature, the hyper-surfaces of constant time has an extrinsic curvature $K_{ij}$ arising from its embedding in four-dimensional spacetime.
$$
K_{ij} = \frac1{2\alpha}(\nabla_i\beta_j + \nabla_j \beta_i - \partial_t\gamma_{ij})
$$
The full spacetime (4 dimensional) curvature is then obviously related to the extrinsic and intrinsic curvature of the hyper-surfaces. This relation has been established using the [[Gauss-Codazzi Equations]].
$$
\!^{(4)}R^0_{\cdot jkl} = -\frac1\alpha \left(\nabla_k K_{jl} -\nabla_l K_{jk}\right)
$$
and 
$$
\!^{(4)}R^i_{\cdot jkl} = \!^{(3)}R^i_{\cdot jkl} - \!^{(4)}R^0_{\cdot jkl}\beta^i + K^i_{\cdot k} K_{jl} - K^i_{\cdot l}K_{jk}. 
$$
The other components of the four-dimensional Riemann tensor are:
$$
\!^{(4)}R^0_{\cdot i0j} = -\frac1\alpha \partial_t K_{ij} - K_i^{\cdot k}K_{jk}-\frac1\alpha\nabla_i\nabla_j\alpha\left[\nabla_j\left(\beta^kK_{ik}\right) + K_{jk}\nabla_i\beta^k\right]
$$
and
$$

\!^{(4)}R^k_{\cdot i0j} = 
\!^{(4)}R^k_{\cdot ilj}\beta^l + \left[
\!^{(4)}R^0_{\cdot ilj}\beta^l - 
\!^{(4)}R^0_{\cdot i0j}
\right]\beta^k + \alpha\left[
\nabla^k K_{ij} -\nabla_i K^k_{\cdot j}
\right]
$$
These equations can be combined to give:
$$
\begin{align}
\!^{(4)}R_{ijkl} =\  &\!^{(3)}R_{ijkl} + K_{ik}K_{jl} - K_{jl}K_{jk},\\
\!^{(4)}R_{0jkl} =\ &\!^{(4)}R_{ijkl}\beta^i +\alpha (\nabla_k K_{jl} - \nabla_lK_{jk})= \!^{(4)}R_{kl0j}\\
\!^{(4)}R^{0j}_{\cdot \cdot kl} = \ &\alpha\partial_t K_{ij}+\alpha^2 K_i^{\cdot k}K_{jk} + \alpha\nabla_i\nabla_j \alpha + \\ &\alpha\beta^k\nabla_kK_{ij} - \alpha\nabla_i(\beta^k K_{jk})- \alpha\nabla_j(\beta^kK_{ik})+ \!^{(4)}R_{kilj}\beta^k\beta^l\end{align}
$$
and 
$$
\begin{align}
\!^{(4)}R^{ij}_{kl} &= \!^{(3)}R^{ij}_{kl}+ K^i_{\cdot k} K^j_{\cdot l} - K^{i}_{\cdot l}K^j_{\cdot k} + \frac4\alpha \beta^{[u}\nabla_{[k}K^{j]}_{\cdot l]},\\
\!^{(4)}R^{0j}_{\cdot \cdot kl} &= -\frac1\alpha (\nabla_k K^j_{\cdot l} - \nabla_l K^j_{\cdot k})\\
\!^{(4)}R^{0i}_{0j} &= -\frac 1\alpha\partial_t K^i_{\cdot j} +K^i_{\cdot k} K^k_{\cdot j} - \frac1\alpha \nabla^i\nabla_j \alpha + \frac1\alpha \nabla_j (\beta^k K^i_{\cdot k}) - \frac1\alpha (\nabla_k\beta^i)K^k_{\cdot j}
\end{align}
$$




---
#introductory #theory 

###### References:
[Hamiltonian Formulation of General Relativity](https://web.mit.edu/edbert/GR/gr11.pdf)
[Gauss Codazzi Equations](https://en.wikipedia.org/wiki/Gauss%E2%80%93Codazzi_equations)

