# Introduction
Relativistic astrophysics studies the most energetic and violent astrophysical processes that are characterized by very high speeds, strong gravitational fields, very high temperatures and ultra-intense magnetic fields.

These conditions are met normally near neutron stars and black holes. Therefore, a truly general relativistic approach if necessary for an accurate description of the physical conditions.

In this context relativistic magnetohydrodynamics represents a very effective framework to describe the dynamics of macroscopic plasma in a relativistic regime, also when considering non-astrophysical scenarios, such as the collision of heavy ions.

In magnetohydrodynamic approximation, the plasma is treated as a macroscopic fluid that is coupled with electromagnetic fields that it produces with its dynamics or that may be present from external sources. In addition, the collision times between particles in astrophysical fluids are usually smaller than the other relevant timescales of the systems. Because of this scale difference, astrophysical plasma can be though of as a continuous medium with well-defined macroscopic average quantities such as velocity, density, and pressure.

This theory is often used to study the dynamics of systems such as the collision of two magnetized neutron stars and the ensuing gamma-ray burst, or the accretion and outflows onto a central compact object. Under these conditions, general relativistic effects become important if not dominant and, hence, the solution of the full set of General Relativistic Magnetohydrodynamics equations represent the only avenue to obtain an accurate and physically consistent description of the system.

# Covariant General Relativistic Magnetohydrodynamic equations
The GRMHD equations can be derived after imposing the [[continuity equation]] and [[energy momentum tensor conservation]], namely the [[Bianchi identities]]:

$$
\begin{align}
\nabla_\mu \tilde J^\mu = \nabla_\mu \rho u^\mu = 0 \\ \nabla_\mu T^{\mu\nu} =0
\end{align}
$$

Here $\nabla_\mu$ is the covariant derivative of the metric $g_{\mu\nu}$, $\rho$ is the rest-mass density and $u^\mu$ is the four-velocity of the fluid. We continue by introducing $h_{\mu\nu}$, the projector operator orthogonal to $u^\mu$.

$$
h_{\mu\nu} := u_\mu u_\nu + g_{\mu\nu}
$$
and such that
$$
h_{\mu\nu} u^\mu = 0
$$
So then it is realizable (How?) that the energy momentum tensor conservation is four equations representing the conservation of energy:
$$
u_\nu \nabla_\mu T^{\mu\nu}= 0
$$
and four equations representing the conservation of four-momentum:
$$
h_{\alpha\nu} = \nabla_\mu T^{\mu\nu}
$$
Note that the total energy-momentum is the linear combination of the contributions coming from the matter and from the electromagnetic fields:
$$
T^{\mu\nu} = T^{\mu\nu}_{\text{matter}} + T^{\mu\nu}_{\text{fields}}
$$
where:
$$
\begin{align}
T^{\mu\nu}_{\text{matter}} &:= \rho h u^\mu u^\nu + p g^{\mu\nu}
\\
T^{\mu\nu}_{\text{fields}} &:= F^\mu_\lambda F^{\lambda\nu} - \frac14(F^{\lambda\delta}F_{\lambda\delta})g^{\mu\nu}
\end{align}
$$
Here $h$ is the specific enthalpy which is:
$$
h:= 1 + \epsilon + \frac p\rho
$$
where $\epsilon$ is the specific internal energy and $p$ is the pressure of the fluid. Also $F^{\mu\nu}$ is the [[Faraday tensor]]. Because electric and magnetic fields are available within the system we must have  the corresponding conservation laws, namely, the Maxwell Equations:
$$
\nabla_\mu F^{\mu\nu} = \mathcal{J}^\mu
$$
and
$$
\nabla_\mu \!^*F^{\mu\nu}= 0
$$
Where $\mathcal{J}^\mu$ is the charge current density, and $\!^*F^{\mu\nu}$ is the dual of the Faraday tensor. The Faraday tensor $F^{\mu\nu}$ is constructed from the electric and magnetic fields, as measured in the generic frame having $U^\alpha$ as tangent vector:
$$
F^{\mu\nu} = U^\mu E^\nu - U^\nu E^\mu - \sqrt{-g} \eta^{\mu\nu\lambda\delta} U_\lambda B_\delta
$$
where $\eta^{\mu\nu\lambda\delta}$ is the fully anti-symmetric symbol and $g$ is the determinant of the spacetime four-metric. The dual Faraday tensor is then:
$$
\!^*F^{\mu\nu} := \sqrt{-g}\eta^{\mu\nu\lambda\delta}F_{\lambda\delta}
$$
which is written as:
$$
\!^* F^{\mu\nu} = U^\mu B^\nu - U^\nu B^\mu + \sqrt{-g}\eta^{\mu\nu\lambda\delta}U_\lambda E_\delta
$$
Most of GRMHD simulations to date have explored scenarios within the so-called "ideal GRMHD" where the electrical conductivity is assumed to be infinite.

This can be a very good approximation because in astrophysical plasmas, the conductivity is actually very large. Under these conditions, the electric charges are infinitely effective in canceling any electric flied, which are therefore, zero in the frame co-moving with the fluid:
$$
F^{\mu\nu}u_\nu = 0
$$
The main consequence of the condition is that the electric fields cease to be independent vector fields and can be obtained from simple algebraic expressions involving the fluid four-velocity and the magnetic fields. 

Specially after defining the electric and magnetic four-vectors in the fluid frame as:
$$\begin{align}
e^\mu &:= F^{\mu\nu}u_\nu\\
b^\mu &:= F^{\mu\nu}u_\nu 
\end{align}
$$
with the constraints that
$$
e^\mu = 0
$$
and that the co-moving magnetic field is fully spatial:
$$
u_\mu b^\mu = 0
$$
Under these conditions we can write the Faraday and its dual tensors as:
$$
\begin{align}
F^{\mu\nu} &= \sqrt{-g}\eta^{\mu\nu\lambda\delta} u_\lambda b_\delta\\
\!^*F^{\mu\nu} &= b^\mu u^\nu - b^\nu u^\mu
\end{align}
$$
We can write the total energy-momentum tensor in terms of the vectors $u^\mu$ and $b^\mu$ as:
$$
T^{\mu\nu} = \rho h_{\text{tot}} u^\mu u^\nu + p_{\text{tot}} g^{\mu\nu} - b^\mu b^\nu
$$
where $p_{\text{tot}} := p + \frac{b^2}{2}$.
###### A List Of Books for GRMHD Numerical Methods:
1. Dinshaw S. Balsara. Higher-order accurate space-time schemes for computational astrophysics—Part I: finite volume methods. Living Reviews in Computational Astrophysics, 3(1):2, December 2017.

2. R. J. Leveque. Finite Volume Methods for Hyperbolic Problems. Cambridge University Press, New York, 2002.

3. J. M. Mart´ı and E. M¨uller. Grid-based Methods in Relativistic Hydrodynamics and Magnetohydrodynamics. Living Reviews in Computational Astrophysics, 1, December 2015.

4. L. Rezzolla and O. Zanotti. Relativistic Hydrodynamics. Oxford University Press, Oxford, UK, 2013.

5. E. F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics. Springer-Verlag, third edition, 2009.


###### A List of papers and codes of GRMHD:
1. M. Anderson, E. W. Hirschmann, S. L. Liebling, and D. Neilsen. Relativistic MHD with adaptive mesh refinement. Class. Quantum Grav., 23:6503–6524, November 2006.

2. Peter Anninos, P. Chris Fragile, and Jay D Salmonson. Cosmos++: Relativistic magnetohydrodynamics on unstructured grids with local adaptive refinement. Astrophys. J., 635:723, 2005.

3. L. Ant´on, O. Zanotti, J. A. Miralles, J. M. Mart´ı, J. M. Iba´ ˜nez, J. A. Font, and J. A. Pons. Numerical 3+1 General Relativistic Magnetohydrodynamics: A Local Characteristic Approach. Astrophys. J., 637:296–312, January 2006.

4. L. Baiotti, I. Hawke, P. J. Montero, F. L¨offler, L. Rezzolla, N. Stergioulas, J. A. Font, and E. Seidel. Three-dimensional relativistic simulations of rotating neutron-star collapse to a Kerr black hole. Phys. Rev. D, 71(2):024035, January 2005.

5. D. Be´gue´, A. Pe’er, G. Q. Zhang, B. B. Zhang, and B. Pevzner. cuHARM: A New GPU-accelerated GRMHD Code and Its Application to ADAF Disks. Astrophys. J., Supp., 264(2):32, February 2023.

6. P. Cerda´-Dura´n, J. A. Font, L. Ant´on, and E. M¨uller. A new general relativistic magnetohydrodynamics code for dynamical spacetimes. Astron. Astrophys., 492:937–953, December 2008.

7. Patrick Chi-Kit Cheong, Alan Tsz-Lok Lam, Harry Ho-Yin Ng, and Tjonnie Guang Feng Li. Gmunu: paralleled, grid-adaptive, general-relativistic magnetohydrodynamics in curvilinear geometries in dynamical space-times. Monthly Notices of the Royal Astronomical Society, 508(2):2279–2301, 09 2021.

8. Federico Cipolletta, Jay Vijay Kalinani, Bruno Giacomazzo, and Riccardo Ciolfi. Spritz: a new fully general-relativistic magnetohydrodynamic code. Class. Quant. Grav., 37(13):135010, 2020.

9. J.-P. De Villiers and J. F. Hawley. A Numerical Method for General Relativistic Magnetohydrodynamics. Astrophys. J., 589:458–480, May 2003.

10. L. Del Zanna, O. Zanotti, N. Bucciantini, and P. Londrillo. ECHO: a Eulerian conservative high-order scheme for general relativistic magnetohydrodynamics and magnetodynamics. Astron. Astrophys., 473:11–30, October 2007.

11. M. D. Duez, Y. T. Liu, S. L. Shapiro, and B. C. Stephens. Relativistic magnetohydrodynamics in dynamical spacetimes: Numerical methods and tests. Phys. Rev. D, 72(2):024028, July 2005.

12. Z. B. Etienne, V. Paschalidis, R. Haas, P. Mo¨sta, and S. L. Shapiro. IllinoisGRMHD: an open-source, user-friendly GRMHD code for dynamical spacetimes. Class. Quantum Grav., 32(17):175009, September 2015.

13. Charles F. Gammie, Jonathan C. McKinney, and Ga´bor T´oth. HARM: A Numerical Scheme for General Relativistic Magnetohydrodynamics. Astrophys. J., 589(1):444–457, May 2003.

14. B. Giacomazzo and L. Rezzolla. WhiskyMHD: a new numerical code for general relativistic magnetohydrodynamics. Class. Quantum Grav., 24:235, June 2007.

15. J. F. Hawley, L. L. Smarr, and J. R. Wilson. A numerical study of nonspherical black hole accretion. I Equations and test problems. Astrophys. J., 277:296–311, February 1984.

16. S. Koide D.L. Meier K. Shibata T. Kudoh. General relativistic simulation of early jet formation in a rapidly rotating black hole magnetosphere. Astrophys. J, 536:668–674, 2000.

17. Matthew Liska, Koushik Chatterjee, Alexand er Tchekhovskoy, Doosoo Yoon, David van Eijnatten, Casper Hesp, Sera Markoff, Adam Ingram, and Michiel van der Klis. H-AMR: A New GPU-accelerated GRMHD Code for Exascale Computing With 3D Adaptive Mesh Refinement and Local Adaptive Time-stepping. arXiv e-prints, page arXiv:1912.10192, December 2019.

18. Y. Mizuno, K.-I. Nishikawa, S. Koide, P. Hardee, and G. J. Fishman. RAISHIN: A HighResolution Three-Dimensional General Relativistic Magnetohydrodynamics Code. ArXiv Astrophysics e-prints, August 2006.

19. Hector Olivares, Oliver Porth, Jordy Davelaar, Elias R. Most, Christian M. Fromm, Yosuke Mizuno, Ziri Younsi, and Luciano Rezzolla. Constrained transport and adaptive mesh refinement in the Black Hole Accretion Code. Astron. Astrophys., 629:A61, September 2019.

20. O. Porth, H. Olivares, Y. Mizuno, Z. Younsi, L. Rezzolla, M. Moscibrodzka, H. Falcke, and M. Kramer. The black hole accretion code. Computational Astrophysics and Cosmology, 4:1, May 2017.

21. M. Shibata and Y.-I. Sekiguchi. Magnetohydrodynamics in full general relativity: Formulation and tests. Phys. Rev. D, 72(4):044014, August 2005.

22. C. J. White, J. M. Stone, and C. F. Gammie. An Extension of the Athena++ Code Framework for GRMHD Based on Advanced Riemann Solvers and Staggered-mesh Constrained Transport. Astrophys. J.s, 225:22, August 2016.

23. O. Zanotti and M. Dumbser. A high order special relativistic hydrodynamic and magnetohydrodynamic code with space-time adaptive mesh refinement. Computer Physics Communications, 188:110–127, March 2015.
24. 29. J. A. Font. Numerical hydrodynamics in general relativity. Living Rev. Relativ., 6:4, 2003.

25. Bruno Giacomazzo and Luciano Rezzolla. The Exact Solution of the Riemann Problem in Relativistic MHD. Journal of Fluid Mechanics, 562:223–259, 2006.

26. J. M. Mart´ı and E. M¨uller. Grid-based Methods in Relativistic Hydrodynamics and Magnetohydrodynamics. Living Reviews in Computational Astrophysics, 1, December 2015.

27. L. Rezzolla and K. Takami. Black-hole production from ultrarelativistic collisions. Class. Quantum Grav., 30(1):012001, January 2013.