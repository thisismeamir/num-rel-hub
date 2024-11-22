from manim import *
from functions import *
from wrapper import *
from title import TitleWithImage
class Presentation(Scene):
    def construct(self):
        TitleWithImage.construct(self)
        self.wait(4)
        self.clear()
        
        slide2 = create_content_scene("OVER\nVIEW","3+1 decomposition is an alternative view on the equations of general relativity. This view is the basics of numerical relativity and hence it's the most important basis of understanding the algorithms, methods, and simulations that we do. ").construct(self)
        self.wait(4)
        self.clear()
        
        slide3 = create_content_scene("INTRO-\nDUCTION","The intrinsically \"covariant view\" of Einstein's theory of general relativity is based on the concept that all coordinates are equivalent and, hence, the distinction between spatial and time coordinates is more an organizational matter than a strict requirement of the theory. Yet, our experience, and the laws of physics on sufficiently large scales, do suggest that a distinction between the time coordinate and from the spatial ones is the most natural one in describing physical processes.").construct(self)
        self.wait(4)
        self.clear()

        slide4 = create_latex_content_scene("DECOM-\nPOSITION", r"""n_\mu =  A\Omega_\mu &= \nabla_\mu t \\ \Omega &:=\nabla_t""", "Given such constant-time hypersurface, we can introduce a timelike four-vector normal to the hypersurface at each event in spacetime and such that its dual one-form is parallel to the gradient of the time coordinate.").construct(self)
        self.wait(4)
        self.clear()

        slide5 = create_latex_content_scene("",r"""
            n^\mu n_\mu &= g^{\mu\nu} n_\mu n_\nu \\
            &= g^{00} A^2 \\
            &= -\frac1{\alpha^2}A^2,\\
            -\frac{1}{\alpha^2}A^2 &= -1\\
            A &= \pm\alpha\\
            """, "If we now require that the four-vector defines an observer and thus that it measures the corresponding four-velocity, then from the normalization condition on timelike four-vectors we find the value of the constant.").construct(self)
        
        self.wait(4)
        self.clear()

        slide6 = create_latex_content_scene("",r"""
\left\{
\begin{matrix}+ & \text{ Past Direction} \\
- & \text{ Future Direction}
\end{matrix}
\right.
""", "And thus we choose the negative sign to showcase the future.").construct(self)
        
        self.wait(4)
        self.clear()

        slide7 = create_latex_content_scene("SPATIAL\nMETRIC", r"""
\begin{matrix}
\gamma_{\mu\nu} = g_{\mu\nu} + n_\mu n_\nu , & \gamma^{\mu\nu} = g^{\mu\nu} + n^\mu n^\nu
\end{matrix}
""","Using this normal vector we can define the spatial metric for the hypersurface:").construct(self)
        self.wait(4)
        self.clear()

        slide8 = create_content_scene("DECOMPOSING\nTOOLKIT","The spatial metric and the hypersurface vector provide us a good toolkit to decompose any tensors into a purely spatial part and a purely timelike part.").construct(self)
        self.wait(2)
        self.clear()
        slide9 = create_latex_content_scene("SPATIAL\nPART", r"""
\gamma^{\mu}_{\cdot\nu} := g^{\mu\alpha}\gamma_{\alpha\nu}=g^\mu_{\cdot\nu} + n^\mu n_\nu = \delta^{\mu}_\nu + n^\mu n_\nu
""","To obtain the spatial part we contract with the spatial projection operator:").construct(self)
        
        self.wait(2)
        self.clear()
        slide10 = create_latex_content_scene("TIMELIKE\nPART", r"""
N^\mu_{\cdot \nu} := -n^\mu n_\nu
""","To obtain the spatial part we contract with the spatial projection operator:").construct(self)
        self.wait(2)
        self.clear()
        slide11 = create_latex_content_scene("GENERAL\nCASE",r"""U^\mu = \gamma^\mu_{\cdot \nu}U^\nu + N^{\mu}_{\cdot \nu}U^\nu""","Note that the spatial part is still a four-vector but now whereas it has the covariant time component, which is non-zero in general. Analogous considerations can be done about tensors of any rank.").construct(self)
        self.wait(4)
        self.clear()
        slide12 = create_latex_content_scene("TIME\nVECTOR", r"""n^\mu \Omega_\mu &= \frac{1}{A} n^\mu n_\mu\\
        &= \frac{1}{\alpha} \not= 1""", "We already know that the unit normal $n$ to a spacelike hypersurface  does not represent the direction along which the time coordinate changes, that is, it is not the direction of the time derivative. Indeed, if we compute the contraction of the two tensors we get a non-unit value:").construct(self)
        self.wait(4)
        self.clear()

        slide13 = create_latex_content_scene("TIME\nVECTOR", r"""t = e_t = \partial_t := \alpha n + \beta""", "Thus we introduce a time-vector, along which to carry out the time evolutions and that is dual to the surface one-form. Such a vector is just the time coordinate basis vector and is defined as the linear superposition of a purely temporal part and of a purely spatial one, namely:").construct(self)
        self.wait(4)
        self.clear()

        slide14 = create_content_scene("SHIFT\nVECTOR", "Here the vector is a purely spatial vector and it's usually referred to as the shift vector and will be another building block of the metric in 3+1 decomposition.").construct(self)
        self.wait(4)
        self.clear()
        slide9 = create_latex_content_scene("COORDINATE\nBASIS", r"""t^\mu \Omega_\mu &= \alpha n^\mu \Omega_\mu + \beta^\mu \Omega_\mu \\
        &= \frac{\alpha}{\alpha} = 1""", "We can check that $t$ is a coordinate basis vector by verifying that:").construct(self)
        self.wait(4)
        self.clear()

        slide10 = create_latex_content_scene("COMPONENTS", r"""\begin{matrix} 
        n_\mu= (-\alpha, 0,0,0), & n^\mu = \frac{1}{\alpha} (1, -\beta)
        \end{matrix}""", "Using the components of $n$, we can express the generic line element in the 3+1 decomposition:").construct(self)
        self.wait(4)
        self.clear()

        slide11 = create_latex_content_scene("LINE\nELEMENT", r"""ds^2 = -(\alpha^2 -\beta_i\beta^i)dt^2 + 2\beta_i dx^i dt + \gamma_{ij}dx^i dx^j""", "The line element clearly shows that to measure the proper time we just need $\\beta^i = 0$.").construct(self)
        self.wait(4)
        self.clear()

        slide12 = create_latex_content_scene("PROPER\nTIME", r"""d\tau^2 = \alpha^2(t, x^j)dt^2""", "While the shift vector measures the change of coordinates of a point from one to another hypersurface:").construct(self)
        self.wait(4)
        self.clear()

        slide13 = create_latex_content_scene("COORDINATE\nCHANGE", r"""x^i_{t+dt}= x^i_t - \beta^i(t,x^j)dt""", "This can be related to the metric covariant and contravariant components:").construct(self)
        self.wait(4)
        self.clear()

        slide14 = create_latex_content_scene("METRIC\nCOMPONENTS", r"""\begin{matrix}
        g_{\mu\nu} = \begin{pmatrix}-\alpha^2 + \beta_i \beta^i & \beta_i\\ 
        \beta_i & \gamma_{ij}\end{pmatrix}, & 
        g^{\mu\nu} = \begin{pmatrix}-\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
        \frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
        \end{pmatrix}
        \end{matrix}""", "An important identity is then derivable from this:").construct(self)
        self.wait(4)
        self.clear()

        slide15 = create_latex_content_scene("DETERMINANT\nIDENTITY", r"""\sqrt{-g} = \alpha\sqrt{\gamma}""", "The unit timelike normal $n$ can be associated to the four-velocity of a special class of observers, which are referred to as normal or Eulerian Observers.").construct(self)
        self.wait(4)
        self.clear()

        slide16 = create_latex_content_scene("METRIC\nDECOMPOSITION", r"""
        g_{00} &= -\alpha^2 + \gamma^{ij}\beta_i\beta_j\\
        g_{0i} &= \beta_i\\
        g_{ij} &= \gamma_{ij} 
        """, "We start off by decomposing the metric as below. This decomposition replaces the 10 independent metric components by the lapse function $\\alpha(x)$, the shift vector $\\beta_i(x)$, and the symmetric spatial metric $\\gamma_{ij}(x)$.").construct(self)
        self.wait(4)
        self.clear()

        slide17 = create_latex_content_scene("INVERSE\nMETRIC", r"""
        g^{00} &= -\frac{1}{\alpha^2}\\
        g^{0i} &= \frac{\beta^i}{\alpha^2}\\
        g^{ij} &= \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
        """, "The inverse spacetime metric components are:").construct(self)
        self.wait(4)
        self.clear()

        slide18 = create_latex_content_scene("3+1\nGRADIENTS", r"""\nabla_j A^i = \partial_j A ^i + \gamma^i_{jk}A^{k}""", "The 3+1 decomposition separates the treatment of time and space coordinates. In place of four dimensional gradients, we use time derivatives and three-dimensional gradients. Thus, as an example:").construct(self)
        self.wait(4)
        self.clear()

        slide19 = create_latex_content_scene("CHRISTOFFEL\nSYMBOLS", r"""\gamma^i_{jk} \equiv \frac{1}{2} \gamma^{il} \left(\partial_j\gamma_{kl} + \partial_k \gamma_{jl} + \partial_l \gamma_{jk}\right)""", "Where the three-dimensional Christoffel symbols are defined as:").construct(self)
        self.wait(4)
        self.clear()