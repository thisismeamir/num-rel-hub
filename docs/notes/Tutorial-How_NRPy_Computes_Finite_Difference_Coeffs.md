## Example: Centered Finite Difference Representation of $u'(x_0) = u'_0$ Accurate to Fourth-order in $\Delta x$

As an illustration, let's first derive for uniform grids the centered, first-order, finite-difference coefficients accurate to fourth-order in $\Delta x$. The fourth-order-accurate, uniformly sampled, centered finite-difference derivative $u'(x_0)$ is equivalent to the derivative of the unique polynomial that passes through $u(x)$ at sampled points $\left\{ x_{-2},x_{-1},x_{0},x_{1},x_{2} \right\}$, where $x_i=x_0 + i \Delta x$.

The Taylor series expansion of a function $u(x)$ about a point $x_0$ is given by

$$u(x) = \sum_{n=0}^\infty \frac{u^{(n)}(x_0)}{n!} (x-x_0)^n,$$

where $u^{(n)}(x_0)$ is the $n$th derivative of $u$ evaluated at point $x_0$. Based on this, we can immediately write the Taylor expansion of $f$ at a point $x=x_0+j\Delta x$. In this case,

\begin{align}
u(x_0+j\Delta x) &= \sum_{n=0}^\infty \frac{u^{(n)}(x_0)}{n!} (j\Delta x)^n,\text{ or equivalently:} \\
u_j &= \sum_{n=0}^\infty \frac{u^{(n)}_0}{n!} (j\Delta x)^n.
\end{align}

Our goal is to compute $u^{(1)}(x_0)=u'_0$ at some point $x_0$, with a dominant error term proportional to $(\Delta x)^4$. We accomplish this by Taylor expanding $u(x_j)$ about $x=x_0$ for $j\in \left\{-2,-1,0,1,2\right\}$, each up to the $n=5$ term.

\begin{align}
u_{-2} &= u_0 - (2 \Delta x) u'_0 + \frac{(2 \Delta x)^2}{2} u''_0 - \frac{(2 \Delta x)^3}{3!} u'''_0 + \frac{(2 \Delta x)^4}{4!} u^{(4)}_0 +\mathcal{O}\left((\Delta x)^5\right) \\
u_{-1} &= u_0 - (\Delta x) u'_0 + \frac{(\Delta x)^2}{2} u''_0 - \frac{(\Delta x)^3}{3!} u'''_0 + \frac{(\Delta x)^4}{4!} u^{(4)}_0 +\mathcal{O}\left((\Delta x)^5\right)\\
u_{0} &= u_0 \\
u_{1} &= u_0 + (\Delta x) u'_0 + \frac{(\Delta x)^2}{2} u''_0 + \frac{(\Delta x)^3}{3!} u'''_0 + \frac{(\Delta x)^4}{4!} u^{(4)}_0 +\mathcal{O}\left((\Delta x)^5\right)\\
u_{2} &= u_0 + (2 \Delta x) u'_0 + \frac{(2 \Delta x)^2}{2} u''_0 + \frac{(2 \Delta x)^3}{3!} u'''_0 + \frac{(2 \Delta x)^4}{4!} u^{(4)}_0 +\mathcal{O}\left((\Delta x)^5\right)\\
\end{align}

Let's combine the above equations to find coefficients $a_j$ such that $(a_{-2} u_{-2} + a_{-1} u_{-1}...)/(\Delta x) = u'_0 + \mathcal{O}\left((\Delta x)^4\right)$.

\begin{align}
& (a_{-2} u_{-2} + a_{-1} u_{-1} + a_0 u_0 + a_{1} u_{1} +a_{2} u_{2})/(\Delta x) \\
= & \left( u_0 - (2 \Delta x) u'_0 + \frac{(2 \Delta x)^2}{2} u''_0 -\frac{(2 \Delta x)^3}{3!} u'''_0+\frac{(2 \Delta x)^4}{4!} u^{(4)}_0 \right) a_{-2} \\
& + \left( u_0 - (\Delta x) u'_0 + \frac{(\Delta x)^2}{2} u''_0 - \frac{(\Delta x)^3}{3!} u'''_0+\frac{(\Delta x)^4}{4!} u^{(4)}_0 \right) a_{-1} \\
& + \left( u_0 \right) a_{0} \\
& + \left( u_0 + (\Delta x) u'_0 + \frac{(\Delta x)^2}{2} u''_0 + \frac{(\Delta x)^3}{3!} u'''_0+\frac{(\Delta x)^4}{4!} u^{(4)}_0 \right) a_{1} \\
& + \left( u_0 + (2 \Delta x) u'_0 + \frac{(2 \Delta x)^2}{2} u''_0 + \frac{(2 \Delta x)^3}{3!} u'''_0+\frac{(2 \Delta x)^4}{4!} u^{(4)}_0 \right) a_{2}
\end{align}

First notice that each time we take a derivative in the Taylor
expansion, we multiply by a $\Delta x$. Notice that this helps to keep
the units consistent (e.g., if $x$ were in units of meters). Let's
just **absorb those $\Delta x$'s into the derivatives (we will extract them again later)** and rearrange terms.

\begin{align}
& a_{-2} u_{-2} + a_{-1} u_{-1} + a_0 u_0 + a_{1} u_{1} + a_{2} u_{2} \\
& = \left( a_{-2} + a_{-1} + a_0 + a_{1} + a_{2} \right) \times u_0 \\
& + \left( -2 a_{-2} - a_{-1} + a_{1} + 2 a_{2} \right) \times u'_0 \\
& + \left( 2^2 a_{-2} + a_{-1} + a_{1} + 2^2 a_{2} \right)/2! \times u''_0 \\
& + \left( -2^3 a_{-2} - a_{-1} + a_{1} + 2^3 a_{2} \right)/3! \times u'''_0 \\
= & u'_0
\end{align}

In order for the above to hold true for any nonzero values of
$\left\{ u_0,u'_0,u''_0,u'''_0,u^{(4)}_0\right\}$, the following set
of equations must also hold:
\begin{align}
0 &= a_{-2} + a_{-1} + a_0 + a_{1} + a_{2}\\
1 &= -2 a_{-2} - a_{-1} + a_{1} + 2 a_{2}\\
0 \times 2! &= 2^2 a_{-2} + a_{-1} + a_{1} + 2^2 a_{2}\\
0 \times 3! &= -2^3 a_{-2} - a_{-1} + a_{1} + 2^3 a_{2} \\
0 \times 4! &= 2^4 a_{-2} + a_{-1} + a_{1} + 2^3 a_{2}.
\end{align}

Now we write this expression in matrix form (note that $0!=1$).
\begin{equation}
\left(
\begin{array}{c}
0\times 0! \\
1\times 1! \\
0\times 2! \\
0\times 3! \\
0\times 4! \\
\end{array}
\right)
=
\left(
\begin{array}{ccccc}
 1 &  1 & 1 & 1 & 1 \\
(-2)^1 &(-1)^1 & 0 & 1 & 2 \\
(-2)^2 &(-1)^2 & 0 & 1 & 2^2 \\
(-2)^3 &(-1)^3 & 0 & 1 & 2^3 \\
(-2)^4 &(-1)^4 & 0 & 1 & 2^4 \\
\end{array}
\right)
\left(
\begin{array}{c}
a_{-2} \\
a_{-1} \\
a_{0} \\
a_{1} \\
a_{2} \\
\end{array}
\right)
\end{equation}

So we have reduced the computation of finite difference coefficients
to the inversion of an $N\times N$ matrix equation. Notice that the
elements of the matrix will vary from the one given above if the grid
spacing is not constant, but are otherwise invariant to $\Delta x$.

The inverted matrix reads
\begin{equation}
\left(
\begin{array}{ccccc}
0 & 1/12 & -1/24 & -1/12 & 1/24 \\
0 & -2/3 & 2/3 & 1/6 & -1/6 \\
1 & 0 & -5/4 & 0 & 1/4 \\
0 & 2/3 & 2/3 & -1/6 & -1/6 \\
0 & -1/12 & -1/24 & 1/12 & 1/24 \\
\end{array}
\right)
\label{fourthorder_inv_matrix}.
\end{equation}

The coefficients for the $M$th derivative can be immediately read by
multiplying the $(M+1)$st column by $M!/(\Delta x)^M$. For example, the zeroth derivative at point $x_0$ is given by 

$$\frac{0!}{(\Delta x)^0} \times (0 u_{-2} + 0 u_{-1} + u_0 + 0 u_{1} + 0 u_{2}) = u_0,$$

which is exact. The first derivative finite difference approximation at point $x_0$ is given by

$$\frac{1!}{(\Delta x)^1} \times \left(\frac{1}{12}( u_{-2} - u_{2}) + \frac{2}{3}( -u_{-1} + u_{1})\right) \approx (\partial_x u)_0,$$

and the second derivative finite difference approximation at point $x_0$ is given by 

$$\frac{2!}{(\Delta x)^2} \times \left(-\frac{1}{24}(u_{-2} + u_{2}) + \frac{2}{3}(u_{-1} + u_{1}) - \frac{5}{4} u_0 \right) \approx (\partial_x^2 u)_0.$$

In short, this matrix yields the finite difference derivative coefficients with the lowest possible error given a stencil size of 5 gridpoints. It can be shown by analyzing cancellations in higher order terms of the Taylor series expansions that the first and second derivative coefficients are correct to $(\Delta x)^4$ and the third and fourth derivatives are correct to $(\Delta x)^2$. 

### Exercise 1: Find the exact expressions for the dominant error term on all derivatives that can be computed from this matrix (zeroth through fourth derivatives).

### Exercise 2: Construct the matrix whose inverse yields the 5-point stencil *upwinded* derivative coefficients (i.e., the stencil includes points $\{u_{-4},u_{-3},u_{-2},u_{-1},u_{0}\}$).

NRPy+ implements this simple matrix inversion strategy to evaluate finite difference coefficients.

<!-- possibly unnecessary. added in case. -->
## Output this notebook to $\LaTeX$-formatted PDF file

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs.pdf](Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs")
```
