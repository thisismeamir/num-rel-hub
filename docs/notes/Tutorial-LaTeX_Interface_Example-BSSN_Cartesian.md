<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# LaTeX Interface Example: BSSN (Cartesian)

## Author: Ken Sible

## This notebook demonstrates parsing BSSN (Cartesian) and validating against NRPy+.
[comment]: <> (**Notebook Status:** <font color='red'><b> Not Validated </b></font>)

### Covariant Formulation of BSSN

\begin{align}
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
    \partial_t \phi &= \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
    \partial_t K &= \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j} \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] + \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_j \beta^j + \frac{1}{3} \bar{D}^i \bar{D}_j \beta^j \\
    &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \partial_j \phi \right) + 2 \alpha \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k} \bar{A}^k{}_j + \alpha \bar{A}_{i j} K \\
    &\qquad + e^{-4 \phi} \left\{-2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \right\}^{\text{TF}} \\
    \partial_t \alpha &= \left[\beta^k \partial_k \alpha \right] - 2 \alpha K \\
    \partial_t \beta^i &= \left[\beta^k \bar{D}_k \beta^i \right] + B^i \\
    \partial_t B^i &= \left[\beta^k \bar{D}_k B^i \right] + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^k \bar{D}_k \bar{\Lambda}^i \right] \right) - \eta B^i \\
    \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi \right) \\
    \mathcal{M}^i &= e^{-4\phi} \left(\hat{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij}\partial_j \phi - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K + \bar{A}^{jk} \Delta^i_{jk} + \bar{A}^{ik} \Delta^j_{jk} \right) \\
    \bar{R}_{i j} &= - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_k \hat{D}_l \bar{\gamma}_{i j} + \bar{\gamma}_{k(i} \hat{D}_{j)} \bar{\Lambda}^k + \Delta^k \Delta_{(i j) k} + \bar{\gamma}^{k l} \left(2 \Delta_{k(i}^m \Delta_{j) m l} + \Delta_{i k}^m \Delta_{m j l} \right)
\end{align}

<a id='top'></a>

# Table of Contents
$$\label{toc}$$

- [Step 1](#step_1): Evolution Equation for $\partial_t \bar{\gamma}_{ij}$
- [Step 2](#step_2): Evolution Equation for $\partial_t \phi$
- [Step 3](#step_3): Evolution Equation for $\partial_t K$
- [Step 4](#step_4): Evolution Equation for $\partial_t \bar{\Lambda}^i$
- [Step 5](#step_5): Evolution Equation for $\partial_t \bar{A}_{ij}$
- [Step 6](#step_6): Gauge Evolution Equation(s)
- [Step 7](#step_7): Constraint Equation(s)
- [Step 8](#step_8): Output Notebook to PDF


```python
!pip install nrpylatex==1.0.8 > /dev/null
!pip freeze | grep nrpylatex
%load_ext nrpylatex.extension
```

    nrpylatex==1.0.8



```python
from UnitTesting.assert_equal import assert_equal

import NRPy_param_funcs as par, reference_metric as rfm
import BSSN.BSSN_quantities as Bq
import BSSN.BSSN_RHSs as Brhs
import BSSN.BSSN_gauge_RHSs as gaugerhs
import BSSN.BSSN_constraints as bssncon

par.set_parval_from_str('reference_metric::CoordSystem', 'Cartesian')
par.set_parval_from_str('BSSN.BSSN_quantities::LeaveRicciSymbolic', 'True')
rfm.reference_metric()
Brhs.BSSN_RHSs()
gaugerhs.BSSN_gauge_RHSs()
bssncon.BSSN_constraints()
par.set_parval_from_str('BSSN.BSSN_quantities::LeaveRicciSymbolic', 'False')
Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
```

We formatted the notebook in the following pedagogical scheme, shown below, to encourage a *learned* approach to parsing tensorial LaTeX.
```
(1) attempt to parse equation
(2) if ParseError or TensorError
        debug error and goto (1)
(3) change equation and goto (1)
```
[YouTube Tutorial](https://www.youtube.com/watch?v=A9IWYpQe6Fo) (Notebook Walkthrough)

<a id='step_1'></a>

## Step 1: Evolution Equation for $\partial_t \bar{\gamma}_{ij}$ [ [^](#top) ]

\begin{equation}
    \partial_t \bar{\gamma}_{i j} = \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j}
\end{equation}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```

    ParseError: \partial_t \bar{\gamma}_{i j} = [\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
                           ^
    cannot index undefined tensor 'gammabarDD' at position 29


First, attempt to parse the equation without any modification, and notice the `ParseError` for indexing an `undefined tensor 'gammabarDD'`.

`gammabarDD :=` $ \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij} $ where `hDD` is the deviation from the flat 3D Cartesian metric `gammahatDD :=` $ \hat{\gamma}_{ij} = \delta_{ij} $.

`vardef`: define a variable, including a derivative option and a (anti)symmetry option.<br>
`assign`: assign a `vardef` option to an already existing variable in the namespace.<br>
`parse`:&nbsp;&nbsp; parse an equation without rendering that equation in a `.tex` document or a Jupyter Notebook.

`-metric` can assign the symmetry `sym01` and automatically generate the metric inverse and determinant.

1. Define the Kronecker delta `deltaDD` using the `vardef` macro.
1. Parse the equation $\hat{\gamma}_{ij} = \delta_{ij}$ using the `parse` macro.
1. Assign the `-metric` option to `gammahatDD` using the `assign` macro.

```
% vardef -kron 'deltaDD'
% parse \hat{\gamma}_{ij} = \delta_{ij}
% assign -diff_type=symbolic -metric 'gammahatDD'
```

**Remark:** If a derivative option should also apply to the metric inverse and determinant, then that option should precede `-metric`. Moreover, you need not define a scalar quantity. However, should a scalar quantity require a `vardef` option, you may assign one using either the `vardef` or `assign` macro.

1. Define the deviation `hDD` using the `vardef` macro.
1. Parse the equation $\bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}$ using the `parse` macro.
1. Assign the `-metric` option to `gammabarDD` using the `assign` macro.

```
% vardef -diff_type=dD -symmetry=sym01 'hDD'
% parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
% assign -diff_type=dD -metric 'gammabarDD'
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```

    ParseError: \partial_t \bar{\gamma}_{i j} = [\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
                                                                                                    ^
    cannot index undefined tensor 'betaU' at position 386


Next, replace every instance of `betaU` with `vetU`, a different variable name chosen for the SymPy output, using the `srepl` macro, and then define the variable `vetU` using the `vardef` macro. The `srepl` macro can preserve the original LaTeX representation while changing the variable name internally.

**Remark:** If you require a variable name longer than one Latin or Greek letter, you should surround that variable name in the `\text` command.

`srepl`: string replacement given the assumption that string `A` and string `B` are equal whenever they are equivalent syntactically

Therefore, for example, `S_{a b} T^{b c}` and `S_{ab}T^{bc}` are considered equal according to the `srepl` macro.

```
% srepl "\beta" -> "\text{vet}"
% vardef -diff_type=dD 'vetU'
```

**Remark**: You can check "syntactic equivalence" for the `srepl` macro using the code snippet provided below.

```
lexer = Lexer()
lexer.initialize(A)
syntax_A = [(lexer.lexeme, token) for token in lexer.tokenize()]
lexer.initialize(B)
syntax_B = [(lexer.lexeme, token) for token in lexer.tokenize()]
assert syntax_A == syntax_B
```

`-persist` will apply the replacement rule `"..." -> "..."` to every subsequent input (internal or external) of the `parse` function.

**Remark:** The `srepl` macro will *not* replace any text inside of a future `srepl` or `ignore` macro.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```

    ParseError: \bar{A}^k_k - \bar{D}_k \text{vet}^k) - 2 \alpha \bar{A}_{i j} \\
                ^
    cannot index undefined tensor 'AbarUD' at position 598


1. Replace every instance of `AbarDD` with `aDD` using the `srepl` macro.
1. Define the variable `aDD` (with symmetry `sym01`) using the `vardef` macro.
1. Generate the variable `aUD` by assigning the metric `gammabar` to `aDD`.

\begin{equation}
    \bar{A}^i{}_j = \bar{\gamma}^{ik} \bar{A}_{kj}
\end{equation}

```
% srepl "\bar{A}" -> "\text{a}"
% vardef -diff_type=dD -symmetry=sym01 'aDD'
% assign -metric='gammabar' 'aDD'
```

**Remark**: The default metric for a variable is the metric associated with the diacritic (or lack thereof) in the variable name. In the current situation, we need to explicitly associate that metric since we removed the `bar` diacritic from the variable `AbarDD` using the `srepl` macro.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in gammabarDD_dD


Next, replace the partial derivative $\partial_t \bar{\gamma}$ on the LHS with the variable name `h_rhs` using the `srepl` macro.

```
% srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
```

**Remark**: We received a `TensorError` that resulted from interpreting `t` as an index in the partial derivative $\partial_t \bar{\gamma}$ on the LHS. However, if the `basis` included the symbol `t`, then numeric indexing would have occurred, i.e. `t -> 0` provided that `basis = [t, x, y, z]` for example.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```




\[
\begin{align}
    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
\]



`keydef`: define a keyword, i.e. define a `basis` or coordinate system, or define an `index` range (for looping or summation).

```
% keydef basis [x, y, z]
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
\]




```python
try:
    assert_equal(h_rhsDD, Brhs.h_rhsDD)
except AssertionError:
    print('Assertion Failed!')
```

    Assertion Failed!


Finally, we upwind every partial derivative in each instance of the pattern `\beta^{...} \partial_{...}`. To enforce upwinding inside of the bracketed term on the RHS where that pattern occurs, use the `vphantom` command to dynamically change the derivative type to upwinded.

```
\partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k ...
```

**Remark**: You could, alternatively, replace the bracketed term with the Lie derivative `\mathcal{L}_\beta \bar{\gamma}_{i j}`, allow automatic generation of that Lie derivative, and then use an `srepl` macro with the `-persist` option and (advanced) capture group syntax to upwind the aforementioned pattern.

```
% srepl -persist "\text{vet}^<1> \partial_<1>" -> "\text{vet}^<1> \vphantom{dupD} \partial_<1>"
```

To simplify the calculation of $\bar{D}_k \beta^k$, substitute a useful tensor identity for that contraction and apply the `-persist` option.
$$ \bar{D}_k \beta^k = \partial_k \beta^k + \beta^k \frac{\partial_k \hat{\gamma}}{2 \hat{\gamma}} $$
```
%% replace '\bar{D}_k \beta^k' with contraction identity
% srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\
\end{align}
\]




```python
assert_equal(h_rhsDD, Brhs.h_rhsDD)
```

    Assertion Passed!


<a id='step_2'></a>

## Step 2: Evolution Equation for $\partial_t \phi$ [ [^](#top) ]

\begin{equation}
    \partial_t \phi = \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right)
\end{equation}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    \partial_t \phi &= \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in phi_dD


First, replace every instance of `phi` with the conformal factor `cf = W = e^{-2\phi}`, including every derivative (partial or covariant) of `phi`.

\begin{equation}
    \partial_t W = \partial_t e^{-2 \phi} = -2 e^{-2\phi} \partial_t \phi = -2 W \partial_t \phi \Rightarrow \partial_t \phi = - \frac{1}{2W} \partial_t W
\end{equation}

```
%% replace 'phi' with conformal factor cf = W = e^{-2\phi}
% srepl "e^{-4\phi}" -> "\text{cf}^2"
% srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
% srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
% srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"
```

**Remark**: We assigned the `-persist` option to the third `srepl` macro to replace the partial derivative inside of each covariant derivative of `phi`.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    \partial_t \phi &= \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    \partial_t \phi &= \left[\beta^k \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
\end{align}
\]




```python
try:
    assert_equal(cf_rhs, Brhs.cf_rhs)
except AssertionError:
    print('Assertion Failed!')
```

    Assertion Failed!


Next, enforce upwinding on the contracted partial derivative inside of the bracketed term using the `vphantom` command.

```
\partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k ...
```

Finally, replace every instance of `K` with `trK` using the `srepl` macro.

```
% srepl "K" -> "\text{trK}"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j} \left(\alpha
        \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\
\end{align}
\]




```python
assert_equal(cf_rhs, Brhs.cf_rhs)
```

    Assertion Passed!


<a id='step_3'></a>

## Step 3: Evolution Equation for $\partial_t K$ [ [^](#top) ]

\begin{equation}
    \partial_t K = \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right)
\end{equation}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    \partial_t K &= \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in trK_dD


First, replace the partial derivative $\partial_t K$ on the LHS with the variable name `trK_rhs` using the `srepl` macro.

```
% srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
\end{align}
\]




```python
try:
    assert_equal(trK_rhs, Brhs.trK_rhs)
except AssertionError:
    print('Assertion Failed!')
```

    Assertion Failed!


Next, enforce upwinding on the contracted partial derivative inside of the bracketed term using the `vphantom` command.

```
\partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K ...
```

Finally, assign the property `-diff_type=dD` to `cf`, `trK`, and `alpha` using the `assign` macro.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\
\end{align}
\]




```python
assert_equal(trK_rhs, Brhs.trK_rhs)
```

    Assertion Passed!


<a id='step_4'></a>

## Step 4: Evolution Equation for $\partial_t \bar{\Lambda}^i$ [ [^](#top) ]

\begin{align}
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] + \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_j \beta^j + \frac{1}{3} \bar{D}^i \bar{D}_j \beta^j \\
    &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \partial_j \phi \right) + 2 \alpha \bar{A}^{j k} \Delta_{j k}^i - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K
\end{align}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```

    ParseError: \partial_t \bar{\Lambda}^i = [\text{vet}^k \partial_k \bar{\Lambda}^i - \partial_k \text{vet}^i \bar{\Lambda}^k ] +
                           ^
    cannot index undefined tensor 'LambdabarU' at position 2099


First, replace every instance of `LambdabarU` with `lambdaU` using the `srepl` macro, and then define the variable `lambdaU` using the `vardef` macro.

```
% srepl "\bar{\Lambda}" -> "\text{lambda}"
% vardef -diff_type=dD 'lambdaU'
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in lambdaU_dD


1. Parse $ \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij} $ using the `parse` macro to define `DeltaUDD`.
1. Parse $ \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij} $ (contraction of `DeltaUDD`) using the `parse` macro to define `DeltaU`.

```
% parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
% parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij} \\
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in lambdaU_dD


Next, replace the partial derivative $\partial_t \bar{\Lambda}$ on the LHS with the variable name `Lambdabar_rhs` using the `srepl` macro.

```
% srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```

    ParseError: \qquad - 2 \text{a}^{i j} (\partial_j \alpha - 6 \alpha \partial_j \text{cf} \frac{-1}{2 \text{cf}} ) + 2 \alpha
                ^
    unsupported operator '\qquad' at position 2793


Next, internally remove the formatting command `\qquad` using the `ignore` macro.

`ignore`: remove a LaTeX command or substring; equivalent to `srepl "..." -> ""` (empty replacement).

```
% ignore "\qquad"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```

    ParseError: - 2 \text{a}^{i j} (\partial_j \alpha - 6 \alpha \partial_j \text{cf} \frac{-1}{2 \text{cf}} ) + 2 \alpha
                ^
    unsupported operator '-' at position 2816


Next, append a percent symbol `%` to the end of the line break preceding `\qquad`, and then remove that custom line break `\\%` using the `ignore` macro.

```
% ignore "\\%"
```

**Remark**: You cannot split an equation across a line break, and hence we removed that line break internally using the `ignore` macro.


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
\]




```python
try:
    assert_equal(Lambdabar_rhsU, Brhs.Lambdabar_rhsU)
except AssertionError:
    print('Assertion Failed!')
```

    Assertion Failed!


Finally, enforce upwinding on the contracted partial derivative inside of the bracketed term using the `vphantom` command.

```
\partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k ...
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\
\end{align}
\]




```python
assert_equal(Lambdabar_rhsU, Brhs.Lambdabar_rhsU)
```

    Assertion Passed!


<a id='step_5'></a>

## Step 5: Evolution Equation for $\partial_t \bar{A}_{ij}$ [ [^](#top) ]

\begin{align}
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k} \bar{A}^k{}_j + \alpha \bar{A}_{i j} K \\
    &\qquad + e^{-4 \phi} \left\{-2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 4 \bar{D}_{(i} \alpha \bar{D}_{j)} \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \right\}^{\text{TF}}
\end{align}

**Remark**: We must manually expand the trace-free term $ \{\ldots\}^{TF} $ and the symmetric term $ \bar{D}_{(i} \alpha \bar{D}_{j)} $ since neither notation is currenlty supported.

\begin{gather}
    \bar{D}_{(i} \alpha \bar{D}_{j)} = \frac{1}{2} \left(\bar{D}_i \alpha \bar{D}_j + \bar{D}_j \alpha \bar{D}_i \right) \\
    X_{i j} = -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} = X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    \partial_t \bar{A}_{i j} = \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k} \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j}
\end{gather}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
```

    ParseError: \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
                                                                                                         ^
    cannot index undefined tensor 'RbarDD' at position 3255


First, define the variable `RbarDD` (with symmetry `sym01`) using the `vardef` macro.

```
% vardef -diff_type=dD -symmetry=sym01 'RbarDD'
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
```

    TensorError: unbalanced free index {'t'} in aDD_dD


Next, replace the partial derivative $\partial_t \bar{A}_{i j}$ on the LHS with the variable name `a_rhs` using the `srepl` macro.

```
% srepl "\partial_t \bar{A}_{i j}" -> "\text{a_rhs}"
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
\]




```python
try:
    assert_equal(a_rhsDD, Brhs.a_rhsDD)
except AssertionError:
    print('Assertion Failed!')
```

    Assertion Failed!


Finally, enforce upwinding on the contracted partial derivative inside of the bracketed term using the `vphantom` command.

```
\partial_t \bar{A}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k ...
```


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\
\end{align}
\]




```python
assert_equal(a_rhsDD, Brhs.a_rhsDD)
```

    Assertion Passed!


<a id='step_6'></a>

## Step 6: Gauge Evolution Equation(s) [ [^](#top) ]

\begin{align}
    \partial_t \alpha &= \left[\beta^k \partial_k \alpha \right] - 2 \alpha K \\
    \partial_t \beta^i &= \left[\beta^k \bar{D}_k \beta^i \right] + B^i \\
    \partial_t B^i &= \left[\beta^k \bar{D}_k B^i \right] + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^k \bar{D}_k \bar{\Lambda}^i \right] \right) - \eta B^i
\end{align}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\

    % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
    \partial_t \alpha &= \left[\beta^k \vphantom{dupD} \partial_k \alpha \right] - 2 \alpha K \\

    % srepl "B" -> "\text{bet}"
    % vardef -diff_type=dD 'betU'
    % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
    \partial_t \beta^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j \beta^i \right] + B^i \\

    % vardef -const 'eta'
    % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
    \partial_t B^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j B^i \right] + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j \vphantom{dupD} \bar{D}_j \bar{\Lambda}^i \right] \right) - \eta B^i \\
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % srepl "\beta" -> "\text{vet}"
    % vardef -diff_type=dD 'vetU'
    %% replace '\bar{D}_k \beta^k' with contraction identity
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % srepl "\bar{A}" -> "\text{a}"
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % assign -metric='gammabar' 'aDD'

    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
    \partial_t \bar{\gamma}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\gamma}_{i j} + \partial_i \beta^k
        \bar{\gamma}_{k j} + \partial_j \beta^k \bar{\gamma}_{i k} \right] + \frac{2}{3} \bar{\gamma}_{i j}
        \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{i j} \\

    % assign -diff_type=dD 'cf'
    %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
    % srepl "e^{-4 \phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    % srepl "K" -> "\text{trK}"
    % assign -diff_type=dD 'trK'
    \partial_t \phi &= \left[\beta^k \vphantom{dupD} \partial_k \phi \right] + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % assign -diff_type=dD 'alpha'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \left[\beta^k \vphantom{dupD} \partial_k K \right] + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{i j}
        \bar{A}^{i j} - e^{-4 \phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi \right) \\

    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % vardef -diff_type=dD 'lambdaU'

    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \left[\beta^k \vphantom{dupD} \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k \right] +
        \bar{\gamma}^{j k} \hat{D}_j \hat{D}_k \beta^i + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
        &\qquad - 2 \bar{A}^{i j} \left(\partial_j \alpha - 6 \alpha \partial_j \phi \right) + 2 \alpha
        \bar{A}^{j k} \Delta^i_{j k} - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{i j} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi + 2 \bar{D}_i \alpha
        \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{i j} \\
    \hat{X}_{i j} &= X_{i j} - \frac{1}{3} \bar{\gamma}_{i j} \bar{\gamma}^{k l} X_{k l} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{i j} &= \left[\beta^k \vphantom{dupD} \partial_k \bar{A}_{i j} + \partial_i \beta^k \bar{A}_{k j} + \partial_j \beta^k
        \bar{A}_{i k} \right] - \frac{2}{3} \bar{A}_{i j} \bar{D}_k \beta^k - 2 \alpha \bar{A}_{i k}
        \bar{A}^k{}_j + \alpha \bar{A}_{i j} K + e^{-4 \phi} \hat{X}_{i j} \\

    % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
    \partial_t \alpha &= \left[\beta^k \vphantom{dupD} \partial_k \alpha \right] - 2 \alpha K \\

    % srepl "B" -> "\text{bet}"
    % vardef -diff_type=dD 'betU'
    % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
    \partial_t \beta^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j \beta^i \right] + B^i \\

    % vardef -const 'eta'
    % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
    \partial_t B^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j B^i \right] + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j \vphantom{dupD} \bar{D}_j \bar{\Lambda}^i \right] \right) - \eta B^i \\
\end{align}
\]




```python
assert_equal(alpha_rhs, gaugerhs.alpha_rhs)
```

    Assertion Passed!



```python
assert_equal(vet_rhsU, gaugerhs.vet_rhsU)
```

    Assertion Passed!



```python
assert_equal(bet_rhsU, gaugerhs.bet_rhsU)
```

    Assertion Passed!


<a id='step_7'></a>

## Step 7: Constraint Equation(s) [ [^](#top) ]

\begin{align}
    \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi \right) \\
    \mathcal{M}^i &= e^{-4\phi} \left(\hat{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij}\partial_j \phi - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K + \bar{A}^{jk} \Delta^i_{jk} + \bar{A}^{ik} \Delta^j_{jk} \right)
\end{align}


```python
%%parse_latex --reset --ignore-warning

\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % vardef -diff_type=dD 'vetU'
    % srepl "\beta" -> "\text{vet}"
    %% upwind pattern inside Lie derivative expansion
    % srepl -persist "\text{vet}^<1> \partial_<1>" -> "\text{vet}^<1> \vphantom{dupD} \partial_<1>"
    %% substitute tensor identity (see appropriate BSSN notebook)
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % vardef -diff_type=dD 'alpha'
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % srepl "\bar{A}" -> "\text{a}"
    % parse \bar{A}^i_j = \bar{\gamma}^{ik} \bar{A}_{kj}
    % assign -diff_type=dD 'aUD'
    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"

    \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij}
        + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right)
        - 2 \alpha \bar{A}_{ij} \\

    % vardef -diff_type=dD 'cf', 'trK'
    % srepl "K" -> "\text{trK}"
    %% replace 'phi' with conformal factor cf = W = e^{{-2\phi}}
    % srepl "e^{-4\phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    \partial_t \phi &= \mathcal{L}_\beta \phi
            + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % parse \bar{A}^{ij} = \bar{\gamma}^{jk} \bar{A}^i_k
    % assign -diff_type=dD 'aUU'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \mathcal{L}_\beta K
        + \frac{1}{3} \alpha K^2
        + \alpha \bar{A}_{ij} \bar{A}^{ij}
        - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

    % vardef -diff_type=dD 'lambdaU'
    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta_{ijk}  = \bar{\gamma}_{il} \Delta^l_{jk}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
            + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
            &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
            + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
        + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
        - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
    \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij}
            - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
            - 2 \alpha \bar{A}_{ik} \bar{A}^k_j
            + \alpha \bar{A}_{ij} K
            + e^{-4\phi} \hat{X}_{ij} \\

    % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
    \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

    % vardef -diff_type=dD 'betU'
    % srepl "B" -> "\text{bet}"
    % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
    \partial_t \beta^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j \beta^i\right] + B^i \\

    % vardef -const 'eta'
    % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
    \partial_t B^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j B^i\right]
        + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j \vphantom{dupD} \bar{D}_j \bar{\Lambda}^i\right]\right)
        - \eta B^i \\

    % parse \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
    % srepl "\bar{D}^2" -> "\bar{D}^i \bar{D}_i", "\mathcal{<1>}" -> "<1>"
    \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi
        \bar{D}_i \phi - 8 \bar{D}^2 \phi \right) \\

    \mathcal{M}^i &= e^{-4\phi} \left(\hat{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij}\partial_j \phi - \frac{2}{3} \bar{\gamma}^{ij}
        \partial_j K + \bar{A}^{jk} \Delta^i_{jk} + \bar{A}^{ik} \Delta^j_{jk} \right) \\

    \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij} + \frac{1}{2} \left(\bar{\gamma}_{ki}
        \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right) + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
        &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
\end{align}
```




\[
\begin{align}
    % keydef basis [x, y, z]
    % ignore "\\%", "\qquad"

    % vardef -kron 'deltaDD'
    % parse \hat{\gamma}_{ij} = \delta_{ij}
    % assign -diff_type=symbolic -metric 'gammahatDD'
    % vardef -diff_type=dD -symmetry=sym01 'hDD'
    % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
    % assign -diff_type=dD -metric 'gammabarDD'

    % vardef -diff_type=dD 'vetU'
    % srepl "\beta" -> "\text{vet}"
    %% upwind pattern inside Lie derivative expansion
    % srepl -persist "\text{vet}^<1> \partial_<1>" -> "\text{vet}^<1> \vphantom{dupD} \partial_<1>"
    %% substitute tensor identity (see appropriate BSSN notebook)
    % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

    % vardef -diff_type=dD 'alpha'
    % vardef -diff_type=dD -symmetry=sym01 'aDD'
    % srepl "\bar{A}" -> "\text{a}"
    % parse \bar{A}^i_j = \bar{\gamma}^{ik} \bar{A}_{kj}
    % assign -diff_type=dD 'aUD'
    % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"

    \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij}
        + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right)
        - 2 \alpha \bar{A}_{ij} \\

    % vardef -diff_type=dD 'cf', 'trK'
    % srepl "K" -> "\text{trK}"
    %% replace 'phi' with conformal factor cf = W = e^{{-2\phi}}
    % srepl "e^{-4\phi}" -> "\text{cf}^2"
    % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
    % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
    % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"

    \partial_t \phi &= \mathcal{L}_\beta \phi
            + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

    % parse \bar{A}^{ij} = \bar{\gamma}^{jk} \bar{A}^i_k
    % assign -diff_type=dD 'aUU'
    % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
    \partial_t K &= \mathcal{L}_\beta K
        + \frac{1}{3} \alpha K^2
        + \alpha \bar{A}_{ij} \bar{A}^{ij}
        - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

    % vardef -diff_type=dD 'lambdaU'
    % srepl "\bar{\Lambda}" -> "\text{lambda}"
    % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
    % parse \Delta_{ijk}  = \bar{\gamma}_{il} \Delta^l_{jk}
    % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
    % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
    \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
            + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
            &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
            + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

    % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
    X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
        + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
        - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
    \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl} \\
    % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
    \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij}
            - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
            - 2 \alpha \bar{A}_{ik} \bar{A}^k_j
            + \alpha \bar{A}_{ij} K
            + e^{-4\phi} \hat{X}_{ij} \\

    % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
    \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

    % vardef -diff_type=dD 'betU'
    % srepl "B" -> "\text{bet}"
    % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
    \partial_t \beta^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j \beta^i\right] + B^i \\

    % vardef -const 'eta'
    % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
    \partial_t B^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j B^i\right]
        + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j \vphantom{dupD} \bar{D}_j \bar{\Lambda}^i\right]\right)
        - \eta B^i \\

    % parse \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
    % srepl "\bar{D}^2" -> "\bar{D}^i \bar{D}_i", "\mathcal{<1>}" -> "<1>"
    \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi
        \bar{D}_i \phi - 8 \bar{D}^2 \phi \right) \\

    \mathcal{M}^i &= e^{-4\phi} \left(\hat{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij}\partial_j \phi - \frac{2}{3} \bar{\gamma}^{ij}
        \partial_j K + \bar{A}^{jk} \Delta^i_{jk} + \bar{A}^{ik} \Delta^j_{jk} \right) \\

    \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij} + \frac{1}{2} \left(\bar{\gamma}_{ki}
        \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right) + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
        &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
\end{align}
\]




```python
assert_equal(H, bssncon.H)
```

    Assertion Passed!



```python
assert_equal(MU, bssncon.MU)
```

    Assertion Passed!



```python
assert_equal(RbarDD, Bq.RbarDD)
```

    Assertion Passed!


<a id='step_8'></a>

## Step 8: Output Notebook to PDF [ [^](#top) ]

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-LaTeX_Interface_Example-BSSN_Cartesian.pdf](Tutorial-LaTeX_Interface_Example-BSSN_Cartesian.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-LaTeX_Interface_Example-BSSN_Cartesian")
```

    Created Tutorial-LaTeX_Interface_Example-BSSN_Cartesian.tex, and compiled
        LaTeX file to PDF file Tutorial-LaTeX_Interface_Example-
        BSSN_Cartesian.pdf

