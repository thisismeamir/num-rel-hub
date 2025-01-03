<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# NRPy+ SymPy LaTeX Interface (NRPyLaTeX)

## Author: Ken Sible

### Formatting Updates by Gabriel M Steward

## The following notebook demonstrates the conversion of LaTeX to SymPy, including support for tensor operations and [Einstein notation](https://en.wikipedia.org/wiki/Einstein_notation).

<a id='top'></a>

# Table of Contents
$$\label{toc}$$

- [Step 1](#step_1): Lexical Analysis and Syntax Analysis
- [Step 2](#step_2): Grammar Demonstration and Sandbox
- [Step 3](#step_3): Tensor Support with Einstein Notation
    - [Example 1](#example_1): Tensor Contraction
    - [Example 2](#example_2): Index Raising
    - [Example 3](#example_3): Cross Product
    - [Example 4](#example_4): Covariant Derivative
    - [Example 5 (1)](#example_5_1): Schwarzschild Metric
    - [Example 5 (2)](#example_5_2): Kretschmann Scalar
    - [Example 6 (1)](#example_6_1): Extrinsic Curvature (ADM Formalism)
    - [Example 6 (2)](#example_6_2): Hamiltonian/Momentum Constraint
- [Step 4](#step_4): Exception Handling and Index Checking
- [Step 5](#step_5): Output Notebook to PDF

Further Reading: [Parsing BSSN (Cartesian) Notebook](Tutorial-LaTeX_Interface_Example-BSSN_Cartesian.ipynb)

<a id='step_1'></a>
## Step 1: Lexical Analysis and Syntax Analysis [ [^](#top) ]

In the following section, we discuss [lexical analysis](https://en.wikipedia.org/wiki/Lexical_analysis) (lexing) and [syntax analysis](https://en.wikipedia.org/wiki/Parsing) (parsing). In lexical analysis, a lexical analyzer (or scanner) can tokenize a character string, called a sentence, using substring pattern matching. In syntax analysis, a syntax analyzer (or parser) can construct a parse tree, containing all syntactic information of the language (specified by a [formal grammar](https://en.wikipedia.org/wiki/Formal_grammar)), after receiving a token iterator from the lexical analyzer.

For LaTeX to SymPy conversion, we implemented a [recursive descent parser](https://en.wikipedia.org/wiki/Recursive_descent_parser) that can construct a parse tree in [preorder traversal](https://en.wikipedia.org/wiki/Tree_traversal#Pre-order_(NLR)), starting from the root [nonterminal](https://en.wikipedia.org/wiki/Terminal_and_nonterminal_symbols), using a [right recursive](https://en.wikipedia.org/wiki/Left_recursion) grammar (partially shown below in the canonical (extended) [BNF](https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form) notation).

```
<EXPRESSION>    -> <TERM> { ( '+' | '-' ) <TERM> }*
<TERM>          -> <FACTOR> { [ '/' ] <FACTOR> }*
<FACTOR>        -> <BASE> { '^' <EXPONENT> }*
<BASE>          -> [ '-' ] ( <ATOM> | <SUBEXPR> )
<EXPONENT>      -> <BASE> | '{' <BASE> '}' | '{' '{' <BASE> '}' '}'
<ATOM>          -> <COMMAND> | <OPERATOR> | <NUMBER> | <TENSOR>
<SUBEXPR>       -> '(' <EXPRESSION> ')' | '[' <EXPRESSION> ']' | '\' '{' <EXPRESSION> '\' '}'
<COMMAND>       -> <FUNC> | <FRAC> | <SQRT> | <NLOG> | <TRIG>
    ⋮            ⋮
```

<small>**Source**: Robert W. Sebesta. Concepts of Programming Languages. Pearson Education Limited, 2016.</small>


```python
import sympy as sp
!pip install nrpylatex~=1.2 > /dev/null
!pip freeze | grep nrpylatex
from nrpylatex import *
```

    nrpylatex==1.2.3



```python
scanner = Scanner(); scanner.initialize(r'(1 + x/n)^n')
print(', '.join(token for token in scanner.tokenize()))
```

    LPAREN, INTEGER, PLUS, LETTER, DIVIDE, LETTER, RPAREN, CARET, LETTER



```python
expr = parse_latex(r'(1 + x/n)^n')
print(expr, '\n  >>', sp.srepr(expr))
```

    (1 + x/n)**n 
      >> Pow(Add(Integer(1), Mul(Pow(Symbol('n', real=True), Integer(-1)), Symbol('x', real=True))), Symbol('n', real=True))


`Grammar Derivation: (1 + x/n)^n`
```
<EXPRESSION> -> <TERM>
             -> <FACTOR>
             -> <BASE>^<EXPONENT>
             -> <SUBEXPR>^<EXPONENT>
             -> (<EXPRESSION>)^<EXPONENT>
             -> (<TERM> + <TERM>)^<EXPONENT>
             -> (<FACTOR> + <TERM>)^<EXPONENT>
             -> (<BASE> + <TERM>)^<EXPONENT>
             -> (<ATOM> + <TERM>)^<EXPONENT>
             -> (<NUMBER> + <TERM>)^<EXPONENT>
             -> (<INTEGER> + <TERM>)^<EXPONENT>
             -> (1 + <TERM>)^<EXPONENT>
             -> (1 + <FACTOR> / <FACTOR>)^<EXPONENT>
             -> ...
```

<a id='step_2'></a>
## Step 2: Grammar Demonstration and Sandbox [ [^](#top) ]

In the following section, we demonstrate the process for extending the parsing module to include a (previously) unsupported LaTeX command.

1. Update the `grammar` dictionary in the `Scanner` class with the mapping `regex` $\mapsto$ `token`.
1. Write a grammar abstraction in BNF notation (similar to a regular expression) for the command.
1. Implement a private method for the nonterminal (command name) to parse the grammar abstraction.

```<SQRT> -> <SQRT_CMD> [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'```
```
def _sqrt(self):
    self.expect('SQRT_CMD')
    if self.accept('LBRACK'):
        integer = self.scanner.lexeme
        self.expect('INTEGER')
        root = Rational(1, integer)
        self.expect('RBRACK')
    else: root = Rational(1, 2)
    self.expect('LBRACE')
    expr = self._expression()
    self.expect('RBRACE')
    if root == Rational(1, 2):
        return sqrt(expr)
    return Pow(expr, root)
```

In addition to expression parsing, we included support for equation parsing, which can produce a dictionary mapping `LHS` $\mapsto$ `RHS`, where `LHS` must be a symbol, and insert that mapping into the global namespace of the previous stack frame, as demonstrated below.

$$ \mathit{s_n} = \left(1 + \frac{1}{n}\right)^n $$


```python
parse_latex(r'\text{s_n} = \left(1 + \frac{1}{n}\right)^n')
```




    ('s_n',)




```python
print('s_n =', s_n)
```

    s_n = (1 + 1/n)**n


Furthermore, we implemented robust error messaging using the custom `ParseError` exception, which should handle every conceivable case to identify, as detailed as possible, invalid syntax inside of a LaTeX sentence. The following are some runnable examples of possible error messages.


```python
try: parse_latex(r'5x^{{4$}}')
except ScanError as e:
    print(type(e).__name__ + ': ' + str(e))
```

    ScanError: 5x^{{4$}}
                      ^
    unexpected '$' at position 6



```python
try: parse_latex(r'\sqrt[0.1]{5x^{{4}}}')
except ParseError as e:
    print(type(e).__name__ + ': ' + str(e))
```

    ParseError: \sqrt[0.1]{5x^{{4}}}
                      ^
    expected token INTEGER at position 6



```python
try: parse_latex(r'\int_0^5 5x^{{4}}dx')
except ParseError as e:
    print(type(e).__name__ + ': ' + str(e))
```

    ParseError: \int_0^5 5x^{{4}}dx
                ^
    unsupported command '\int' at position 0


In the sandbox code cell below, you can experiment with converting LaTeX to SymPy using the wrapper function `parse(sentence)`, where `sentence` must be a Python [raw string](https://docs.python.org/3/reference/lexical_analysis.html) to interpret a backslash as a literal character rather than an [escape sequence](https://en.wikipedia.org/wiki/Escape_sequence). You could, alternatively, use the supported cell magic `%%parse_latex` to automatically escape every backslash and parse the cell (more convenient than `parse(sentence)` in a notebook format).


```python
# Write Sandbox Code Here
```

<a id='step_3'></a>
## Step 3: Tensor Support with Einstein Notation [ [^](#top) ]

In the following section, we demonstrate parsing tensor notation using the Einstein summation convention. In each example, every tensor should appear either on the LHS of an equation or on the RHS of a `vardef` macro before appearing on the RHS of an equation. Furthermore, an exception will be raised upon violation of the Einstein summation convention, i.e. the occurrence of an invalid free or bound index.

**Configuration Grammar**

```
<MACRO>  -> <PARSE> | <SREPL> | <VARDEF> | <ATTRIB> | <ASSIGN> | <IGNORE>
<PARSE>  -> <PARSE_MACRO> <ASSIGNMENT> { ',' <ASSIGNMENT> }* '\\'
<SREPL>  -> <SREPL_MACRO> [ '-' <PERSIST> ] <STRING> <ARROW> <STRING> { ',' <STRING> <ARROW> <STRING> }*
<VARDEF> -> <VARDEF_MACRO> { '-' ( <OPTION> | <ZERO> ) }* <VARIABLE> [ '::' <DIMENSION> ]
            { ',' <VARIABLE> [ '::' <DIMENSION> ] }*
<ATTRIB> -> <ATTRIB_MACRO> ( <COORD_KWRD> ( <COORD> | <DEFAULT> ) | <INDEX_KWRD> ( <INDEX> | <DEFAULT> ) )
<ASSIGN> -> <ASSIGN_MACRO> { '-' <OPTION> }* <VARIABLE> { ',' <VARIABLE> }*
<IGNORE> -> <IGNORE_MACRO> <STRING> { ',' <STRING> }*
<OPTION> -> <CONSTANT> | <KRONECKER> | <METRIC> [ '=' <VARIABLE> ] | <WEIGHT> '=' <NUMBER>
            | <DIFF_TYPE> '=' <DIFF_OPT> | <SYMMETRY> '=' <SYM_OPT>
<COORD>  -> <COORD_KWRD> <LBRACK> <SYMBOL> [ ',' <SYMBOL> ]* <RBRACK>
<INDEX>  -> ( <LETTER> | '[' <LETTER> '-' <LETTER> ']' ) '::' <DIMENSION>
```

<a id='example_1'></a>
### Example 1 [Tensor Contraction](https://en.wikipedia.org/wiki/Tensor_contraction) [ [^](#top) ]


```python
parse_latex(r"""
    % define hUD --dim 4
    h = h^\mu{}_\mu
""", reset=True)
```




    ('hUD', 'h')




```python
print('h =', h)
```

    h = hUD00 + hUD11 + hUD22 + hUD33


<a id='example_2'></a>
### Example 2 [Index Raising](https://en.wikipedia.org/wiki/Raising_and_lowering_indices) [ [^](#top) ]


```python
parse_latex(r"""
    % define gUU --dim 3 --metric
    % define vD --dim 3
    v^a = g^{ab} v_b
""", reset=True)
```




    ('vU', 'gUU', 'gDD', 'GammaUDD', 'epsilonDDD', 'vD', 'gdet')




```python
print('vU =', vU)
```

    vU = [gUU00*vD0 + gUU01*vD1 + gUU02*vD2, gUU01*vD0 + gUU11*vD1 + gUU12*vD2, gUU02*vD0 + gUU12*vD1 + gUU22*vD2]


<a id='example_3'></a>
### Example 3 [Cross Product](https://en.wikipedia.org/wiki/Cross_product) [ [^](#top) ]


```python
parse_latex(r"""
    % define vU wU --dim 3
    u_i = \epsilon_{ijk} v^j w^k
""", reset=True)
```




    ('uD', 'vU', 'epsilonDDD', 'wU')




```python
print('uD =', uD)
```

    uD = [vU1*wU2 - vU2*wU1, -vU0*wU2 + vU2*wU0, vU0*wU1 - vU1*wU0]


<a id='example_4'></a>
### Example 4 [Covariant Derivative](https://en.wikipedia.org/wiki/Covariant_derivative) [ [^](#top) ]

The following are contextually inferred, dynamically generated, and injected into the global namespace for expansion of the covariant derivative $\nabla_\nu F^{\mu\nu}$
$$
\begin{align*}
    \Gamma^\mu_{ba} &= \frac{1}{2} g^{\mu c}(\partial_b\,g_{a c} + \partial_a\,g_{c b} - \partial_c\,g_{b a}) \\
    \Gamma^\nu_{ba} &= \frac{1}{2} g^{\nu c}(\partial_b\,g_{a c} + \partial_a\,g_{c b} - \partial_c\,g_{b a}) \\
    \nabla_a F^{\mu \nu} &= \partial_a F^{\mu \nu} + \Gamma^\mu_{b a} F^{b \nu} + \Gamma^\nu_{b a} F^{\mu b}
\end{align*}
$$


```python
parse_latex(r"""
    % define FUU --dim 4 --deriv dD --sym anti01
    % define gDD --dim 4 --deriv dD --metric
    % define k --const
    J^\mu = (4\pi k)^{-1} \nabla_\nu F^{\mu\nu}
""", reset=True)
```




    ('gUU',
     'FUU',
     'gDD',
     'JU',
     'epsilonUUUU',
     'GammaUDD',
     'FUU_dD',
     'k',
     'gDD_dD',
     'FUU_cdD',
     'gdet')




```python
parse_latex(r"""
    % define FUU --dim 4 --deriv dD --sym anti01
    % define ghatDD --dim 4 --deriv dD --metric
    % define k --const
    J^\mu = (4\pi k)^{-1} \hat{\nabla}_\nu F^{\mu\nu}
""", reset=True)
```




    ('ghatDD',
     'JU',
     'FUU',
     'ghatDD_dD',
     'epsilonUUUU',
     'FUU_dD',
     'k',
     'FUU_cdhatD',
     'ghatUU',
     'ghatdet',
     'GammahatUDD')



<a id='example_5_1'></a>
### Example 5 (1) [Schwarzschild Metric](https://en.wikipedia.org/wiki/Schwarzschild_metric) [ [^](#top) ]


```python
%load_ext nrpylatex.extension
```


```python
%%parse_latex --reset --ignore-warning

% coord [t, r, \theta, \phi]
% define gDD --dim 4 --zero
% define G M --const

% ignore "\begin{align}" "\end{align}"

\begin{align}
    g_{t t} &= -\left(1 - \frac{2GM}{r}\right) \\
    g_{r r} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
    g_{\theta \theta} &= r^2 \\
    g_{\phi \phi} &= r^2 \sin^2\theta
\end{align}

% assign gDD --metric
```




\[
% coord [t, r, \theta, \phi]
% define gDD --dim 4 --zero
% define G M --const

% ignore "\begin{align}" "\end{align}"

\begin{align}
    g_{t t} &= -\left(1 - \frac{2GM}{r}\right) \\
    g_{r r} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
    g_{\theta \theta} &= r^2 \\
    g_{\phi \phi} &= r^2 \sin^2\theta
\end{align}

% assign gDD --metric
\]




```python
sp.Matrix(gDD)
```




$\displaystyle \left[\begin{matrix}\frac{2 G M}{r} - 1 & 0 & 0 & 0\\0 & \frac{1}{- \frac{2 G M}{r} + 1} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]$



<a id='example_5_2'></a>
### Example 5 (2) [Kretschmann Scalar](https://en.wikipedia.org/wiki/Kretschmann_scalar) [ [^](#top) ]


```python
%%parse_latex

% ignore "\begin{align}" "\end{align}"

\begin{align}
    R^\alpha{}_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu}
        + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
    K &= R^{\alpha\beta\mu\nu} R_{\alpha\beta\mu\nu} \\
    R_{\beta\nu} &= R^\alpha{}_{\beta\alpha\nu} \\
    R &= g^{\beta\nu} R_{\beta\nu} \\
    G_{\beta\nu} &= R_{\beta\nu} - \frac{1}{2}g_{\beta\nu}R
\end{align}

```




\[
% ignore "\begin{align}" "\end{align}"

\begin{align}
    R^\alpha{}_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu}
        + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
    K &= R^{\alpha\beta\mu\nu} R_{\alpha\beta\mu\nu} \\
    R_{\beta\nu} &= R^\alpha{}_{\beta\alpha\nu} \\
    R &= g^{\beta\nu} R_{\beta\nu} \\
    G_{\beta\nu} &= R_{\beta\nu} - \frac{1}{2}g_{\beta\nu}R
\end{align}
\]




```python
sp.simplify(sp.Matrix(RDD))
```




$\displaystyle \left[\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right]$




```python
display(sp.Matrix(GammaUDD[0][:][:]))
```


$\displaystyle \left[\begin{matrix}0 & - \frac{G M}{r^{2} \cdot \left(\frac{2 G M}{r} - 1\right)} & 0 & 0\\- \frac{G M}{r^{2} \cdot \left(\frac{2 G M}{r} - 1\right)} & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right]$



```python
display(sp.Matrix(GammaUDD[1][:][:]))
```


$\displaystyle \left[\begin{matrix}\frac{G M \left(- \frac{2 G M}{r} + 1\right)}{r^{2}} & 0 & 0 & 0\\0 & - \frac{G M}{r^{2} \left(- \frac{2 G M}{r} + 1\right)} & 0 & 0\\0 & 0 & - r \left(- \frac{2 G M}{r} + 1\right) & 0\\0 & 0 & 0 & - r \left(- \frac{2 G M}{r} + 1\right) \sin^{2}{\left(\theta \right)}\end{matrix}\right]$



```python
display(sp.Matrix(GammaUDD[2][:][:]))
```


$\displaystyle \left[\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & \frac{1}{r} & 0\\0 & \frac{1}{r} & 0 & 0\\0 & 0 & 0 & - \sin{\left(\theta \right)} \cos{\left(\theta \right)}\end{matrix}\right]$



```python
display(sp.Matrix(GammaUDD[3][:][:]))
```


$\displaystyle \left[\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{r}\\0 & 0 & 0 & \frac{\cos{\left(\theta \right)}}{\sin{\left(\theta \right)}}\\0 & \frac{1}{r} & \frac{\cos{\left(\theta \right)}}{\sin{\left(\theta \right)}} & 0\end{matrix}\right]$


For the Schwarzschild metric, the Kretschmann scalar $K$ has the property that $K\to\infty$ as $r\to 0$, and hence the metric and spacetime itself are undefined at the point of infinite curvature $r=0$, indicating the presence of a physical singularity since the Kretschmann scalar is an [invariant quantity](https://en.wikipedia.org/wiki/Curvature_invariant_(general_relativity)) in general relativity.


```python
display(sp.simplify(K))
```


$\displaystyle \frac{48 G^{2} M^{2}}{r^{6}}$


In a [vacuum region](https://en.wikipedia.org/wiki/Vacuum_solution_(general_relativity)#:~:text=In%20general%20relativity%2C%20a%20vacuum,non%2Dgravitational%20fields%20are%20present.), such as the spacetime described by the Schwarzschild metric, $T_{\mu\nu}=0$ and hence $G_{\mu\nu}=0$ since $G_{\mu\nu}=8\pi G\,T_{\mu\nu}$ ([Einstein Equations](https://en.wikipedia.org/wiki/Einstein_field_equations)).


```python
sp.simplify(sp.Matrix(GDD))
```




$\displaystyle \left[\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right]$



<a id='example_6_1'></a>
### Example 6 (1) [Extrinsic Curvature](https://en.wikipedia.org/wiki/Curvature) ([ADM Formalism](https://en.wikipedia.org/wiki/ADM_formalism)) [ [^](#top) ]


```python
%%parse_latex --ignore-warning

% coord [r, \theta, \phi]
% ignore "\begin{align}" "\end{align}"
\begin{align}
    \gamma_{ij} &= g_{ij} \\
    % assign gammaDD --metric
    \beta_i &= g_{r i} \\
    \alpha &= \sqrt{\gamma^{ij}\beta_i\beta_j - g_{r r}} \\
    K_{ij} &= \frac{1}{2\alpha}\left(\nabla_i \beta_j + \nabla_j \beta_i\right) \\
    K &= \gamma^{ij} K_{ij}
\end{align}

```




\[
% coord [r, \theta, \phi]
% ignore "\begin{align}" "\end{align}"
\begin{align}
    \gamma_{ij} &= g_{ij} \\
    % assign gammaDD --metric
    \beta_i &= g_{r i} \\
    \alpha &= \sqrt{\gamma^{ij}\beta_i\beta_j - g_{r r}} \\
    K_{ij} &= \frac{1}{2\alpha}\left(\nabla_i \beta_j + \nabla_j \beta_i\right) \\
    K &= \gamma^{ij} K_{ij}
\end{align}
\]



For the Schwarzschild metric (defined in the previous example), the extrinsic curvature in the ADM formalism should evaluate to zero.


```python
display(sp.Matrix(KDD))
```


$\displaystyle \left[\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right]$


<a id='example_6_2'></a>
### Example 6 (2) [Hamiltonian/Momentum Constraint](https://en.wikipedia.org/wiki/Hamiltonian_constraint) [ [^](#top) ]


```python
%%parse_latex --ignore-warning

% ignore "\begin{align}" "\end{align}"

\begin{align}
    R_{ij} &= \partial_k \Gamma^k_{ij} - \partial_j \Gamma^k_{ik}
        + \Gamma^k_{ij}\Gamma^l_{kl} - \Gamma^l_{ik}\Gamma^k_{lj} \\
    R &= \gamma^{ij} R_{ij} \\
    E &= \frac{1}{16\pi}\left(R + K^{{2}} - K_{ij}K^{ij}\right) \\
    p_i &= \frac{1}{8\pi}\left(D_j \gamma^{jk} K_{ki} - D_i K\right)
\end{align}
```




\[
% ignore "\begin{align}" "\end{align}"

\begin{align}
    R_{ij} &= \partial_k \Gamma^k_{ij} - \partial_j \Gamma^k_{ik}
        + \Gamma^k_{ij}\Gamma^l_{kl} - \Gamma^l_{ik}\Gamma^k_{lj} \\
    R &= \gamma^{ij} R_{ij} \\
    E &= \frac{1}{16\pi}\left(R + K^{{2}} - K_{ij}K^{ij}\right) \\
    p_i &= \frac{1}{8\pi}\left(D_j \gamma^{jk} K_{ki} - D_i K\right)
\end{align}
\]



Every solution to the Einstein Equations, including Schwarzschild, must satisfy the Hamiltonian constraint ($E=0$) and the Momentum constraint ($p_i=0$).


```python
print('E = %s, pD = %s' % (sp.simplify(E), pD))
```

    E = 0, pD = [0, 0, 0]


<a id='step_4'></a>
## Step 4: Exception Handling and Index Checking ( [^](#top) )

We extended our robust error messaging using the custom `TensorError` exception, which should handle any inconsistent tensor dimension and any violation of the Einstein summation convention, specifically that a bound index must appear exactly once as a superscript and exactly once as a subscript in any single term and that a free index must appear in every term with the same position and cannot be summed over in any term.


```python
%%parse_latex --reset

% define TUD uD --dim 4
v^\mu = T^\mu_\nu u_\nu
```

    TensorError: illegal bound index 'nu' in vU



```python
%%parse_latex --reset

% define TUD uD --dim 4
v^\mu = T^\mu_\nu u_\mu
```

    TensorError: unbalanced free indices {'nu', 'mu'} in vU



```python
%%parse_latex --reset

% define TUD --dim 4
% define uD --dim 3

v_\nu = T^\mu_\nu u_\mu
```

    ParseError: index out of range; change loop/summation range



```python
%%parse_latex --reset

% define vD --dim 4
T_{\mu\nu} = v_\mu w_\nu
```

    ParseError: T_{\mu\nu} = v_\mu w_\nu
                                   ^
    cannot index undefined variable 'wD' at position 39



```python
%%parse_latex --reset

% define FUU --dim 4 --sym anti01
% define k --const
J^\mu = (4\pi k)^{-1} \nabla_\nu F^{\mu\nu}
```

    ParseError: J^\mu = (4\pi k)^{-1} \nabla_\nu F^{\mu\nu}
                                      ^
    cannot generate covariant derivative without defined metric 'g'


<a id='step_5'></a>
## Step 5: Output Notebook to PDF ( [^](#top) )

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-SymPy_LaTeX_Interface.pdf](Tutorial-SymPy_LaTeX_Interface.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means).


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-SymPy_LaTeX_Interface")
```

    Created Tutorial-SymPy_LaTeX_Interface.tex, and compiled LaTeX file to PDF
        file Tutorial-SymPy_LaTeX_Interface.pdf

