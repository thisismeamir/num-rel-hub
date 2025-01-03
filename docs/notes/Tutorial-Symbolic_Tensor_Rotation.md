<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Symbolic Tensor (Quaternion) Rotation

## Author: Ken Sible

## The following module presents an algorithm for symbolic vector or tensor rotation using SymPy. It validates the algorithm against known results and demonstrates common subexpression elimination for optimized symbolic computation.

### NRPy+ Source Code for this module:
1. [tensor_rotation.py](../edit/tensor_rotation.py); [\[**tutorial**\]](Tutorial-Symbolic_Tensor_Rotation.ipynb) The tensor_rotation.py script will perform symbolic tensor rotation using the following function: `rotate(tensor, axis, angle)`. 

<a id='toc'></a>

# Table of Contents
$$\label{toc}$$

0. [Step 0](#prelim): Derivation of quaternion rotation from matrix rotation using [linear algebra](https://en.wikipedia.org/wiki/Linear_algebra) and [ring theory](https://en.wikipedia.org/wiki/Ring_theory)
1. [Step 1](#algorithm): Discussion of the tensor rotation algorithm using the [SymPy](https://www.sympy.org) package for symbolic manipulation </header_section>
1. [Step 2](#validation): Validation and demonstration of the tensor rotation algorithm, including common subexpression elimination
    1. [Step 2.a](#axisrotation): Verification of three-dimensional rotation about a coordinate axis using an equivalent rotation matrix
    1. [Step 2.b](#symtensor): Verification of symbolic tensor rotation and the invariance of argument permutation for a [symmetric tensor](https://en.wikipedia.org/wiki/Symmetric_tensor)
    1. [Step 2.c](#cse): Demonstration of [common subexpression elimination](https://en.wikipedia.org/wiki/Common_subexpression_elimination) applied to a rotated symbolic matrix
1. [Step 3](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

<a id='prelim'></a>

# Mathematical Background (Optional): Quaternion Rotation \[Back to [top](#toc)\]
$$\label{prelim}$$ 

Let $\vec{v}$ denote a vector in $\mathbb{R}^2$. We recall from linear algebra that $\vec{v}$ rotated through an angle $\theta$ about the x-axis, denoted $\vec{v}'$, has the following matrix formula

$$\vec{v}'=\begin{bmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{bmatrix}\vec{v}.$$

Let $z=a+bi\in\mathbb{C}$ for some $a,b\in\mathbb{R}$. Consider the corresponding (or [isomorphic](https://en.wikipedia.org/wiki/Isomorphism)) vector $\vec{z}=(a,b)\in\mathbb{R}^2$. We observe from the rotation formula that $\vec{v}'=a(\cos\theta,\sin\theta)+b(-\sin\theta,\cos\theta)$ after expanding the matrix product. Let $w=c+di\in\mathbb{C}$ for some $c,d\in\mathbb{R}$. We recall the definition of the complex product as $zw=(ac-bd)+i(ad+bc)$ where $i^2=-1$. Hence, $z'=(a+bi)(\cos\theta+i\sin\theta)=(a+bi)e^{i\theta}$ after comparing the rotated vector $\vec{v}'$ with the complex product as defined above. Therefore, the following compact formula arises for vector rotation (given the correspondence between $\mathbb{C}$ and $\mathbb{R}^2$):

$$z'=e^{i\theta}z.$$

However, for vector rotation in three-dimensional space, we curiously require a four-dimensional, non-commutative extension of the complex number system. The following discussion will provide an overview of the connection with vector rotation. For further reading, we suggest the document ([quaternion_rotation](http://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf)).

$$\mathcal{H}=\{a+bi+cj+dk:a,b,c,d\in\mathbb{R}\text{ and }i^2=j^2=k^2=ijk=-1\}$$

Consider the special case of rotating a vector $\vec{v}\in\mathbb{R}^3$ through an angle $\theta$ about a normalized rotation axis $\vec{n}$ perpendicular to the vector $\vec{v}$. We typically decompose a quaternion into a scalar and vector component whenever performing vector rotation, such as the decomposition of $q=a+bi+cj+dk$ into $q=(a,\vec{w})$ where $\vec{w}=(b,c,d)$. For future convenience, we define the following quaternions for vector rotation: $v=(0,\vec{v})$, $v'=(0,\vec{v}')$, and $n=(0,\vec{n})$.

From the fundamental formula of quaternion algebra, $i^2=j^2=k^2=ijk=-1$, we could derive quaternion multiplication, known as the [Hamilton product](https://en.wikipedia.org/wiki/Quaternion#Hamilton_product). Let $q_1,q_2\in\mathcal{H}$ where both $q_1$ and $q_2$ have zero scalar component. The Hamilton product of $q_1$ and $q_2$, after some straightforward verification, has the form $q_1q_2=(-\vec{q}_1\cdot\vec{q}_2,\vec{q}_1\times\vec{q}_2)$. Hence, $nv=(0,\vec{n}\times\vec{v})$ since $\vec{n}$ and $\vec{v}$ are orthogonal. We observe that the projection of the rotated vector $\vec{v}'$ onto $\vec{v}$ is $\cos\theta\,\vec{v}$ and the projection of $\vec{v}'$ onto $\vec{n}\times\vec{v}$ is $\sin\theta\,(\vec{n}\times\vec{v})$, and hence $\vec{v}'=\cos\theta\,\vec{v}+\sin\theta\,(\vec{n}\times\vec{v})$. We define the quaternion exponential as $e^{n\theta}=\cos\theta+n\sin\theta$, analogous to the complex exponential. Finally, we identify the formula for vector rotation in $\mathbb{R}^3$ by comparing the rotated vector $\vec{v}'$ with the Hamilton product:

$$v'=(\cos\theta+n\sin\theta)v=e^{n\theta}v$$

The tensor rotation algorithm defined by the function `rotate(tensor, axis, angle)` in `tensor_rotation.py` does general vector and tensor rotation in $\mathbb{R}^3$ about an arbitrary rotation axis, not necessarily orthogonal to the original vector or tensor. For general three-dimensional rotation, we define the rotation quaternion as $q=e^{n(\theta/2)}$ and the conjugate of $q$ as $q^*=e^{-n(\theta/2)}$. The quaternion rotation operator has the following form for general vector and tensor rotation.  

$$\mathcal{L}[v]=qvq^*,\,\,\mathcal{L}[M]=(q(qMq^*)^\text{T}q^*)^\text{T}$$

We should remark that quaternion-matrix multiplication is defined as column-wise quaternion multiplication [(source)](https://people.dsv.su.se/~miko1432/rb/Rotations%20of%20Tensors%20using%20Quaternions%20v0.3.pdf). Furthermore, we claim that $vq^*=qv$ whenever the rotation axis $\vec{n}$ and the vector $\vec{v}$ are perpendicular (straightforward verification), which will recover the rotation formula for the special case.

<a id='algorithm'></a>

# Step 1: The Tensor Rotation Algorithm \[Back to [top](#toc)\]
$$\label{algorithm}$$

### Algorithm Pseudocode
```
function rotate(tensor, axis, angle):
    initialize quaternion q from rotation axis and angle
    if tensor is a nested list (i.e. matrix)
        convert tensor to a symbolic matrix
        if not (size of symbolic matrix is (3, 3))
            throw error for invalid matrix size
        end if
        initialize empty vector M
        for column in tensor
            append quaternion(column) onto M
        end for
        M = q * M *(conjugate of q) // rotate each column
        replace each column of tensor with the vector part
            of the associated column quaternion in M
        initialize empty vector M
        for row in tensor
            append quaternion(row) onto M
        end for
        M = q * M *(conjugate of q) // rotate each row
        replace each row of tensor with the vector part
            of the associated row quaternion in M
        convert tensor to a nested list
        return tensor
    else if tensor is a list (i.e. vector)
        if not (length of vector is 3):
            throw error for invalid vector length
        v = q * quaternion(tensor) * (conjugate of q)
        replace each element of tensor with the vector part
            of the associated vector quaternion v
        return tensor
    end if
    throw error for unsupported tensor type
end function
```


```python
""" Symbolic Tensor (Quaternion) Rotation

The following script will perform symbolic tensor rotation using quaternions.
"""

from sympy import Quaternion as quat
from sympy import Matrix
from sympy.functions import transpose

# Input:  tensor = 3-vector or (3x3)-matrix
#         axis   = rotation axis (normal 3-vector)
#         angle  = rotation angle (in radians)
# Output: rotated tensor (of original type)
def rotate(tensor, axis, angle):
    # Quaternion-Matrix Multiplication
    def mul(*args):
        if isinstance(args[0], list):
            q, M = args[1], args[0]
            for i, col in enumerate(M):
                M[i] = col * q
        else:
            q, M = args[0], args[1]
            for i, col in enumerate(M):
                M[i] = q * col
        return M
    # Rotation Quaternion (Axis, Angle)
    q = quat.from_axis_angle(axis, angle)
    if isinstance(tensor[0], list):
        tensor = Matrix(tensor)
        if tensor.shape != (3, 3):
            raise Exception('Invalid Matrix Size')
        # Rotation Formula: M' = (q.(q.M.q*)^T.q*)^T
        M = [quat(0, *tensor[:, i]) for i in range(tensor.shape[1])]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(tensor.shape[1]):
            tensor[:, i] = [M[i].b, M[i].c, M[i].d]
        M = [quat(0, *tensor[i, :]) for i in range(tensor.shape[0])]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(tensor.shape[0]):
            tensor[i, :] = [[M[i].b, M[i].c, M[i].d]]
        return tensor.tolist()
    if isinstance(tensor, list):
        if len(tensor) != 3:
            raise Exception('Invalid Vector Length')
        # Rotation Formula: v' = q.v.q*
        v = q * quat(0, *tensor) * q.conjugate()
        return [v.b, v.c, v.d]
    raise Exception('Unsupported Tensor Type')
```

<a id='validation'></a>

# Step 2: Validation and Demonstration \[Back to [top](#toc)\]
$$\label{validation}$$

In the following section, we demonstrate the rotation algorithm, specifically for symbolic tensor rotation, and validate that the numeric or symbolic output from the algorithm does agree with a known result, usually obtained from a rotation matrix. The format of each code cell has the structure: generate expected result from a rotation matrix or NRPy+, generate received result from the rotation algorithm, and assert that these are both equivalent.

<a id='axisrotation'></a>

## Step 2.a: Verifying Coordinate Axis Rotation \[Back to [top](#toc)\]
$$\label{axisrotation}$$

We recall that any general three-dimensional rotation can be expressed as the composition of a rotation about each coordinate axis (see [rotation theorem](https://en.wikipedia.org/wiki/Euler%27s_rotation_theorem)). Therefore, we validate our rotation algorithm using only rotations about each coordinate axis rather than about an arbitrary axis in three-dimensional space. Consider the following vector $\vec{v}$ and matrix $M$ defined below using the SymPy package.


```python
from sympy.matrices import rot_axis1, rot_axis2, rot_axis3
from sympy import latex, pi
from IPython.display import Math

v, angle = [1, 0, 1], pi/2
M = [[1, 2, 1], [0, 1, 0], [2, 1, 2]]
Math('\\vec{v}=%s,\\,M=%s' % (latex(Matrix(v)), latex(Matrix(M))))
```




$\displaystyle \vec{v}=\left[\begin{matrix}1\\0\\1\end{matrix}\right],\,M=\left[\begin{matrix}1 & 2 & 1\\0 & 1 & 0\\2 & 1 & 2\end{matrix}\right]$



We further recall that for any rotation matrix $R$ and vector $\vec{v}$, the rotated vector $\vec{v}'$ has the formula $\vec{v}'=R\vec{v}$ and the rotated matrix $M'$ has the formula $M'=RMR^\text{-1}$, where $R^{-1}=R^\text{T}$ since every rotation matrix is an [orthogonal matrix](https://en.wikipedia.org/wiki/Orthogonal_matrix).


```python
# vector rotation about x-axis
expected = rot_axis1(-angle) * Matrix(v)
received = Matrix(rotate(v, [1, 0, 0], angle))
assert expected == received; v_ = received

# matrix rotation about x-axis
expected = rot_axis1(-angle) * Matrix(M) * transpose(rot_axis1(-angle))
received = Matrix(rotate(M, [1, 0, 0], angle))
assert expected == received; M_ = received

Math('\\vec{v}\'=%s,\\,M\'=%s' % (latex(Matrix(v_)), latex(Matrix(M_))))
```




$\displaystyle \vec{v}'=\left[\begin{matrix}1\\-1\\0\end{matrix}\right],\,M'=\left[\begin{matrix}1 & -1 & 2\\-2 & 2 & -1\\0 & 0 & 1\end{matrix}\right]$




```python
# vector rotation about y-axis
expected = rot_axis2(-angle) * Matrix(v)
received = Matrix(rotate(v, [0, 1, 0], angle))
assert expected == received; v_ = received

# matrix rotation about y-axis
expected = rot_axis2(-angle) * Matrix(M) * transpose(rot_axis2(-angle))
received = Matrix(rotate(M, [0, 1, 0], angle))
assert expected == received; M_ = received

Math('\\vec{v}\'=%s,\\,M\'=%s' % (latex(Matrix(v_)), latex(Matrix(M_))))
```




$\displaystyle \vec{v}'=\left[\begin{matrix}1\\0\\-1\end{matrix}\right],\,M'=\left[\begin{matrix}2 & 1 & -2\\0 & 1 & 0\\-1 & -2 & 1\end{matrix}\right]$




```python
# vector rotation about z-axis
expected = rot_axis3(-angle) * Matrix(v)
received = Matrix(rotate(v, [0, 0, 1], angle))
assert expected == received; v_ = received

# matrix rotation about z-axis
expected = rot_axis3(-angle) * Matrix(M) * transpose(rot_axis3(-angle))
received = Matrix(rotate(M, [0, 0, 1], angle))
assert expected == received; M_ = received

Math('\\vec{v}\'=%s,\\,M\'=%s' % (latex(Matrix(v_)), latex(Matrix(M_))))
```




$\displaystyle \vec{v}'=\left[\begin{matrix}0\\1\\1\end{matrix}\right],\,M'=\left[\begin{matrix}1 & 0 & 0\\-2 & 1 & 1\\-1 & 2 & 2\end{matrix}\right]$



<a id='symtensor'></a>

## Step 2.b: Verifying Symbolic Tensor Rotation \[Back to [top](#toc)\]
$$\label{symtensor}$$

The rotation algorithm does support symbolic rotation, as shown below with the second rank, symmetric tensor $h^{\mu\nu}$ rotated about the 4-vector $v^\mu$.


```python
from sympy import symbols, simplify, cse
import indexedexp as ixp # NRPy+: symbolic indexed expressions

angle = symbols('theta', real=True)
vU    = ixp.declarerank1("vU")
hUU   = ixp.declarerank2("hUU", "sym01")
rotatedhUU, N = rotate(hUU, vU, angle), len(hUU)

Math('hUU=%s' % latex(Matrix(hUU)))
```




$\displaystyle hUU=\left[\begin{matrix}hUU_{00} & hUU_{01} & hUU_{02}\\hUU_{01} & hUU_{11} & hUU_{12}\\hUU_{02} & hUU_{12} & hUU_{22}\end{matrix}\right]$



We demonstrate that a completely symbolic rotation applied to the tensor $h^{\mu\nu}$ does preserve the index symmetry ($h^{\mu\nu}=h^{\nu\mu}$).


```python
for i in range(N):
    for j in range(N):
        if j >= i: continue
        assert simplify(rotatedhUU[i][j] - rotatedhUU[j][i]) == 0
        print('Assertion Passed: rotatedhUU[{i}][{j}] == rotatedhUU[{j}][{i}]'.format(i=i, j=j))
```

    Assertion Passed: rotatedhUU[1][0] == rotatedhUU[0][1]
    Assertion Passed: rotatedhUU[2][0] == rotatedhUU[0][2]
    Assertion Passed: rotatedhUU[2][1] == rotatedhUU[1][2]


<a id='cse'></a>

## Step 2.c: Common Subexpression Elimination \[Back to [top](#toc)\]
$$\label{cse}$$

If the rotation algorithm is given any symbolic input, then the resulting expression will support common subexpression elimination.


```python
cse_rotatedhUU = cse(Matrix(rotatedhUU))[1][0]
Math('\\text{cse}(hUU)=%s' % latex(Matrix(cse_rotatedhUU)))
```




$\displaystyle \text{cse}(hUU)=\left[\begin{matrix}x_{1} x_{33} + x_{3} x_{34} - x_{36} x_{5} + x_{37} x_{9} & x_{1} x_{36} + x_{3} x_{37} + x_{33} x_{5} - x_{34} x_{9} & x_{1} x_{34} - x_{3} x_{33} + x_{36} x_{9} + x_{37} x_{5}\\x_{1} x_{41} + x_{3} x_{42} - x_{44} x_{5} + x_{45} x_{9} & x_{1} x_{44} + x_{3} x_{45} + x_{41} x_{5} - x_{42} x_{9} & x_{1} x_{42} - x_{3} x_{41} + x_{44} x_{9} + x_{45} x_{5}\\x_{1} x_{49} + x_{3} x_{50} - x_{5} x_{52} + x_{53} x_{9} & x_{1} x_{52} + x_{3} x_{53} + x_{49} x_{5} - x_{50} x_{9} & x_{1} x_{50} - x_{3} x_{49} + x_{5} x_{53} + x_{52} x_{9}\end{matrix}\right]$



<a id='latex_pdf_output'></a>

# Step 3: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
$$\label{latex_pdf_output}$$

The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
[Tutorial-Symbolic_Tensor_Rotation.pdf](Tutorial-Symbolic_Tensor_Rotation.pdf). (Note that clicking on this link may not work; you may need to open the PDF file through another means.)


```python
import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Symbolic_Tensor_Rotation")
```

    Created Tutorial-Symbolic_Tensor_Rotation.tex, and compiled LaTeX file to
        PDF file Tutorial-Symbolic_Tensor_Rotation.pdf

