<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Creating the Finite-Difference Table for Centered First and Second  Derivatives, from 2nd through 10th-Order Accuracy

<!-- abstract? -->

## *Courtesy Brandon Clark*


```python
print("Installing astropy, needed for creating the output table. Please wait a few seconds...")
!pip install -U pip astropy > /dev/null
print("astropy installed.")
# Step 0: Import needed modules
import numpy as np
import finite_difference as fin
from astropy.table import Table

# Step 1: Set the maximum finite-difference accuracy order computed in the table
max_fdorder = 10

# Step 2: Set up table parameters
#    One column for deriv order, one for deriv accuracy, and max_fdorder+1
numcols = 2 + max_fdorder + 1
#    8 rows: max_fdorder accuracy orders per derivative order, times 2 derivative orders (first & second derivative)
numrows = int(max_fdorder/2 * 2)
#    Center column index of table will be at 2 + max_fdorder/2  (zero-offset indexing)
column_corresponding_to_zero_fd_point = 2 + int(max_fdorder/2)
#    The table is initialized as a matrix of zeroes in numpy...
numpy_matrix = np.zeros((numrows, numcols), dtype=object)
#    Then we replace all elements with the empty string to match the Wikipedia article.
for row in range(numrows):
    for col in range(numcols):
        numpy_matrix[row,col] = ""

# Step 3: Construct the first-order derivative finite difference coefficients
rowcount = 0
for fdorder in range(2, max_fdorder+1, 2): # loop runs from 2 to max_fdorder inclusive, skipping odd orders.
    numpy_matrix[rowcount, 0] = "1st"
    numpy_matrix[rowcount, 1] = fdorder
    fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl("D0", fdorder)
    for i in range(fdorder):
        numpy_matrix[rowcount, column_corresponding_to_zero_fd_point + fdstencl[i][0]] = fdcoeffs[i]
    rowcount += 1

# Step 4: Construct the second-order derivative finite difference coefficients
for fdorder in range(2, max_fdorder+1, 2): # loop runs from 2 to max_fdorder inclusive, skipping odd orders.
    numpy_matrix[rowcount, 0] = "2nd"
    numpy_matrix[rowcount, 1] = fdorder
    fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl("DD00", fdorder)
    for i in range(fdorder+1):
        numpy_matrix[rowcount, column_corresponding_to_zero_fd_point + fdstencl[i][0]] = fdcoeffs[i]
    rowcount += 1

# Step 5: Construct an astropy table from the numpy matrix with the following header info, and then print it:
colnames = ['Derivative','Accuracy']
for i in range(-int(max_fdorder/2),int(max_fdorder/2)+1):
    colnames.append(str(i))
table = Table(numpy_matrix, names=colnames)
table.pprint(max_width=-1)
```

    Installing astropy, needed for creating the output table. Please wait a few seconds...
    [33mWARNING: Retrying (Retry(total=4, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2436baa30>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/pip/[0m
    [33mWARNING: Retrying (Retry(total=3, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2436bac40>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/pip/[0m
    [33mWARNING: Retrying (Retry(total=2, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2436b58b0>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/pip/[0m
    [33mWARNING: Retrying (Retry(total=1, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2436b5dc0>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/pip/[0m
    [33mWARNING: Retrying (Retry(total=0, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2436b5cd0>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/pip/[0m
    [33mWARNING: Retrying (Retry(total=4, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2429bdd90>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/astropy/[0m
    [33mWARNING: Retrying (Retry(total=3, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2429bdf40>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/astropy/[0m
    [33mWARNING: Retrying (Retry(total=2, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2429c9130>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/astropy/[0m
    [33mWARNING: Retrying (Retry(total=1, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2429c92e0>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/astropy/[0m
    [33mWARNING: Retrying (Retry(total=0, connect=None, read=None, redirect=None, status=None)) after connection broken by 'NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7ff2429c9490>: Failed to establish a new connection: [Errno -2] Name or service not known')': /simple/astropy/[0m
    astropy installed.
    Derivative Accuracy    -5      -4     -3     -2   -1      0       1    2     3      4      5   
    ---------- -------- ------- ------- ------ ----- ---- ---------- --- ----- ----- ------- ------
           1st        2                              -1/2            1/2                           
           1st        4                         1/12 -2/3            2/3 -1/12                     
           1st        6                  -1/60  3/20 -3/4            3/4 -3/20  1/60               
           1st        8           1/280 -4/105   1/5 -4/5            4/5  -1/5 4/105  -1/280       
           1st       10 -1/1260   5/504  -5/84  5/21 -5/6            5/6 -5/21  5/84  -5/504 1/1260
           2nd        2                                 1         -2   1                           
           2nd        4                        -1/12  4/3       -5/2 4/3 -1/12                     
           2nd        6                   1/90 -3/20  3/2     -49/18 3/2 -3/20  1/90               
           2nd        8          -1/560  8/315  -1/5  8/5    -205/72 8/5  -1/5 8/315  -1/560       
           2nd       10  1/3150 -5/1008  5/126 -5/21  5/3 -5269/1800 5/3 -5/21 5/126 -5/1008 1/3150

