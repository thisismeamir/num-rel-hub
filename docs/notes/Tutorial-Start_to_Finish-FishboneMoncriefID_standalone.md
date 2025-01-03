<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Standalone Fishbone-Moncrief C Code

We start with the NRPy+ expressions generated in the [Tutorial-FishboneMoncriefID](Tutorial-FishboneMoncriefID.ipynb), and output them to the C file "FishboneMoncriefID/FMstandalone.h".

Further, $\Gamma = \alpha u^0$ is given by (as shown [here](Tutorial-u0_smallb_Poynting-Cartesian.ipynb)):
$$
\Gamma = \alpha u^0 = \sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}.
$$


```python
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
from outputC import lhrh         # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import FishboneMoncriefID.FishboneMoncriefID as fmid

# Step 1: Set up the Fishbone-Moncrief initial data. This sets all the ID gridfunctions.
fmid.FishboneMoncriefID("Spherical")

gammaDD = ixp.zerorank2()

DIM = 3
for i in range(DIM):
    for j in range(DIM):
        if i<=j:
            gammaDD[i][j] = fmid.IDgammaDD[i][j]
        else:
            gammaDD[i][j] = fmid.IDgammaDD[j][i]

# gamma_{ij} v^i_{(n)} v^j_{(n)}
Gammacontraction = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        Gammacontraction += gammaDD[i][j] * fmid.IDValencia3velocityU[i] * fmid.IDValencia3velocityU[j]

Gammafactor = sp.sqrt(1 / (1 - Gammacontraction))

# -={ F-M quantities: Generate C code from expressions and output to file }=-
FishboneMoncrief_to_print = [\
                     lhrh(lhs="alpha",rhs=fmid.IDalpha),\
                     lhrh(lhs="betaU0",rhs=fmid.IDbetaU[0]),\
                     lhrh(lhs="betaU1",rhs=fmid.IDbetaU[1]),\
                     lhrh(lhs="betaU2",rhs=fmid.IDbetaU[2]),\
                     lhrh(lhs="Gammafactor",rhs=Gammafactor),\
                     lhrh(lhs="Gamma_times_ValenciavU0",rhs=Gammafactor*fmid.IDValencia3velocityU[0]),\
                     lhrh(lhs="Gamma_times_ValenciavU1",rhs=Gammafactor*fmid.IDValencia3velocityU[1]),\
                     lhrh(lhs="Gamma_times_ValenciavU2",rhs=Gammafactor*fmid.IDValencia3velocityU[2]),\
                     lhrh(lhs="uKS4U1",rhs=fmid.uKS4U[1]),\
                     lhrh(lhs="uKS4U2",rhs=fmid.uKS4U[2]),\
                     lhrh(lhs="uKS4U3",rhs=fmid.uKS4U[3]),\
                     lhrh(lhs="uBL4U1",rhs=fmid.uBL4U[1]),\
                     lhrh(lhs="uBL4U2",rhs=fmid.uBL4U[2]),\
                     lhrh(lhs="uBL4U3",rhs=fmid.uBL4U[3])
                     ]
# print(fmid.uKS4U[3])
fin.FD_outputC("FishboneMoncriefID/FM_standalone.h",FishboneMoncrief_to_print,params="outCverbose=False,CSE_enable=False")
```

    Wrote to file "FishboneMoncriefID/FM_standalone.h"



```python
%%writefile FishboneMoncriefID/FM_standalone.c

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

const double a = 0.9375;
const double M = 1.0;
const double r_at_max_density = 12.0;
const double r_in = 6.0;

int main(int argc, const char *argv[]) {

    // Step 0a: Read command-line input, error out if nonconformant
    double xx0,xx1,xx2;
/*
    if(argc != 4) {
        printf("Error: Expected three command-line arguments: ./FM_standalone r theta phi\n");
        exit(1);
    }
    xx0 = strtod(argv[1],NULL);
    xx1 = strtod(argv[2],NULL);
    xx2 = strtod(argv[3],NULL);
*/

//    printf("# Output: r,th,ph, alpha, betaU0, betaU1, betaU2, Gamma, Gamma*vValenciaU0, Gamma*vValenciaU1, Gamma*vValenciaU2\n");
    for(double xx0=1.6;xx0<50.0;xx0+=0.2) {
        xx1 = 1.56463634120e0; //M_PI/2.0;
        xx2 = 0.0;
        double alpha,betaU0,betaU1,betaU2,Gammafactor,Gamma_times_ValenciavU0,Gamma_times_ValenciavU1,Gamma_times_ValenciavU2;
        double uKS4U1,uKS4U2,uKS4U3,uBL4U1,uBL4U2,uBL4U3;
#include "FM_standalone.h"
        if(xx0 < r_in) {
            Gammafactor = 1.0;
            Gamma_times_ValenciavU0 = Gamma_times_ValenciavU1 = Gamma_times_ValenciavU2 = 0.0;
            uKS4U1 = uKS4U2 = uKS4U3 = 0.0;
            uBL4U1 = uBL4U2 = uBL4U3 = 0.0;
        }
        printf("%e %e %e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
               xx0,xx1,xx2,
               alpha,betaU0,betaU1,betaU2,
               Gammafactor,
               Gamma_times_ValenciavU0, // util1(1) in FMtorus.f90; util(1,i,j,k) near the write statement
               Gamma_times_ValenciavU1, // util1(3) in FMtorus.f90.
               Gamma_times_ValenciavU2, // util1(2) in FMtorus.f90.
               uKS4U1,uKS4U2,uKS4U3,
               uBL4U1,uBL4U2,uBL4U3);
    }
    return 0;
}
```

    Writing FishboneMoncriefID/FM_standalone.c



```python
!gcc -O2 FishboneMoncriefID/FM_standalone.c -o FM_standalone -lm
```


```python
!./FM_standalone > out.txt
```


```python
%matplotlib inline

import matplotlib.pyplot as plt
import mpmath as mp
import csv

# Download torus_cuts.csv:
URL     = "https://etienneresearch.com/torus_cuts.csv"
outfile = "torus_cuts.csv"
try:
    with open(outfile,"w") as file:
        file.write(urllib.request.urlopen(URL).read().decode("utf-8"))
except:
    try:
        with open(outfile,"w") as file:
            file.write(urllib.urlopen(URL).read().decode("utf-8"))
    except:
        # If all else fails, hope wget does the job
        !wget -O $outfile $URL

def file_reader(filename,list_of_cols,delim=" "):
    with open(filename) as file:
        reader = csv.reader(file, delimiter=delim)
        data  = list(zip(*reader))
#         print(data)
        # data is a tuple of strings. Tuples are immutable, and we need to perform math on
        #   the data, so here we convert tuple to lists of floats:
#        data_output = [[sp.sympify(0) for i in range(len(list_of_cols))] for j in range(len(data[0]))]
        data_output = [[sp.sympify(0) for i in range(len(data[0]))] for j in range(len(list_of_cols))]
        for i in range(len(data[0])):
            for j in range(len(list_of_cols)):
#                 print(i,j,data[list_of_cols[j]][i])
                data_output[j][i] = float(data[list_of_cols[j]][i])
        return data_output

NRPy_data_output = file_reader('out.txt',       [0,7,8,9,10])
std_data_output  = file_reader('torus_cuts.csv',[0,4,1,3,2])


ylabels = ['Lorentz Gamma_{KS}=G','G*v^r_{KS,Val.}','G*v^{\\theta}_{KS,Val.}','G*v^{\phi}_{KS,Val.}']

for i in range(len(ylabels)):
    # https://matplotlib.org/gallery/text_labels_and_annotations/legend.html#sphx-glr-gallery-text-labels-and-annotations-legend-py
    fig, ax = plt.subplots()
    plt.title("NRPy's FM solve with FMtorus.f90: "+ylabels[i])
    plt.xlabel("r/M")
    plt.ylabel(ylabels[i])
    ax.plot(NRPy_data_output[0], NRPy_data_output[i+1], 'k--', label='NRPyFMSolve')
    ax.plot(std_data_output[0],   std_data_output[i+1], 'k-',  label='FMtorus.f90')
    legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
    legend.get_frame().set_facecolor('C1')
    plt.show()
```

    --2020-12-16 21:14:25--  http://astro.phys.wvu.edu/zetienne/torus_cuts.csv
    Resolving astro.phys.wvu.edu (astro.phys.wvu.edu)... 157.182.3.45
    Connecting to astro.phys.wvu.edu (astro.phys.wvu.edu)|157.182.3.45|:80... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 16445 (16K) [text/csv]
    Saving to: ‘torus_cuts.csv’
    
    torus_cuts.csv      100%[===================>]  16.06K  --.-KB/s    in 0.01s   
    
    2020-12-16 21:14:25 (1.40 MB/s) - ‘torus_cuts.csv’ saved [16445/16445]
    



    
![png](output_5_1.png)
    



    
![png](output_5_2.png)
    



    
![png](output_5_3.png)
    



    
![png](output_5_4.png)
    

