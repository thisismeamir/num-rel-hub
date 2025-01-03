<script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-59152712-8');
</script>

# Exercise (Loop Generation) Solution

## *Courtesy Ken Sible*


```python
from loop import loop # Import NRPy+ module for loop generation

boundary = 'u[n][0] = u[n][Nx] = 0;\n'
inner_1 = loop('k', '1', '(Nx - 1)', '1', '', interior='u[n + 1][k] = u[n][k] + r*(u[n][k + 1] - 2*u[n][k] + u[n][k - 1]);')
inner_2 = loop('k', '0', 'Nx', '1', '', interior='u[n][k] = u[n + 1][k];')
print(loop('n', '0', '(Nt - 1)', '1', '', interior=(boundary + inner_1 + inner_2[:-1])))
```

    for (int n = 0; n < (Nt - 1); n++) {
        u[n][0] = u[n][Nx] = 0;
        for (int k = 1; k < (Nx - 1); k++) {
            u[n + 1][k] = u[n][k] + r*(u[n][k + 1] - 2*u[n][k] + u[n][k - 1]);
        } // END LOOP: for (int k = 1; k < (Nx - 1); k++)
        for (int k = 0; k < Nx; k++) {
            u[n][k] = u[n + 1][k];
        } // END LOOP: for (int k = 0; k < Nx; k++)
    } // END LOOP: for (int n = 0; n < (Nt - 1); n++)
    

