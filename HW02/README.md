# HW 2: Galaxy deblur

d is the blured image of the galaxy via the transformation

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=d(x,&space;y)&space;=&space;\int_{x'_{min}}^{x'_{max}}\int_{y'_{min}}^{y'_{max}}&space;k(x-x',&space;y-y')m(x',&space;y')dx'dy'" target="_blank"><img src="https://latex.codecogs.com/svg.latex?d(x,&space;y)&space;=&space;\int_{x'_{min}}^{x'_{max}}\int_{y'_{min}}^{y'_{max}}&space;k(x-x',&space;y-y')m(x',&space;y')dx'dy'" title="d(x, y) = \int_{x'_{min}}^{x'_{max}}\int_{y'_{min}}^{y'_{max}} k(x-x', y-y')m(x', y')dx'dy'" /></a>
</p>

where the kernel is given by

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=k(x-x',&space;y-y')&space;=&space;\frac{1}{\sqrt{\pi}}\exp\left(-\frac{(x-x')^2&plus;(y-y')^2)}{2}\right)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?k(x-x',&space;y-y')&space;=&space;\frac{1}{\sqrt{\pi}}\exp\left(-\frac{(x-x')^2&plus;(y-y')^2)}{2}\right)" title="k(x-x', y-y') = \frac{1}{\sqrt{\pi}}\exp\left(-\frac{(x-x')^2+(y-y')^2)}{2}\right)" /></a>
</p>

Consider the discrete integration then

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=d(x,&space;y)&space;=&space;\sum_{x'_{min}}^{x'_{max}}\sum_{y'_{min}}^{y'_{max}}&space;k(x-x',&space;y-y')m(x',&space;y')" target="_blank"><img src="https://latex.codecogs.com/svg.latex?d(x,&space;y)&space;=&space;\sum_{x'_{min}}^{x'_{max}}\sum_{y'_{min}}^{y'_{max}}&space;k(x-x',&space;y-y')m(x',&space;y')" title="d(x, y) = \sum_{x'_{min}}^{x'_{max}}\sum_{y'_{min}}^{y'_{max}} k(x-x', y-y')m(x', y')" /></a>
</p>

Since governing equation is

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=d&space;=&space;K(m)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?d&space;=&space;K(m)" title="d = K(m)" /></a>
</p>

The Kernel will be a tensor rank 4 where it could leading to a complicated operation. Then we need to convert the indices from matrix to vector and the kernel will become a tensor rank 2 as

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=d(i,&space;i')&space;=&space;\sum_{i'}&space;k_{ii'}m(i')" target="_blank"><img src="https://latex.codecogs.com/svg.latex?d(i,&space;i')&space;=&space;\sum_{i'}&space;k_{ii'}m(i')" title="d(i, i') = \sum_{i'} k_{ii'}m(i')" /></a>
</p>

The original image could be calculated by using relation

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=m&space;=&space;(A^TA&plus;\epsilon&space;I)^{-1}A^Td" target="_blank"><img src="https://latex.codecogs.com/svg.latex?m&space;=&space;(A^TA&plus;\epsilon&space;I)^{-1}A^Td" title="m = (A^TA+\epsilon I)^{-1}A^Td" /></a>
</p>

which A is equivalent to kernel (K) in this case.

#### Result

<p align="center">
    <img src="blured_galaxy.png" width="500px" >
    <br>
    <em>Blured Galaxy</em>
</p>

<p align="center">
    <img src="deblured_galaxy_ep1e-08.png" width="500px" >
    <br>
    <em>Deblured Galaxy</em>
</p>