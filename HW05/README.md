# HW 5: Curved Ray Traveltimg Tomography

## Fast Sweeping Algorithm
Let
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=|\vec{\nabla}u(x)|&space;=&space;f(x)\&space;\Leftrightarrow\&space;|\vec{\nabla}T|&space;=&space;s" target="_blank"><img src="https://latex.codecogs.com/svg.latex?|\vec{\nabla}u(x)|&space;=&space;f(x)\&space;\Leftrightarrow\&space;|\vec{\nabla}T|&space;=&space;s" title="|\vec{\nabla}u(x)| = f(x)\ \Leftrightarrow\ |\vec{\nabla}T| = s" /></a>
</p>

### Discretization
Since
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;|\vec{\nabla}T|^2&space;&=&space;s^2&space;\\&space;\left(\frac{\partial&space;T}{\partial&space;x}\right)^2&space;&plus;&space;\left(\frac{\partial&space;T}{\partial&space;y}\right)^2&space;&=&space;s^2&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;|\vec{\nabla}T|^2&space;&=&space;s^2&space;\\&space;\left(\frac{\partial&space;T}{\partial&space;x}\right)^2&space;&plus;&space;\left(\frac{\partial&space;T}{\partial&space;y}\right)^2&space;&=&space;s^2&space;\end{align*}" title="\begin{align*} |\vec{\nabla}T|^2 &= s^2 \\ \left(\frac{\partial T}{\partial x}\right)^2 + \left(\frac{\partial T}{\partial y}\right)^2 &= s^2 \end{align*}" /></a>
</p>

Consider discrete form with Godunov upwind different scheme
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=[H(T^h_{ij}-T^h_{x,\text{min}})]^2&space;&plus;&space;[H(T^h_{ij}-T^h_{y,\text{min}})]^2&space;=&space;s^2_{ij}h^2\&space;:\&space;h&space;\equiv&space;dx" target="_blank"><img src="https://latex.codecogs.com/svg.latex?[H(T^h_{ij}-T^h_{x,\text{min}})]^2&space;&plus;&space;[H(T^h_{ij}-T^h_{y,\text{min}})]^2&space;=&space;s^2_{ij}h^2\&space;:\&space;h&space;\equiv&space;dx" title="[H(T^h_{ij}-T^h_{x,\text{min}})]^2 + [H(T^h_{ij}-T^h_{y,\text{min}})]^2 = s^2_{ij}h^2\ :\ h \equiv dx" /></a>
</p>

Where
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;T^h_{x,\text{min}}&space;&=&space;\min(T^h_{i-1,j},&space;T^h_{i&plus;1,j})&space;\\&space;T^h_{y,\text{min}}&space;&=&space;\min(T^h_{i,j-1},&space;T^h_{i,j&plus;1})&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;T^h_{x,\text{min}}&space;&=&space;\min(T^h_{i-1,j},&space;T^h_{i&plus;1,j})&space;\\&space;T^h_{y,\text{min}}&space;&=&space;\min(T^h_{i,j-1},&space;T^h_{i,j&plus;1})&space;\end{align*}" title="\begin{align*} T^h_{x,\text{min}} &= \min(T^h_{i-1,j}, T^h_{i+1,j}) \\ T^h_{y,\text{min}} &= \min(T^h_{i,j-1}, T^h_{i,j+1}) \end{align*}" /></a>
</p>

### Gauss-Seidei iterations 
The new travel time could be updated via

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\begin{align*}&space;\bar{T}_{ij}=&space;\begin{cases}&space;\min(T^h_{x,\text{min}},&space;T^h_{y,\text{min}})&space;&plus;&space;s_{ij}h\&space;&:\&space;|T^h_{x,\text{min}}-&space;T^h_{y,\text{min}}|&space;\geq&space;s_{ij}h&space;\\&space;\frac{T^h_{x,\text{min}}&plus;&space;T^h_{y,\text{min}}&plus;\sqrt{2s^2_{ij}h^2-(T^h_{x,\text{min}}-&space;T^h_{y,\text{min}})^2}}{2}&space;\&space;&:\&space;|T^h_{x,\text{min}}-&space;T^h_{y,\text{min}}|&space;<&space;s_{ij}h&space;\end{cases}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;\bar{T}_{ij}=&space;\begin{cases}&space;\min(T^h_{x,\text{min}},&space;T^h_{y,\text{min}})&space;&plus;&space;s_{ij}h\&space;&:\&space;|T^h_{x,\text{min}}-&space;T^h_{y,\text{min}}|&space;\geq&space;s_{ij}h&space;\\&space;\frac{T^h_{x,\text{min}}&plus;&space;T^h_{y,\text{min}}&plus;\sqrt{2s^2_{ij}h^2-(T^h_{x,\text{min}}-&space;T^h_{y,\text{min}})^2}}{2}&space;\&space;&:\&space;|T^h_{x,\text{min}}-&space;T^h_{y,\text{min}}|&space;<&space;s_{ij}h&space;\end{cases}&space;\end{align*}" title="\begin{align*} \bar{T}_{ij}= \begin{cases} \min(T^h_{x,\text{min}}, T^h_{y,\text{min}}) + s_{ij}h\ &:\ |T^h_{x,\text{min}}- T^h_{y,\text{min}}| \geq s_{ij}h \\ \frac{T^h_{x,\text{min}}+ T^h_{y,\text{min}}+\sqrt{2s^2_{ij}h^2-(T^h_{x,\text{min}}- T^h_{y,\text{min}})^2}}{2} \ &:\ |T^h_{x,\text{min}}- T^h_{y,\text{min}}| < s_{ij}h \end{cases} \end{align*}" /></a>
</p>

Sweep the travel time 2D field with for four quadrant of the field and update the new travel time via
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=T^\text{new}_{ij}&space;=&space;\min(T^\text{old}_{ij},&space;\bar{T}_{ij})" target="_blank"><img src="https://latex.codecogs.com/svg.latex?T^\text{new}_{ij}&space;=&space;\min(T^\text{old}_{ij},&space;\bar{T}_{ij})" title="T^\text{new}_{ij} = \min(T^\text{old}_{ij}, \bar{T}_{ij})" /></a>
</p>
The amount of iteration could be either determined from a given maximum number or the while condiiton where a variation between the previous and the current state is so small.

### Example of Travel Time Field from a single source
<p align="center">
    <img src="unit_test/travel_time_1.png" width="400px" >
</p>
<p align="center">
    <img src="unit_test/travel_time_2.png" width="400px" >
</p>
<p align="center">
    <img src="unit_test/travel_time_3.png" width="400px" >
</p>

### Ray path
A ray path be calculated from the travel time field from a source and a given regarded receiver position. The process is simply start with the point of receiver on the travel time field and compute a gradient to update the previous position of the ray until it reach to the source.

<p align="center">
    <img src="unit_test/ray_path1.png" width="400px" >
    <img src="unit_test/ray_path2.png" width="400px" >
    <br>
</p>

### Real Velocity Map 
<p align="center">
    <img src="real_v.png" width="500px" >
    <br>
</p>

### Initial Velocity Model
<p align="center">
    <img src="pb1/model_v0.png" width="500px" >
    <br>
    <em>Initial Velocity Model</em>
</p>

<p align="center">
    <img src="unit_test/trial_T_v0_1.png" width="400px" >
    <img src="unit_test/trial_Ray_v0_1.png" width="400px" >
    <br>
    <em>An example of the one pair of source-receiver from an initial velocity model</em>
</p>


## Pb 1 Steepest descent method with a back-tracking line search


## Pb 2 Nonlinear conjugate gradient method with the Wolfe conditions and quadratic/cubic interpolation for the line search

## Pb 3 Marquardt-Levenberg method