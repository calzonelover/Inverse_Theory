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

[Source code](utility.py) of the travel time.

### Ray path
A ray path be calculated from the travel time field from a source and a given regarded receiver position. The process is simply start with the point of receiver on the travel time field and compute a gradient to update the previous position of the ray until it reach to the source.

<p align="center">
    <img src="unit_test/ray_path1.png" width="400px" >
    <img src="unit_test/ray_path2.png" width="400px" >
    <br>
</p>

[Source code](ray.py) of a curved ray tracing.

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

The ray legnths in each pair of source-receiver could be parallelize into multiple threads that could easily implement with the built-in Python's multiprocessing to reduce the calculation time in each ray tracing that could taking the .

Note: All of the configuration is in this [file](settings.py).

## Pb 1 Steepest descent method with a back-tracking line search ([Source](pb1/sd.py))

<p align="center">
    <img src="pb1/pk_k0.png" width="400px" >
    <img src="pb1/pk_k9.png" width="400px" >
    <br>
    <em>The blurred gradient fields after applied a uniform filter kernel size of (10, 40)</em>
</p>

<p align="center">
    <img src="pb1/model_v.png" width="400px" >
    <img src="pb1/res.png" width="400px" >
    <br>
    <em>Results from hyper-parameters: alpha0=0.0001 and alpha_decay=0.5</em>
</p>

SD method with backtracing step length yield a reasonable result where we could see a remarkable band around x = 6000 and y ~ 8000 which agree to the real velocity map.

## Pb 2 Nonlinear conjugate gradient method with the Wolfe conditions and quadratic/cubic interpolation for the line search ([Source](pb2/cg.py))

<p align="center">
    <img src="pb2/pk_k0.png" width="400px" >
    <img src="pb2/pk_k19.png" width="400px" >
    <br>
    <em>The blurred gradient fields after applied a uniform filter kernel size of (10, 40)</em>
</p>

<p align="center">
    <img src="pb2/model_v.png" width="400px" >
    <img src="pb2/res.png" width="400px" >
</p>

The iteration that CG method with quad/cube interpolation step length reach to the minimum point much faster than SD method with backtracing. Nevertheless, a gradient field looks similar for both methods in an early iteration except for the latter iteration where it almost plateau in the loss figure.

## Pb 3 Marquardt-Levenberg method

The line search direction from Newton's method are governed by
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;(J^TJ)p^\text{N}_k&space;&\equiv&space;-J^Tr&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;(J^TJ)p^\text{N}_k&space;&\equiv&space;-J^Tr&space;\end{align*}" title="\begin{align*} (J^TJ)p^\text{N}_k &\equiv -J^Tr \end{align*}" /></a>
</p>

The method of Gauss-Newton line search direction in non-linear least square could also be regularized by adding the constrain lambda in the hessian approximation term as
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;(J^TJ&plus;\lambda&space;I)p^\text{LM}_k&space;&\equiv&space;-J^Tr&space;\\&space;p^\text{LM}_k&space;&=&space;-(J^TJ&plus;\lambda&space;I)^{-1}J^Tr&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;(J^TJ&plus;\lambda&space;I)p^\text{LM}_k&space;&\equiv&space;-J^Tr&space;\\&space;p^\text{LM}_k&space;&=&space;-(J^TJ&plus;\lambda&space;I)^{-1}J^Tr&space;\end{align*}" title="\begin{align*} (J^TJ+\lambda I)p^\text{LM}_k &\equiv -J^Tr \\ p^\text{LM}_k &= -(J^TJ+\lambda I)^{-1}J^Tr \end{align*}" /></a>
</p>

The lambda factor is essentially one of the hyper-parameters that has to tune to the proper value.
One way to determine the proper value lambda is to try various quantities and visual the L-Curve which the bottom is the proper hyper-parameter where it does not make the model overfitting the problem and vice versa.

An overwhelming of memory allocation occurs in this problem. The easiest way to fix the problem without any modified the algorithm code is to reduce the size of velocity field as in the [settings](settings.py) file.

### Results from various alpha
Splitting the alpha parameters in the log-scale from xxx to xxx. The amount of trial alpha is 12. 

<p align="center">
    <img src="pb3/alpha0/model_v_alpha0.png" width="400px" >
    <img src="pb3/alpha0/res.png" width="400px" >
</p>
<p align="center">
    <img src="pb3/alpha1/model_v_alpha1.png" width="400px" >
    <img src="pb3/alpha1/res.png" width="400px" >
</p>
<p align="center">
    <img src="pb3/alpha6/model_v_alpha6.png" width="400px" >
    <img src="pb3/alpha6/res.png" width="400px" >
</p>
<p align="center">
    <img src="pb3/alpha10/model_v_alpha10.png" width="400px" >
    <img src="pb3/alpha10/res.png" width="400px" >
</p>
<p align="center">
    <img src="pb3/alpha11/model_v_alpha11.png" width="400px" >
    <img src="pb3/alpha11/res.png" width="400px" >
</p>

### L-Curve
<p align="center">
    <img src="l_curve_unlog.png" width="400px" >
    <img src="l_curve.png" width="400px" >
    <br>
    <em>L-Curve from various trial of alpha with log and linear scale</em>
</p>