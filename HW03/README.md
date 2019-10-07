# HW 3

## Problem 1: Tomography

<p align="center">
    <img src="pb1/real_v.png" width="500px" >
    <br>
    <em>Real velocity map of the stone layer</em>
</p>

### Finding a ray's path length

* pic with reference
* eq

<p align="center">
    <img src="pb1/img/out.gif" width="500px" >
    <br>
    <em>Fix source and shifting receiver along y-axis</em>
</p>

<p align="center">
    <img src="pb1/img/s.gif" width="500px" >
    <br>
    <em>Measured wave velocity from various alpha</em>
</p>

<p align="center">
    <img src="pb1/l_curve_origin.png" width="400px" >
    <img src="pb1/l_curve_ref.png" width="400px" >
    <br>
    <em>L-Curve of getting proper alpha from origin and shifted reference</em>
</p>

#### Corner compare

<p align="center">
    <img src="pb1/img/model_v_alpha024.png" width="250px" >
    <img src="pb1/real_v.png" width="250px" >
    <img src="pb1/img/model_v_alpha015.png" width="250px" >
    <br>
    <em>Real velocity map of the stone layer</em>
</p>

### 1.1 Least Square with Tikhonov Regularization

### 1.2 Applying linear SD and CG inverse problem

#### SD

<p align="center">
    <img src="pb1/sd_r.png" width="500px" >
    <br>
    <em>Measured wave velocity from SD</em>
</p>

<p align="center">
    <img src="pb1/v_sd.png" width="500px" >
    <br>
    <em>Measured wave velocity from SD</em>
</p>


#### CGLS

since bla

<p align="center">
    <img src="pb1/cg_r.png" width="500px" >
    <br>
    <em>Measured wave velocity from CGLS</em>
</p>

<p align="center">
    <img src="pb1/v_cg.png" width="500px" >
    <br>
    <em>Measured wave velocity from CGLS</em>
</p>


## Problem 2: Optimize Rosenblock function

* Rosenblock

### 2.1 Fix step-length

<p align="center">
    <img src="pb2/1_fixalpha/sd_r.png" width="400px" >
    <img src="pb2/1_fixalpha/sd.png" width="400px" >
    <br>
    <em>SD</em>
</p>

<p align="center">
    <img src="pb2/1_fixalpha/cg_r.png" width="400px" >
    <img src="pb2/1_fixalpha/cg.png" width="400px" >
    <br>
    <em>CG</em>
</p>

### 2.2 Backtracking

<p align="center">
    <img src="pb2/2_backtracking/sd_r.png" width="400px" >
    <img src="pb2/2_backtracking/sd.png" width="400px" >
    <br>
    <em>SD</em>
</p>

<p align="center">
    <img src="pb2/2_backtracking/cg_r.png" width="400px" >
    <img src="pb2/2_backtracking/cg.png" width="400px" >
    <br>
    <em>CG</em>
</p>

### 2.3 Quadratic/Cubic

<p align="center">
    <img src="pb2/3_quad/sd_r.png" width="400px" >
    <img src="pb2/3_quad/sd.png" width="400px" >
    <br>
    <em>SD</em>
</p>

* report cg too small alpha

<p align="center">
    <img src="pb2/3_quad/cg_r.png" width="400px" >
    <img src="pb2/3_quad/cg.png" width="400px" >
    <br>
    <em>CG</em>
</p>

### 2.4 Newton

<p align="center">
    <img src="pb2/4_newton/newton_r.png" width="400px" >
    <img src="pb2/4_newton/newton.png" width="400px" >
    <br>
    <em>Newton</em>
</p>

