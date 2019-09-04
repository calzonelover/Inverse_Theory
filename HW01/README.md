# HW 1 

## Problem 1

The result
```bash
Standard least squares method
Parameter a = 1.0087216227203744, b = -0.5104779304736144
Weight least squares method
Parameter a = 1.0087216227203744, b = -0.5104779304736127
```

## Problem 2
The result
```bash
Standard least squares method
Parameter a = 0.9990946986294076, b = -0.5074347563623611
Weight least squares method
Parameter a = 1.0015782107490603, b = -0.4968111847386467
```

## Problem 3

### Pesudocode

* Notation for the relations
    * <a href="https://www.codecogs.com/eqnedit.php?latex=R_{ki}&space;\leftarrow&space;[I(d-Am_{ki}))]^{-1}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R_{ki}&space;\leftarrow&space;[I(d-Am_{ki}))]^{-1}" title="R_{ki} \leftarrow [I(d-Am_{ki}))]^{-1}" /></a>
    * <a href="https://www.codecogs.com/eqnedit.php?latex=m_{kf}&space;\leftarrow&space;(A^TR_{ki}A)^{-1}A^TR_{ki}d" target="_blank"><img src="https://latex.codecogs.com/svg.latex?m_{kf}&space;\leftarrow&space;(A^TR_{ki}A)^{-1}A^TR_{ki}d" title="m_{kf} \leftarrow (A^TR_{ki}A)^{-1}A^TR_{ki}d" /></a>
    * <a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon_{\text{now}}&space;\leftarrow&space;\frac{||&space;m_{kf}&space;-&space;m_{ki}&space;||}{1&plus;||m_{kf}||}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\epsilon_{\text{now}}&space;\leftarrow&space;\frac{||&space;m_{kf}&space;-&space;m_{ki}&space;||}{1&plus;||m_{kf}||}" title="\epsilon_{\text{now}} \leftarrow \frac{|| m_{kf} - m_{ki} ||}{1+||m_{kf}||}" /></a>

```
func get_model_by_l1r(A, d, m_0, epsilon)
    m_ki = m_0
    while true do
        R_ki <- get_R(A, d, m_ki)
        m_kf <- get_m(A, d, R_ki)
        e_now <- get_e(m_ki, m_kf)
        if e_now < epsilon then
            break
        m_ki <- m_kf
    return m_kf
end func
```