---
title: "Transformer"
linkTitle: "Transformer"
date: 2023-06-09
---

## 2-Winding Power Transformer
The transformer model is composed of an RL-segment and an ideal transformer.
The single line diagram is depicted in the figure below.

Diagram

An RL- segment in time domain is described by

$$
\frac{\mathrm{d} i_L(t)}{\mathrm{d}t} = \frac{v_j(t) - v_k(t) - R \cdot i_L(t)}{L}
$$

Transforming it to DP domain , the RL segment is described by 

$$
  \frac{\mathrm{d} i_L(t)}{\mathrm{d}t} + j\omega \cdot i_L(t)= \frac{v_j(t) - v_k(t) - R \cdot i_L(t)}{L}
$$

Applying the trapezoidal method leads to the finite difference equation:
```math
i_L(t) = \frac{(1-b^2 - (a*R)^2 -j*2b) }{(1+a*R)^2+b^2} \cdot i_L(t - \Delta t) + \frac{a+a^2*R - j *ab}{(1+aR)^2 + b^2} \cdot (v(t) + v(t + \Delta t))
```
with
```math
a = \frac{\Delta t}{2L}, \qquad b = \frac{\Delta t \omega}{2} , \qquad v(t) = v_j(t) - v_k(t)
```

Applying this to the transformer model 

Diagram 

```math
i_{n0}(t) = \frac{(1-b^{2} - (aR)^{2} -j2b) }{(1+aR)^{2}+b^{2}} * i_{n0}(t - \Delta t) + \frac{a+a^{2}R - jab}{(1+aR)^{2} + b^{2}} \cdot (v(t) + v(t + \Delta t))
```
with
```math
v(t) = v_n_0(t) - v'_n_1(t), \qquad v'_n_1(t)= n * v_n_1(t) , \qquad i_n_1(t) = n * i_n_0(t)
```

Substituting for v(t) , we get 
```math
i_{n0}(t) = \frac{a+a^{2}R - jab}{(1+aR)^{2} + b^{2}} \cdot (v_{n0}(t) - n v_{n1}(t)) + \frac{(1-b^{2} - (aR)^{2} -j2b) }{(1+aR)^{2}+b^{2}} \cdot i_{n0}(t - \Delta t) + \frac{a+a^{2} R - j *ab}{(1+aR)^{2} + b^{2}} \cdot (v_{n0}(t-\Delta t) - n v_{n1}(t-\Delta t))
```

This is equivalent to 
```math
i_n_0(t) = Y_e_q * (v_n_0(t) - n*v_n_1(t)) + I_e_q 
```
```math
i_n_1(t) = - n * i_n_0(t)
```

So the admittance matrix becomes 

```math
\begin{pmatrix}
i_{n0} \cr
i_{n1} \cr
\end{pmatrix}
=
\begin{pmatrix}
  Y_{eq} & -n*Y_{eq}\cr
  -n*Y_{eq} & n^2* Y_{eq} \cr
\end{pmatrix} 
* 
\begin{pmatrix}
v_n_0 \cr
v_n_1 \cr
\end{pmatrix} + 
\begin{pmatrix}
I_e_q \cr
-I_e_q \cr
\end{pmatrix}
```
