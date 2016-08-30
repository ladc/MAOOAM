#  Modular arbitrary-order ocean-atmosphere model: The Tangent Linear and Adjoint model #

## Description ##

 The Tangent Linear and Adjoint models are implemented in the same way as the nonlinear model, with a tensor storing the different terms. The Tangent Linear (TL) tensor `TL[i][j][k]` is defined as:

    TL[i][j][k] = T[i][k][j] + T[i][j][k]

while the Adjoint (AD) tensor `AD[i][j][k]` is defined as:

    AD[i][j][k] = T[j][k][i] + T[j][i][k]

where `T[i][j][k]` is the tensor of the nonlinear model.

These two tensors are used to compute the trajectories of the models, with the equations

      d Dy[i] / dt =  TL[i][j][k]  y*[k] Dy[j]

    - d Dy[i] / dt =  AD[i][j][k]  y*[k] Dy[j]

with an implicit summation over `j=1..ndim` and `k=0...ndim`. Here, `y*` is the point where the Tangent Linear model is defined, with `y*[0]=1`

## Implementation ##

The two tensors are implemented in the module @{tl_ad_tensor} and are obtained with the functions @{tl_ad_tensor.get_tltensor|get_tltensor} and @{tl_ad_tensor.get_adtensor|get_adtensor}.

These tensors are initialized once when loading the @{tl_ad_tensor} module. The tendencies for the two models are then given by the functions @{tl_ad_tensor.tl_traj|tl_traj(t,ystar,deltay,buf)} and @{tl_ad_tensor.ad_traj|ad_traj(t,ystar,deltay,buf)}. 

A fourth-order Runge-Kutta integrator is provided in @{rk4_tl_ad}. 
