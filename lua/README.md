# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Lua implementation

------------------------------------------------------------------------

## About ##

(c) 2013-2016 Lesley De Cruz and Jonathan Demaeyer

See LICENSE.txt for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
[doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

**Please cite this article if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <lesley.decruz@meteo.be>, <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the MAOOAM [code repository](http://www.github.com/Climdyn/MAOOAM)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

------------------------------------------------------------------------

## Installation ##

Unpack the archive in a folder. There is no compilation required to run the
model once you have LuaJIT. 

Installation instructions for LuaJIT can be found on the [LuaJIT
website](http://luajit.org/install.html). In short, you can install it with the
following commands (on Linux):

```
git clone http://luajit.org/git/luajit-2.0.git
cd luajit-2.0
make && sudo make install
```

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called aotensor which
includes all the coefficients. This tensor is computed once at the program
initialization.

* `maooam.lua` : Main program.
* `aotensor.lua` : Tensor aotensor computation module.
* `inprod_analytic.lua` : Inner products computation module.
* `restore.lua` : Function to restore a previous model state.
* `write_IC.lua` : Function to write an initial condition template.
* `array.lua` : Module for arrays and their mathematical operations.
* `rk2.lua` : Module with the Heun (second order Runge-Kutta) ODE integrator.
* `rk4.lua` : Module with the fourth-order Runge-Kutta integrator.
* `tensor.lua` : Sparse tensor tools.
* `gz.lua` : Module with basic bindings to zlib.
* `stat.lua` : A module for statistic accumulation.
* `tl_ad_tensor.lua` : Module which defines the Tangent Linear and Adjoint models.
* `README.md` : The present file.
* `TL_AD.md` : Information on the tangent linear and adjoint models.
* `LICENSE.txt` : The license text of the program.
* `params.lua` : A configuration file to specify the integration and model parameters.
* `modeselection.lua` : A configuration file to specify which spectral decomposition will be used.
* `test_inprod_analytic.lua` : Program to write out the inner products. For testing purposes.
* `test_aotensor.lua` : Program to write out the inner products. For testing purposes.
* `test_tl_ad.lua` : Program which tests the Tangent Linear and Adjoint models.

------------------------------------------------------------------------

## Usage ##

Model and integration parameters can be specified in @{params|params.lua}.

The resolution of the model can be modified by specifying which modes should be
included in the file @{modeselection|modeselection.lua}.

The initial conditions for each of the variables are defined in IC.lua.
To obtain an example of this file corresponding to the model you have defined
in @{modeselection|modeselection.lua}, simply delete the current IC.lua file (if
it exists) and run the program. It will generate a new one and exit. Just fill
the newly generated one.

To run a simulation, just run:
    luajit maooam.lua

In Windows:
    luajit.exe maooam.lua

This will generate several files:

* `output_maooam_meanfields.txt` : the @{stat|mean and variance} (climatology) of the variables
* `output_maooam_snapshot.txt` : snapshots of the dynamical state of the model
* `output_maooam_trajectory.txt_NNNN.gz` : the recorded time evolution of the
  variables (gzipped if compression is enabled).

To continue a previous run, use the "continue" argument:
    luajit maooam.lua continue

The tangent linear and adjoint models of MAOOAM are provided in
@{tl_ad_tensor|tl_ad_tensor.lua}.

The LuaJIT source code can be obtained [here](http://luajit.org/download.html).

------------------------------------------------------------------------

## Implementation notes ##

As the system of differential equations is at most bilinear in `y[j] (j=1..n)`, `y`
being the `array` of variables, it can be expressed as a `tensor` contraction
(written using Einstein convention, i.e. indices that occur twice on one side
of an equation are summed over):

    dy  / dt =  T        y   y      (y  == 1)
      i          i,j,k    j   k       0

The tensor T that encodes the differential equations is composed so that:

* `T[i][j][k]` contains the contribution of `dy[i]/dt` proportional to `y[j]*y[k]`.
* Furthermore, `y[0]` is always equal to 1, so that `T[i][0][0]` is the constant
contribution to var `dy[i]/dt`.
* `T[i][j][0] + T[i][0][j]` is the contribution to `dy[i]/dt` which is linear in
`y[j]`.

Ideally, the tensor is composed as an upper triangular matrix (in the last two
coordinates).

The tensor for this model is composed in @{aotensor|aotensor.lua} and uses the
inner products defined in @{inprod_analytic|inprod_analytic.lua}.

