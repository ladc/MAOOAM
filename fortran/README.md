# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Fortran implementation
## Lyapunov vectors computation version

------------------------------------------------------------------------

## About ##

(c) 2013-2016 Lesley De Cruz, Jonathan Demaeyer and Sebastian Schubert

See [LICENSE.txt](../LICENSE.txt) for license information.

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

## Note on the Lyapunov vectors computation

This version of the code allows for the computation of the backward and 
forward Lyapunov vectors. The method used is described in

* Benettin, G., Galgani, L., Giorgilli, A., and Strelcyn, J. M. : Lyapunov
characteristic exponents for smooth dynamical systems; a method for computing
all of them. Part 2: Numerical application. \a Meccanica \a, 15, 21-30,
[doi:10.1007/BF02128237](http://dx.doi.org/10.1007/BF02128237), 1980.

------------------------------------------------------------------------

## Installation ##

The program can be installed with Makefile. We provide configuration files for 
two compilers : gfortran and ifort.

By default, gfortran is selected. To select one or the other, simply modify the 
Makefile accordingly. If gfortran is selected, the code should be compiled 
with gfortran 4.7+ (allows for allocatable arrays in namelists). 
If ifort is selected, the code has been tested with the version 14.0.2 and we 
do not guarantee compatibility with older compiler version.

To install, unpack the archive in a folder, and run:
     make
 
 Remark: The command "make clean" removes the compiled files.

For Windows users, a minimalistic GNU development environment
 (including gfortran and make) is available at [www.mingw.org](http://www.mingw.org) .

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called aotensor which
includes all the coefficients. This tensor is computed once at the program
initialization.

* maooam.f90 : Main program.
* maooam_lyap.f90 : Version of the main program with the computation of the Lyapunov vectors and exponents.
* maooam_lyap_div.f90 : Version of the main program with the computation of the first Lyapunov exponents with the divergence method.
* aotensor_def.f90 : Tensor aotensor computation module.
* lyap_vectors.f90 : Module for the computation of Lyapunov exponents and vectors
* IC_def.f90 : A module which loads the user specified initial condition.
* inprod_analytic.f90 : Inner products computation module.
* rk2_integrator.f90 : A module which contains the Heun integrator for the model equations.
* rk4_integrator.f90 : A module which contains the RK4 integrator for the model equations.
* Makefile : The Makefile.
* gfortran.mk : Gfortran compiler options file.
* ifort.mk : Ifort compiler options file.
* params.f90 : The model parameters module.
* tl_ad_tensor.f90 : Tangent Linear (TL) and Adjoint (AD) model tensors definition module
* rk2_tl_ad_integrator.f90 : Heun Tangent Linear (TL) and Adjoint (AD) model integrators module
* rk4_tl_ad_integrator.f90 : RK4 Tangent Linear (TL) and Adjoint (AD) model integrators module
* test_tl_ad.f90 : Tests for the Tangent Linear (TL) and Adjoint (AD) model versions
* README.md : The present file.
* LICENSE.txt : The license text of the program.
* util.f90 : A module with various useful functions.
* tensor.f90 : Tensor utility module.
* stat.f90 : A module for statistic accumulation.
* params.nml : A namelist to specify the model parameters.
* int_params.nml : A namelist to specify the integration parameters.
* modeselection.nml : A namelist to specify which spectral decomposition will be used.

A documentation is available [here](./doc/html/index.html) (html) and [here](./doc/latex/Reference_manual.pdf) (pdf).
 
------------------------------------------------------------------------

## Usage ##

The user first has to fill the params.nml and int_params.nml namelist files according to their needs.
Indeed, model and integration parameters can be specified respectively in the params.nml and int_params.nml namelist files. Some examples related to already published article are available in the params folder.

The parameters related to the Lyapunov vectors and exponents computation are located in the int_params.nml namelist file.

The modeselection.nml namelist can then be filled : 
* NBOC and NBATM specify the number of blocks that will be used in respectively the ocean and
  the atmosphere. Each block corresponds to a given x and y wavenumber.
* The OMS and AMS arrays are integer arrays which specify which wavenumbers of
  the spectral decomposition will be used in respectively the ocean and the
  atmosphere. Their shapes are OMS(NBOC,2) and AMS(NBATM,2).
* The first dimension specifies the number attributed by the user to the block and the second
  dimension specifies the x and the y wavenumbers.
* The VDDG model, described in Vannitsem et al. (2015) is given as an example
  in the archive.
* Note that the variables of the model are numbered according to the chosen
  order of the blocks.


Finally, the IC.nml file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current IC.nml file (if it exists)
and run the main program :

    ./maooam

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart :

    ./maooam

It will generate two files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)

The tangent linear and adjoint models of MAOOAM are provided in the
tl_ad_tensor, rk2_tl_ad_integrator and rk4_tl_ad_integrator modules. It is
documented [here](./doc/html/md_tl_ad_doc.html).

The computation of the Lyapunov vectors can be done by running the program :

    ./maooam_lyap

It will generate four files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)
 * lyapunov_exponents.dat : the recorded time evolution of the Lyapunov exponents
 * mean_lyapunov.dat : the mean Lyapunov exponents as well as their variance.

Alternatively, the first Lyapunov exponent can be computed by the divergence method, 
by running :

   ./maooam_lyap_div

It will generate three files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)
 * lyapunov_exponents_div.dat : the recorded time evolution of the first Lyapunov exponent

------------------------------------------------------------------------

## Implementation notes ##

As the system of differential equations is at most bilinear in y[j] (j=1..n), y
being the array of variables, it can be expressed as a tensor contraction
(written using Einstein convention, i.e. indices that occur twice on one side
of an equation are summed over):

    dy  / dt =  T        y   y      (y  == 1)
      i          i,j,k    j   k       0

The tensor T that encodes the differential equations is composed so that:

* T[i][j][k] contains the contribution of dy[i]/dt proportional to y[j]*y[k].
* Furthermore, y[0] is always equal to 1, so that T[i][0][0] is the constant
contribution to var dy[i]/dt.
* T[i][j][0] + T[i][0][j] is the contribution to  dy[i]/dt which is linear in
y[j].

Ideally, the tensor is composed as an upper triangular matrix 
(in the last two coordinates).

The tensor for this model is composed in the aotensor_def module and uses the
inner products defined in the inprod_analytic module.


------------------------------------------------------------------------

## Final Remarks ##

The authors would like to thank Kris for help with the lua2fortran project. It
has greatly reduced the amount of (error-prone) work.

  No animals were harmed during the coding process.
