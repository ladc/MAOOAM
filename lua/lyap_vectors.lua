-- lyap_vectors.lua
-- (C) 2016 Sebastian Schubert, Jonathan Demaeyer & Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- This module contains the necessary tools to perform the Bennettin
-- steps to compute the lyapunov exponents. (Ginelli for CLV will be added
-- later) 
------------------------------------------------------------------------

local m = require("lgsl.matrix")
local gsl = require("lgsl.gsl")
local gsl_check = require("lgsl.gsl-check")
local ffi = require("ffi")

local rand, log, abs = math.random, math.log, math.abs
local NT=gsl.CblasNoTrans

-- fill a matrix mat with uniform random numbers in [-1, 1]
-- @param mat
local function fillrand(mat)
  local n1,n2=mat:dim()
  for i=1,n1 do
    for j=1,n2 do
      mat:set(i,j, 2*rand() - 1)
    end
  end
end

return function(n)
  -- Initialize buffers
  local prop_buf, prop = m.new(n-1,n-1), m.unit(n-1)
  local ensemble = m.new(n-1,n-1)
  local tau = ffi.gc(gsl.gsl_vector_alloc(n-1), gsl.gsl_vector_free)
  local lyapunov, loclyap = m.new(n-1,1), m.new(n-1,1)
  fillrand(ensemble)
  gsl_check(gsl.gsl_linalg_QR_decomp(ensemble,tau))

  -- Multiply prop with the next propagator
  -- prop = prop_mul x (prop = prop_buf)
  -- @param prop_mul propagator to multiply with prop
  local function multiply_prop(prop_mul)
    m.set(prop_buf,prop) -- prop_buf <- prop
    -- prop <- prop_mul * prop_buf
    gsl_check(gsl.gsl_blas_dgemm(NT,NT,1,prop_mul,prop_buf,0,prop))
  end

  -- Performs the Benettin step in integration. Multiplies the aggregated
  -- propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  -- orthogonalization gives Q and upper triangular matrix R). Computes also the
  -- Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  -- the subroutine and restored to a unit matrix.
  -- In short:
  -- compute ensemble = prop * Q(stored in ensemble)
  -- compute QR decomposition of ensemble, store in ensemble and tau
  -- reset prop to unit mat
  local function benettin_step(rescaling_time)
    -- transpose prop
    gsl_check(gsl.gsl_matrix_transpose(prop))
    -- compute Q^T prop^T
    gsl_check(gsl.gsl_linalg_QR_QTmat(ensemble, tau, prop))
    -- transpose to obtain prop = prop * Q
    gsl_check(gsl.gsl_matrix_transpose(prop))
    -- switch prop and ensemble: now ensemble = prop * Q
    ensemble, prop = prop, ensemble
    -- do a new QR decomposition of ensemble
    gsl_check(gsl.gsl_linalg_QR_decomp(ensemble,tau))
    -- reinitialize prop
    gsl.gsl_matrix_set_identity(prop)
    for i=1,n-1 do
      loclyap:set(i,1,log(abs(ensemble:get(i,i)))/rescaling_time)
    end
    return loclyap
  end

  -- Return the current state of the Lyapunov computation:
  -- accumulated propagator and the ensemble of (orthonormal) Lyapunov vectors.
  local function get_lyap_state()
    return ensemble, prop
  end
  
  -- Set the current state of the Lyapunov computation:
  -- accumulated propagator and the ensemble of (orthonormal) Lyapunov vectors.
  local function set_lyap_state(ens,pr)
    ensemble = ens or ensemble
    prop = pr or prop
  end

  return {
    benettin_step = benettin_step,
    multiply_prop = multiply_prop,
    get_lyap_state = get_lyap_state,
    set_lyap_state = set_lyap_state,
  }
end

