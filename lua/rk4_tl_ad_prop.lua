-- rk4_tl_ad_prop.lua
-- (C) 2015-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- A classical fourth order Runge-Kutta scheme. 
-- Butcher tableau:
--
--    0 |
--  1/2 | 1/2
--  1/2 | 0    1/2
--    1 | 0    0    1
--      +-------------------
--        1/6  1/3  1/3  1/6
------------------------------------------------------------------------

local m = require("lgsl.matrix")

--- Create a TL propagator integrator for n-arrays.
-- y is also evolved to calculate the correct TL.
-- Note that the function f_tl has more arguments than a regular integrator.
-- You need to pass both has both the trajectory y (also ystar) and the
-- perturbation deltay:
-- f_tl(t,y,deltay,buf) instead of f(t,y,buf)
-- @function rk4_tl_ad
-- @param n number of variables (length of the array)
-- @return Integrator for n-arrays
return function(n)
  local array = require("array")(n)
  local buf_kA, buf_kB, buf_y1 = array(), array(), array()
  local buf_j1, buf_j2, buf_j3, buf_j4 = m.new(n-1,n-1), m.new(n-1,n-1),
    m.new(n-1,n-1), m.new(n-1,n-1)
  local one = m.unit(n-1)
  --- Integrator for n-arrays
  -- @function integrator
  -- @param y variables at time t
  -- @param f function to calculate the time derivatives of the variables
  -- @param J Jacobian function
  -- @param t time
  -- @param dt time integration step
  -- @param ynew n-array (buffer) to store the new y value 
  -- @return t+dt (incremented time)
  -- @return ynew array with variables at time t+dt
  return function(y,f,J,t,dt,ynew,propagator)
    f(t,y,buf_kA)                               -- k_1 (A)  = f(t, y)
    J(y,buf_j1)

    buf_kA:nmul(0.5*dt,buf_y1):add(y,buf_y1)    -- y_1       = y + 0.5*dt*k_1 (A)

    f(t+0.5*dt,buf_y1,buf_kB)                       -- k_2 (B) = f(t + 0.5*dt, y_1)
    J(buf_y1,buf_j2)

    buf_kB:nmul(0.5*dt,buf_y1):add(y,buf_y1) -- y_1'    = y + 0.5*dt*k_2

    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_1 (A) + 2 * k_2 (B)

    f(t+0.5*dt,buf_y1,buf_kB)                -- k_3 (B) = f(t + 0.5*dt, y_1')
    J(buf_y1,buf_j3)

    buf_kB:nmul(dt,buf_y1):add(y,buf_y1)     -- y_1''   = y + dt*k_3 (B)

    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_S (A) + 2 * k_3 (B)

    f(t+dt,buf_y1,buf_kB)                    -- k_4 (B) = f(t + dt, y_1'')
    J(buf_y1,buf_j4)

    buf_kB:add(buf_kA,buf_kA)                -- k_S (A) = k_S (A) + k_4 (B)

    -- propagator = one + dt/6 *(buf_j1 + 2*buf_j2 + 2*buf_j3 + buf_j4) -- wasteful: will create many matrices

    -- explicitly add the matrices to avoid creating too many intermediate objects
    for i=1,n-1 do
      for j=1,n-1 do
        propagator:set(i,j,one:get(i,j) + dt/6* (buf_j1:get(i,j) + 2*buf_j2:get(i,j) +
          2*buf_j3:get(i,j) + buf_j4:get(i,j)))
      end
    end

     -- dynew = dy + dt/6 * k_S
    return t+dt, 
           buf_kA:nmul(dt/6,buf_kA):add(y,ynew),
           propagator
  end
end

