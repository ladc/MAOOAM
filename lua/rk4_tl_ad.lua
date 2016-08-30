-- rk4_tl_ad.lua
-- (C) 2015-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- A classical fourth order Runge-Kutta scheme for the TL/AD model. 
-- Butcher tableau:
--
--    0   |
--    1/2 | 1/2
--    1/2 | 0    1/2
--    1   | 0    0    1
--        +-------------------
--          1/6  1/3  1/3  1/6
------------------------------------------------------------------------

--- Create a TL/AD integrator for n-arrays.
-- y is also evolved to calculate the correct TL.
-- Note that the function `f_tl` has more arguments than a regular integrator.
-- You need to pass both has both the trajectory y (also ystar) and the
-- perturbation deltay:
-- `f_tl(t,y,deltay,buf)` instead of `f(t,y,buf)`
-- @function rk4_tl_ad
-- @param n number of variables (length of the array)
-- @return Integrator for n-arrays
return function(n)
  local array = require("array")(n)
  local buf_kA, buf_kB, buf_y1 = array(), array(), array()
  local buf_kAA, buf_kBB, buf_y11 = array(), array(), array()
  --- Integrator for n-arrays
  -- @function integrator
  -- @param y variables at time t
  -- @param deltay perturbation at time t
  -- @param f function to calculate the time derivatives of the variables
  -- @param f_tl function to calculate the time derivatives of the tangent linear model
  -- @param t time
  -- @param dt time integration step
  -- @param ynew n-array (buffer) to store the new y value 
  -- @param deltaynew n-array (buffer) to store the new perturbation value 
  -- @return t+dt (incremented time)
  -- @return ynew array with variables at time t+dt
  -- @return deltaynew array with perturbation at time t+dt
  return function(y,deltay,f,f_tl,t,dt,ynew,deltaynew)
    f(t,y,buf_kA)                               -- k_1 (A)  = f(t, y)
    f_tl(t,y,deltay,buf_kAA)                    -- k_11 (AA) = f_tl(t, y, dy)

    buf_kA:nmul(0.5*dt,buf_y1):add(y,buf_y1)    -- y_1       = y + 0.5*dt*k_1 (A)
    buf_kAA:nmul(0.5*dt,buf_y11):add(deltay,buf_y11)-- y_11      = dy + 0.5*dt*k_11 (AA)

    f(t+0.5*dt,buf_y1,buf_kB)                       -- k_2 (B) = f(t + 0.5*dt, y_1)
    f_tl(t+0.5*dt,buf_y1,buf_y11,buf_kBB)    -- k_22 (BB) = f_tl(t + 0.5*dt, y_1, y_11)

    buf_kB:nmul(0.5*dt,buf_y1):add(y,buf_y1) -- y_1'    = y + 0.5*dt*k_2
    buf_kBB:nmul(0.5*dt,buf_y11):add(deltay,buf_y11) -- y_11'  = dy + 0.5*dt*k_22

    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_1 (A) + 2 * k_2 (B)
    buf_kBB:nmul(2,buf_kBB):add(buf_kAA,buf_kAA) -- k_SS (AA) = k_11 (AA) + 2 * k_22 (BB)

    f(t+0.5*dt,buf_y1,buf_kB)                -- k_3 (B) = f(t + 0.5*dt, y_1')
    f_tl(t+0.5*dt,buf_y1,buf_y11,buf_kBB)    -- k_33 (BB) = f_tl(t + 0.5*dt, y_1', y_11')

    buf_kB:nmul(dt,buf_y1):add(y,buf_y1)     -- y_1''   = y + dt*k_3 (B)
    buf_kBB:nmul(dt,buf_y11):add(deltay,buf_y11)    -- y_11''   = dy + dt*k_33 (BB)

    buf_kB:nmul(2,buf_kB):add(buf_kA,buf_kA) -- k_S (A) = k_S (A) + 2 * k_3 (B)
    buf_kBB:nmul(2,buf_kBB):add(buf_kAA,buf_kAA) -- k_SS (AA) = k_SS (AA) + 2 * k_33 (BB)

    f(t+dt,buf_y1,buf_kB)                    -- k_4 (B) = f(t + dt, y_1'')
    f_tl(t+dt,buf_y1,buf_y11,buf_kBB)        -- k_44 (BB) = f_tl(t + dt, y_1'',y_11'')

    buf_kB:add(buf_kA,buf_kA)                -- k_S (A) = k_S (A) + k_4 (B)
    buf_kBB:add(buf_kAA,buf_kAA)             -- k_SS (AA) = k_SS (AA) + k_44 (BB)

     -- dynew = dy + dt/6 * k_S
    return t+dt,
           buf_kA:nmul(dt/6,buf_kA):add(y,ynew),
           buf_kAA:nmul(dt/6,buf_kAA):add(deltay,deltaynew)
  end
end

