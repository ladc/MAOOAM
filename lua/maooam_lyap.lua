#!/usr/bin/env luajit
-- maooam_lyap.lua
-- (C) 2013-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Lua implementation of the modular arbitrary-order ocean-atmosphere
-- model MAOOAM allowing for the simultaneous integration of both the nonlinear
-- and tangent linear model. This version of the model also computes the
-- Lyapunov spectrum.
--
-- Run as follows:
-- luajit maooam_lyap.lua
--
-- To continue a previous run, use the "continue" argument:
-- luajit maooam_lyap.lua continue
--
-- Optimised for LuaJIT (http://luajit.org).
--
-- See @{readme.md} for more details.
------------------------------------------------------------------------

local params = require("params")

-- Tensor with the coefficients of the nonlinear (polynomial) system of
-- diffeqs. N.B.: the model tensor is fixed, so changes to model parameters
-- will only be effective up to this point!
local n, aotensor, sparse_mul3
do
  local inprod = require("inprod_analytic")
  n = inprod.natm*2+inprod.noc*2 + 1
  local tensor = require("tensor")
  aotensor = tensor.coo_to_fficoo(tensor.simplify_coo(require("aotensor")))
  sparse_mul3 = tensor.sparse_mul3
  require("write_IC")(inprod)
end

-- Utilities
local array = require("array")(n) -- n-array constructor
local evolve = require(params.i.integrator)(n) -- integrator
local evolve_prop = require("rk4_tl_ad_prop")(n) -- combined integrator for the model and the propagator matrix (TL/AD)
local lyap_vectors = require("lyap_vectors")(n) -- Utilities: benettin_step multiply_prop

-- checks if the value x, if truthy, is a number, and print a message if falsy.
local function check_number(x,name)
  name = name or ""
  if x then
    return assert(tonumber(x),"Error in params: "..name.." should be a number if truthy.")
  else
    io.write("* ",name," disabled. \n")
    return x
  end
end

-- IO functions
local format = string.format
local function fprintf(f,...) return f:write(format(...)) end

local function write(firstcol,y,f,fmt)
  f = f or io.stdout
  fmt = fmt or "\t%.12f"
  f:write(firstcol)
  for i=1,n-1 do fprintf(f,fmt,y[i]) end
  f:write("\n")
end

local clock0, clockprev = os.clock(), os.clock()
local walltime = check_number(params.i.walltime,"walltime")

-- Take a full-precision snapshot of the current state (y) and statistics (stats)
-- Return false when time's up.
local function take_snapshot(t,y,f,stats)
  write(t,y,f," \t%13a")
  if stats then
    write("mean", stats.mean(), f, "\t%13a")
    write("var ", stats.var() , f, "\t%13a")
    fprintf(f,"%13a\n",stats.iter())
  end
  f:flush()
  -- We still have time if walltime is still more than clockstep(+10%) away
  if walltime then
    local clock = os.clock()
    local clockstep = clock - clockprev
    clockprev = clock
    return clock - clock0 + clockstep*1.5 < walltime
  else
    return true
  end
end

local function save_matrix(matrix,f)
  for i=1,n-1 do write("",matrix[i],f,"\t%13a") end
end

local function close(f) if f and f~=io.stdout then f:close() end end

------------------------------------------------------------------------
-- Integrate coupled ocean-atmosphere model.
------------------------------------------------------------------------

--- Function that calculates the time derivative of the n variables.
-- The function reduces to a sparse tensor contraction due to the bilinear
-- nature of the equations.
-- @param t time
-- @param y array with variables at time t
-- @param buf n-arrayi (buffer) to store derivatives.
local function ao(t,y,buf)
  return sparse_mul3(aotensor,y,y,buf)
end

-- Check some input values.
local writeout = check_number(params.i.writeout,"writeout")
local snapshot = check_number(params.i.snapshot,"snapshot")
if walltime then
  assert(snapshot,"If a walltime is defined, you must also specify the snapshot interval. Please edit your parameters.")
end
local statistics = check_number(params.i.statistics,"statistics")
local st = statistics and require("stat")(n) -- statistics object
local stlyap = statistics and require("stat")(n) -- statistics object

-- Get initial conditions and statistics.
local t_init, y
do
  local initialconditions
  -- If it's a continuation of a previous run, restore the t, statistics and IC
  -- from the snapshot file.
  if arg[1]=="continue" then
    t_init, initialconditions = require("restore")(params, st)
    if t_init > params.i.t_trans then
      -- only restore statistics for lyap
      require("restore")(params, stlyap, params.i.getoutfn("_snapshot_lyap"))
      -- restore the ensemble of lyap vectors
      lyap_vectors.set_lyap_state(require("read_matrix")(params.i.getoutfn("_lyap_state")))
    end
  end
  -- No continuation or no snapshot file: read IC file.
  if not initialconditions then
    t_init = 0
    initialconditions = require("IC")
  end
  assert(#initialconditions==n-1, "Dimension of initial conditions ".. #initialconditions.." ~= expected number "..n-1)
  initialconditions[0] = 1
  y = array({data=initialconditions})
end

-- Set up snapshot files
local snapf, snaplyapf, snapensf
local mode = arg[1]=="continue" and "a" or "w" -- append or (over)write
if snapshot then
  local snapfn = params.i.getoutfn("_snapshot")
  local snaplyapfn = params.i.getoutfn("_snapshot_lyap")
  local snapensfn = params.i.getoutfn("_snapshot_ens")
  snapf = io.open(snapfn,mode)
  snaplyapf = io.open(snaplyapfn,mode)
  snapensf = io.open(snapensfn,mode)
  fprintf(io.stdout,"* Writing snapshots to: %s, %s, %s\n", snapfn, snaplyapfn, snapensfn)
end

local dt = params.i.dt
-- Perform integration for transitory period (no writeout).
-- Only if t_init does not exceed the transitory period
if t_init <= params.i.t_trans then 
  fprintf(io.stdout,"* Starting transient period of %e time units from time %e.\n",params.i.t_trans,t_init)
  for t=t_init,params.i.t_trans,dt do
    if snapshot and t%snapshot<dt then take_snapshot(t,y,snapf,st) end
    evolve(y,ao,t,dt,y)
  end
  fprintf(io.stdout,"* Finished transient period, starting run of %e time units.\n",params.i.t_run)
  t_init = 0
else -- start counting from zero at the end of the transitory period.
  t_init = t_init - params.i.t_trans
end

for i=1,n do
  assert(y[i]==y[i], "NaN encountered; the system is unstable. Please modify your initial conditions, forcing or coupling parameter.")
end

-- Setup output files.
local logf
if writeout then -- setup trajectory file
  local logfn_base = params.i.getoutfn("_trajectory")
  local logfn = logfn_base
  if params.i.compression then
    -- If compression is enabled: always create a new file to avoid corrupting the archive.
    local i = 0
    repeat
      logfn = string.format("%s_%04d.gz",logfn_base,i)
      logf = io.open(logfn,"r") -- check if file exists
      i = i+1
    until not (logf and logf:close()) -- always close it if it exists
    logf = require("gz").open(logfn,"wb") -- binary mode
  else
    logf = io.open(logfn,mode)
  end
  if logf then fprintf(io.stdout,"* Writing trajectory to: %s\n",logfn) end
end

logf = logf or io.stdout

-- Jacobian matrix, propagator and lyapunov vector buffers
local m = require("lgsl.matrix")
local jmul = require("jmul").arr
local function J(v,buf_J)
  return jmul(aotensor,v,buf_J)
end
local rescaling_time = params.i.rescaling_time
local prop_buf = m.new(n-1,n-1)

local copy=require("ffi").copy

-- Utilities for computation of Lyapunov exponents
local benettin_step,multiply_prop,get_lyap_state = 
  lyap_vectors.benettin_step,lyap_vectors.multiply_prop,lyap_vectors.get_lyap_state


-- Perform integration (and write out).
for t=t_init,params.i.t_run,dt do
  if snapshot   and t%snapshot<dt then
    -- Save the lyapunov exponents spectrum: mean and variance
    take_snapshot(t+params.i.t_trans,y,snaplyapf,stlyap)
    -- save the matrix containing the ensemble of lyap vectors
    save_matrix(get_lyap_state(),snapensf)
    -- take_snapshot returns false if our walltime is up.
    if not take_snapshot(t+params.i.t_trans,y,snapf,st) then
      close(logf)
      close(snapf)
      close(snaplyapf)
      io.write("- Clean exit due to walltime restrictions.\n")
      os.exit(0)
    end
  end
  if writeout   and t%writeout<dt then write(t,y,logf) end
  if statistics and t%statistics<dt then st.acc(y) end
  -- Evolve both y and the propagator; accumulate the internal propagator.
  local _, loclyap_arr = 0, array()
  _,y,prop_buf= evolve_prop(y,ao,J,t,dt,y,prop_buf)
  multiply_prop(prop_buf)
  -- perform the Benettin step every rescaling_time.
  if rescaling_time and t%rescaling_time<dt then
    local loclyap = benettin_step(rescaling_time)
    copy(loclyap_arr.data+1,loclyap.data,(n-1)*8)
    -- Compute statistics of the local lyapunov exponents.
    stlyap.acc(loclyap_arr)
  end
end

-- Finalize output files
if writeout then write(params.i.t_run+dt,y,logf) end
close(logf)
if snapshot then
  take_snapshot(params.i.t_run+params.i.t_trans+dt,y,snapf,st)
  take_snapshot(params.i.t_run+params.i.t_trans+dt,y,snaplyapf,stlyap)
end
close(snapf)
close(snaplyapf)

-- Write out statistics
if statistics then
  local logfn = params.i.getoutfn("_meanfields")
  logf = io.open(logfn,mode) or io.stdout
  fprintf(io.stdout,"* Writing statistics to: %s\n", logfn)
  write("mean", st.mean(),logf)
  write("var ", st.var(),logf)
  fprintf(logf,"samples\t %d\n", st.iter())
  close(logf)
  logfn = params.i.getoutfn("_meanlyap")
  logf = io.open(logfn,mode) or io.stdout
  fprintf(io.stdout,"* Writing lyapunov statistics to: %s\n", logfn)
  write("mean", stlyap.mean(),logf)
  write("var ", stlyap.var(),logf)
  fprintf(logf,"samples\t %d\n", stlyap:iter())
  close(logf)
end
