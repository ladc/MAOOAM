-- restore.lua
-- (C) 2014-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- The restore function restores the dynamical state and the statistics
-- of a previous run from a snapshot file.
------------------------------------------------------------------------

local inprod = require("inprod_analytic")
local n = inprod.natm*2+inprod.noc*2 + 1
local array = require("array")(n) -- n-array constructor
local format = string.format

--- Restore the dynamical state and the statistics of a previous run
-- from a snapshot file.
-- @function restore
-- @param params parameter table as loaded from "params.lua"
-- @param st statistics object to be initialized; if not given, none will be
-- initialized.
-- @param snapfn Snapshot filename.
-- @return time step at which the snapshot was restored
-- @return state of the system from the snapshot
local function restore(params, st, snapfn)
  snapfn = snapfn or params.i.getoutfn("_snapshot")
  -- Check if snapshot file exists
  local snapf = io.open(snapfn,"r")
  if not snapf or snapf:read("*a")=="" then return end
  snapf:close()
  -- number of lines depends on whether statistics are also in the snapshot file
  local lines = params.i.statistics and 4 or 1
  -- Read in the state from the last lines.
  snapf = io.popen(format("tail -n%d %s",lines,snapfn),"r")
  local t_init = assert(tonumber(snapf:read("*n")), "Could not read t_init from "..snapfn..
  ". Please verify that the previous run was done with the same parameters.")
  io.write(format("* Restored from %s:\n",snapfn))
  io.write(format("* Continuing run of %e time units from time %e\n",params.i.t_run+params.i.t_trans,t_init))
  local initialconditions = {}
  for i=1,n-1 do
    initialconditions[i] = assert(tonumber(snapf:read("*n")), "Could not read initial conditions from "..snapfn..
  ". Please verify that the previous run was done with the same model and parameters.")
  end
  if st and t_init > params.i.t_trans then
  -- Restore the statistics state as well.
    local stats = {mean={[0]=1}, var={[0]=0}}
    for _=1,2 do
      local stat_type = snapf:read(5):match("%S+") or ""
      local stat = assert(stats[stat_type], "Could not read statistics from "..snapfn)
      for i=1,n-1 do
        stat[i] = assert(tonumber(snapf:read("*n")), "Could not read "..stat_type.." from "..snapfn)
      end
    end
    local iter = assert(tonumber(snapf:read("*n")), "Could not read number of iterations from "..snapfn)
    io.write(format("* Restored statistics from %e samples\n",iter))
    st.set(array({data=stats.mean}),array({data=stats.var}),iter)
  end
  snapf:close()
  return t_init, initialconditions
end

-- @export
return restore
