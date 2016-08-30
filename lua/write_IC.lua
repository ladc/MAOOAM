-- write_IC.lua
-- (C) 2014-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- The function write_IC writes an initial condition template in IC.lua if the
-- file does not yet exist. It then also updates the file called
-- translation.txt, which contains a translation of the linear array of
-- variables into the specific ocean and atmosphere variable names.
------------------------------------------------------------------------

local rand = math.random

local function tabstr(t)
  local str={}
  for k,v in pairs(t) do str[#str+1] = k.."="..v end
  table.sort(str)
  return str
end

--- Write IC template if necessary, return true if it is written.
-- @param inprod Inner products table; contains info about variables.
-- @param initvec Initial condition vector: OVERRIDES the other options, will
-- be written to IC file, along with the option init_type="read".
-- @param fn String which contains the filename; default: IC.lua
-- @return Nothing. Its effect is to write the new IC file, if needed.
local function write_IC(inprod, initvec, fn)
  fn = fn or "IC.lua"
  local IC = io.open(fn,"r")
  local rewrite = true
  local init_type, size_of_random_noise = "rand", 1e-2
  if IC then 
    IC:close()
    local IC_options = loadfile(fn)()
    init_type, size_of_random_noise = IC_options.init_type or init_type, IC_options.size_of_random_noise or size_of_random_noise
    if init_type=="read" and not initvec then
      rewrite = false
    end
  end
  if init_type=="rand" or initvec=="random" then
    assert(type(size_of_random_noise)=="number",
    "Please specify the table entry size_of_random_noise (number) in "..fn)
    initvec = setmetatable({},{__index=function() return size_of_random_noise*(rand()-0.5)*2 end})
  end
  initvec = initvec or setmetatable({},{__index=function() return 0 end})

  if not rewrite then return end
  IC=io.open(fn,"w")
  assert(IC,fn.." is not present and not writeable.")
  IC:write("-- ",fn,"\n\n-- Initial conditions for MAOOAM.\n\nreturn {\n")
  local i_vec = 1
  for _,var in ipairs{"psi","theta"} do
    IC:write("-- ",var,"\n")
    for _,wavenums in ipairs(inprod.get_awavenum) do
      IC:write(" "..initvec[i_vec].. ", -- ",var,"(",table.concat(tabstr(wavenums),", "),")\n")
      i_vec = i_vec + 1
    end
  end
  for _,var in ipairs{"A","T"} do
    IC:write("-- ",var,"\n")
    for _,wavenums in ipairs(inprod.get_owavenum) do
      IC:write(" "..initvec[i_vec]..", -- ",var,"(",table.concat(tabstr(wavenums),", "),")\n")
      i_vec = i_vec + 1
    end
  end
  IC:write("  init_type = \"",init_type,"\",\n  size_of_random_noise = ",size_of_random_noise,",\n}")
  IC:close()
  print("* Wrote template in "..fn)

  local tab=io.open("translation.txt","w")  -- always update translation.txt
  assert(tab,"Translation file is not present and not writeable.")
  tab:write("-- Translation table for MAOOAM :\n\n")
  for _,var in ipairs{"psi","theta"} do
    tab:write("-- ",var,"\n")
    local i=0
    for _,wavenums in ipairs(inprod.get_awavenum) do
      i=i+1
      tab:write(var,"_",string.format("%d",i),"  <==>  ",table.concat(tabstr(wavenums),", "),"\n")
    end
    tab:write("\n")
  end
  for _,var in ipairs{"A","T"} do
    tab:write("-- ",var,"\n")
    local i=0
    for _,wavenums in ipairs(inprod.get_owavenum) do
      i=i+1
      tab:write(var,"_",string.format("%d",i),"  <==>  ",table.concat(tabstr(wavenums),", "),"\n")
    end
    tab:write("\n")
  end
  tab:write("\n")
  tab:close()
  print("* Wrote translation table in translation.txt.")
end

return write_IC
