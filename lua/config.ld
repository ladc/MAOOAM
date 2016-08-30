-- To build the MAOOAM documentation, run "ldoc . " in this directory.
file= {"maooam.lua", "tl_ad_tensor.lua", "IC.lua",
  "params.lua","modeselection.lua", "aotensor.lua",
  "inprod_analytic.lua","write_IC.lua","restore.lua", "array.lua",
  "tensor.lua", "rk2.lua", "rk4.lua", "rk4_tl_ad.lua", "stat.lua","rand.lua",
  "gz.lua", }
project = "<img style='margin:auto' alt='MAOOAM' src='https://raw.githubusercontent.com/Climdyn/MAOOAM/master/fortran/doc/MAOOAM_small.png' />"
description = "Lua implementation of the modular arbitrary-order ocean-atmosphere model: MAOOAM"
readme = {"README.md", "TL_AD.md"}
format = "markdown"
dir = "doc"
all=true
