local m = require("lgsl.matrix")

local function read_matrix(fn,size)
  local f = io.open(fn,"r")
  if not f then return end
  local mat = m.new()
  for i=1,size do
    for j=1,size do
      local v = f:read("*n")
      mat:set(i,j,v)
    end
  end
end

return read_matrix
