-- rand.lua
-- (C) 2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Random number generation utility.

--- Set the randomseed to seed.
-- @param seed Number to be passed to math.randomseed. If nil, a number is
-- taken from /dev/urandom or os.time() is used.
-- @param verbose Boolean: if true, print the seed.
local function setrandomseed(seed, verbose)
  if not seed then
    local devurandom = io.open("/dev/urandom","rb")
    if devurandom then
      local b1,b2,b3,b4 = devurandom:read(4):byte(1,4)
      seed = b1 + (256*(b2 + (256 *(b3 + 256 * b4))))
    else
      seed = os.time()
    end
  end
  math.randomseed(seed)
  print("* Set random seed to: ",seed)
end

return {
  setrandomseed = setrandomseed
}
