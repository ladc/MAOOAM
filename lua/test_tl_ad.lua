-- test_tl_ad.lua
-- (C) 2016 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Tests for the Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
------------------------------------------------------------------------

local inprod = require("inprod_analytic")
local params =require("params")
local n = inprod.natm*2+inprod.noc*2 + 1
local array = require("array")(n) -- n-array constructor
local rk = require("rk2")(n)
local tensor = require("tensor")
local aotensor_tab = tensor.simplify_coo(require("aotensor"))

local tl_ad = require("tl_ad_tensor")

-- Test Taylor property for the Tangent Linear.
-- lim(\lambda->0) M(x+\lambda dx) - M(x) / M'(\lambda dx) = 1

local IC = {[0]=1}

-- Set all values to random.
math.randomseed(1000)

for i=1,n-1 do IC[i]=math.random() end
local y0_IC = array({data=IC})

local ao, tl, ad = tl_ad.model, tl_ad.get_tl(y0_IC), tl_ad.get_ad(y0_IC)

-- Evolve during transient period
local dt = params.i.dt
for t=0,params.i.t_trans,dt do
  rk(y0_IC,ao,t,dt,y0_IC)
end

print("Initial values:\n",y0_IC)

-- Test 1: Taylor test
-- Integrate the original model by one step, integrate with perturbed
-- IC, and test if the difference approximates the TL
local t = 0

for N=0,6 do
  -- Small perturbation.
  local pert = {[0]=0}
  for i=1,n-1 do pert[i]=(2^-N)/math.sqrt(n-1) end
  local dy = array({data=pert})
  print("Perturbation size:\n",dy%dy)

  local y0 = y0_IC*1
  local y0prime = y0 + dy
  local y1 = array()
  local y1prime = array()
  rk(y0,ao,t,dt,y1)
  rk(y0prime,ao,t,dt,y1prime)

  local dy1 = y1prime - y1

  local dy0 = dy*1
  local dy1_tl = array()
  rk(dy0,tl,t,dt,dy1_tl)

  -- Don't forget to set 0'th component to 0...
  dy1[0]=0
  dy1_tl[0]=0

  print("Resulting difference in trajectory: (epsilon ~ 2^-"..N..")")
  print("diff:  ",dy1%dy1)
  print("tl:    ",dy1_tl%dy1_tl)
  print("ratio: ",(dy1%dy1)/(dy1_tl%dy1_tl))
end

-- Test 2: Adjoint Identity: <M(TL).x,y> = <x,M(AD).y>
for j=1,100 do
  -- Any perturbation.
  local pert = {[0]=0}
  for i=1,n-1 do pert[i]=math.random()/math.sqrt(n-1) end
  local dy = array({data=pert})

  for i=1,n-1 do pert[i]=math.random()/math.sqrt(n-1) end
  local dy_bis = array({data=pert})

  -- Calculate M(TL).x in dy1_tl
  local dy0 = dy*1
  local dy1_tl = array() 
  rk(dy0,tl,t,dt,dy1_tl)

  -- Calculate M(AD).x in dy1_ad
  local dy1_ad = array() 
  rk(dy0,ad,t,dt,dy1_ad)

  -- Calculate M(TL).y in dy1_bis_tl
  local dy0_bis= dy_bis*1
  local dy1_bis_tl = array() 
  rk(dy0_bis,tl,t,dt,dy1_bis_tl)

  -- Calculate M(AD).y in dy1_bis_ad
  local dy1_bis_ad = array() 
  rk(dy0_bis,ad,t,dt,dy1_bis_ad)

  -- Calculate norm <M(TL).x,y>
  local norm1 = dy1_tl%dy0_bis
  -- Calculate norm <x,M(AD).y>
  local norm2 = dy0%dy1_bis_ad

  print("<M(TL).x,y> = ",norm1)
  print("<x,M(AD).y> = ",norm2)
  print("Ratio       = ",norm1/norm2)

  -- Calculate norm <M(TL).y,x>
  norm1 = dy1_bis_tl%dy0
  -- Calculate norm <y,M(AD).x>
  norm2 = dy0_bis%dy1_ad

  print("<M(TL).y,x> = ",norm1)
  print("<y,M(AD).x> = ",norm2)
  print("Ratio       = ",norm1/norm2)
end
