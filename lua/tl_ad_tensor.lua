-- tl_ad_tensor.lua
-- (C) 2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
------------------------------------------------------------------------

local tensor = require("tensor")
local coo_to_fficoo = tensor.coo_to_fficoo
local aotensor_tab = require("aotensor")
local aotensor = coo_to_fficoo(tensor.simplify_coo(aotensor_tab))

local jsparse_mul = tensor.jsparse_mul
local sparse_mul3 = tensor.sparse_mul3

--- Compute the Jacobian of MAOOAM in point ystar.
-- @param ystar array with variables in which the jacobian should be evaluated.
-- @return Jacobian in coolist-form (table of tuples {i,j,value})
local function jacobian(ystar)
  return jsparse_mul(aotensor,ystar)
end

--- Compute the TL tensor from the original MAOOAM one
-- @param aotensor_tab model tensor coolist (table form)
-- @return tangent linear model tensor (table form)
local function get_tltensor(aotensor_tab)
  local tltensor = {}
  for _,entry in pairs(aotensor_tab) do
    if entry[2]~=0 then
      tltensor[#tltensor+1] = {entry[1],entry[2],entry[3],entry[4]}
    end
    if entry[3]~=0 then
      tltensor[#tltensor+1] = {entry[1],entry[3],entry[2],entry[4]}
    end
  end
  return tltensor
end

--- Compute the AD tensor from the TL tensor (method 1)
-- @param tltensor_tab model TL tensor coolist (not yet in fficoo form).
-- @return adjoint model tensor (table form)
local function  adtensor_tab_ref(tltensor_tab)
  local adtensor = {}
  for _,entry in pairs(tltensor_tab) do
    adtensor[#adtensor+1] = {entry[2],entry[1],entry[3],entry[4]}
  end
  return adtensor
end

-- Should yield the same tensor as the previous one (test!)

--- Compute the AD tensor from the original MAOOAM one (method 2)
-- @param aotensor_tab model tensor coolist (not yet in fficoo form).
-- @return adjoint model tensor (table form)
local function get_adtensor(aotensor_tab)
  local adtensor = {}
  for _,entry in pairs(aotensor_tab) do
    if entry[1]~=0 then
      if entry[3]~=0 then
        adtensor[#adtensor+1] = {entry[3],entry[1],entry[2],entry[4]}
      end
      if entry[2]~=0 then
        adtensor[#adtensor+1] = {entry[2],entry[1],entry[3],entry[4]}
      end
    end
  end
  return adtensor
end

--- Tendencies for MAOOAM.
-- @param t time
-- @param y array with variables at time t
-- @param buf n-array (buffer) to store derivatives.
local function model(t,y,buf)
  return sparse_mul3(aotensor,y,y,buf)
end

local tltensor = coo_to_fficoo(tensor.simplify_coo(get_tltensor(aotensor_tab)))

--- Tendencies for the TL of MAOOAM in point ystar for perturbation deltay.
-- @param t time
-- @param ystar array with the variables (current point in trajectory)
-- @param deltay array with the perturbation of the variables at time t
-- @param buf n-array (buffer) to store derivatives.
local function tl_traj(t,ystar,deltay,buf)
  return sparse_mul3(tltensor,deltay,ystar,buf)
end

--- Get a function that computes the tendencies for the TL of MAOOAM in point
-- ystar in a form that can be integrated (same function signature as model)
-- @param ystar array with the variables (current point in trajectory)
-- @return function that computes the tendencies for the TL of MAOOAM in point
-- ystar.
local function get_tl(ystar)
  --- Tendencies for the TL of MAOOAM in point ystar.
  -- @function tl
  -- @param t time
  -- @param deltay array with the perturbation of the variables at time t
  -- @param buf n-array (buffer) to store derivatives.
  -- @return buf n-array with tendencies for TL model
  return function(t,deltay,buf)
    return tl_traj(t,ystar,deltay,buf)
  end
end

local adtensor = coo_to_fficoo(tensor.simplify_coo(get_adtensor(aotensor_tab)))

--- Tendencies for the adjoint of MAOOAM in point ystar for perturbation deltay.
-- @param t time
-- @param ystar array with the variables (current point in trajectory)
-- @param deltay array with the perturbation of the variables at time t
-- @param buf n-array (buffer) to store derivatives.
local function ad_traj(t,ystar,deltay,buf)
  return sparse_mul3(adtensor,deltay,ystar,buf)
end

--- Get a function that computes the tendencies for the adjoint of MAOOAM in
-- point ystar in a form that can be integrated (same function signature as
-- model)
-- @param ystar array with the variables (current point in trajectory)
-- @return function that computes the tendencies for the adjoint of MAOOAM in
-- point ystar.
local function get_ad(ystar)
  --- Tendencies for the adjoint of MAOOAM in point ystar.
  -- @function ad
  -- @param t time
  -- @param deltay array with the perturbation of the variables at time t
  -- @param buf n-array (buffer) to store derivatives.
  -- @return buf n-array with tendencies for adjoint model
  return function(t,deltay,buf)
    return ad_traj(t,ystar,deltay,buf)
  end
end

-- @export 
return {
  model = model,
  tl_traj = tl_traj,
  get_tltensor = get_tltensor,
  get_tl = get_tl,
  ad_traj = ad_traj,
  get_adtensor = get_adtensor,
  get_ad = get_ad,
  jacobian = jacobian,
}


