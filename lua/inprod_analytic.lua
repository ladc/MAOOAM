-- inprod_analytic.lua
-- (C) 2013-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Inner products between the truncated set of basis functions for the
-- ocean and atmosphere streamfunction fields.
-- These are partly calculated using the analytical expressions from 
-- Cehelsky, P., & Tung, K. K. 87). Theories of multiple equilibria and
-- weather regimes-A critical reexamination. Part II: Baroclinic two-layer
-- models. Journal of the atmospheric sciences, 44(21), 3282-3303, 1987.
------------------------------------------------------------------------

local params = require("params")

local n = params.m.n
local pi, sqrt = math.pi, math.sqrt
-- Convert sparse arrays stored as {coordinate1, coordinate2, ..., value} to sparse tables:
local new_sparse_table = require("tensor").new_sparse_table
local atmos,ocean = {},{}

------------------------------------------------------------------------

-- Build up translation table in both directions
local modeselection = require("modeselection")
table.sort(modeselection.atmos,function(tab1,tab2) return tab1[1]<tab2[1] or (tab1[1]==tab2[1] and tab1[2]<tab2[2]) end)
table.sort(modeselection.ocean,function(tab1,tab2) return tab1[1]<tab2[1] or (tab1[1]==tab2[1] and tab1[2]<tab2[2]) end)

local get_awavenum = {}
for _,wavenumbers in ipairs(modeselection.atmos) do
  if wavenumbers[1]==1 then -- extra mode f_{A(Pi)}
    get_awavenum[#get_awavenum+1]= {P=wavenumbers[2],typ="A",Nx=0,Ny=wavenumbers[2]}
  end
  get_awavenum[#get_awavenum+1]= {M=wavenumbers[1],P=wavenumbers[2],typ="K",Nx=wavenumbers[1], Ny=wavenumbers[2]}
  get_awavenum[#get_awavenum+1]= {H=wavenumbers[1],P=wavenumbers[2],typ="L",Nx=wavenumbers[1], Ny=wavenumbers[2]}
end
local get_a_index={}
for i,triplet in pairs(get_awavenum) do 
  get_a_index[triplet]=i
end

local get_owavenum = {}
for _,wavenumbers in ipairs(modeselection.ocean) do
  get_owavenum[#get_owavenum+1]= {H=wavenumbers[1],P=wavenumbers[2],Nx=wavenumbers[1]/2, Ny=wavenumbers[2]} -- no need of a function type !!
end
local get_o_index={}
for i,doublet in pairs(get_owavenum) do
  get_o_index[doublet]=i
end


------------------------------------------------------------------------

-- Helper functions from Cehelsky & Tung
local function B1(Pi,Pj,Pk) return (Pk+Pj)/Pi end
local function B2(Pi,Pj,Pk) return (Pk-Pj)/Pi end
local function delta(r) return r==0 and 1 or 0 end
local function lambda(r) return r%2==0 and 0 or 1 end
local function S1(Pj,Pk,Mj,Hk) return -(Pk*Mj + Pj*Hk)/2 end
local function S2(Pj,Pk,Mj,Hk) return (Pk*Mj - Pj*Hk)/2 end
local function S3(Pj,Pk,Hj,Hk) return (Pk*Hj + Pj*Hk)/2 end
local function S4(Pj,Pk,Hj,Hk) return (Pk*Hj - Pj*Hk)/2 end

------------------------------------------------------------------------
-- Generate permutations, yielding the in-place permutation and the sign of the
-- permutation.
local function permute(a, n, sign)
  sign = sign or 1
  n = n or #a
  if n == 0 then coroutine.yield(a, sign)
  else
    for i=1,n do
      -- put i-th element as the last one
      a[n], a[i] = a[i], a[n]
      local flip = (n==i) and 1 or -1
      -- generate all permutations of the other elements
      permute(a, n-1, sign*flip)
      -- restore i-th element
      a[n], a[i] = a[i], a[n]
    end
  end
end
local function getpermuter(a)
  return coroutine.wrap(function() permute(a) end)
end

------------------------------------------------------------------------
-- Fill the (sparse) table `tofill` using the permutations of `indices` with
-- `value` multiplied by the sign of the permutation.
local function fill_permutations(indices,sparse_t,value)
  if not value then return end
  local perm = getpermuter(indices)
  local p,sign = perm()
  while p do
    sparse_t:assign(sign*value,unpack(p))
    p,sign = perm()
  end
end

----------------------------------------------------------------
-- Inner products in the equations for the atmosphere
----------------------------------------------------------------

--- `a_{i,j} = (F_i, \nabla^2 F_j)`.
atmos.a = new_sparse_table()
local function calculate_a(i,j)
  local Ti = get_awavenum[i]
  local value = (i==j) and -n^2*Ti.Nx^2 - Ti.Ny^2 or 0
  atmos.a:assign(value,i,j)
end

--- `b_{i,j,k} = (F_i, J(F_j, \nabla^2 F_k))`.
atmos.b = new_sparse_table()
local function calculate_b(i,j,k)
  local value = atmos.a[k][k]*atmos.g[i][j][k]
  if value ~=0 then atmos.b:assign(value,i,j,k) end
end

--- `c_{i,j} = (F_i, \partial_x F_j)`.
-- Beta term for the atmosphere
atmos.c = new_sparse_table()
-- Strict function !! Only accepts KL type.
-- For any other combination, it will not calculate anything
local function calculate_c(i,j)
  local Ti,Tj = get_awavenum[i],get_awavenum[j]
  local value = (Ti.typ=="K") and (Tj.typ=="L") and n*Ti.M*delta(Ti.M-Tj.H)*delta(Ti.P-Tj.P) or nil
  fill_permutations({i,j},atmos.c,value)
end

--- `d_{i,j} = (F_i, \nabla^2 \eta_j)`.
-- Forcing of the ocean on the atmosphere.
atmos.d = new_sparse_table()
-- Atmospheric s tensor and oceanic M tensor must be computed before
--  calling this function !
local function calculate_d(i,j)
  local value = atmos.s[i][j] * ocean.M[j][j]
  atmos.d:assign(value,i,j)
end

--- `g_{i,j,k} = (F_i, J(F_j, F_k))`.
atmos.g = new_sparse_table()
-- This is a strict function: it only accepts AKL KKL and LLL types.
-- For any other combination, it will not calculate anything.
local function calculate_g(i,j,k)
--  io.write("calculate_g(",i,",",j,",",k,")...")
  local value = tonumber(atmos.g[i][j][k])
  if value then return end
  local Ti, Tj, Tk = get_awavenum[i], get_awavenum[j], get_awavenum[k]
  if Ti.typ=="A" and Tj.typ=="K" and Tk.typ=="L" then
    local b1,b2 = B1(Ti.P,Tj.P,Tk.P), B2(Ti.P,Tj.P,Tk.P)
    value = -2*sqrt(2)/pi * Tj.M  * delta(Tj.M-Tk.H) * lambda(Ti.P + Tj.P + Tk.P)
    -- avoid NaNs by checking this first
    if value~=0 then value = value * (b1^2/(b1^2-1) - b2^2/(b2^2-1)) end
  elseif Ti.typ=="K" and Tj.typ=="K" and Tk.typ=="L" then
    local s1,s2 = S1(Tj.P,Tk.P,Tj.M,Tk.H), S2(Tj.P,Tk.P,Tj.M,Tk.H)
    value = s1*(delta(Ti.M-Tk.H-Tj.M)*delta(Ti.P-Tk.P+Tj.P) 
              - delta(Ti.M-Tk.H-Tj.M)*delta(Ti.P+Tk.P-Tj.P)
             + (delta(Tk.H-Tj.M+Ti.M)+delta(Tk.H-Tj.M-Ti.M))*delta(Tk.P+Tj.P-Ti.P))
          + s2*(delta(Ti.M-Tk.H-Tj.M)*delta(Ti.P-Tk.P-Tj.P)
             + (delta(Tk.H-Tj.M-Ti.M) + delta(Ti.M+Tk.H-Tj.M))
              *(delta(Ti.P-Tk.P+Tj.P) - delta(Tk.P-Tj.P+Ti.P)))
  elseif Ti.typ=="L" and Tj.typ=="L" and Tk.typ=="L" then
    local s3,s4 = S3(Tj.P,Tk.P,Tj.H,Tk.H), S4(Tj.P,Tk.P,Tj.H,Tk.H)
    value = s3 * ((delta(Tk.H-Tj.H-Ti.H) - delta(Tk.H-Tj.H+Ti.H))*delta(Tk.P+Tj.P-Ti.P)
		+ delta(Tk.H+Tj.H-Ti.H)*(delta(Tk.P-Tj.P+Ti.P) - delta(Tk.P-Tj.P-Ti.P)))
	 + s4 *((delta(Tk.H+Tj.H-Ti.H)*delta(Tk.P-Tj.P-Ti.P))
	 + (delta(Tk.H-Tj.H+Ti.H) - delta(Tk.H-Tj.H-Ti.H))
	 * (delta(Tk.P-Tj.P-Ti.P) - delta(Tk.P-Tj.P+Ti.P)))
  else return end
  value = value * n
  -- Store the value at all the permuted indices.
  fill_permutations({i,j,k},atmos.g,value)
end

--- `s_{i,j} = (F_i, \eta_j)`.
-- Forcing (thermal) of the ocean on the atmosphere.
atmos.s = new_sparse_table()
local function calculate_s(i,j)
  local Ti,Dj = get_awavenum[i],get_owavenum[j]
  local value
  if Ti.typ=="A" then
    value = lambda(Dj.H)*lambda(Dj.P+Ti.P)
    -- avoiding NaNs by checking this first
    if value~=0 then value = value * 8*sqrt(2)*Dj.P/(pi^2*(Dj.P^2-Ti.P^2)*Dj.H) end
  elseif Ti.typ=="K" then
    value = lambda(2*Ti.M+Dj.H)*delta(Dj.P-Ti.P)
    -- avoiding NaNs by checking this first
    if value~=0 then value = value * 4*Dj.H/(pi*(-4*Ti.M^2+Dj.H^2)) end
  elseif Ti.typ=="L" then
    value = delta(Dj.P-Ti.P) * delta(2*Ti.H-Dj.H)
  else return end
  atmos.s:assign(value,i,j)
end


----------------------------------------------------------------
-- Inner products in the equations for the ocean
----------------------------------------------------------------

--- `K_{i,j} = (\eta_i, \nabla^2 F_j)`.
-- Forcing of the atmosphere on the ocean.
ocean.K = new_sparse_table()
-- atmospheric a and s tensors must be computed before calling
-- this function !
local function calculate_K(i,j)
  local value = atmos.s[j][i] * atmos.a[j][j]
  ocean.K:assign(value,i,j)
end 

--- `M_{i,j} = (eta_i, \nabla^2 \eta_j)`.
-- Forcing of the ocean fields on the ocean.
ocean.M = new_sparse_table()
local function calculate_M(i,j)
  local Di = get_owavenum[i]
  local value = (i==j) and -n^2*Di.Nx^2 - Di.Ny^2 or 0
  ocean.M:assign(value,i,j)
end

--- `N_{i,j} = (eta_i, \partial_x \eta_j)`.
-- Beta term for the ocean
ocean.N = new_sparse_table()
local function calculate_N(i,j)
  local Di,Dj = get_owavenum[i],get_owavenum[j]
  local value = delta(Di.P-Dj.P)*lambda(Di.H+Dj.H)
  value = value~=0 and value*(-2)*Dj.H*Di.H*n/((Dj.H^2-Di.H^2)*pi) or 0
  ocean.N:assign(value,i,j)
end

--- `O_{i,j,k} = (eta_i, J(\eta_j, \eta_k))`.
-- Temperature advection term (passive scalar)
ocean.O = new_sparse_table()
local function calculate_O(i,j,k)
  local value = tonumber(ocean.O[i][j][k])
  if value then return end
  local Di,Dj,Dk = get_owavenum[i],get_owavenum[j],get_owavenum[k]
  local s3,s4 = S3(Dj.P,Dk.P,Dj.H,Dk.H), S4(Dj.P,Dk.P,Dj.H,Dk.H)
  value = s3 * ((delta(Dk.H-Dj.H-Di.H) - delta(Dk.H-Dj.H+Di.H))*delta(Dk.P+Dj.P-Di.P)
                + delta(Dk.H+Dj.H-Di.H)*(delta(Dk.P-Dj.P+Di.P) - delta(Dk.P-Dj.P-Di.P)))
         + s4 *((delta(Dk.H+Dj.H-Di.H)*delta(Dk.P-Dj.P-Di.P))
         + (delta(Dk.H-Dj.H+Di.H) - delta(Dk.H-Dj.H-Di.H))
         * (delta(Dk.P-Dj.P-Di.P) - delta(Dk.P-Dj.P+Di.P)))
  value = value * n/2
  fill_permutations({i,j,k},ocean.O,value)
end
  
--- `C_{i,j,k} = (\eta_i, J(\eta_j,\nabla^2 \eta_k))`.
ocean.C = new_sparse_table()
-- Requires O_{i,j,k} and M_{i,j} to be calculated beforehand.
local function calculate_C(i,j,k)
  local value = ocean.M[k][k]*ocean.O[i][j][k]
  if value ~=0 then ocean.C:assign(value,i,j,k) end
end

--- `W_{i,j} = (\eta_i, F_j)`.
-- Short-wave radiative forcing of the ocean.
ocean.W = new_sparse_table()
-- atmospheric s tensor must be computed before calling
-- this function !
local function calculate_W(i,j)
  local value = atmos.s[j][i]
  ocean.W:assign(value,i,j)
end 

local natm = #get_awavenum
local noc = #get_owavenum

for i=1,natm do
  for j=1,natm do
    calculate_a(i,j)
    calculate_c(i,j)
    for k=1,natm do
      calculate_g(i,j,k)
    end
  end
  for j=1,noc do
    calculate_s(i,j)
  end
end

for i=1,natm do
  for j=1,natm do
    for k=1,natm do
      calculate_b(i,j,k)
    end
  end
end


for i=1,noc do
  for j=1,noc do
    calculate_M(i,j)
    calculate_N(i,j)
    for k=1,noc do
      calculate_O(i,j,k)
    end
  end
end

for i=1,noc do
  for j=1,noc do
    for k=1,noc do
      calculate_C(i,j,k)
    end
  end
end

for i=1,natm do
  for j=1,noc do
    calculate_d(i,j)
    calculate_K(j,i)
    calculate_W(j,i)
  end
end

return {atmos=atmos, ocean=ocean,
        natm=natm, noc=noc, 
        get_awavenum=get_awavenum, get_owavenum=get_owavenum}
