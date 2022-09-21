using Parameters, Plots, Interpolations, Optim
include("PS3_interpolation_func.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 0.42)

model_solver(prim, res, param)
@elapsed vf_test = cons_opt(prim, res, param)

F_finder(prim, res, param)

cap_supply(prim, res, param)
lab_supply(prim, res, param)
@unpack val_func, cap_pf, labor_pf = res #unpack the value and policy functions
@unpack A = prim #unpack the capital grid for plotting


#part III
#social security
prim, res = Initialize()
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 0.42) #benchmark model
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("Social security benchmark: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")

#eliminate social security
prim, res = Initialize()
param = init_param(θ=0.0, w=1.05, r = 0.05, b = 0.0, Z = [3.0, 0.5], γ = 0.42) 
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("No social security: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")

#no idiosyncratic risk, with social security
prim, res = Initialize()
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [0.5, 0.5], γ = 0.42)
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("No idiosyncratic risk, with SS: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")

#no idiosyncratic risk, without social security
prim, res = Initialize()
param = init_param(θ=0.0, w=1.05, r = 0.05, b = 0.0, Z = [0.5, 0.5], γ = 0.42)
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("No idiosyncratic risk, without SS: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")

#exogenous labor supply, with social security
prim, res = Initialize()
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0)
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("Exogenous labor supply, with SS: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")

#exogenous labor supply, without social security
prim, res = Initialize()
param = init_param(θ=0.0, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0)
k_e, l_e, w_e, r_e, b_e = kl_search(prim, res, param)
welf_e = welfare(prim, res, param)
cv_e = coeff_of_var(prim, res, param)
println("Exogenous labor supply, without SS: K = $(k_e), L = $(l_e), w = $(w_e), r = $(r_e), b = $(b_e), welfare = $(welf_e), cv = $(cv_e)")
