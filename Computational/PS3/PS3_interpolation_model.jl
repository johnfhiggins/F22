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
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 0.42) #benchmark model
kl_search(prim, res, param)
welfare(prim, res, param)