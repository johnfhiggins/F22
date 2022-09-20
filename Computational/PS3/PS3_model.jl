using Parameters, Plots
include("PS3_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
param = init_param(θ=0.11, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 0.42)

@profile vf_test = cons_opt(prim, res, param)
Profile.print()
F_finder(prim, res, param)


@unpack val_func, pol_func, pol_func_ind = res #unpack the value and policy functions
@unpack A = prim #unpack the capital grid for plotting
