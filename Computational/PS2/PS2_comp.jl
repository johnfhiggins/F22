using Parameters, Plots
include("PS2_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
@elapsed solve_model(prim,res, 0.99425) #solve the model using functions in the included file

using Profile

Profile.clear()
@profile solve_model(prim, res, 0.99425)
Profile.print()



@unpack val_func, pol_func = res #unpack the value and policy functions
@unpack k_grid = prim #unpack the capital grid for plotting
