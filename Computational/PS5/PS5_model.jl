#using Distributed
#addprocs(4)
using Distributions, Parameters, Plots, Interpolations, Optim

include("PS5_func.jl")

prim, res = Initialize()

iterate(prim, res)
K_test, V_test = panel_sim(prim, res)