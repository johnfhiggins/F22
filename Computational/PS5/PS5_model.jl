using Distributed
addprocs(4)
@everywhere using Distributions, Parameters, SharedArrays, Plots, Interpolations, Optim

@everywhere include("PS5_func.jl")

@everywhere prim, res = Initialize()

iterate(prim, res)