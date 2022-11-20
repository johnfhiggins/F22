using Parameters, CSV, StatFiles, DataFrames, Kronecker, LinearAlgebra, Plots, Optim

include("PS3_func.jl")
data = Data()

β_IV =  IV(data.X, data.Z, inv(data.Z' * data.Z), data.δ_iia)
@elapsed test = BLP_contraction(data, data.shares, β_IV, 0.6)

β_IV_2 = IV(data.X, data.Z, inv(data.Z' * data.Z), test)
test2 =  BLP_contraction(data, data.shares, β_IV_2, 0.6)


test = λ_grid_search(data, inv(data.Z' * data.Z))

opt = optimize(λ -> GMM_objective(data, λ, method="two-step", λ_hat = 0.62), [0.62], BFGS(), Optim.Options(show_trace = true, show_every = 5))