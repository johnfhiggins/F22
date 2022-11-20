using Parameters, CSV, StatFiles, DataFrames, Kronecker, LinearAlgebra, Plots, Optim

include("PS3_func.jl")
data = Data()

β_IV =  IV(data.X, data.Z, inv(data.Z' * data.Z), data.δ_iia)
@elapsed test = BLP_contraction(data, data.shares, β_IV, 0.6)

β_IV_2 = IV(data.X, data.Z, inv(data.Z' * data.Z), test)
test2 =  BLP_contraction(data, data.shares, β_IV_2, 0.6)


test = λ_grid_search(data, inv(data.Z' * data.Z))
λ_ind = argmin(test[:,1])
gmm_val = test[λ_ind, 1] #234.7611
λ_opt = test[λ_ind, 2] # 0.62

opt = optimize(λ -> GMM_objective(data, λ, method="two-step", λ_hat = 0.62), [0.62], LBFGS(), Optim.Options(show_trace = true, show_every = 1))
opt.minimizer
opt.minimum
#
#minimizer: 0.5639
#minimum: 163.7396

