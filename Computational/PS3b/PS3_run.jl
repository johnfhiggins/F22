using Parameters, CSV, StatFiles, DataFrames, Kronecker, LinearAlgebra, Plots

include("PS3_func.jl")
data = Data()

β_IV =  IV(data.X, data.Z, inv(data.Z' * data.Z), data.δ_iia)
@elapsed test = BLP_contraction(data, data.shares, β_IV, 0.6)

β_IV_2 = IV(data.X, data.Z, inv(data.Z' * data.Z), test)
test2 =  BLP_contraction(data, data.shares, β_IV_2, 0.6)