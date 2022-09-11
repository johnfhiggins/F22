using Parameters, Plots

@with_kw struct Primitives
    β :: Float64 = 0.9932 #discount rate
    α :: Float64 = 1.5 #crra parameter
    S :: Array{Float64,1} = [1, 0.5] #possible earning levels
    Π :: Array{Float64} = [0.97 0.03; 0.5 0.5] #transition matrix
    A :: Array{Float64, 1} = collect(range(-2.0, length = 1000, stop= 5.0)) #asset grid
end
