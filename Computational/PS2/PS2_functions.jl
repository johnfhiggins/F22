using Parameters, Plots

@with_kw struct Primitives
    β :: Float64 = 0.9932 #discount rate
    α :: Float64 = 1.5 #crra parameter
    S :: Array{Float64,1} = [1, 0.5] #possible earning levels
    Π :: Array{Float64} = [0.97 0.03; 0.5 0.5] #transition matrix
    A :: Array{Float64, 1} = collect(range(-2.0, length = 1000, stop= 5.0)) #asset grid
end

mutable struct Results
    val_func :: Array{Float64} #value function struct
    pol_func :: Array{Float64} #policy function struct
end

function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nk, 3) #initial value function guess
    pol_func = zeros(prim.nk, 3) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return structs
end

function Bellman(prim::Primitives, res::Results, q::Float64) #bellman operator which takes primitives, results struct, employment status, and the market price as given
    @unpack β, α, S, Π, A = prim
     
end