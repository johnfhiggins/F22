using Random, Distributions, LinearAlgebra, Parameters, Plots
Random.seed!(123)

@with_kw struct Primitives
    size_dist = DiscreteUniform(1,10)
    city_N ::Vector{Int64} = rand(size_dist, 1000) #randomly draw number of firms in each city
    active_indicator :: Vector{Int64} = vcat(ones(500), zeros(500)) #indicator for whether there is antitrust enforcement in city
    city_mat ::Array{Int64, 2} = hcat(city_N, active_indicator) #create matrix indicating the number of firms in each city and whether antitrust is active
end

@with_kw struct Parameters_1
    a_0 ::Float64 = 3.0 
    a_1 ::Float64= 0.1
    ν ::Float64= 0.0
    F ::Float64= 1.0
    b_0 ::Float64= 1.0
    b_1 ::Float64= 0.0
    η ::Float64= 0.0
end

@with_kw struct Parameters_2
    a_0 ::Float64 = 5.0 
    a_1 ::Float64= 1.0
    F ::Float64= 1.0
    b_0 ::Float64= 1.0
end

function cournot_results_1(param_1::Parameters_1, N::Int64)
    @unpack a_0, a_1, ν, F, b_0, b_1, η= param_1
    elas = (a_0 + ν + N*(b_0 - η))/(N*(a_0 + ν - (b_0 - η)))
    L_I = 1/(N*elas)
    HHI = 1/N
    L_I, HHI, elas
end


function collusion_results_1(param_1::Parameters_1, N::Int64)
    @unpack a_0, a_1, ν, F, b_0, b_1, η= param_1
    elas = (a_0 + ν + (b_0 - η))/(a_0 + ν - (b_0 - η))
    L_I = 1/(N*elas)
    HHI = 1/N
    L_I, HHI, elas
end

function data_gen_1(prim::Primitives, param_1::Parameters_1)
    @unpack city_mat = prim
    L = zeros(1000)
    HHI = zeros(1000)
    elas = zeros(1000)
    for i=1:1000
        if city_mat[i,2] == 1 || city_mat[i, 1] > 8 #if antitrust enforcement is active or there are too many firms for collusion
            L[i], HHI[i], elas[i] = cournot_results_1(param_1, city_mat[i,1]) #set lerner index, HHI, and elasticity array values equal to the equilibrium values under cournot competition
        else #if collusion happens
            L[i], HHI[i], elas[i] = collusion_results_1(param_1, city_mat[i,1]) #set lerner index, HHI, and elasticity array values equal to the equilibrium values under collusion
        end
    end
    L, HHI, elas
end

function add_noise(data, lower, upper)
    dist = Uniform(lower, upper)
    n = length(data)
    ε = rand(dist, n)
    new_data = data + ε
    new_data
end

function OLS_est(Y_m::Array{Float64}, X_m::Array{Float64})
    β = inv(X_m'X_m)X_m'Y_m
    sse = sum((Y_m - X_m*β).^2)
    covar = inv(X_m'X_m)/sse
    se_0 = sqrt.(covar[1,1])
    se_1 = sqrt.(covar[2,2])
    β[1], β[2], se_0, se_1
end

prim = Primitives()

function find_estimates_1(prim::Primitives)
    param_1 = Parameters_1()
    L_test, HHI_test, elas_test = data_gen_1(prim, param_1)
    Y = add_noise(log.(L_test), -0.05, 0.05)
    X = hcat(ones(1000), log.(HHI_test)) #construct matrix for regression; will consist of a column of ones and a column of the HHI observations

    β_0, β_1, se_beta_0, se_beta_1 =  OLS_est(Y, X)
    println("Pooled estimate:  beta_0 = $(β_0), standard error = $(se_beta_0)")
    println("beta_1 = $(β_1), standard error = $(se_beta_1)")

    X_act = X[1:500, :]
    Y_act = Y[1:500, :]
    β_0_act, β_1_act, se_beta_0_act, se_beta_1_act =  OLS_est(Y_act, X_act)
    println("No collusion estimate: beta_0_act = $(β_0_act), standard error = $(se_beta_0_act)")
    println("beta_1_act = $(β_1_act), standard error = $(se_beta_1_act)")


    X_inact = X[501:1000, :]
    Y_inact = Y[501:1000, :]
    β_0_inact, β_1_inact, se_beta_0_inact, se_beta_1_inact =  OLS_est(Y_inact, X_inact)
    println("Possible collusion estimate: beta_0_inact = $(β_0_inact), standard error = $(se_beta_0_inact)")
    println("beta_1_inact = $(β_1_inact), standard error = $(se_beta_1_inact)")
end

find_estimates_1(prim)

### problem 3

function cournot_results_2(param_2::Parameters_2, ν::Float64, η::Float64)
    @unpack a_0, a_1, F, b_0 = param_2
    N =  (a_0 + ν - (b_0 - η))/sqrt(a_1*F)-1
    elas = (a_0 + ν + N*(b_0 - η))/(N*(a_0 + ν - (b_0 - η)))
    L_I = 1/(N*elas)
    HHI = 1/N
    L_I, HHI, elas
end

#=
function collusion_results_2(param_2::Parameters_2, N::Int64, ν::Float64, η::Float64)
    @unpack a_0, a_1, F, b_0 = param_2
    elas = (a_0 + ν + (b_0 - η))/(a_0 + ν - (b_0 - η))
    L_I = 1/(N*elas)
    HHI = 1/N
    L_I, HHI, elas
end
=#

function data_gen_2(prim::Primitives, param_2::Parameters_2, ν_vec::Vector{Float64}, η_vec::Vector{Float64})
    @unpack city_mat = prim
    L = zeros(1000)
    HHI = zeros(1000)
    elas = zeros(1000)
    for i=1:1000
        L[i], HHI[i], elas[i] = cournot_results_2(param_2, ν_vec[i], η_vec[i])
    end
    L, HHI, elas
end

function find_estimates_2(prim::Primitives,ν_draw::Vector{Float64}, η_draw::Vector{Float64})
    param_2 = Parameters_2()
    L_test, HHI_test, elas_test = data_gen_2(prim,param_2, ν_draw, η_draw)
    Y = add_noise(log.(L_test), -0.05, 0.05)
    X = hcat(ones(1000), log.(HHI_test)) #construct matrix for regression; will consist of a column of ones and a column of the HHI observations

    β_0, β_1, se_beta_0, se_beta_1 =  OLS_est(Y, X)
    println("beta_0 = $(β_0), standard error = $(se_beta_0)")
    println("beta_1 = $(β_1), standard error = $(se_beta_1)")
end

ν_dist = Uniform(-1,1)
ν_draw = rand(ν_dist, 1000)
η_draw = zeros(1000)
find_estimates_2(prim, ν_draw, η_draw)


η_dist = Uniform(-1,1)
η_draw = rand(η_dist, 1000)
ν_draw = zeros(1000)
find_estimates_2(prim, ν_draw, η_draw)


