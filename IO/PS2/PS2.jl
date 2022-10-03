using Random, Distributions, LinearAlgebra, Parameters, Plots
Random.seed!(125)

#struct to hold primitives for both parts of the problem
@with_kw struct Primitives
    size_dist = DiscreteUniform(1,10) #uniform distribution for industry size
    city_N ::Vector{Int64} = rand(size_dist, 1000) #randomly draw number of firms in each city
    active_indicator :: Vector{Int64} = vcat(ones(500), zeros(500)) #indicator for whether there is antitrust enforcement in city
    city_mat ::Array{Int64, 2} = hcat(city_N, active_indicator) #create matrix indicating the number of firms in each city and whether antitrust is active
end

#create struct to hold parameters for the first part of the problem
@with_kw struct Parameters_1
    a_0 ::Float64 = 3.0 
    a_1 ::Float64= 0.1
    ν ::Float64= 0.0
    F ::Float64= 1.0
    b_0 ::Float64= 1.0
    b_1 ::Float64= 0.0
    η ::Float64= 0.0
end

#struct to hold parameters for the second pard of the problem 
@with_kw struct Parameters_2
    a_0 ::Float64 = 5.0 
    a_1 ::Float64= 1.0
    F ::Float64= 1.0
    b_0 ::Float64= 1.0
end

#this function takes a struct of parameters and an industry size and outputs the equilibrium Lerner index, HHI, and elasticity of demand assuming cournot competition
function cournot_results_1(param_1::Parameters_1, N::Int64)
    @unpack a_0, a_1, ν, F, b_0, b_1, η= param_1 #unpack relevant parameters from struct
    elas = (a_0 + ν + N*(b_0 - η))/(N*(a_0 + ν - (b_0 - η))) #compute the elasticity of demand in the city
    L_I = 1/(N*elas) #compute lerner index
    HHI = 1/N #compute HHI
    L_I, HHI, elas #return computed variables
end

#this function takes a struct of parameters and an industry size and outputs the equilibrium Lerner index, HHI, and elasticity of demand assuming perfect collusion between firms
function collusion_results_1(param_1::Parameters_1, N::Int64)
    @unpack a_0, a_1, ν, F, b_0, b_1, η= param_1 #unpack relevant parameters
    elas = (a_0 + ν + (b_0 - η))/(a_0 + ν - (b_0 - η)) #compute elasticity of demand
    L_I = 1/(elas) #find lerner index
    HHI = 1/N #find HHI
    L_I, HHI, elas #return quantities
end

function data_gen_1(prim::Primitives, param_1::Parameters_1)
    @unpack city_mat = prim
    #create empty arrays to store equilibrium results
    L = zeros(1000)
    HHI = zeros(1000)
    elas = zeros(1000)
    for i=1:1000 #iterate over cities
        if city_mat[i,2] == 1 || city_mat[i, 1] > 8 #if antitrust enforcement is active or there are too many firms for collusion
            L[i], HHI[i], elas[i] = cournot_results_1(param_1, city_mat[i,1]) #set lerner index, HHI, and elasticity array values equal to the equilibrium values under cournot competition
        else #if collusion happens
            L[i], HHI[i], elas[i] = collusion_results_1(param_1, city_mat[i,1]) #set lerner index, HHI, and elasticity array values equal to the equilibrium values under collusion
        end
    end
    L, HHI, elas
end

#this function takes an input array and bounds for the uniform distribution and adds noise to it 
function add_noise(data, lower, upper)
    dist = Uniform(lower, upper) #create uniform distribution
    n = length(data) #find length of shock vector
    ε = rand(dist, n) #create vector of shocks
    new_data = data + ε #perturb observed data by adding shocks
    new_data #return new data
end

#this function takes input arrays Y_m and X_m and finds the OLS estimate beta as well as the standard error of each coefficient
function OLS_est(Y_m::Array{Float64}, X_m::Array{Float64})
    β = inv(X_m'X_m)X_m'Y_m #find projection coefficient
    sse = sum((Y_m - X_m*β).^2) #find sum of squared errors
    covar = inv(X_m'X_m)/sse #construct covariance matrix
    se_0 = sqrt.(covar[1,1]) #se of beta_0
    se_1 = sqrt.(covar[2,2]) #se of beta_1
    β[1], β[2], se_0, se_1 #return values
end

prim = Primitives()

function find_estimates_1(prim::Primitives)
    param_1 = Parameters_1()
    L_test, HHI_test, elas_test = data_gen_1(prim, param_1) #find equilibrium results for each city 
    Y = add_noise(log.(L_test), -0.05, 0.05) #disturb the observations of the Lerner index
    X = hcat(ones(1000), log.(HHI_test)) #construct matrix for regression; will consist of a column of ones and a column of the HHI observations

    #find OLS estimates for full sample
    β_0, β_1, se_beta_0, se_beta_1 =  OLS_est(Y, X)
    println("Pooled estimate:  beta_0 = $(β_0), standard error = $(se_beta_0)")
    println("beta_1 = $(β_1), standard error = $(se_beta_1)")

    #OLS estimates for sample with active antitrust 
    X_act = X[1:500, :]
    Y_act = Y[1:500, :]
    β_0_act, β_1_act, se_beta_0_act, se_beta_1_act =  OLS_est(Y_act, X_act)
    println("No collusion estimate: beta_0_act = $(β_0_act), standard error = $(se_beta_0_act)")
    println("beta_1_act = $(β_1_act), standard error = $(se_beta_1_act)")

    #OLS estimates for sample with inactive antitrust
    X_inact = X[501:1000, :]
    Y_inact = Y[501:1000, :]
    β_0_inact, β_1_inact, se_beta_0_inact, se_beta_1_inact =  OLS_est(Y_inact, X_inact)
    println("Possible collusion estimate: beta_0_inact = $(β_0_inact), standard error = $(se_beta_0_inact)")
    println("beta_1_inact = $(β_1_inact), standard error = $(se_beta_1_inact)")
end

#find model estimates
find_estimates_1(prim)

### problem 3

#computes cournot equilibrium results given specified parameters nu and eta. Note that N is endogenously determined by nu and eta
function cournot_results_2(param_2::Parameters_2, ν::Float64, η::Float64)
    @unpack a_0, a_1, F, b_0 = param_2
    N =  (a_0 + ν - (b_0 - η))/sqrt(a_1*F)-1 #compute the total number of firms which choose to enter
    elas = (a_0 + ν + N*(b_0 - η))/(N*(a_0 + ν - (b_0 - η))) #compute elasticity
    L_I = 1/(N*elas) #compute lerner
    HHI = 1/N #compute HHI
    L_I, HHI, elas
end

#generates the data; takes primitives and parameter structs, as well as vectors of the nu and eta parameters for each city
function data_gen_2(prim::Primitives, param_2::Parameters_2, ν_vec::Vector{Float64}, η_vec::Vector{Float64})
    @unpack city_mat = prim
    #initialize empty arrays for eq outcomes
    L = zeros(1000)
    HHI = zeros(1000)
    elas = zeros(1000)
    for i=1:1000 #iterate over each city
        L[i], HHI[i], elas[i] = cournot_results_2(param_2, ν_vec[i], η_vec[i]) #for each city, find eq results given specified parameters
    end
    L, HHI, elas
end

#find model estimates given specified draws of nu and eta
function find_estimates_2(prim::Primitives,ν_draw::Vector{Float64}, η_draw::Vector{Float64})
    param_2 = Parameters_2()
    L_test, HHI_test, elas_test = data_gen_2(prim,param_2, ν_draw, η_draw) #find eq results in each city
    Y = add_noise(log.(L_test), -0.05, 0.05) #perturb the observed lerner index
    X = hcat(ones(1000), log.(HHI_test)) #construct matrix for regression; will consist of a column of ones and a column of the HHI observations

    #find and report OLS results
    β_0, β_1, se_beta_0, se_beta_1 =  OLS_est(Y, X)
    println("beta_0 = $(β_0), standard error = $(se_beta_0)")
    println("beta_1 = $(β_1), standard error = $(se_beta_1)")
end

#part b, where nu is uniformly distributed and eta = 0
ν_dist = Uniform(-1,1)
ν_draw = rand(ν_dist, 1000)
η_draw = zeros(1000)
find_estimates_2(prim, ν_draw, η_draw)

#part c, where nu is 0 and eta is uniformly distributed
η_dist = Uniform(-1,1)
η_draw = rand(η_dist, 1000)
ν_draw = zeros(1000)
find_estimates_2(prim, ν_draw, η_draw)


