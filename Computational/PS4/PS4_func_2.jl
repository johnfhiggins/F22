@with_kw struct Primitives
    β::Float64 = 0.97 #discount factor
    Π::Array{Float64,2} = [0.9261 0.0739; 0.0189 0.9811 ] #transition matrix
    Z_freq::Vector{Float64} = [0.2037, 0.7963] #ergodic distribution of initial productivity
    η_arr ::Vector{Float64} = [0.59923239,0.63885106 ,0.67846973 ,0.71808840 ,0.75699959 ,0.79591079 ,0.83482198 ,0.87373318 ,0.91264437 ,0.95155556 ,0.99046676 ,0.99872065 ,1.0069745 ,1.0152284 ,1.0234823 ,1.0317362 ,1.0399901 ,1.0482440 ,1.0564979 ,1.0647518 , 1.0730057 ,1.0787834 ,1.0845611 ,1.0903388 , 1.0961165 , 1.1018943 , 1.1076720 , 1.1134497 , 1.1192274 , 1.1250052 , 1.1307829 ,1.1233544 ,1.1159259 ,1.1084974 ,1.1010689 ,1.0936404 ,1.0862119 ,1.0787834 ,1.0713549 ,1.0639264,1.0519200,1.0430000,1.0363000,1.0200000,1.0110000] #age efficiency profile (sorry!)
    n::Float64 = 0.011 #population growth rate
    N::Int64 = 66 #age of death :(
    jr::Int64 = 46 #retirement age
    σ :: Float64 = 2.0 #CRRA parameter
    α::Float64  = 0.36  #capital share in production function
    δ::Float64 = 0.06 #depreciation rate of capital
    na ::Int64 = 1000
    a_max = 100.0
    a_min = 0.0001
    a_range = range(a_min, length=na, stop=a_max)
    A :: Array{Float64} = collect(a_range) #capital grid
end

mutable struct Param
    θ_seq :: Vector{Float64} #sequence of SS taxes 
    T :: Vector{Float64} #number of time periods until new steady state
    
end


mutable struct res_0  
    r_0_SS :: Float64 #r at initial steady state
    r_T_SS :: Float64 #r at new steady state
    w_0_SS ::Float64 #w at initial steady state
    w_T_SS :: Float64 #w at new steady state
    b_0_SS :: Float64 #pension benefits at initial steady state
    Γ_0 :: Array{Float64,3} #steady state distribution at t=0
    Γ_T :: Array{Float64,3} #ss dist at t = T
    K_0_SS :: Float64 #initial SS capital
    L_0_SS :: Float64 #initial SS labor 
    K_T_SS :: Float64 #SS capital at new steady state
    L_T_SS :: Float64 #SS labor at new steady state
    vf_SS :: Array{Float64, 3} #value function at time 0 with Social Security
    vf_0_tran :: Array{Float64,3} #value function at time 0 with transition path
    vf_T :: Array{Float64, 3} #value function at new steady state K_T
    k_path_guess :: Vector{Float64} 
    k_path :: Vector{Float64} #transition path
end

mutable struct res_t
    val_func_next_t :: Array{Float64, 3} #for a given period t, this is next period's value function V_{t+1}
    val_func :: Array{Float64,3} #value function for state t; i.e. V_t
    cap_pf :: Array{Float64,3} #capital policy function struct
    labor_pf :: Array{Float64,3} #labor supply policy function struct
    μ :: Vector{Float64} #size of cohorts by age
    Γ :: Array{Float64,3} # dist of agents over age, prod, and assets from previous period
    Γ_0 :: Array{Float64,3} #initial ss dist of agents over age, prod, and assets

    
    cap_pf_path :: Array{Float64, 4} #capital policy functions indexed by time t
    lab_pf_path :: Array{Float64, 4} #labor policy functions indexed by time t
    k_path_g :: Array{Float64,1} #current guessed capital path
    l_path_g :: Array{Float64,1} #current guessed labor path




mutable struct forward_solve 
    k_seq
    r_seq
    w_seq
    Γ #previous distribution
end

#find the mass of agents in each age based on population growth
function mu_finder(prim::Primitives,mu_1::Float64)
    @unpack n, N = prim
    μ = zeros(N) #empty array
    μ[1] = mu_1
    for i=2:N
        μ[i] = μ[i-1]/(1+n)
    end
    μ = μ/sum(μ)
    μ
end

#based on agent's productivity and chosen level of a_prime, find the optimal labor supply
function opt_l(ap::Float64,a::Float64, prod::Float64, param::Param)
    @unpack 
    @unpack γ, θ_seq = param
    θ = θ_seq[t+1] #based on the time index t (for exercise 2)
    val = max(0,min(1,(γ*(1-θ)*prod*w - (1-γ)*((1+r)*a - ap))/((1-θ)*w*prod))) #the agent's labor supply must be in [0,1]. If it is inside, it will be given by the following function
    val 
end

#based on a consumer's chosen level of a_prime, current assets, and productivity, find their income
function work_budget(ap::Float64, a::Float64, prod::Float64, param::Param, θ::Float64)
    @unpack w, r ,t = param
    θ = θ_seq[t+1]
    #find income of agent
    budg = w*(1-θ)*prod*opt_l(ap, a, prod, param) + (1+r)*a
    budg
end

function utility_w(ap::Float64, a::Float64, prod::Float64, prim::Primitives, param::Param )
    @unpack σ = prim
    @unpack γ = param
    u = 0.0
    if work_budget(ap, a, prod, param) > ap && (opt_l(ap, a, prod, param) < 1 || γ == 1.0)
            u = (((work_budget(ap, a, prod, param)-ap)^(γ) *(1-opt_l(ap, a, prod, param))^(1-γ))^(1-σ))/(1-σ)
    else
        u = -Inf
    end
    u
end

#at a given time t, given aggregate capital and price paths in res_bw, update vfs and pfs for each agent
function cons_opt(t)
    #unpack theta[t+1] based on t

#accept arbitrary theta sequence from param, use it to compute b_t
#given arbitrary capital path, work backwards from T to solve for agents' policy fns
function backward_solver()

#given initial distribution and policy functions from backward_solver, find implied paths of K, r, and w
function forward_solver(t)
    #accept time index as parameter
    #given previous period's asset distribution and policy, find aggregate capital and Labor
    #this will be K_{t}
    #given aggregate capital, and L, solve for implied r and w as r_t, w_t (also solve for b!)
    #do vf iteration given r and w
    #store policy function and find this period's asset distribution
end

#compare the guessed path to the resulting path
function path_compare()

#iterate over guesses of paths until convergence 
function path_finder()
    #init starting guess
    #run path_compare()
    #if error > tol:
        #update guess
        #run path_compare() again
    #end
    #when convergence is reached, the correct paths of K, L, r, and w, b will be in forward_solve. Unpack these and return them (store in trans_path struct)
end

function wf_comp()
    #unpack V_0^SS from res_0
    #also unpack V_0^T from res_0 (the value at new steady state)
    #unpack V_0(j, a, z, theta^T) from res_bw
    #find lambda