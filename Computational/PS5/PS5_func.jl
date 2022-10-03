@everywhere begin

#first, construct pi transition matrix
#find sequence of z_t's given transition probability


#guess coefficients and store them
#solve dynamic programming problem (given h_1, which is determined by coefficients)
    #intuitively, agents forecast future prices given the conjectured law of motion for capital
    #have DP problem take the law of motion as a parameter
# construct N x T matrix. This is a panel of N households over time T
    #each row will have the agent's k_t+1 choice in the given state they're in

#using the panel data, run regression of aggregate K_{t+1} on K_{t} for each t. Use this to update the conjectured coefficients for h_1 

#if coefficients are close to previous guess, stop. If not, run again


##### STRUCTS #####
#primitives struct (discount factor, utility parameters, etc)

@with_kw struct Primitives
    β :: Float64 = 0.99 #discount rate
    α :: Float64 = 0.36 #production parameter
    z :: Vector{Float64} = [1.01, 0.99] #technology shock vector
    δ :: Float64 = 0.025 #depreciation rate
    ε_grid :: Vector{Float64} = [0, 1] #employment state vector
    e :: Float64 = 0.3271 #efficiency wage
    N::Int64 = 5000 #number of households in panel
    T::Int64 = 11000 #number of time periods
    Π ::Array{Float64, 2} = pi_init()[1]
    Π_z::Array{Float64, 2} = pi_init()[2]
    Π_st ::Array{Float64,2} = Π^10000
    #Π_z = pi_z_gen(Π)
    Z ::Array{Float64,1} = z_generator(N,T, Π_z)
    E :: Array{Float64, 2} = E_generator()
end

#current results struct (mutable) which stores coefficients, policy functions, and the matrix of observed capital choices for each HH
function Initialize()
    prim = Primitives() #initialize primtiives
    #=val_func_next = zeros(prim.N, prim.na, 2) #initial value function guess for next period
    val_func = zeros(prim.N, prim.na, 2) #initial value function guess; depends on age, asset, and productivity 
    cap_pf = zeros(prim.N, prim.na, 2)#initial policy function guess
    labor_pf = ones(prim.N, prim.na,2)#initial policy function guess index
    Γ = ones(prim.N, prim.na, 2) #create initial guess for F, will change later
    Γ_0 = ones(prim.N, prim.na, 2) #create initial guess for Gamma_0
    μ = mu_finder(prim, 1.0) #find size of cohorts, will sum to one
    K = 3.36#3.6239 #aggregate capital supply guess informed by model output with starting parameters
    L = 0.343#0.3249 #aggregate labor supply guess informed by model output with starting parameters 
    ##Note: not sure if I need these or not ###
    cap_pf_path = zeros(30, prim.N, prim.na, 2) #empty array to store capital pf over time - will change T length later
    lab_pf_path = zeros(30, prim.N, prim.na, 2) #empty array to store labor pf over time - will change T length later
    k_path_guess = zeros(30)
    l_path_guess = zeros(30)
    
    res = Results(val_func_next, val_func, cap_pf, labor_pf, μ, Γ, Γ_0, K, L, cap_pf_path, lab_pf_path, k_path_guess, l_path_guess) #initialize results struct
    prim, res #return structs
    =#
end

function pi_init()
    #construct pi transition matrix
    Π = zeros(4,4)
    Π_z = zeros(2,2)
    durug=1.5
    unempg=0.04
    durgd=8.0
    unempb=0.1
    durbd=8.0
    durub=2.5
    #transition probabilities
    pgg00 = (durug-1)/durug
    pbb00 = (durub-1)/durub
    pbg00 = 1.25*pbb00
    pgb00 = 0.75*pgg00
    pgg01 = (unempg - unempg*pgg00)/(1-unempg)
    pbb01 = (unempb - unempb*pbb00)/(1-unempb)
    pbg01 = (unempb - unempg*pbg00)/(1-unempg)
    pgb01 = (unempg - unempb*pgb00)/(1-unempb)
    pgg = (durgd-1)/durgd
    pgb = 1 - (durbd-1)/durbd
    pgg10 = 1 - (durug-1)/durug
    pbb10 = 1 - (durub-1)/durub
    pbg10 = 1 - 1.25*pbb00
    pgb10 = 1 - 0.75*pgg00
    pgg11 = 1 - (unempg - unempg*pgg00)/(1-unempg)
    pbb11 = 1 - (unempb - unempb*pbb00)/(1-unempb)
    pbg11 = 1 - (unempb - unempg*pbg00)/(1-unempg)
    pgb11 = 1 - (unempg - unempb*pgb00)/(1-unempb)
    pbg = 1 - (durgd-1)/durgd
    pbb = (durbd-1)/durbd
    #matrix
    Π[1,1] = pgg*pgg11
    Π[2,1] = pbg*pbg11
    Π[3,1] = pgg*pgg01
    Π[4,1] = pbg*pbg01
    Π[1,2] = pgb*pgb11
    Π[2,2] = pbb*pbb11
    Π[3,2] = pgb*pgb01
    Π[4,2] = pbb*pbb01
    Π[1,3] = pgg*pgg10
    Π[2,3] = pbg*pbg10
    Π[3,3] = pgg*pgg00
    Π[4,3] = pbg*pbg00
    Π[1,4] = pgb*pgb10
    Π[2,4] = pbb*pbb10
    Π[3,4] = pgb*pgb00
    Π[4,4] = pbb*pbb00

    Π_z[1,1] = pgg
    Π_z[1,2] = pgb
    Π_z[2,1] = pbg
    Π_z[2,2] = pbb
    Π, Π_z
end

function z_generator(N::Int64, T::Int64, Π_z::Array{Float64})
    #given transition probabilities, simulate N paths of shocks of length T
    Z = zeros(N)
    Z[1] = 1
    for t = 2:T
        p = Π_z[Z[t-1],1] #probability of z_b given z_{t-1} in previous period
        z_dist = Bernoulli(p) #probability distribution of whether z_g is observed in next period
        z_next = rand(z_dist, 1) + 1 #draw a random value from dist. Will be 0 if z_g, 1 if z_b. We add 1 to get the correct index for the next value of z (index is 1 for z_g and 2 for z_b)
        Z[t] = z_next
    end
    Z
end

function E_generator(N, T, Z)
    ε_mat = zeros(N,T)
    ε_mat[:,1] = 1
    for t = 2:T, n = 1:N
        z_new = Z[t]
        prob_z_1 = Π[Z[t-1], Z[t], ε_mat[n, t-1], 2] #probability of going to eps = 1 and z_t given we're coming from z_{t-1} and eps_t-1
        #we went from Z[t-1] to Z[t], the prob of going from e to e' is pi_{z z' e e'}
        ε_next = rand(Bernoulli(prob_z_1)) + 1 #adding one to get the index right; when we draw zero, this corresponds to e = 0 which has index 1. When we have e = 1, this corresponds to index of 2
        ε_mat[t, n] = ε_next
    end
end

function coeff_finder()
    pi_init()
    z_generator()

    #initialize pi matrix and zequence of z_t's for each household in panel. 
    #start with initial guess for h_1 coefficients
    #given initial guess, solve HH DP problem
    #given policy functions, 
end
end