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
    z_grid :: Vector{Float64} = [1.01, 0.99] #technology shock vector
    δ :: Float64 = 0.025 #depreciation rate
    ε_grid :: Vector{Float64} = [0, 1] #employment state vector
    e :: Float64 = 0.3271 #efficiency wage
    N::Int64 = 5000 #number of households in panel
    T::Int64 = 11000 #number of time periods
    Π ::Array{Float64, 2} = pi_init()[1]
    Π_z::Array{Float64, 2} = pi_init()[2]
    Π_st ::Array{Float64,2} = Π^10000
    #Π_z = pi_z_gen(Π)
    Z ::Array{Float64,1} = z_generator(T, Π_z)
    E :: Array{Float64, 2} = E_generator(N, T, Z, Π)
    nk::Int64 = 40
    nkb::Int64 = 30
    k_min ::Float64 = 0.0001
    k_max :: Float64 = 15.0
    kb_min :: Float64 = 11.0
    k_grid :: Vector{Float64} = collect(range(start = k_min, length=nk, stop = k_max))
    k_bar_gr ::Vector{Float64} = collect(range(start = kb_min, length=nkb, stop = k_max))
    u_g :: Float64 = 0.04
    u_b :: Float64 = 0.10
end

#current results struct (mutable) which stores coefficients, policy functions, and the matrix of observed capital choices for each HH

mutable struct Results
    g_coeff :: Vector{Float64} #coefficient guess for good state
    b_coeff :: Vector{Float64} #coefficient guess for bad state
    vf_i :: Array{Float64, 4} #value function array
    pf_i :: SharedArray{Float64, 4} #policy function array
    K_mat :: Array{Float64,2} #matrix containing K_t' for each household and time index
end

function Initialize()
    prim = Primitives() #initialize primitives
    g_coeff  = [0.095, 0.999] #coefficient guess for good state
    b_coeff  = [0.085, 0.999] #coefficient guess for bad state
    vf_i = zeros(prim.nk, 2, prim.nkb,2) #value function array
    pf_i = zeros(prim.nk, 2, prim.nkb, 2) #policy function array
    K_mat = zeros(prim.N, prim.T) #matrix containing K_t' for each household
    res = Results(g_coeff, b_coeff, vf_i, pf_i, K_mat)
    prim, res
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

function z_generator(T::Int64, Π_z::Array{Float64})
    #given transition probabilities, simulate path of productivity shocks of length T
    Z = zeros(Int64, T)
    Z[1] = 1
    for t = 2:T
        p = Π_z[Z[t-1],1] #probability of z_g given z_{t-1} in previous period
        z_dist = Bernoulli(p) #probability distribution of whether z_g is observed in next period
        z_next = Int64(rand(z_dist) + 1) #draw a random value from dist. Will be 0 if z_g, 1 if z_b. We add 1 to get the correct index for the next value of z (index is 1 for z_g and 2 for z_b)
        Z[t] = z_next
    end
    Z
end

function pi_index(z, zp, e, ep)
    ind = []
    if z==1 && zp == 1 && e == 2 && ep==2
        ind = [1,1]   
    elseif z==2 && zp == 1 && e == 2 && ep==2
        ind = [2,1]    
    elseif z==1 && zp == 1 && e == 1 && ep==2
        ind = [3,1]  
    elseif z==2 && zp == 1 && e == 1 && ep==2
        ind = [4,1]   
    elseif z==1 && zp == 2 && e == 2 && ep==2
        ind = [1,2]  
    elseif z==2 && zp == 2 && e == 2 && ep==2
        ind = [2,2]  
    elseif z==1 && zp == 2 && e == 1 && ep==2
        ind = [3,2]
    elseif z==2 && zp == 2 && e == 1 && ep==2
        ind = [4,2]  
    elseif z==1 && zp == 1 && e == 2 && ep==1
        ind = [1,3]    
    elseif z==2 && zp == 1 && e == 2 && ep==1
        ind = [2,3]  
    elseif z==1 && zp == 1 && e == 1 && ep==1
        ind = [3,3]  
    elseif z==2 && zp == 1 && e == 1 && ep==1
        ind = [4,3]  
    elseif z==1 && zp == 2 && e == 2 && ep==1
        ind = [1,4]  
    elseif z==2 && zp == 2 && e == 2 && ep==1
        ind = [2,4]  
    elseif z==1 && zp == 2 && e == 1 && ep==1
        ind = [3,4]  
    elseif z==2 && zp == 2 && e == 1 && ep==1
        ind = [4,4]
    end
    ind
end

function E_generator(N, T, Z, Π)
    ε_mat = zeros(Int64, N,T)
    ε_mat[:,1] .= 2
    for t = 2:T
        if mod(t, 1000) ==0
            print(t)
        end
        for n = 1:N
            ind = pi_index(Z[t-1], Z[t], ε_mat[n, t-1], 2) #find index in Pi matrix corresponding to going from z_t-1 to z_t and eps_t-1 to eps-t = 1
            #print(ind)
            prob_z_1 = Π[ind[1], ind[2]] #find corresponding probability in pi matrix
            ε_dist = Bernoulli(prob_z_1)
            ε_next = Int64(rand(ε_dist) + 1) #draw random variable based on probability in Pi matrix
            ε_mat[n, t] = ε_next
        end
    end
    ε_mat
end

function prices(prim::Primitives, K::Float64, L::Float64, z::Float64)
    @unpack α = prim
    r = α*z*(K/L)^(α-1)
    w = (1-α)*z*(K/L)^(α)
    r, w
end

function utility(c::Float64)
    if c > 0
        u = log(c)
    else
        u = -Inf
    end
    u
end

function ev(prim::Primitives, z_i::Int64, e_i::Int64, vf, kp::Float64, kb::Float64)
    @unpack Π = prim
    @unpack g_coeff, b_coeff = res
    val = 0.0
    for zp_i =1:2
        for ep_i = 1:2
            p = pi_index(z_i, zp_i, e_i, ep_i) 
            pi_transition = Π[p[1], p[2]]
            val += pi_transition*vf[ep_i, zp_i](kp, kb)
        end
    end
    val
end

function Bellman(prim::Primitives, res::Results)
    @unpack k_grid, k_bar_gr, nk, nkb, k_min, k_max, kb_min, β, δ, z_grid, e, u_g, u_b = prim
    @unpack b_coeff, g_coeff = res
    v_next = SharedArray{Float64, 4}(zeros(nk, 2, nkb, 2))

    vf_int_11 = interpolate((k_grid, k_bar_gr),res.vf_i[:, 1, :, 1], Gridded(Linear()))
    vf_int_12 = interpolate((k_grid, k_bar_gr),res.vf_i[:, 1, :, 2], Gridded(Linear()))
    vf_int_21 = interpolate((k_grid, k_bar_gr),res.vf_i[:, 2, :, 1], Gridded(Linear()))
    vf_int_22 = interpolate((k_grid, k_bar_gr),res.vf_i[:, 2, :, 2], Gridded(Linear()))
    vf_int = [[vf_int_11, vf_int_21] [vf_int_12,  vf_int_22]]
    for K_i = 1:nkb
        K = k_bar_gr[K_i]
        for z_i = 1:2
            z = z_grid[z_i]
            #forecast next K_bar
            if z_i == 1
                kb_pred = b_coeff[1] + b_coeff[2]*K #next period's expected K_bar in bad state
                L = e*(1-u_g)
            else
                kb_pred = g_coeff[1] + g_coeff[2]*K #next period's expected K_bar in good state
                L = e*(1-u_b)
            end
            kb_i = findmin(abs.(k_bar_gr .- kb_pred))[2] #find the closest grid point to predicted k_bar
            #kb = k_bar_gr[kb_i] #set k_bar equal to corresponding value
            kb = min(kb_pred, k_max)
            r, w = prices(prim, K, L, z)
            @sync @distributed for k_i = 1:nk
                k = k_grid[k_i]
                for e_i = 1:2
                    budget = r*k + w*(e_i - 1) + (1-δ)*k
                    val_opt = optimize(kp ->  -utility(budget- kp) - β*ev(prim, z_i, e_i, vf_int, kp, kb) , k_min, min(budget, k_max))
                    opt_kp = val_opt.minimizer
                    #of course, we searched over an interval - this means the optimal k_prime may not be on our grid. Here, we find the closest level of capital which lies on the grid
                    closest_kp_i = findmin(abs.(k_grid .- opt_kp ))[2]
                    #find the corresponding level of k in the grid
                    closest_kp = k_grid[closest_kp_i]
                    #update the value function evaluated at the closest point in the grid to the optimal k_prime
                    v_next[k_i, e_i, kb_i, z_i] = -val_opt.minimum #utility(budget - closest_kp) + β*ev(prim, z_i, e_i, vf_int, closest_kp, kb)
                    #update policy function for capital
                    res.pf_i[k_i, e_i, kb_i, z_i] = opt_kp
                end
            end
        end
    end
    v_next
end

function iterate(prim::Primitives, res::Results)
    @unpack nk, k_grid, k_bar_gr = prim
    tol = 1e-4
    n = 0
    error = 100
    while error > tol
        n +=1
        v_next = Bellman(prim, res)
        error = maximum(abs.(res.vf_i - v_next))/maximum(abs.(v_next))
        println("n = $(n), error = $(error)")
        res.vf_i = v_next
    end
    print("Convergence!")
end

function coeff_finder(prim::Primitives)
    

    #initialize pi matrix and zequence of z_t's for each household in panel. 
    #start with initial guess for h_1 coefficients
    #given initial guess, solve HH DP problem
    #given policy functions, 
end
end