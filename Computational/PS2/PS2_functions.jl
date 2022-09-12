using Parameters, Plots

@with_kw struct Primitives
    β :: Float64 = 0.9932 #discount rate
    α :: Float64 = 1.5 #crra parameter
    S :: Array{Float64,1} = [1, 0.5] #possible earning levels
    Π :: Array{Float64} = [0.97 0.03; 0.5 0.5] #transition matrix
    Π_st :: Array{Float64} = [0.94, 0.06]
    A :: Array{Float64, 1} = collect(range(-2.0, length = 1000, stop= 5.0)) #asset grid
    na :: Int64 = length(A) #number of grid points
end

mutable struct Results
    val_func :: Array{Float64} #value function struct
    pol_func :: Array{Float64} #policy function struct
    pol_func_ind :: Array{Int64} #policy function index struct
end

function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.na, 2) #initial value function guess
    pol_func = zeros(prim.na, 2) #initial policy function guess
    pol_func_ind = ones(Int64, prim.na,2) #initial policy function guess index
    res = Results(val_func, pol_func, pol_func_ind) #initialize results struct
    prim, res #return structs
end

function Bellman(prim::Primitives, res::Results, q::Float64) #bellman operator which takes primitives, results struct, employment status, and the market price as given
    @unpack β, α, S, Π, A, na = prim
    v_next = zeros(na, 2)
    for s_i=1:2
        s = S[s_i] #find current income given employment status index
        for a_i=1:na
            a = A[a_i] #given current asset index, find asset holding in grid
            val_prev = -1e16 #v bad number to start as candidate maximum
            pol_cand = 1 #initial policy choice
            pol_cand_ind = 1 #corresponding index
            budget = s + a #find the agent's budget given income and asset holdings
            for ap_i=1:na #iterate over index of next period asset holdings
                ap = A[ap_i] #find corresponding asset level
                c = budget - q*ap #given the agent's budget and choice of a prime, find resulting consumption
                val = -1e16
                if c > 0
                    val = -2*(1/sqrt(c) -1) #set value equal to the utility of consumption
                    for s_i_n=1:2 #iterate over index of employment status in next period
                        val += β*Π[s_i, s_i_n]*res.val_func[ap_i, s_i_n] #add the probability of state s_i_n times the value of assets ap and status s_i_n, discounted by beta
                    end
                end
                if val > val_prev #check if the value from choosing ap_i is higher than the value from choosing ap_i-1 
                    val_prev = val
                    pol_cand = ap
                    pol_cand_ind = ap_i
                else
                    break
                end
            end
            res.pol_func[a_i, s_i] = pol_cand
            res.pol_func_ind[a_i, s_i] = pol_cand_ind
            v_next[a_i, s_i] = val_prev
        end
    end
    v_next
end

function solve_model(prim::Primitives, res::Results, q::Float64)
    error = 100
    tol = 1e-4
    n=0
    #v_next= zeros(prim.na, 2)
    while error > tol
        n+=1
        v_next = Bellman(prim, res, q)
        error = maximum(abs.(res.val_func .- v_next))/maximum(abs.(v_next))
        res.val_func = v_next #update value function held in results vector
        println(n, "  ",  error) #iteration number and error level
    end
    println("Convergence in $(n) iterations!")
end

function Q_finder(prim::Primitives, res::Results) #find Q
    @unpack A, na, Π = prim
    pf_ind = res.pol_func_ind #policy index matrix
    Q = zeros(na, 2, na, 2)
    for sp_i = 1:2
        for ap_i = 1:na
            ap_choosers = findall(==(ap_i), pf_ind) #find all indices a_i, s_i which lead to choice of ap 
            for x in ap_choosers #iterate over each element
                ai = x[1]
                si = x[2]
                Q[ai, si, ap_i, sp_i] = Π[si, sp_i]
            end
        end
    end
    Q
end
            


function next_mu_finder(prim::Primitives, res::Results, Q::Array{Float64}, μ_0::Array{Float64})
    @unpack A, na, Π, Π_st = prim
    μ_1 = zeros(na, 2)
    for ap_i =1:na
        for sp_i = 1:2
            val = 0.0
            for a_i = 1:na
                for s_i = 1:2
                    val += Q[a_i, s_i, ap_i, sp_i]*μ_0[a_i, s_i]
                end
            end
            μ_1[ap_i, sp_i] = val
        end
    end
    μ_1
end

function invar_dist(prim::Primitives, res::Results)
    @unpack A, na, Π, Π_st = prim
    Q = Q_finder(prim::Primitives, res::Results)
    μ_0 = zeros(na, 2)
    for j=1:na
        for si = 1:2
            μ_0[j, si] = ((A[j] + 2)/(32)) * Π_st[si]
        end
    end
    error = 100
    tol = 1e-4
    n=0
    while error > tol
        n+=1
        μ_1 = next_mu_finder(prim, res, Q, μ_0)
        error = maximum(abs.(μ_1 - μ_0))/maximum(abs.(μ_0)) 
        μ_0 = copy(μ_1)
        if mod(n, 100) == 0
            println(n, ": ", error)
        end
    end
    print("Found it in $(n) tries!")
    μ_0
end
