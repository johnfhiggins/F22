#using Parameters, Plots #read in necessary packages

@with_kw struct Primitives
    β::Float64 = 0.99 #discount factor
    θ::Float64 = 0.36 #production
    δ::Float64 = 0.025 #depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 1000, stop= 75.0)) #capital grid
    nk::Int64 = length(k_grid) #number of capital grid states
    Pi::Array{Float64,2} = [0.977  0.023; 0.074  0.926] #initialize productivity transition grid
    prod_mat::Array{Float64,1} = [1.25, 0.2]#productivity matrix
    zk::Int64 = 2#length of productivity vector
end

mutable struct Results
    val_func::Array{Float64,2}  #value function
    pol_func::Array{Float64,2} #policy function
end

function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.nk, 2) #initial value function guess
    pol_func = zeros(prim.nk, 2) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return structs
end

function solve_model(prim::Primitives, res::Results)
    error = 100
    n = 0
    tol = 1e-4
    while error>tol
        n+=1
        v_next=zeros(prim.nk,2) #initialize empty array for updated value function
        for i in 1:2 #iterate over productivity shocks and update value function for both
            v_next[:,i] .= Bellman(prim,res, i) 
        end
        error = maximum(abs.(res.val_func .- v_next)) #reset error term
        res.val_func = v_next #update value function held in results vector
        println(n, "  ",  error) #iteration number and error level

        #if either we reach the max number of iterations or our error is less than the tolerance parameter, stop and output the iteration number and max error
        if mod(n, 5000) == 0 || error <tol
            println(" ")
            println("*************************************************")
            println("AT ITERATION = ", n)
            println("MAX DIFFERENCE = ", error)
            println("*************************************************")
        end
    end
end

#Bellman operator (for specified productivity). Takes the primitives struct, results struct, and the productivity index as arguments
function Bellman(prim::Primitives, res::Results, prod_ind)
    @unpack β,δ,θ,nk, k_grid, Pi, prod_mat,zk = prim
    v_next = zeros(nk,2)

    for i_k = 1:nk #iterate over capital grid
        candidate_max = -1e10 #very low number 
        k = k_grid[i_k]#convert state indices to state values
        prod = prod_mat[prod_ind] #convert current productivity index to productivity value
        budget = prod*k^θ + (1-δ)*k #budget given current state

        for i_kp = 1:nk #loop over choice of k_prime
            kp = k_grid[i_kp]
            c = budget - kp #consumption
            if c>0 #check to make sure that consumption is positive
                val = log(c) #utility from consumption choice
                for i in 1:zk #find expectation of value function in the next period
                    val += β*Pi[prod_ind,i]*res.val_func[i_kp, i]
                end
                if val>candidate_max #check for new max value
                    candidate_max = val
                    res.pol_func[i_k,prod_ind] = kp #update policy function
                end
            end
        end
        v_next[i_k,prod_ind] = candidate_max #update next guess of value function
    end
    v_next[:,prod_ind]
end

#############
