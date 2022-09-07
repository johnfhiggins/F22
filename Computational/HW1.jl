using Parameters, Plots #read in necessary packages


@with_kw struct Primitives
    β::Float64 = 0.99 #discount factor
    θ::Float64 = 0.36 #production
    δ::Float64 = 0.025 #depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 1000, stop= 45.0)) #capital grid
    nk::Int64 = length(k_grid) #number of capital grid states
    Pi::Array{Float64,2} = [0.977  0.023; 0.074  0.926] #initialize productivity transition grid
    prod_mat::Array{Float64,1} = [1.25, 0.2]#productivity matrix
    zk::Int64 = 2#length of productivity vector
end

mutable struct Results
    val_func::Array{Float64,2}  #value function
    pol_func::Array{Float64,2} #policy function
end

function solve_model()
    prim = Primitives()
    val_func, pol_func = zeros(prim.nk,2), zeros(prim.nk,2)#initialize policy function vectors
    res = Results(val_func, pol_func)#results

    error = 100
    n = 0
    tol = 1e-4
    while error>tol
        n+=1
        v_next=zeros(prim.nk,2)
        for i in 1:2
            v_next[:,i] .= Bellman(prim,res, i)
        end
        error = maximum(abs.(res.val_func .- v_next)) #reset error term
        res.val_func = v_next
        #update value function held in results vector
        println(n, "  ",  error)

        if mod(n, 5000) == 0 || error <tol
            println(" ")
            println("*************************************************")
            println("AT ITERATION = ", n)
            println("MAX DIFFERENCE = ", error)
            println("*************************************************")
        end
    end
    
    println("Value function converged in ", n, " iterations.")
    println(res.val_func)
    vfplot = Plots.plot(prim.k_grid, res.val_func, title="Value Functions", legend=:bottomright, labels=["High Productivity" "Low Productivity"]) #plot value function
    pfplot = Plots.plot(prim.k_grid, res.pol_func, title="Policy Functions", legend=:bottomright, labels=["High Productivity" "Low Productivity"]) #plot value function
    
    net_k_pf = res.pol_func .- prim.k_grid
    net_pfplot = Plots.plot(prim.k_grid, net_k_pf, title="Change in decision rule, by K", legend=:bottomleft, labels=["High Productivity" "Low Productivity"])
    
    
    display(vfplot)
    display(pfplot)
    display(net_pfplot)
    savefig(vfplot, "vfplot.png")
    savefig(pfplot, "pfplot.png")
    savefig(net_pfplot, "netpfplot.png")
end

#Bellman operator. Note the lack of type declarations inthe function -- another exaple of sub-optimal coding
function Bellman(prim::Primitives, res::Results, prod_ind)
    @unpack β,δ,θ,nk, k_grid, Pi, prod_mat,zk = prim
    v_next = zeros(nk,2)

    for i_k = 1:nk #loop over state space
        candidate_max = -1e10 #something crappy
        k = k_grid[i_k]#convert state indices to state values
        prod = prod_mat[prod_ind] #convert current productivity index to productivity value
        budget = prod*k^θ + (1-δ)*k #budget given current state. Doesn't this look nice?

        for i_kp = 1:nk #loop over choice of k_prime
            kp = k_grid[i_kp]
            c = budget - kp #consumption
            if c>0 #check to make sure that consumption is positive
                val = log(c)
                for i in 1:zk
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

solve_model() #solve the model.





using Profile

Profile.clear()
@profile solve_model()
Profile.print()

#############
