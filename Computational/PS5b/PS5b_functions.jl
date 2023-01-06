@with_kw struct Primitives 
    δ :: Float64 = 0.1
    β :: Float64 = 1/1.05
    α :: Float64 = 0.06
    a :: Float64 = 40.0
    b :: Float64 = 10.0
    q_min :: Int64 = 1
    q_max :: Int64 = 10
    cap_grid :: Vector{Float64} = collect(q_min-1:5.0:(q_max-1)*5)
    nq :: Int64 = length(cap_grid)
    state_space :: Vector{Tuple{Int64, Int64}} = vec(collect(Base.product(1:nq, 1:nq)))
    profits :: Vector{Float64} = profits(state_space, cap_grid, a, b)
end

mutable struct Results
    #value function
    val_func :: Vector{Float64}
    #W_func :: Vector{Float64}
    x_k :: Vector{Float64}
    x_k_old :: Vector{Float64}
    Q_mat :: Array{Float64, 2}
end

function Initialize()
    prim = Primitives()
    N = prim.nq * prim.nq
    val_func = zeros(N)
    x_k = zeros(N)
    x_k_old = zeros(N)
    Q_mat = zeros(N, N)
    res = Results(val_func, x_k, x_k_old, Q_mat)
    prim, res
end


#this function computes the probability of transitioning from q to q_prime for each q_prime (given investment decision x_i)
#returns a vector of probabilities
function investment_transition(prim, res, q_ind, x_i)
    @unpack cap_grid, α, δ, q_min, q_max, nq = prim
    probs = zeros(nq)
    if q_ind == 1
        probs[1] = 1/(1 + α*x_i) 
        probs[2] = (α*x_i)/(1 + α*x_i)
    elseif q_ind == nq
        probs[nq] = ((1-δ)+ α*x_i)/(1 + α*x_i)
        probs[nq-1] = δ/(1 + α*x_i)
    else
        probs[q_ind - 1] = δ/(1 + α*x_i)
        probs[q_ind] = (1 - δ)/(1 + α*x_i) + (δ*α*x_i)/(1 + α*x_i)
        probs[q_ind + 1] = ((1-δ)*α*x_i)/(1 + α*x_i)
    end
    probs
end

function continuation_value(prim, res, q_i, q_j, x_j)
    @unpack nq, state_space=prim
    #There is certainly a more elegant way of doing this, but I don't care to find it. This picks the indices of the state space which correspond to capacity q_i for player i. 
    ω_indices = vec(collect(Int64(q_i):nq:(nq*nq)))
    probs = investment_transition(prim, res, q_j, x_j)
    V = res.val_func[ω_indices,:]
    W = sum(probs .* V)
    W
end

function best_response(prim, res, ω_ind, x_j)
    @unpack β, δ, α, state_space, q_max, q_min = prim
    q_i = state_space[ω_ind][1]
    q_j = state_space[ω_ind][2]
    
    #x_j = res.pol_func[ω_j]
    q_plus = min(q_max, q_i + 1)
    q_minus = max(q_min, q_i - 1)
    W_plus = continuation_value(prim, res, q_plus, q_j, x_j)
    W  = continuation_value(prim, res, q_i, q_j, x_j)
    W_minus = continuation_value(prim, res, q_minus, q_j, x_j)
    inside = β*α*((1-δ)*(W_plus - W) + δ*(W - W_minus))
    if inside > 0
        br = max(0, (-1 + sqrt(inside))/(α))
    else
        br = 0
    end
    br
end

function br(q_bar_i, q_j, a)
    q = max(min(q_bar_i, (a - q_j)/2), 0.0)
    q
end

function nash_eq(q_bar_i, q_bar_j, a, b)
    br_i = br(q_bar_i, q_bar_j, a)
    br_j = br(q_bar_j, q_bar_i, a)
    if min(q_bar_i, q_bar_j) >= a/3
        q_i = a/3
        q_j = a/3
    elseif br_i <= q_bar_i && br_j > q_bar_j 
        q_j = br(q_bar_j, br_i, a)
        q_i = br(q_bar_i, q_j, a)
    elseif br_j <= q_bar_j && br_i > q_bar_i
        q_i = br(q_bar_i, br_j, a)
        q_j = br(q_bar_j, q_i, a)
    else
        q_i = q_bar_i
        q_j = q_bar_j
    end
    profit_i = max(0.0, (a/b - (q_i + q_j)/b)*q_i)
    profit_i
end

#gauss-jacobi
function nash_2(q_bar_i, q_bar_j, a, b)
    br_i_0 = br(q_bar_i, a/3, a)
    br_j_0 = br(q_bar_j, a/3, a)
    br_i_1 = br(q_bar_i, br_j_0, a)
    br_j_1 = br(q_bar_j, br_i_0, a)
    error = maximum(abs.([br_i_1- br_i_0, br_j_1 - br_j_0]))
    while error > 1e-6
        br_i_0 = br_i_1
        br_j_0 = br_j_1
        br_i_1 = br(q_bar_i, br_j_0, a)
        br_j_1 = br(q_bar_j, br_i_0, a)
        error = maximum(abs.([br_i_1- br_i_0, br_j_1 - br_j_0]))
    end
    profit = max(0.0, (a/b - (br_i_1 + br_j_1)/b)*br_i_1)
    profit
end


function profits(state_space, cap_grid, a, b)
    nq = length(cap_grid)
    profs = zeros(nq*nq)
    for ω_i=1:(nq*nq)
        q_bar_i= cap_grid[state_space[ω_i][1]]
        q_bar_j = cap_grid[state_space[ω_i][2]]
        profs[ω_i] = nash_2(q_bar_i, q_bar_j, a, b)
    end
    profs
end

function bellman(prim, res)
    @unpack nq, profits, state_space, β = prim
    N = nq * nq
    val = zeros(N)
    for ω_i= 1:N
        q_i =state_space[ω_i][1]
        q_j = state_space[ω_i][2]
        ω_j = findall((state_space .== [(q_j, q_i)]))
        x_j = res.x_k_old[ω_j[1]]
        res.x_k[ω_i] = best_response(prim, res, ω_i, x_j)
        val[ω_i] = profits[ω_i] - res.x_k[ω_i]
        q = investment_transition(prim, res, q_i, res.x_k[ω_i])
        for q_i in 1:nq
            val[ω_i] += β*continuation_value(prim, res, q_i, q_j, x_j)*q[q_i]
        end
    end
    #update policy function
    res.x_k_old = res.x_k
    val
end

function solve_model(prim, res)
    error = 100
    tol = 1e-12
    v_new = bellman(prim, res)
    error = maximum(abs.(v_new - res.val_func))
    n = 0
    while error > tol
        n +=1
        println(n, error)
        res.val_func = v_new
        v_new = bellman(prim, res)
        error = maximum(abs.(v_new - res.val_func))
    end
    print("Convergence!")
    res.val_func = v_new
end


function transition_matrix(prim, res)
    @unpack nq, state_space = prim
    N = nq * nq
    Q = zeros(N, N)
    for ω_i in 1:N
        q_i =state_space[ω_i][1]
        q_j = state_space[ω_i][2]
        ω_j = findall((state_space .== [(q_j, q_i)]))[1]
        x_i = res.x_k[ω_i]
        x_j = res.x_k[ω_j]
        q_i_transition = investment_transition(prim, res, q_i, x_i)
        q_j_transition = investment_transition(prim, res, q_j, x_j)
        for q_i_p = 1:nq
            for q_j_p = 1:nq
                ω_p = findall((state_space .== [(q_i_p, q_j_p)]))[1]
                Q[ω_i,ω_p] = q_i_transition[q_i_p]*q_j_transition[q_j_p]
            end
        end
    end
    Q
end

function recompose(prim, input)
    @unpack state_space, nq=prim
    output = zeros(nq, nq)
    for q_i = 1:nq
        for q_j = 1:nq
            ω = findall((state_space .== [(q_i, q_j)]))[1]
            x_i = input[ω]
            output[q_i, q_j] = x_i
        end
    end
    output
end

function market_sim(ω_init, Q, T, B)
    states = zeros(Int64, B)
    for b=1:B
        ω = ω_init
        for t =1:T
            probs = Q[ω, :]
            cdf = cumsum(probs)/sum(probs)
            shock = rand(Uniform(0,1))
            ω = findfirst((cdf .> shock))
        end
        states[b]= ω
    end
    states
end

function counter(prim, data)
    @unpack state_space, nq = prim
    N = nq*nq
    freq = zeros(N)
    for ω = 1:N
        freq[ω] = sum(data .== ω )/length(data)
    end
    probs = recompose(prim, freq)
    probs
end
