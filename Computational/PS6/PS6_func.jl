@with_kw struct Primitives
    β::Float64 = 0.8
    θ::Float64 = 0.64
    ns::Int64 = 5 #number of grid pts
    s_grid ::Vector{Float64} = [3.98e-4, 3.58, 6.82, 12.18, 18.79]
    n_grid::Vector{Float64} = [1.3e-9, 10.0, 60.0, 300.0, 1000.0]
    F::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;
                           0.1997 0.7201 0.0420 0.0326 0.0056;
                           0.2000 0.2000 0.5555 0.0344 0.0101;
                           0.2000 0.2000 0.2502 0.3397 0.0101; 
                           0.2000 0.2000 0.2500 0.3400 0.0100]
    ent_dist::Vector{Float64} = [0.37, 0.4631, 0.1102, 0.0504, 0.0063]
    A::Float64 = 0.005
    cf::Float64 = 10.0
    ce::Float64 = 5.0
end

mutable struct Res_standard
    W::Vector{Float64}
    pol_func::Vector{Float64}
    p::Float64
end

function Init_standard()
    prim = Primitives()
    W = zeros(prim.ns)
    pol_func = zeros(prim.ns)
    p = 1.0
    res = Res_standard(W, pol_func, p)
    prim, res
end

function W_comp(prim::Primitives, res::Res_standard, p::Float64)
    @unpack s_grid,ns, n_grid,β, θ, F, A, cf, ce = prim
    vf_new = zeros(length(s_grid))
    for s_i = 1:ns #iterate over firm size index
        s = s_grid[s_i] #find corresponding firm size
        n = n_grid[s_i] #find optimal labor demand for firm
        val = p*s*n^θ - n - p*cf
        cont_val = 0.0
        for sp_i = 1:ns #iterate over future states
            cont_val += β*F[s_i, sp_i]*res.W[sp_i]
        end
        if cont_val < 0
            res.pol_func[s_i] = 1
            vf_new[s_i] = val
        else
            res.pol_func[s_i] = 0
            vf_new[s_i] = val + cont_val
        end
    end
    vf_new
end

function Bellman(prim::Primitives, res::Res_standard, p::Float64)
    n = 0
    error = 100
    tol = 1e-9
    while error > tol
        w_next = W_comp(prim, res, p)
        error = maximum(abs.(w_next - res.W))/(maximum(abs.(w_next)))
        res.W = w_next
    end
    println("Convergence!")
end

function free_entry(prim::Primitives, res::Res_standard, p::Float64)
    @unpack ns, ent_dist, ce = prim
    Bellman(prim, res, p)
    ec = 0.0
    for s_i = 1:ns
        ec += ent_dist[s_i]*res.W[s_i]/p
    end
    ec -= ce
    ec
end

function price_finder(prim::Primitives, res::Res_standard)
    p_low = 0.0
    p_high = 1.0
    p_guess = 0.5
    tol = 1e-6
    EC = free_entry(prim, res, p_guess)
    while abs(EC) > tol
        println("Guessed p: $(p_guess), error: $(EC)")
        if EC >= 0
            p_high = p_guess
        else
            p_low = p_guess
        end
        p_guess = (p_low + p_high)/2
        EC = free_entry(prim, res, p_guess)
    end
    p_guess
end
