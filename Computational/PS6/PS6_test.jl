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
    exit_pol::Vector{Float64}
    labor_pol::Vector{Float64}
    p::Float64
    μ::Vector{Float64}
    M::Float64
    M_exit::Float64
    L_I::Float64
    L_E::Float64
end

function Init_standard()
    prim = Primitives()
    W = zeros(prim.ns)
    exit_pol = zeros(prim.ns)
    labor_pol = zeros(prim.ns)
    p = 1.0
    μ = zeros(prim.ns)
    M = 1.0
    M_exit = 0.1
    L_I = 100.0
    L_E = 100.0
    res = Res_standard(W, exit_pol, labor_pol, p, μ, M, M_exit, L_I, L_E)
    prim, res
end

function opt_n(prim::Primitives, p::Float64,s::Float64)
    @unpack θ = prim
    n_opt = max(0, (θ*p*s)^(1/(1-θ)))
    n_opt
end

function W_comp(prim::Primitives, res::Res_standard, p::Float64)
    @unpack s_grid,ns, n_grid,β, θ, F, A, cf, ce = prim
    vf_new = zeros(length(s_grid))
    for s_i = 1:ns #iterate over firm size index
        s = s_grid[s_i] #find corresponding firm size
        n = opt_n(prim, p, s) #find optimal labor demand for firm
        res.labor_pol[s_i] = n
        val = p*s*n^θ - n - p*cf
        cont_val = 0.0
        for sp_i = 1:ns #iterate over future states
            cont_val += β*F[s_i, sp_i]*res.W[sp_i]
        end
        if cont_val < 0
            res.exit_pol[s_i] = 1
            vf_new[s_i] = val
        else
            res.exit_pol[s_i] = 0
            vf_new[s_i] = val + cont_val
        end
    end
    vf_new
end

function Bellman(prim::Primitives, res::Res_standard, p::Float64)
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
    res.p = p_guess
    p_guess
end

function invar_mu(prim::Primitives, exit_pf::Vector{Float64}, M::Float64)
    @unpack F, ent_dist = prim
    F_hat = (exit_pf .== 0.0) .* F  #multiply the transition matrix by an indicator for whether the firm stays in the market or not. Only keep transition probabilities for firms which remain
    μ = M*inv(I - F_hat)*F_hat*ent_dist
    μ
end

function labor_demand(prim::Primitives, res::Res_standard, M::Float64)
    @unpack ns, ent_dist, ce = prim
    μ = invar_mu(prim,  res.exit_pol, M)
    ld = 0.0
    ld_I = 0.0
    ld_E = 0.0
    for s_i = 1:ns
        ld_I += res.labor_pol[s_i]*μ[s_i]
        ld_E +=  M*res.labor_pol[s_i]*ent_dist[s_i]
        ld += ld_I + ld_E
    end
    res.L_I = ld_I
    res.L_E = ld_E
    ld
end

function profits(prim::Primitives, res::Res_standard, μ::Vector{Float64}, p::Float64, M::Float64)
    @unpack ns, s_grid, θ, cf, ent_dist =prim
    Π = 0.0
    for s_i=1:ns
        Π += (p*s_grid[s_i]*res.labor_pol[s_i]^θ - res.labor_pol[s_i] - p*cf)*μ[s_i] + M*(p*s_grid[s_i]*res.labor_pol[s_i]^θ - res.labor_pol[s_i] - p*cf)*ent_dist[s_i]
    end
    Π
end
    

function cons_opt(prim::Primitives, Π::Float64)
    @unpack A = prim
    n_supp = 1/A - Π
    n_supp
end

function market_clearing(prim::Primitives, res::Res_standard, p::Float64)
    M_low = 0.0
    M_high = 100.0
    M_guess = 50.0
    tol = 1e-4
    μ = invar_mu(prim, res.exit_pol, M_guess)
    Π = profits(prim, res, μ, p, M_guess)
    L_supp = cons_opt(prim, Π)
    L_d = labor_demand(prim, res, M_guess)
    error = L_d - L_supp
    n = 0
    while abs(error) > tol
        n += 1
        println("Iteration $(n), error: $(error)")
        if error >= 0
            M_high = M_guess
        else
            M_low = M_guess
        end
        M_guess = (M_high + M_low)/2
        μ = invar_mu(prim, res.exit_pol, M_guess)
        Π = profits(prim, res, μ, p, M_guess)
        L_supp = cons_opt(prim, Π)
        L_d = labor_demand(prim, res, M_guess)
        error = L_d - L_supp
    end
    println("Convergence!")
    print("Agg labor: $(L_d)")
    res.M = M_guess
    M_guess
end

function model_solve(prim::Primitives, res::Res_standard)
    res.p = price_finder(prim, res)
    market_clearing(prim, res, res.p)
    res.μ = invar_mu(prim, res.exit_pol, res.M)
end
