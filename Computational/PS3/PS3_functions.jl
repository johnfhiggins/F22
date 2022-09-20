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
    A :: Array{Float64} = collect(range(0.0, length=1000, stop=100.0)) #capital grid
    na = Int64 = length(A) #number of grid points in capital grid
end

mutable struct Param
    θ :: Float64 #social security tax rate
    w :: Float64 #wage rate
    r ::Float64 #interest rate
    b :: Float64 #retirement benefits
    Z :: Vector{Float64} #idiosyncratic productivity shock vector
    γ  :: Float64  #weight on consumption
end

mutable struct Results
    val_func :: Array{Float64,3} #value function struct
    cap_pf :: Array{Float64,3} #capital policy function struct
    labor_pf :: Array{Float64,3} #labor supply policy function struct
    μ :: Vector{Float64} #size of cohorts by age
    F :: Array{Float64,3} #ss dist of agents over age, prod, and assets
end

function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = zeros(prim.N, prim.na, 2) #initial value function guess; depends on age, asset, and productivity 
    cap_pf = zeros(prim.N, prim.na, 2) #initial policy function guess
    labor_pf = ones(prim.N, prim.na,2) #initial policy function guess index
    F = ones(prim.N, prim.na, 2) #create initial guess for F, will change later
    μ = mu_finder(prim, 1.0) #find size of cohorts given mu_1 = 1 (we will rescale later)
    res = Results(val_func, cap_pf, labor_pf, μ, F) #initialize results struct
    prim, res #return structs
end

#this function allows us to change the relevant model parameters in a compact, systematic manner
function init_param(;θ::Float64, w::Float64, r::Float64, b::Float64, Z::Vector{Float64}, γ::Float64)
    param = Param(θ, w, r, b, Z, γ)
    param
end


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



#solve consumer's problem given w, r, and b fixed.
function cons_opt(prim::Primitives, res::Results, param::Param)
    @unpack na, A, N, jr, σ, δ, β, Π, η_arr = prim
    θ = param.θ
    w = param.w
    r = param.r
    b = param.b
    Z = param.Z
    γ = param.γ

    for a_i=1:na #eat everything! 
        a = A[a_i]
        c = (1+r)*a + b
        res.val_func[N, a_i, :] .= (c)^((1-σ)*γ)/(1-σ)
        res.cap_pf[N,a_i, :] .= 0.0
    end
    for age in reverse(jr:N-1)#loop backwards over retirement ages
        print(age)
        for a_i=1:na
            budget = (1+r)*A[a_i] + b
            val_prev = -1e10
            cand_max = 0.0
            for ap_i=1:na
                c = budget - A[ap_i] #amount left over for consumption
                if c >= 0
                    val = (c)^((1-σ)*γ)/(1-σ) + β*res.val_func[age+1, ap_i, 1]
                    if val < val_prev #due to concavity we know that once the value starts decreasing, we can stop iterating
                        break
                    else
                        val_prev = val
                        res.cap_pf[age, a_i, :] .= A[ ap_i]
                    end
                end
            end
            res.val_func[age, a_i, :] .= val_prev
            #res.cap_pf[age, a_i, :] .= cand_max
        end
    end

    for age=reverse(1:jr-1)#loop backwards over working ages
        print(age)
        η = η_arr[age]
        for a_i=1:na
            a = A[a_i]
            for z_i = 1:2
                val_prev = -1e10
                cand_max = 0.0
                for ap_i= 1:na
                    ap = A[ap_i]
                    prod = η*Z[z_i]
                    opt_l = max(0,min(1,(γ*(1-θ)*prod*w - (1-γ)*((1+r)*a - ap))/((1-θ)*w*prod)))
                    res.labor_pf[age, a_i, z_i] = opt_l
                    budget = w*(1-θ)*prod*opt_l + (1+r)*a
                    c = budget - ap #amount left over for consumption
                    if c >= 0
                        val = (((c)^γ *(1-opt_l)^(1-γ))^(1-σ))/(1-σ) 
                        for zp_i=1:2
                            val += β*res.val_func[age+1, ap_i, zp_i]*Π[z_i, zp_i]
                        end
                        if val < val_prev
                            break
                        else
                            val_prev = val
                            res.cap_pf[age, a_i, z_i ] = A[ap_i]
                        end
                    end
                end
                res.val_func[age, a_i, z_i] = val_prev
                res.cap_pf[age, a_i, z_i] = cand_max
            end
        end
    end
    res.val_func
end


function F_finder(prim::Primitives, res::Results, param::Param)
    @unpack N, n,na, A, Z_freq, Π = prim
    μ = res.μ
    F = zeros(N, na, 2)
    F[1, 1, 1] = Z_freq[1]*μ[1] #mass of agents of age one with zero assets and shock z_h 
    F[1, 1, 2] = Z_freq[2]*μ[1] #above, but with shock z_l
    for age=2:N
        for a_i=1:na
            a = A[a_i] #the selected capital level corresponding to a_i
            a_i_sel = findall(az -> az == a, res.cap_pf[age-1, :, :]) #find all (a,z) pairs in the previous age which select a
            for az in a_i_sel #iterate over all pairs in previous age which select a_i
                a_pre_i = az[1] #extract the previous a index
                z_pre_i = az[2] #extract previous z index
                #find mass of agents who chose a_i for each productivity type
                F[age, a_i, 1] += F[age-1, a_pre_i, z_pre_i]*Π[z_pre_i, 1]/(1+n) #all agents who chose a_i and got z_h
                F[age, a_i, 2] += F[age-1, a_pre_i, z_pre_i]*Π[z_pre_i, 2]/(1+n) #all agents who chose a_i and got z_l
            end
        end
    end
    res.F = F
end

#function agg_supply(prim::Primitives, res::Results, param::Param)
