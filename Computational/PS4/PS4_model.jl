using Parameters, Plots, Interpolations, Optim
include("PS4_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
θ_seq_SS = fill(0.11, 30)
param = init_param(θ_seq = θ_seq_SS, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 0.42, t =0, T = 30)
k_0, l_0, w_0, r_0, b_0 = kl_search(prim, res, param) #search for equilibrium quantities/parameters
welf_0 = welfare(prim, res, param) #find total welfare in equilibrium
cv_0 = coeff_of_var(prim, res, param) #compute coefficient of variation
Γ_0 = res.Γ
println("Social security benchmark: K = $(k_0), L = $(l_0), w = $(w_0), r = $(r_0), b = $(b_0), welfare = $(welf_0), cv = $(cv_0)")



#solve for eq with SS 
#solve for eq without SS
#guess transition path K_t^i

#loop over iterations i
#given guessed K_t^i, for each t=N-1:1 (i.e. work backwards), find V_t given V_{t+1} and thus the optimal policy function
#given policy function in time t, compute aggregate capital
#store as K_t^{i+1}
#once we get to t=0, we can then check |K_t^i - K_t^{i-1}| < tol_i
#then, if this is satisfied, check if |K_T^i - K_T^SS| < tol_K
#if so, we're done! If not, choose a higher T
typeof(θ_seq_SS)