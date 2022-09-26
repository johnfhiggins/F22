using Parameters, Plots, Interpolations, Optim
include("PS4_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
θ_seq_SS = fill(0.11, 31)
param = init_param(θ_seq = θ_seq_SS, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0, t =0, T = 30)
k_0, l_0, w_0, r_0, b_0 = kl_search(prim, res, param) #search for equilibrium quantities/parameters
welf_0 = welfare(prim, res, param) #find total welfare in equilibrium
cv_0 = coeff_of_var(prim, res, param) #compute coefficient of variation
Γ_0 = res.Γ
println("Social security benchmark: K = $(k_0), L = $(l_0), w = $(w_0), r = $(r_0), b = $(b_0), welfare = $(welf_0), cv = $(cv_0)")

prim, res = Initialize()
θ_seq_T = vcat([0.11],  fill(0.0, 30))
param = init_param(θ_seq = θ_seq_T, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0, t =30, T = 30)
k_T, l_T, w_T, r_T, b_T = kl_search(prim, res, param) #search for equilibrium quantities/parameters
welf_T = welfare(prim, res, param) #find total welfare in equilibrium
cv_T = coeff_of_var(prim, res, param) #compute coefficient of variation
Γ_T = res.Γ
println("Social security benchmark: K = $(k_T), L = $(l_T), w = $(w_T), r = $(r_T), b = $(b_T), welfare = $(welf_T), cv = $(cv_T)")

k_path_30 =path_finder(prim, res, param, k_0, 30)

k_path_50 =path_finder(prim, res, 80)
plot(k_path_50)


#NOTE: clean up variable names

k_path_50_test, l_path_50_test = path_finder(prim, res, 80)

r_path, w_path = price_paths(k_path_50, k_path_50_test[2])
rplot = plot(r_path, title="Evolution of interest rate", legend=false)
wplot = plot(w_path, title="Evolution of wage", legend=false)
kplot = plot(k_path_50_test, title="Evolution of capital", legend=false)
savefig(kplot, "kplot.png")
savefig(rplot, "rplot.png")
savefig(wplot, "wplot.png")

k_path_2, l_path_2 = path_finder(prim, res, 110, vcat(fill(0.11, 21) , fill(0.0, 90)))

#T = 30: 10.30622 vs goal of 10.94903
#T = 50: 10.48655 vs goal of 10.49403
#T = 80: 10.49403 vs goal of 10.49403 :) 


#pt 2: T = 50: 10.35369 vs goal of 10.49403
#T = 80: 10.49378 vs 10.49403
#T = 110: 10.49403 :)