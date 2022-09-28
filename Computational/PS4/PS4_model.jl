using Parameters, Plots, Interpolations, Optim
include("PS4_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
#θ_seq_SS = fill(0.11, 31)
#param = init_param(θ_seq = θ_seq_SS, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0, t =0, T = 30)
#k_0, l_0, w_0, r_0, b_0 = kl_search(prim, res, param) #search for equilibrium quantities/parameters
#welf_0 = welfare(prim, res, param) #find total welfare in equilibrium
#cv_0 = coeff_of_var(prim, res, param) #compute coefficient of variation
#Γ_0 = res.Γ
#println("Social security benchmark: K = $(k_0), L = $(l_0), w = $(w_0), r = $(r_0), b = $(b_0), welfare = $(welf_0), cv = $(cv_0)")

#prim, res = Initialize()
#θ_seq_T = vcat([0.11],  fill(0.0, 30))
#param = init_param(θ_seq = θ_seq_T, w=1.05, r = 0.05, b = 0.2, Z = [3.0, 0.5], γ = 1.0, t =30, T = 30)
#k_T, l_T, w_T, r_T, b_T = kl_search(prim, res, param) #search for equilibrium quantities/parameters
#welf_T = welfare(prim, res, param) #find total welfare in equilibrium
#cv_T = coeff_of_var(prim, res, param) #compute coefficient of variation
#Γ_T = res.Γ
#println("Social security benchmark: K = $(k_T), L = $(l_T), w = $(w_T), r = $(r_T), b = $(b_T), welfare = $(welf_T), cv = $(cv_T)")

#k_path_30 =path_finder(prim, res, param, k_0, 30)

#k_path_50 =path_finder(prim, res, 30)
#plot(k_path_50[1])


#NOTE: clean up variable names

#k_path_50_test, l_path_50_test = path_finder(prim, res, 80)

r_path, w_path = price_paths(k_path_50[1], fill(0.754, 81))
rplot = plot(r_path, title="Evolution of interest rate", legend=false)
wplot = plot(w_path, title="Evolution of wage", legend=false)
kplot = plot(k_path_50[1], title="Evolution of capital", legend=false)
savefig(kplot, "kplot.png")
savefig(rplot, "rplot.png")
savefig(wplot, "wplot.png")

tg = 30
k_path, l_path = path_finder(prim, res, tg)

#T = 30: 10.30622 vs goal of 10.94903
#T = 50: 10.48655 vs goal of 10.49403
#T = 80: 10.49403 vs goal of 10.49403 :) 
k_path = res.k_path_guess


lambd, gamm = equivalent_variation_bench(vcat([0.11],  fill(0.0, tg)), k_path_50[1])
ce_age = ce_avg(prim, lambd, gamm)
ceavg = plot(ce_age, title="Average CE by model age", legend=false)
savefig(ceavg, "ceavg2.png")
prop_favor(prim, lambd .- 1 , gamm)
#
#7.97 percent of the population votes in favor w transition path

pk_path_2, l_path_2 = path_finder(prim, res, 110, vcat(fill(0.11, 21) , fill(0.0, 90)))

r_path_2, w_path_2 = price_paths(k_path_2, l_path_2)
rplot2 = plot(r_path_2, title="Evolution of interest rate (t = 21 phaseout)", legend=false)
wplot2 = plot(w_path_2, title="Evolution of wage (t = 21 phaseout)", legend=false)
kplot2 = plot(k_path_2, title="Evolution of capital (t = 21 phaseout)", legend=false)
savefig(kplot2, "kplot2.png")
savefig(rplot2, "rplot2.png")
savefig(wplot2, "wplot2.png")
#pt 2: T = 50: 10.35369 vs goal of 10.49403
#T = 80: 10.49378 vs 10.49403
#T = 110: 10.49403 :)

