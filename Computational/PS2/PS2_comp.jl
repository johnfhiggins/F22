using Parameters, Plots
include("PS2_functions.jl")

#initialize the Primitives and Results structs, respectively
prim, res = Initialize()
res.q = 0.994273
@elapsed solve_vf(prim,res) #solve the model using functions in the included file

market_clearing(prim, res)


@unpack val_func, pol_func, pol_func_ind = res #unpack the value and policy functions
@unpack A = prim #unpack the capital grid for plotting

#part a
pfplot = plot(A, [pol_func, A], legend =:bottomright, labels=["g(a,e)" "g(a,u)" "45-degree line"], title="Policy function g(a,s) across a,s")
a_diff = pol_func .- A
pfdiffplot = plot(A, [a_diff, zeros(prim.na)], legend=:bottomright, labels=["g(a,e) - a" "g(a,u) - a" "g(a,s) - a = 0"], title="Difference between g(a, s) - a across a,s")
savefig(pfplot, "pfplot.png")
savefig(pfdiffplot, "pfdiffplot.png")

#part b
#eq price found above using market_clearing function is given by 0.99427
wealth_d = wealth_dist(prim, res)
wealth_plot = plot(A, wealth_d, legend=:topright, labels=["Employed" "Unemployed"], title="Wealth distribution by a,s")
savefig(wealth_plot, "wealthplot.png")