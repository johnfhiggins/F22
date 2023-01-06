using Parameters, Plots, Distributions, GMT

include("PS5b_functions.jl")

prim, res = Initialize()

solve_model(prim, res)

pi_plot = plot(1:10, 1:10, recompose(prim, prim.profits)', st=:surface, camera=(30,30), xlab="q_1_bar", ylab="q_2_bar")
savefig(pi_plot, "piplot.png")

x_surface = recompose(prim, res.x_k)
xplot = plot(0.0:5.0:45.0, 0.0:5.0:45.0, xtest', st=:surface, camera=(345, 15))
savefig(xplot, "xplot.png")

Q = transition_matrix(prim, res)
Q25 = Q^25
q_00 = Q25[1, :]
Q_25 = recompose(prim, q_00)
plot(0.0:5.0:45.0, 0.0:5.0:45.0, Q_25, st=:surface, c=:deep,camera=(30, 15), xlab="q_1_bar", ylab="q_2_bar")

qtest = market_sim(1, Q, 25, 1000)
sim_freq = counter(prim, qtest)
sim_plot = plot(0.0:5.0:45.0, 0.0:5.0:45.0, sim_freq, st=:surface, camera=(30,15), c=:deep,  xlab="q_1_bar", ylab="q_2_bar")
savefig(sim_plot, "simplot.png")