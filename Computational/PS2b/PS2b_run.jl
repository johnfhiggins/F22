using StatFiles, CSV,DataFrames, ForwardDiff, Optim, LatexPrint, BenchmarkTools, Plots
df = DataFrame(load("Mortgage_performance_data.dta"))
