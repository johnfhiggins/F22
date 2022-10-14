using Distributions, Random, Plots, Optim, Statistics
Random.seed!(17) 


N = 500

#this function simulates p(x_i | theta) using a vector of iid normals 
function p_S(xi::Float64, theta_g::Vector{Float64}, z_seq::Vector{Float64})
    val = 0.0 #starting value!
    s = length(z_seq) #find the length of the sample of z's 
    d = Normal() #initialize distribution
    for i=1:s #loop over sample of z's
        val += cdf(d, theta_g[1] + (sqrt(theta_g[3])*z_seq[i] + theta_g[2])*xi)/s #iteratively add the cdf divided by the sample size to p_S - this will give us the average probability of observing Y_i given x_i and the guessed theta
    end
    val #return value
end

#this function computes the simulated log likelihood of Y given x and the guessed theta_g. It takes a sequence of normals to be used for the Monte Carlo simulation
function SLL(X::Vector{Float64}, Y, theta_g::Vector{Float64}, z_seq::Vector{Float64})
    val = 0.0 #starting value
    n = length(X) #find length of sample 
    for i in 1:n #iterate over all observations in the sample
        p = p_S(X[i], theta_g, z_seq) #find simulated probability of X given theta_g using Monte Carlo
        arg = p^(Y[i])*(1-p)^(1-Y[i]) #find the input to the log in the likelihood function.  
        if arg > 0 #possibly unnecessary, but I check if the argument of the log will be positive. If not, I make the likelihood incredibly small as a penalty
            val += log(arg)/n #add the likelihood of this observation to the log likelihood (divide by n so we get the average at the end)
        else
            val -= Inf
        end
    end
    val
end
#this function is identical to the above function except that it takes an array of shocks, one for each observation. 
function SLL2(X::Vector{Float64}, Y, theta_g::Vector{Float64}, z_seq::Array{Float64})
    val = 0.0
    n = length(X)
    for i in 1:n
        p = p_S(X[i], theta_g, z_seq[:,i]) #note z_seq[:,i], as opposed to z_seq[i] - we take the shock vector for observation i (which is independent of the shock vector for each other observation)
        arg = p^(Y[i])*(1-p)^(1-Y[i])
        if arg > 0
            val += log(arg)/n
        else
            val -= Inf
        end
    end
    val
end

#this function finds the theta vector which maximizes the log likelihood of the sample
function max_like(X::Vector{Float64},Y, S::Int64, norm::Int64)
    d2 = Normal() #specify distribution of z's
    if norm==0
        z_seq = rand(d2, S)
    else
        z_pre = rand(d2, S)
        z_seq = (z_pre .- mean(z_pre))./std(z_pre)
    end
    sll = optimize(theta_g ->  -SLL(X, Y, theta_g, z_seq), [-Inf, -Inf, 0.00001], [Inf, Inf, 10], [0.5,0.5,0.5]) #find theta_hat which maximizes log likelihood (by minimizing negative log likelihood)
    theta_hat = sll.minimizer #find argmin
    theta_hat_trans = theta_hat 
    #theta_hat_trans[3] = theta_hat[3]^2
    theta_hat_trans
end

function max_like_2(X::Vector{Float64}, Y, S::Int64, norm::Int64)
    d2 = Normal()
    if norm==0
        z_seq = rand(d2, S, length(X))
    else
        z_pre = rand(d2, S, length(X))
        z_seq = zeros( S, length(X))
        for i=1:length(X)
            z_seq[:,i] = (z_pre[:,i] .- mean(z_pre[:,i]))./std(z_pre[:,i])
        end
    end
    sll = optimize(theta_g ->  -SLL2(X, Y, theta_g, z_seq), [-Inf, -Inf, 0.00001], [Inf, Inf, 10], [0.5,0.5,0.5]) #find theta_hat which maximizes log likelihood (by minimizing negative log likelihood)
    theta_hat = sll.minimizer #find argmin
    theta_hat_trans = theta_hat
    #theta_hat_trans[3] = theta_hat[3]^2
    theta_hat_trans
end


function gen_theta_S(N::Int64, M::Int64, S::Int64, norm::Int64)
    θ_S = zeros(M,3)
    for m=1:M
        if mod(m, 100)==0
            print(m/M)
        end
        dist = Normal()
        X_m = rand(dist, N)
        eps = rand(dist,N)
        β = rand(Normal(1,1), N)
        Y_pre = β .* X_m + eps
        Y_m = (Y_pre .> 0.0)
        θ_S[m,:] .= max_like(X_m, Y_m, S, norm)
    end
    θ_S
end

function gen_theta_S2(N::Int64, M::Int64, S::Int64, norm::Int64)
    θ_S = zeros(M,3)
    for m=1:M
        if mod(m, 100)==0
            print(m/M)
        end
        dist = Normal()
        X_m = rand(dist, N)
        eps = rand(dist,N)
        β = rand(Normal(1,1), N)
        Y_pre = β .* X_m + eps
        Y_m = (Y_pre .> 0.0)
        θ_S[m,:] .= max_like_2(X_m, Y_m, S, norm)
    end
    θ_S
end
theta_S1 = gen_theta_S(500, 10,50,0)
theta_S2 = gen_theta_S2(500, 10, 50,0)

theta_S1_n = gen_theta_S(500, 100,25, 1)
theta_S2_n = gen_theta_S2(500, 100, 25,1)

function error_stats(est::Vector{Float64}, true_v::Float64)
    bias = 0.0
    n = length(est)
    for i=1:n
        bias += (est[i] - true_v)/n
    end
    sd = std(est)
    rmse = sqrt(bias^2 + sd^2)
    bias, sd, rmse
end


function result_finder(N::Int64, S::Int64, M::Int64, norm::Int64)
    theta = [0.0, 1.0, 1.0]
    theta_S1 = gen_theta_S(N, M, S, norm)
    theta_S2 = gen_theta_S2(N,M,S, norm)
    if norm==1
        normalized = "normalized"
    else
        normalized = "not normalized"
    end
    bias_1 = zeros(3)
    sd_1 = zeros(3)
    rmse_1 = zeros(3)
    bias_2 = zeros(3)
    sd_2 = zeros(3)
    rmse_2 = zeros(3)
    for i=1:3
        bias_1[i], sd_1[i], rmse_1[i] = error_stats(theta_S1[:,i], theta[i])
        bias_2[i], sd_2[i], rmse_2[i] = error_stats(theta_S2[:,i], theta[i])
    end
    println("Same draw simulated log-likelihood, S = $(S), M = $(M), $(normalized):")
    println("Bias: $(bias_1[1]), $(bias_1[2]), $(bias_1[3])")
    println("Standard deviation: $(sd_1[1]), $(sd_1[2]), $(sd_1[3])")
    println("RMSE: $(rmse_1[1]), $(rmse_1[2]), $(rmse_1[3])")
    println("Independent draw simulated log-likelihood, S = $(S), M = $(M), $(normalized):")
    println("Bias: $(bias_2[1]), $(bias_2[2]), $(bias_2[3])")
    println("Standard deviation: $(sd_2[1]), $(sd_2[2]), $(sd_2[3])")
    println("RMSE: $(rmse_2[1]), $(rmse_2[2]), $(rmse_2[3])")
end


s_grid = [50, 100, 25]
n_grid = [0,1]
for s in s_grid, n in n_grid
    result_finder(500, s, 1000, n)
end