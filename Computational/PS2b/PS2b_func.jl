using StatFiles, CSV, Parameters, DataFrames, ForwardDiff, Optim, LatexPrint, BenchmarkTools, Plots, Random, Distributions


@with_kw struct Data 
    df :: DataFrame = DataFrame(load("Mortgage_performance_data.dta"))
    #select the relevant columns from the data and convert into array 
    X ::Array{Float64} = identity.(Array(select(df, [:score_0, :rate_spread, :i_large_loan, :i_medium_loan, :i_refinance, :age_r, :cltv, :dti, :cu, :first_mort_r, :i_FHA, :i_open_year2, :i_open_year3, :i_open_year4, :i_open_year5])))
    #add a column of ones to handle the intercept term
    #X ::Array{Float64}= hcat(ones(size(X_pre,1)), X_pre)
    #construct time-varying characteristics and convert to array
    Z :: Array{Float64} = identity.(Array(select(df,[:score_0, :score_1, :score_2])))
    #extract the column of Y variables, convert into array
    Y ::Array{Float64} = identity.(Array(select(df,[:i_close_0, :i_close_1, :i_close_2])))
    #construct the variable T based on Y
    T :: Array{Float64} = T_maker(Y) 
    #load in our quadrature weights for 1 and 2 dimensional integration
    d1_quad :: Array{Float64} = identity.(Array(DataFrame(CSV.File("KPU_d1_l20.asc", header=false))))
    d2_quad :: Array{Float64} = identity.(Array(DataFrame(CSV.File("KPU_d2_l20.asc", header=false))))
end

#constructs the variable T based on Y
#I know if statements are passe - but I don't care!
function T_maker(Y)
    n = size(Y,1)
    T = zeros(n)
    for i=1:size(Y,1)
        if Y[i, 1] == 1
            T[i] =1
        elseif Y[i,1]==0 && Y[i,2]==1
            T[i] = 2
        elseif Y[i,1]==0 && Y[i,2]==0 && Y[i,3]==1
            T[i]=3
        else 
            T[i]=4
        end
    end
    T
end

function quadrature_d1(data::Data, inside, bound, ρ, σ)
    @unpack d1_quad = data
    dist = Normal(0,1)
    n=size(d1_quad,1)
    val = 0.0
    #transform initial quadrature nodes to lie in (-inf, bound)
    nodes = log.(d1_quad[:,1]) .+ bound
    #find jacobian of transformation 
    jacob = 1 ./(nodes)
    weight = d1_quad[:,2]
    #for i=1:n
    #    val += weight[i]*cdf(dist, -inside - ρ*nodes[i])*(pdf(dist,nodes[i]/σ)/σ)*jacob[i]
    #end
    val = sum(weight .* cdf.(dist, -inside .- ρ.*nodes) .* (pdf.(dist, nodes ./ σ) ./ σ) .* jacob)
    val
end

function quadrature_d2(data::Data, inside, bound_0, bound_1, ρ, σ)
    @unpack d2_quad = data
    dist = Normal(0,1)
    n=size(d2_quad,1)
    val = 0.0
    nodes_0 = log.(d2_quad[:,1]) .+ bound_0
    nodes_1 = log.(d2_quad[:,2]) .+ bound_1
    jacob_0 = 1 ./(nodes_0)
    jacob_1 = 1 ./(nodes_1)
    weights = d2_quad[:,3]
    #for i=1:n
    #    val += weights[i]*cdf(dist, -inside - ρ*nodes_1[i])*(pdf(dist, nodes_1[i]-ρ*nodes_0[i])*pdf(dist, nodes_0[i]/σ)/σ)*jacob_0[i]*jacob_1[i]
    #end
    val = sum( weights .* cdf.(dist, -inside .- ρ .* nodes_1) .* (pdf.(dist, nodes_1 .- ρ .* nodes_0) .* pdf.(dist, nodes_0 ./ σ) ./ σ) .* jacob_0 .* jacob_1)
    val
end

function likelihood(data::Data, T, X, Z, α, β, γ, ρ)
    dist = Normal(0,1)
    σ = 1/(1-ρ)
    ind_0 = α[1] + X*β + Z[:,1]*γ
    ind_1 = α[2] + X*β + Z[:, 2]*γ
    ind_2 = α[3] + X*β + Z[:, 3]*γ
    val = 0.0
    if T ==1 
        val = cdf(dist, ind_0/σ)
    elseif T == 2
        val = quadrature_d1(data, -ind_1, ind_0, ρ, σ)
    elseif T == 3
        val = quadrature_d2(data, -ind_2, ind_0, ind_1, ρ, σ)
    else
        val = quadrature_d2(data, ind_2, ind_0, ind_1, ρ, σ)
    end
    val
end

function log_like(data::Data, α, β, γ, ρ)
    for i=1:length(T)


