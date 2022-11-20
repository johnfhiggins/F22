
@with_kw struct Data 
    char_df :: DataFrame = DataFrame(load("Car_demand_characteristics_spec1.dta"))
    IV_df :: DataFrame = DataFrame(load("Car_demand_iv_spec1.dta")) 
    type_df :: DataFrame = DataFrame(load("Simulated_type_distribution.dta")) 
    X :: Array{Float64} = identity.(Array(char_df[!, Not([:Model_id, :Year, :Version_id, :share, :delta_iia])]))
    shares :: Array{Float64} = identity.(Array(char_df[!, "share"]))
    price :: Array{Float64} = identity.(Array(char_df[!, "price"]))
    #select the relevant columns from the data and convert into array 
    char ::Array{Float64} = Array(char_df) 
    δ_iia ::Vector{Float64} = identity.(Array(char_df[!, "delta_iia"]))
    Z ::Array{Float64} = Array(IV_df[!, Not([:Model_id, :Year])]) 
    income :: Array{Float64} = Array(type_df)
    min_year :: Int64 = Int64(minimum(char_df[!, "Year"]))
    max_year :: Int64 = Int64(maximum(char_df[!, "Year"]))
    year_range = collect(min_year:1:max_year)
    mkt_index :: BitMatrix = mkt_indices(char_df, year_range, max_year, min_year)
end

#finds the indices of all products sold in each year, stores as matrix
function mkt_indices(char_df, year_range, max_year, min_year)
    year_indices = zeros(size(char_df,1), max_year- min_year + 1)
    for (index, year) in enumerate(year_range)
        year_indices[:,index] .= (char_df[!, "Year"] .== year)
    end
    val = (year_indices .== 1)
    val
end


function μ_comp(income, price, λ_p)
    μ = zeros(length(price), length(income))
    for i=1:length(income)
        μ[:,i] .= λ_p .* income[i] .* price
    end
    μ
end

function δ_comp(X, price, α, β)
    δ = X * β .+ α .* price
    δ
end 

function ind_share_pred(μ, δ)
    delta_mu = δ .+  μ 
    share_num = exp.(delta_mu)
    share_denom = 1 .+ sum(share_num, dims=1)
    share = share_num ./ share_denom
    share
end

function share_pred(income, μ, δ)
    ind_share = ind_share_pred(μ, δ)
    R = length(income)
    share = sum(ind_share, dims=2) ./ R
    share
end

#to do: jacobian is messed up
function jacobian(σ)
    J_t = size(σ, 1)
    R = size(σ, 2)
    Δ = (1/R).* diagm(ones(J_t)) .* (σ*(1 .- σ)') - (1/R) .* (1 .- diagm(ones(J_t))) .* (σ*σ') 
    Δ
end

function BLP_contraction_mkt( X, income, price, shares, δ_t, α, β, λ_p)
    ε = 1e-12
    μ = μ_comp(income, price, λ_p)
    δ_0 = δ_t
    σ_pred = share_pred(income, μ, δ_0)
    δ_1 = δ_0 + log.(shares) - log.(σ_pred)
    error = maximum(abs.(δ_1 - δ_0))
    δ_0 = δ_1
    ε_newt = 1
    n =1
    while error > ε
        #n +=1
        #println(n, error)
        if error > ε_newt
            σ_pred = share_pred(income, μ, δ_0)
            δ_1 = δ_0 + log.(shares) - log.(σ_pred)
            error = maximum(abs.(δ_1 - δ_0))
            δ_0 = δ_1
        else
            σ_pred = share_pred(income, μ, δ_0)
            Δ = jacobian(σ_pred)
            δ_1 = δ_0 +(inv(Δ ./ σ_pred))*(log.(shares) - log.(σ_pred))
            error = maximum(abs.(δ_1 - δ_0))
            δ_0 = δ_1
        end
    end
    #println("Convergence!")
    δ_1
end 

function BLP_contraction(data::Data, β_iv, λ_p)
    @unpack char, X, shares, δ_iia, income, year_range, mkt_index= data
    δ = zeros(size(char,1))
    α = β_iv[1]
    β = β_iv[2:length(β_iv)]
    for (index, year) in enumerate(year_range) 
        mkt_indicator = mkt_index[:, index]
        X_t = X[mkt_indicator,:]
        price_t = X_t[:,1]
        X_t = X_t[:,2:size(X_t,2)]
        δ_t = δ_iia[mkt_indicator]
        shares_t = shares[mkt_indicator]
        δ[mkt_indicator] = BLP_contraction_mkt(X_t, income, price_t, shares_t, δ_t, α, β, λ_p)
    end
    δ
end

function IV(X, Z, W, δ)
    β =inv((X' * Z) * W * (Z' * X)) * (X' * Z) * W * Z' * δ
    β
end

#function GMM

function λ_grid_search(data, β_IV)
    @unpack X, Z = data
    λ_grid = collect(0.0:0.1:1.0)
    val_array = zeros(length(λ_grid))
    exog = X * β_IV
    W = inv(Z' * Z)
    inner = Z * W * Z'
    for (i,λ_p) in enumerate(λ_grid)
        δ = BLP_contraction(data, β_IV, λ_p)
        ρ = δ - exog
        println(i/10)
        val_array[i] = ρ' * inner * ρ
    end
    val_array
end