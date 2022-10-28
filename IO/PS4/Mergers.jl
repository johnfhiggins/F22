

function own_dum_comp(firm_own)
    own_mat = zeros(24,24)
    for i=1:24 #loop over each product
        for j = 1:24
            own_mat[i, j] = (firm_own[i] .== firm_own[j]) #determine which brands are owned by the same firm as product i, store as 1 if yes and 0 of not
        end
    end
    own_mat
end

function omega1(shares, alpha, own_dum)
    n = 24
    Ω_pre = zeros(n,n)
    for i=1:n
        for j =1:n
            if i == j
                Ω_pre[i,j] = -alpha[i]*shares[i]*(1-shares[i])
            else
                Ω_pre[i,j] = alpha[i]*shares[i]*shares[j]
            end
        end
    end
    Ω = Ω_pre .* own_dum
    Ω
end

function omega1_blp(shares,alpha_i,own_dum)
    n=24
    Ω_pre = zeros(n,n)
    for i=1:n
        for j =1:n
            if i == j
                for k=1:ns
                    Ω_pre[i,j] += -(alpha_i[i,k]*shares[i,k]*(1 .-shares[i,k]))/ns
                end 
            else
                for k=1:ns
                    Ω_pre[i,j] += alpha_i[i,k]*shares[i,k]*shares[j,k]/ns
                end
            end
        end
    end
    Ω = Ω_pre .* own_dum
    Ω
end

function mkup_i( s_jt, alpha_i, own_mat_d) #data is market-level share data; only shares from products in the same market
    mkups = zeros(length(s_jt))
    omeg = omega1(s_jt, alpha_i, own_mat_d)
    mkups = inv(omeg)*s_jt
end

function mkup_i_blp( s_ijt, alpha_i, own_mat_d) 
    n = size(s_ijt,1)
    s_jt = sum(s_ijt, dims=2)/ns
    mkups = zeros(n,1)
    omeg = omega1_blp(s_ijt, alpha_i, own_mat_d)
    mkups = inv(omeg)*s_jt
end

function mkup_finder_pre(alpha_i, own_dum)
    nm = 94
    nbr = 24
    N = nm*nbr
    mkup = zeros(N)
    mc = zeros(N)
    marg = zeros(N)
    s_jt = cereal[:,8]
    prices = cereal[:,9]
    for m = 1:nm
        mkt_range = (m-1)*nbr + 1:m*nbr
        mkup[mkt_range] .= mkup_i(s_jt[mkt_range], alpha_i[mkt_range], own_dum)
        mc[mkt_range] .= prices[mkt_range] - mkup[mkt_range]
        marg[mkt_range] .= mkup[mkt_range]./(prices[mkt_range])
    end
    mkup, mc, marg, prices, s_jt
end

function blp_mkup_pre(alpha_i, own_dum)
    nm = 94
    nbr = 24
    N = nm*nbr
    mkup = zeros(N)
    mc = zeros(N)
    marg = zeros(N)
    delta = meanval(Data(), res)
    expmu = exp.(mufunc(x2,theta2w, Data()))
    s_ijt = ind_sh(exp.(delta), expmu, Data())
    #s_jt = sum(s_ijt, dims=2)/ns
    prices = cereal[:,9]
    for m = 1:nm
        mkt_range = (m-1)*nbr + 1:m*nbr
        mkup[mkt_range] .= mkup_i_blp(s_ijt[mkt_range,:], alpha_i[mkt_range,:], own_dum)
        mc[mkt_range] .= prices[mkt_range] - mkup[mkt_range]
        marg[mkt_range] .= mkup[mkt_range]./(prices[mkt_range])
    end
    mkup, mc, marg, prices, s_ijt
end
#mkblp, mcblp, margblp, pblp, sblp = blp_mkup_pre(a_mat, dum)

function share_finder(expdeltas, cdindex, cdid)
    temp = cumsum(expdeltas, dims=1)
    sum1 = temp[cdindex]
    sum1[2:length(sum1[:,1]),:] = diff(sum1, dims=1)
    denom1 = 1 ./(1 .+ sum1)
    denom = denom1[cdid]
    s_jt = expdeltas .*denom
    s_jt
end


function mkup_finder_merge(own_dum_m, prev_shares, mc, beta, reg_type )
    nm = 94
    nbr = 24
    N = nm*nbr
    mkup = zeros(N)
    marg = zeros(N)
    s_jt = zeros(N)
    prices = zeros(N)
    alph_vec = fill(beta[1], 94*24)
    for m = 1:nm
        mkt_range = (m-1)*nbr + 1:m*nbr
        mkup[mkt_range] .= mkup_i(prev_shares[mkt_range], alph_vec[mkt_range], own_dum_m)
        prices[mkt_range] .= mkup[mkt_range] + mc[mkt_range]
        marg[mkt_range] .= mkup[mkt_range]./(prices[mkt_range])
    end
    #find the shares at the new prices:

    XFE = hcat(prices, brand_dummies)
    if reg_type==1 #if we want regression only on prices
        pred_delta = prices .*beta #beta will be just alpha
    elseif reg_type ==2 #if we want regression on prices with fixed effects #beta will have alpha in the first column and the individual fixed effects in the second
        pred_delta = XFE*beta
    elseif reg_type ==3 #if we want IV with only prices, no fixed effects (beta comes from IV regression)
        pred_delta = prices .*beta
    else #if we want IV with prices and fixed effects (beta comes from IV reg with fixed effects)
        pred_delta = XFE*beta
    end
    #find the numerator of the expected share
    raw_shares = exp.(pred_delta)
    #call the share finder function to find the actual share. It'll be exp(δ_j)/(1 + sum(exp(δ_i))) applied to each element
    s_jt = share_finder(raw_shares, cdindex, cdid)
    #return all of the desired quantities
    mkup, marg, s_jt, prices
end

function mkup_finder_merge_blp(own_dum_m, prev_shares, prev_prices, mc, alpha_i)
    nm = 94
    nbr = 24
    N = nm*nbr
    mkup = zeros(N)
    marg = zeros(N)
    s_ijt = zeros(N)
    prices = zeros(N)
    #delta = meanval(Data(), res, prev_prices)
    #expmu = exp.(mufunc(x2,theta2w, Data()))
    #s_ijt = ind_sh(exp.(delta), expmu, Data())
    for m = 1:nm
        mkt_range = (m-1)*nbr + 1:m*nbr
        mkup[mkt_range] .= mkup_i_blp(prev_shares[mkt_range,:], alpha_i[mkt_range,:], own_dum_m)
        prices[mkt_range] .= mkup[mkt_range] + mc[mkt_range]
        marg[mkt_range] .= mkup[mkt_range]./(prices[mkt_range])
    end
    #find the shares at the new prices:

    new_delta = meanval(Data(), res, prices)
    expmu = exp.(mufunc(x2,theta2w, Data()))
    s_ijt = ind_sh(exp.(new_delta), expmu, Data())

    #find the numerator of the expected share
    #raw_shares = exp.(pred_delta)
    #call the share finder function to find the actual share. It'll be exp(δ_j)/(1 + sum(exp(δ_i))) applied to each element
    #s_jt = share_finder(raw_shares, cdindex, cdid)
    #return all of the desired quantities
    mkup, marg, s_ijt, prices
end


#estimate the four coefficients below:
@unpack y, brand_dummies, ns= Data()

Z = cereal[:,12:21]
p = cereal[:, 9]
X_FE = hcat(p, brand_dummies)

dum = own_dum_comp(cereal[:,3])

function reg_results()
    reg_types = ["OLS", "OLS with FE", "IV", "IV with FE"]
    println("Regression estimates by regression specification:")

    b_ols, seb_ols, R2_ols, adjR2_ols = ols(y, p, 1)
    b_ols_FE, seb_ols_FE, R2_ols_FE, adjR2_ols_FE = ols(y, X_FE, 1)
    b_IV, seb_IV, R2_IV, adjR2_IV = regressIV(y, p, Z, 1)
    b_IV_FE, seb_IV_FE, R2_IV_FE, adjR2_IV_FE = regressIV_FE(y, p, brand_dummies ,Z, 1)
    
    println("OLS: α = $(b_ols[1]), se_α = $(seb_ols), R2 = $(R2_ols), adjR2 =$(adjR2_ols)")
    println("OLS w/ FE: α = $(b_ols_FE[1]), se_α = $(seb_ols_FE[1]), R2 = $(R2_ols_FE), adjR2 =$(adjR2_ols_FE)")
    println("IV: α = $(b_IV[1]), se_α = $(seb_IV), R2 = $(R2_IV), adjR2 =$(adjR2_IV)")
    println("IV w/ FE: α = $(b_IV_FE[1]), se_α = $(seb_IV_FE), R2 = $(R2_IV_FE), adjR2 =$(adjR2_IV_FE)")
end

###merger simulation section!!!###
cereal_pn = copy(cereal) #copy data 
cereal_pn[cereal_pn[:,3] .== 6,3 ] .= 3 #merge post and nabisco by setting each Nabisco firm entry to 3 (the post firm ID)

cereal_gmq = copy(cereal)
cereal_gmq[cereal_gmq[:,3] .== 4,3] .= 2 #merge general mills and quaker by setting each quaker entry to 2 (the General mills firm ID)

#mkup_pre, mc_pre, marg_pre , p_0, s_pre = mkup_finder_pre(fill(-30.6011, 94*24), dum)
#dum_test = own_dum_comp(cereal_pn[1:24,3])
#mkt, mgt, st, pt = mkup_finder_merge(dum_test, s_pre, mc_pre, b_ols, 1) #fill(abs(α_ols), 94*24), 1)
function merger_sim(dataset, reg_type)
    if reg_type == 1 #OLS coefficient vector on just price
        β = ols(y, p, 0)
    elseif reg_type == 2 #OLS coefficient vector for reg on price with brand FE
        β = ols(y, X_FE, 0)
    elseif reg_type == 3 #IV coefficient vector for reg on price
        β = regressIV(y, p, Z, 0)
    else #IV coefficient vector with price and brand FE
        β = regressIV_FE(y, p, brand_dummies, Z, 0)
    end
    α = β[1]
    #find initial mkups, marginal costs, margins, prices, and shares given regression type
    mkup_pre, mc_pre, marg_pre , p_0, s_pre = mkup_finder_pre(fill(α, 94*24), dum)
    #println(mean(p_0))
    p_pre = p_0 
    own_dum_m = own_dum_comp(dataset[1:24,3])
    mk_t, marg_t, s_t, p_1 = mkup_finder_merge(own_dum_m, s_pre, mc_pre, β, reg_type)
    tol = 1e-9
    error = maximum(abs.(p_1 - p_0))
    n = 0
    maxiter = 200
    while error > tol && n < maxiter
        n +=1
        p_0 = p_1
        mk_t, marg_t, s_t, p_1 = mkup_finder_merge(own_dum_m, s_t, mc_pre, β, reg_type)
        error = maximum(abs.(p_1 - p_0))
        #println(mean(p_1), ", ", mean(marg_t), ", ", mean(s_t))
    end
    #print("Convergence!")
    mkup_pre, mc_pre, marg_pre, p_pre, s_pre, mk_t, marg_t, p_1, s_t
end

function logit_merger_results()
    merger_types = ["Post-Nabisco", "GM-Quaker"]
    merger_datasets = [cereal_pn, cereal_gmq]
    reg_types = ["OLS", "OLS with FE", "IV", "IV with FE"]
    for j=1:2
        print("$(merger_types[j]) merger results:")
        for i=1:4
            mkup_pre, mc_pre, marg_pre, p_0, s_pre, mkup_m, marg_m, p_m, s_m = merger_sim(merger_datasets[j], i)
            println(reg_types[i])
            println("***************************")
            println("Mean initial price: $(mean(p_0)), median initial price: $(median(p_0)), standard deviation: $(std(p_0))")
            println("Mean initial markup: $(mean(mkup_pre)), median initial markup: $(median(mkup_pre)), standard deviation: $(std(mkup_pre))")
            println("Mean initial margins: $(mean(marg_pre)), median initial margins: $(median(marg_pre)), standard deviation: $(std(marg_pre))")
            println("Mean marginal costs: $(mean(mc_pre)), median marginal costs: $(median(mc_pre)), standard deviation: $(std(mc_pre))")
            println("Mean initial market share: $(mean(s_pre)), median initial market share: $(median(s_pre)), standard deviation: $(std(s_pre))")
            println("***************************")
            println("Mean post-merger price change: $(mean(p_m)-mean(p_0)), median post-merger price change: $(median(p_m)-median(p_0)), standard deviation change: $(std(p_m))")
            println("Mean post-merger markup change: $(mean(mkup_m) - mean(mkup_pre)), median post-merger markup change: $(median(mkup_m) - median(mkup_pre)), standard deviation change: $(std(mkup_m))")
            println("Mean post-merger margin change: $(mean(marg_m)-mean(marg_pre)), median post-merger margin change: $(median(marg_m) - median(marg_pre)), standard deviation: $(std(marg_m))")
            println("Mean post-merger market share change: $(mean(s_m)-mean(s_pre)), median post-merger market share change: $(median(s_m)-median(s_pre)), standard deviation: $(std(s_m))")
        end
    end
end

function merger_sim_blp(dataset)
    ##
    #find alpha_i
    a_mat = alpha_finder(Data(), res, -32.4374)

    #find initial mkups, marginal costs, margins, prices, and shares 
    mkblp, mcblp, margblp, p_pre, sblp = blp_mkup_pre(a_mat, dum)

    p_0 = p_pre
    own_dum_m = own_dum_comp(dataset[1:24,3])

    mkblpm, margblpm, sblpm, p_1 = mkup_finder_merge_blp(own_dum_m, sblp, p_pre, mcblp, a_mat)
    tol = 1e-9
    error = maximum(abs.(p_1 - p_0))
    n = 0
    maxiter = 200
    while error > tol && n < maxiter
        n +=1
        p_0 = p_1
        mkblpm, margblpm, sblpm, p_1 = mkup_finder_merge_blp(own_dum_m, sblp, p_pre, mcblp, a_mat)
        error = maximum(abs.(p_1 - p_0))
        #println(mean(p_1), ", ", mean(marg_t), ", ", mean(s_t))
    end
    #print("Convergence!")
    mkblp, mcblp, margblp, p_pre, sblp, mkblpm, margblpm, p_1, sblpm
end

function blp_merger_results()
    for j=1:2
        merger_types = ["Post-Nabisco", "GM-Quaker"]
        merger_datasets = [cereal_pn, cereal_gmq]
        println("$(merger_types[j]) merger results, BLP:")
            mkup_pre, mc_pre, marg_pre, p_0, s_pre, mkup_m, marg_m, p_m, s_m = merger_sim_blp(merger_datasets[j])
            println("***************************")
            println("Mean initial price: $(mean(p_0)), median initial price: $(median(p_0)), standard deviation: $(std(p_0))")
            println("Mean initial mc: $(mean(mc_pre)), median initial mc: $(median(mc_pre)), standard deviation: $(std(mc_pre))")
            println("Mean initial markup: $(mean(mkup_pre)), median initial markup: $(median(mkup_pre)), standard deviation: $(std(mkup_pre))")
            println("Mean initial margins: $(mean(marg_pre)), median initial margins: $(median(marg_pre)), standard deviation: $(std(marg_pre))")
            println("Mean initial market share: $(mean(s_pre)), median initial market share: $(median(s_pre)), standard deviation: $(std(s_pre))")
            println("***************************")
            println("Mean post-merger price change: $(mean(p_m)-mean(p_0)), median post-merger price change: $(median(p_m)-median(p_0)), standard deviation: $(std(p_m))")
            println("Mean post-merger markup change: $(mean(mkup_m)-mean(mkup_pre)), median post-merger markup change: $(median(mkup_m)- median(mkup_pre)), standard deviation: $(std(mkup_m))")
            println("Mean post-merger margin change: $(mean(marg_m)-mean(marg_pre)), median post-merger margin change: $(median(marg_m)-median(marg_pre)), standard deviation: $(std(marg_m))")
            println("Mean post-merger market share change: $(mean(s_m)-mean(s_pre)), median post-merger market share change: $(median(s_m)-median(s_pre)), standard deviation: $(std(s_m))")
    end
end