

function own_dum_comp(firm_own)
    own_mat = zeros(24,24)
    for i=1:24 #loop over each product
        for j = 1:24
            own_mat[i, j] = (firm_own[i] .== firm_own[j]) #determine which brands are owned by the same firm as product i, store as 1 if yes and 0 of not
        end
    end
    own_mat
end

function mkup_i( s_jt, alpha_i, own_mat_d) #data is market-level share data; only shares from products in the same market
    mkups = zeros(length(s_jt))
    omeg = omega1(s_jt, alpha_i, own_mat_d)
    mkups = inv(omeg)*s_jt
end

function omega1(shar_i,alpha_i,own_dum)
    #
    #inputs:		p--price(over which we will perform a search);
    #				shar_i--individual purchase probabilities
    #				alpha_i--individual price coeff
    #				mc_hat1--marginal costs;
    #				own_dum--ownership structure
    n=24
        o1=(ones(n).*alpha_i') * shar_i*shar_i'
        o2=diagm(vec(sum(((ones(n)*alpha_i').*shar_i)', dims=2)))
        omega=(o2-o1) #need to do a random coefficients version where we aggregate across consumers within each market? This is why the original code divides by 20
        f = omega.*(own_dum)
        f
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
    s_jt = cereal[:,8]
    prices = zeros(N)
    alph_vec = fill(abs(beta[1]), 94*24)
    for m = 1:nm
        mkt_range = (m-1)*nbr + 1:m*nbr
        mkup[mkt_range] .= mkup_i(prev_shares[mkt_range], alph_vec[mkt_range], own_dum_m)
        prices[mkt_range] .= mkup[mkt_range] + mc[mkt_range]
        marg[mkt_range] .= mkup[mkt_range]./(prices[mkt_range])
    end
    #find the shares at the new prices:

    if reg_type==1 #if we want regression only on prices
        pred_delta = prices .*beta #beta will be just alpha
    elseif reg_type ==2 #if we want regression on prices with fixed effects #beta will have alpha in the first column and the individual fixed effects in the second
        pred_delta = X_FE*beta
    elseif reg_type ==3 #if we want IV with only prices, no fixed effects (beta comes from IV regression)
        pred_delta = prices .*beta
    else #if we want IV with prices and fixed effects (beta comes from IV reg with fixed effects)
        pred_delta = X_FE*beta
    end
    #find the numerator of the expected share
    raw_shares = exp.(pred_delta)
    #call the share finder function to find the actual share. It'll be exp(δ_j)/(1 + sum(exp(δ_i))) applied to each element
    s_jt = share_finder(raw_shares, cdindex, cdid)
    #return all of the desired quantities
    mkup, marg, s_jt, prices
end

#estimate the four coefficients below:
@unpack y, brand_dummies= Data()

Z = cereal[:,12:21]
p = cereal[:, 9]
X_FE = hcat(p, brand_dummies)

b_ols = ols(y, p, 0)
b_ols_FE = ols(y, X_FE, 0)
b_IV = regressIV(y, p, Z, 0)
b_IV_FE = regressIV_FE(y, p, brand_dummies ,Z, 0)

α_ols = b_ols[1]
α_ols_FE = b_ols_FE[1]
α_IV = b_IV[1]
α_IV_FE = b_IV_FE[1]

dum = own_dum_comp(cereal[:,3])

alpha_list = [abs(α_ols), abs(α_ols_FE), abs(α_IV), abs(α_IV_FE)]
for alph in alpha_list
    mkup, mc, marg , p_new, s_mat = mkup_finder_pre(fill(alph, 94*24), dum)
    println("α = $(alph):")
    println("Mean markup = $(mean(mkup)), median markup = $(median(mkup)), standard deviation = $(std(mkup))")
    println("Mean mc = $(mean(mc)), median mc = $(median(mc)), standard deviation = $(std(mc))")
    println("Mean margin = $(mean(marg)), median margin = $(median(marg)), standard deviation = $(std(marg))")
end

mkup_ols_pre, mc_ols_pre, marg_ols_pre , p_ols_pre, s_ols_pre = mkup_finder_pre(fill(29.454, 94*24), dum)

###merger simulation section!!!###
cereal_pn = copy(cereal) #copy data 
cereal_pn[cereal_pn[:,3] .== 6,3 ] .= 3 #merge post and nabisco by setting each Nabisco firm entry to 3 (the post firm ID)

cereal_gmq = copy(cereal)
cereal_gmq[cereal_gmq[:,3] .== 4,3] .= 2 #merge general mills and quaker by setting each quaker entry to 2 (the General mills firm ID)

mkup_pre, mc_pre, marg_pre , p_0, s_pre = mkup_finder_pre(fill(18.6011, 94*24), dum)
dum_test = own_dum_comp(cereal_pn[1:24,3])
mkt, mgt, st, pt = mkup_finder_merge(dum_test, st, mc_pre, b_IV_FE, 4) #fill(abs(α_ols), 94*24), 1)
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
    α = abs.(β[1])
    #find initial mkups, marginal costs, margins, prices, and shares given regression type
    mkup_pre, mc_pre, marg_pre , p_0, s_pre = mkup_finder_pre(fill(α, 94*24), dum)
    println(mean(p_0))
    own_dum_m = own_dum_comp(dataset[1:24,3])
    mk_t, marg_t, s_t, p_1 = mkup_finder_merge(own_dum_m, s_pre, mc_pre, β, reg_type)
    tol = 1e-9
    error = maximum(abs.(p_1 - p_0))
    while error > tol
        p_0 = p_1
        mk_t, marg_t, s_t, p_1 = mkup_finder_merge(own_dum_m, s_t, mc_pre, β, reg_type)
        error = maximum(abs.(p_1 - p_0))
        println(mean(p_1), ", ", mean(marg_t), ", ", mean(s_t))
    end
    print("Convergence!")
    #mkup_pre, mc_pre, marg_pre, p_0, s_pre, mk_t, marg_t, p_1, s_t
end

merger_sim(cereal_pn, 4)

#=

mk_t, mc_t, marg_t, p_m_t, s_t = merger_sim_iter(cereal_pn, b_ols, mc_ols_pre, s_ols_pre,1)

fig = plot([s_ols_pre[1:300], s_t[1:24]])

for i=1:100
    mk_t, mc_t, marg_t, p_m_t, s_t = merger_sim_iter(cereal_pn, b_ols, mc_ols_pre, s_t,1)
    println(mean(p_m_t))
    plot!(fig, s_t[1:24])
end
display(fig)


mkt, mgt, st, pt = mkup_finder_merge(fill(abs(α_ols), 94*24), dum, Data().s_jt, fill(0.12,24*94), b_ols, 1)
histogram(pt)
=#