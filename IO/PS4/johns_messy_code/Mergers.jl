

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
        o1=(ones(n)*alpha_i') .* shar_i*shar_i'
        o2=diagm(sum(((ones(n)*alpha_i').*shar_i)'))
        omega=(o2-o1)/ns
        f = omega.*(own_dum)
        f
end

function mkup_finder_pre(alpha_i, own_dum, cereal)
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

function mkup_finder_merge(alpha_i, own_dum, cereal, mc, x1; x2)


