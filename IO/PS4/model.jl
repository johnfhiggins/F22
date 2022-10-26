using CSV, DataFrames, StatsBase, LinearAlgebra, SparseArrays, Kronecker, Parameters, Optim, Plots

cereal = Array(DataFrame(CSV.File("cereal.csv")))
cdid = kronecker(1:94,ones(Int64,24,1))
cdindex = 24:24:24*94   

temp = cumsum(cereal[:,8])
sum1 = temp[cdindex,:]
sum1[2:size(sum1,1),:] = diff(sum1, dims=1)
outshr = 1.0 .- sum1
outshr_copied = kronecker(outshr, ones(24))

@with_kw struct Data
    y::Matrix{Float64} = log.(cereal[:,8]) - log.(outshr_copied)
    demog::Array{Float64,2} = Array(DataFrame(CSV.File("demog.csv", header=false)))
    v::Array{Float64,2} = Array(DataFrame(CSV.File("v.csv")))[:,4:83]
    nmkt::Int64 = 94
    nbrn ::Int64 = 24
    s_jt ::Vector{Float64} = cereal[:, 8]
    brand_dummies ::Array{Float64,2} = sparse(kronecker(ones(nmkt), diagm(ones(nbrn))))
    x1 ::Array{Float64, 2} = hcat(cereal[:,9], brand_dummies)
    x2 ::Array{Float64,2} = hcat(ones(nmkt*nbrn), cereal[:, 9:11]) #x2 consists of an intercept term, prices, sugar content, and mushy dummy

    IV ::Array{Float64,2} = hcat(cereal[:,12:21], brand_dummies)
    invA::Array{Float64,2} = inv(Array(IV)'*Array(IV))

    cdid::Kronecker.KroneckerProduct{Int64} = kronecker(1:nmkt,ones(Int64,nbrn,1))
    cdindex::StepRange{Int64,Int64} = nbrn:nbrn:nbrn*nmkt   

    ns::Int64 = 20
    vfull ::Array{Float64,2} = repeat(v, inner=(nbrn,1), outer=(1,1))
    dfull ::Array{Float64,2} = repeat(demog, inner=(nbrn,1), outer=(1,1))

    theta2w::Array{Float64,2} =    [0.3    5.0      0.0    0.2         0.0;
             2.2   13.0    -1.0         0.0   2.5;
            0.1  -0.2         0.0    0.03         0;
             0.2    1.3         0   -0.8         0]

    thet::Vector{CartesianIndex{2}} = findall(!=(0), theta2w)
    theti::Vector{Int64} = getindex.(thet, [1])
    thetj::Vector{Int64} = getindex.(thet, [2])
    theta2::Vector{Float64} = theta2w[thet]
end
Data()
Data().vfull

mutable struct Results
    mvalold ::Array{Float64}
    oldt2 :: Vector{Float64}
    theti ::Vector{Float64}
    thetj ::Vector{Float64}
    theta2 ::Vector{Float64}
    gmmresid::Vector{Float64}
end

function Initialize()
    data = Data()
    mvalold = zeros(data.nmkt*data.nbrn)
    oldt2 = zeros(13)
    theti = data.theti
    thetj = data.thetj
    theta2 = data.theta2
    gmmresid = zeros(13)
    res = Results(mvalold, oldt2, theti, thetj, theta2, gmmresid)
    res
end

res = Initialize()
#I think this is all of them lol
include("mufunc.jl")
include("meanval.jl")
include("gmm.jl")
include("jacob.jl")
include("ind_sh.jl")
include("mktsh.jl")
include("var_cov.jl")
include("ols.jl")
include("IV.jl")
include("Mergers.jl")


function iv_initial(data::Data, res::Results)
    @unpack x1, IV, invA,s_jt,y = data
    mid = x1'*IV*invA*IV'
    t = inv(mid*x1)*mid*y
    mvalold = x1*t
    res.mvalold = exp.(mvalold)
end
#need to run this to get the initial predicted market shares 
iv_initial(Data(), res)
#sh_test = mktsh(res.mvalold, exp.(mufunc(Data().x2, Data().theta2w, Data())), Data())
#mv = meanval(Data(), res)

#f, df = gmmobj(Data(), res)

@unpack theti, thetj,x1, x2, IV, invA= Data()

opt = optimize(theta2 ->  -gmmobj(theta2, Data(), res)[1][1], res.theta2)
theta2 = opt.minimizer
vcov = var_cov(Data(), res)
###I'm cheating rn; remove the absolute value!
se = sqrt.(diag(abs.(vcov)))
theta2w = Array(sparse(theti, thetj,theta2 ))
t = size(se,1) - size(theta2,1)
se2w = Array(sparse(theti, thetj, se[t+1:size(se,1)]))

omega = inv(vcov[2:25, 2:25])
xmd = [x2[1:24,1] x2[1:24, 3:4]]

delta = meanval(Data(),res)
temp1 = x1'*IV;
temp2 = delta'*IV;
theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2';
ymd = theta1[2:25]
beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd
resmd = ymd - xmd*beta
semd = sqrt.(diag(inv(xmd'*omega*xmd)))
mcoef = [beta[1]; theta1[1]; beta[2:3]]
semcoef = vcat(semd[1], se[1], semd)
Rsq = 1-((resmd .-mean(resmd))'*(resmd .-mean(resmd)))./((ymd .-mean(ymd))'*(ymd .-mean(ymd)))
Rsq_G = 1-(resmd'*omega*resmd)./((ymd .-mean(ymd))'*omega*(ymd.-mean(ymd)))
Chisq = size(x1,1)*resmd'*omega*resmd

for i=1:size(theta2w,1)
    println(mcoef[i], theta2w[i,:])
    println(semcoef[i], se2w[i,:])
end


##### ignore the code below this, this was my initial pass at the merger stuff and is not relevant for the above parts ###
###I am including my new merger functions in a separate Mergers.jl file###

#this function takes the data and returns the ownership matrix
function own_dum_comp(firm_own)
    #total number of products in each market
    own_mat = zeros(24,24)
    for i=1:24 #loop over each product
        for j = 1:24
            own_mat[i, j] = (firm_own[i] .== firm_own[j]) #determine which brands are owned by the same firm as product i, store as 1 if yes and 0 of not
        end
    end
    #own_mat[24,24] = 1 #outside option 
    own_mat
end

function mkup_i(shares_t, alpha_i, own_mat_d) #data is market-level share data; only shares from products in the same market
    mkups = zeros(length(shares_t))
    #shares_0 = vcat(shares_t, 1-sum(shares_t))
    omeg = omega1(29.45, shares_t, alpha_i, own_mat_d)
    mkups = inv(omeg)*shares_t
end

own_dum = own_dum_comp(cereal_a[1:24,3])

function mkup_finder(alpha, own_dum)
    mkup_mat = zeros(94, 24)
    mc_mat = zeros(94, 24)
    marg_mat = zeros(94,24)
    p_mat = zeros(94, 24)
    s_mat = zeros(94, 24)
    for nm = 1:94
        shares_i = cereal_a[:,8][(nm-1)*24 + 1: nm*24]
        mkup_mat[nm, :] = mkup_i(shares_i, alpha, own_dum)
        mc_mat[nm, :] = cereal_a[:,9][(nm-1)*24 + 1: nm*24] - mkup_mat[nm,:]
        marg_mat[nm, :] = mkup_mat[nm, :]./(cereal_a[:,9][(nm-1)*24 + 1: nm*24] )
        p_mat[nm, :] = cereal_a[:,9][(nm-1)*24 + 1: nm*24]
        s_mat[nm,:] = shares_i
    end
    mkup_mat, mc_mat, marg_mat, p_mat, s_mat
end

function mkup_finder_merger(alpha, own_dum_m, mc_mat, shares_t)
    mkup_mat = zeros(94, 24)
    marg_mat = zeros(94,24)
    p_mat = zeros(94, 24)
    share_mat = zeros(94, 24)
    for nm = 1:94
        #shares_i = cereal_a[:,8][(nm-1)*24 + 1: nm*24] ##Need to accept passed shares vector, then recompute shares at end of function based on prices
        mkup_mat[nm, :] = mkup_i(shares_t[nm, :], alpha, own_dum_m)
        #mc_mat[nm, :] = cereal_a[:,9][(nm-1)*24 + 1: nm*24] - mkup_mat[nm,:]
        p_mat[nm, :] = mc_mat[nm,:] + mkup_mat[nm,:]
        marg_mat[nm, :] = mkup_mat[nm, :]./p_mat[nm,:]
        share_mat[nm,:] = mkt_shares(p_mat[nm,:], alpha)
    end
    mkup_mat, mc_mat, marg_mat, p_mat, share_mat
end

alpha_list = [29.454, 29.037, 29.5266, 30.1856]
for alph in alpha_list
    mkup, mc, marg , p_new, s_mat = mkup_finder(alph, own_dum)
    println("α = $(alph):")
    println("Mean markup = $(mean(mkup)), median markup = $(median(mkup)), standard deviation = $(std(mkup))")
    println("Mean mc = $(mean(mc)), median mc = $(median(mc)), standard deviation = $(std(mc))")
    println("Mean margin = $(mean(marg)), median margin = $(median(marg)), standard deviation = $(std(marg))")
end

mkup_pre, mc_pre, marg_pre , p_pre, s_pre= mkup_finder(30, own_dum)
mean(p_pre)

###Merger simulation
cereal_pn = copy(cereal_a) #copy data 
cereal_pn[cereal_pn[:,3] .== 6,3 ] .= 3 #merge post and nabisco by setting each Nabisco firm entry to 3 (the post firm ID)

cereal_gmq = copy(cereal_a)
cereal_gmq[cereal_gmq[:,3] .== 4,3] .= 2 #merge general mills and quaker by setting each quaker entry to 2 (the General mills firm ID)

function mkt_shares(p, alpha) #start simple: only price coef
    shares_vec = zeros(length(p))
    shares_vec .= exp.(-alpha.*p)./(1 + sum(exp.(-alpha.*p)))
    shares_vec
end


function merger_sim(dataset, alpha, mc_pre, prev_share)
    own_dum_m = own_dum_comp(dataset[1:24,3])
    mkup, mc, marg, p_m, s_mat = mkup_finder_merger(alpha, own_dum_m, mc_pre, prev_share )
    mkup, mc_pre, marg, p_m, s_mat
end
    
mk_t, mc_t, marg_t, p_m_t, s_t = merger_sim(cereal_pn, 30, mc_pre, s_pre)
fig = plot([cumsum(s_pre, dims=1)[94,:]/94, cumsum(s_t, dims=1)[94,:]/94])

for i=1:100
    mk_t, mc_t, marg_t, p_m_t, s_t = merger_sim(cereal_pn, 30, mc_pre, s_t)
    println(mean(p_m_t))
    plot!(fig, cumsum(s_t, dims=1)[94,:]/94)
end
display(fig)

fig = plot([cumsum(s_pre, dims=1)[94,:]/94, cumsum(s_t, dims=1)[94,:]/94])


b_OLS, se_OLS = ols(y, cereal_a[:,9],1)
b_IV, se_IV = regressIV(y, cereal_a[:,9], zeros(24*94), cereal_a[:,12:21],1)

brand_dummies = sparse(kronecker(ones(nmkt), diagm(ones(nbrn))))

x1 = hcat(cereal_a[:,9], brand_dummies)
IV = hcat(cereal_a[:,12:21], brand_dummies)
invA = inv(Array(IV)'*Array(IV))




