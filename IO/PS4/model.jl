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

    IV ::Array{Float64,2} = hcat(cereal[:,12:31], brand_dummies)
    invA::Array{Float64,2} = inv(Array(IV)'*Array(IV))

    cdid::Kronecker.KroneckerProduct{Int64} = kronecker(1:nmkt,ones(Int64,nbrn,1))
    cdindex::StepRange{Int64,Int64} = nbrn:nbrn:nbrn*nmkt   

    ns::Int64 = 20
    vfull ::Array{Float64,2} = repeat(v, inner=(nbrn,1), outer=(1,1))
    dfull ::Array{Float64,2} = repeat(demog, inner=(nbrn,1), outer=(1,1))

    theta2w::Array{Float64,2} =    [0.3772    3.0888         0    1.1859         0;
    1.8480   16.5980    -.6590         0   11.6245;
   -0.0035   -0.1925         0    0.0296         0;
    0.0810    1.4684         0   -1.5143         0];

    thet::Vector{CartesianIndex{2}} = findall(!=(0), theta2w)
    theti::Vector{Int64} = getindex.(thet, [1])
    thetj::Vector{Int64} = getindex.(thet, [2])
    theta2::Vector{Float64} = theta2w[thet]
end
Data()

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
mv = meanval(Data(), res)

f= gmmobj(res.theta2,Data(), res,1)[2]

function gmm_wrapper(theta2)
    val = -gmmobj(theta2, Data(), res, 1)[1][1]
    val
end
function g!(storage, theta2)
    storage .= gmmobj(theta2, Data(), res, 1)[2]
    #val = vec(val)
    #val
end

@unpack theti, thetj,x1, x2, IV, invA= Data()
function gmm_estimation()
    opt = optimize(gmm_wrapper, g!, zeros(13), BFGS(),  Optim.Options(g_tol = 1e-20, iterations = 100, store_trace = true, show_trace = true))
    theta2 = opt.minimizer
    #theta2 = res.theta2
    varcov = var_cov(Data(), res)
    se = sqrt.(Complex.((diag(varcov))))
    theta2w = Array(sparse(theti, thetj,theta2 ))
    t = size(se,1) - size(theta2,1)
    se2w = Array(sparse(theti, thetj, se[t+1:size(se,1)]))

    omega = inv(varcov[2:25, 2:25])
    xmd = [x2[1:24,1] x2[1:24, 3:4]]

    delta = meanval(Data(),res)
     #println("***", delta)
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
        println("R^2: $(Rsq), adj R^2: $(Rsq_G), Chisq: $(Chisq)")
    end
end
gmm_estimation()
###setting theta 



#res.theta2 = [0.3772
#1.8480
#-0.0035
#0.0810
#3.0888
##16.5980
#-0.1926
#1.4684
#-0.6590
#1.1859
#0.0295
#-1.5143
#11.6245]

function alpha_finder(data::Data, res::Results, α_mean)
    @unpack ns, nmkt, nbrn, x1, x2, theti, thetj, vfull, dfull = data
    @unpack theta2 = res
    theta2w = Array(sparse(theti, thetj, theta2))
    n,k = size(x2)
    j = size(theta2w,2)-1
    alpha_mat = zeros(nmkt*nbrn, ns)
    for i=1:ns
        v_i = vfull[:, i:ns:ns*k]
        d_i = dfull[:, i:ns:ns*j]
        alpha_mat[:,i] = α_mean .+ x2[:,2] .* v_i[:,1] * theta2w[2,1] + (x2[:,2] .* d_i[:,1] * theta2w[2,2:j+1]')*ones(k)
    end
    alpha_mat
end


theta2w = Array(sparse(theti, thetj,res.theta2))
delta = meanval(Data(), res)
a_mat = alpha_finder(Data(), res, -32.497106)
expmu = exp.(mufunc(x2, theta2w, Data()))
pred_shares = ind_sh(exp.(delta), expmu, Data())

reg_results()
logit_merger_results()
blp_merger_results()

