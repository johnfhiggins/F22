
function bert_eq(p::Float64, expall_i::Array{Float64}, α_i::Array{Float64}, mc_hat1::Float64, own_dum::Array{Float64})
    eg = expall_i .*exp(fill(p, ns) .* kronecker(α_i, ones(24)))
    shar_i = eg ./(ones(length(eg[:,1]))*(1 + sum(eg)))
    f = omega1(p, shar_i, α_i', own_dum)*(p - mc_hat1) + (mean(shar_i')')
    f
end