function mufunc(x2::Array{Float64,2}, theta2w::Array{Float64,2}, data::Data)
    @unpack ns, vfull, dfull = data
    n,k = size(x2)
    j = size(theta2w,2)-1
    mu = zeros(n, ns)
    for i=1:ns
        v_i = vfull[:, i:ns:ns*k]
        d_i = dfull[:, i:ns:ns*j]
        mu[:,i] = (x2 .* v_i * theta2w[:,1]) + x2 .*(d_i *theta2w[:, 2:j+1]') * ones(k)
    end
    mu
end
