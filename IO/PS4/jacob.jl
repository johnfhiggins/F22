function jacob(data::Data, res::Results)
    @unpack ns, cdid, cdindex, x1, x2, vfull, dfull = data
    @unpack theti, thetj, theta2, mvalold = res
    theta2w = Array(sparse(theti, thetj, theta2))
    ###CHANGE THIS BACK
    expmu = fill(1.0, 94*24,20) #exp.(mufunc(x2, theta2w, data)) #
    shares = ind_sh(mvalold, expmu, data)[:,:,1]

    n, K = size(x2)
    J = size(theta2w, 2)-1
    f1 = zeros(size(cdid, 1), K*(J+1))

    for i=1:K
        xv = (x2[:,i]*ones(1,ns)) .* vfull[:, ns*(i-1)+1:ns*i]
        temp = cumsum(xv .*shares, dims=1)
        sum1 = temp[cdindex,:]
        sum1[2:size(sum1,1),:] = diff(sum1, dims=1)
        f1[:,i] .= mean((shares .*(xv .- sum1[cdid])))'
    end
    for j=1:J
        d = dfull[:, ns*(j-1)+1:ns*j]
        temp1 = zeros(size(cdid, 1), K)
        for i=1:K
            xd = (x2[:,i]*ones(1,ns)).*d
            temp = cumsum(xd .*shares, dims=1)
            sum1 = temp[cdindex,:]
            sum1[2:size(sum1,1), :] = diff(sum1, dims=1)
            temp1[:,i] .= mean((shares .*(xd .- sum1[cdid,:])))'
        end
        f1[:, K*j + 1:K*(j+1)] = temp1
    end

    rel = Int64.(theti + (thetj .- 1) *maximum(theti))
    f = zeros(size(cdid,1), size(rel,1))
    n=1
    for i=1:size(cdindex, 1)
        temp = shares[n:cdindex[i],:]
        H1 = temp*temp'
        H = (diagm(vec(sum(temp, dims=2))) - H1) ./ns
        f[n:cdindex[i],:] = -inv(H)*f1[n:cdindex[i], rel]
        n = cdindex[i] + 1
    end
    f
end