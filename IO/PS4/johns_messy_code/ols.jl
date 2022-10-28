function ols(y, X, fid)
    b = inv(X'*X)*X'*y
    if fid==1
        e = y - X*b
        N = size(X,1)
        K = size(X, 2)
        se2 = (e'*e)./(N - K)
        varb = se2.*inv(X'*X)
        seb = varb[1].^(1/2)
        R2 = 1.0 .- sum(e'*e)/sum((y .- mean(y)).^2)
        adjR2 = 1 .- ((N-1)/(N - K))*(1 .-R2)
        return b, seb, R2, adjR2
    else
        return b
    end
end