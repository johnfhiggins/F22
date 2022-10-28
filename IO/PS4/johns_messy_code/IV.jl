function regressIV(y, X, Z, fid)
    #Xhat = Z*inv(Z'*Z)*Z'*XI
    #Xhat = [Xhat X]
    #X = [XI X]
    invA = inv(Z'*Z)
    mid = X'*Z*invA*Z'
    b = inv(mid*X)*mid*y
    #b = inv(Xhat'*Xhat)*Xhat'y
    e = y - X*b
    #println((e'* e)[1,1])
    N = size(X,1)
    K = size(X,2)
    se2 = (e'*e)[1,1]/(N - K)
    varb = se2 .*[inv(X'*X)]
    seb = varb[1,1].^(1/2)
    R2 = 1 .- (e'*e)/sum((y .- mean(y)).^2)
    adjR2 = 1 .- ((N-1)/(N - K))*(1 .-R2)
    if fid==1
        return b,  seb,R2, adjR2
    else
        return b
    end
end

function regressIV_FE(y, XI, X, Z, fid)
    Xhat = Z*inv(Z'*Z)*Z'*XI
    Xhat = [Xhat X]
    X = [XI X]
    b = inv(Xhat'*Xhat)*Xhat'y
    e = y - X*b
    println((e'* e)[1,1])
    N = size(X,1)
    K = size(X,2)
    se2 = (e'*e)[1,1]/(N - K)
    varb = se2 .*[inv(Xhat'*Xhat)]
    seb = (abs.(varb[1,1][1])).^(1/2)
    R2 = 1 .- (e'*e)/sum((y .- mean(y)).^2)
    adjR2 = 1 .- ((N-1)/(N - K))*(1 .-R2)
    if fid==1
        return b, seb,R2, adjR2
    else
        return b
    end
end

