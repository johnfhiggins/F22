function meanval(data::Data, res::Results)
    @unpack theti, thetj, theta2, oldt2, mvalold= res
    @unpack s_jt, x2 = data
    if maximum(abs.(theta2 - oldt2 )) < 0.01
        tol = 1e-9
        flag = 0
    else
        tol = 1e-6
        flag = 1
    end

    theta2w = Array(sparse(theti, thetj, theta2))
    expmu = exp.(mufunc(x2,theta2w, data))
    norm = 1
    avgnorm = 1

    i = 0
    mval=zeros(3)
    while norm>tol*10^(flag*floor(i/50)) && avgnorm > 1e-3*tol*10^(flag*floor(i/50))
        #println(i)
        ###NOTE: CHANGE THIS BACK######
        mval = mvalold .*s_jt ./mktsh(mvalold, fill(1.0, 94*24,20), data) #log.(s_jt) .- log.(mktsh(mvalold, expmu, data)) #
        #println(mktsh(mvalold, expmu, data))
        t = abs.(mval-mvalold)
        norm = maximum(t)
        avgnorm = mean(t)
   
            println(i, ", norm: ", norm, ", avgnorm: ", avgnorm)
        
        mvalold = mval
        i += 1
    end


    if flag ==1 && maximum(isnan.(mval)) < 1
        print("Updated")
        mvalold = mval
        oldt2 = theta2
        res.mvalold = mvalold
        res.oldt2 = oldt2
    end
    
    log.(mval)
end