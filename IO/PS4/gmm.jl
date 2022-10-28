function gmmobj(theta2::Vector{Float64},data::Data, res::Results, df_out)
    @unpack invA, x1, IV = data
    res.theta2 = theta2
    delta = meanval(data, res)
    df = fill(100.0, 13)
    if maximum(isnan.(delta)) ==1
        f = 1e10
    else
        temp1 = x1'*IV
        temp2 = delta'*IV
        theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2'
        gmmresid = delta - x1*theta1
        res.gmmresid = vec(gmmresid)
        temp1 = gmmresid'*IV
        f1 = temp1*invA*temp1'
        f = f1
        temp = jacob(data,res)'
        df = 2*temp*IV*invA*IV'*gmmresid
    end
    if df_out==0
        return f
    else
        return f, vec(df)
    end
end