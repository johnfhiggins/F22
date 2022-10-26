using Kronecker
function ind_sh(expmval, expmu, data::Data)
    @unpack ns, cdindex, cdid = data 
    eg = expmu .* kronecker(ones(1, ns), expmval)
    temp = cumsum(eg, dims=1) #rather unfortunate function name
    #i'm assuming we want the sum along the first dimension (i.e. add up each column)
    sum1 = temp[cdindex, :] #find the last entry in each column (which will be the cumulative sum of all of the entries in the column)
    sum1[2:length(sum1[:,1]), :] = diff(sum1, dims=1) #not entirely sure my syntax is right here

    denom1 = 1 ./(1 .+ sum1)
    denom = denom1[cdid,:]
    f = eg .*denom
    f
end