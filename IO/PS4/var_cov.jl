###need to import mvalold, ps2, gmmresid from structs
###To do: we need to create results struct to hold these things! 

#also: need to have a data struct holding x1 and IV (and potentially other stuff)!

#function accepts theta2 and results struct
function var_cov(data_s::Data, res::Results)
    @unpack mvalold, gmmresid = res
    @unpack x1, IV, invA = data_s
    N = size(x1,1) #guessing we want the number of observations in x1
    Z = size(IV,2) #I'm guessing we want the number of instruments
    temp = jacob(data_s, res) #call the jacobian function and see what's up
    a = [x1 temp]' * IV
    IVres = IV .* (gmmresid .* ones(1,Z))
    b = IVres'*IVres 

    f = inv(a*invA*a')*a*invA*b*invA*a'*inv(a*invA*a') #I am assuming we want this boi returned
    f #return :)
end