function omega1(p,shar_i,alpha_i,own_dum)
#
#inputs:		p--price(over which we will perform a search);
#				shar_i--individual purchase probabilities
#				alpha_i--individual price coeff
#				mc_hat1--marginal costs;
#				own_dum--ownership structure

    n=24
    o1=(ones(n)*alpha_i') .* shar_i*shar_i'
    o2=diagm(sum(((ones(n)*alpha_i').*shar_i)'))
    omega=(o2-o1)/ns
    f = omega.*(own_dum)
    f
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%9.1.1  Calculates the ownership share matrix x the (partial share)/(partial price)
#%





