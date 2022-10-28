function mktsh(mval, expmu, data::Data)
    @unpack ns = data
    val = sum(ind_sh(mval, expmu, data),dims=2)./ns
    val
end


#sh_test = mktsh(res.mvalold, fill(1.0, 94*24,20), Data())