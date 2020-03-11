# Define # of binding sites for clusters with n monomers up to imax
function Ns(imax, prefactor)

n  = range(1, stop=imax)
Ns = round.(Int64, prefactor.*n.^(2/3))

return Ns

end
