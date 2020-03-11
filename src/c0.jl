# Initialize 2D matrix w/rows corresponding to monomers
# and columns corresponding to ligands
# Indices are as follows:
# c[1,1] = M+, c[3,2] = M, c[2,3] = L, c[3,3] = ML
# Julia has 1-based indices therefore # of monomers = row index + 1 and
# # of ligands = column index + 1
function c0(imax, Ns, Mplus, L)

rows    = imax + 3
columns = Ns[imax] + 3

c0 = zeros(Float64, rows, columns)
c0[1,1] = Mplus
c0[2,3] = L

return c0

end
