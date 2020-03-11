# kminus,i,j = k*(j/Ns,i)
function kminusij(imax, Ns)

rows    = imax + 3
columns = Ns[imax] + 3
istart  = 3
idimer  = 4
iend    = rows - 1
j0      = 2
jmax    = columns - 1

kmij = zeros(Float64, rows, columns)

for i = idimer:iend
    monomers = i - 2
    for j = j0:Ns[monomers] + 2
        ligands = j - 2
        kmij[i,j] = (ligands)/Ns[monomers]
    end
end

return kmij

end
