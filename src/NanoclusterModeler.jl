module NanoclusterModeler

greet() = print("Hello World!")

export odes!

# Define max # of monomers imax
imax = 40

# Define # of binding sites for clusters with n monomers up to imax
n  = range(1, stop=imax)
Ns = round.(Int64, 2.08.*n.^(2/3))

# Define important indices for 2D array
rows    = imax + 3
columns = Ns[imax] + 3
istart  = 3
idimer  = 4
iend    = rows - 1
j0      = 2
jmax    = columns - 1

# Initialize 2D matrix w/rows corresponding to monomers
# and columns corresponding to ligands
# Indices are as follows:
# c[1,1] = M+, c[3,2] = M, c[2,3] = L, c[3,3] = ML
# Julia has 1-based indices therefore # of monomers = row index + 1 and
# # of ligands = column index + 1
c0 = zeros(Float64, rows, columns)
c0[1,1] = 0.05
c0[2,3] = 6.00

# Define parameters
# kplus,i,j = k*(1.0 - j/Ns,i)
kpij = zeros(Float64, rows, columns)
for i = idimer:iend
    monomers = i - 2
    for j = j0:Ns[monomers] + 2
        ligands = j - 2
        kpij[i,j] = 1.0 - (ligands)/Ns[monomers]
    end
end

# kminus,i,j = k*(j/Ns,i)
kmij = zeros(Float64, rows, columns)
for i = idimer:iend
    monomers = i - 2
    for j = j0:Ns[monomers] + 2
        ligands = j - 2
        kmij[i,j] = (ligands)/Ns[monomers]
    end
end

# p = (kp,     kb,     kub,    kn,     kg1,    kd1,    kg2,    kd2,    ka,     ke,     kc)
p   = (1.0e+3, 1.0e+1, 1.0e-7, 1.0e+1, 1.0e+4, 1.0e-9, 1.0e+4, 1.0e-9, 1.0e1, 1.0e-3, 1.0e+3,
        imax, idimer, iend, j0, jmax, Ns, kpij, kmij)

function odes!(dc, c, p, t)
    # unpack parameters
    kp, kb, kub, kn, kg1, kd1, kg2, kd2, ka, ke, kc, imax, idimer, iend, j0, jmax, Ns, kpij, kmij = p

    # array element (1,1) is assigned to M+ or c1,0+
    dc[1,1] = -kp*c[1,1]

    # array element (3,2) is assigned to M or c1,0
    c_kg2 = 0.0
    for j = j0:jmax
        for i = idimer:iend-1
            c_kg2 += kpij[i,j]*c[i,j]
        end
    end
    dc[3,2] = kp*c[1,1] - kb*c[3,2]*c[2,3] + kub*c[3,3] -
        2*kn*c[3,2]*c[3,2] - kn*c[3,2]*c[3,3] - kg2*c[3,2]*c_kg2

    # array element (3,3) is assigned to ML or c1,1
    c_kg1 = 0.0
    for j = j0:jmax
        for i = idimer:iend-1
            c_kg1 += kpij[i,j]*c[i,j]
        end
    end
    c_kd1 = 0.0
    for j = j0:jmax
        for i = idimer+1:iend
            c_kd1 += kmij[i,j]*c[i,j]
        end
    end
    dc[3,3] = kb*c[3,2]*c[2,3] - kub*c[3,3] - 2*kn*c[3,3]*c[3,3] - kn*c[3,3]*c[3,2] -
        kg1*c[3,3]*c_kg1 + kd1*c_kd1

    # array element (2,3) is assigned to L or c0,1
    c_ka = 0.0
    for j = j0:jmax
        for i = idimer:iend
            c_ka += kpij[i,j]*c[i,j]
        end
    end
    c_ke = 0.0
    for j = j0:jmax
        for i = idimer:iend
            c_ke += kmij[i,j]*c[i,j]
        end
    end
    dc[2,3] = -kb*c[3,2]*c[2,3] + kub*c[3,3] - ka*c[2,3]*c_ka + ke*c_ke

    # nested loop provides expression for ci,j starting at c2,0
    # ligands are the outer loop for memory efficiency in julia
    for j = j0:Ns[imax]+2
        for i = idimer:iend
            # # of monomers is index i minus 2
            # # of ligands on NC is the index j minus 2
            monomers = i - 2
            ligands  = j - 2
            # don't access array elements that correspond to non-existent clusters
            if ( ligands <= Ns[monomers] )

                # sum up clusters consumed by coalescence c_cns = kc * kpij *ci,j * sum_k,l kpkl * ck,l
                # k, l are dummy indices
                # conditions: i+k <= imax, j+l <= Ns[i+k]
                c_cns = 0.0
                if i < iend - 1
                for l = j0:jmax
                    for k = idimer:iend-monomers+2
                        if ( (i-2)+(k-2) <= imax)
                            if ( (j-2)+(l-2) <= Ns[(i-2)+(k-2)] )
                                c_cns += kpij[i,j]*c[i,j]*kpij[k,l]*c[k,l]
                            end
                        end
                    end
                end
                end

                # sum up clusters produced by coalescence c_prd = kc * sum_w,y sum_x,z kpwx * cw,x * kpyz * cy,z
                # w, x, y, z are dummy indices
                # conditions: i=w+y, j=x+z
                # sum over-counts by factor of 4, so divide by 4 later
                c_prd = 0.0
                if i > idimer + 1
                for x = j0:jmax
                    for z = j0:jmax
                        for w = idimer:i
                            for y = idimer:i
                                if ( (i-2) == (w-2)+(y-2) )
                                    if ( (j-2) == (x-2)+(z-2) )
                                        c_prd += kpij[w,x]*c[w,x]*kpij[y,z]*c[y,z]
                                    end
                                end
                            end
                        end
                    end
                end
                end

                # if clause provides special cases for ci,j
                if i == idimer
                    if j == j0
                        dc[i,j] = kn*c[3,2]*c[3,2] -
                            kg1*kpij[i,j]*c[3,3]*c[i,j] + kd1*kmij[i+1,j+1]*c[i+1,j+1] -
                            kg2*kpij[i,j]*c[3,2]*c[i,j] -
                            ka*kpij[i,j]*c[2,3]*c[i,j] +
                            ke*kmij[i,j+1]*c[i,j+1] -
                            kc*c_cns
                    elseif j == j0 + 1
                        dc[i,j] = kn*c[3,3]*c[3,2] -
                            kg1*kpij[i,j]*c[3,3]*c[i,j] + kd1*kmij[i+1,j+1]*c[i+1,j+1] -
                            kg2*kpij[i,j]*c[3,2]*c[i,j] +
                            ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                            ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] -
                            kc*c_cns
                    elseif j == j0 + 2
                        dc[i,j] = kn*c[3,3]*c[3,3] -
                            kg1*kpij[i,j]*c[3,3]*c[i,j] + kd1*kmij[i+1,j+1]*c[i+1,j+1] -
                            kg2*kpij[i,j]*c[3,2]*c[i,j] +
                            ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                            ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] -
                            kc*c_cns
                    elseif j == j0 + 3
                        dc[i,j] = -kg1*kpij[i,j]*c[3,3]*c[i,j] + kd1*kmij[i+1,j+1]*c[i+1,j+1] -
                            kg2*kpij[i,j]*c[3,2]*c[i,j] +
                            ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                            ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] -
                            kc*c_cns
                    end
                elseif i == idimer + 1
                    dc[i,j] = kg1*kpij[i-1,j-1]*c[3,3]*c[i-1,j-1] - kg1*kpij[i,j]*c[3,3]*c[i,j] +
                        kd1*kmij[i+1,j+1]*c[i+1,j+1] - kd1*kmij[i,j]*c[i,j] +
                        kg2*kpij[i-1,j]*c[3,2]*c[i-1,j] - kg2*kpij[i,j]*c[3,2]*c[i,j] +
                        ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                        ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] -
                        kc*c_cns
                elseif i > idimer + 1 && i < iend - 1
                    dc[i,j] = kg1*kpij[i-1,j-1]*c[3,3]*c[i-1,j-1] - kg1*kpij[i,j]*c[3,3]*c[i,j] +
                        kd1*kmij[i+1,j+1]*c[i+1,j+1] - kd1*kmij[i,j]*c[i,j] +
                        kg2*kpij[i-1,j]*c[3,2]*c[i-1,j] - kg2*kpij[i,j]*c[3,2]*c[i,j] +
                        ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                        ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] -
                        kc*c_cns +
                        kc*c_prd/4.0
                elseif i == iend - 1
                    dc[i,j] = kg1*kpij[i-1,j-1]*c[3,3]*c[i-1,j-1] - kg1*kpij[i,j]*c[3,3]*c[i,j] +
                        kd1*kmij[i+1,j+1]*c[i+1,j+1] - kd1*kmij[i,j]*c[i,j] +
                        kg2*kpij[i-1,j]*c[3,2]*c[i-1,j] - kg2*kpij[i,j]*c[3,2]*c[i,j] +
                        ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                        ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] +
                        kc*c_prd/4.0
                elseif i == iend
                    dc[i,j] = kg1*kpij[i-1,j-1]*c[3,3]*c[i-1,j-1] +
                        kd1*kmij[i+1,j+1]*c[i+1,j+1] - kd1*kmij[i,j]*c[i,j] +
                        kg2*kpij[i-1,j]*c[3,2]*c[i-1,j] +
                        ka*kpij[i,j-1]*c[2,3]*c[i,j-1] - ka*kpij[i,j]*c[2,3]*c[i,j] +
                        ke*kmij[i,j+1]*c[i,j+1] - ke*kmij[i,j]*c[i,j] +
                        kc*c_prd/4.0
                end
            end
        end
    end
end

end # module
