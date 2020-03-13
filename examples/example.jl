using NanoclusterModeler
using DifferentialEquations
using ODEInterfaceDiffEq
using Queryverse
using DataFrames
using JLD2, FileIO

# Set parameters for model
imax = 68 # Maximum size cluster considered in model
prefactor = 2.08 # Scaling factor in Ns
Ns = NanoclusterModeler.Ns(imax, prefactor) # Number of surface sites on the cluster
Mplus = 0.05 # Concentration in mM
L     = 6.00 # Concentration in mM
c0 = NanoclusterModeler.c0(imax, Ns, Mplus, L) # Initialize 2D matrix containing intial concentrations of M+ and L
kpij = NanoclusterModeler.kplusij(imax, Ns) # Define 2D matrix containing fraction of free sites on surface of NCs 
kmij = NanoclusterModeler.kminusij(imax, Ns) # Define 2D matrix containing fraction of occupied sites on surface of NCs
kp  = 1.0e+3
kb  = 1.0e+5
kub = 1.0e-7
kn  = 1.0e+1
kg1 = 1.0e+4
kd1 = 1.0e-9
kg2 = 1.0e+4
kd2 = 1.0e-9
ka  = 1.0e-3
ke  = 1.0e+3
kc  = 1.0e+3
odes = NanoclusterModeler.odes!

# Set file prefix for output
file = "kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e-3_kc1e3_M+005_Ns208"

# Time is seconds
tspan = (0.0, 1000000.0)

# Parameters being passed to ode function
p   = (imax, Ns, kpij, kmij,
    kp, kb, kub, kn, kg1, kd1, kg2, kd2, ka, ke, kc)

# Define ode problem and solve
prob = ODEProblem(odes, c0, tspan, p)
sol  = solve(prob, radau(), reltol=1e-10, abstol=1e-10)

# Post processing of data
rows, columns, tlast = size(sol)
t = sol.t

# If solution is negative, it is out of the tolerance. Set to zero.
for t = 1:tlast
    for j = 1:columns
        for i = 1:rows
            if sol[i,j,t] < 0.0
                sol[i,j,t] = 0.0
            end
        end
    end
end

Mplus = sol[1,1,:]
M     = sol[3,2,:]
L     = sol[2,3,:]
ML    = sol[3,3,:]
C_20  = sol[4,2,:]
C_21  = sol[4,3,:]
C_22  = sol[4,4,:]
C_23  = sol[4,5,:]
C_2   = hcat(C_20, C_21, C_22, C_23)

C_tot = zeros(Float64, tlast)
C_bar = zeros(Float64, imax, tlast)
for t = 1:tlast
    C_tot[t]   = sum(sol[3:rows-1, 2:columns-1, t])
    C_bar[:,t] = sum(sol[3:rows-1, 2:columns-1, t], dims=2)
end

file_time = string(file, "_time.csv")
df_t = DataFrame([sol.t, Mplus, M, L, ML, C_20, C_21, C_22, C_23, C_tot])
rename!(df_t, [:t, :Mplus, :M, :L, :ML, :C_20, :C_21, :C_22, :C_23, :C_tot])
save(file_time, df_t)

file_Cbar = string(file, "_cbar.csv")
C_bar_T = transpose(C_bar)
df_Cbar = DataFrame(C_bar_T)
save(file_Cbar, df_Cbar)

P_j = zeros(Float64, imax, Ns[imax]+1, tlast)
for t = 1:tlast
    for j = 2:columns
        for i = 4:rows
            if sol[i,j,t] >= 1e-10
                P_j[i-2,j-1,t] = sol[i,j,t]/C_bar[i-2,t]
            end
        end
    end
end

file_P2D_t = string(file, "_prob_j_t.jld2")
@save file_P2D_t P_j

file_sol = string(file, "_full_sol.jld2")
@save file_sol sol
