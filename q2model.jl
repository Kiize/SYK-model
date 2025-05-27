using ITensors, ITensorMPS
using Random
using DelimitedFiles

ITensors.op(::OpName"M_even", ::SiteType"Fermion") = im*[0 -1;1 0]

ITensors.op(::OpName"M_odd", ::SiteType"Fermion") = [0 1; 1 0]

function syk_model_q2(N, J) #N_fermions, q = 2
    sites = siteinds("Fermion", N)

    os_j = Op("M_even", 1) + Op("M_odd", 1)
    os_i = Op("M_even", 1) + Op("M_odd", 1)
    
    for i=2:N
        for j=2:N
            os_j += Op("M_even", j) + Op("M_odd", j)
        end
        os_i += Op("M_even", i) + Op("M_odd", i)
    end
    os = os_i * os_j

    ope_local = OpSum()
    for i=1:N_majo, j=1:N_majo
        idx_i = [if (os[1][i].which_op == "M_even") 2*os[1][i].sites[1] else 2*os[1][i].sites[1] - 1 end]
        idx_j = [if (os[2][j].which_op == "M_even") 2*os[2][j].sites[1] else 2*os[2][j].sites[1] - 1 end]
        ope_local +=  J[idx_i[1], idx_j[1]], os[1][i].which_op, os[1][i].sites[1], os[2][j].which_op, os[2][j].sites[1]
    end

    H = MPO(ope_local, sites)

    state0 = [isodd(n) ? "1" : "0" for n=1:N]
    psi0_i = MPS(sites, state0)

    return H, psi0_i
end

N_majo = 8
N_fermions = N_majo ÷ 2 

# Gaussianly distributed coefficients (randn) with mean 0 and variance J^2/N^(q-1), J^2 = 1.
# The matrix must be hermitian to represent a physical system.

variance = 1/N_majo
J = randn(N_majo, N_majo) .* variance

J = (J + J')/2

# SYK Hamiltonian

H, psi0_i = syk_model_q2(N_fermions, J)

# Run the DMRG algorithm

sweeps = Sweeps(10)
setmaxdim!(sweeps,10,20,100,200,400,800)
setcutoff!(sweeps,1E-8)

energy0,psi0 = dmrg(H,psi0_i,sweeps)
energy1, psi1 = dmrg(H, [psi0], psi0_i, sweeps)

states_arr = [psi0, psi1]
energies_arr = [energy0, energy1]
N_levels = 10
for ind=1:(N_levels-2)
    local energy, psi = dmrg(H, states_arr, psi0_i, sweeps)
    global states_arr = [states_arr..., psi]
    append!(energies_arr, energy)
end

open("energies_q2.txt", "a") do io
    writedlm(io, (N_majo, energies_arr), "\n")
end