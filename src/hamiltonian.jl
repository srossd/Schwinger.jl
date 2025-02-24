abstract type SchwingerHamiltonian{N,F} end

struct EDHamiltonian{N,F} <: SchwingerHamiltonian{N,F}
    lattice::SchwingerLattice{N,F}
    matrix::SparseMatrixCSC{ComplexF64,Int64}
    L_max::Int64

    function EDHamiltonian(lattice::SchwingerLattice{N,F}, matrix::SparseMatrixCSC{ComplexF64,Int64}, L_max::Int) where {N,F}
        new{N,F}(lattice, matrix, L_max)
    end
end

"""
`EDHamiltonian(lattice)`
Computes the Hamiltonian for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDHamiltonian(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += q
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    states = basis(lattice; L_max = L_max, q = lattice.q, universe = universe)
    positions = Dict{Tuple{BitMatrix,Int},Int}((states[i].occupations, states[i].L₀) => i for i in eachindex(states))

    I = Vector{Int}(undef, length(states)*N*F*2)
    J = Vector{Int}(undef, length(states)*N*F*2)
    V = Vector{ComplexF64}(undef, length(states)*N*F*2)
    idx = 1
    for i in eachindex(states)
        I[idx] = i
        J[idx] = i
        V[idx] = energy(states[i])
        idx += 1

        state = states[i]
        parities = [(-1)^sum(state.occupations[:,k]) for k in 1:F]
        for j in 1:N
            for k in 1:F
                for dir in [-1,1]
                    !(1 ≤ j + dir ≤ N) && !(lattice.periodic) && continue
                    j2 = j + dir - (j + dir > N ? N : 0) + (j + dir < 1 ? N : 0)
                    if state.occupations[j,k] && !state.occupations[j2,k]
                        sign = 1

                        hopoccupations=copy(state.occupations)
                        hopL₀=L₀(state)
                        hopoccupations[j,k] = false
                        hopoccupations[j2,k] = true
                        if j==1 && dir==-1
                            hopL₀ += lattice.q
                            sign *= parities[k]
                        elseif j==N && dir==1
                            hopL₀ -= lattice.q
                            sign *= parities[k]
                        end

                        if haskey(positions, (hopoccupations,hopL₀))
                            I[idx] = i
                            J[idx] = positions[(hopoccupations,hopL₀)]
                            V[idx] = sign*(-(dir*1im)/(2*lattice.a) + (-1)^j*1im*lattice.mprime[j][k]/2)
                            idx += 1
                        end
                    end
                end
            end
        end
    end

    matrix = sparse(I[1:idx-1], J[1:idx-1], V[1:idx-1], length(states), length(states))
    return EDHamiltonian(lattice, matrix, L_max)
end

struct MPOHamiltonian{N,F} <: SchwingerHamiltonian{N,F}
    lattice::SchwingerLattice{N,F}
    mpo::MPO
    L_max::Int64

    function MPOHamiltonian(lattice::SchwingerLattice{N,F}, mpo::MPO, L_max::Int) where {N,F}
        new{N,F}(lattice, mpo, L_max)
    end
end

"""
`MPOhamiltonian(lattice)`

Computes the MPO Hamiltonian for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOHamiltonian(lattice::SchwingerLattice{N,F}; L_max::Int = 3, universe::Int = 0) where {N,F}
    @unpack q, periodic, a, θ2π, mlat, mprime = lattice

    universe = mod(universe, q)
    if universe < 0
        universe += q
    end
    θ2πu = θ2π .+ universe

    hamiltonian = OpSum()

    # Hopping term
    
    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += 1/(2a),"S+",ind,"S-",ind + F
            hamiltonian += 1/(2a),"S-",ind,"S+",ind + F
        end
    end

    if periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                hamiltonian += 1/(2a),"S+",ind,"S-",k,"raise",N * F + 1
                hamiltonian += 1/(2a),"S-",ind,"S+",k,"lower",N * F + 1
            else
                hamiltonian += 1/(2a),"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                hamiltonian += 1/(2a),"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    end

    # Real mass term

    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += mlat[j][k] * ((-1) ^ (j - 1)),"Sz",ind
        end
    end

    # Imaginary mass term

    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += (mprime[j][k]/2)*(-1)^(j-1),"S+",ind,"S-",ind + F
            hamiltonian += (mprime[j][k]/2)*(-1)^(j-1),"S-",ind,"S+",ind + F
        end
    end

    if periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                hamiltonian += -(mprime[N][k]/2),"S+",ind,"S-",k,"raise",N * F + 1
                hamiltonian += -(mprime[N][k]/2),"S-",ind,"S+",k,"lower",N * F + 1
            else
                hamiltonian += -(mprime[N][k]/2),"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                hamiltonian += -(mprime[N][k]/2),"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    end

    # Gauge kinetic term

    if periodic
        hamiltonian += q^2 * a * N / 2,"L0",N * F + 1,"L0",N * F + 1
        hamiltonian += q * a * sum(θ2πu),"L0",N * F + 1

        hamiltonian += q^2 * a * N * F / 4,"L0",N * F + 1
        
        for j in 1:N
            for k in 1:F
                ind = F*(j - 1) + k
                hamiltonian += q^2 * a * (N - j),"L0",N * F + 1,"Sz",ind
            end
        end
    end
    hamiltonian += a / 2 * sum(θ2πu .^ 2),"Id",1
    hamiltonian += q * a * F / 4 * sum(θ2πu),"Id",1
    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += q * a * sum(θ2πu[j+1:N]),"Sz",ind
        end
    end

    hamiltonian += q^2 * a * N * F * F / 16,"Id",1
    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += q^2 * a * F * ((N + 1 - j) ÷ 2) / 2,"Sz",ind
        end
    end

    for j1 in 1:N
        for k1 in 1:F
            ind1 = F * (j1 - 1) + k1
            for j2 in j1:N
                for k2 in (j2 == j1 ? (k1:F) : 1:F)
                    ind2 = F*(j2 - 1) + k2
                    hamiltonian += q^2*a*(N - j2)/(j1 == j2 && k1 == k2 ? 2 : 1),"Sz",ind1,"Sz",ind2
                end
            end
        end
    end

    mpo = MPO(hamiltonian, sites(lattice; L_max = L_max))
    return MPOHamiltonian(lattice, mpo, L_max)
end