"""
`EDMass(lattice)`
Computes the mass operator ∑ (-1)ⁿ χ†ₙχₙ for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDMass(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function massenergy(state::SchwingerBasisState{N,F})
        occs = occupations(state)
        return [(occs, L₀(state)) => sum((-1)^(j) * (bare ? 1 : state.lattice.mlat[j][k]) * occs[j,k] for j=1:N, k=1:F)]
    end
    return constructoperator(lattice, massenergy; L_max = L_max, universe = universe)
end

"""
`EDGaugeKinetic(lattice)`
Computes the gauge kinetic operator ∑(Lₙ+θ/2π)² for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDGaugeKinetic(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function electricenergy(state::SchwingerBasisState{N,F})
        efs = electricfields(state)
        return [(occupations(state), L₀(state)) => (bare ? 1 : lattice.a/2) * sum(efs .^ 2)]
    end
    return constructoperator(lattice, electricenergy; L_max = L_max, universe = universe)
end

"""
`EDHopping(lattice)`
Computes the hopping term -i ∑(χ†ₙ χₙ₊₁ - χ†ₙ₊₁ χₙ) for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDHopping(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function hopping(state::SchwingerBasisState{N,F})
        parities = [(-1)^sum(state.occupations[:,k]) for k in 1:F]
        hops = Dict{Tuple{BitMatrix,Int},ComplexF64}()
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

                       hops[(hopoccupations, hopL₀)] = sign*dir*1im*(bare ? 1 : 1/(2*lattice.a))
                    end
                end
            end
        end
        return hops
    end
    return constructoperator(lattice, hopping; L_max = L_max, universe = universe)
end

@memoize function EDHoppingMass(lattice::SchwingerLattice{N,F}, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function hoppingmass(state::SchwingerBasisState{N,F})
        parities = [(-1)^sum(state.occupations[:,k]) for k in 1:F]
        hops = Dict{Tuple{BitMatrix,Int},ComplexF64}()
        for (j, dir) in [(site, 1), (mod(site, N) + 1, -1)]
            for k in 1:F
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

                    hops[(hopoccupations, hopL₀)] = sign*(-1)^(j+1)*1im*(bare ? 1 : lattice.mprime[j][k])
                end
            end
        end
        return hops
    end
    return constructoperator(lattice, hoppingmass; L_max = L_max, universe = universe)
end

"""
`EDHoppingMass(lattice)`
Computes the hopping-type mass term i/2 ∑(-1)^n (χ†ₙ₊₁ χₙ + χ†ₙ₋₁ χₙ) for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDHoppingMass(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDHoppingMass.(Ref(lattice), 1:N; L_max = L_max, universe = universe, bare = bare))
end

"""
`EDHamiltonian(lattice)`
Computes the Hamiltonian for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function EDHamiltonian(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return  EDGaugeKinetic(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDMass(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDHopping(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDHoppingMass(lattice; L_max = L_max, universe = universe, bare = false)
end

function opsum_gaugekinetic(lattice::SchwingerLattice{N,F}; universe::Int = 0, bare::Bool = true) where {N,F}
    @unpack q, periodic, a, θ2π = lattice

    θ2πu = θ2π .+ universe

    term = OpSum()

    if periodic
        term += q^2 * N,"L0",N * F + 1,"L0",N * F + 1
        term += 2q * sum(θ2πu),"L0",N * F + 1

        term += -q^2 * N * F / 2,"L0",N * F + 1
        
        for j in 1:N
            for k in 1:F
                ind = F*(j - 1) + k
                term += 2q^2 * (N + 1 - j),"L0",N * F + 1,"Sz",ind
            end
        end
    end

    term += sum(θ2πu .^ 2),"Id",1
    term += -q * F / 2 * sum(θ2πu),"Id",1
    term += q^2 * N * F * F / 8,"Id",1

    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            term += 2q * sum(θ2πu[j+1:N]),"Sz",ind
        end
    end

    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            term += -q^2 * F * ((N + 1 - j) ÷ 2),"Sz",ind
        end
    end

    for j1 in 1:N
        for k1 in 1:F
            ind1 = F * (j1 - 1) + k1
            for j2 in j1:N
                for k2 in (j2 == j1 ? (k1:F) : 1:F)
                    ind2 = F*(j2 - 1) + k2
                    term += 2q^2*(N - j2)/(j1 == j2 && k1 == k2 ? 2 : 1),"Sz",ind1,"Sz",ind2
                end
            end
        end
    end

    return (bare ? 1 : a/2)*term
end

"""
`MPOGaugeKinetic(lattice)`

Computes the MPO gauge kinetic operator for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOGaugeKinetic(lattice::SchwingerLattice{N,F}; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_gaugekinetic(lattice; universe = universe, bare = bare)
    mpo = MPO(opsum, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

function opsum_mass(lattice::SchwingerLattice{N,F}; bare::Bool = true) where {N,F}
    @unpack mlat = lattice
    term = OpSum()

    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            term += (bare ? 1 : mlat[j][k]) * ((-1) ^ j),"Sz",ind
        end
    end

    return term
end

"""
`MPOMass(lattice)`

Computes the MPO mass operator for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOMass(lattice::SchwingerLattice{N,F}; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_mass(lattice; bare = bare)
    mpo = MPO(opsum, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

function opsum_hopping(lattice::SchwingerLattice{N,F}; bare::Bool = true) where {N,F}
    @unpack a, periodic = lattice
    term = OpSum()

    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            term += "S+",ind,"S-",ind + F
            term += "S-",ind,"S+",ind + F
        end
    end

    if periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                term += "S+",ind,"S-",k,"raise",N * F + 1
                term += "S-",ind,"S+",k,"lower",N * F + 1
            else
                term += 2^N,"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                term += 2^N,"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    end

    return (bare ? 1 : 1/(2a)) * term
end

"""
`MPOHopping(lattice)`

Computes the MPO hopping operator for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOHopping(lattice::SchwingerLattice{N,F}; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hopping(lattice; bare = bare)
    mpo = MPO(opsum, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

function opsum_hoppingmass(lattice::SchwingerLattice{N,F}, site::Int; bare::Bool = false) where {N,F}
    if !(1 ≤ site ≤ N)
        throw(ArgumentError("site must be between 1 and N"))
    end

    term = OpSum()

    if site < N
        for k in 1:F
            ind = F*(site - 1) + k
            term += (bare ? 1 : lattice.mprime[site][k])*(-1)^(site+1),"S+",ind,"S-",ind + F
            term += (bare ? 1 : lattice.mprime[site][k])*(-1)^(site+1),"S-",ind,"S+",ind + F
        end
    elseif lattice.periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                term += -(bare ? 1 : lattice.mprime[N][k]),"S+",ind,"S-",k,"raise",N * F + 1
                term += -(bare ? 1 : lattice.mprime[N][k]),"S-",ind,"S+",k,"lower",N * F + 1
            else
                term += -(2^N)*(bare ? 1 : lattice.mprime[N][k]),"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                term += -(2^N)*(bare ? 1 : lattice.mprime[N][k]),"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    else
        term += 0,"Id",1
    end

    return term
end

function opsum_hoppingmass(lattice::SchwingerLattice{N,F}; bare::Bool = false) where {N,F}
    sum(opsum_hoppingmass.(Ref(lattice), 1:N; bare = bare))
end

@memoize function MPOHoppingMass(lattice::SchwingerLattice{N,F}, site::Int; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hoppingmass(lattice, site; bare = bare)
    mpo = MPO(opsum, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

"""
`MPOHoppingMass(lattice)`

Computes the MPO hopping-mass operator for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOHoppingMass(lattice::SchwingerLattice{N,F}; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hoppingmass(lattice; bare = bare)
    mpo = MPO(opsum, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

"""
`MPOHamiltonian(lattice)`

Computes the MPO Hamiltonian for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function MPOHamiltonian(lattice::SchwingerLattice{N,F}; L_max::Union{Int,Nothing} = nothing, universe::Int = 0) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    hamiltonian =   opsum_gaugekinetic(lattice; bare = false, universe = universe) + 
                    opsum_mass(lattice; bare = false) + 
                    opsum_hopping(lattice; bare = false) + 
                    opsum_hoppingmass(lattice; bare = false)

    mpo = MPO(hamiltonian, sites(lattice; L_max = L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end