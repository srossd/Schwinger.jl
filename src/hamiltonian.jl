# =============================================================================
# ED Backend
# =============================================================================

"""
`EDMass(lattice)`
Computes the mass operator ∑ (-1)ⁿ χ†ₙχₙ for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDMass(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true)
    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function massenergy(state::BasisState)
        occs = occupations(state)
        return [(occs, L₀(state)) => sum((-1)^(j) * (bare ? 1 : state.lattice.mlat[j][k]) * occs[j,k] for j=1:N, k=1:F)]
    end
    return constructoperator(lattice, massenergy; L_max = L_max, universe = universe)
end

"""
`EDGaugeKinetic(lattice)`
Computes the gauge kinetic operator ∑(Lₙ+θ/2π)² for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDGaugeKinetic(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function electricenergy(state::BasisState)
        efs = electricfields(state)
        return [(occupations(state), L₀(state)) => (bare ? 1 : lattice.a/2) * sum(efs .^ 2)]
    end
    return constructoperator(lattice, electricenergy; L_max = L_max, universe = universe)
end

"""
`EDHopping(lattice)`
Computes the hopping term -i ∑(χ†ₙ χₙ₊₁ - χ†ₙ₊₁ χₙ) for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHopping(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true)
    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function hopping(state::BasisState)
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

@memoize function EDHoppingMass(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true)
    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function hoppingmass(state::BasisState)
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
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHoppingMass(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDHoppingMass.(Ref(lattice), 1:Int(lattice.N); L_max = L_max, universe = universe, bare = bare))
end

"""
`EDHamiltonian(lattice)`
Computes the Hamiltonian for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHamiltonian(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return  EDGaugeKinetic(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDMass(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDHopping(lattice; L_max = L_max, universe = universe, bare = false) + 
            EDHoppingMass(lattice; L_max = L_max, universe = universe, bare = false)
end

# =============================================================================
# ITensors Backend
# =============================================================================

function opsum_gaugekinetic(lattice::Lattice; universe::Int = 0, bare::Bool = true)
    N, F = Int(lattice.N), lattice.F
    q, periodic, a, θ2π = lattice.q, lattice.periodic, lattice.a, lattice.θ2π

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

    term += sum(θ2πu[1:N] .^ 2),"Id",1
    term += -q * F * sum(θ2πu[1:N]) / 2,"Id",1
    term += q^2 * N * F^2 / 8,"Id",1

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
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function ITensorGaugeKinetic(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    N, F = lattice.N, lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_gaugekinetic(lattice; universe = universe, bare = bare)
    mpo = ITensorMPS.MPO(opsum, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

function opsum_mass(lattice::Lattice; bare::Bool = true)
    N, F, mlat = Int(lattice.N), lattice.F, lattice.mlat
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
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function ITensorMass(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    N, F = lattice.N, lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_mass(lattice; bare = bare)
    mpo = MPO(opsum, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

function opsum_hopping(lattice::Lattice; bare::Bool = true)
    N, F, a, periodic = Int(lattice.N), lattice.F, lattice.a, lattice.periodic
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
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function ITensorHopping(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    N, F = lattice.N, lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hopping(lattice; bare = bare)
    mpo = ITensorMPS.MPO(opsum, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

function opsum_hoppingmass(lattice::Lattice, site::Int; bare::Bool = false)
    N, F = Int(lattice.N), lattice.F
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

function opsum_hoppingmass(lattice::Lattice; bare::Bool = false)
    N, F = Int(lattice.N), lattice.F
    sum(opsum_hoppingmass.(Ref(lattice), 1:N; bare = bare))
end

@memoize function ITensorHoppingMass(lattice::Lattice, site::Int; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hoppingmass(lattice, site; bare = bare)
    mpo = ITensorMPS.MPO(opsum, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

"""
`MPOHoppingMass(lattice)`

Computes the MPO hopping-mass operator for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function ITensorHoppingMass(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    opsum = opsum_hoppingmass(lattice; bare = bare)
    mpo = ITensorMPS.MPO(opsum, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

"""
`MPOHamiltonian(lattice)`

Computes the MPO Hamiltonian for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function ITensorHamiltonian(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, sites::Union{Vector{Index{T}},Nothing} = nothing) where {T}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    hamiltonian =   opsum_gaugekinetic(lattice; bare = false, universe = universe) +
                    opsum_mass(lattice; bare = false) +
                    opsum_hopping(lattice; bare = false) +
                    opsum_hoppingmass(lattice; bare = false)

    mpo = ITensorMPS.MPO(hamiltonian, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end

# =============================================================================
# MPSKit Backend
# =============================================================================

"""
`MPSKitGaugeKinetic(lattice)`

Computes the MPSKit gauge kinetic operator for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitGaugeKinetic(lattice::Lattice; universe::Int = 0)
    if lattice.periodic
        throw(ArgumentError("MPSKit backend not implemented for periodic lattices"))
    end
    N, F = lattice.N, lattice.F
    _, universe = process_L_max_universe(lattice, nothing, universe)
    θ2πu = copy(lattice.θ2π) .+ universe

    link_fcts = [r::U1Irrep -> (lattice.a/2) * (r.charge + θ2πu[n])^2 for n in 1:(isinf(N) ? 2 : Int(N))]
    if isinf(N)
        mpo = InfiniteLEMPOHamiltonian(get_mpskit_spaces(lattice), vcat(fill(missing, F - 1), [link_fcts[1]], fill(missing, F - 1), [link_fcts[2]]))
    else
        all_link_fcts = vcat([vcat(fill(missing, F - 1), [link_fcts[n]]) for n in 1:Int(N)]...)
        mpo = FiniteLEMPOHamiltonian(get_mpskit_spaces(lattice), all_link_fcts)
    end
    return MPSKitOperator(lattice, mpo, universe)
end

function matrices_mass(lattice::Lattice)
    if lattice.periodic
        throw(ArgumentError("MPSKit backend not implemented for periodic lattices"))
    end
    
    N, F = lattice.N, lattice.F

    spaces = get_mpskit_spaces(lattice)
    oddspace = spaces[1]
    evenspace = spaces[F + 1]

    masseven = ones(U1Space(0 => 1) ⊗ evenspace ← evenspace ⊗ U1Space(0 => 1))
    block(masseven, U1Irrep(0)) .= -0.5
    block(masseven, U1Irrep(lattice.q)) .= 0.5
    massodd  = ones(U1Space(0 => 1) ⊗ oddspace  ← oddspace  ⊗ U1Space(0 => 1))
    block(massodd, U1Irrep(0)) .= -0.5
    block(massodd, U1Irrep(-lattice.q)) .= 0.5

    Elt = Union{Missing, typeof(masseven), scalartype(masseven)}
    A = Vector{Matrix{Elt}}(undef, (isinf(N) ? 2 : Int(N))*F)

    for n in 1:(isinf(N) ? 2 : Int(N))*F
        W = Matrix{Elt}(missing, 2, 2)
        W[1, 1] = 1.0
        W[end, end] = 1.0

        site_ind = 1 + ((n - 1) ÷ F)
        flavor_ind = mod(n - 1, F) + 1

        if site_ind % 2 == 1
            W[1, 2] = lattice.mlat[site_ind][flavor_ind] * massodd
        else
            W[1, 2] = lattice.mlat[site_ind][flavor_ind] * masseven
        end

        A[n] = W
    end

    return A
end

"""
`MPSKitMass(lattice)`

Computes the MPSKit mass operator for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitMass(lattice::Lattice; universe::Int = 0, bare::Bool = true)
    A = matrices_mass(lattice)

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        A[1] = A[1][1:1, :]
        A[end] = A[end][:, end:end]
        mpo = FiniteMPOHamiltonian(A)
    end

    if !bare
        mpo *= lattice.mlat
    end

    return MPSKitOperator(lattice, mpo, universe)
end

function matrices_hopping(lattice::Lattice)
    if lattice.periodic
        throw(ArgumentError("MPSKit backend not implemented for periodic lattices"))
    end
    
    N, F = lattice.N, lattice.F

    spaces = get_mpskit_spaces(lattice)
    oddspace = spaces[1]
    evenspace = spaces[F + 1]

    q = lattice.q
    hopleftevenL =  ones(U1Space(0 => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(q => 1))
    hopleftevenR =  ones(U1Space(q => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(0 => 1))
    hopleftoddL =   ones(U1Space(0 => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(q => 1))
    hopleftoddR =   ones(U1Space(q => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(0 => 1))
    hoprightevenL = ones(U1Space(0 => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(-q => 1))
    hoprightevenR = ones(U1Space(-q => 1) ⊗ evenspace ← evenspace ⊗ U1Space(0 => 1))
    hoprightoddL =  ones(U1Space(0 => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(-q => 1))
    hoprightoddR =  ones(U1Space(-q => 1) ⊗ oddspace  ← oddspace  ⊗ U1Space(0 => 1))

    Elt = Union{Missing, typeof(hopleftevenL), scalartype(hopleftevenL)}
    A = Vector{Matrix{Elt}}(undef, (isinf(N) ? 2 : Int(N))*F)

    for n in 1:(isinf(N) ? 2 : Int(N))*F
        W = Matrix{Elt}(missing, 2 + 2*F, 2 + 2*F)
        W[1, 1] = 1.0
        W[end, end] = 1.0

        site_ind = 1 + ((n - 1) ÷ F)
        flavor_ind = mod(n - 1, F) + 1

        if isodd(site_ind)
            W[1, 1 + 2*(flavor_ind - 1) + 1] = 1/(2*lattice.a) * hopleftoddL
            W[1, 1 + 2*(flavor_ind - 1) + 2] = 1/(2*lattice.a) * hoprightoddL
            W[1 + 2*(flavor_ind - 1) + 1, end] = hopleftoddR
            W[1 + 2*(flavor_ind - 1) + 2, end] = hoprightoddR

            for j=1:F
                j == flavor_ind && continue
                W[1 + 2*(j - 1) + 1, 1 + 2*(j - 1) + 1] = 1.0
                W[1 + 2*(j - 1) + 2, 1 + 2*(j - 1) + 2] = 1.0
            end
        else
            W[1, 1 + 2*(flavor_ind - 1) + 1] = 1/(2*lattice.a) * hopleftevenL
            W[1, 1 + 2*(flavor_ind - 1) + 2] = 1/(2*lattice.a) * hoprightevenL
            W[1 + 2*(flavor_ind - 1) + 1, end] = hopleftevenR
            W[1 + 2*(flavor_ind - 1) + 2, end] = hoprightevenR

            for j=1:F
                j == flavor_ind && continue
                W[1 + 2*(j - 1) + 1, 1 + 2*(j - 1) + 1] = 1.0
                W[1 + 2*(j - 1) + 2, 1 + 2*(j - 1) + 2] = 1.0
            end
        end

        A[n] = W
    end

    return A
end

"""
`MPSKitHopping(lattice)`

Computes the MPSKit hopping operator for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitHopping(lattice::Lattice; universe::Int = 0)
    A = matrices_hopping(lattice)

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        A[1] = A[1][1:1, :]
        A[end] = A[end][:, end:end]
        mpo = FiniteMPOHamiltonian(A)
    end

    return MPSKitOperator(lattice, mpo, universe)
end

function matrices_hoppingmass(lattice::Lattice; bare::Bool = false)
    if lattice.periodic
        throw(ArgumentError("MPSKit backend not implemented for periodic lattices"))
    end
    
    N, F = lattice.N, lattice.F

    spaces = get_mpskit_spaces(lattice)
    oddspace = spaces[1]
    evenspace = spaces[F + 1]

    q = lattice.q
    hopleftevenL =  ones(U1Space(0 => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(q => 1))
    hopleftevenR =  ones(U1Space(q => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(0 => 1))
    hopleftoddL =   ones(U1Space(0 => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(q => 1))
    hopleftoddR =   ones(U1Space(q => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(0 => 1))
    hoprightevenL = ones(U1Space(0 => 1)  ⊗ evenspace ← evenspace ⊗ U1Space(-q => 1))
    hoprightevenR = ones(U1Space(-q => 1) ⊗ evenspace ← evenspace ⊗ U1Space(0 => 1))
    hoprightoddL =  ones(U1Space(0 => 1)  ⊗ oddspace  ← oddspace  ⊗ U1Space(-q => 1))
    hoprightoddR =  ones(U1Space(-q => 1) ⊗ oddspace  ← oddspace  ⊗ U1Space(0 => 1))

    Elt = Union{Missing, typeof(hopleftevenL), scalartype(hopleftevenL)}
    A = Vector{Matrix{Elt}}(undef, (isinf(N) ? 2 : Int(N))*F)

    for n in 1:(isinf(N) ? 2 : Int(N))*F
        W = Matrix{Elt}(missing, 2 + 2*F, 2 + 2*F)
        W[1, 1] = 1.0
        W[end, end] = 1.0

        site_ind = 1 + ((n - 1) ÷ F)
        flavor_ind = mod(n - 1, F) + 1

        if isodd(site_ind)
            W[1, 1 + 2*(flavor_ind - 1) + 1] = (bare ? 1 : lattice.mprime[site_ind][flavor_ind]) * hopleftoddL
            W[1, 1 + 2*(flavor_ind - 1) + 2] = (bare ? 1 : lattice.mprime[site_ind][flavor_ind]) * hoprightoddL
            W[1 + 2*(flavor_ind - 1) + 1, end] = hopleftoddR
            W[1 + 2*(flavor_ind - 1) + 2, end] = hoprightoddR

            for j=1:F
                j == flavor_ind && continue
                W[1 + 2*(j - 1) + 1, 1 + 2*(j - 1) + 1] = 1.0
                W[1 + 2*(j - 1) + 2, 1 + 2*(j - 1) + 2] = 1.0
            end
        else
            W[1, 1 + 2*(flavor_ind - 1) + 1] = -(bare ? 1 : lattice.mprime[site_ind][flavor_ind]) * hopleftevenL
            W[1, 1 + 2*(flavor_ind - 1) + 2] = -(bare ? 1 : lattice.mprime[site_ind][flavor_ind]) * hoprightevenL
            W[1 + 2*(flavor_ind - 1) + 1, end] = hopleftevenR
            W[1 + 2*(flavor_ind - 1) + 2, end] = hoprightevenR

            for j=1:F
                j == flavor_ind && continue
                W[1 + 2*(j - 1) + 1, 1 + 2*(j - 1) + 1] = 1.0
                W[1 + 2*(j - 1) + 2, 1 + 2*(j - 1) + 2] = 1.0
            end
        end

        A[n] = W
    end

    return A
end

@memoize function MPSKitHoppingMass(lattice::Lattice, site::Int; bare::Bool = true, universe::Int = 0)
    A = matrices_hoppingmass(lattice; bare = bare) # TODO: optimize to only build needed tensors

    for n in eachindex(A)
        site_ind = 1 + ((n - 1) ÷ lattice.F)
        if site_ind != site && site_ind != site + 1
            for m in 2:size(A[n])[1]-1
                A[n][1, m] *= 0.0
                A[n][m, end] *= 0.0
            end
        end
    end

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        A[1] = A[1][1:1, :]
        for j=2:lattice.F
            for k=2:(2*lattice.F + 1)
                A[j][k, end] = missing
            end
        end
        for j=(Int(lattice.N) - 1)*lattice.F + 1:Int(lattice.N)*lattice.F - 1
            for k=2:(2*lattice.F + 1)
                A[j][1, k] = missing
            end
        end
        A[end] = A[end][:, end:end]
        mpo = FiniteMPOHamiltonian(A)
    end

    return MPSKitOperator(lattice, mpo, universe)
end

"""
`MPSKitHoppingMass(lattice)`

Computes the MPSKit hopping-mass operator for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitHoppingMass(lattice::Lattice; universe::Int = 0)
    A = matrices_hoppingmass(lattice)

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        A[1] = A[1][1:1, :]
        for j=2:lattice.F
            for k=2:(2*lattice.F + 1)
                A[j][k, end] = missing
            end
        end
        for j=(Int(lattice.N) - 1)*lattice.F + 1:Int(lattice.N)*lattice.F - 1
            for k=2:(2*lattice.F + 1)
                A[j][1, k] = missing
            end
        end
        A[end] = A[end][:, end:end]
        mpo = FiniteMPOHamiltonian(A)
    end

    return MPSKitOperator(lattice, mpo, universe)
end

"""
`MPSKitHamiltonian(lattice)`

Computes the MPSKit Hamiltonian for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitHamiltonian(lattice::Lattice; universe::Int = 0)
    Amass = matrices_mass(lattice)
    Ahoppingmass = matrices_hoppingmass(lattice)
    A = matrices_hopping(lattice)
    gauge = MPSKitGaugeKinetic(lattice; universe = universe)

    for n in eachindex(A)
        A[n][1, end] = Amass[n][1, end]
        for k in 1:lattice.F
            A[n][1, 1 + 2*(k - 1) + 1] += Ahoppingmass[n][1, 1 + 2*(k - 1) + 1]
            A[n][1, 1 + 2*(k - 1) + 2] += Ahoppingmass[n][1, 1 + 2*(k - 1) + 2]
        end
    end

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        A[1] = A[1][1:1, :]
        for j=2:lattice.F
            for k=2:(2*lattice.F + 1)
                A[j][k, end] = missing
            end
        end
        for j=length(A) - lattice.F + 1:length(A) - 1
            for k=2:(2*lattice.F + 1)
                A[j][1, k] = missing
            end
        end
        A[end] = A[end][:, end:end]
        mpo = FiniteMPOHamiltonian(A)
    end

    lempo = if isinf(lattice.N)
        InfiniteLEMPOHamiltonian(mpo, gauge.lempo.link_fcts)
    else
        FiniteLEMPOHamiltonian(mpo, gauge.lempo.link_fcts)
    end

    return MPSKitOperator(lattice, lempo, universe)
end