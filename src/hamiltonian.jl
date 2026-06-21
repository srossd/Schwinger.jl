# =============================================================================
# ED Backend
# =============================================================================

"""
`EDMass(lattice)`
Computes the mass operator ∑ (-1)ⁿ χ†ₙχₙ for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDMass(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    F = lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function massenergy(state::BasisState)
        occs = occupations(state)
        return [(occs, L₀(state)) => sum((-1)^(site) * (bare ? 1 : state.lattice.mlat[site][k]) * occs[site,k] for k=1:F)]
    end
    return constructoperator(lattice, massenergy; L_max = L_max, universe = universe, in_charge = charge, out_charge = charge)
end

@memoize function EDMass(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDMass.(Ref(lattice), 1:Int(lattice.N); L_max = L_max, universe = universe, bare = bare, charge = charge))
end

"""
`EDGaugeKinetic(lattice)`
Computes the gauge kinetic operator ∑(Lₙ+θ/2π)² for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDGaugeKinetic(lattice::Lattice, link::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function electricenergy(state::BasisState)
        return [(occupations(state), L₀(state)) => (bare ? 1 : lattice.a/2) * electricfields(state)[link]^2]
    end
    return constructoperator(lattice, electricenergy; L_max = L_max, universe = universe, in_charge = charge, out_charge = charge)
end

@memoize function EDGaugeKinetic(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDGaugeKinetic.(Ref(lattice), 1:Int(lattice.N); L_max = L_max, universe = universe, bare = bare, charge = charge))
end

"""
`EDHopping(lattice)`
Computes the hopping term -i ∑(χ†ₙ χₙ₊₁ - χ†ₙ₊₁ χₙ) for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHopping(lattice::Lattice, bond::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    function hopping(state::BasisState)
        parities = [(-1)^sum(state.occupations[:,k]) for k in 1:F]
        hops = Dict{Tuple{BitMatrix,Int},ComplexF64}()
        for (j, dir) in [(bond, 1), (mod(bond, N) + 1, -1)]   # the two directed hops on bond (bond, bond+1)
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

                    hops[(hopoccupations, hopL₀)] = sign*dir*1im*(bare ? 1 : 1/(2*lattice.a))
                end
            end
        end
        return hops
    end
    return constructoperator(lattice, hopping; L_max = L_max, universe = universe, in_charge = charge, out_charge = charge)
end

@memoize function EDHopping(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDHopping.(Ref(lattice), 1:Int(lattice.N); L_max = L_max, universe = universe, bare = bare, charge = charge))
end

@memoize function EDHoppingMass(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
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
    return constructoperator(lattice, hoppingmass; L_max = L_max, universe = universe, in_charge = charge, out_charge = charge)
end

"""
`EDHoppingMass(lattice)`
Computes the hopping-type mass term i/2 ∑(-1)^n (χ†ₙ₊₁ χₙ + χ†ₙ₋₁ χₙ) for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHoppingMass(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0, bare::Bool = true, charge::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return sum(EDHoppingMass.(Ref(lattice), 1:Int(lattice.N); L_max = L_max, universe = universe, bare = bare, charge = charge))
end

"""
`EDHamiltonian(lattice)`
Computes the Hamiltonian for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function EDHamiltonian(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0,
                                charge::Int = 0, defects::Vector{DefectCharge} = DefectCharge[])
    if !isempty(defects)   # a defect charge is a θ2π step internally, but we keep the
        # original (unshifted) lattice on the operator and record the defect list. The
        # defect's charge lives in the flux (θ shift), NOT in the matter filling, so the
        # matter sector is `charge` (default 0, neutral) — the same one acted on by the
        # charge-changing `FermionField` when building a Wilson line.
        op = EDHamiltonian(_lattice_with_defects(lattice, defects); L_max = L_max, universe = universe, charge = charge)
        return EDOperator(lattice, op.matrix, op.L_max, op.universe, charge, charge, defects)
    end
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    return  EDGaugeKinetic(lattice; L_max = L_max, universe = universe, bare = false, charge = charge) +
            EDMass(lattice; L_max = L_max, universe = universe, bare = false, charge = charge) +
            EDHopping(lattice; L_max = L_max, universe = universe, bare = false, charge = charge) +
            EDHoppingMass(lattice; L_max = L_max, universe = universe, bare = false, charge = charge)
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
    term += -q * F * sum(θ2πu[1:2:N]),"Id",1   # 2·Σₙ θₙ Gₙ  (Gₙ ∝ Σ_{m≤n}(-1)ᵐ)
    term += q^2 * N * F^2 / 8,"Id",1

    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            term += 2q * sum(θ2πu[j:N]),"Sz",ind   # Sz_j couples to θ on links n ≥ j
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

# ---- per-link / per-site / per-bond ITensor term operators (used by `energy_density`) ----

# electric energy on a single `link`: (a/2)⟨E_link²⟩, reusing AverageElectricField (power 2)
@memoize function ITensorGaugeKinetic(lattice::Lattice, link::Int; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true)
    op = ITensorAverageElectricField(lattice; power = 2, sitelist = [link], L_max = L_max, universe = universe)
    return ITensorOperator(lattice, (bare ? 1.0 : lattice.a / 2) * op.mpo, op.L_max, op.universe)
end

function opsum_mass(lattice::Lattice, site::Int; bare::Bool = true)
    term = OpSum()
    for k in 1:lattice.F
        ind = lattice.F*(site - 1) + k
        term += (bare ? 1 : lattice.mlat[site][k]) * ((-1)^site),"Sz",ind
    end
    return term
end
@memoize function ITensorMass(lattice::Lattice, site::Int; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    mpo = ITensorMPS.MPO(opsum_mass(lattice, site; bare = bare), get_sites(lattice; L_max = L_max))
    return ITensorOperator(lattice, mpo, L_max, universe)
end

function opsum_hopping(lattice::Lattice, bond::Int; bare::Bool = true)
    N, F, a = Int(lattice.N), lattice.F, lattice.a
    term = OpSum()
    if bond < N
        for k in 1:F
            ind = F*(bond - 1) + k
            term += "S+",ind,"S-",ind + F
            term += "S-",ind,"S+",ind + F
        end
    elseif lattice.periodic
        for k in 1:F
            ind = F*(N - 1) + k
            if F ≤ 2
                term += "S+",ind,"S-",k,"raise",N*F + 1
                term += "S-",ind,"S+",k,"lower",N*F + 1
            else
                term += 2^N,"S+",ind,"S-",k,"raise",N*F + 1, wigner_string(N, F, k)...
                term += 2^N,"S-",ind,"S+",k,"lower",N*F + 1, wigner_string(N, F, k)...
            end
        end
    else
        term += 0,"Id",1
    end
    return (bare ? 1 : 1/(2a)) * term
end
@memoize function ITensorHopping(lattice::Lattice, bond::Int; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, bare::Bool = true)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    mpo = ITensorMPS.MPO(opsum_hopping(lattice, bond; bare = bare), get_sites(lattice; L_max = L_max))
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
@memoize function ITensorHamiltonian(lattice::Lattice; L_max::Union{Int,Nothing} = nothing, universe::Int = 0, sites::Union{AbstractVector{<:Index},Nothing} = nothing, defects::Vector{DefectCharge} = DefectCharge[])
    if !isempty(defects)   # θ2π step internally; keep the original lattice + defect list
        L_max_p, _ = process_L_max_universe(lattice, L_max, universe)
        # build the (θ-shifted) MPO on the ORIGINAL lattice's site indices so that the
        # operator stays compatible with states built from the original lattice
        sites_orig = isnothing(sites) ? get_sites(lattice; L_max = L_max_p) : sites
        op = ITensorHamiltonian(_lattice_with_defects(lattice, defects); L_max = L_max, universe = universe, sites = sites_orig)
        return ITensorOperator(lattice, op.mpo, op.L_max, op.universe, defects)
    end
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

# ---- per-site / per-link / per-bond MPSKit term operators (used by `energy_density`) ----

"""
`MPSKitGaugeKinetic(lattice, link)`

Electric energy `(a/2)(Lₙ+θ/2π)²` on a single `link`, as a LEMPO whose link function is
nonzero only on `link` (the charge still accumulates in the virtual bond).
"""
@memoize function MPSKitGaugeKinetic(lattice::Lattice, link::Int; universe::Int = 0)
    isinf(lattice.N) && throw(ArgumentError("per-link MPSKitGaugeKinetic requires a finite lattice"))
    N, F = Int(lattice.N), lattice.F
    _, universe = process_L_max_universe(lattice, nothing, universe)
    θ2πu = Base.convert(Vector{Float64}, lattice.θ2π[1:N]) .+ universe
    fct = r::U1Irrep -> (lattice.a / 2) * (r.charge + θ2πu[link])^2
    fcts = Union{Missing,Function}[k == link * F ? fct : missing for k in 1:N*F]
    return MPSKitOperator(lattice, FiniteLEMPOHamiltonian(get_mpskit_spaces(lattice), fcts), universe)
end

"""
`MPSKitMass(lattice, site)`

Staggered mass term on a single `site` (summed over flavors).
"""
@memoize function MPSKitMass(lattice::Lattice, site::Int; universe::Int = 0, bare::Bool = true)
    isinf(lattice.N) && throw(ArgumentError("per-site MPSKitMass requires a finite lattice"))
    F = lattice.F
    _, universe = process_L_max_universe(lattice, nothing, universe)
    sp = collect(get_mpskit_spaces(lattice))
    function massop(P)                                  # diagonal: empty −½, occupied +½
        T = zeros(ComplexF64, P ← P)
        block(T, U1Irrep(0)) .= -0.5
        block(T, U1Irrep(isodd(site) ? -lattice.q : lattice.q)) .= 0.5
        return T
    end
    terms = [(site - 1)*F + f => ComplexF64(bare ? 1 : lattice.mlat[site][f]) * massop(sp[(site - 1)*F + f]) for f in 1:F]
    return MPSKitOperator(lattice, FiniteMPOHamiltonian(sp, terms...), universe)
end

"""
`MPSKitHopping(lattice, bond)`

Hopping term `1/(2a)(χ†ₙχₙ₊₁ + h.c.)` on a single `bond` `(bond, bond+1)`, summed over
flavors and the two charge-transfer directions.
"""
@memoize function MPSKitHopping(lattice::Lattice, bond::Int; universe::Int = 0, bare::Bool = true)
    isinf(lattice.N) && throw(ArgumentError("per-bond MPSKitHopping requires a finite lattice"))
    N, F = Int(lattice.N), lattice.F
    _, universe = process_L_max_universe(lattice, nothing, universe)
    sp = collect(get_mpskit_spaces(lattice))
    idop(P) = isomorphism(ComplexF64, U1Space(0 => 1) ⊗ P, P ⊗ U1Space(0 => 1))
    coeff = bare ? 1.0 : 1 / (2 * lattice.a)
    ops = MPSKitOperator[]
    for f in 1:F, qs in (lattice.q, -lattice.q)
        i1, i2 = (bond - 1)*F + f, bond*F + f
        ts = Any[idop(sp[k]) for k in 1:N*F]
        for k in (i1 + 1):(i2 - 1)                       # carry the charge across intervening sites
            ts[k] = Base.convert(TensorMap, BraidingTensor(sp[k], U1Space(qs => 1)))
        end
        ts[i1] = coeff * ones(ComplexF64, U1Space(0 => 1) ⊗ sp[i1] ← sp[i1] ⊗ U1Space(qs => 1))
        ts[i2] = ones(ComplexF64, U1Space(qs => 1) ⊗ sp[i2] ← sp[i2] ⊗ U1Space(0 => 1))
        push!(ops, MPSKitOperator(lattice, FiniteMPO([t for t in ts]), universe))
    end
    return ops
end

# Combined matter MPO channel-matrices (hopping + mass + hopping-mass), before any defect
# handling or boundary trimming, together with the gauge operator (for its link functions).
# Shared by `MPSKitHamiltonian` and the fused-defect LEMPO builder.
function _matter_mpo_matrices(lattice::Lattice; universe::Int = 0)
    Amass = matrices_mass(lattice)
    Ahoppingmass = matrices_hoppingmass(lattice)
    A = matrices_hopping(lattice)
    for n in eachindex(A)
        A[n][1, end] = Amass[n][1, end]
        for k in 1:lattice.F
            A[n][1, 1 + 2*(k - 1) + 1] += Ahoppingmass[n][1, 1 + 2*(k - 1) + 1]
            A[n][1, 1 + 2*(k - 1) + 2] += Ahoppingmass[n][1, 1 + 2*(k - 1) + 2]
        end
    end
    return A, MPSKitGaugeKinetic(lattice; universe = universe)
end

# Trim the boundary channels of a finite matter MPO channel-matrix list (in place).
function _trim_finite_mpo!(A, F::Int)
    A[1] = A[1][1:1, :]
    for j in 2:F, k in 2:(2F + 1); A[j][k, end] = missing; end
    for j in (length(A) - F + 1):(length(A) - 1), k in 2:(2F + 1); A[j][1, k] = missing; end
    A[end] = A[end][:, end:end]
    return A
end

"""
`MPSKitHamiltonian(lattice)`

Computes the MPSKit Hamiltonian for the Schwinger model.

# Arguments
- `lattice::Lattice`: Schwinger model lattice.
"""
@memoize function MPSKitHamiltonian(lattice::Lattice; universe::Int = 0,
                                    defects::Vector{DefectCharge} = DefectCharge[])
    isinf(lattice.N) && !isempty(defects) &&
        throw(ArgumentError("defects are not supported for infinite lattices"))

    A, gauge = _matter_mpo_matrices(lattice; universe = universe)

    # Insert a static defect as a genuine extra lattice site: an inert MPO site
    # that passes every channel through as identity (so the matter Hamiltonian
    # acts across it unchanged) on a 1-D `U1Space(charge=>1)` physical leg, plus a
    # zero-width (`missing`) gauge link after it.  The charge enters the LEMPO's
    # Gauss-law accumulation automatically, shifting all links to its right.
    link_fcts = isempty(defects) ? gauge.lempo.link_fcts : collect(Any, gauge.lempo.link_fcts)
    if !isempty(defects)
        Wsize = size(A[1], 1)                       # 2 + 2F MPO channels
        # virtual (bond) space carried by each MPO channel: identity-in/out and
        # the Hamiltonian accumulator carry charge 0; hop-left/right carry ±q.
        levelspace = Vector{Any}(undef, Wsize)
        levelspace[1] = U1Space(0 => 1); levelspace[Wsize] = U1Space(0 => 1)
        for k in 1:lattice.F
            levelspace[1 + 2*(k-1) + 1] = U1Space( lattice.q => 1)
            levelspace[1 + 2*(k-1) + 2] = U1Space(-lattice.q => 1)
        end
        for d in sort(collect(defects); by = x -> x.link, rev = true)   # right-to-left
            (2 ≤ d.link ≤ Int(lattice.N)) ||
                throw(ArgumentError("defect link $(d.link) must lie in 2:$(Int(lattice.N))"))
            idx = (d.link - 1) * lattice.F + 1   # between physical sites link-1 and link
            Pdef = U1Space(d.charge => 1)
            Wdef = Matrix{eltype(A[1])}(missing, Wsize, Wsize)
            # corners (identity-in / accumulator channels): scalar 1.0 → braiding;
            # middle channels (hop-left/right): explicit identity on the 1-D leg.
            Wdef[1, 1] = 1.0
            Wdef[Wsize, Wsize] = 1.0
            for i in 2:Wsize-1
                Wdef[i, i] = ones(levelspace[i] ⊗ Pdef ← Pdef ⊗ levelspace[i])
            end
            insert!(A, idx, Wdef)
            insert!(link_fcts, idx, missing)
        end
    end

    if isinf(lattice.N)
        mpo = InfiniteMPOHamiltonian(A)
    else
        _trim_finite_mpo!(A, lattice.F)
        mpo = FiniteMPOHamiltonian(A)
    end

    lempo = if isinf(lattice.N)
        InfiniteLEMPOHamiltonian(mpo, link_fcts)
    else
        FiniteLEMPOHamiltonian(mpo, link_fcts)
    end

    return MPSKitOperator(lattice, lempo, universe, defects)
end