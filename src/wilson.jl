"""
`ITensorWilsonLoop(lattice, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function ITensorWilsonLoop(lattice::Lattice, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) 
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * F + 1
    
    mpo = ITensorMPS.MPO(holonomy, get_sites(lattice; L_max=L_max))
    return ITensorOperator(lattice, mpo, L_max, universe)
end

"""
`wilsonloop(hamiltonian, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function EDWilsonLoop(lattice::Lattice, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) 
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    function shiftL₀(state::BasisState)
        occs = occupations(state)
        return [(occs, L₀(state) + (conjugate ? -1 : 1)) => 1]
    end
    return constructoperator(lattice, shiftL₀; L_max = L_max, universe = universe)
end

"""
`ITensorWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)`

Returns the spatial Wilson line operator for `lattice`.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson line.
- `flavor::Int`: Flavor of the Wilson line.
- `start::Int`: Starting site of the Wilson line.
- `finish::Int`: Finishing site of the Wilson line.
"""
function ITensorWilsonLine(lattice::Lattice, conjugate::Bool = false, flavor::Int = 1, start::Int = 1, finish::Int = Int(lattice.N); L_max::Union{Nothing,Int} = nothing, universe::Int = 0) 
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    N, F = Int(lattice.N), lattice.F
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    line = OpSum()

    ind1 = F*(start-1) + flavor
    ind2 = F*(finish-1) + flavor
    if start == finish
        line += 1,"S+",ind1,"S-",ind1          # zero-length line = number operator χ†ₙχₙ
    elseif conjugate
        line += 1,"S-",ind1,"S+",ind2
    else
        line += 1,"S+",ind1,"S-",ind2
    end

    mpo = ITensorMPS.MPO(line, get_sites(lattice; L_max=L_max))
    return ITensorOperator(lattice, mpo, L_max, universe)
end

"""
`EDWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)`

Returns the spatial Wilson line operator for `lattice`.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson line.
- `flavor::Int`: Flavor of the Wilson line.
- `start::Int`: Starting site of the Wilson line.
- `finish::Int`: Finishing site of the Wilson line.
"""
function EDWilsonLine(lattice::Lattice, conjugate::Bool = false, flavor::Int = 1, start::Int = 1, finish::Int = Int(lattice.N); L_max::Union{Nothing,Int} = nothing, universe::Int = 0) 
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += lattice.q
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    if start == finish
        function numop(state::BasisState)      # zero-length line = number operator χ†ₙχₙ
            occs = occupations(state)
            return occs[start, flavor] == 1 ? [(copy(occs), L₀(state)) => 1] : []
        end
        return constructoperator(lattice, numop; L_max = L_max, universe = universe)
    end

    function actW(state::BasisState)
        occs = copy(occupations(state))
        if occs[start, flavor] == (conjugate ? 1 : 0) && occs[finish, flavor] == (conjugate ? 0 : 1)
            occs[start, flavor] = (conjugate ? 0 : 1)
            occs[finish, flavor] = (conjugate ? 1 : 0)
            return [(occs, L₀(state)) => 1]
        end
        return []
    end
    return constructoperator(lattice, actW; L_max = L_max, universe = universe)
end

# =============================================================================
# MPSKit Backend
# =============================================================================

"""
`MPSKitWilsonLoop(lattice, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice` using MPSKit.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function MPSKitWilsonLoop(lattice::Lattice, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) 
    error("Periodic lattices are not currently supported in MPSKit backend.")
end

"""
`MPSKitWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)`

Returns the spatial Wilson line operator for `lattice` using MPSKit.

# Arguments
- `lattice::Lattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson line.
- `flavor::Int`: Flavor of the Wilson line.
- `start::Int`: Starting site of the Wilson line.
- `finish::Int`: Finishing site of the Wilson line.
"""
function MPSKitWilsonLine(lattice::Lattice, conjugate::Bool = false, flavor::Int = 1, start::Int = 1, finish::Int = Int(lattice.N); universe::Int = 0) 
    if lattice.periodic
        throw(ArgumentError("MPSKit backend not implemented for periodic lattices"))
    end
    _, universe = process_L_max_universe(lattice, 0, universe)
    N, F = lattice.N, lattice.F
    if isinf(N)
        throw(ArgumentError("Wilson lines not implemented for infinite lattices"))
    end
    N = Int(N)

    start_ind = F*(start-1) + flavor
    finish_ind = F*(finish-1) + flavor

    spaces = get_mpskit_spaces(lattice)

    if start_ind == finish_ind
        # zero-length line = number operator χ†ₙχₙ: an on-site, charge-conserving MPO.
        numop(P) = begin
            T = zeros(ComplexF64, U1Space(0 => 1) ⊗ P ← P ⊗ U1Space(0 => 1))
            for (i, c) in enumerate(sort(collect(sectors(P)); by = c -> c.charge))
                block(T, c) .= (i == 2 ? 1.0 : 0.0)   # occupied → 1, empty → 0
            end
            T
        end
        idop(P) = isomorphism(ComplexF64, U1Space(0 => 1) ⊗ P, P ⊗ U1Space(0 => 1))
        MPOT = typeof(idop(spaces[1]))
        ts = MPOT[k == start_ind ? numop(spaces[k]) : idop(spaces[k]) for k in 1:N*F]
        return MPSKitOperator(lattice, FiniteMPO(ts), universe)
    end

    start_space = spaces[start_ind]
    finish_space = spaces[finish_ind]

    charge = conjugate ? -lattice.q : lattice.q

    startop = ones(U1Space(0 => 1) ⊗ start_space ← start_space ⊗ U1Space(charge => 1))
    finishop = ones(U1Space(charge => 1) ⊗ finish_space ← finish_space ⊗ U1Space(0 => 1))

    Elt = Union{Missing, typeof(startop), scalartype(startop)}
    A = Vector{Matrix{Elt}}(undef, N*F)

    for n in 1:N*F
        W = Matrix{Elt}(missing, 3, 3)
        W[1, 1] = 1.0
        W[end, end] = 1.0

        if n == start_ind
            W[1, 2] = startop
        elseif n == finish_ind
            W[2, end] = finishop
        elseif start_ind < n < finish_ind
            W[2, 2] = convert(TensorMap, BraidingTensor(spaces[n], U1Space(charge => 1)))
        else
            W[2, 2] = convert(TensorMap, BraidingTensor(spaces[n], U1Space(0 => 1)))
        end

        A[n] = W
    end
    A[1] = A[1][1:1, :]
    if start != 1
        A[1][1,2] = convert(TensorMap, BraidingTensor(spaces[1], U1Space(0 => 1)))
    end
    A[end] = A[end][:, end:end]
    if finish != N
        A[end][2,1] = convert(TensorMap, BraidingTensor(spaces[end], U1Space(0 => 1)))
    end

    mpo = FiniteMPOHamiltonian(A)
    return MPSKitOperator(lattice, mpo, universe)
end
# =============================================================================
# Single staggered-fermion field  φ̃_s = ∏_{k<s}(±iσ^z_k) σ∓_s   (MPSKit)
# =============================================================================

"""
`MPSKitFermionField(lattice, site; dagger=false, flavor=1, universe=0)`

The single staggered-fermion field `φ̃_site = ∏_{k<site}(-iσ^z_k) σ⁻_site`
(annihilation; `dagger=true` gives the creation field with `σ⁺` and `+iσ^z`).

Unlike `MPSKitWilsonLine` (a charge-conserving bilinear), this is a *single* fermion
field: it changes the total U(1) charge by `∓q`, so acting with it sends the MPS into
a different symmetry sector. Returned as an `MPSKitOperator`; apply with `act`/`*`.
"""
function MPSKitFermionField(lattice::Lattice, site::Int; dagger::Bool = false, flavor::Int = 1, universe::Int = 0)
    if lattice.periodic
        throw(ArgumentError("MPSKitFermionField not implemented for periodic lattices"))
    end
    if isinf(lattice.N)
        throw(ArgumentError("MPSKitFermionField not implemented for infinite lattices"))
    end
    _, universe = process_L_max_universe(lattice, 0, universe)
    N, F = Int(lattice.N), lattice.F
    q = lattice.q
    spaces = get_mpskit_spaces(lattice)
    ind = F * (site - 1) + flavor          # MPSKit index carrying the fermion operator
    Δq    = dagger ? q : -q                # charge injected into the bond (to the right)
    phase = dagger ? 1im : -1im            # JW string factor (±i)·σ^z to the left

    siteop(P) = ones(ComplexF64, U1Space(0 => 1) ⊗ P ← P ⊗ U1Space(Δq => 1))   # σ∓
    function jwstr(P)                       # (±i)·σ^z  (diagonal: empty −1, occupied +1)
        T = zeros(ComplexF64, U1Space(0 => 1) ⊗ P ← P ⊗ U1Space(0 => 1))
        for (i, c) in enumerate(sort(collect(sectors(P)); by = c -> c.charge))
            block(T, c) .= phase * (i == 2 ? 1.0 : -1.0)
        end
        return T
    end
    passthru(P) = isomorphism(ComplexF64, U1Space(Δq => 1) ⊗ P, P ⊗ U1Space(Δq => 1))

    MPOT = typeof(siteop(spaces[ind]))
    ts = MPOT[k < ind ? jwstr(spaces[k]) : k == ind ? siteop(spaces[k]) : passthru(spaces[k]) for k in 1:N*F]
    return MPSKitOperator(lattice, FiniteMPO(ts), universe)
end

"""
`EDFermionField(lattice, site; dagger=false, flavor=1)`

ED version of the staggered fermion field `φ̃_site` (annihilation; `dagger=true` for
creation), as a sparse operator with the Jordan–Wigner string. It maps the
charge-`c` sector to the charge-`c ∓ q` sector (annihilation/creation), so acting
with it on the (charge-neutral) vacuum lands in the charge-`∓q` basis.

!!! note
    The ED Hamiltonian uses the imaginary staggered hopping `−i(χ†ₙχₙ₊₁ − h.c.)`, while
    the Jordan–Wigner string here is in the real-hopping gauge of the ITensors/MPSKit
    backends. The two differ by `χₙ → (∓i)ⁿ χₙ`, so a staggered factor `(∓i)^site` keeps
    the field consistent with the ED ground state.
"""
function EDFermionField(lattice::Lattice, site::Int; dagger::Bool = false, flavor::Int = 1,
                        L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    F = lattice.F
    ind = F * (site - 1) + flavor
    phase = dagger ? 1im : -1im
    stagger = phase^site                       # gauge factor for ED's imaginary hopping
    function action(state::BasisState)
        occs = occupations(state)
        occs[site, flavor] == (dagger ? 0 : 1) || return Pair{Tuple{BitMatrix,Int},ComplexF64}[]
        new = copy(occs)
        new[site, flavor] = dagger ? 1 : 0
        sgn = ComplexF64(phase^(ind - 1)) * stagger
        for k in 1:ind-1                       # Jordan–Wigner σ^z string (σ^z = 2n−1)
            sgn *= (2 * occs[(k - 1) ÷ F + 1, (k - 1) % F + 1] - 1)
        end
        return [(BitMatrix(new), L₀(state)) => sgn]
    end
    # σ⁻ removes a charge-q fermion (out_charge = −q); the creation field adds one (+q).
    out_charge = dagger ? lattice.q : -lattice.q
    return constructoperator(lattice, action; L_max = L_max, universe = universe, in_charge = 0, out_charge = out_charge)
end

"""
`ITensorFermionField(lattice, site; dagger=false, flavor=1)`

ITensors version of the staggered fermion field `φ̃_site` (annihilation; `dagger=true`
for creation): the single-site `S∓` operator dressed with the Jordan–Wigner `σ^z`
string, as an MPO. Acting with it changes the total quantum number of the MPS.
"""
function ITensorFermionField(lattice::Lattice, site::Int; dagger::Bool = false, flavor::Int = 1,
                             L_max::Union{Nothing,Int} = nothing, universe::Int = 0,
                             sites::Union{AbstractVector{<:Index},Nothing} = nothing)
    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    F = lattice.F
    ind = F * (site - 1) + flavor
    phase = dagger ? 1im : -1im
    term = Any[(2 * phase)^(ind - 1)]          # σ^z = 2·Sz  ⇒  (±2i)^(ind−1)
    for k in 1:ind-1
        push!(term, "Sz", k)
    end
    push!(term, dagger ? "S+" : "S-", ind)
    os = OpSum()
    os += Tuple(term)
    mpo = ITensorMPS.MPO(os, isnothing(sites) ? get_sites(lattice; L_max = L_max) : sites)
    return ITensorOperator(lattice, mpo, L_max, universe)
end
