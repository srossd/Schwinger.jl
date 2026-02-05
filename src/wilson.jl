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
    if conjugate
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