"""
`wilsonloop(hamiltonian, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function MPOWilsonLoop(lattice::SchwingerLattice{N,F}, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

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

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * F + 1
    
    mpo = MPO(holonomy, sites(lattice; L_max=L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

"""
`wilsonloop(hamiltonian, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function EDWilsonLoop(lattice::SchwingerLattice{N,F}, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

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

    function shiftL₀(state::SchwingerBasisState{N,F})
        occs = occupations(state)
        return [(occs, L₀(state) + (conjugate ? -1 : 1)) => 1]
    end
    return constructoperator(lattice, shiftL₀; L_max = L_max, universe = universe)
end

"""
`MPOWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)`

Returns the spatial Wilson line operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson line.
- `flavor::Int`: Flavor of the Wilson line.
- `start::Int`: Starting site of the Wilson line.
- `finish::Int`: Finishing site of the Wilson line.
"""
function MPOWilsonLine(lattice::SchwingerLattice{N,F}, conjugate::Bool = false, flavor::Int = 1, start::Union{Nothing,Int} = nothing, finish::Union{Nothing,Int} = nothing; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
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

    if isnothing(start)
        start = 1
    end
    if isnothing(finish)
        finish = N
    end

    line = OpSum()

    ind1 = F*(start-1) + flavor
    ind2 = F*(finish-1) + flavor
    if conjugate
        line += 1,"S-",ind1,"S+",ind2
    else
        line += 1,"S+",ind1,"S-",ind2
    end

    mpo = MPO(line, sites(lattice; L_max=L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
end

"""
`EDWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)`

Returns the spatial Wilson line operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson line.
- `flavor::Int`: Flavor of the Wilson line.
- `start::Int`: Starting site of the Wilson line.
- `finish::Int`: Finishing site of the Wilson line.
"""
function EDWilsonLine(lattice::SchwingerLattice{N,F}, conjugate::Bool = false, flavor::Int = 1, start::Union{Nothing,Int} = nothing, finish::Union{Nothing,Int} = nothing; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
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
    
    if isnothing(start)
        start = 1
    end
    if isnothing(finish)
        finish = N
    end

    function actW(state::SchwingerBasisState{N,F})
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