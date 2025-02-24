"""
`holonomy(hamiltonian, conjugate = false)`

Returns the holonomy operator corresponding to `hamiltonian`.

# Arguments
- `hamiltonian::EDHamiltonian`: Schwinger model Hamiltonian.
"""
function holonomy(hamiltonian::EDHamiltonian{N,F}, conjugate::Bool = false) where {N,F}
    if !hamiltonian.lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * F + 1

    return MPO(holonomy, sites(hamiltonian))
end

"""
`holonomy(hamiltonian, conjugate = false)`

Returns the holonomy operator corresponding to `hamiltonian`.

# Arguments
- `hamiltonian::EDHamiltonian`: Schwinger model Hamiltonian.
"""
function holonomy(lattice::SchwingerLattice{N,F}, conjugate::Bool = false) where {N,F}
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * F + 1

    return MPO(holonomy, sites(lattice; L_max = L_max))
end

"""
`wilsonline(lattice, i, j, flavor = 1, conjugate = false)`

Returns the Wilson line operator acting on `lattice` between sites `i` and `j`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `i::Integer`: Starting site
- `j::Integer`: Ending site
- `flavor::Integer`: Fermion flavor at endpoints
"""
function wilsonline(lattice::SchwingerLattice{N,F}, i::Integer, j::Integer; flavor::Integer = 1, L_max::Int = 3, conjugate::Bool = false) where {N,F}
    if i < 1 || i > N || j < 1 || j > N || j > i
        throw(DomainError(self, "Site indices must satisfy 1 ≤ i < j ≤ N"))
    end

    line = OpSum()
    line += 1,(conjugate ? "S-" : "S+"),(i - 1)*F + flavor,(conjugate ? "S+" : "S-"),(j - 1)*F + flavor

    return MPO(line, sites(lattice; L_max = L_max))
end