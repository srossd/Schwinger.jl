function ITensors.op(::OpName"L0", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d
      mat[k, k] = (k - 1) - (d ÷ 2)
    end
    return mat
end
function ITensors.op(::OpName"raise", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d-1
      mat[k + 1, k] = 1
    end
    return mat
end
function ITensors.op(::OpName"lower", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d-1
      mat[k, k + 1] = 1
    end
    return mat
end

tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = tuplejoin(x, tuplejoin(y, z))

function wigner_string(N::Int, F::Int, flavor::Int)
    return tuplejoin((("Sz",j*N + k) for j in 1:N)...)
end

@memoize function sites(lattice::SchwingerLattice{N,F}; L_max::Int = 3) where {N,F}
    total_sites = lattice.periodic ? N * F + 1 : N * F
    sites = Vector{Index}(undef, total_sites)

    # Create each site in the SchwingerLattice
    for i in 1:total_sites
        if lattice.periodic && i == total_sites
            sites[i] = Index([QN() => 2*L_max + 1]; tags="Boson,Site,n=$i")
        else
            sites[i] = siteind("S=1/2", i; conserve_qns = true)
        end
    end

    return sites
end

"""
`holonomy(lattice, conjugate = false)`

Returns the holonomy operator acting on `lattice`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function holonomy(lattice::SchwingerLattice{N,F}, conjugate::Bool = false; L_max::Int = 3) where {N,F}
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