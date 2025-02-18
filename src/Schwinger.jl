module Schwinger

export SchwingerLattice
export findstates!
export energies, occupations, charges, electricfields, entanglements
export setmass!, setspacing!, setθ2π!

using ITensors, ITensorMPS
using Parameters

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

"""
`SchwingerLattice(N, Nf; periodic=false, L_max=0)`

Constructs a SchwingerLattice for the Schwinger model.

# Arguments
- `N::Int`: Number of sites.
- `Nf::Int`: Number of flavors.
- `periodic::Bool=false`: Whether the lattice is periodic.
- `L_max::Int=0`: Maximum absolute value for L_0 when periodic.

# Returns
A `SchwingerLattice` object.

"""
@with_kw mutable struct SchwingerLattice
    sites::Vector{Index}
    N::Int
    Nf::Int
    periodic::Bool
    L_max::Int

    q::Int
    a::Real
    m::Vector{Real}
    mlat::Vector{Real}
    θ2π::Real
    hamiltonian_updated::Bool

    psis::Vector{ITensorMPS.MPS}
    eigenstates::Bool

    H::ITensorMPS.MPO

    function SchwingerLattice(N::Int, Nf::Int; periodic::Bool=false, L_max::Int=0, θ2π::Real=0, m::Union{Real,Vector{Float64}}=0, a::Real=1, q::Int=1)
        if periodic && L_max == 0
            L_max = 3  # Default value for L_max when periodic is true
        elseif !periodic && L_max > 0
            periodic = true  # Setting L_max overrides setting periodic
        end

        total_sites = periodic ? N * Nf + 1 : N * Nf
        sites = Vector{Index}(undef, total_sites)

        # Create each site in the SchwingerLattice
        for i in 1:total_sites
            if periodic && i == total_sites
                sites[i] = Index([QN() => 2*L_max + 1]; tags="Boson,Site,n=$i")
            else
                sites[i] = siteind("S=1/2", i; conserve_qns = true)
            end
        end

        if typeof(m) <: Real
            m = [m for _=1:Nf]
        end
        new(sites, N, Nf, periodic, L_max, q, a, m, m .- q^2*Nf*a/8, θ2π, false, [], true)
    end
end



"""
`update_hamiltonian()`

Computes the Hamiltonian MPO using the current lattice parameters.
"""
function update_hamiltonian(self::SchwingerLattice)
    @unpack m, q, a, θ2π, N, Nf, periodic, sites, hamiltonian_updated = self

    if hamiltonian_updated
        return
    end

    update_mlat(self)
    @unpack mlat = self
    self.psis = []

    hamiltonian = OpSum()

    # Hopping term
    
    for j in 1:N-1
        for k in 1:Nf
            ind = Nf*(j - 1) + k
            hamiltonian += 1/(2*a),"S+",ind,"S-",ind + Nf
            hamiltonian += 1/(2*a),"S-",ind,"S+",ind + Nf
        end
    end

    if periodic
        for k in 1:Nf
            ind = Nf * (N - 1) + k
            hamiltonian += 1/(2*a),"S+",ind,"S-",k,"raise",N * Nf + 1
            hamiltonian += 1/(2*a),"S-",ind,"S+",k,"lower",N * Nf + 1
        end
    end

    # Mass term

    for j in 1:N
        for k in 1:Nf
            ind = Nf*(j - 1) + k
            hamiltonian += mlat[k] * ((-1) ^ (j + 1)),"Sz",ind
        end
    end

    # Gauge kinetic term

    if periodic
        hamiltonian += q^2 * a * N / 2,"L0",N * Nf + 1,"L0",N * Nf + 1
        hamiltonian += q * a * N * θ2π,"L0",N * Nf + 1

        hamiltonian += q^2 * a * N * Nf / 4,"L0",N * Nf + 1
        
        for j in 1:N
            for k in 1:Nf
                ind = Nf*(j - 1) + k
                hamiltonian += q^2 * a * (N - j),"L0",N * Nf + 1,"Sz",ind
            end
        end
    end
    hamiltonian += a * N / 2 * (θ2π ^ 2),"Id",1
    hamiltonian += q * a * N * Nf / 4 * θ2π,"Id",1
    for j in 1:N
        for k in 1:Nf
            ind = Nf*(j - 1) + k
            hamiltonian += q * a * (N - j) * θ2π,"Sz",ind
        end
    end

    hamiltonian += q^2 * a * N * Nf * Nf / 16,"Id",1
    for j in 1:N
        for k in 1:Nf
            ind = Nf*(j - 1) + k
            hamiltonian += q^2 * a * Nf * ((N + 1 - j) ÷ 2) / 2,"Sz",ind
        end
    end

    for j1 in 1:N
        for k1 in 1:Nf
            ind1 = Nf * (j1 - 1) + k1
            for j2 in j1:N
                for k2 in (j2 == j1 ? (k1:Nf) : 1:Nf)
                    ind2 = Nf*(j2 - 1) + k2
                    hamiltonian += q^2*a*(N - j2)/(j1 == j2 && k1 == k2 ? 2 : 1),"Sz",ind1,"Sz",ind2
                end
            end
        end
    end

    self.H = MPO(hamiltonian, sites)

    self.hamiltonian_updated = true
end

"""
`findstates(lattice, nstates)`

Use DMRG to find the lowest few eigenstates.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `nstates::Int`:: number of states to determine.
"""

function findstates!(self::SchwingerLattice, nstates::Int; maxiters::Int = 500, maxlinkdim::Int = 600, energy_tol::Real = 1E-6, weight::Real = 100, outputlevel::Int = 1, minsweeps::Int = 10)
    update_hamiltonian(self)
    @unpack N, Nf, L_max, periodic, sites, hamiltonian_updated = self

    for idx in (length(self.psis) + 1):nstates
        state = [n == N * Nf + 1 ? L_max + 1 : isodd(floor((n-1)/Nf)) ? "Up" : "Dn" for n=1:(N * Nf + (periodic ? 1 : 0))]
        psi = randomMPS(sites, state; linkdims = 4)

        sweeps = Sweeps(maxiters)
        maxdim!(sweeps, maxlinkdim)

        if outputlevel >= 1
            println("Running DMRG for state $idx")
        end
        _, psi = ITensorMPS.dmrg(self.H, self.psis, psi, sweeps; weight = weight, observer = DMRGObserver(;energy_tol = energy_tol,minsweeps=minsweeps), outputlevel = outputlevel)
        if outputlevel >= 1
            println()
        end

        push!(self.psis, psi)
    end

    self.eigenstates = true
end

"""
`energies(lattice)`

Return the Hamiltonian eigenvalues of all known eigenstates.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function energies(self::SchwingerLattice)
    @unpack H, psis = self

    return [inner(psi', H, psi) for psi=psis]
end

"""
`occupations(lattice)`

Return a list of the expectations of χ†χ operators on each site and for each known eigenstate.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function occupations(self::SchwingerLattice)
    @unpack psis, N, Nf = self

    return [expect(psi, "Sz", sites=1:N*Nf) + [1/2 for n=1:N*Nf] for psi=psis]
end

"""
`charges(lattice)`

Return a list of the expectations of Q operators on each site and for each known eigenstate.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function charges(self::SchwingerLattice)
    @unpack N, Nf = self
    occs = occupations(self)

    return [occ + [-1/2 + (-1)^((n - 1) ÷ Nf)/2 for n=1:N*Nf] for occ=occs]
end

"""
`electricfields(lattice)`

Return a list of the expectations of (L + θ/2π) operators on each site and for each known eigenstate.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function electricfields(self::SchwingerLattice)
    @unpack N, Nf, θ2π, psis = self
    Qs = charges(self)

    L0s = [self.periodic ? expect(psi, "L0", sites=N*Nf + 1) : 0 for psi=psis]

    return [accumulate(+, Qs[i]) .+ (L0s[i] + θ2π) for i=eachindex(psis)]
end

"""
`entanglements(lattice)`

Return a list of the entanglement entropies for each bisection of the lattice and for each known eigenstate.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
function entanglements(self::SchwingerLattice)
    @unpack psis, N, Nf = self

    data = zeros(length(psis), N*Nf)
    for i=eachindex(psis)
        psi = psis[i]
        for j=1:N*Nf-1
            orthogonalize!(psi, j)
            _,S,_ = svd(psi[j], (linkinds(psi, j-1)..., siteinds(psi,j)...))
            for n=1:dim(S, 1)
                p = S[n,n]^2
                data[i,j] -= p * log(p)
            end
        end
    end

    data = [data[i,:] for i in axes(data)[1]]

    return data
end

"""
`holonomy(lattice, conjugate = false)`

Returns the holonomy operator acting on `lattice`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""

function holonomy(self::SchwingerLattice, conjugate::Bool = false)
    @unpack sites, N, Nf, periodic = self

    if !periodic
        throw(DomainError(self, "Lattice must be periodic."))
    end

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * Nf + 1

    return MPO(holonomy, sites)
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

function wilsonline(self::SchwingerLattice, i::Integer, j::Integer, flavor::Integer = 1, conjugate::Bool = false)
    @unpack sites, N, Nf, periodic = self

    if i < 1 || i > N || j < 1 || j > N || j > i
        throw(DomainError(self, "Site indices must satisfy 1 ≤ i < j ≤ N"))
    end

    line = OpSum()
    line += 1,(conjugate ? "S-" : "S+"),(i - 1)*Nf + flavor,(conjugate ? "S+" : "S-"),(j - 1)*Nf + flavor

    return MPO(line, sites)
end

"""
`setspacing!(lattice, a)`

Sets the lattice spacing in coupling units to `a`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `a::Real`:: lattice spacing.
"""
function setspacing!(self::SchwingerLattice, a::Real)
    if a != self.a
        self.hamiltonian_updated = false
    end
    self.a = a
    update_mlat(self)
end

"""
`setmass!(lattice, m)`

Sets the physical mass in coupling units to `m`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `m::Real`:: Physical mass.
"""
function setmass!(self::SchwingerLattice, m::Real)
    if [m for k=1:self.Nf] != self.m
        self.hamiltonian_updated = false
    end
    self.m = [m for k=1:self.Nf]
    update_mlat(self)
end

"""
`setmass!(lattice, m)`

Sets the physical masses (in coupling units) for each flavor to be the elements of `m`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `m::Vector{Real}`:: Physical masses for each flavor.
"""
function setmass!(self::SchwingerLattice, m::Vector{Real})
    if m != self.m
        self.hamiltonian_updated = false
    end
    self.m = m
    update_mlat(self)
end

"""
`setθ2π!(lattice, θ2π)`

Sets the background electric field to `θ2π`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `θ2π::Real`:: Theta parameter.
"""
function setθ2π!(self::SchwingerLattice, θ2π::Real)
    if θ2π != self.θ2π
        self.hamiltonian_updated = false
    end
    self.θ2π = θ2π
end

"""
`setcharge!(lattice, θ2π)`

Sets the background electric field to `θ2π`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `θ2π::Real`:: Theta parameter.
"""
function setcharge!(self::SchwingerLattice, q::Int)
    if q != self.q
        self.hamiltonian_updated = false
    end
    self.q = q
    update_mlat(self)
end

function update_mlat(self::SchwingerLattice)
    self.mlat = self.m .- self.q^2*self.Nf*self.a/8
end

end
