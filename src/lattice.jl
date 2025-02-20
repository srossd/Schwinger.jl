"""
`SchwingerLattice(N, F; periodic=false, L_max=0)`

Constructs a SchwingerLattice for the Schwinger model.

# Arguments
- `N::Int`: Number of sites.
- `F::Int`: Number of flavors.
- `periodic::Bool=false`: Whether the lattice is periodic.
- `L_max::Int=0`: Maximum absolute value for L_0 when periodic.

# Returns
A `SchwingerLattice` object.

"""
struct SchwingerLattice{N,F}
    # global parameters
    q::Int
    periodic::Bool
    a::Float64

    # local parameters (typically position-independent)
    θ2π::NTuple{N,Float64}
    mlat::NTuple{N,NTuple{F,Float64}}
    mprime::NTuple{N,NTuple{F,Float64}}

    function SchwingerLattice{N,F}(;
        q::Int = 1, periodic::Bool=false, 
        L::Union{Nothing,Real}=nothing, a::Union{Nothing,Real}=nothing, 
        m::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}}=0., 
        mlat::Union{Nothing,Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}}=nothing, 
        mprime::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}}=0., 
        θ2π::Union{Real,NTuple{N,Real}}=0.) where {N,F}

        if !iseven(N) || N < 2
            throw(DomainError(N, "Number of sites must be even and ≥ 2."))
        end
        if F < 1
            throw(DomainError(F, "Number of flavors must be ≥ 1."))
        end

        if isnothing(a) && isnothing(L)
            a = 1.
        end

        if !isnothing(L)
            if !isnothing(a)
                throw(ArgumentError("Cannot specify both lattice spacing and length."))
            end
            a = L/N
        end

        if isa(θ2π, Float64)
            θ2π = NTuple{N, Float64}((θ2π for _ in 1:N))
        end

        if isnothing(mlat)
            if isa(m, Float64)
                mlat = NTuple{N, NTuple{F, Float64}}((NTuple{F,Float64}((m - F*q^2*a/8 for _ in 1:F)) for _ in 1:N))
            elseif isa(m, NTuple{F, Float64})
                mlat = NTuple{N, NTuple{F, Float64}}((m .- F*q^2*a/8 for _ in 1:N))
            end
        end

        if isa(mprime, Float64)
            mprime = NTuple{N, NTuple{F, Float64}}((NTuple{F,Float64}((mprime for _ in 1:F)) for _ in 1:N))
        elseif isa(mprime, NTuple{F, Float64})
            mprime = NTuple{N, NTuple{F, Float64}}((mprime for _ in 1:N))
        end
        
        new(q, periodic, a, θ2π, mlat, mprime)
    end
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
`setlength!(lattice, L)`

Sets the lattice length in coupling units to `L`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `L::Real`:: lattice length.
"""
function setlength!(self::SchwingerLattice, L::Real)
    setspacing!(self, L/self.N)
end

"""
`setmass!(lattice, m)`

Sets the physical mass in coupling units to `m`.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `m::Real`:: Physical mass.
"""
function setmass!(self::SchwingerLattice, m::Real)
    if [m for k=1:self.F] != self.m
        self.hamiltonian_updated = false
    end
    self.m = [m for k=1:self.F]
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
    self.mlat = self.m .- self.q^2*self.F*self.a/8
    return
end