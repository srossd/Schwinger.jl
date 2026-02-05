# =============================================================================
# LatticeSize: Represents finite or infinite lattice sizes
# =============================================================================

import Base: isfinite, isinf

"""
    LatticeSize

Represents the size of a lattice, which can be either finite or infinite.

# Constructors
- `LatticeSize(n::Integer)`: Create a finite lattice of size `n`
- `LatticeSize(::Type{Inf})` or `LatticeSize(Inf)`: Create an infinite lattice

# Examples
```julia
finite = LatticeSize(10)
infinite = LatticeSize(Inf)
```
"""
struct LatticeSize
    value::Int
    isinfinite::Bool
    
    LatticeSize(n::Integer) = new(Int(n), false)
    LatticeSize(::Type{T}) where T<:Real = new(0, true)
end

# Constructor for Inf literal
LatticeSize(::typeof(Inf)) = LatticeSize(Float64)

# Comparison operations
Base.:(==)(a::LatticeSize, b::LatticeSize) = (a.isinfinite == b.isinfinite) && (a.isinfinite || a.value == b.value)
Base.:(==)(a::LatticeSize, b::Integer) = !a.isinfinite && a.value == b
Base.:(==)(a::Integer, b::LatticeSize) = b == a
Base.:(==)(a::LatticeSize, ::typeof(Inf)) = a.isinfinite
Base.:(==)(::typeof(Inf), a::LatticeSize) = a.isinfinite

# Property checks
Base.isfinite(n::LatticeSize) = !n.isinfinite
Base.isinf(n::LatticeSize) = n.isinfinite

# Arithmetic operations (for finite sizes)
Base.:(*)(a::LatticeSize, b::Real) = a.isinfinite ? Inf : a.value * b
Base.:(*)(a::Real, b::LatticeSize) = b * a
Base.:(/)(a::LatticeSize, b::Real) = a.isinfinite ? Inf : a.value / b

# Convert to Int (throws error if infinite)
Base.Int(n::LatticeSize) = n.isinfinite ? error("Cannot convert infinite LatticeSize to Int") : n.value

# Show method
Base.show(io::IO, n::LatticeSize) = print(io, n.isinfinite ? "∞" : string(n.value))

# =============================================================================
# Lattice Structure
# =============================================================================

"""
`Lattice(;kwargs...)`

Constructs a Lattice for the Schwinger model.

# Arguments
- `N::Union{Integer,Inf}`: Number of sites (use `Inf` for infinite lattices).
- `F::Int=1`: Number of flavors.
- `periodic::Bool=false`: Whether the lattice is periodic.
- `q::Int=1`: Charge.
- `L::Union{Nothing,Real}=nothing`: Length of the lattice.
- `a::Union{Nothing,Real}=nothing`: Lattice spacing.
- `m::Union{Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=0.`: Mass parameter.
- `mlat::Union{Nothing,Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=nothing`: Lattice mass parameter.
- `mprime::Union{Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=0.`: Hopping mass parameter.
- `θ2π::Union{Real,NTuple{N,Real}}=0.`: Theta angle.

# Returns
A `Lattice` object.

"""
struct Lattice
    # global parameters
    N::LatticeSize
    F::Int
    q::Int
    periodic::Bool
    a::Float64
    L::Float64

    # local parameters (typically position-independent)
    θ2π::AbstractVector{Float64}
    m::AbstractVector{AbstractVector{Float64}}
    mlat::AbstractVector{AbstractVector{Float64}}
    mprime::AbstractVector{AbstractVector{Float64}}

    function Lattice(N::Union{Real,Type{<:Real}}; F::Int = 1,
        q::Int = 1, periodic::Bool=false, 
        L::Union{Nothing,Real}=nothing, a::Union{Nothing,Real}=nothing, 
        m::Union{Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=0., 
        mlat::Union{Nothing,Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=nothing, 
        mprime::Union{Real,AbstractVector{<:Real},AbstractArray{<:Real,2}}=0., 
        θ2π::Union{Real,AbstractVector{<:Real}}=0.)

        # Convert N to LatticeSize
        N = LatticeSize(N)

        if isfinite(N) && (!iseven(Int(N)) || Int(N) < 2)
            throw(DomainError(Int(N), "Number of sites must be even and ≥ 2."))
        end
        if isinf(N) && periodic
            throw(ArgumentError("Infinite lattices cannot be periodic."))
        end
        if F < 1
            throw(DomainError(F, "Number of flavors must be ≥ 1."))
        end
        if q < 1
            throw(DomainError(q, "Charge must be ≥ 1."))
        end

        if isnothing(a) && isnothing(L)
            a = 1.
        end

        if !isnothing(L) && L != Inf
            if isinf(N)
                throw(ArgumentError("Finite lattice length L cannot be specified for infinite lattices."))
            end
            if !isnothing(a) && !isapprox(L, a * Int(N); atol=1e-10, rtol=1e-10)
                throw(ArgumentError("Lattice spacing $L not equal to $a * $N."))
            end
            a = L/Int(N)
        else
            L = isinf(N) ? Inf : Int(N)*a
        end

        if θ2π isa Real
            θ2π = if isinf(N)
                PeriodicArrays.PeriodicVector([θ2π])
            else
                fill(θ2π, Int(N))
            end
        end

        # Helper function to convert mass parameters to proper format
        function process_mass_param(param, param_name)
            if isa(param, Real)
                # Scalar: broadcast to all sites and flavors
                if isinf(N)
                    return PeriodicArrays.PeriodicVector([fill(param, F)])
                else
                    return [fill(param, F) for _=1:Int(N)]
                end
            elseif isa(param, AbstractVector{<:Real})
                if length(param) == F
                    # Length-F vector: same for all sites
                    if isinf(N)
                        return PeriodicArrays.PeriodicVector([collect(param)])
                    else
                        return [collect(param) for _=1:Int(N)]
                    end
                elseif isfinite(N) && length(param) == Int(N)
                    # Length-N vector: broadcast across flavors at each site
                    return [fill(param[i], F) for i=1:Int(N)]
                else
                    throw(DomainError(param, "If $param_name is a vector, must have length F=$F or N=$N."))
                end
            else
                # Array: must be N×F
                if isinf(N)
                    throw(DomainError(param, "Cannot specify matrix $param_name for infinite lattices."))
                end
                if size(param) != (Int(N), F)
                    throw(DomainError(param, "If $param_name is a matrix, must have dimensions N×F = $N×$F."))
                end
                return [param[i, :] for i=1:Int(N)]
            end
        end

        # Process m and mlat, ensuring mlat = m - F*q^2*a/8
        mass_shift = F * q^2 * a / 8
        
        if isnothing(mlat)
            # mlat not provided: compute from m
            m_processed = process_mass_param(m, "m")
            mlat = if isinf(N)
                PeriodicArrays.PeriodicVector([m_processed[1] .- mass_shift])
            else
                [m_processed[i] .- mass_shift for i=1:Int(N)]
            end
            m = m_processed
        else
            # mlat provided: use it and compute m from mlat
            mlat = process_mass_param(mlat, "mlat")
            m = if isinf(N)
                PeriodicArrays.PeriodicVector([mlat[1] .+ mass_shift])
            else
                [mlat[i] .+ mass_shift for i=1:Int(N)]
            end
        end

        # Process mprime
        mprime = process_mass_param(mprime, "mprime")
        new(N, F, q, periodic, a, L, θ2π, m, mlat, mprime)
    end
end

Base.isfinite(lattice::Lattice) = isfinite(lattice.N)
Base.isinf(lattice::Lattice) = isinf(lattice.N)

