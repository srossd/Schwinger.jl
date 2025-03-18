"""
`SchwingerLattice{N,F}(;kwargs...)`

Constructs a SchwingerLattice for the Schwinger model.

# Arguments
- `periodic::Bool=false`: Whether the lattice is periodic.
- `q::Int=1`: Charge.
- `L::Union{Nothing,Real}=nothing`: Length of the lattice.
- `a::Union{Nothing,Real}=nothing`: Lattice spacing.
- `m::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=0.`: Mass parameter.
- `mlat::Union{Nothing,Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=nothing`: Local mass parameter.
- `mprime::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=0.`: Prime mass parameter.
- `θ2π::Union{Real,NTuple{N,Real}}=0.`: Theta angle.

# Returns
A `SchwingerLattice` object.

"""
struct SchwingerLattice{N,F}
    # global parameters
    N::Int
    F::Int
    q::Int
    periodic::Bool
    a::Float64
    L::Float64

    # local parameters (typically position-independent)
    θ2π::Vector{Float64}
    m::Array{Float64,2}
    mlat::Array{Float64,2}
    mprime::Array{Float64,2}

    function SchwingerLattice{N,F}(;
        q::Int = 1, periodic::Bool=false, 
        L::Union{Nothing,Real}=nothing, a::Union{Nothing,Real}=nothing, 
        m::Union{Real,NTuple{F,<:Real},AbstractVector{<:Real},AbstractArray{<:Real,2}}=0., 
        mlat::Union{Nothing,Real,NTuple{F,<:Real},AbstractVector{<:Real},AbstractArray{<:Real,2}}=nothing, 
        mprime::Union{Real,NTuple{F,<:Real},AbstractVector{<:Real},AbstractArray{<:Real,2}}=0., 
        θ2π::Union{Real,AbstractVector{<:Real}}=0.) where {N,F}

        if !iseven(N) || N < 2
            throw(DomainError(N, "Number of sites must be even and ≥ 2."))
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

        if !isnothing(L)
            if !isnothing(a)
                throw(ArgumentError("Cannot specify both lattice spacing and length."))
            end
            a = L/N
        end

        L = N*a

        if isa(θ2π, Real)
            θ2π = [θ2π for _ in 1:N]
        end

        if isa(m, NTuple{F,<:Real})
            m = [m[i] for i in 1:F]
        end

        if isnothing(mlat)
            mlat = m .- F*q^2*a/8 # the mass shift
        end
            
        if isa(mlat, Real)
            mlat = [mlat for _ in 1:N, _ in 1:F]
        elseif isa(mlat, AbstractVector{<:Real})
            if length(mlat) == F
            mlat = stack([mlat for _ in 1:N], dims=1)
            elseif length(mlat) == N
            mlat = stack([mlat for _ in 1:F], dims=2)
            else
            throw(DomainError(mlat, "If mass parameter is a vector, must have length F = $F or N = $N."))
            end
        else
            if size(mlat) != (N,F)
            throw(DomainError(m, "If mass parameter is a matrix, must have dimensions N×F = $N×$F."))
            end
        end

        m = mlat .+ F*q^2*a/8

        if isa(mprime, NTuple{F,<:Real})
            mprime = [mprime[i] for i in 1:F]
        end
        if isa(mprime, Real)
            mprime = [mprime for _ in 1:N, _ in 1:F]
        elseif isa(mprime, AbstractVector{<:Real})
            if length(mprime) == F
            mprime = stack([mprime for _ in 1:N], dims=1)
            elseif length(mprime) == N
            mprime = stack([mprime for _ in 1:F], dims=2)
            else
            throw(DomainError(mprime, "If mass parameter is a vector, must have length F = $F or N = $N."))
            end
        else
            if size(mprime) != (N,F)
            throw(DomainError(mprime, "If mass parameter is a matrix, must have dimensions N×F = $N×$F."))
            end
        end

        new(N, F, q, periodic, a, L, θ2π, m, mlat, mprime)
    end
end