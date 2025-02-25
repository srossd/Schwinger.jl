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
    L::Float64

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

        L = N*a

        if isa(θ2π, Real)
            θ2π = NTuple{N, Float64}((θ2π for _ in 1:N))
        end

        if isnothing(mlat)
            if isa(m, Real)
                mlat = NTuple{N, NTuple{F, Float64}}((NTuple{F,Float64}((m - F*q^2*a/8 for _ in 1:F)) for _ in 1:N))
            elseif isa(m, NTuple{F, Real})
                mlat = NTuple{N, NTuple{F, Float64}}((m .- F*q^2*a/8 for _ in 1:N))
            end
        elseif isa(mlat, Real)
            mlat = NTuple{N, NTuple{F, Float64}}((NTuple{F,Float64}((mlat for _ in 1:F)) for _ in 1:N))
        elseif isa(mlat, NTuple{F, Real})
            mlat = NTuple{N, NTuple{F, Float64}}((mlat for _ in 1:N))
        end

        if isa(mprime, Real)
            mprime = NTuple{N, NTuple{F, Float64}}((NTuple{F,Float64}((mprime for _ in 1:F)) for _ in 1:N))
        elseif isa(mprime, NTuple{F, Real})
            mprime = NTuple{N, NTuple{F, Float64}}((mprime for _ in 1:N))
        end
        
        new(q, periodic, a, L, θ2π, mlat, mprime)
    end
end