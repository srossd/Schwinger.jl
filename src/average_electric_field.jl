"""
`EDAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)`

Construct an `EDOperator` that computes the average electric field (raised to some power).

# Arguments
- `lattice::SchwingerLattice{N,F}`: The lattice to compute the average electric field on.
- `power::Int=1`: The power to raise the electric field to.
- `L_max::Union{Nothing,Int}=nothing`: The maximum absolute value of L₀.
- `universe::Int=0`: The universe to compute the average electric field in.
- `sitelist::Union{Nothing,Vector{Int}}=nothing`: List of sites to average over.
"""
function EDAverageElectricField(lattice::SchwingerLattice{N,F}; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0, sitelist::Union{Nothing,Vector{Int}} = nothing) where {N,F}
    if isnothing(sitelist)
        sitelist = collect(1:N)
    else
        if !all(x -> 1 <= x <= N, sitelist) || length(unique(sitelist)) != length(sitelist) && length(sitelist) < 1
            throw(ArgumentError("sitelist must be between 1 and $(N), unique and non-empty"))
        end
    end
    nsites = length(sitelist)

    function averagefield(state::SchwingerBasisState{N,F})
        return [(occupations(state), L₀(state)) => (sum(electricfields(state)[sitelist])/nsites)^power]
    end

    return constructoperator(lattice, averagefield; L_max = L_max, universe = universe)
end

"""
`MPOAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)`

Construct an `MPOOperator` that computes the average electric field (raised to some power).

# Arguments
- `lattice::SchwingerLattice{N,F}`: The lattice to compute the average electric field on.
- `power::Int=1`: The power to raise the electric field to.
- `L_max::Union{Nothing,Int}=nothing`: The maximum absolute value of L₀.
- `universe::Int=0`: The universe to compute the average electric field in.
- `sitelist::Union{Nothing,Vector{Int}}=nothing`: List of sites to average over.
"""
function MPOAverageElectricField(lattice::SchwingerLattice{N,F}; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0, sitelist::Union{Nothing,Vector{Int}} = nothing) where {N,F}
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end
    if isnothing(sitelist)
        sitelist = collect(1:N)
    else
        if !all(x -> 1 <= x <= N, sitelist) || length(unique(sitelist)) != length(sitelist) && length(sitelist) < 1
            throw(ArgumentError("sitelist must be between 1 and $(N), unique and non-empty"))
        end
    end
    nsites = length(sitelist)

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += q
    end
    θ2πu = lattice.θ2π .+ universe

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max
    
    opsum = OpSum()
    
    if lattice.periodic
        opsum += 1,"L0",N*F+1
    end
    opsum += sum(θ2πu[s] - lattice.q * F/2 * isodd(s) for s = sitelist) / nsites,"Id",1 

    for j in 1:N
        for f in 1:F
            ind = F*(j-1)+f
            opsum += lattice.q*sum(j <= s for s = sitelist) / nsites,"Sz",ind
        end
    end

    mpo1 = MPO(opsum, sites(lattice; L_max=L_max))
    mpo = foldl(ITensors.apply, (mpo1 for _ in 1:power))

    return MPOOperator(lattice, mpo, L_max, universe)
end