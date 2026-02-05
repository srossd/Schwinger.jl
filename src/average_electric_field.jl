"""
`EDAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)`

Construct an `EDOperator` that computes the average electric field (raised to some power).

# Arguments
- `lattice::Lattice`: The lattice to compute the average electric field on.
- `power::Int=1`: The power to raise the electric field to.
- `L_max::Union{Nothing,Int}=nothing`: The maximum absolute value of L₀.
- `universe::Int=0`: The universe to compute the average electric field in.
- `sitelist::Union{Nothing,Vector{Int}}=nothing`: List of sites to average over.
"""
function EDAverageElectricField(lattice::Lattice; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0, sitelist::Union{Nothing,Vector{Int}} = nothing) 
    N = Int(lattice.N)
    if isnothing(sitelist)
        sitelist = collect(1:N)
    else
        if !all(x -> 1 <= x <= N, sitelist) || length(unique(sitelist)) != length(sitelist) && length(sitelist) < 1
            throw(ArgumentError("sitelist must be between 1 and $(N), unique and non-empty"))
        end
    end
    nsites = length(sitelist)

    function averagefield(state::BasisState)
        return [(occupations(state), L₀(state)) => (sum(electricfields(state)[sitelist])/nsites)^power]
    end

    return constructoperator(lattice, averagefield; L_max = L_max, universe = universe)
end

"""
`ITensorAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)`

Construct an `ITensorOperator` that computes the average electric field (raised to some power).

# Arguments
- `lattice::Lattice`: The lattice to compute the average electric field on.
- `power::Int=1`: The power to raise the electric field to.
- `L_max::Union{Nothing,Int}=nothing`: The maximum absolute value of L₀.
- `universe::Int=0`: The universe to compute the average electric field in.
- `sitelist::Union{Nothing,Vector{Int}}=nothing`: List of sites to average over.
"""
function ITensorAverageElectricField(lattice::Lattice; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0, sitelist::Union{Nothing,Vector{Int}} = nothing) 
    N, F = Int(lattice.N), lattice.F
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

    L_max, universe = process_L_max_universe(lattice, L_max, universe)
    θ2πu = lattice.θ2π .+ universe
    
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

    mpo1 = ITensorMPS.MPO(opsum, get_sites(lattice; L_max=L_max))
    mpo = foldl(ITensors.apply, (mpo1 for _ in 1:power))

    return ITensorOperator(lattice, mpo, L_max, universe)
end

# =============================================================================
# MPSKit Backend
# =============================================================================

"""
`MPSKitAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)`

Construct an `MPSKitOperator` that computes the average electric field (raised to some power).

# Arguments
- `lattice::Lattice`: The lattice to compute the average electric field on.
- `power::Int=1`: The power to raise the electric field to.
- `universe::Int=0`: The universe to compute the average electric field in.
- `sitelist::Union{Nothing,Vector{Int}}=nothing`: List of sites to average over.
"""
function MPSKitAverageElectricField(lattice::Lattice; power::Int = 1, universe::Int = 0, sitelist::Union{Nothing,Vector{Int}} = nothing) 
    if power != 1
        throw(ArgumentError("MPSKitAverageElectricField currently only supports power=1"))
    end
    _, universe = process_L_max_universe(lattice, 0, universe)

    N, F = lattice.N, lattice.F
    if isnothing(sitelist)
        if isinf(N)
            throw(ArgumentError("sitelist must be specified for infinite lattices"))
        end
        sitelist = collect(1:Int(N))
    else
        if !all(x -> 1 <= x <= (isinf(N) ? 2 : Int(N)), sitelist) || length(unique(sitelist)) != length(sitelist) && length(sitelist) < 1
            throw(ArgumentError("sitelist must be between 1 and $(isinf(N) ? 2 : Int(N)), unique and non-empty"))
        end
    end

    # Build the average electric field operator
    θ2πu = lattice.θ2π .+ universe
    nsites = length(sitelist)

    link_fcts = [n in sitelist ? (r::U1Irrep -> (r.charge + θ2πu[n])/nsites) : missing for n in 1:(isinf(N) ? 2 : Int(N))]
    if isinf(N)
        mpo = InfiniteLEMPOHamiltonian(get_mpskit_spaces(lattice), vcat(fill(missing, F - 1), [link_fcts[1]], fill(missing, F - 1), [link_fcts[2]]))
    else
        all_link_fcts = vcat([vcat(fill(missing, F - 1), [link_fcts[n]]) for n in 1:Int(N)]...)
        mpo = FiniteLEMPOHamiltonian(get_mpskit_spaces(lattice), all_link_fcts)
    end
    return MPSKitOperator(lattice, mpo, universe)
end