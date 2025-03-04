function EDAverageElectricField(lattice::SchwingerLattice{N,F}; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    function averagefield(state::SchwingerBasisState{N,F})
        return [(occupations(state), L₀(state)) => mean(electricfields(state))^power]
    end

    return constructoperator(lattice, averagefield; L_max = L_max, universe = universe)
end

function MPOAverageElectricField(lattice::SchwingerLattice{N,F}; power::Int = 1, L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
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
    θ2πu = lattice.θ2π .+ universe

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max
    
    opsum = OpSum()
    
    if lattice.periodic
        opsum += 1,"L0",N*F+1
    end
    opsum += mean(θ2πu)+1/4,"Id",1 # 1/4 to account for shift terms in the Gauss law not included below

    for j in 1:N
        for f in 1:F
            ind = F*(j-1)+f
            opsum += (N-j)/N,"Sz",ind
        end
    end

    mpo1 = MPO(opsum, sites(lattice; L_max=L_max))
    mpo = foldl(ITensors.apply, (mpo1 for _ in 1:power))

    return MPOOperator(lattice, mpo, L_max, universe)
end