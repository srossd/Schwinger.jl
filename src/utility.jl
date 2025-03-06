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
tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

function wigner_string(N::Int, F::Int, flavor::Int)
    return tuplejoin((("Sz",(j-1)*F + flavor) for j in 1:N)...)
end

@memoize function sites(lattice::SchwingerLattice{N,F}; L_max::Int = 3) where {N,F}
    total_sites = lattice.periodic ? N * F + 1 : N * F
    sites = Vector{Index}(undef, total_sites)

    # Create each site in the SchwingerLattice
    for i in 1:total_sites
        if lattice.periodic && i == total_sites
            sites[i] = Index([QN() => 2L_max + 1]; tags="Boson,Site,n=$i")
        else
            sites[i] = siteind("S=1/2", i; conserve_qns = true)
        end
    end

    return sites
end

function sites(hamiltonian::MPOOperator{N,F}) where {N,F}
    return sites(hamiltonian.lattice; L_max = hamiltonian.L_max)
end

function constructoperator(lattice::SchwingerLattice{N,F}, action::Function; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    L_max, universe = process_L_max_universe(lattice, L_max, universe)

    states = schwingerbasis(lattice; L_max = L_max, universe = universe)
    positions = positionindex(lattice; L_max = L_max, universe = universe)

    I = Vector{Int}()
    J = Vector{Int}()
    V = Vector{ComplexF64}()
    for i in eachindex(states)
        for (key, val) in action(states[i])
            if haskey(positions,key)
                push!(I, positions[key])
                push!(J, i)
                push!(V, val)
            end
        end
    end

    matrix = sparse(I, J, V, length(states), length(states))
    return EDOperator(lattice, matrix, L_max, universe)
end

function process_L_max_universe(lattice::SchwingerLattice, L_max::Union{Nothing,Int}, universe::Int)
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += lattice.q
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    return L_max, universe
end