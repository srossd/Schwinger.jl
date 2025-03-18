function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, state::SchwingerMPS)
    g = create_group(parent, name)
    attributes(g)["type"] = "SchwingerMPS"
    attributes(g)["version"] = 1
    write(g, "hamiltonian", state.hamiltonian)
    write(g, "psi", state.psi)
    return g
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, state::SchwingerEDState)
    g = create_group(parent, name)
    attributes(g)["type"] = "SchwingerEDState"
    attributes(g)["version"] = 1
    write(g, "hamiltonian", state.hamiltonian)
    write(g, "coeffs", state.coeffs)
    return g
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, op::MPOOperator)
    g = create_group(parent, name)
    attributes(g)["type"] = "MPOOperator"
    attributes(g)["version"] = 1
    write(g, "lattice", op.lattice)
    write(g, "mpo", op.mpo)
    write(g, "L_max", op.L_max)
    write(g, "universe", op.universe)
    return g
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, op::EDOperator)
    g = create_group(parent, name)
    attributes(g)["type"] = "EDOperator"
    attributes(g)["version"] = 1
    write(g, "lattice", op.lattice)
    write(g, "matrix", op.matrix)
    write(g, "L_max", op.L_max)
    write(g, "universe", op.universe)
    return g
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, lattice::SchwingerLattice{N,F}) where {N,F}
    g = create_group(parent, name)
    attributes(g)["type"] = "SchwingerLattice"
    attributes(g)["version"] = 1
    write(g, "N", N)
    write(g, "F", F)
    write(g, "q", lattice.q)
    write(g, "periodic", lattice.periodic)
    write(g, "a", lattice.a)
    write(g, "theta2pi", lattice.θ2π)
    write(g, "mlat", lattice.mlat)
    write(g, "mprime", lattice.mprime)
    return g
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ::Type{SchwingerMPS})
    g = open_group(parent, name)
    hamiltonian = read(g, "hamiltonian", MPOOperator)
    psi = read(g, "psi", MPS)
    return SchwingerMPS(hamiltonian, psi)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ::Type{SchwingerEDState})
    g = open_group(parent, name)
    hamiltonian = read(g, "hamiltonian", EDOperator)
    coeffs = read(g, "coeffs")
    return SchwingerEDState(hamiltonian, coeffs)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ::Type{MPOOperator})
    g = open_group(parent, name)
    lattice = read(g, "lattice", SchwingerLattice)
    mpo = read(g, "mpo", MPO)
    L_max = read(g, "L_max")
    universe = read(g, "universe")
    return MPOOperator(lattice, mpo, L_max, universe)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ::Type{EDOperator})
    g = open_group(parent, name)
    lattice = read(g, "lattice", SchwingerLattice)
    matrix = read(g, "matrix")
    L_max = read(g, "L_max")
    universe = read(g, "universe")
    return EDOperator(lattice, matrix, L_max, universe)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ::Type{SchwingerLattice})
    g = open_group(parent, name)
    N = read(g, "N")
    F = read(g, "F")
    q = read(g, "q")
    periodic = read(g, "periodic")
    a = read(g, "a")
    θ2π = read(g, "theta2pi")
    mlat = read(g, "mlat")
    mprime = read(g, "mprime")
    return SchwingerLattice{N,F}(q = q, periodic = periodic, a = a, θ2π = θ2π, mlat = mlat, mprime = mprime)
end