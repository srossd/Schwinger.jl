"""
    evolve(state::EDState, t::Real; nsteps::Int = 1, tol = 1E-12, observable = nothing, kwargs...)

Evolve an exact diagonalization state by time `t` using matrix exponentiation.

# Arguments
- `state::EDState`: The state to evolve.
- `t::Real`: The time to evolve by.
- `nsteps::Int = 1`: Number of time steps to divide the evolution into.
- `tol = 1E-12`: Tolerance for the matrix exponentiation.
- `observable::Union{Nothing,Function,Dict} = nothing`: Observable(s) to monitor during evolution.
  Can be a single function `(state, t) -> value` or a dictionary of name => function pairs.

# Returns
- `EDState`: The evolved state.
- `obs`: Observer object containing the history of observables and metadata.

# Examples
```julia
evolved_state, obs = evolve(state, 1.0; nsteps=10, observable=s->energy(s))
```
"""
function evolve(state::EDState, t::Real; nsteps::Int = 1, tol = 1E-12, observable::Union{Nothing,Function,Dict} = nothing) 
    H = state.hamiltonian.matrix
    ψ = copy(state.coeffs)

    step(; step) = step
    current_time(; current_time) = current_time

    if observable isa Function
        observable = Dict("observable" => observable)
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(EDState(state.hamiltonian, kwargs[:state]), kwargs[:current_time])
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    for (i, tnow) in enumerate(t/nsteps .* (1:nsteps))
        ψ = exponentiate(H, -1im*t/nsteps, ψ; ishermitian = true, tol = tol/nsteps)[1]
        Observers.update!(obs; step = i, current_time = tnow, state = ψ)
    end

    return EDState(state.hamiltonian, ψ), obs
end

"""
    evolve(state::ITensorState, t::Real; nsteps::Int = 1, observable = nothing, kwargs...)

Evolve an ITensor MPS state by imaginary time `t` using TDVP.

# Arguments
- `state::ITensorState`: The state to evolve.
- `t::Real`: The time to evolve by (will be multiplied by -im for imaginary time evolution).
- `nsteps::Int = 1`: Number of time steps for the TDVP algorithm.
- `observable::Union{Nothing,Function,Dict} = nothing`: Observable(s) to monitor during evolution.
  Can be a single function `(state, t) -> value` or a dictionary of name => function pairs.
- `kwargs...`: Additional keyword arguments passed to the TDVP algorithm.

# Returns
- `ITensorState`: The evolved state.
- `obs`: Observer object containing the history of observables and metadata.

# Notes
- Uses ITensors.jl TDVP algorithm.
- Time is tracked as real values (divided by -im) in the observer.
"""
function evolve(state::ITensorState, t::Real; nsteps::Int = 1, observable::Union{Nothing,Function,Dict} = nothing, kwargs...) 
    H = state.hamiltonian.mpo

    step(; sweep) = sweep
    current_time(; current_time) = real(current_time/(-1im))

    if observable isa Function
        observable = Dict("observable" => observable)
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(ITensorState(state.hamiltonian, kwargs[:state]), kwargs[:current_time]/(-1im))
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    ψ0 = state.psi
    ψt = tdvp(H, -1im*t, ψ0; nsteps = nsteps, (step_observer!) = obs, kwargs...,)

    return ITensorState(state.hamiltonian, ψt), obs
end

"""
    evolve(state::MPSKitState, t::Real; nsteps::Int = 1, dt = nothing, observable = nothing, kwargs...)

Evolve an MPSKit MPS state by time `t` using MPSKit's TDVP algorithm.

# Arguments
- `state::MPSKitState`: The state to evolve.
- `t::Real`: The time to evolve by.
- `nsteps::Int = 1`: Number of time steps to divide the evolution into.
- `dt::Union{Nothing,Real} = nothing`: Time step size (currently unused, evolution uses t/nsteps).
- `observable::Union{Nothing,Function,Dict} = nothing`: Observable(s) to monitor during evolution.
  Can be a single function `(state, t) -> value` or a dictionary of name => function pairs.
- `kwargs...`: Additional keyword arguments.

# Returns
- `MPSKitState`: The evolved state.
- `obs`: Observer object containing the history of observables and metadata.

# Notes
- Uses MPSKit.jl TDVP algorithm.
- Evolves with imaginary time (-im*t/nsteps).
"""
function evolve(state::MPSKitState, t::Real; nsteps::Int = 1, observable::Union{Nothing,Function,Dict} = nothing, kwargs...) 
    H = state.hamiltonian.lempo

    ψ = copy(state.psi)

    step(; step) = step
    current_time(; current_time) = current_time

    if observable isa Function
        observable = Dict("observable" => observable)
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(MPSKitState(state.hamiltonian, kwargs[:state]), kwargs[:current_time])
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    alg = MPSKit.TDVP()
    for (i, tnow) in enumerate(t/nsteps .* (1:nsteps))
        ψ = MPSKit.timestep(ψ, H, tnow, t/nsteps, alg; kwargs...)[1]
        Observers.update!(obs; step = i, current_time = tnow, state = ψ)
    end

    return MPSKitState(state.hamiltonian, ψ), obs

end