"""
evolve(state::SchwingerEDState{N,F}, t::Real, observable::Function)

Evolve the state by a time `t`, monitoring an observable.

# Arguments
- `state::SchwingerEDState{N,F}`: The state to evolve.
- `t::Real`: The time to evolve by.
- `observable::Function`: The observable `(ψ, t) -> obs` to monitor.
"""
function evolve(state::SchwingerEDState{N,F}, t::Real; nsteps::Int = 1, tol = 1E-12, observable::Union{Nothing,Function,Dict} = nothing, kwargs...) where {N,F}
    H = state.hamiltonian.matrix
    ψ = copy(state.coeffs)

    step(; step) = step
    current_time(; current_time) = current_time

    if observable isa Function
        observable = Dict("observable" => observable)
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(SchwingerEDState(state.hamiltonian, kwargs[:state]), kwargs[:current_time])
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    for (i, tnow) in enumerate(t/nsteps .* (1:nsteps))
        ψ = exponentiate(H, -1im*t/nsteps, ψ; ishermitian = true, tol = tol/nsteps)[1]
        Observers.update!(obs; step = i, current_time = tnow, state = ψ)
    end

    return SchwingerEDState(state.hamiltonian, ψ), obs
end

"""
evolve(state::SchwingerMPS{N,F}, t::Real)

Evolve the state by a time `t`.

# Arguments
- `state::SchwingerMPS{N,F}`: The state to evolve.
- `t::Real`: The time to evolve by.
"""
function evolve(state::SchwingerMPS{N,F}, t::Real; nsteps::Int = 1, observable::Union{Nothing,Function,Dict} = nothing, kwargs...) where {N,F}
    H = state.hamiltonian.mpo

    step(; sweep) = sweep
    current_time(; current_time) = real(current_time/(-1im))

    if observable isa Function
        observable = Dict("observable" => observable)
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(SchwingerMPS(state.hamiltonian, kwargs[:state]), kwargs[:current_time]/(-1im))
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )
    
    ψ0 = state.psi
    ψt = tdvp(H, -1im*t, ψ0; nsteps = nsteps, (step_observer!) = obs, kwargs...,)

    return SchwingerMPS(state.hamiltonian, ψt), obs
end