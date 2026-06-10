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
    if !isnothing(observable)
        for (name, obs_fn) in observable
            observables[name] = (; kwargs...) -> obs_fn(EDState(state.hamiltonian, kwargs[:state], state.defects, state.net_charge), kwargs[:current_time])
        end
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    for (i, tnow) in enumerate(t/nsteps .* (1:nsteps))
        ψ = exponentiate(H, -1im*t/nsteps, ψ; ishermitian = true, tol = tol/nsteps)[1]
        Observers.update!(obs; step = i, current_time = tnow, state = ψ)
    end

    return EDState(state.hamiltonian, ψ, state.defects, state.net_charge), obs
end

"""
    evolve(state::ITensorState, t::Real; nsteps::Int = 1, observable = nothing, kwargs...)

Evolve an ITensor MPS state by imaginary time `t` using TDVP.

# Arguments
- `state::ITensorState`: The state to evolve.
- `t::Real`: The time to evolve by (will be multiplied by -im for imaginary time evolution).
- `nsteps::Int = 1`: Number of time steps for the TDVP algorithm.
- `maxlinkdim::Union{Nothing,Int} = nothing`: Maximum bond dimension to keep
  during the evolution (forwarded to ITensorMPS `tdvp` as `maxdim`).
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
function evolve(state::ITensorState, t::Real; nsteps::Int = 1,
                maxlinkdim::Union{Nothing,Int} = nothing,
                observable::Union{Nothing,Function,Dict} = nothing, kwargs...)
    H = state.hamiltonian.mpo

    step(; sweep) = sweep
    current_time(; current_time) = real(current_time/(-1im))

    if observable isa Function
        observable = Dict("observable" => observable)
    elseif isnothing(observable)
        observable = Dict{String,Function}()
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(ITensorState(state.hamiltonian, kwargs[:state], state.defects), kwargs[:current_time]/(-1im))
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    maxdim_kw = isnothing(maxlinkdim) ? (;) : (; maxdim = maxlinkdim)
    ψ0 = state.psi
    ψt = tdvp(H, -1im*t, ψ0; nsteps = nsteps, (step_observer!) = obs, maxdim_kw..., kwargs...,)

    return ITensorState(state.hamiltonian, ψt, state.defects), obs
end

"""
    evolve(state::MPSKitState, t::Real; nsteps::Int = 1, two_site = false,
           maxlinkdim = nothing, trscheme = nothing, observable = nothing, kwargs...)

Evolve an MPSKit MPS state by real time `t` using MPSKit's TDVP algorithm.

# Arguments
- `state::MPSKitState`: The state to evolve.
- `t::Real`: The time to evolve by.
- `nsteps::Int = 1`: Number of time steps to divide the evolution into.
- `two_site::Bool = false`: Use the two-site integrator (`TDVP2`) instead of the
  single-site `TDVP`. Two-site updates let the bond dimension grow and be
  truncated during the evolution (needed e.g. when the initial state has a small
  bond dimension or when entanglement grows); single-site `TDVP` preserves the
  bond dimension of the input state exactly.
- `maxlinkdim::Union{Nothing,Int} = nothing`: Maximum bond dimension to keep at
  each two-site truncation. Only used when `two_site = true`. If both `maxlinkdim`
  and `trscheme` are `nothing`, no truncation is applied.
- `trscheme = nothing`: A full TensorKit truncation scheme (e.g.
  `trunctol(; rtol = 1e-10)`), taking precedence over `maxlinkdim`. Only used when
  `two_site = true`.
- `observable::Union{Nothing,Function,Dict} = nothing`: Observable(s) to monitor.
  A single function `(state, t) -> value`, or a dictionary of name => function.
- `kwargs...`: Additional keyword arguments forwarded to `MPSKit.timestep`.

# Returns
- `MPSKitState`: The evolved state.
- `obs`: Observer object containing the history of observables and metadata.
"""
function evolve(state::MPSKitState, t::Real; nsteps::Int = 1, two_site::Bool = false,
                maxlinkdim::Union{Nothing,Int} = nothing, trscheme = nothing,
                observable::Union{Nothing,Function,Dict} = nothing, kwargs...)
    # With static defects, evolve in the absorbed representation (defect sites fused into
    # their matter neighbours) so the bond dimension can grow, then split them back out.
    if !isempty(state.defects)
        lat = state.hamiltonian.lattice
        fop = MPSKitOperator(lat, _fused_defect_lempo(lat, state.defects, state.hamiltonian.universe),
                             state.hamiltonian.universe, DefectCharge[])
        fstate = MPSKitState(fop, _fuse_defect_mps(state.psi, lat, state.defects), DefectCharge[])
        fevolved, obs = evolve(fstate, t; nsteps = nsteps, two_site = two_site,
                               maxlinkdim = maxlinkdim, trscheme = trscheme,
                               observable = observable, kwargs...)
        return MPSKitState(state.hamiltonian, _split_defect_mps(fevolved.psi, lat, state.defects),
                           state.defects), obs
    end

    H = state.hamiltonian.lempo

    ψ = copy(state.psi)

    step(; step) = step
    current_time(; current_time) = current_time

    if observable isa Function
        observable = Dict("observable" => observable)
    elseif isnothing(observable)
        observable = Dict{String,Function}()
    end

    observables = Dict{String,Function}()
    for (name, obs_fn) in observable
        observables[name] = (; kwargs...) -> obs_fn(MPSKitState(state.hamiltonian, kwargs[:state], state.defects), kwargs[:current_time])
    end

    obs = Observers.observer(
        merge(Dict("step" => step, "time" => current_time), observables)...
    )

    # Choose the integrator. TDVP2 (two-site) allows the bond dimension to grow
    # and be truncated; cap it with `maxlinkdim` (or a full `trscheme`).
    alg = if two_site
        scheme = !isnothing(trscheme)   ? trscheme :
                 !isnothing(maxlinkdim) ? truncrank(maxlinkdim) : notrunc()
        MPSKit.TDVP2(; trscheme = scheme)
    else
        (isnothing(maxlinkdim) && isnothing(trscheme)) ||
            @warn "maxlinkdim/trscheme are ignored by single-site TDVP; pass two_site=true to bound the bond dimension"
        MPSKit.TDVP()
    end

    for (i, tnow) in enumerate(t/nsteps .* (1:nsteps))
        @info "TDVP step $i / $nsteps (t = $tnow / $t)"
        ψ = MPSKit.timestep(ψ, H, tnow, t/nsteps, alg; kwargs...)[1]
        Observers.update!(obs; step = i, current_time = tnow, state = ψ)
    end

    return MPSKitState(state.hamiltonian, ψ, state.defects), obs

end