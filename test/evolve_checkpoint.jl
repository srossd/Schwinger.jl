# The `evolve` checkpoint hook fires every `checkpoint_every` steps and receives
# (state, current_time, step, observer-with-history-so-far).
@testset "evolve checkpoint hook" begin
    lat = Lattice(6; F = 1, a = 1, m = 0.5)
    gs  = groundstate(Hamiltonian(lat, EDBackend()))

    steps = Int[]; times = Float64[]; cp_energy = Float64[]; cp_histlen = Int[]
    _, obs = evolve(gs, 1.0; nsteps = 10,
        observable = Dict("e" => (ψ, t) -> real(energy(ψ))),
        checkpoint = (st, t, step, o) -> begin
            push!(steps, step); push!(times, t)
            push!(cp_energy, real(energy(st))); push!(cp_histlen, length(o.e))
        end,
        checkpoint_every = 5)

    @test steps == [5, 10]                                # fired every 5 steps
    @test times[end] ≈ 1.0                                # final checkpoint at the total time
    @test cp_histlen == [5, 10]                           # observer carries the history up to each checkpoint
    @test cp_energy ≈ [obs.e[5], obs.e[10]]               # the checkpointed state matches the trajectory

    # with no checkpoint given, evolution runs normally (no firing, no error)
    nfired = Ref(0)
    evolve(gs, 0.5; nsteps = 4, checkpoint = (args...) -> (nfired[] += 1), checkpoint_every = 10)
    @test nfired[] == 0                                   # checkpoint_every > nsteps => never fires
end
