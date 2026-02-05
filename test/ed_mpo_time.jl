using ProgressMeter

@testset "Wilson line/loop at t = 0" begin
    @showprogress desc="ED vs MPO" for (lat, t) in [
        (Lattice(10; F = 1, a = .7, m = 0.5, θ2π = 0.3, periodic = true), .8), 
        (Lattice(10; F = 1, a = .5, m = 1, θ2π = 0.5, periodic = true), 1.2), 
        (Lattice(12; F = 1, a = .7, m = 0.5, θ2π = 0.3, q = 4, mprime = 0.8), .7), 
        (Lattice(12; F = 1, a = .7, m = 0.5, q = 4, mprime = -0.2), .8)]

        ed_gs = groundstate(Hamiltonian(lat, EDBackend()))
        itensor_gs = groundstate(Hamiltonian(lat, ITensorsBackend()); energy_tol=1E-12)
        mpskit_gs = lat.periodic ? missing : groundstate(Hamiltonian(lat, MPSKitBackend()); bonddim=100)

        ed_gsW = if lat.periodic
            act(WilsonLoop(lat; backend = EDBackend()), ed_gs)
        else
            act(WilsonLine(lat; backend = EDBackend()), ed_gs)
        end
        itensor_gsW = if lat.periodic
            act(WilsonLoop(lat; backend = ITensorsBackend()), itensor_gs)
        else
            act(WilsonLine(lat; backend = ITensorsBackend()), itensor_gs)
        end
        mpskit_gsW = lat.periodic ? missing : act(WilsonLine(lat; backend = MPSKitBackend()), mpskit_gs)

        ed_gsW_t, ed_obs = evolve(ed_gsW, t; nsteps=20, observable = Dict("energy" => (ψ, t) -> energy(ψ), "scalar" => (ψ, t) -> scalar(ψ), "electricfields" => (ψ, t) -> electricfields(ψ)))
        itensor_gsW_t, itensor_obs = evolve(itensor_gsW, t; nsteps=20, observable = Dict("energy" => (ψ, t) -> energy(ψ), "scalar" => (ψ, t) -> scalar(ψ), "electricfields" => (ψ, t) -> electricfields(ψ)))
        mpskit_gsW_t, mpskit_obs = lat.periodic ? (missing, missing) : evolve(mpskit_gsW, t; nsteps=20, observable = Dict("energy" => (ψ, t) -> energy(ψ), "scalar" => (ψ, t) -> scalar(ψ), "electricfields" => (ψ, t) -> electricfields(ψ)))

        @test energy(ed_gsW_t) == ed_obs.energy[20]
        @test energy(itensor_gsW_t) == itensor_obs.energy[20]
        if !lat.periodic
            @test energy(mpskit_gsW_t) == mpskit_obs.energy[20]
        end

        @test scalar(ed_gsW_t) == ed_obs.scalar[20]
        @test scalar(itensor_gsW_t) == itensor_obs.scalar[20]
        if !lat.periodic
            @test scalar(mpskit_gsW_t) == mpskit_obs.scalar[20]
        end

        # Energy conservation
        @test ed_obs.energy[1] ≈ ed_obs.energy[20] rtol=1E-10
        @test itensor_obs.energy[1] ≈ itensor_obs.energy[20] rtol=1E-10
        if !lat.periodic
            @test mpskit_obs.energy[1] ≈ mpskit_obs.energy[20] rtol=1E-10
        end

        #ED vs MPO
        @test ed_obs.energy ≈ itensor_obs.energy rtol=1E-4
        @test ed_obs.scalar ≈ itensor_obs.scalar rtol=1E-4
        @test ed_obs.electricfields ≈ itensor_obs.electricfields rtol=1E-4
        if !lat.periodic
            @test ed_obs.energy ≈ mpskit_obs.energy rtol=1E-4
            @test ed_obs.scalar ≈ mpskit_obs.scalar rtol=1E-4
            @test ed_obs.electricfields ≈ mpskit_obs.electricfields rtol=1E-4
        end
    end
end