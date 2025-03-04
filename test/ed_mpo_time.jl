using ProgressMeter

@testset "Wilson line/loop at t = 0" begin
    @showprogress desc="ED vs MPO" for (lat, t) in [
        (SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3, periodic = true), .8), 
        (SchwingerLattice{10,1}(a = .5, m = 1, θ2π = 0.5, periodic = true), 1.2), 
        (SchwingerLattice{12,1}(a = .7, m = 0.5, θ2π = 0.3, q = 4, mprime = 0.8), .7), 
        (SchwingerLattice{12,1}(a = .7, m = 0.5, q = 4, mprime = -0.2), .8)]

        ed_gs = groundstate(EDHamiltonian(lat))
        mpo_gs = groundstate(MPOHamiltonian(lat); energy_tol=1E-12)

        ed_gsW = if lat.periodic
            act(EDWilsonLoop(lat), ed_gs)
        else
            act(EDWilsonLine(lat), ed_gs)
        end
        mpo_gsW = if lat.periodic
            act(MPOWilsonLoop(lat), mpo_gs)
        else
            act(MPOWilsonLine(lat), mpo_gs)
        end

        ed_gsW_t, ed_obs = evolve(ed_gsW, t; nsteps=20, observable = Dict("energy" => (ψ, t) -> energy(ψ), "scalarvev" => (ψ, t) -> scalarvev(ψ), "electricfields" => (ψ, t) -> electricfields(ψ)))
        mpo_gsW_t, mpo_obs = evolve(mpo_gsW, t; nsteps=20, observable = Dict("energy" => (ψ, t) -> energy(ψ), "scalarvev" => (ψ, t) -> scalarvev(ψ), "electricfields" => (ψ, t) -> electricfields(ψ)))

        @test energy(ed_gsW_t) == ed_obs.energy[20]
        @test energy(mpo_gsW_t) == mpo_obs.energy[20]

        @test scalarvev(ed_gsW_t) == ed_obs.scalarvev[20]
        @test scalarvev(mpo_gsW_t) == mpo_obs.scalarvev[20]

        # Energy conservation
        @test ed_obs.energy[1] ≈ ed_obs.energy[20] rtol=1E-10
        @test mpo_obs.energy[1] ≈ mpo_obs.energy[20] rtol=1E-10

        #ED vs MPO
        @test ed_obs.energy ≈ mpo_obs.energy rtol=1E-4
        @test ed_obs.scalarvev ≈ mpo_obs.scalarvev rtol=1E-4
        @test ed_obs.electricfields ≈ mpo_obs.electricfields rtol=1E-4
    end
end