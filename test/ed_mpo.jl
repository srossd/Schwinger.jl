using ProgressMeter, StatsBase

@testset "Ground state and gap" begin
    @showprogress desc="ED vs MPO" for lat in [
        Lattice(10; F = 1, a = .7), 
        Lattice(10; F = 1, a = .7, m = 0.5), 
        Lattice(10; F = 1, a = .7, m = 0.5, θ2π = 0.3), 
        Lattice(10; F = 1, a = .7, m = 0.5, θ2π = 0.3, q = 4),
        Lattice(10; F = 1, a = .7, m = 0.5, θ2π = 0.3, q = 4, mprime = 0.8),
        Lattice(6;  F = 2, a = .7, m = 0.5, θ2π = 0.3), 
        Lattice(4;  F = 3, a = .7, m = 0.5, θ2π = 0.3),
        Lattice(4;  F = 3, a = .7, m = 0.5, θ2π = 0.3, periodic = true), # test wigner string, wilson loop
        Lattice(2;  F = 4, a = .7, m = 0.5, θ2π = 0.3)]

        for universe in 0:lat.q-1
            ed_gs = groundstate(Hamiltonian(lat, EDBackend(); universe=universe))
            itensors_gs = groundstate(Hamiltonian(lat, ITensorsBackend(); universe=universe), energy_tol=1E-10)
            mpskit_gs = lat.periodic ? missing : groundstate(Hamiltonian(lat, MPSKitBackend(); universe=universe), energy_tol=1E-10)

            ed_energy = energy(ed_gs)
            itensors_energy = energy(itensors_gs)
            @test ed_energy ≈ itensors_energy rtol=1E-6
            if !lat.periodic
                mpskit_energy = energy(mpskit_gs)
                @test ed_energy ≈ mpskit_energy rtol=1E-6
            end

            ed_scalar = scalar(ed_gs)
            itensors_scalar = scalar(itensors_gs)
            @test ed_scalar ≈ itensors_scalar rtol=1E-4
            if !lat.periodic
                mpskit_scalar = scalar(mpskit_gs)
                @test ed_scalar ≈ mpskit_scalar rtol=1E-4
            end

            ed_scalardensities = scalardensities(ed_gs)
            itensors_scalardensities = scalardensities(itensors_gs)
            @test ed_scalardensities ≈ itensors_scalardensities rtol=1E-4
            if !lat.periodic
                mpskit_scalardensities = scalardensities(mpskit_gs)
                @test ed_scalardensities ≈ mpskit_scalardensities rtol=1E-4
            end

            ed_pseudoscalar = pseudoscalar(ed_gs)
            itensors_pseudoscalar = pseudoscalar(itensors_gs)
            @test ed_pseudoscalar ≈ itensors_pseudoscalar rtol=1E-4
            if !lat.periodic
                mpskit_pseudoscalar = pseudoscalar(mpskit_gs)
                @test ed_pseudoscalar ≈ mpskit_pseudoscalar rtol=1E-4
            end

            ed_pseudoscalardensities = pseudoscalardensities(ed_gs)
            itensors_pseudoscalardensities = pseudoscalardensities(itensors_gs)
            @test ed_pseudoscalardensities ≈ itensors_pseudoscalardensities rtol=1E-4
            if !lat.periodic
                mpskit_pseudoscalardensities = pseudoscalardensities(mpskit_gs)
                @test ed_pseudoscalardensities ≈ mpskit_pseudoscalardensities rtol=1E-4
            end

            ed_efs = electricfields(ed_gs)
            itensors_efs = electricfields(itensors_gs)
            @test ed_efs ≈ itensors_efs rtol=1E-4
            if !lat.periodic
                mpskit_efs = electricfields(mpskit_gs)
                @test ed_efs ≈ mpskit_efs rtol=1E-4
            end

            ed_ents = entanglements(ed_gs)
            itensors_ents = entanglements(itensors_gs)
            @test ed_ents ≈ itensors_ents rtol=1E-4
            if !lat.periodic
                mpskit_ents = entanglements(mpskit_gs)
                @test ed_ents ≈ mpskit_ents rtol=1E-4
            end

            ed_avgE = expectation(AverageElectricField(lat, EDBackend(); power=1, universe=universe), ed_gs)
            itensors_avgE = expectation(AverageElectricField(lat, ITensorsBackend(); power=1, universe=universe), itensors_gs)
            @test ed_avgE ≈ itensors_avgE rtol=1E-4
            if !lat.periodic
                mpskit_avgE = expectation(AverageElectricField(lat, MPSKitBackend(); power=1, universe=universe), mpskit_gs)
                @test ed_avgE ≈ mpskit_avgE rtol=1E-4
            end

            sitelist = sample(1:Int(lat.N), rand(1:Int(lat.N)), replace=false)
            ed_avgE2_sub = real(expectation(AverageElectricField(lat, EDBackend(); power=2, universe=universe, sitelist=sitelist), ed_gs))
            itensors_avgE2_sub = real(expectation(AverageElectricField(lat, ITensorsBackend(); power=2, universe=universe, sitelist=sitelist), itensors_gs))
            @test ed_avgE2_sub ≈ itensors_avgE2_sub rtol=1E-3

            if lat.periodic
                ed_wilson = expectation(WilsonLoop(lat; backend = EDBackend(), universe=universe), ed_gs)
                itensors_wilson = expectation(WilsonLoop(lat; backend = ITensorsBackend(), universe=universe), itensors_gs)
                @test ed_wilson ≈ itensors_wilson rtol=1E-4
            end
            
            ed_gap = energygap(Hamiltonian(lat, EDBackend(); universe=universe))
            itensors_gap = energygap(Hamiltonian(lat, ITensorsBackend(); universe=universe))
            @test ed_gap ≈ itensors_gap rtol=1E-3
            if !lat.periodic
                mpskit_gap = energygap(Hamiltonian(lat, MPSKitBackend(); universe=universe))
                @test ed_gap ≈ mpskit_gap rtol=1E-3
            end
        end
    end
end