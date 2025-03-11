using ProgressMeter, StatsBase

@testset "Ground state and gap" begin
    @showprogress desc="ED vs MPO" for lat in [
        SchwingerLattice{10,1}(a = .7), 
        SchwingerLattice{10,1}(a = .7, m = 0.5), 
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3), 
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3, q = 4),
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3, q = 4, mprime = 0.8),
        SchwingerLattice{6,2}(a = .7, m = 0.5, θ2π = 0.3), 
        SchwingerLattice{4,3}(a = .7, m = 0.5, θ2π = 0.3),
        SchwingerLattice{4,3}(a = .7, m = 0.5, θ2π = 0.3, periodic = true), # test wigner string, wilson loop
        SchwingerLattice{2,4}(a = .7, m = 0.5, θ2π = 0.3)]

        for universe in 0:lat.q-1
            ed_gs = groundstate(EDHamiltonian(lat; universe=universe))
            mpo_gs = groundstate(MPOHamiltonian(lat; universe=universe))

            ed_energy = energy(ed_gs)
            mpo_energy = energy(mpo_gs)
            @test ed_energy ≈ mpo_energy rtol=1E-6

            ed_scalar = scalar(ed_gs)
            mpo_scalar = scalar(mpo_gs)
            @test ed_scalar ≈ mpo_scalar rtol=1E-4

            ed_scalardensities = scalardensities(ed_gs)
            mpo_scalardensities = scalardensities(mpo_gs)
            @test ed_scalardensities ≈ mpo_scalardensities rtol=1E-4

            ed_pseudoscalar = pseudoscalar(ed_gs)
            mpo_pseudoscalar = pseudoscalar(mpo_gs)
            @test ed_pseudoscalar ≈ mpo_pseudoscalar rtol=1E-4

            ed_pseudoscalardensities = pseudoscalardensities(ed_gs)
            mpo_pseudoscalardensities = pseudoscalardensities(mpo_gs)
            @test ed_pseudoscalardensities ≈ mpo_pseudoscalardensities rtol=1E-4

            ed_efs = electricfields(ed_gs)
            mpo_efs = electricfields(mpo_gs)
            @test ed_efs ≈ mpo_efs rtol=1E-4

            ed_ents = entanglements(ed_gs)
            mpo_ents = entanglements(mpo_gs)
            @test ed_ents ≈ mpo_ents rtol=1E-4

            ed_avgE2 = expectation(EDAverageElectricField(lat; power=2, universe=universe), ed_gs)
            mpo_avgE2 = expectation(MPOAverageElectricField(lat; power=2, universe=universe), mpo_gs)
            @test ed_avgE2 ≈ mpo_avgE2 rtol=1E-4

            sitelist = sample(1:lat.N, rand(1:lat.N), replace=false)
            ed_avgE2_sub = real(expectation(EDAverageElectricField(lat; power=2, universe=universe, sitelist=sitelist), ed_gs))
            mpo_avgE2_sub = real(expectation(MPOAverageElectricField(lat; power=2, universe=universe, sitelist=sitelist), mpo_gs))
            @test ed_avgE2_sub ≈ mpo_avgE2_sub rtol=1E-4

            if lat.periodic
                ed_wilson = expectation(EDWilsonLoop(lat; universe=universe), ed_gs)
                mpo_wilson = expectation(MPOWilsonLoop(lat; universe=universe), mpo_gs)
                @test ed_wilson ≈ mpo_wilson rtol=1E-4
            end
            
            ed_gap = energygap(EDHamiltonian(lat; universe=universe))
            mpo_gap = energygap(MPOHamiltonian(lat; universe=universe))

            @test ed_gap ≈ mpo_gap rtol=1E-3
        end
    end
end