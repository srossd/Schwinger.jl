using ProgressMeter

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

        ed_gs = groundstate(EDHamiltonian(lat))
        mpo_gs = groundstate(MPOHamiltonian(lat))

        ed_energy = energy(ed_gs)
        mpo_energy = energy(mpo_gs)
        @test ed_energy ≈ mpo_energy rtol=1E-6

        ed_scalarvev = scalarvev(ed_gs)
        mpo_scalarvev = scalarvev(mpo_gs)
        @test ed_scalarvev ≈ mpo_scalarvev rtol=1E-4

        ed_pseudoscalarvev = pseudoscalarvev(ed_gs)
        mpo_pseudoscalarvev = pseudoscalarvev(mpo_gs)
        @test ed_pseudoscalarvev ≈ mpo_pseudoscalarvev rtol=1E-4

        ed_efs = electricfields(ed_gs)
        mpo_efs = electricfields(mpo_gs)
        @test ed_efs ≈ mpo_efs rtol=1E-4

        if lat.periodic
            ed_wilson = expectation(EDWilsonLoop(lat), ed_gs)
            mpo_wilson = expectation(MPOWilsonLoop(lat), mpo_gs)
            @test ed_wilson ≈ mpo_wilson rtol=1E-4
        end
        
        ed_gap = energygap(EDHamiltonian(lat))
        mpo_gap = energygap(MPOHamiltonian(lat))

        @test ed_gap ≈ mpo_gap rtol=1E-4
    end
end