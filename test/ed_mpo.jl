using ProgressMeter

@testset "Ground state and gap" begin
    @showprogress desc="ED vs MPO" for lat in [
        SchwingerLattice{10,1}(a = .7), 
        SchwingerLattice{10,1}(a = .7, m = 0.5), 
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3), 
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3, q = 4),
        SchwingerLattice{10,1}(a = .7, m = 0.5, θ2π = 0.3, q = 4, mprime = 1.2),
        SchwingerLattice{6,2}(a = .7, m = 0.5, θ2π = 0.3), 
        SchwingerLattice{4,3}(a = .7, m = 0.5, θ2π = 0.3), 
        SchwingerLattice{2,4}(a = .7, m = 0.5, θ2π = 0.3)]
        ed_gs = groundstate(EDHamiltonian(lat))
        mpo_gs = groundstate(MPOHamiltonian(lat))

        ed_energy = energy(ed_gs)
        mpo_energy = energy(mpo_gs)
        @test ed_energy ≈ mpo_energy rtol=1E-6

        ed_efs = electricfields(ed_gs)
        mpo_efs = electricfields(mpo_gs)
        @test ed_efs ≈ mpo_efs rtol=1E-6
        
        ed_gap = energygap(EDHamiltonian(lat))
        mpo_gap = energygap(MPOHamiltonian(lat); energy_tol=1E-10)

        @test ed_gap ≈ mpo_gap rtol=1E-6
    end
end