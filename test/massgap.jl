using ProgressMeter

@testset "Schwinger mass" begin
    gap = energygap(Hamiltonian(Lattice(10; periodic=true), EDBackend()))

    @test gap ≈ 1/√(π) rtol = 0.01

    @showprogress desc="OBC gap" for L in [10,15,20]
        itensors_gap = energygap(Hamiltonian(Lattice(20; L=L), ITensorsBackend()))
        # mpskit_gap = energygap(Hamiltonian(Lattice(20; L=L), MPSKitBackend()))

        # @test itensors_gap ≈ mpskit_gap rtol = 1E-6
        @test itensors_gap ≈ √(1/π + (π/L)^2) rtol = 0.02
    end
end