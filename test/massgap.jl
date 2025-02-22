using ProgressMeter

@testset "Schwinger mass" begin
    gap = energygap(EDHamiltonian(SchwingerLattice{10,1}(periodic=true)))

    @test gap ≈ 1/√(π) rtol = 0.01

    @showprogress desc="OBC gap" for L in [10,15,20]
        gap = energygap(MPOHamiltonian(SchwingerLattice{20,1}(L=L)))

        @test gap ≈ √(1/π + (π/L)^2) rtol = 0.02
    end
end