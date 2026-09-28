# §2.2 (two-point correlator) and §3.3 (chargeprofile).
using Schwinger, LinearAlgebra, Test

@testset "correlator2pt + chargeprofile" begin
    N = 8
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.5)
    H   = Hamiltonian(lat, EDBackend())
    gs  = groundstate(H)
    E0  = energy(gs)
    ts  = [0.0, 0.3, 0.7, 1.1]

    @testset "H,H correlator is the constant E0^2 (phase cancellation)" begin
        # H|0⟩ = E0|0⟩, so A(t)=H(t)=H and C(t) = ⟨0|H e^{-iHt} H|0⟩ e^{iE0 t} = E0^2 for all t.
        C = correlator2pt(gs, H, H, ts)
        @test all(c -> isapprox(real(c), E0^2; atol = 1e-6), C)
        @test all(c -> abs(imag(c)) < 1e-6, C)
        # connected subtracts ⟨H⟩² = E0² ⇒ identically zero
        Cc = correlator2pt(gs, H, H, ts; connected = true)
        @test all(c -> abs(c) < 1e-6, Cc)
    end

    @testset "C(0) equals the equal-time operator product ⟨A B⟩" begin
        M = Mass(lat; backend = EDBackend(), bare = false)   # Hermitian
        C0 = correlator2pt(gs, M, M, [0.0])[1]
        @test isapprox(C0, expectation(M * M, gs); atol = 1e-8)   # independent path (operator product)
        # non-trivial time dependence: a Hermitian non-conserved M does evolve
        Ct = correlator2pt(gs, M, M, ts)
        @test abs(Ct[end] - Ct[1]) > 1e-4
    end

    @testset "argument checks" begin
        @test correlator2pt(gs, H, H, Float64[]) == ComplexF64[]
        @test_throws ArgumentError correlator2pt(gs, H, H, [0.5, 0.2])   # not sorted
        @test_throws ArgumentError correlator2pt(gs, H, H, [-0.1, 0.2])  # negative time
    end

    @testset "chargeprofile" begin
        cp = chargeprofile(gs)
        @test cp isa Vector{Float64}
        @test length(cp) == N
        @test cp ≈ real.(vec(charges(gs))) atol = 1e-12
    end
end
