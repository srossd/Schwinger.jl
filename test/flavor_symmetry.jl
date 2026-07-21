using TensorKit: dim

# SU(F) flavor-symmetric MPSKit construction (Lattice(...; flavor_sym=true)).
# Validated against Schwinger.jl's independent exact diagonalization: the two use
# entirely different representations (flavor_sym fuses F flavors into one site as SU(F)
# multiplets; ED uses the occupation basis), so equal ground-state energies confirm the
# construction and its conventions. Kept light for CI: N=4, F=2, exact ED reference.

@testset "SU(2) flavor_sym ground state == ED" begin
    # N=4, F=2, m=0.5 (exercises mass + hopping + electric together). Open boundaries
    # fix the electric field by Gauss law, so ED (L_max=0) is exact.
    N, a, m = 4, 1.0, 0.5
    Eed = energy(loweststates(Hamiltonian(Lattice(N; F=2, m=m, a=a); backend=:ED), 1)[1])
    latfs = Lattice(N; F=2, m=m, a=a, flavor_sym=true)
    Efs = energy(loweststates(Hamiltonian(latfs; backend=:MPSKit), 1; energy_tol=1e-10, bonddim=24)[1])
    @test Efs ≈ Eed atol=1e-6
end

@testset "SU(3) flavor_sym Hamiltonian builds" begin
    lat = Lattice(4; F=3, m=0.3, flavor_sym=true)
    @test Hamiltonian(lat; backend=:MPSKit) isa MPSKitOperator
    @test length(Schwinger.get_mpskit_spaces(lat)) == 4      # one physical site per staggered site
end

@testset "flavor_sym argument validation" begin
    @test_throws ArgumentError Lattice(4; F=1, flavor_sym=true)               # needs F ≥ 2
    @test_throws ArgumentError Lattice(4; F=2, periodic=true, flavor_sym=true) # no periodic
    # unequal masses across flavors break the SU(F) symmetry -> error at build
    latbad = Lattice(4; F=2, mlat=[0.1 0.2; 0.1 0.2; 0.1 0.2; 0.1 0.2], flavor_sym=true)
    @test_throws ArgumentError Hamiltonian(latbad; backend=:MPSKit)
end

@testset "flavor-aware observables == ED" begin
    N, a, m = 4, 1.0, 0.5
    sed = loweststates(Hamiltonian(Lattice(N; F=2, m=m, a=a); backend=:ED), 1)[1]
    sfs = loweststates(Hamiltonian(Lattice(N; F=2, m=m, a=a, flavor_sym=true); backend=:MPSKit), 1;
                       energy_tol=1e-10, bonddim=24)[1]
    # occupations collapse to one number per site (the common per-flavor value)
    @test occupations(sfs) isa AbstractVector
    @test length(occupations(sfs)) == N
    @test occupations(sfs) ≈ occupations(sed)[:, 1] atol=1e-4        # per-flavor value
    # derived per-site observables agree with the (F-site) ED computation
    @test vec(charges(sfs)) ≈ vec(charges(sed)) atol=1e-4
    @test electricfields(sfs) ≈ electricfields(sed) atol=1e-4
    @test scalar(sfs) ≈ scalar(sed) atol=1e-4
    @test length(entanglements(sfs)) == N - 1
    # per-site energy_density decomposition is not implemented for flavor_sym
    @test_throws ArgumentError energy_density(sfs, 1)
end

@testset "flavor_irrep helpers & validation" begin
    lat = Lattice(4; F=2, flavor_sym=true)
    @test flavor_singlet(lat) == Schwinger._flavor_irrep(2, 0)              # SU(2) singlet
    @test dim(flavor_adjoint(lat)) == 3                           # SU(2) adjoint = triplet
    @test dim(flavor_adjoint(Lattice(4; F=3, flavor_sym=true))) == 8   # SU(3) adjoint = 8
    # flavor_irrep-resolved excitations require flavor_sym (here the lattice is not flavor_sym)
    @test_throws ArgumentError loweststates(Hamiltonian(Lattice(4; F=2, m=0.5); backend=:MPSKit), 2;
                                            flavor_irrep=flavor_adjoint(lat))
end

@testset "flavor-resolved excitations split == ED (finite)" begin
    # N=4, m=1: ED has a 3-fold degenerate isotriplet (adjoint) below a non-degenerate
    # isosinglet. The flavor_sym excitation, targeted by flavor_irrep, must land the adjoint
    # channel on the (lighter) triplet level and the singlet channel on the singlet level —
    # confirming the SU(2) flavor labelling of excited states, not just the ground state.
    N, a, m = 4, 1.0, 1.0
    Eed = sort!([real(energy(s)) for s in loweststates(Hamiltonian(Lattice(N; F=2, m=m, a=a); backend=:ED), 6)])
    E_triplet, E_singlet = Eed[2], Eed[5]      # Eed[2]==Eed[3]==Eed[4] (triplet), Eed[5] singlet
    @test Eed[2] ≈ Eed[4] atol=1e-6            # isotriplet is 3-fold degenerate in ED
    latf = Lattice(N; F=2, m=m, a=a, flavor_sym=true)
    Hf = Hamiltonian(latf; backend=:MPSKit)
    E_adj = real(energy(loweststates(Hf, 2; flavor_irrep=flavor_adjoint(latf), energy_tol=1e-10, bonddim=24)[2]))
    E_sng = real(energy(loweststates(Hf, 2; flavor_irrep=flavor_singlet(latf), energy_tol=1e-10, bonddim=24)[2]))
    @test E_adj ≈ E_triplet atol=1e-5
    @test E_sng ≈ E_singlet atol=1e-5
    @test E_adj < E_sng                        # triplet lighter than singlet at m/g = 1
end
