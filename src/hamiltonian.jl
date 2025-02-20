"""
`hamiltonian(lattice)`

Computes the MPO Hamiltonian for the Schwinger model.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""
@memoize function hamiltonian(lattice::SchwingerLattice{N,F}; L_max::Int = 3) where {N,F}
    @unpack q, periodic, a, θ2π, mlat, mprime = lattice

    hamiltonian = OpSum()

    # Hopping term
    
    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += 1/(2a),"S+",ind,"S-",ind + F
            hamiltonian += 1/(2a),"S-",ind,"S+",ind + F
        end
    end

    if periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                hamiltonian += 1/(2a),"S+",ind,"S-",k,"raise",N * F + 1
                hamiltonian += 1/(2a),"S-",ind,"S+",k,"lower",N * F + 1
            else
                hamiltonian += 1/(2a),"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                hamiltonian += 1/(2a),"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    end

    # Real mass term

    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += mlat[j][k] * ((-1) ^ (j + 1)),"Sz",ind
        end
    end

    # Imaginary mass term

    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += (mprime[j][k]/2)*(-1)^j,"S+",ind,"S-",ind + F
            hamiltonian += (mprime[j][k]/2)*(-1)^j,"S-",ind,"S+",ind + F
        end
    end

    if periodic
        for k in 1:F
            ind = F * (N - 1) + k

            if F ≤ 2 # safe to ignore fermion parity factors
                hamiltonian += (mprime[N][k]/2),"S+",ind,"S-",k,"raise",N * F + 1
                hamiltonian += (mprime[N][k]/2),"S-",ind,"S+",k,"lower",N * F + 1
            else
                hamiltonian += (mprime[N][k]/2),"S+",ind,"S-",k,"raise",N * F + 1, wigner_string(N, F, k)...
                hamiltonian += (mprime[N][k]/2),"S-",ind,"S+",k,"lower",N * F + 1, wigner_string(N, F, k)...
            end
        end
    end

    # Gauge kinetic term

    if periodic
        hamiltonian += q^2 * a * N / 2,"L0",N * F + 1,"L0",N * F + 1
        hamiltonian += q * a * sum(θ2π),"L0",N * F + 1

        hamiltonian += q^2 * a * N * F / 4,"L0",N * F + 1
        
        for j in 1:N
            for k in 1:F
                ind = F*(j - 1) + k
                hamiltonian += q^2 * a * (N - j),"L0",N * F + 1,"Sz",ind
            end
        end
    end
    hamiltonian += a / 2 * sum(θ2π .^ 2),"Id",1
    hamiltonian += q * a * F / 4 * sum(θ2π),"Id",1
    for j in 1:N-1
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += q * a * sum(θ2π[j+1:N]),"Sz",ind
        end
    end

    hamiltonian += q^2 * a * N * F * F / 16,"Id",1
    for j in 1:N
        for k in 1:F
            ind = F*(j - 1) + k
            hamiltonian += q^2 * a * F * ((N + 1 - j) ÷ 2) / 2,"Sz",ind
        end
    end

    for j1 in 1:N
        for k1 in 1:F
            ind1 = F * (j1 - 1) + k1
            for j2 in j1:N
                for k2 in (j2 == j1 ? (k1:F) : 1:F)
                    ind2 = F*(j2 - 1) + k2
                    hamiltonian += q^2*a*(N - j2)/(j1 == j2 && k1 == k2 ? 2 : 1),"Sz",ind1,"Sz",ind2
                end
            end
        end
    end

    return MPO(hamiltonian, sites(lattice; L_max = L_max))
end