import libscf
import numpy as np
import pyscf

h = 1e-6

# H2O

mol = pyscf.gto.M(atom='''
                    O   -0.0000000   -0.1113512    0.0000000
                    H    0.0000000    0.4454047   -0.7830363
                    H   -0.0000000    0.4454047    0.7830363''',
                    basis='sto-3g')
natm, nbas, nelec, nshells = 3, 5, 10, 7

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

P = libscf.RHF(natm, nbas, nelec, nshells, atm, bas, env)
print("P:")
libscf.pmat(P)

print("denv:")
for i in range(32, 56):
    env[i] += h
    E1 = libscf.energy(natm, nbas, nshells, atm, bas, env, P)

    env[i] -= 2.0*h
    E2 = libscf.energy(natm, nbas, nshells, atm, bas, env, P)

    print(round(-(E2 - E1)/(2.0*h), 6), end=' ')

    env[i] += h

print()