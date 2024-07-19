import libscf
import numpy as np
import pyscf

mol = pyscf.gto.M(atom='''
                    H 0 0 -0.8
                    H 0 0 0.8''',
                    basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

natm, nbas, nelec, nshells = 2, 2, 2, 2

alpha = 1e-4

max_i = 500
epsilon = 1e-6

delta = 1

i = 0

while not (delta < epsilon or i == max_i):
    P = libscf.RHF(natm, nbas, nelec, nshells, atm, bas, env)
    E = libscf.energy(natm, nbas, nshells, atm, bas, env, P)
    denv = libscf.grad(natm, nbas, nshells, atm, bas, env, P)
    
    env[28:34] = env[28:34] - alpha * denv[28:34]

    delta = np.linalg.norm(denv[28:34])

    # if i % 100 == 0:
    print(i)
    print("E:   ", E)
    # print("grad:", denv[28:34])

    i += 1

if i == max_i:
    print("did not converge")
else:
    print("converged")

print("env:   ", env[28:34])
print("P:   ", P)
print("E:   ", E)
print("grad:", denv[28:34])

S = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')
print("S: ", S)