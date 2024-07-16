import libscf
import numpy as np
import pyscf

mol = pyscf.gto.M(atom='''
                    O   -0.0000000   -0.1113512    0.0000000
                    H    0.0000000    0.4454047   -0.7830363
                    H   -0.0000000    0.4454047    0.7830363''',
                    basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

natm = 3
nbas = 5
    
nelec = 10
nshells = 7

P = libscf.RHF(natm, nbas, nelec, nshells, atm, bas, env)
print(P)

E = libscf.energy(natm, nbas, nshells, atm, bas, env, P)
print(E)