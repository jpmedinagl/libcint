import libscf
import numpy as np
import pyscf

def printf(A):
    for i in range(len(A)):
      for j in range(len(A[0])):
        if (A[i,j] <= 0):
            print("{:f} ".format(A[i,j]), end='')
        else:
            print(" {:f} ".format(A[i,j]), end='')
      print()


h = 1e-6

# H2

# mol = pyscf.gto.M(atom='''
#                     H 0 0 -0.8
#                     H 0 0 0.8''',
#                     basis='sto-3g')

# for i in range(6):
#     mol._env[-6 + i] += h
#     S1 = mol.intor("int1e_ovlp")

#     mol._env[-6 + i] -= 2.0*h
#     S2 = mol.intor("int1e_ovlp")

#     print(-(S2 - S1)/(2.0*h))

#     mol._env[-6 + i] += h

# atm = np.asarray(mol._atm, dtype=np.int32, order='C')
# bas = np.asarray(mol._bas, dtype=np.int32, order='C')
# env = np.asarray(mol._env, dtype=np.double, order='C')

# natm, nbas, nelec, nshells = 2, 2, 2, 2

# for i in range(6):
#     env[-6 + i] += h
#     S1 = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')

#     env[-6 + i] -= 2.0*h
#     S2 = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')

#     print(-(S2 - S1)/(2.0*h))

#     env[-6 + i] += h

# H2O

# mol = pyscf.gto.M(atom='''
#                     O   -0.0000000   -0.0000000    0.0000000''',
#                     basis='sto-3g')

mol = pyscf.gto.M(atom='''
                    O   -0.0000000   -0.1113512    0.0000000
                    H    0.0000000    0.4454047   -0.7830363
                    H   -0.0000000    0.4454047    0.7830363''',
                    basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

natm, nbas, nshells = 3, 5, 7

env[-1] += h
S1 = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')

env[-1] -= 2.0*h
S2 = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')

printf(-(S2 - S1)/(2.0*h))

env[-1] += h