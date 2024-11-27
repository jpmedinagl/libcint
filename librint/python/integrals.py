import numpy as np
import pyscf

import libcscf
import utils

# mol = pyscf.gto.M(atom='''
#                     H 0 0 -0.4
#                     H 0 0 0.4''',
#                     basis='sto-3g')

def pmat4(A: np.ndarray):
    for i in range(len(A)):
        for j in range(len(A[0])):
            for k in range(len(A[0][0])):
                for l in range(len(A[0][0][0])):
                    if (A[i,j,k,l] < 0):
                        print("{:f} ".format(A[i,j,k,l]), end='')
                    else:
                        print(" {:f} ".format(A[i,j,k,l]), end='')
                print()
            print()
        print()

mol = pyscf.gto.M(atom='''
                        O   -0.0000000   -0.1113512    0.0000000
                        H    0.0000000    0.4454047   -0.7830363
                        H   -0.0000000    0.4454047    0.7830363''',
                        basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

S = libcscf.int1e(atm, bas, env, 'ovlp')
T = libcscf.int1e(atm, bas, env, 'kin')
V = libcscf.int1e(atm, bas, env, 'nuc')
R = libcscf.int2e(atm, bas, env)

print("overlap")
utils.pmat(S)
# utils.pmat(T)
# utils.pmat(V)
print("repulsion")
pmat4(R)
