import numpy as np
import pyscf

import time
import utils
import libcscf

SAVE = True

# mol = pyscf.gto.M(atom='''
#                     O   -0.0000000   -0.1113512    0.0000000
#                     H    0.0000000    0.4454047   -0.7830363
#                     H   -0.0000000    0.4454047    0.7830363''',
#                     basis='sto-3g')

# atm = np.asarray(mol._atm, dtype=np.int32, order='C')
# bas = np.asarray(mol._bas, dtype=np.int32, order='C')
# env = np.asarray(mol._env, dtype=np.double, order='C')

# nelec = 10

mol = pyscf.gto.M(atom='''
                    H 0 0 -0.4
                    H 0 0 0.4''',
                    basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

nelec = 2

P = libcscf.RHF(atm, bas, env, nelec)

start = time.time()
E = libcscf.energy(atm, bas, env, P)
end = time.time()
print("E time:", (end - start) * 1000000)

start = time.time()
denv = libcscf.grad(atm, bas, env, P)
end = time.time()
print("dE time:", (end - start) * 1000000)

a, b = utils.split(bas)
print(a, b)

max_i = 100
i = 0

total = 0

while not (i == max_i):
    start = time.perf_counter()
    utils.basis_atm(bas, env)
    
    P = libcscf.RHF(atm, bas, env, nelec)
    E = libcscf.energy(atm, bas, env, P)
    denv = libcscf.grad(atm, bas, env, P)
    end = time.perf_counter()
    
    total += (end - start) * 1000000
    i += 1

print("nparams time")
print(b-a, total/max_i)
