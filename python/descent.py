import libscf
import numpy as np
import pyscf

from libscf import pmat

NORMALIZE_GTO = True

def norm(basis_add):
    _env = []
    for b in basis_add:
        angl = b[0]
        b_coeff = np.array(sorted(list(b[1:]), reverse=True))
        es = b_coeff[:,0]
        cs = b_coeff[:,1:]
        nprim, nctr = cs.shape
        # cs = np.einsum('pi,p->pi', cs, pyscf.gto.gto_norm(angl, es))
        if NORMALIZE_GTO:
            cs = pyscf.gto.mole._nomalize_contracted_ao(angl, es, cs)

        _env.append(es)
        _env.append(cs.T.reshape(-1))

    return np.array(_env).flatten()

mol = pyscf.gto.M(atom='''
                    H 0 0 -0.4
                    H 0 0  0.4''',
                    basis='sto-3g')

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

natm, nbas, nelec, nshells = 2, 2, 2, 2

S = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')
print("S: ")
pmat(S)
print(env[28:34])

alpha = 1e-6
max_i = 10000
epsilon = 1e-5
delta = 1
i = 0

norm_arr = []
E_arr = []

while not (delta < epsilon or i == max_i):
    # NORMALIZE ENV
    tmp = np.array([env[28:31], env[31:34]]).T.tolist()
    basis_add = [[0, tmp[0], tmp[1], tmp[2]]]
    env[28:34] = norm(basis_add)

    P = libscf.RHF(natm, nbas, nelec, nshells, atm, bas, env)
    E = libscf.energy(natm, nbas, nshells, atm, bas, env, P)
    denv = libscf.grad(natm, nbas, nshells, atm, bas, env, P)
    
    env[28:34] = env[28:34] - alpha * denv[28:34] #/np.linalg.norm(denv[28:34]))

    delta = np.linalg.norm(denv[28:34])
    norm_arr.append(delta)

    # if i % 100 == 0:
    print(i)
    print("E: ", E)
    E_arr.append(E)

    print("grad: ", delta)

    print("env: ", env[28:34])
    S = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')
    print("S: ")
    pmat(S)

    i += 1

if i == max_i:
    print("did not converge")
else:
    print("converged")

print("env:   ", env[28:34])
# print("P:   ", P)
print("E:   ", E)
print("grad:", denv[28:34])

S = libscf.int1e(natm, nbas, nshells, atm, bas, env, 'ovlp')
print("S:")
pmat(S)

f1 = open("files/norm.txt", "w")
for norm in norm_arr:
    f1.write(str(norm) + "\n")
f1.close()

f2 = open("files/E.txt", "w")
for EN in E_arr:
    f2.write(str(EN) + "\n")
f2.close()