import numpy as np
import pyscf

import utils
import libcscf

import time

SAVE = True

# mol = pyscf.gto.M(atom='''
#                     O   -0.0000000   -0.1113512    0.0000000
#                     H    0.0000000    0.4454047   -0.7830363
#                     H   -0.0000000    0.4454047    0.7830363''',
#                     basis='sto-3g')

# nelec = 10

mol = pyscf.gto.M(atom='''
                    H 0 0 -0.4
                    H 0 0 0.4''',
                    basis='sto-3g')

nelec = 2

atm = np.asarray(mol._atm, dtype=np.int32, order='C')
bas = np.asarray(mol._bas, dtype=np.int32, order='C')
env = np.asarray(mol._env, dtype=np.double, order='C')

a, b = utils.split(bas)
print(a, b)

alpha = 1e-2
max_i = 10000
epsilon = 1e-2
delta = 1
i = 0

h = 1e-6

norm_arr = []
E_arr = []

total = 0

while not (delta < epsilon or i == max_i):
    # NORMALIZE ENV
    start = time.perf_counter()
    utils.basis_atm(bas, env)
    
    P = libcscf.density(atm, bas, env, nelec)
    denergy = libcscf.denergyf(atm, bas, env, P)
    danalytical = libcscf.danalyticalf(atm, bas, env, P)

    E = libcscf.energy(atm, bas, env, P)

    fd = np.zeros(b-a)
    for j in range(a, b):
        env[j] -= h
        P1 = libcscf.density(atm, bas, env, nelec)
        E1 = libcscf.energy(atm, bas, env, P1)
        env[j] += 2.0*h
        P2 = libcscf.density(atm, bas, env, nelec)
        E2 = libcscf.energy(atm, bas, env, P2)

        fd[j-a] = (E2 - E1)/(2.0*h)
        env[j] -= h
    
    # print(fd)
    # print(denergy)
    # print(danalytical)

    env[a:b] = env[a:b] - alpha * denergy
    delta = np.linalg.norm(denergy)
    end = time.perf_counter()


    # exit(1)

    norm_arr.append(delta)
    E_arr.append(E)
    
    print(i)
    print("E: ", E)
    print("delta: ", delta)
    # print("denv:  ", denergy)
    # if (np.linalg.norm(denv - fd) > 1e-4):
    #     print("grad and fd diff i: {0}".format(i))
    #     print("{0}".format(denv)) 
    #     print("{0}".format(fd))
    
    # exit(1)

    total += (end - start) * 1000000

    i += 1

if i == max_i:
    print("did not converge")
else:
    print("converged")

print("Param, time/iter")
print(b-a, total/max_i)

print("env:   ", env[a:b])
# print("P:   ", P)
print("E:   ", E)
# print("grad:", grad)

# S = libcscf.int1e(atm, bas, env, 'ovlp')
# print("S:")
# utils.pmat(S)

if SAVE == True:
    from datetime import datetime

    current_time = datetime.now()
    formatted_time = current_time.strftime("%H%M-%d%m%y")

    f1 = open("files/" + formatted_time + "_norm.txt", "w")
    for norm in norm_arr:
        f1.write(str(norm) + "\n")
    f1.close()

    f2 = open("files/" + formatted_time + "_E.txt", "w")
    for EN in E_arr:
        f2.write(str(EN) + "\n")
    f2.close()