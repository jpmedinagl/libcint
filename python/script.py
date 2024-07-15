import libscf
import numpy as np

natm = 1
nbas = 3
    
nelec = 6
nshells = 5

atm = np.array([[6, 20, 1, 23, 0, 0]], dtype=np.int32)
bas = np.array([[0, 0, 3, 1, 0, 24, 27, 0],
                [0, 0, 3, 1, 0, 30, 33, 0],
                [0, 1, 3, 1, 0, 36, 39, 0]], dtype=np.int32)
env = np.array([
0.,         0.,         0.,         0.,         0.,         0.,
0.,         0.,         0.,         0.,         0.,         0.,
0.,         0.,         0.,         0.,         0.,         0.,
0.,         0.,         0.,         0.,         0.,         0.,
71.616837,  13.045096,  3.5305122,  9.59895231, 9.28368844, 2.89332026,
2.9412494,  0.6834831,  0.2222899, -0.56724617, 0.75873757, 0.57263053,
2.9412494,  0.6834831,  0.2222899,  1.75202641, 1.10172125, 0.17453114])

P = libscf.RHF(natm, nbas, nelec, nshells, atm, bas, env)
print(P)

E = libscf.energy(natm, nbas, nshells, atm, bas, env, P)
print(E)
