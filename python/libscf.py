import os
import ctypes
import numpy as np

path = '/u/jpmedina/libcint/python/libgrad.so'
# path = os.path.abspath("libgrad.so")
# print(path)

libc = ctypes.CDLL(path)

libc.RHF.argtypes = (
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_double
)
libc.RHF.restype = ctypes.POINTER(ctypes.c_double)

libc.energy.argtypes =(
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
)
libc.energy.restype = ctypes.c_double

libc.grad.argtypes = (
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
)
libc.grad.restype = ctypes.POINTER(ctypes.c_double)

def RHF(natm: int, nbas: int, nelec: int, nshells: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, imax: int = 200, conv: float = 1e-6) -> np.ndarray:    
    atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    P_c = libc.RHF(natm, nbas, nelec, nshells, atm_ctypes, bas_ctypes, env_ctypes, imax, conv)
    P = np.ctypeslib.as_array(P_c, shape=(nshells, nshells))
    return P

def energy(natm: int, nbas: int, nshells: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, P: np.ndarray) -> float:
    atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    P_ctypes = P.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return libc.energy(natm, nbas, nshells, atm_ctypes, bas_ctypes, env_ctypes, P_ctypes)

def grad(natm: int, nbas: int, nshells: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, P: np.ndarray) -> np.ndarray:
    atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    P_ctypes = P.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    denv_c = libc.grad(natm, nbas, nshells, atm_ctypes, bas_ctypes, env_ctypes, P_ctypes)
    denv = np.ctypeslib.as_array(denv_c, shape=(1, len(env)))
    return denv.flatten()

if __name__ == '__main__':
    print("Hello world")