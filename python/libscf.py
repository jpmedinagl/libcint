import os
import ctypes
import numpy as np

# path = '/u/jpmedina/libcint/python/libgrad.a'
path = os.path.abspath("libgrad.so")
print(path)

libc = ctypes.CDLL(path)

libc.RHF.argtypes = (
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_double
)
libc.RHF.restype = ctypes.c_double

# libc.grad.argtypes = (
#     ctypes.c_int,
#     ctypes.c_int,
#     ctypes.c_int,
#     ctypes.POINTER(ctypes.c_int),
#     ctypes.POINTER(ctypes.c_int),
#     ctypes.POINTER(ctypes.c_double),
#     ctypes.POINTER(ctypes.c_double)
# )
# libc.grad.restype = ctypes.POINTER(ctypes.c_double)

libc.grad.argtypes =(
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
)
libc.energy.restype = ctypes.c_double

def RHF(natm: int, nbas: int, nelec: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, imax: int = 200, conv: float = 1e-6) -> float:
    atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return libc.RHF(natm, nbas, nelec, atm_ctypes, bas_ctypes, env_ctypes, imax, conv)

# def grad(natm: int, nbas: int, nelec: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, P: np.ndarray) -> np.ndarray:
#     atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
#     bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
#     env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
#     P_ctypes = P.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

#     grad = libc.grad(natm, nbas, atm, bas, env, P)
#     # convert to numpy array
    
#     return grad

def energy(natm: int, nbas: int, nelec: int, atm: np.ndarray, bas: np.ndarray, env: np.ndarray, P: np.ndarray) -> float:
    atm_ctypes = atm.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    bas_ctypes = bas.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    env_ctypes = env.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    P_ctypes = P.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return libc.energy(natm, nbas, atm, bas, env, P)

if __name__ == '__main__':
    print("hello world")
