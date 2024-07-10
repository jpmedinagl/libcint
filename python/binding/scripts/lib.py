import os
import ctypes
import numpy as np

path = os.path.abspath("libscf.so")
print(path)

libgrad = ctypes.CDLL(path)

libgrad.add.argtypes = (ctypes.c_int, ctypes.c_int)
libgrad.add.restype = ctypes.c_int

def subtract(n1, n2) -> float:
    return libgrad.subtract(n1, n2)