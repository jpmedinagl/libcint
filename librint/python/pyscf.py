import numpy as np

import libcscf

def int1e(mol, typei: str = 'ovlp', coord: str = 'cart') -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.int2e(atm, bas, env, typei, coord)


def int2e(mol, coord: str = 'cart') -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.int2e(atm, bas, env, coord)


def density(mol, nelec: int, imax: int = 200, conv: float = 1e-6) -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.density(atm, bas, env, nelec, imax, conv)


def energy(mol, P: np.ndarray) -> float:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.energy_c(atm, bas, env, P)


def scf(mol, nelec: int, imax: int = 200, conv: float = 1e-6) -> float:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.scf(atm, bas, env, nelec, imax, conv)


def grad(mol, P: np.ndarray) -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.grad(atm, bas, env, P)


def dSf(mol, P: np.ndarray) -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.dSf(atm, bas, env, P)


def dHcoref(mol, P: np.ndarray) -> np.ndarray:    
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.dHcoref(atm, bas, env, P)


def dRf(mol, P: np.ndarray) -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.dRf(atm, bas, env, P)


def danalyticalf(mol, P: np.ndarray) -> np.ndarray:
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.danalyticalf(atm, bas, env, P)


def denergyf(mol, P: np.ndarray) -> np.ndarray:    
    atm = np.asarray(mol._atm, dtype=np.int32, order='C')
    bas = np.asarray(mol._bas, dtype=np.int32, order='C')
    env = np.asarray(mol._env, dtype=np.double, order='C')

    return libcscf.denergyf(atm, bas, env, P)


if __name__ == '__main__':
    print("Hello world")