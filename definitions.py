import numpy as np
gamma = {                                                                   # gamma matrices copied from github.com/lehner/gpt/blob/master/lib/gpt/core/gamma.py
    0: np.array(
        [[0, 0, 0, 1j], [0, 0, 1j, 0], [0, -1j, 0, 0], [-1j, 0, 0, 0]],
        dtype=np.complex128,
    ),
    1: np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]], dtype=np.complex128),
    2: np.array(
        [[0, 0, 1j, 0], [0, 0, 0, -1j], [-1j, 0, 0, 0], [0, 1j, 0, 0]],
        dtype=np.complex128,
    ),
    3: np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.complex128),
    4: np.diagflat([1, 1, -1, -1]).astype(dtype=np.complex128),
    11: np.diagflat([1, 1, 1, 1]).astype(dtype=np.complex128)
    }
def comm(A,B):
    return A@B - B@A
def anticomm(A,B):
    return A@B + B@A
S = {(mu,nu) : 0.25 * comm(gamma[mu],gamma[nu]) for mu in range(3) for nu in range(3)}         # generators of spinor representation, S_{mu,nu} = 1/4[gamma_mu,gamma_nu], see e.g. Tong QFT Lecture, ch. 4

## check Clifford algebra of gamma matrices under Euclidian metric:
for mu in range(3):
    for nu in range(3):
        assert (anticomm(gamma[mu],gamma[nu]) - 2*np.eye(4)[mu][nu]*np.eye(4) == 0).all()

## check Lorentz algebra of S_{mu,nu} under Euclidian metric, see Tong ch.4 eq. 4.19
for mu in range(3):
    for nu in range(3):
        for rho in range(3):
            for sigma in range(3):
                assert (comm(S[(mu,nu)],S[(rho,sigma)]) == np.eye(4)[nu][rho]*S[(mu,sigma)] - np.eye(4)[mu][rho]*S[(nu,sigma)]
                         + np.eye(4)[mu][sigma]*S[(nu,rho)] - np.eye(4)[nu][sigma]*S[(mu,rho)]).all()
