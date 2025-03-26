import numpy as np
from scipy import linalg

num_tol = 1e-8
num_rtol = 1e-5

gen_actions_O_h = ["Rot0","Rot1","Rot2","Inv"]
gen_actions_O = ["Rot0","Rot1","Rot2"]

gamma = {                                                                   # gamma matrices copied from github.com/lehner/gpt/blob/master/lib/gpt/core/gamma.py
    0: np.array(
        [[0, 0, 0, 1j], [0, 0, 1j, 0], [0, -1j, 0, 0], [-1j, 0, 0, 0]],         #x
        dtype=np.complex128,
    ),
    1: np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]], dtype=np.complex128),   #y
    2: np.array(                                                                                    #z
        [[0, 0, 1j, 0], [0, 0, 0, -1j], [-1j, 0, 0, 0], [0, 1j, 0, 0]],
        dtype=np.complex128,
    ),
    3: np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.complex128),     #t
    4: np.diagflat([1, 1, -1, -1]).astype(dtype=np.complex128),                                     #5
    11: np.diagflat([1, 1, 1, 1]).astype(dtype=np.complex128)                                       #id
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

## check Lorentz algebra of S_{mu,nu} under Euclidian metric, see e.g. Tong QFT Lecture ch.4 eq. 4.19
for mu in range(3):
    for nu in range(3):
        for rho in range(3):
            for sigma in range(3):
                assert (comm(S[(mu,nu)],S[(rho,sigma)]) == np.eye(4)[nu][rho]*S[(mu,sigma)] - np.eye(4)[mu][rho]*S[(nu,sigma)]
                         + np.eye(4)[mu][sigma]*S[(nu,rho)] - np.eye(4)[nu][sigma]*S[(mu,rho)]).all()
                
## compute rotations by 90 degrees from generators

R_z = linalg.expm(1/2*2*(-np.pi)/2*S[(0,1)])
R_y = linalg.expm(1/2*2*(-np.pi)/2*S[(0,2)])
R_x = linalg.expm(1/2*2*(-np.pi)/2*S[(1,2)])

# Pauli matrices - convention: 0 = x, 1 = y and 2 = z
sigma = {}
sigma[0] = np.array([[0,1],[1,0]], dtype=np.complex128)
sigma[1] = np.array([[0,-1j],[1j,0]], dtype=np.complex128)
sigma[2] = np.array([[1,0],[0,-1]], dtype=np.complex128)

# generating set of G1 irrep - convention: 0 = x, 1 = y and 2 = z
A = {}
A[0]= 1/np.sqrt(2)*(np.eye(2) -1j*sigma[0])
A[1] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[1])
A[2] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[2])

A_inv = {}
for i in range(3):
    A_inv[i] = np.linalg.inv(A[i])



