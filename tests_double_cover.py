import numpy as np
import groups as g
import representations as r
import objects as o

np.set_printoptions(precision = 8, suppress = True)

gen_actions = ["Rot0","Rot1","Rot2"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O = o.generate_group(gen_actions,v)

# Pauli matrices - convention: 0 = x, 1 = y and 2 = z
sigma = {}
sigma[0] = np.array([[0,1],[1,0]], dtype=np.complex128)
sigma[1] = np.array([[0,-1j],[1j,0]], dtype=np.complex128)
sigma[2] = np.array([[1,0],[0,-1]], dtype=np.complex128)

# generators of double cover of O -  convention: 0 = x, 1 = y and 2 = z
A = {}
A_inv = {}

A[0]= 1/np.sqrt(2)*(np.eye(2) -1j*sigma[0])
A[1] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[1])
A[2] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[2])

for i in range(3):
    A_inv[i] = np.linalg.inv(A[i])

# test against generators of rotation group O - see if double cover holds

def R_entry(A,j,k):
    return 1/2*np.trace(sigma[j]@A@sigma[k]@np.linalg.inv(A))
def R(A):
    M = np.zeros((3,3),dtype = np.complex128)
    for j in range(3):
        for k in range(3):
            M[j][k] = R_entry(A,j,k)
    return M

def group_hom_via_reps(Rep1,Rep2,F):        #F:function relating Rep1 matrices to Rep2 matrices, e.g. Rep1:Double cover, Rep2: O, F: R(A); returns dict{groupelement1:groupelement2}
    group_hom = {}
    for g1,m1 in Rep1.hom.items():
        n = R(m1)
        matches = []
        for g2,m2 in Rep2.hom.items():
            if np.allclose(n,m2):
                matches.append(g2)
        assert len(matches) == 1
        group_hom[g1] = matches[0]
    return group_hom

b = o.generate_basis([o.Vector([1,0,0])],O)
T1 = r.rep_from_action(O,b,"T1")
T1.check_if_homomorphism()


# print(T1.hom["Rot0"]-R(A[0]))
# print(T1.hom["Rot1"]-R(A[1]))
# print(T1.hom["Rot2"]-R(A[2]))

# print(T1.hom["Rot0"]-R(-A[0]))
# print(T1.hom["Rot1"]-R(-A[1]))
# print(T1.hom["Rot2"]-R(-A[2]))

# R(A_i) = R(-A_i) = Rot_i, thus 2-to-1 homomorphism as expected

#test find_closure

gen_matrices = {}
gen_matrices["A0"] = A[0]
gen_matrices["A1"] = A[1]
gen_matrices["A2"] = A[2]

dcover = r.find_closure_from_matrices(gen_matrices)
mtable = r.mult_table_from_matrices(dcover)
# print("#:", len(dcover))

DC_O = g.Group(list(dcover.keys()),mtable)
G1 = r.Representation(DC_O,dcover,"G1")
# provisory assignment of basis
s1 = o.L_Spinor(1,0)
s2 = o.L_Spinor(0,1)
G1.basis = [s1,s2]
G1.check_if_homomorphism()
assert not G1.is_reducible()

DC_O_to_O = group_hom_via_reps(G1,T1,R)
assert g.is_group_homomorphism(DC_O_to_O,DC_O,O)

a = [o.Scalar(1)]
A1_dc = r.rep_from_action(DC_O,a,"A1_dc")
A1_dc.check_if_homomorphism()

b = o.generate_basis([o.Vector([1,0,0])],O)
T1_dc = r.rep_from_action(DC_O,b,"T1_dc",group_homomorphism=DC_O_to_O)
T1_dc.check_if_homomorphism()


T1_x_T1_dc =  r.product_rep(T1_dc,T1_dc)           #no irrep
T1_x_T1_dc.check_if_homomorphism()

T2_dc = T1_x_T1_dc.copy("T2_dc")
r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2_dc)
T2_dc.check_if_homomorphism()

E_dc = T1_x_T1_dc.copy("E_dc")
r.apply_projectors([r.traceless_projector,r.diagonal_projector],E_dc)
E_dc.check_if_homomorphism()

A2_dc = r.product_rep(T1_dc,r.product_rep(T1_dc,T1_dc))
A2_dc.name = "A2_dc"
r.project_out_irreps(A2_dc, [A1_dc,T1_dc,T2_dc,E_dc])
A2_dc.check_if_homomorphism()

G2 = r.product_rep(A2_dc,G1)
G2.name = "G2"
G2.check_if_homomorphism()

H = r.product_rep(G1,E_dc)
H.name = "H"
H.check_if_homomorphism()


DC_O.set_char_table([A1_dc,T1_dc,A2_dc,T2_dc,E_dc,G1,G2,H])
print("Character table of Double Cover of O:")
for irrep in DC_O.char_table.keys():
    print(irrep, ": " , DC_O.char_table[irrep])

s1 = o.L_Spinor(1,0)
s2 = o.L_Spinor(0,1)
G1_from_actions = r.rep_from_action(DC_O,[s1,s2],"G1_again")
# G1_from_actions.check_if_homomorphism()

A = "A0"
print("G1 matrix for", A)
print(G1.hom[A])
print("again but from actions")
print(G1_from_actions.hom[A])

print("characters of G1 as from actions")
print(G1_from_actions.characters)
