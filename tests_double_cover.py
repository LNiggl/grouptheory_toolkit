import numpy as np
from base import groups as g
from base import objects as o
from base import representations as r
from base import testing as t
from base import O_h,gen_actions_O_h                    # import group list of actions used in generation of O_h

np.set_printoptions(precision = 8, suppress = True)

gen_actions = ["Rot0","Rot1","Rot2"]
# v = o.Vector(["a","b","c"])
# ac = o.find_closure(gen_actions,v)
# mt = o.mult_table_from_actions(ac)
# O = o.generate_group(gen_actions,v)

# # Pauli matrices - convention: 0 = x, 1 = y and 2 = z
# sigma = {}
# sigma[0] = np.array([[0,1],[1,0]], dtype=np.complex128)
# sigma[1] = np.array([[0,-1j],[1j,0]], dtype=np.complex128)
# sigma[2] = np.array([[1,0],[0,-1]], dtype=np.complex128)

# # generators of double cover of O -  convention: 0 = x, 1 = y and 2 = z
# A = {}
# A_inv = {}

# A[0]= 1/np.sqrt(2)*(np.eye(2) -1j*sigma[0])
# A[1] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[1])
# A[2] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[2])

# for i in range(3):
#     A_inv[i] = np.linalg.inv(A[i])

# # test against generators of rotation group O - see if double cover holds

# def R_entry(A,j,k):
#     return 1/2*np.trace(sigma[j]@A@sigma[k]@np.linalg.inv(A))
# def R(A):
#     M = np.zeros((3,3),dtype = np.complex128)
#     for j in range(3):
#         for k in range(3):
#             M[j][k] = R_entry(A,j,k)
#     return M

# def group_hom_via_reps(Rep1,Rep2,F):        #F:function relating Rep1 matrices to Rep2 matrices, e.g. Rep1: Double cover of O, Rep2: O, F: R(A); returns dict{groupelement1:groupelement2}
#     group_hom = {}
#     for g1,m1 in Rep1.hom.items():
#         n = R(m1)
#         matches = []
#         for g2,m2 in Rep2.hom.items():
#             if np.allclose(n,m2):
#                 matches.append(g2)
#         assert len(matches) == 1            # any-to-one homomorphism
#         group_hom[g1] = matches[0]
#     return group_hom

# b = o.generate_basis([o.Vector([1,0,0])],O)
# T1 = r.rep_from_action(O,b,"T1")
# T1.check_if_homomorphism()

# #test find_closure

# gen_matrices = {}
# gen_matrices["A0"] = A[0]
# gen_matrices["A1"] = A[1]
# gen_matrices["A2"] = A[2]

# dcover = r.find_closure_from_matrices(gen_matrices)
# mtable = r.mult_table_from_matrices(dcover)
# # print("#:", len(dcover))

# DC_O = g.Group(list(dcover.keys()),mtable)
# G1 = r.Representation(DC_O,dcover,"G1")
# G1.direction_action = "left"
# # provisory assignment of basis
# s1 = o.WeylSpinor(1,0)
# s2 = o.WeylSpinor(0,1)
# G1.basis = [s1,s2]
# G1.check_if_homomorphism()
# assert not G1.is_reducible()

# DC_O_to_O = group_hom_via_reps(G1,T1,R)
# assert g.is_group_homomorphism(DC_O_to_O,DC_O,O)

# a = [o.Scalar(1)]
# A1_dc = r.rep_from_action(DC_O,a,"A1_dc")
# A1_dc.check_if_homomorphism()

# b = o.generate_basis([o.Vector([1,0,0])],O)
# T1_dc = r.rep_from_action(DC_O,b,"T1_dc",group_homomorphism=DC_O_to_O)
# T1_dc.check_if_homomorphism()

# ## test ##
# for A in DC_O.elements:
#     print(T1_dc.hom[A]-T1.hom[DC_O_to_O[A]])


# ############

# T1_x_T1_dc =  r.product_rep(T1_dc,T1_dc)           #no irrep
# T1_x_T1_dc.check_if_homomorphism()

# T2_dc = T1_x_T1_dc.copy("T2_dc")
# r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2_dc)
# T2_dc.check_if_homomorphism()

# E_dc = T1_x_T1_dc.copy("E_dc")
# r.apply_projectors([r.traceless_projector,r.diagonal_projector],E_dc)
# E_dc.check_if_homomorphism()

# A2_dc = r.product_rep(T1_dc,r.product_rep(T1_dc,T1_dc))
# A2_dc.name = "A2_dc"
# r.project_out_irreps(A2_dc, [A1_dc,T1_dc,T2_dc,E_dc])
# A2_dc.check_if_homomorphism()

# G2 = r.product_rep(A2_dc,G1)
# G2.name = "G2"
# G2.check_if_homomorphism()

# H = r.product_rep(G1,E_dc)
# H.name = "H"
# H.check_if_homomorphism()


# DC_O.set_char_table([A1_dc,T1_dc,A2_dc,T2_dc,E_dc,G1,G2,H])
# f = open("D:/Master/Masterarbeit/tests/double_cover/test.txt", "w")
# f.write("Character table of Double Cover of O:\n")
# for irrep in DC_O.char_table.keys():
#     f.write(irrep+ ": ")
#     for  key,val in DC_O.char_table[irrep].items():
#         f.write(str(key)+ ": " +str(val) +". ")
#     f.write("\n")


## group and matrices from weyl spinor action
s_ab = o.WeylSpinor(1/np.sqrt(2),1/np.sqrt(2))
O_from_spinor = o.generate_group(gen_actions,s_ab)
print("haha")
G1_from_spinor = r.rep_from_action(O_from_spinor,[s1,s2],"From_spinor_trafo")
f.write("\ncharacters G1 from spinor trafo")
f.write(str(G1_from_spinor.characters))
f.write("\nclasses of group:\n")
for c,m in O_from_spinor.classes.items():
    f.write(str(c)+ ": "+str(m) +"\n")

#irreps
a = [o.Scalar(1)]
A1_from_spinor = r.rep_from_action(O_from_spinor,a,"A1_from_spinor")
A1_from_spinor.check_if_homomorphism()
print("char table A1:")
print(A1_from_spinor.characters)
b = o.generate_basis([o.Vector([1,0,0])],O_from_spinor)
# o.print_all(b)
O_from_spinor_to_O = group_hom_via_reps(G1_from_spinor,T1,R)
print(O_from_spinor_to_O)
assert g.is_group_homomorphism(O_from_spinor_to_O,O_from_spinor,O)

T1_from_spinor = r.rep_from_action(O_from_spinor,b,"T1_from_spinor",group_homomorphism=O_from_spinor_to_O)
T1_from_spinor.check_if_homomorphism()
print("char table T1:")
print(T1_from_spinor.characters)

T1_x_T1_from_spinor =  r.product_rep(T1_from_spinor,T1_from_spinor)           #no irrep
T1_x_T1_from_spinor.check_if_homomorphism()

T2_from_spinor = T1_x_T1_from_spinor.copy("T2_from_spinor")
r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2_from_spinor)
T2_from_spinor.check_if_homomorphism()

E_from_spinor = T1_x_T1_from_spinor.copy("E_from_spinor")
r.apply_projectors([r.traceless_projector,r.diagonal_projector],E_from_spinor)
E_from_spinor.check_if_homomorphism()

A2_from_spinor = r.product_rep(T1_from_spinor,r.product_rep(T1_from_spinor,T1_from_spinor))
A2_from_spinor.name = "A2_from_spinor"
r.project_out_irreps(A2_from_spinor, [A1_from_spinor,T1_from_spinor,T2_from_spinor,E_from_spinor])
A2_from_spinor.check_if_homomorphism()

## test ##

# for A in DC_O.elements:#["A0","A1","A2"]:
#     print(T1_dc.hom[A]-T1_from_spinor.hom[DC_O_to_O[A]])
#     print(T2_dc.hom[A]-T2_from_spinor.hom[DC_O_to_O[A]])
#     print(E_dc.hom[A]-E_from_spinor.hom[DC_O_to_O[A]])
#     print(A2_dc.hom[A]-A2_from_spinor.hom[DC_O_to_O[A]])

## irreps continued ##

G2_from_spinor = r.product_rep(G1_from_spinor,A2_from_spinor)
G2_from_spinor.name = "G2_from_spinor"
# G2_from_spinor.check_if_homomorphism()


H_from_spinor = r.product_rep(E_from_spinor,G1_from_spinor)
H_from_spinor.name = "H_from_spinor"
# H_from_spinor.check_if_homomorphism()

# test
# for A in DC_O.elements:#["A0","A1","A2"]:
#     print(T1_dc.hom[A]-T1_from_spinor.hom[DC_O_to_O[A]])
#     print(T2_dc.hom[A]-T2_from_spinor.hom[DC_O_to_O[A]])
#     print(E_dc.hom[A]-E_from_spinor.hom[DC_O_to_O[A]])
#     print(A2_dc.hom[A]-A2_from_spinor.hom[DC_O_to_O[A]])
    # print(G2.hom[A]-G2_from_spinor.hom[DC_O_to_O[A]])     ## all equal except those for H (maybe incorrect comparison?)
    # print(H.hom[A]-H_from_spinor.hom[DC_O_to_O[A]])


# O_from_spinor.set_char_table([A1_from_spinor,T1_from_spinor,A2_from_spinor,T2_from_spinor,E_from_spinor,G1_from_spinor,G2_from_spinor,H_from_spinor])
# # f = open("D:/Master/Masterarbeit/tests/double_cover/test.txt", "w")
# f.write("Character table of Double Cover of O:\n")
# for irrep in O_from_spinor.char_table.keys():
#     f.write(irrep+ ": ")
#     for  key,val in O_from_spinor.char_table[irrep].items():
#         f.write(str(key)+ ": " +str(val) +". ")
#     f.write("\n")
