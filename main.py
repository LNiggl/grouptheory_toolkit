import numpy as np
import groups as g
import representations as r
import objects as o
np.set_printoptions(precision = 6, suppress = True)
# group O_h

gen_actions = ["Rot0","Rot1","Rot2","Inv"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O_h = o.generate_group(gen_actions,v)

# irreps of O_h

A1p = r.Representation(O_h,r.rep_trivial(O_h),"A1p")
A1p.check_if_homomorphism()

b = o.generate_basis([o.Vector([1,0,0])],O_h)
T1m = r.rep_from_action(O_h,b,"T1m")
T1m.check_if_homomorphism()

A1m = r.Representation(O_h,r.rep_determinant(T1m.hom),"A1m")
A1m.check_if_homomorphism()

T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
r.apply_projectors([r.antisymmetric_projector],T1p)
T1p.check_if_homomorphism()

T2p = r.Representation(O_h,T1m_T1m_reps,"T2p")
r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2p)
T2p.check_if_homomorphism()

T2m = r.product_rep(T2p,A1m)
T2m.name = "T2m"
T2m.check_if_homomorphism()

Ep = r.Representation(O_h,T1m_T1m_reps,"Ep")
r.apply_projectors([r.traceless_projector,r.diagonal_projector],Ep)
Ep.check_if_homomorphism()
Ep.round_chars()

Em = r.product_rep(Ep,A1m)
Em.name = "Em"
Em.check_if_homomorphism()
Em.round_chars()

A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
A2m.name = "A2m"
r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
A2m.check_if_homomorphism()
A2m.round_chars()

A2p = r.product_rep(A2m,A1m)
A2p.name = "A2p"
A2p.check_if_homomorphism()
A2p.round_chars()

O_h.set_char_table([A1m,A1p,T1m,T1p,T2m,T2p,Em,Ep,A2m,A2p])
print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": " , O_h.char_table[irrep])

### rho meson (vector particle): rep and irreps

# Rho1 = o.Rho([0,0,0])
# b_r = o.generate_basis([Rho1],O_h) 
# # print("Rho Basis")
# # o.print_all(b_r)
# # print("# ", len(b_r))
# Rho_Rep = r.rep_from_action(O_h,b_r,"Rho_Rep")

# rho_red = r.find_irreps(Rho_Rep,O_h)

# # vectors that transform like (x,y,z): 

# P_Rho_T1m = rho_red["T1m"]
# Rho_vector_components = r.T1_identify_components(P_Rho_T1m,Rho_Rep)
# print("(x,y,z)-like eigenvectors in Rho_Rep:")
# print(Rho_vector_components)

# # sanity check: compare to Vector rep
# v = o.Vector([1,0,0])
# b_v = o.generate_basis([v],O_h) 
# V_Rep = r.rep_from_action(O_h,b_v,"V_Rep")
# print("General Vector Basis")
# o.print_all(b_v)
# print("# ",len(b_v))

# V_red = r.find_irreps(V_Rep,O_h)
# P_A1m_V = V_red["A1m"][0]
# print(P_A1m_V)
# import sys
# sys.exit(0)
# same result as for the rho Vector particle

### Two-Pion systems: reps and irreps

##(1,0,0) type single momenta

p1 = o.Pion([1,0,0],"+")
p2 = o.Pion([-1,0,0],"+")

# Two Pions, distinguishable
tp1 = o.Two_Pion(p1,p2)
b_tp = o.generate_basis([tp1],O_h) 
print("TP1 Basis: (1,0,0)-type single momenta")
o.print_all(b_tp)
print("# ", len(b_tp))
TP_Rep1 = r.rep_from_action(O_h,b_tp,"TP_Rep1")

TP_red1 = r.find_irreps(TP_Rep1,O_h)

# # A1p: 1
# # T1m: 1
# # Ep: 1
# # total dim: 6


# vectors that span A1p,Ep:
P_TP1_A1p = TP_red1["A1p"][0]
print(P_TP1_A1p)
P_TP1_Ep = TP_red1["Ep"][0]
TP1_A1p_vecs = r.list_nonzero_eigvecs(P_TP1_A1p)
TP1_Ep_vecs = r.list_nonzero_eigvecs(P_TP1_Ep)
print("A1p subspace:") 
print(TP1_A1p_vecs)
print("Em subspace:")
print(TP1_Ep_vecs)
# vectors that span T1m:
P_TP1_T1m = TP_red1["T1m"]
TP1_vector_components = r.T1_identify_components(P_TP1_T1m,TP_Rep1)
print("(x,y,z)-like eigenvectors in TP1:")
print(TP1_vector_components)

# ## (1,1,0)-type momenta

p3 = o.Pion([1,1,0],"+")
p4 = o.Pion([-1,-1,0],"+")

# Two Pions, distinguishable
tp2 = o.Two_Pion(p3,p4)
b_tp2 = o.generate_basis([tp2],O_h) 
print("TP2 Basis: (1,1,0)-type single momenta")
o.print_all(b_tp2)
print("# ", len(b_tp2))
TP_Rep2 = r.rep_from_action(O_h,b_tp2,"TP_Rep2")

TP_red2 = r.find_irreps(TP_Rep2,O_h)

# # A1p: 1
# # T1m: 1
# # T2m: 1
# # T2p: 1
# # Ep: 1

## vectors spanning subspaces:

P_TP2_A1p = TP_red2["A1p"][0]
P_TP2_T1m = TP_red2["T1m"]
P_TP2_T2m = TP_red2["T2m"][0]

TP2_A1p_vecs = r.list_nonzero_eigvecs(P_TP2_A1p)
print("A1p subspace:") 
print(TP2_A1p_vecs)

TP2_vector_components = r.T1_identify_components(P_TP2_T1m,TP_Rep2)
print("(x,y,z)-like eigenvectors in TP1:")
print(TP2_vector_components)

#(1,1,1) type single momenta

p5 = o.Pion([1,1,1],"+")
p6 = o.Pion([-1,-1,-1],"+")

# Two Pions, distinguishable
tp3 = o.Two_Pion(p5,p6)
b_tp3 = o.generate_basis([tp3],O_h) 
print("TP3 Basis: (1,1,1)-type single momenta")
o.print_all(b_tp3)
print("# ", len(b_tp3))
TP_Rep3 = r.rep_from_action(O_h,b_tp3,"TP_Rep3")

TP_red3 = r.find_irreps(TP_Rep3,O_h)

P_TP3_A1p = TP_red3["A1p"][0]
P_TP3_T1m = TP_red3["T1m"]

TP3_A1p_vecs = r.list_nonzero_eigvecs(P_TP3_A1p)
print("A1p subspace:") 
print(TP3_A1p_vecs)

TP3_vector_components = r.T1_identify_components(P_TP3_T1m,TP_Rep3)
print("(x,y,z)-like eigenvectors in TP3:")
print(TP3_vector_components)
#(2,1,0) type single momenta

p7 = o.Pion([2,1,0],"+")
p8 = o.Pion([-2,-1,0],"+")

# Two Pions, distinguishable
tp4 = o.Two_Pion(p7,p8)
b_tp4 = o.generate_basis([tp4],O_h) 
print("TP3 Basis: (2,1,0)-type single momenta")
o.print_all(b_tp4)
print("# ", len(b_tp4))
TP_Rep4 = r.rep_from_action(O_h,b_tp4,"TP_Rep3")

TP_red4 = r.find_irreps(TP_Rep4,O_h)

# A1p : 1
# T1m: 2
# T1p : 1
# T2m : 2
# T2p : 1
# Ep : 2
# A2p : 1

P_TP4_A1p = TP_red4["A1p"][0]
P_TP4_T1m = TP_red4["T1m"]

TP4_A1p_vecs = r.list_nonzero_eigvecs(P_TP4_A1p)
print("A1p subspace:") 
print(TP4_A1p_vecs)

TP4_vector_components = r.T1_identify_components(P_TP4_T1m,TP_Rep4)
print("(x,y,z)-like eigenvectors in TP4:")
print(TP4_vector_components)

#(2,1,1) type single momenta

p9 = o.Pion([2,1,1],"+")
p10 = o.Pion([-2,-1,-1],"+")

# Two Pions, distinguishable
tp5 = o.Two_Pion(p9,p10)
b_tp5 = o.generate_basis([tp5],O_h) 
print("TP4 Basis: (2,1,1)-type single momenta")
o.print_all(b_tp5)
print("# ", len(b_tp5))
TP_Rep5 = r.rep_from_action(O_h,b_tp5,"TP_Rep5")

TP_red5 = r.find_irreps(TP_Rep5,O_h)

# A1p : 1
# T1m : 2
# T1p : 1
# T2m : 1
# T2p : 2
# Em : 1
# Ep : 1
# A2m : 1

P_TP5_A1p = TP_red5["A1p"][0]
P_TP5_T1m = TP_red5["T1m"]

TP5_A1p_vecs = r.list_nonzero_eigvecs(P_TP5_A1p)
print("A1p subspace:") 
print(TP5_A1p_vecs)

TP5_vector_components = r.T1_identify_components(P_TP5_T1m,TP_Rep5)
print("(x,y,z)-like eigenvectors in TP5:")
print(TP5_vector_components)

#shared with rho: T1m, T2m, Em, A2m


# #tests

# lc = o.linear_combination(b_tp,[0,1,1,0,1,1])





