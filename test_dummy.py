import numpy as np
import groups as g
import representations as r
import objects as o

# group O_h

gen_actions = ["Rot0","Rot1","Rot2","Inv"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O_h = o.generate_group(gen_actions,v)

# group O

# gen_actions_O = ["Rot0","Rot1","Rot2"]
# v = o.Vector(["a","b","c"])
# ac_O = o.find_closure(gen_actions_O,v)
# mt_O = o.mult_table_from_actions(ac_O)
# O = o.generate_group(gen_actions_O,v)

b = o.generate_basis([o.Vector([1,0,0])],O_h)

A1p = r.Representation(O_h,r.rep_trivial(O_h),"A1p")
A1p.check_if_homomorphism()

T1m = r.rep_from_action(O_h,b,"T1m")
# T1m = r.Representation(O_h,Rep_T1m,"T1m")
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

## tests Pion and Two_Pion reps and finding their irreps

# op = "InvRot2"

p1 = o.Pion([1,0,0],"+")
p2 = o.Pion([-1,0,0],"+")
# b_p = o.generate_basis([p1],O_h)
# print("P Basis")
# o.print_all(b_p)
# P_Rep = r.rep_from_action(O_h,b_p,"P_Rep")
# print("P_rep ",op, ":")
# print(P_Rep.hom[op])
# r.find_irreps(P_Rep,O_h)

# Two Pion distinguishable
tp1 = o.Two_Pion(p1,p2)
b_tp = o.generate_basis([tp1],O_h) 
print("TP Basis")
o.print_all(b_tp)
TP_Rep = r.rep_from_action(O_h,b_tp,"TP_Rep")
# print("TP_rep ",op, ":")
# print(TP_Rep.hom[op])

s = r.find_irreps(TP_Rep,O_h)
# print("invariant subspaces:")
# for irr in s.keys():
#     print(irr , ":" , s[irr][0])
#     print("EV decomp: " , r.list_nonzero_eigvecs(s[irr][0]))


# print("test distinct_eigvals:")
# m = np.diag([1,2,3])
# n = np.matrix([[1,1,3],[9,4,5],[2,3,5]])
# p = np.diag([1,2,1,3])
# d = np.diag([0,0,2,0,1])
# l = [m,n,p,d]

# print("test make_subspaces_comparable")
# paired_evecs = {}
# temp1 = {2: 3 , 1.0000000002 : 2, 3 : 4 , 4 : np.array([1,2])}
# temp2 = {1 : 2 , 2: 6 , 3 : 5 , 4:3 }  
# print(r.all_evs_different(temp1))
# print(r.all_evs_different(temp2))
# if r.all_evs_different(temp2) and r.all_evs_different(temp1):                   # Evals all different
#             ev_temp1 = sorted(temp1.keys())
#             ev_temp2 = sorted(temp2.keys())
#             if np.allclose(ev_temp1,ev_temp2,1e-5,1e-10):                                                   # same Evals
#                 for i in range(len(ev_temp1)):
#                     paired_evecs[ev_temp1[i]] = [temp1[ev_temp1[i]],temp2[ev_temp2[i]]]
#             print(paired_evecs)

# Two Pion indistinguishable
# tp2 = o.Two_Pion(p1,p2,distinguishable = False)
# tp3 = o.Two_Pion(p2,p1,distinguishable = False)
# # o.print_all([tp2,tp3])
# b_tp2 = o.generate_basis([tp2],O_h) 
# print("TP Basis, indistinguishable")
# o.print_all(b_tp2)
# TP_Rep_indist = r.rep_from_action(O_h,b_tp2,"TP_Rep_indist")
# # print("Two Pion, insdistinguishable ",op, ":")
# # print(TP_Rep.hom[op])
# s1 = r.find_irreps(TP_Rep_indist,O_h)
# print("invariant subspaces:")
# for irr in s1.keys():
#     print(irr , ":" , s1[irr])

# product reps to compare to Andreas' notes

Ep_x_Ep = r.product_rep(Ep,Ep)
Ep_x_T1p = r.product_rep(Ep,T1p)
Ep_x_T2p = r.product_rep(Ep,T2p)
T1p_x_T1p = r.product_rep(T1p,T1p)
T1p_x_T2p = r.product_rep(T1p,T2p)
T2p_x_T2p = r.product_rep(T2p,T2p)

# r.find_irreps(Ep_x_Ep,O_h)
# r.find_irreps(Ep_x_T1p,O_h)
sol1 = r.find_irreps(Ep_x_T2p,O_h)
sol2 = r.find_irreps(T1p_x_T1p,O_h)
# r.find_irreps(T1p_x_T2p,O_h)
# r.find_irreps(T2p_x_T2p,O_h)

## all in accordance with note (assuming no issue with: "subduction" SW_3 -> O and splitting e.g. E -> Ep,Em means Ep = E, Em = E_x_A1m)
print("test: finding same subspace decomp")
P1_T1p = sol1["T1p"][0]
P2_T1p = sol2["T1p"][0]

P1_T2p = sol1["T2p"][0]
P2_T2p = sol2["T2p"][0]
# print(P1_T1p,P2_T1p,P1_T2p,P2_T2p)
# print((abs(P2_T1p) < 1e-2).all())
# print((abs(P2_T2p) < 1e-2).all())           # not zero projectors
rel_T1p = r.make_subspaces_comparable(P1_T1p,P2_T1p,Ep_x_T2p,T1p_x_T1p)
rel_T2p = r.make_subspaces_comparable(P1_T2p,P2_T2p,Ep_x_T2p,T1p_x_T1p)
print(rel_T1p)
print(rel_T2p)
P_T1m_TP = s["T1m"][0]
P_T1m_T1m = r.find_irreps(T1m,O_h)["T1m"][0]
TP_T1m = r.make_subspaces_comparable(P_T1m_TP,P_T1m_T1m,TP_Rep,T1m)


Rho1 = o.Rho([0,0,0])
b_r = o.generate_basis([Rho1],O_h) 
print("Rho Basis")
o.print_all(b_r)
print(len(b_r))
Rho_Rep = r.rep_from_action(O_h,b_r,"Rho_Rep")

# Ve = o.Vector([1,2,3])
# b_v = o.generate_basis([Ve],O_h)
# print("compare vector:" , len(b_v))


rh = r.find_irreps(Rho_Rep,O_h)


# regular representation

# Reg = r.rep_regular(O_h,"Reg")
# Reg.check_if_homomorphism()
# r.find_irreps(Reg,O_h)

## every irrep occurs as many times as its dimension






# current:

 

# next: 

# write project_irreps function


################### MILESTONE ####################

######## double cover



