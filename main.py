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
print("invariant subspaces:")
for irr in s.keys():
    print(irr , ":" , s[irr])

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

# Ep_x_Ep = r.product_rep(Ep,Ep)
# Ep_x_T1p = r.product_rep(Ep,T1p)
# Ep_x_T2p = r.product_rep(Ep,T2p)
# T1p_x_T1p = r.product_rep(T1p,T1p)
# T1p_x_T2p = r.product_rep(T1p,T2p)
# T2p_x_T2p = r.product_rep(T2p,T2p)

# r.find_irreps(Ep_x_Ep,O_h)
# r.find_irreps(Ep_x_T1p,O_h)
# r.find_irreps(Ep_x_T2p,O_h)
# r.find_irreps(T1p_x_T1p,O_h)
# r.find_irreps(T1p_x_T2p,O_h)
# r.find_irreps(T2p_x_T2p,O_h)

## all in accordance with note (assuming no issue with: "subduction" SW_3 -> O and splitting e.g. E -> Ep,Em means Ep = E, Em = E_x_A1m)

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



