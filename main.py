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

# Two-Pion system: rep and irreps

p1 = o.Pion([1,0,0],"+")
p2 = o.Pion([-1,0,0],"+")

# Two Pions, distinguishable
tp1 = o.Two_Pion(p1,p2)
b_tp = o.generate_basis([tp1],O_h) 
print("TP Basis")
o.print_all(b_tp)
print("# ", len(b_tp))
TP_Rep = r.rep_from_action(O_h,b_tp,"TP_Rep")

TP_red = r.find_irreps(TP_Rep,O_h)

# rho meson: rep and irreps

Rho1 = o.Rho([0,0,0])
b_r = o.generate_basis([Rho1],O_h) 
print("Rho Basis")
o.print_all(b_r)
print("# ", len(b_r))
Rho_Rep = r.rep_from_action(O_h,b_r,"Rho_Rep")

rho_red = r.find_irreps(Rho_Rep,O_h)

# sanity check: compare to general Vector rep
v = o.Vector([1,2,3])
b_v = o.generate_basis([v],O_h) 
print("General Vector Basis")
o.print_all(b_v)
print("# ",len(b_v))
V_Rep = r.rep_from_action(O_h,b_r,"V_Rep")

V_red = r.find_irreps(V_Rep,O_h)
# same result as for Vector particle 

#only overlap in irreps of Two_Pion and Rho Meson: T1m

P_TP_T1m = TP_red["T1m"][0]
P_Rho_T1m = rho_red["T1m"][0]
r.make_subspaces_comparable(P_TP_T1m,P_Rho_T1m,TP_Rep,Rho_Rep)
