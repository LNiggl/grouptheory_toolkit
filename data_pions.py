import numpy as np
import groups as g
import representations as r
import objects as o

np.set_printoptions(precision = 8, suppress = True)

gen_actions = ["Rot0","Rot1","Rot2","Inv"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O_h = o.generate_group(gen_actions,v)

# irreps of O_h

a = [o.Scalar(1)]
A1p = r.rep_from_action(O_h,a,"A1p")
A1p.check_if_homomorphism()

b = o.generate_basis([o.Vector([1,0,0])],O_h)
T1m = r.rep_from_action(O_h,b,"T1m")
T1m.check_if_homomorphism()

ps = o.PseudoScalar(1)
b_ps = o.generate_basis([ps],O_h)
A1m = r.rep_from_action(O_h,b_ps,"A1m")

T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
r.apply_projectors([r.antisymmetric_projector],T1p)
T1p.check_if_homomorphism()


T2p = T1m_x_T1m.copy("T2p")
r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2p)
T2p.check_if_homomorphism()

T2m = r.product_rep(T2p,A1m)
T2m.name = "T2m"
T2m.check_if_homomorphism()

Ep = T1m_x_T1m.copy("Ep")
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

O_h.set_char_table([A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep])
print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": " , O_h.char_table[irrep])

#(1,0,0) - type 2Pion    

p1 = o.PseudoScalarField([1,0,0])
p2 = o.PseudoScalarField([-1,0,0])
pi1 = o.TensorProduct(p1,p2)

b1 = o.generate_basis(pi1,O_h)
pi100 = r.rep_from_action(O_h,b1,"pi100")
r.study_irreps(pi100,O_h,"../results/twopi100_irreps.txt")

#(1,1,0) - type 2Pion    

p3 = o.PseudoScalarField([1,1,0])
p4 = o.PseudoScalarField([-1,-1,0])
pi2 = o.TensorProduct(p3,p4)

b2 = o.generate_basis(pi2,O_h)
pi110 = r.rep_from_action(O_h,b2,"pi110")
r.study_irreps(pi110,O_h,"../results/twopi110_irreps.txt")

#(1,1,1) - type 2Pion    

p5 = o.PseudoScalarField([1,1,1])
p6 = o.PseudoScalarField([-1,-1,-1])
pi3 = o.TensorProduct(p5,p6)

b3 = o.generate_basis(pi3,O_h)
pi111 = r.rep_from_action(O_h,b3,"pi111")
r.study_irreps(pi111,O_h,"../results/twopi111_irreps.txt")

#(2,1,0) - type 2Pion    

p7 = o.PseudoScalarField([2,1,0])
p8 = o.PseudoScalarField([-2,-1,0])
pi4 = o.TensorProduct(p7,p8)

b4 = o.generate_basis(pi4,O_h)
pi210 = r.rep_from_action(O_h,b4,"pi210")
r.study_irreps(pi210,O_h,"../results/twopi210_irreps.txt")

#(2,1,1) - type 2Pion    

p9 = o.PseudoScalarField([2,1,1])
p10 = o.PseudoScalarField([-2,-1,-1])
pi5 = o.TensorProduct(p9,p10)

b5 = o.generate_basis(pi5,O_h)
pi211 = r.rep_from_action(O_h,b5,"pi210")
r.study_irreps(pi211,O_h,"../results/twopi211_irreps.txt")