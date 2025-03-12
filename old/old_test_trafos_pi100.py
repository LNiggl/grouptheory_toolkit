import numpy as np
import groups as g
import objects as o
import representations as r
import testing as t
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
# T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

# T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
T1p = T1m_x_T1m.copy("T1p")
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
# Ep.round_chars()

Em = r.product_rep(Ep,A1m)
Em.name = "Em"
Em.check_if_homomorphism()
# Em.round_chars()

A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
A2m.name = "A2m"
r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
A2m.check_if_homomorphism()
# A2m.round_chars()

A2p = r.product_rep(A2m,A1m)
A2p.name = "A2p"
A2p.check_if_homomorphism()


list_irreps = [A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep]
O_h.set_char_table(list_irreps)

#(2,1,0) - type 2Pion    
f = open("../tests/twopi100_test_trafos.txt","w")

p1 = o.PseudoScalarField([1,0,0])
p2 = o.PseudoScalarField([-1,0,0])
pi4 = o.TensorProduct(p1,p2)
print(pi4.direction_action)

b4 = o.generate_basis(pi4,O_h)
pi100 = r.rep_from_action(O_h,b4,"pi100")
pi100.check_if_homomorphism()
for A in gen_actions:
    f.write(A + ":\n")
    f.write(str(pi100.hom[A]) + "\n")
weights_pi100 = []
for i in range(len(pi100.basis)):
    basisv = [0 for j in range(len(pi100.basis))]
    basisv[i] = 1
    weights_pi100.append(basisv)

for i in range(len(weights_pi100)):
    print(i)
    assert t.matrices_agree_with_LinComb_objects(pi100,weights_pi100[i])
    
subspaces = r.study_irreps(pi100,O_h,"../results/twopi210_irreps.txt")
subspaces_reordered = t.reorder_subspaces(subspaces)
all_disjoint = t.subspaces_disjoint(subspaces_reordered,pi100)
print(subspaces_reordered)
x = np.array(subspaces_reordered["T1m"][0][0])
y = np.array(subspaces_reordered["T1m"][0][1])
z = np.array(subspaces_reordered["T1m"][0][2])

LC_x = o.LinearCombination(pi100.basis,x,label = "x")
LC_y = o.LinearCombination(pi100.basis,y,label = "y")
LC_z = o.LinearCombination(pi100.basis,z,label = "z")

Rot2_x = np.matmul(pi100.hom["Rot1"],x)
LC_Rot2_x = o.LinearCombination(pi100.basis,Rot2_x,label = "Rot2_x")
for coord in [LC_x,LC_y,LC_z]:
    print("factor between ", LC_Rot2_x.label, " and ", coord.label, ": " , LC_Rot2_x.lin_factor(coord))
o.print_all(pi100.basis)
print(pi100.hom["Rot2"])




