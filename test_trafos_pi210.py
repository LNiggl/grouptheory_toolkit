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
f = open("../tests/twopi210_test_trafos.txt","w")

p7 = o.PseudoScalarField([2,1,0])
p8 = o.PseudoScalarField([-2,-1,0])
pi4 = o.TensorProduct(p7,p8)
p7.set_name_gpt(p7.namestring_gpt_pion("minus"))
p8.set_name_gpt(p8.namestring_gpt_pion("plus"))
pi4.set_name_gpt(pi4.namestring_gpt_twopion())

b4 = o.generate_basis(pi4,O_h)
pi210 = r.rep_from_action(O_h,b4,"pi210")
weights_pi210 = []

for i in range(len(pi210.basis)):
    basisv = [0 for j in range(len(pi210.basis))]
    basisv[i] = 1
    weights_pi210.append(basisv)

f.write("Test: LinearCombination action agrees with matrix multiplication:   ")
for i in range(len(weights_pi210)):
    assert t.matrices_agree_with_LinComb_objects(pi210,weights_pi210[i])
f.write("Test successful.\n")

subspaces = r.study_irreps(pi210,O_h,"../results/twopi210_irreps.txt")

# # check approach with LinerCombination objects
# for irrep_name,eval_evecs in subspaces.items():
#     print("Checking equality of approaches in ", irrep_name)
#     for evecs in eval_evecs.values():
#         for evec in evecs:
#             assert t.matrices_agree_with_LinComb_objects(pi210,evec)

### check if invariant irred. subspaces are independent of each other

subspaces_ordered_by_space = t.reorder_subspaces(subspaces)
all_disjoint = t.subspaces_disjoint(subspaces_ordered_by_space,pi210)
f.write("Test: All subspaces disjoint:   ")
assert all_disjoint
f.write("Test successful.\n")

f.write("Basis of the representation:\n")
for i in range(len(pi210.basis)):
    f.write(str(i+1)+":\n")
    f.write(pi210.basis[i].name_gpt)
    f.write("\n")
f.write("Matrices of rotations and parity:\n")
for A in gen_actions:
    f.write(A + ":\n")
    f.write(str(pi210.hom[A]) + "\n")

#######################################

f.write("\nInvariant irreducible subspaces:\n\n")

f.write("T1m\n")

LC_T1m_space_1 = t.T1_labelling(subspaces_ordered_by_space["T1m"][0],pi210)
LC_T1m_space_2 = t.T1_labelling(subspaces_ordered_by_space["T1m"][1],pi210)
for lc in LC_T1m_space_1:
    lc.set_name_gpt()
for lc in LC_T1m_space_2:
    lc.set_name_gpt()
T1m_trafos_1 = t.test_trafo_behavior(LC_T1m_space_1,gen_actions)
T1m_trafos_2 = t.test_trafo_behavior(LC_T1m_space_2,gen_actions)

f.write("Subspace 1/2. Basis:\n")
for lc in LC_T1m_space_1:
    f.write(lc.label + ": " + lc.name_gpt + "\n")
f.write("\nTransformations:\n")
f.write(str(T1m_trafos_1))
f.write("\n\n")

f.write("Subspace 2/2. Basis:\n")
for lc in LC_T1m_space_2:
    f.write(lc.label+ ": " + lc.name_gpt + "\n")
f.write("\nTransformations:\n")
f.write(str(T1m_trafos_2))
f.write("\n\n\n")


f.write("T1p\n")

LC_T1p_space = t.T1_labelling(subspaces_ordered_by_space["T1p"][0],pi210)
T1p_trafos = t.test_trafo_behavior(LC_T1p_space,gen_actions)
for lc in LC_T1p_space:
    lc.set_name_gpt()

f.write("Subspace 1/1. Basis:\n")
for lc in LC_T1p_space:
    f.write(lc.label + ": " + lc.name_gpt + "\n")
f.write("Transformations:\n")
f.write(str(T1p_trafos))
f.write("\n\n\n")

f.write("T2m\n\n")

LC_T2m_space_1 = t.T2_labelling(subspaces_ordered_by_space["T2m"][0],pi210)
LC_T2m_space_2 = t.T2_labelling(subspaces_ordered_by_space["T2m"][1],pi210)

T2m_trafos_1 = t.test_trafo_behavior(LC_T2m_space_1,gen_actions)
T2m_trafos_2 = t.test_trafo_behavior(LC_T2m_space_2,gen_actions)

for lc in LC_T2m_space_1:
    lc.set_name_gpt()
for lc in LC_T2m_space_2:
    lc.set_name_gpt()
f.write("Subspace 1/2. Basis:\n")
for lc in LC_T2m_space_1:
    f.write(lc.label + ": " + lc.name_gpt + "\n")
f.write("\nTransformations:\n")
f.write(str(T2m_trafos_1))
f.write("\n\n")

f.write("Subspace 2/2. Basis:\n")
for lc in LC_T2m_space_2:
    f.write(lc.label+ ": " + lc.name_gpt + "\n")
f.write("\nTransformations:\n")
f.write(str(T2m_trafos_2))
f.write("\n\n\n")

f.write("T2p\n")

LC_T2p_space = t.T2_labelling(subspaces_ordered_by_space["T2p"][0],pi210)
T2p_trafos = t.test_trafo_behavior(LC_T2p_space,gen_actions)
for lc in LC_T2p_space:
    lc.set_name_gpt()
f.write("Subspace 1/1. Basis:\n")
for lc in LC_T2p_space:
    f.write(lc.label + ": " + lc.name_gpt + "\n")
f.write("Transformations:\n")
f.write(str(T2p_trafos))
f.write("\n\n\n")

