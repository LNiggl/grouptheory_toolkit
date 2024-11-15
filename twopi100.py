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

#(1,0,0) - type 2Pion    

p1 = o.PseudoScalarField([1,0,0],explicit_momentum_trafo=True)
p2 = o.PseudoScalarField([-1,0,0],explicit_momentum_trafo=True)
pi1 = o.TensorProduct(p1,p2)
p1.set_name_gpt(p1.namestring_gpt_pion("minus"))
p2.set_name_gpt(p2.namestring_gpt_pion("plus"))
pi1.set_name_gpt(pi1.namestring_gpt_twopion())

b1 = o.generate_basis(pi1,O_h)
pi100 = r.rep_from_action(O_h,b1,"pi100")

## tests ##
f = open("../tests/twopi100_test_trafos.txt","w")
f.write("Test: LinearCombination action agrees with matrix multiplication:   ")
assert t.test_matrices_against_LinearCombinations(pi100)
f.write("Test successful.\n")

subspaces = r.study_irreps(pi100,O_h,"../results/twopi100_irreps.txt")
subspaces_ordered_by_space = t.reorder_subspaces(subspaces)
all_disjoint = t.subspaces_disjoint(subspaces_ordered_by_space,pi100)
subspaces_LC_labelled = t.label_all(subspaces_ordered_by_space,pi100)

f.write("Test: All subspaces disjoint:   ")
assert all_disjoint
f.write("Test successful.\n")

## Invariant subspaces ##

f.write("Basis of the representation:\n")
for i in range(len(pi100.basis)):
    f.write(str(i+1)+":\n")
    f.write(pi100.basis[i].name_gpt)
    f.write("\n")
f.write("Matrices of rotations and parity:\n")
for A in gen_actions:
    f.write(A + ":\n")
    f.write(str(pi100.hom[A]) + "\n")

f.write("\nInvariant irreducible subspaces:\n\n")

for irrep,spaces in subspaces_LC_labelled.items():
    f.write("Irrep " + irrep + ".\n")
    for i in range(len(spaces)):
        f.write("Subspace " + str(i+1) + "/" + str(len(spaces)) + ". Basis:\n")
        for lc in spaces[i]:
            f.write(lc.label + ":\n" + lc.name_gpt + "\n")
        f.write("\n")
        add_candidates = None
        if "E" in irrep:
            difference_vec_of_eps = spaces[i][0].copy()
            difference_vec_of_eps.add(spaces[i][1].negative())
            difference_vec_of_eps.set_label("(eps1 - eps2)")
            add_candidates = [difference_vec_of_eps]

        trafos = t.test_trafo_behavior(spaces[i],gen_actions,add_outcome_candidates=add_candidates)
        f.write("Transformations:\n")
        f.write(str(trafos))
        f.write("\n")
        f.write("Trafo as expected: ")
        b = t.compare_string_to_file(str(trafos),"D:/Master/Masterarbeit/tests/expected_trafos/" + irrep + "_expected.txt")
        f.write(str(b) + "\n")
    f.write("\n")

## create files for operators in gpt convention ##
master_filepath = "D:/Master/Masterarbeit/tests/gpt_folder_structure/"
irrep_folder_prefix = "I1_"
operator_name = "2pi.g5.0.0.1"
# to be put together for a total filepath of: master_filepath/irrep_folder_prefix + <irrep_name> + /<int as basis vector enumerator>/operator_name
for irrep,spaces in subspaces_LC_labelled.items():
    for i in range(len(spaces)):
        if i == 0: 
            folder_name = irrep_folder_prefix + irrep.upper() + "/"
            t.create_operator_files(spaces[i],master_filepath=master_filepath,irrep_folder=folder_name,filename=operator_name)
        else: 
            folder_name = irrep_folder_prefix + irrep.upper() + "/"
            t.create_operator_files(spaces[i],master_filepath=master_filepath,irrep_folder=folder_name,filename=operator_name+".v"+str(i+1))