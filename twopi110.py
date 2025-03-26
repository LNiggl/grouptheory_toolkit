from base import objects as o
from base import representations as r
from base import tools as t
from base import O_h,gen_actions_O_h                    # import group list of actions used in generation of O_h


#(1,1,0) - type 2Pion    
filepath = "D:/Master/Masterarbeit/results/twopi/data/"
name = "twopi110"
p1 = o.PseudoScalarField([1,1,0],modified_momentum_trafo=True)
p2 = o.PseudoScalarField([-1,-1,0],modified_momentum_trafo=True)
pi1 = o.TensorProduct(p1,p2)
p1.set_name_gpt(p1.namestring_gpt_pion("minus"))
p2.set_name_gpt(p2.namestring_gpt_pion("plus"))
pi1.set_name_gpt(pi1.namestring_gpt_twopion())

b1 = o.generate_basis(pi1,O_h)
pi110 = r.rep_from_action(O_h,b1,"twopi110")

## tests ##
f = open(filepath + "tests/" + name + "_trafos.txt", "w")
f.write("Test: LinearCombination action agrees with matrix multiplication:   ")
assert t.test_matrices_against_LinearCombinations(pi110)
f.write("Test successful.\n")

subspaces = r.study_irreps(pi110,O_h,filepath + "summary_irreps/" + name + ".txt")
subspaces_ordered_by_space = t.reorder_subspaces(subspaces)
t.export_vectors(subspaces_ordered_by_space,filepath + name + "_vecdata", real = True)
all_disjoint = t.subspaces_disjoint(subspaces_ordered_by_space,pi110)
subspaces_LC_labelled = t.label_all(subspaces_ordered_by_space,pi110)

f.write("Test: All subspaces disjoint:   ")
assert all_disjoint
f.write("Test successful.\n")

## Invariant subspaces ##

f.write("Basis of the representation:\n")
for i in range(len(pi110.basis)):
    f.write(str(i+1)+":\n")
    f.write(pi110.basis[i].name_gpt)
    f.write("\n")
f.write("Matrices of rotations and parity:\n")
for A in gen_actions_O_h:
    f.write(A + ":\n")
    f.write(str(pi110.hom[A]) + "\n")

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
        trafos = t.test_trafo_behavior(spaces[i],gen_actions_O_h,add_outcome_candidates=add_candidates)
        f.write("Transformations:\n")
        f.write(str(trafos))
        f.write("\n")
        f.write("Trafo as expected: ")
        b = t.compare_string_to_file(str(trafos),"../results/expected_trafos/" + irrep + "_expected.txt")
        f.write(str(b) + "\n")
    f.write("\n")

## create files for operators in gpt convention ##
master_filepath = filepath + "files_operators_gpt/"
irrep_folder_prefix = "I1_"
operator_name = "2pi.g5.0.1.1"
# to be put together for a total filepath of: master_filepath/irrep_folder_prefix + <irrep_name> + /<int as basis vector enumerator>/operator_name
for irrep,spaces in subspaces_LC_labelled.items():
    for space in spaces:
        folder_name = irrep_folder_prefix + irrep.upper() + "/"
        t.create_operator_files(space,master_filepath=master_filepath,irrep_folder=folder_name,filename=operator_name)