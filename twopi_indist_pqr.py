import numpy as np

from base import groups as g
from base import objects as o
from base import representations as r
from base import tools as t
from base import O_h,gen_actions_O_h                    # import group list of actions used in generation of O_h


#(1,0,0) - type 2Pion    
filepath = "D:/Master/Masterarbeit/results/twopi_indist/data/"
name = "twopi_indist_pqr"
p1 = o.PseudoScalarField(["p","q","r"],modified_momentum_trafo=True)
p2 = o.PseudoScalarField(["-p","-q","-r"],modified_momentum_trafo=True)
pi1 = o.TensorProduct(p1,p2,distinguishable = False)
p1.set_name_gpt(p1.namestring_gpt_pion("minus"))
p2.set_name_gpt(p2.namestring_gpt_pion("plus"))
pi1.set_name_gpt(pi1.namestring_gpt_twopion())

b1 = o.generate_basis(pi1,O_h)
pipqr = r.rep_from_action(O_h,b1,"twopi_indist_pqr")

## tests ##
f = open(filepath + "tests/" + name + "_trafos.txt", "w")
f.write("Test: LinearCombination action agrees with matrix multiplication:   ")
assert t.test_matrices_against_LinearCombinations(pipqr)
f.write("Test successful.\n")

subspaces = r.study_irreps(pipqr,O_h,filepath + "summary_irreps/" + name + ".txt")
subspaces_ordered_by_space = t.reorder_subspaces(subspaces)
t.export_vectors(subspaces_ordered_by_space,filepath + name + "_vecdata", real = True) # D:/Master/Masterarbeit/results/test_01/twopi_indist_pqr_vecdata", real = True)
# import sys
# sys.exit()
all_disjoint = t.subspaces_disjoint(subspaces_ordered_by_space,pipqr)
subspaces_LC_labelled = t.label_all(subspaces_ordered_by_space,pipqr)

f.write("Test: All subspaces disjoint:   ")
assert all_disjoint
f.write("Test successful.\n")

## Invariant subspaces ##

f.write("Basis of the representation:\n")
for i in range(len(pipqr.basis)):
    f.write(str(i+1)+":\n")
    f.write(pipqr.basis[i].name_gpt)
    f.write("\n")
# f.write("Matrices of rotations and parity:\n")
# for A in gen_actions_O_h:
#     f.write(A + ":\n")
#     f.write(str(pipqr.hom[A]) + "\n")

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

## all (complex) values are very close to natural square roots (normalization factors). Testing accuracy if these values would be assumed 
nums = [0,1/np.sqrt(2),1/np.sqrt(4),1/np.sqrt(6),1/np.sqrt(8),1/np.sqrt(12),1/np.sqrt(16),1/np.sqrt(24),1/np.sqrt(32),1/np.sqrt(48),1/np.sqrt(96)]
for irrep in subspaces_ordered_by_space:
    for i in range(len(subspaces_ordered_by_space[irrep])):
        d = t.test_largest_deviations_from_set(subspaces_ordered_by_space[irrep][i],nums)
        f.write("For" + irrep + " space " + str(i+1) + ":\nLargest deviation of a any vector component from a value expressed by 1/sqrt(n), n integer: " + str(d[0]) + "; in vector:\n" + str(d[1]) + "\n")
        vecs_modified_values = t.adjust_values(subspaces_ordered_by_space[irrep][i],nums)
        modified_labelled_vecs = t.label(vecs_modified_values,irrep,pipqr)
        add_candidates = None
        if "E" in irrep:
            difference_vec_of_eps = modified_labelled_vecs[0].copy()
            difference_vec_of_eps.add(modified_labelled_vecs[1].negative())
            difference_vec_of_eps.set_label("(eps1 - eps2)")
            add_candidates = [difference_vec_of_eps]
        trafos_mod = t.test_trafo_behavior(modified_labelled_vecs,gen_actions_O_h,add_outcome_candidates=add_candidates)
        b = t.compare_string_to_file(str(trafos_mod),"../results/expected_trafos/" + irrep + "_expected.txt")
        f.write("If vector components are rounded to 1/sqrt(n) values, trafos still as expected: " + str(b)+"\n")
    f.write("\n")

## create files for operators in gpt convention ##
# master_filepath = filepath + "files_operators_gpt/" # D:/Master/Masterarbeit/results/files_operators_gpt/"
# irrep_folder_prefix = "I1_"
# operator_name = "2pi.g5.0.0.1"
# # to be put together for a total filepath of: master_filepath/irrep_folder_prefix + <irrep_name> + /<int as basis vector enumerator>/operator_name
# for irrep,spaces in subspaces_LC_labelled.items():
#     for i in range(len(spaces)):
#         if i == 0: 
#             folder_name = irrep_folder_prefix + irrep.upper() + "/"
#             t.create_operator_files(spaces[i],master_filepath=master_filepath,irrep_folder=folder_name,filename=operator_name)
#         else: 
#             folder_name = irrep_folder_prefix + irrep.upper() + "/"
#             t.create_operator_files(spaces[i],master_filepath=master_filepath,irrep_folder=folder_name,filename=operator_name+".v"+str(i+1))