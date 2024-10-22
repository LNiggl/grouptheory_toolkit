import numpy as np
import groups as g
import objects as o
import representations as r
import numbers
def matrix_equals_LinComb_approach(A,Rep,vec):       # apply transformation two ways: 1. matrix mult of Rep.hom[A]*vec; 2. Via LinearCombination objects -> compare
    #direct trafo of LinearCombination
    weights1 = list(complex(vec[i]) for i in range(len(vec)))
    LinComb1 = o.LinearCombination(Rep.basis,weights1)
    LinComb1.action(A)
    #matrix multiplication
    print("In matrix_equals_LinComb_approach: matrix ", A)
    print(Rep.hom[A])
    if Rep.direction_action == "right":
        vec2 = np.matmul(Rep.hom[A],vec)
    else: 
        vec2 = np.matmul(Rep.hom[A],vec)
    print("In matrix_equals_LinComb_approach: result of matrix mult for:", A)
    print(vec2)
    weights2 = list(complex(vec2[i]) for i in range(len(vec2)))
    LinComb2 = o.LinearCombination(Rep.basis,weights2)
    if not LinComb1.is_equal_to(LinComb2):
        # print("result of LinComb Trafo:")
        l = []
        for i in range(len(LinComb1.lin_comb)):
            l.append(LinComb1.lin_comb[i].obj[0].num)
        # print(l)
        # print("result of matrix mult:")
        # print(weights2)
        print("differences in approaches:")
        for i in range(len(weights2)):
            print(weights2[i]-l[i])
    return LinComb1.is_equal_to(LinComb2)
def matrices_agree_with_LinComb_objects(Rep,vec):
    for A in Rep.group.elements:
        if not matrix_equals_LinComb_approach(A,Rep,vec):
            print("FAIL: in matrices_agree_with_LinComb_objects: no equality for ", A)
            return False
    return True
def subspaces_disjoint(dict_spaces,Rep):                    # dict_spaces: {key: [[vectors of subspace1],[vectors of subspace2], ..]} -> returns true if all spaces are disjoint under all group operations         
    subspaces = []
    for irrep,vecs_list in dict_spaces.items():
        for vecs in vecs_list:
            Orb = []
            for vec in vecs:
                
                weights = list(complex(vec[i]) for i in range(len(vec)))
                LC = o.LinearCombination(Rep.basis,weights)
                result = o.true_orbit(LC,Rep.group)
                Orb.extend(result)
            subspaces.append(Orb)
            print(irrep, ": space appended")
    print("#subspaces for comparison: ", len(subspaces))
    j = 0
    disjoint = True
    while j < len(subspaces):        
        for i in range(j+1,len(subspaces)):
            if o.intersection(subspaces[j],subspaces[i]) == None:
                print(j,i,": disjoint")
            else: 
                print(j,i,": not disjoint")
                disjoint = False
        j += 1
    return disjoint

def reorder_subspaces(subspaces):                   # takes study_irreps output. returns {irrepname : [[vector(s) of 1st subspace],[vector(s) of 2nd subspace], ..]}
    subspaces_reordered = {}
    for irrep,evecs_list in subspaces.items():      # evecs_list: {vector_name: [vectors of such kind]}; vector_name is e.g. "x", and the following ordered list consists of vectors that transform like x
        subspaces_reordered[irrep] = []
        names_vectors = list(evecs_list.keys())
        n_spaces = len(evecs_list[names_vectors[0]])
        for j in range(n_spaces):
            subspace = []                           # append collected vectors of each subspace as a list 
            for component in evecs_list.keys():
                subspace.append(evecs_list[component][j])
            subspaces_reordered[irrep].append(subspace)
    return subspaces_reordered

## functions for labelling subspace vectors -> use after applying reorder subspaces

def T1_labelling(vecs,rep):                     # takes one list of the list of subspaces per each irrep as from reorder_subspaces. returns [LC_objects with .labels x,y,z]
    x = o.LinearCombination(rep.basis,np.array(vecs[0]),label = "x")
    y = o.LinearCombination(rep.basis,np.array(vecs[1]),label = "y")
    z = o.LinearCombination(rep.basis,np.array(vecs[2]),label = "z")
    return [x,y,z]

def T2_labelling(vecs,rep):                     # takes one list of the list of subspaces per each irrep as from reorder_subspaces. returns [LC_objects with .labels tau_{1,2,3}]
    tau_1 = o.LinearCombination(rep.basis,np.array(vecs[0]),label = "tau_1")
    tau_2 = o.LinearCombination(rep.basis,np.array(vecs[1]),label = "tau_2")
    tau_3 = o.LinearCombination(rep.basis,np.array(vecs[2]),label = "tau_3")
    return [tau_1,tau_2,tau_3]
    
## functions applying the tests
def test_trafo_behavior(LCs,actions,outcome_LCs = None):                    # takes [LinearCombination objects] after labelling from one subspace,list of actions. returns dict of results of applied trafos
                                                                            # extend for E: (e_1-e_2) must be added to candidates to check for -> additional linear combs. in outcome_LCs
    trafos = {}
    for A in actions:
        trafos[A] = {}
        for c in LCs:
            temp = c.copy()      
            temp.action(A)
            res = o.match_in_list(temp,LCs)
            if res == None:
                res = o.negative_match_in_list(temp,LCs)
                if res == None: 
                    print("Problem")
                trafos[A][c.label] = o.minus(res.label)
            else:
                trafos[A][c.label] = res.label
    return trafos
################################
# gen_actions = ["Rot0","Rot1","Rot2","Inv"]
# v = o.Vector(["a","b","c"])
# ac = o.find_closure(gen_actions,v)
# mt = o.mult_table_from_actions(ac)
# O_h = o.generate_group(gen_actions,v)

# # irreps of O_h

# a = [o.Scalar(1)]
# A1p = r.rep_from_action(O_h,a,"A1p")
# A1p.check_if_homomorphism()

# b = o.generate_basis([o.Vector([1,0,0])],O_h)
# T1m = r.rep_from_action(O_h,b,"T1m")
# T1m.check_if_homomorphism()

# ps = o.PseudoScalar(1)
# b_ps = o.generate_basis([ps],O_h)
# A1m = r.rep_from_action(O_h,b_ps,"A1m")

# T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
# # T1m_T1m_reps = T1m_x_T1m.hom.copy()
# T1m_x_T1m.check_if_homomorphism()

# # T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
# T1p = T1m_x_T1m.copy("T1p")
# r.apply_projectors([r.antisymmetric_projector],T1p)
# T1p.check_if_homomorphism()

# T2p = T1m_x_T1m.copy("T2p")
# r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2p)
# T2p.check_if_homomorphism()

# T2m = r.product_rep(T2p,A1m)
# T2m.name = "T2m"
# T2m.check_if_homomorphism()

# Ep = T1m_x_T1m.copy("Ep")
# r.apply_projectors([r.traceless_projector,r.diagonal_projector],Ep)
# Ep.check_if_homomorphism()
# # Ep.round_chars()

# Em = r.product_rep(Ep,A1m)
# Em.name = "Em"
# Em.check_if_homomorphism()
# # Em.round_chars()

# A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
# A2m.name = "A2m"
# r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
# A2m.check_if_homomorphism()
# # A2m.round_chars()

# A2p = r.product_rep(A2m,A1m)
# A2p.name = "A2p"
# A2p.check_if_homomorphism()


# list_irreps = [A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep]
# O_h.set_char_table(list_irreps)

# example: (2,1,0) - type 2Pion    

# p7 = o.PseudoScalarField([2,1,0])
# p8 = o.PseudoScalarField([-2,-1,0])
# pi4 = o.TensorProduct(p7,p8)

# b4 = o.generate_basis(pi4,O_h)
# pi210 = r.rep_from_action(O_h,b4,"pi210")
# subspaces = r.study_irreps(pi210,O_h,"../results/twopi210_irreps.txt")


# print(subspaces)

# for irrep_name,eval_evecs in subspaces.items():
#     print("Checking equality of approaches in ", irrep_name)
#     for evecs in eval_evecs.values():
#         for evec in evecs:
#             assert matrices_agree_with_LinComb_objects(pi210,evec)

### check if invariant irred. subspaces are independent of each other


# subspaces_reordered = reorder_subspaces(subspaces)
# # print(subspaces_reordered["T1m"][0])
# subspaces_disjoint(subspaces_reordered,pi210)

###

### test trafo behavior of components of invariant subspaces 
# gen_actions = ["Rot0","Rot1","Rot2","Inv"]
# f = open("../results/trafo_behavior_O_h.txt","w")
# ## T1m
# f.write("T1m\n")
# x = o.LinearCombination(T1m.basis,[1,0,0], label = "x")
# y = o.LinearCombination(T1m.basis,[0,1,0], label = "y")
# z = o.LinearCombination(T1m.basis,[0,0,1], label = "z")
# LC_basis_T1m = [x,y,z]
# f.write("T1m\nBasis: \n")
# for obj in LC_basis_T1m:
#     f.write(obj.label + ": " + obj.name)
#     f.write("\n")
# trafos_T1m = {}
# for A in gen_actions:
#     trafos_T1m[A] = {}
#     for c in LC_basis_T1m:
#         temp = c.copy()       
#         temp.action(A)
#         res = o.match_in_list(temp,LC_basis_T1m)
#         if res == None:
#             res = o.negative_match_in_list(temp,LC_basis_T1m)
#             if res == None: 
#                 print("Problem")
#             trafos_T1m[A][c.label] = o.minus(res.label)
#         else:
#             trafos_T1m[A][c.label] = res.label
# f.write(str(trafos_T1m))
# f.write("\n\n")
# ##T1p test
# b_test = o.generate_basis([o.PseudoVector([1,0,0])],O_h)
# T1p_test= r.rep_from_action(O_h,b_test,"T1p_test")
# T1p_test.check_if_homomorphism()


# # o.print_all(T1p_test.basis)
# x_t = o.LinearCombination(T1p_test.basis,[1,0,0], label = "x")
# y_t = o.LinearCombination(T1p_test.basis,[0,1,0], label = "y")
# z_t = o.LinearCombination(T1p_test.basis,[0,0,1], label = "z")

#+++

# # print(T1p_test.characters)
# # print(O_h.char_table["T1p"])    # same


# ## T1p
# # print("Basis T1p:")
# # o.print_all(T1p.basis)
# # T1p_P = np.array(r.antisymmetric_projector(3))
# # T1p_weights = []
# # print(T1p_P)
# # for i in range(len(T1p_P)):
# #     for j in range(len(T1p_P)):
# #         # print(T1p_P[i][j])
# #         if abs(T1p_P[i][j]) > o.num_tol:
# #     # if any(abs(T1p_P[i][j]) > o.num_tol for j in range(len(T1p_P[i]))):
# #             T1p_weights.append(T1p_P[i])
# #             break
# # # print(T1p_weights)
# # f.write("T1p\n")
# # x = o.LinearCombination(T1p.basis,[1,0,0], label = "x")
# # y = o.LinearCombination(T1p.basis,[0,1,0], label = "y")
# # z = o.LinearCombination(T1p.basis,[0,0,1], label = "z")
# LC_basis_T1p = [x_t,y_t,z_t]
# f.write("T1p\nBasis: \n")
# # for obj in LC_basis_T1p:
# #     f.write(obj.label + ": " + obj.name)
# #     f.write("\n")
# # trafos_T1p = {}
# # for A in gen_actions:
# #     trafos_T1p[A] = {}
# #     for c in LC_basis_T1p:
# #         temp = c.copy()
# #         # print("Start")
# #         # print(temp.name, "label: ", temp.label)        
# #         temp.action(A)
# #         # print("after action " , A)
# #         # print(temp.name, "label: ", temp.label)
# #         res = o.match_in_list(temp,LC_basis_T1p)
# #         if res == None:
# #             res = o.negative_match_in_list(temp,LC_basis_T1p)
# #             if res == None: 
# #                 print("Problem")
# #             trafos_T1p[A][c.label] = o.minus(res.label)
# #         else:
# #             trafos_T1p[A][c.label] = res.label

# trafos_T1p = test_trafo_behavior(LC_basis_T1p,gen_actions)
# f.write(str(trafos_T1p))
# f.write("\n\n")
# ## T2p

# # print("T2p Basis:")
# # o.print_all(T2p.basis)
# T2p_P = np.array(r.symmetric_projector(3)*r.invert_projector(r.diagonal_projector)(3))
# # print(T2p_P)
# ## per reading off
# T2p_weights = [T2p_P[1],T2p_P[2],T2p_P[5]]
# tau_1 = o.LinearCombination(T2p.basis,T2p_weights[0], label = "tau_1")
# tau_2 = o.LinearCombination(T2p.basis,T2p_weights[1], label = "tau_2")
# tau_3 = o.LinearCombination(T2p.basis,T2p_weights[2], label = "tau_3")

# LC_basis_T2p = [tau_1,tau_2,tau_3]
# f.write("T2p\nBasis: \n")
# for obj in LC_basis_T2p:
#     f.write(obj.label + ": " + obj.name)
#     f.write("\n")
# trafos_T2p = {}
# for A in gen_actions:
#     trafos_T2p[A] = {}
#     for c in LC_basis_T2p:
#         temp = c.copy()
#         # print("Start")
#         # print(temp.name, "label: ", temp.label)        
#         temp.action(A)
#         # print("after action " , A)
#         # print(temp.name, "label: ", temp.label)
#         res = o.match_in_list(temp,LC_basis_T2p)
#         if res == None:
#             res = o.negative_match_in_list(temp,LC_basis_T2p)
#             if res == None: 
#                 print("Problem")
#             trafos_T2p[A][c.label] = o.minus(res.label)
#         else:
#             trafos_T2p[A][c.label] = res.label
# f.write(str(trafos_T2p))
# f.write("\n\n")
# ## T2m

# print("T2m Basis:")
# o.print_all(T2m.basis)

# # print(T2m_P)
# ## per reading off
# T2m_weights = [T2p_P[1],T2p_P[2],T2p_P[5]]
# tau_1 = o.LinearCombination(T2m.basis,T2p_weights[0], label = "tau_1")
# tau_2 = o.LinearCombination(T2m.basis,T2p_weights[1], label = "tau_2")
# tau_3 = o.LinearCombination(T2m.basis,T2p_weights[2], label = "tau_3")

# LC_basis_T2m = [tau_1,tau_2,tau_3]
# f.write("T2m\nBasis: \n")
# for obj in LC_basis_T2m:
#     f.write(obj.label + ": " + obj.name)
#     f.write("\n")
# trafos_T2m = {}
# for A in gen_actions:
#     trafos_T2m[A] = {}
#     for c in LC_basis_T2m:
#         temp = c.copy()
#         # print("Start")
#         # print(temp.name, "label: ", temp.label)        
#         temp.action(A)
#         # print("after action " , A)
#         # print(temp.name, "label: ", temp.label)
#         res = o.match_in_list(temp,LC_basis_T2m)
#         if res == None:
#             res = o.negative_match_in_list(temp,LC_basis_T2m)
#             if res == None: 
#                 print("Problem")
#             trafos_T2m[A][c.label] = o.minus(res.label)
#         else:
#             trafos_T2m[A][c.label] = res.label
# f.write(str(trafos_T2m))
# f.write("\n\n")

# ##A2m
# rep = r.product_rep(T1m,r.product_rep(T1m,T1m))
# A2m_P = 0.0 * rep.hom["I"]
# irreps = [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep]
# for R in irreps:
#     A2m_P += r.projector_irrep(rep,R)
# A2m_P = np.matrix(np.eye(len(A2m_P))) - A2m_P
# o.print_all(rep.basis)
# print(A2m_P)

#+++


################################# miscellaneous tests ###############################
# x = o.Vector([1,0,0])
# y = o.Vector([0,1,0])
# z = o.Vector([0,0,1])
# test = o.LinearCombination([x,y,z],[1,2,3])
# test.set_label("test")
# print(test.label)
# test2 = o.LinearCombination([x,y,z],[2,4,7])
# print(test.lin_factor(test2))



