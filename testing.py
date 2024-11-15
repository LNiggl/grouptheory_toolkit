import numpy as np
import groups as g
import objects as o
import representations as r
def matrix_equals_LinComb_approach(A,Rep,vec):       # apply transformation two ways: 1. matrix mult of Rep.hom[A]*vec; 2. Via LinearCombination objects -> compare
    #direct trafo of LinearCombination
    weights1 = list(complex(vec[i]) for i in range(len(vec)))    
    LinComb1 = o.LinearCombination(Rep.basis,weights1)
    LinComb1.action(A)
    if Rep.direction_action == "right":
        vec2 = np.matmul(Rep.hom[A],vec)
    else: 
        vec2 = np.matmul(Rep.hom[A],vec)
    weights2 = list(complex(vec2[i]) for i in range(len(vec2)))
    LinComb2 = o.LinearCombination(Rep.basis,weights2)
    if not LinComb1.is_equal_to(LinComb2):
        l = []
        for i in range(len(LinComb1.lin_comb)):
            l.append(LinComb1.lin_comb[i].obj[0].num)        
        print("WARNING: in matrix_equals_LinComb_approach: differences in approaches:")
        for i in range(len(weights2)):
            print(weights2[i]-l[i])
    return LinComb1.is_equal_to(LinComb2)
def matrices_agree_with_LinComb_objects(Rep,vec):
    for A in Rep.group.elements:
        if not matrix_equals_LinComb_approach(A,Rep,vec):
            print("FAIL: in matrices_agree_with_LinComb_objects: no equality for ", A)
            return False
    return True
def test_matrices_against_LinearCombinations(Rep):
    weights = []
    for i in range(len(Rep.basis)):
        basisv = [0 for j in range(len(Rep.basis))]
        basisv[i] = 1
        weights.append(basisv)
    all_good = True
    for i in range(len(weights)):
        if not matrices_agree_with_LinComb_objects(Rep,weights[i]):  # check approach for all basis vectors
            all_good = False
    return all_good
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
    print("Test if subspaces are disjoint.\n#subspaces for comparison: ", len(subspaces))
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
    print("Test complete.")
    return disjoint

def reorder_subspaces(subspaces):                   # takes study_irreps output. Returns {irrepname : [[vector(s) of 1st subspace],[vector(s) of 2nd subspace],.. (each in order as r.<..>_identify_components functions yield, e.g. [x,y,z])]}
    subspaces_reordered = {}
    for irrep,evecs_list in subspaces.items():      # evecs_list structure: {vector_name: [vectors of such kind]}; vector_name is e.g. "x", and the following ordered list consists of vectors that transform like x
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
def A_labelling(vec,rep):
    a = o.LinearCombination(rep.basis,np.array(vec[0]),label = "a")
    return [a]
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

def E_labelling(vecs,rep):                     # takes one list of the list of subspaces per each irrep as from reorder_subspaces. returns [LC_objects with .labels eps_{1,2}]
    eps_1 = o.LinearCombination(rep.basis,np.array(vecs[0]),label = "eps_1")
    eps_2 = o.LinearCombination(rep.basis,np.array(vecs[1]),label = "eps_2")
    return [eps_1,eps_2]    
def label_all(dict_reps,Rep):                   # takes dictionary of ordered irrep subspace vectors (arrays) from reorder_subspaces, and representation. Returns dict of LinearCombination objects with labels
    for irrep, spaces in dict_reps.items():
        for i in range(len(spaces)):
            if "A" in irrep:
               dict_reps[irrep][i] = A_labelling(spaces[i],Rep)
            if "T1" in irrep:
               dict_reps[irrep][i] = T1_labelling(spaces[i],Rep)
            if "T2" in irrep:
               dict_reps[irrep][i] = T2_labelling(spaces[i],Rep)  
            if "E" in irrep:
               dict_reps[irrep][i] = E_labelling(spaces[i],Rep)  
            if hasattr(Rep.basis[0],"name_gpt"):
                for x in dict_reps[irrep][i]:
                    x.set_name_gpt()
    return dict_reps
## functions applying the tests
def test_trafo_behavior(LCs,actions,add_outcome_candidates = None):                     # takes [LinearCombination objects] after labelling from one subspace,list of actions. returns dict of results of applied trafos
    outcome_LCs = []
    for x in LCs:
        temp = x.copy()
        outcome_LCs.append(temp)
    if hasattr(add_outcome_candidates,"__iter__"):
        for x in add_outcome_candidates:
            temp = x.copy()
            outcome_LCs.append(temp)                                                                                      # for E: (e_1-e_2) must be added to candidates to check for
    trafos = {}
    for A in actions:
        trafos[A] = {}
        for c in LCs:
            temp = c.copy()      
            temp.action(A)
            res = o.match_in_list(temp,outcome_LCs)
            if res == None:
                res = o.negative_match_in_list(temp,outcome_LCs)
                if res == None: 
                    print("Problem")
                trafos[A][c.label] = o.minus(res.label)
            else:
                trafos[A][c.label] = res.label
    return trafos

def create_operator_files(vecs,master_filepath,irrep_folder,filename):
    for v in vecs:
        if "A1" in irrep_folder or "A2" in irrep_folder:
            filepath = master_filepath + irrep_folder + "0/"
        if "T1" in irrep_folder:
            if v.label == "x":
                filepath = master_filepath + irrep_folder + "0/"
            if v.label == "y":
                filepath = master_filepath + irrep_folder + "1/"
            if v.label == "z":
                filepath = master_filepath + irrep_folder + "2/"
        if "T2" in irrep_folder:
            if v.label == "tau_1":
                filepath = master_filepath + irrep_folder + "0/"
            if v.label == "tau_2":
                filepath = master_filepath + irrep_folder + "1/"
            if v.label == "tau_3":
                filepath = master_filepath + irrep_folder + "2/"
        if "E" in irrep_folder:
            if v.label == "eps_1":
                filepath = master_filepath + irrep_folder + "0/"
            if v.label == "eps_2":
                filepath = master_filepath + irrep_folder + "1/"
        import os
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        f = open(filepath+filename + ".txt","w")
        f.write(v.name_gpt)
        print("File written: " + filepath+filename + ".txt")
        f.close()
def compare_strings(file1,file2):                               #returns True if strings are the same
    from difflib import Differ
    with open(file1) as f1, open(file2) as f2:
        differ = Differ()
        diff = differ.compare(f1.readlines(),f2.readlines())
        c = 0               #counts lines   
        e = -1               #marker for first error
        for line in diff:
            if not line.startswith(" "):
                # print(line)
                e = c
                print("first difference in line:", e)
                break
            c += 1
    if e < 0:
        # print("no difference; linecount:" , c)
        return True
    else:   
        return False
    
def compare_string_to_file(string1,file2):                               #returns True if strings are the same
    from difflib import Differ
    with open(file2) as f2:
        differ = Differ()
        diff = differ.compare(string1.readlines(),f2.readlines())
        c = 0               #counts lines   
        e = -1               #marker for first error
        for line in diff:
            if not line.startswith(" "):
                # print(line)
                e = c
                print("first difference in line:", e)
                break
            c += 1
    if e < 0:
        # print("no difference; linecount:" , c)
        return True
    else:   
        return False
def compare_string_to_file(string1,file2):                               #returns True if strings are the same
    with open(file2,"r") as f2:
        string2 = f2.read()
        return string1 == string2




