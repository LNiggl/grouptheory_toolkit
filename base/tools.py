import numpy as np

from . import groups as g
from . import objects as o
from . import representations as r
from .definitions import num_tol


#####   tools for handling the names of group elements  #####    

def remove_I(A):                    #remove excess "I"s to make keys less ambiguous
    if A == "I":
        return A
    B = ""
    for i in range(len(A)):
        if not A[i] == "I":
            B += A[i]
        else: 
            if not i == len(A) - 1:
                if A[i+1] == "n":               # written for the case that the only other operation with "I" in its name is "Inv"
                    B+= A[i]
    return B
def rename_key(dict,old,new):
    dict[new] = dict.pop(old)
def split(O):                   #splits string of actions ABC into array of actions [C,B,A]
    if O == "I":
        return O
    O = remove_I(O)
    P = []
    if "Rot" in O or "Inv" in O:            # works for formalism of O_h element names
        while len(O)>0:
            if O[-1] == "v":
                Op = O[-3:]
                P.append(Op)
                O = O[:-3]
            else:
                Op = O[-4:]
                P.append(Op)
                O = O[:-4]
    else:                                  # works for names with at most one leading letter and one or multiple digits in name
        while len(O)>0:
            n = 1
            while not O[-n].isalpha():
                n +=1
            Op = O[len(O)-n:]
            O = O[:-n]
            P.append(Op)   
    return P



def matrix_equals_LinComb_approach(A,Rep,vec):       # apply transformation two ways: 1. matrix mult of Rep.hom[A]*vec; 2. Via LinearCombination objects. Then compare both
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
            print("FAIL: in matrices_agree_with_LinComb_objects: inequality for ", A)
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
        if not matrices_agree_with_LinComb_objects(Rep,weights[i]):                         # check approach for all basis vectors
            all_good = False
    return all_good
def subspaces_disjoint(dict_spaces,Rep):                    # takes dict_spaces: {key: [[vectors of subspace1],[vectors of subspace2], ..]} -> returns true if all spaces are disjoint under all group operations         
    subspaces = []
    for irrep,vecs_list in dict_spaces.items():
        print(irrep)
        for vecs in vecs_list:
            Orb = []
            for vec in vecs:
                vec = np.array(vec)
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

def export_vectors(subspaces,filename,real = False):                     # takes output of reorder_subspaces; real: only export real part of complex vectors, note of largest imag omission in seq file
                                                            # writes them in a .npy file as 2dim array
    
    all_vectors = []
    irreps_sequence = []
    if real:
        largest_abs_imag = 0
    for irrep in subspaces.keys():
        for i in range(len(subspaces[irrep])):
            for j in range(len(subspaces[irrep][i])):
                temp = []
                for k in range(len(subspaces[irrep][i][j])):
                    if real:
                        x = subspaces[irrep][i][j][k].item()
                        temp.append(x.real)
                        if abs(largest_abs_imag)-abs(subspaces[irrep][i][j][k].item().imag) < 0:
                            largest_abs_imag = subspaces[irrep][i][j][k].item().imag 
                    else:
                        temp.append(subspaces[irrep][i][j][k].item())
                all_vectors.append(temp)
            irreps_sequence.append(irrep)
    all_vectors = np.array(all_vectors)
    if real:
        irreps_sequence.append(largest_abs_imag)
        print("in " + filename + " largest imag: " +  str(largest_abs_imag))
    np.save(filename,all_vectors)
    np.save(filename+"_irreps_seq",irreps_sequence)
def test_largest_deviations_from_set(vecs,nums):
    if type(vecs) == str:
        vecs = np.load(vecs)
    diff_total = [0,np.zeros(np.shape(vecs[0][0]))]
    for v in vecs:   
        diff_vec = [0,np.zeros(np.shape(v[0]))]            
        for c in v:                   
            diff = [1,np.zeros(np.shape(v))]          
            for n in nums:
                if abs(c-n)<diff[0]:
                    diff[0] = abs(c-n)
                    diff[1] = v
                if abs(c+n)<diff[0]:
                    diff[0] = abs(c+n)
                    diff[1] = v
            if diff_vec[0] < diff[0]:
                diff_vec[0] = diff[0]
                diff_vec[1] = diff[1]
        if diff_total[0] < diff_vec[0]:
            diff_total[0] = diff_vec[0]
            diff_total[1] = diff_vec[1]
    return diff_total
def adjust_values(vecs,nums):                           # modifies every value of 2d array vecs to closest value within {+/-nums}
    for i in range(len(vecs)):   
            for j in range(len(vecs[i])):                   
                diff = [1,0]          
                for n in nums:
                    if abs(vecs[i][j]-n)<diff[0]:
                        diff[0] = abs(vecs[i][j]-n)
                        diff[1] = 0+n
                    if abs(vecs[i][j]+n)<diff[0]:
                        diff[0] = abs(vecs[i][j]+n)
                        diff[1] = 0-n                
                vecs[i][j] = diff[1]
    return vecs
def save_group(G,path):
    from pathlib import Path
    Path(path).mkdir(parents = True,exist_ok=True) 
    f = open(path+"/elements.txt", "w")
    # f.write("begin_elements" + str(G.elements)+"end_elements\n")
    f.write(str(G.elements))
    f.close()
    f = open(path+"/mult_table.txt", "w")
    f.write(str(G.mult_table))
    f.close()
    f = open(path+"/inverse_list.txt", "w")
    f.write(str(G.inverse_list))
    f.close()
    f = open(path+"/classes.txt", "w")
    f.write(str(G.classes))
    f.close()
    if hasattr(G,"char_table"):
        f = open(path+"/char_table.txt", "w")
        f.write(str(G.char_table))
        f.close()
def load_group(path,load_char_table = True):
    import ast
    with open(path + "/elements.txt","r") as file:
        elements = ast.literal_eval(file.read())
    with open(path + "/classes.txt","r") as file:
        classes = ast.literal_eval(file.read())
    with open(path + "/mult_table.txt","r") as file:
        mult_table = ast.literal_eval(file.read())
    with open(path + "/inverse_list.txt","r") as file:
        inverse_list = ast.literal_eval(file.read())
    G = g.Group(elements,mult_table,inv_list = inverse_list,class_list = classes)
    if load_char_table:
        with open(path + "/char_table.txt","r") as file:
            char_table = ast.literal_eval(file.read())
        G.set_char_table(char_table,load_dict = True)
    return G
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
def label(vecs,irrep,Rep):
    if "A" in irrep:
        vecs2 = A_labelling(vecs,Rep)
    if "T1" in irrep:
        vecs2 = T1_labelling(vecs,Rep)
    if "T2" in irrep:
        vecs2 = T2_labelling(vecs,Rep)  
    if "E" in irrep:
        vecs2 = E_labelling(vecs,Rep)  
    if hasattr(Rep.basis[0],"name_gpt"):
        for x in vecs2:
            x.set_name_gpt()
    return vecs2
def label_all(dict_Reps,Rep):                   # takes dictionary of ordered irrep subspace vectors (arrays) from reorder_subspaces, and representation. Returns dict of LinearCombination objects with labels
    from copy import deepcopy
    dict_reps = deepcopy(dict_Reps)
    for irrep, spaces in dict_reps.items():
        for i in range(len(spaces)):
            dict_reps[irrep][i] = label(spaces[i],irrep,Rep)
            # if "A" in irrep:
            #    dict_reps[irrep][i] = A_labelling(spaces[i],Rep)
            # if "T1" in irrep:
            #    dict_reps[irrep][i] = T1_labelling(spaces[i],Rep)
            # if "T2" in irrep:
            #    dict_reps[irrep][i] = T2_labelling(spaces[i],Rep)  
            # if "E" in irrep:
            #    dict_reps[irrep][i] = E_labelling(spaces[i],Rep)  
            # if hasattr(Rep.basis[0],"name_gpt"):
            #     for x in dict_reps[irrep][i]:
            #         x.set_name_gpt()
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
                    print("test_trafo_behavior: Problem:")
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
                e = c
                print("first difference in line:", e)
                break
            c += 1
    if e < 0:
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
                e = c
                print("first difference in line:", e)
                break
            c += 1
    if e < 0:
        return True
    else:   
        return False
def compare_string_to_file(string1,file2):                               #returns True if strings are the same
    with open(file2,"r") as f2:
        string2 = f2.read()
        return string1 == string2

def restore_1d_array(M):
    M = np.array(M)
    T = []
    if not hasattr(M,"__len__"):
        return M
    else:
        for i in range(len(M)):
            if hasattr(M[i],"__len__"):
                if len(M[i])>1:
                    return M
                T.append(M[i][0])
            else:
                T.append(M[i])
    T = np.array(T)
    return T

def test_trafos_matmul(dict_vecs,operations,Rep,skip_H = True):                 # Vecs in format as from reorder_subspaces,e.g. "G1p": [[all s1's],[all s2's]]. operations: list of Strings. 
                                                                                # Rep: Representation. skip_H: True -> do not check H (no concrete phase relation implemented yet for H subspace vectors)
    result = {}
    for irrep in dict_vecs.keys():
        result[irrep] = {}
        vecs = dict_vecs[irrep]
        for i in range(len(vecs)):                                                      # loop over multiplicity of inv subspace
            s1 = restore_1d_array(vecs[i][0].copy())                                  # 1 dim. spaces       
            s = []
            s.append(s1)
            if "E" in irrep or "G" in irrep:    # 2 dim. spaces
                s2 = restore_1d_array(vecs[i][1].copy())
                s.append(s2)
            if "T" in irrep:                                        # 3 dim. spaces
                s2 = restore_1d_array(vecs[i][1].copy())
                s3 = restore_1d_array(vecs[i][2].copy())
                s.append(s2)
                s.append(s3)
            if "H" in irrep:                                        # 4 dim. spaces
                if skip_H:
                    continue
                s2 = restore_1d_array(vecs[i][1].copy())
                s3 = restore_1d_array(vecs[i][2].copy())
                s4 = restore_1d_array(vecs[i][3].copy())
                s.append(s2)
                s.append(s3)
                s.append(s4)
            
            result[irrep]["Space {}".format(i+1)] = {}
            for A in operations:
                result[irrep]["Space {}".format(i+1)][A] = {}
                for n in range(len(s)):
                    temp = Rep.hom[A]@s[n]
                    if np.shape(temp) != np.shape(s[n]):
                        temp = temp.T
                    temp = restore_1d_array(temp)
                    c1 = 0
                    c2 = 0
                    c3 = 0
                    c4 = 0   
                    if "A" in irrep:
                        c1 = r.scalar_prod(s1,temp)              
                        lc = c1*s1
                    if "E" in irrep: 
                        c = [0,0,0,0]
                        compare_to = [s1,s2,(s1-s2),(s1+s2)]
                        for l in range(len(compare_to)):
                            temp1 = temp.copy()
                            temp2 = compare_to[l]
                            j = 0
                            while j < len(temp1)-1:        
                                if temp1[j]>num_tol:
                                    break
                                j+= 1
                            k = temp2[j]/temp1[j]                            
                            for u in range(len(temp1)):
                                temp1[u] = k*temp1[u]                    
                            close = True
                            for o in range(len(temp1)):
                                if abs(temp1[o]-temp2[o])>num_tol:
                                    close = False
                            if close:
                                c[l] = k
                        lc = 0
                        for l in range(len(compare_to)):
                            lc += c[l]*compare_to[l]
                    if "G" in irrep:  
                        for k in range(len(temp)): 
                            c1 += np.conj(s1[k])*temp[k]-np.conj(s1[k])*s2[k]
                            c2 += np.conj(s2[k])*temp[k]-np.conj(s2[k])*s1[k]
                        lc =c1*s1+c2*s2
                    if "T" in irrep:
                        for k in range(len(temp)):
                            c1 += np.conj(s1[k])*temp[k]
                            c2 += np.conj(s2[k])*temp[k]                    
                            c3 += np.conj(s3[k])*temp[k]
                        lc =c1*s1+c2*s2+c3*s3
                    if "H" in irrep:
                        for k in range(len(temp)):
                            c1 += np.conj(s1[k])*temp[k]
                            c2 += np.conj(s2[k])*temp[k]                    
                            c3 += np.conj(s3[k])*temp[k]
                            c4 += np.conj(s4[k])*temp[k]
                        lc =c1*s1+c2*s2+c3*s3+c4*s4
                    if not np.allclose(temp-lc,np.zeros(np.shape(lc)),rtol = num_tol,atol = num_tol):
                        print(A)
                    assert np.allclose(temp-lc,np.zeros(np.shape(lc)),rtol = num_tol,atol = num_tol)                              # vector transforms into given lin. comb.
                    assert abs(1-r.scalar_prod(temp,temp))< num_tol
                    if "A" in irrep:
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["a"] = str(c1) + "*a"
                    if "E" in irrep:
                        result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] = ""
                        if abs(c[0]) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] += str(c[0]) + "*eps_1"
                        if abs(c[1]) > num_tol:
                            if abs(c[0]) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] += "+" + str(c[1]) + "*eps_2"
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] += str(c[1]) + "*eps_2"
                        if abs(c[2]) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] += str(c[2]) + "*(eps_1-eps_2)"
                        if abs(c[3]) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["eps_{}".format(n+1)] += str(c[3]) + "*(eps_1+eps_2)"
                    if "G1" in irrep:
                        result[irrep]["Space {}".format(i+1)][A]["s_{}".format(n+1)] = ""
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["s_{}".format(n+1)] += str(c1) + "*s_1"
                        if abs(c2) > num_tol:
                            if abs(c1) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["s_{}".format(n+1)] += "+" + str(c2) + "*s_2"
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["s_{}".format(n+1)] += str(c2) + "*s_2"
                    if "G2" in irrep:
                        result[irrep]["Space {}".format(i+1)][A]["t_{}".format(n+1)] = ""
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["t_{}".format(n+1)] += str(c1) + "*t_1"
                        if abs(c2) > num_tol:
                            if abs(c1) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["t_{}".format(n+1)] += "+" + str(c2) + "*t_2"
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["t_{}".format(n+1)] += str(c2) + "*t_2"
                    if "T1" in irrep:
                        v = ["x","y","z"]
                        result[irrep]["Space {}".format(i+1)][A][v[n]] = ""
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A][v[n]] += str(c1) + "*x"
                        if abs(c2) > num_tol:
                            if abs(c1) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A][v[n]] += "+" + str(c2) + "*y"
                            else:
                                result[irrep]["Space {}".format(i+1)][A][v[n]] += str(c2) + "*y"
                        if abs(c3) > num_tol:
                            if abs(c1) > num_tol or abs(c2) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A][v[n]] += "+" +str(c3) + "*z"                            
                            else:
                                result[irrep]["Space {}".format(i+1)][A][v[n]] +=  str(c3) + "*z"
                    if "T2" in irrep:
                        result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] = ""
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] += str(c1) + "*tau_1"
                        if abs(c2) > num_tol:
                            if abs(c1) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] += "+" + str(c2) + "*tau_2"
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] += str(c2) + "*tau_2"
                        if abs(c3) > num_tol:
                            if abs(c1) > num_tol or abs(c2) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] += "+" +str(c3) + "*tau_3"                            
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["tau_{}".format(n+1)] +=  str(c3) + "*tau_3"
                    if "H" in irrep:
                        result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] = ""
                        if abs(c1) > num_tol:
                            result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += str(c1) + "*h_1"
                        if abs(c2) > num_tol:
                            if abs(c1) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += "+" + str(c2) + "*h_2"
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += str(c2) + "*h_2"
                        if abs(c3) > num_tol:
                            if abs(c1) > num_tol or abs(c2) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += "+" +str(c3) + "*h_3"                            
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] +=  str(c3) + "*h_3"
                        if abs(c4) > num_tol:
                            if abs(c1) > num_tol or abs(c2) > num_tol or abs(c3) > num_tol:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += "+" +str(c4) + "*h_4"                            
                            else:
                                result[irrep]["Space {}".format(i+1)][A]["h_{}".format(n+1)] += str(c4) + "*h_4"
    return result

    

            






