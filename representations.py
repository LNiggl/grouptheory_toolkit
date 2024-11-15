import numpy as np
import groups as g
import objects as o
import numbers
class Representation(g.Group):
    # hom = {}                           #dict for homomorphism: G -> GL(n)
    # characters = {}
    def __init__(self,group,dict_matrix,name):             #initialize with dict: g -> M(g)
        self.group = group                                 # UGLY structure, but atm needed to make copy() work. RETHINK: inheritance and Group and Representation constructors
        self.elements = group.elements
        self.classes = group.classes
        self.mult_table = group.mult_table
        self.inverse_list = group.inverse_list

        self.name = name
        self.hom = dict_matrix.copy()
        self.characters = {}
        for c in self.classes:
            self.characters[c] = np.trace(self.hom[c])
    def copy(self,name):
        hom = self.hom.copy()
        new_rep = Representation(self.group,hom,name)
        if hasattr(self,"basis"):
            new_rep.basis = self.basis
        if hasattr(self,"direction_action"):
            new_rep.direction_action = self.direction_action
        return new_rep
    def update(self):
        for c in self.classes:
            self.characters[c] = np.trace(self.hom[c])
    def check_if_homomorphism(self):
        for g in self.elements:
            for h in self.elements:
                if self.direction_action == "right":
                    eps = np.linalg.norm(np.matmul(self.hom[g].T,self.hom[h].T) - self.hom[self.mult_table[g][h]].T)
                elif self.direction_action == "left":
                    eps = np.linalg.norm(np.matmul(self.hom[g],self.hom[h]) - self.hom[self.mult_table[g][h]])
                else:                       # cases with neither right nor left action -> must be combination, which does not form homomorphism                
                    assert False
                assert eps < 1e-13
                
    def is_reducible(self):
        return g.inner_product(self.characters,self.characters,self.group) - 1 > 1e-5
    def round_chars(self):
        for c,chi in self.characters.items():
            if chi - round(chi) < 1e-12:
                self.characters[c] = round(chi)

def matrix_from_action(A,basis,group_homomorphism = None):
    if not group_homomorphism == None:
        print(A,"->",group_homomorphism[A])
        A = group_homomorphism[A]
    d = len(basis)
    M = np.zeros((d,d),dtype = complex)
    acted_basis = []
    for b in basis:
        x = b.copy()
        x.action(A)
        acted_basis.append(x)
    for i in range(d):
        for j in range(d):
            if hasattr(basis[j],"lin_factor"):
                k = basis[j].lin_factor(basis[i])
                if k == None:
                    M[i][j] = 0
                else:
                    M[i][j] = k
            else:
                if basis[j].is_equal_to(acted_basis[i]):                # "if e_i gets transformed to e_j via action A" -> ge_i = M(A)_{ji} e_j == +/- e_j
                    M[j][i] = 1                                         # --> M(A)_{ji} == +/- 1
                if basis[j].is_negative_of(acted_basis[i]):
                    M[j][i] = -1
    # CAREFUL: with right action: M is technically defined as the transpose of the above, but application of trafo must be by right multiplication of that matrix M.T,                          
    # therefore again by left multiplication of (M.T).T -> as implemented, the difference left/right is NOT respected in the matrices, but only in check_if_homomorphism()
    return M

def hom_from_action(group,basis,group_homomorphism = None):                                 # returns the dict for Rep.hom; 
                                                                                            # group_homomorphism: dictionary, e.g. for mapping double cover to single cover actions
    return {g : matrix_from_action(g,basis,group_homomorphism) for g in group.elements}
def rep_from_action(group,basis,name,group_homomorphism = None):                       # group_homomorphism param: dictionary, e.g. for mapping double cover to single cover actions
    hom = hom_from_action(group,basis,group_homomorphism)
    Gamma = Representation(group,hom,name)
    Gamma.basis = basis
    Gamma.direction_action = basis[0].direction_action
    return Gamma  
def find_closure_from_matrices(gen_matrices):                   # gen_matrices: dict{name:matrix}
    names = ["I"]
    matrices = [np.eye(len(list(gen_matrices)[0]))]
    hom = {}
    mult_table = {}
    for i in range(len(gen_matrices.keys())):
        names.append(list(gen_matrices.keys())[i])
        matrices.append(list(gen_matrices.values())[i])
    for i in range(len(gen_matrices.keys())+1):    
        hom[names[i]] = matrices[i]
    n = 0
    while n != len(hom):
        n = len(hom)
        for i in range(len(matrices)):
            mult_table[names[i]] = {}
            for j in range(len(matrices)):  
                new_m = matrices[i]@matrices[j]         
                New = True
                for x in matrices:
                    if np.allclose(new_m,x):
                        New = False
                if New:
                    new_name = g.remove_I(names[i] + names[j])
                    names.append(new_name)
                    matrices.append(new_m)
                    hom[new_name] = new_m
    return hom
def mult_table_from_matrices(dict_hom):
    mt = {}
    for a,A in dict_hom.items():
        mt[a] = {}
        for b,B in dict_hom.items():
            for c,C in dict_hom.items():
                if np.allclose(A@B,C):
                    mt[a][b]=c
    return mt      
def rep_trivial(group):
    return {g: np.matrix([1]) for g in group.elements}                  #returns dict of the homomorphism
def rep_determinant(rep):                                               #pass dict of homomorphism, returns new dict of homomorphism
    return {g: np.matrix(np.linalg.det(rep[g])) for g in rep.keys()}
def matrix_regular(i,group):          # formula g_i*g_j = sum g_m (delta_i)_{mj}; (delta_i)_{mj} = 1 for m = k and 0 otherwhise; (delta_i) is rep matrix of g_i in regular rep
    el = list(group.elements)
    M = np.zeros((len(el),len(el)))
    for j in range(len(el)):            #one column of regular rep
        for m in range(len(el)):
            if group.mult_table[el[i]][el[j]] == el[m]:
                M[m][j] = 1
    return M
def rep_regular(group,name):
    el = list(group.elements)
    hom = {}
    for i in range(len(el)):
        hom[el[i]] = matrix_regular(i,group)
    r = Representation(group,hom,name)
    return r
        
def product_rep(A, B):                                  #pass two  representation objects (of same group). Returns new representation object
    assert hasattr(A,"basis")
    assert hasattr(B,"basis") 
    hom = {}
    for g in A.elements:
       hom[g] = np.kron(A.hom[g],B.hom[g])
    r = Representation(A,hom,A.name + "_x_" + B.name)
    r.basis = product_basis(A.basis,B.basis)
    r.direction_action = r.basis[0].direction_action
    return r            
def product_basis(A,B):                 # arrays A,B of objects-> basis of product space as done via Kronecker product (np.kron())
    basis = []
    for i in range(len(A)):
        for j in range(len(B)):
            # basis[i*len(A)+j] = o.TensorProduct(A[i],B[j])
            basis.append(o.TensorProduct(A[i],B[j]))
    return basis

# projectors (adapted from the lecture notes: Prof. Dr. Christoph Lehner: Lattice QCD, chapter 11)
def projector_irrep(rep,irrep):                 #takes rep and EITHER Representation object OR line in char_table ( -> dict object) irrep, returns projector as matrix                            
    ret = 0.0 * rep.hom["I"]
    if isinstance(irrep,Representation):
        for c in rep.classes.keys():
            for g in rep.classes[c]:
                ret += rep.hom[g]*np.conj(irrep.characters[c])
        P = np.matrix(ret * (irrep.characters["I"] / len(rep.elements)))
    else:
        for c in rep.classes.keys():
            for g in rep.classes[c]:
                ret += rep.hom[g]*irrep[c]
        P = np.matrix(ret * (irrep["I"] / len(rep.elements)))
    assert np.linalg.norm(np.matmul(P,P) - P) < 1e-14      
    return P

def project_out_irreps(rep,irreps):         #irreps: list of representation objects; manipulate rep object directly
    P = 0.0 * rep.hom["I"]
    for r in irreps:
        P += projector_irrep(rep,r)
    P = np.matrix(np.eye(len(P))) - P
    assert np.linalg.norm(np.matmul(P,P) - P) < 1e-14  
    for g,m in rep.hom.items():
        rep.hom[g] = P*m*P
    rep.update()

def antisymmetric_projector(n):
    ei = [np.zeros(shape=(n,)) for i in range(n)]
    for i in range(n):
        ei[i][i] = 1        
    # 00 01 10 11
    # (|01> - |10>)(<01| - <10|) / 2
    # v_i v_j -> v_i v_j - v_j v_i
    P = sum([
        np.matrix(np.outer(np.kron(ei[i],ei[j]) - np.kron(ei[j],ei[i]),
                           np.kron(ei[i],ei[j]) - np.kron(ei[j],ei[i])))/2
        for i in range(n) for j in range(i)])        
    assert np.linalg.norm(P*P - P) < 1e-14
    return P

def symmetric_projector(n):
    ei = [np.zeros(shape=(n,)) for i in range(n)]
    for i in range(n):
        ei[i][i] = 1        
    # 00 01 10 11
    # (|01> + |10>)(<01| + <10|) / 2 + |00><00| + |11><11|
    # v_i v_j -> v_i v_j + v_j v_i
    P = sum([
        np.matrix(np.outer(np.kron(ei[i],ei[j]) + np.kron(ei[j],ei[i]),
                           np.kron(ei[i],ei[j]) + np.kron(ei[j],ei[i])))/4.
        for i in range(n) for j in range(n)])        
    assert np.linalg.norm(P*P - P) < 1e-14
    return P

def diagonal_projector(n):
    ei = [np.zeros(shape=(n,)) for i in range(n)]
    for i in range(n):
        ei[i][i] = 1        
    # |00><00| + |11><11|
    # v_i v_j -> v_i v_j \delta_{ij}
    P = sum([
        np.matrix(np.outer(np.kron(ei[i],ei[i]),
                           np.kron(ei[i],ei[i])))
     for i in range(n)])        
    assert np.linalg.norm(P*P - P) < 1e-14
    return P

def invert_projector(projector):
    def _proj(n):
        P = projector(n)
        P = np.matrix(np.eye(n*n)) - P
        assert np.linalg.norm(P*P - P) < 1e-14
        return P
    return _proj
        
def traceless_projector(n):
    ei = [np.zeros(shape=(n,)) for i in range(n)]
    for i in range(n):
        ei[i][i] = 1        
    # (|ii>-\sum_l |ll> / ndim)(<ii|-\sum_l <ll| / ndim) + sum_{i!=j} |ij><ij|
    # v_i v_j -> v_i v_j - v_l v_l \delta_{ij} / ndim
    P = sum([
        np.matrix(np.outer(np.kron(ei[i],ei[j]),np.kron(ei[i],ei[j])))
     for i in range(n) for j in range(n) if i != j])    
    sll = sum([np.kron(ei[l],ei[l]) for l in range(n)]) / n    
    P += sum([
        np.matrix(np.outer(np.kron(ei[i],ei[i]) - sll,np.kron(ei[i],ei[i]) - sll))
     for i in range(n)])    
    assert np.linalg.norm(P*P - P) < 1e-14
    return P

def apply_projectors(projectors,rep):       #input: array of projector functions, representation object which gets changed directly
    n = np.sqrt(len(rep.hom["I"]))          # projectors(arg n) are defined to act on product rep of two nxn reps
    assert n.is_integer()
    n = int(n)
    if not hasattr(projectors,"__iter__"):
        P = projectors(n)
    else:
        P = projectors[0](n)
        for p in projectors[1:]:
            P = P * p(n)  
    apply_projector(P,rep)                 
def apply_projector(P,rep):
    for g,v in rep.hom.items():
        rep.hom[g] = P*rep.hom[g]*P
    rep.update()  

def projected_to_zero(rep):             #returns True if the matrix of identity transformation is zero
    return (abs(rep.hom["I"]) < 1e-8).all()
def find_irreps(rep,group):             # Representation and Group objects -> dict: {irrep : [projector, multiplicity]}
    c_dim = 0
    result = {}
    Rep = rep.copy("Rep")               # change Rep, leave the input Repr. rep unchanged
    print("irreps of ", rep.name , ":")
    for irrep,chars in group.char_table.items():         
        temp = Rep.copy("temp")
        P = projector_irrep(temp,chars)
        print("P")
        print(P)
        apply_projector(P,temp)
        if not projected_to_zero(temp):
            result[irrep] = []
            R = np.matrix(np.eye(len(P))) - P
            apply_projector(R,Rep)                 
            if temp.is_reducible():
                print(irrep, "(reducible inv. subspace)-> occurrence: ", temp.characters["I"] / chars["I"])
                c_dim += temp.characters["I"]
            else:
                print(irrep, "-> occurrence: 1") 
                c_dim += temp.characters["I"]
            result[irrep] = [P,temp.characters["I"] / chars["I"]]
    print("sum of dims: ",  c_dim)
    print("end find_irreps ",rep.name)
    return result

def list_eigvecs(M,e):                  # Matrix M and eigenvalue e -> list of eigenvectors for this eigenvalue
    result = []
    eigvals,eigvecs = np.linalg.eig(M)              
    for i in range(len(eigvals)):
        if abs(eigvals[i] - e) < 1e-8:
            result.append(eigvecs[:,i])
    return result

def list_nonzero_eigvecs(M):            # Matrix M -> dict {eigenvalue: [eigenvectors]}; the vectors in np.matrix format 
    result = {}
    c = 0  
    eigvals,eigvecs = np.linalg.eig(M)          
    for i in range(len(eigvals)):
        if not abs(eigvals[i]) < 1e-8:  
            ev = o.match_in_list(eigvals[i],list(result.keys()))   
            if ev == None: 
                result[eigvals[i]] = [eigvecs[:,i]]
            else:
                if hasattr(result[ev],"__len__") and not isinstance(result[ev],np.matrix):
                    result[ev].append(np.matrix(eigvecs[:,i]))
                else:
                    result[ev] = [result[ev]]
                    result[ev].append(np.matrix(eigvecs[:,i]))
            c += 1
    if not c:
        print("no nonzero EVs found")
    return result

def nonzero_evs_different(ev_dict,mult):                     # to be used with list_nonzero_eigvals() -> True if each Eval which is nonzero appears only up to given multiplicity mult
    if abs(1-mult)<1e-8:       
        return all([not isinstance(evecs,list) for evecs in ev_dict.values()]) and nonzeros_all_different(ev_dict.keys())
    else: 
        mult = round(mult)
        if len(list(ev_dict.values())) == mult and nonzeros_all_different(ev_dict.keys()):
            print("evs appear ", mult, "times.")
        return len(list(ev_dict.values())) == mult          #and nonzeros_all_different(ev_dict.keys())


def nonzeros_all_different(nums):                   # nonzero list or similar -> True if all nonzero elements are different, False otherwhise
    nums = list(nums)                               # CAREFUL: in combination with list_nonzero_eigvals(), use all_evs_different() 
    while len(nums) > 0:
        ev = nums.pop()
        if not abs(ev) < o.num_tol:
            for comp in nums:
                if abs(ev - comp) < o.num_tol:
                    return False
    return True
def make_subspaces_comparable(P1,P2,R1,R2):         # projectors with multiplicity of irrep as P1 = [p1,m1], P2 = [p2,m2], Representation R1,R2 of same underlying group -> Evecs of M(g)P1 and M(g)P2 for same g such that all nonzero Evals are different
                                                # NB: this makes these Evecs transform in the same way, e.g., but not literally, "like the x coordinate of T1m" (see Aaron's note)
                                                #     for this literal trafo behavior, use <irrepname>_identify_components functions
                                                # NB2: similar matrices have the same Evals, but not the same Evecs (wiki)
    paired_evecs = {}
    candidates = [c for c in R1.classes.keys() if c not in {"I","Inv"}]
    for c in candidates:
        temp1 = list_nonzero_eigvecs(np.matmul(R1.hom[c],P1[0]))                                # reminder: -> dict {Evals : Evecs}
        temp2 = list_nonzero_eigvecs(np.matmul(R2.hom[c],P2[0]))
        if nonzero_evs_different(temp1,P1[1]) and nonzero_evs_different(temp2,P2[1]):                   # Evals all different
            ev_temp1 = sorted(temp1.keys())
            ev_temp2 = sorted(temp2.keys())
            print("all different evs.")
            if np.allclose(ev_temp1,ev_temp2,1e-5,1e-10):                                                   # same Evals
                print("success for element " , c)
                for i in range(len(ev_temp1)):
                    paired_evecs[ev_temp1[i]] = [temp1[ev_temp1[i]],temp2[ev_temp2[i]]]
                return paired_evecs
            else:
                print("fail: not all matching eigenvalues.")
        else:
            print("fail: not all different eigenvalues for both P1 and P2 and class element ", c)
    return None 
def A_identify_components(P):
    components = list_nonzero_eigvecs(P[0])
    comps = {}
    for ev, evecs in components.items():
        if abs(ev-1) < o.num_tol:
            comps["1"] = gram_schmidt(*evecs)
            for i in range(len(comps["1"])):
                comps["1"][0] = rotate_to_real_valued(comps["1"][0])
                comps["1"][0] = make_first_entry_pos(comps["1"][0])        
    return comps
def T1_identify_components(P,Rep):          # P = [proj,mult] -> dict {"x" : v or [v_1,v_2, ..], "y" : ..} such that same array entries, e.g. dict["x"][0],dict["y"][0],.. form an inv. subspace
    components = {}
    M = np.matmul(Rep.hom["Rot0"],P[0])
    ev_x = list_nonzero_eigvecs(M)            # evecs of Rep(Rot0).P with ev 1 transform like x coordinate 
    for e in list(ev_x.keys()):
        if abs(e-1)< 1e-8:
            components["x"] = gram_schmidt(*ev_x[e])
            for i in range(len(components["x"])):
                # components["x"][i] = rotate_to_real_valued(components["x"][i])
                components["x"][i] = make_first_entry_pos(components["x"][i])
    components["y"] = []
    components["z"] = [] 
    # define such that (x,y,z) form a right system. Depending on direction_action: x->y and y->z done by opposite rotations since f.Rot_i^{right} = Rot_i^{left}^{-1}.f   
    if Rep.direction_action == "right":
        for x in components["x"]:
            y = np.matmul(Rep.hom["Rot2"].T,x)
            components["y"].append(y)
            components["z"].append(np.matmul(Rep.hom["Rot0"].T,y))
    elif Rep.direction_action == "left":
        for x in components["x"]:
            y = np.matmul(Rep.hom["Rot2"],x)
            components["y"].append(y)
            components["z"].append(np.matmul(Rep.hom["Rot0"],y)) 
    else:
        print("PROBLEM: in T1_identify_components: neither left nor right action applicable.")
    return components

def E_identify_components(P,Rep):          # P = [proj,mult] -> dict {"xx-yy": .., "xx-zz": ..} such that same array entries, e.g. dict["x"][0],dict["y"][0],.. form an inv. subspace
    components = {}
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               # Eigvec with EV -1 transforms like (x_1,x_2) - (y_1,y_2) in T1m_x_T1m
    for e in list(P_orientation.keys()):
        if abs(e+1)< 1e-8:
            components["xx-yy"] = gram_schmidt(*P_orientation[e])
            for i in range(len(components["xx-yy"])):
                components["xx-yy"][i] = rotate_to_real_valued(components["xx-yy"][i])
                components["xx-yy"][i] = make_first_entry_pos(components["xx-yy"][i]) 
    components["xx-zz"] = []
    # direction_action no factor here, since this is taken care of in Rep.hom via transposition: there, every right-action problem is transformed in a left-action problem
    for x in components["xx-yy"]:
        components["xx-zz"].append(np.matmul(Rep.hom["Rot0"],x))                # rotate xx-yy to xx-zz
    return components

def T2_identify_components(P,Rep):          # P = [proj,mult] -> dict {"xx-yy": .., "xx-zz": ..} such that same array entries, e.g. dict["x"][0],dict["y"][0],.. form an inv. subspace
    components = {}
    # P_t1 = np.matmul(Rep.hom["Rot2"],P[0])
    # np.save("../results/projectors/"+Rep.name+"_T2m_t1",P_t1)
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               # Eigvec with EV -1 transforms like (x_1,x_2) - (y_1,y_2) in T1m_x_T1m
    for e in list(P_orientation.keys()):
        if abs(e+1)< 1e-8:
            components["xy+yx"] = P_orientation[e]
            for i in range(len(components["xy+yx"])-1):
                components["xy+yx"][0] = subtract_for_zero_entries(components["xy+yx"][0],components["xy+yx"][i+1])     #makes for "nicer" vectors
            components["xy+yx"][0] = normalize(components["xy+yx"][0])
            components["xy+yx"][0] = rotate_to_real_valued(components["xy+yx"][0])
            components["xy+yx"] = gram_schmidt(*components["xy+yx"])                # orthonormalize
            for i in range(len(components["xy+yx"])):
                
                components["xy+yx"][i] = make_first_entry_pos(components["xy+yx"][i])       # introduce sign convention for vectors
    components["xz+zx"] = []
    components["yz+zy"] = [] 
    if Rep.direction_action == "left":
        for a in components["xy+yx"]:
            components["xz+zx"].append(np.matmul(Rep.hom["Rot0"],a))                # rotate xy+yx to xz+zx   
        for b in components["xz+zx"]:
            components["yz+zy"].append(np.matmul(Rep.hom["Rot2"],b))                # rotate xz+zx to yz+zy
    elif Rep.direction_action == "right":
        for a in components["xy+yx"]:
            components["xz+zx"].append(np.matmul(Rep.hom["Rot0"].T,a))                # rotate xy+yx to xz+zx   
        for b in components["xz+zx"]:
            components["yz+zy"].append(np.matmul(Rep.hom["Rot2"].T,b))                # rotate xz+zx to yz+zy
    else:
        print("PROBLEM: in T1_identify_components: neither left nor right action applicable.")
    return components
def study_irreps(rep,group,filename):                               # find all irreps(like find_irreps()), then decompose all subspaces and write results to file; 
                                                                    # filename: absolute or relative path, must end in desired format,e.g. ".txt" 
                                                                    # returns {irrep : {basisvector_name : [such vectors]}}
    c_dim = 0
    P_irrep = {}
    Rep = rep.copy("Rep")               # change Rep, leave the input Repr. rep unchanged
    print("irreps of ", rep.name , ":")
    for irrep,chars in group.char_table.items():         
        temp = Rep.copy("temp")
        P = projector_irrep(temp,chars)
        np.save("../results/projectors/"+rep.name+"_"+irrep,P)
        if irrep == "T2m":
            P_t1 = np.matmul(Rep.hom["Rot2"],P)
            np.save("../results/projectors/"+Rep.name+"_T2m_t1",P_t1)
        apply_projector(P,temp)
        if not projected_to_zero(temp):
            P_irrep[irrep] = []
            R = np.matrix(np.eye(len(P))) - P
            apply_projector(R,Rep)                 
            if temp.is_reducible():
                print(irrep, "(reducible inv. subspace)-> occurrence: ", temp.characters["I"] / chars["I"])
                c_dim += temp.characters["I"]
            else:
                print(irrep, "-> occurrence: 1") 
                c_dim += temp.characters["I"]
            P_irrep[irrep] = [P,temp.characters["I"] / chars["I"]]
    print("sum of dims: ",  c_dim)
    

    inv_subspaces = {}
    f = open(filename, "w")
    f.write("Basis of " + rep.name + ": \n") 
    for b in rep.basis:
        f.write(b.name + "\n")
    f.write("Dimension: " + str(len(rep.basis)) + "\n")
    for key,P_n_pair in P_irrep.items():
        if key == "A1m" or key == "A1p":   
            vecs = A_identify_components(P_n_pair)                
            if key == "A1m":
                f.write("A1m subspace:"+"\n") 
            else:
                f.write("A1p subspace:"+"\n")
            inv_subspaces[key] = vecs
            f.write(str(vecs)+"\n")
        if key == "T1m" or key == "T1p":
            vecs = T1_identify_components(P_n_pair,rep)
            if key == "T1m":
                f.write("T1m subspace:"+"\n")
            else:
                f.write("T1p subspace:"+"\n")
            inv_subspaces[key] = vecs
            f.write(str(vecs)+"\n")
        if key == "A2m" or key == "A2p":  
            vecs = A_identify_components(P_n_pair)                 
            if key == "A2m":
                f.write("A2m subspace:"+"\n") 
            else:
                f.write("A2p subspace:"+"\n")
            inv_subspaces[key] = vecs
            f.write(str(vecs)+"\n")
        if key == "T2m" or key == "T2p":
            vecs = T2_identify_components(P_n_pair,rep)
            if key == "T2m":
                f.write("T2m subspace:"+"\n")
            else:
                f.write("T2p subspace:"+"\n")
            inv_subspaces[key] = vecs
            f.write(str(vecs)+"\n")
        if key == "Em" or key == "Ep":
            vecs = E_identify_components(P_n_pair,rep)
            if key == "Em":
                f.write("Em subspace:"+"\n")
            else:
                f.write("Ep subspace:"+"\n")
            inv_subspaces[key] = vecs
            f.write(str(vecs)+"\n")
    f.close()
    return inv_subspaces

def scalar_prod(x,y):                           #returne (complex or real) scalar product between x and y 
    assert len(x) == len(y)
    temp = 0
    for i in range(len(x)):
        u = complex(x[i])
        v = np.conj(complex(y[i]))
        temp += u * v
    return temp            
def normalize(x):                           # returns normalized vector(array etc.)
        return np.sqrt(scalar_prod(x,x))**(-1)*x
def gram_schmidt(*v):                      # returns array of orthonormal vectors 
    def proj(x,y):
        nom = scalar_prod(y,x)
        denom = scalar_prod(x,x)
        return (nom/denom)*x            
    def sum_proj(x,y,i): 
        if i < 0:
            return 0
        return proj(x[i],y) + sum_proj(x,y,i-1)
    u = []
    e = []
    u.append(v[0])
    e.append(normalize(v[0]))
    for i in range(len(v)-1):
        temp = v[i+1] - sum_proj(u,v[i+1],i)
        u.append(temp)
        e.append(normalize(temp))
    return e

def make_first_entry_pos(v):                # looks for first nonzero entry in v; multiplies v with (-1) if this entry is negative
    assert hasattr(v,"__len__")
    i = 0
    while abs(v[i]) < o.num_tol:
        i += 1
    if np.real(v[i]) < 0:
        for j in range(len(v)):
            v[j] = (-1)*v[j]
    return v    
def subtract_for_zero_entries(u,v):             #subtracts vector v from u such that u obtains more/most zero entries; returns u
    assert len(u) == len(v)
    new_u = []
    i = 0
    while True:
        if abs(u[i]) > o.num_tol:
            if abs(v[i]) > o.num_tol:
                break
        i += 1
        if i > len(u)+1:
            print("In subtract_for_zero_entries: only entries of zero.")
            i = 0
            break
    f = complex(10000000000*u[i]/(100000000000*v[i]))
    for j in range(len(u)):
        new_u.append([complex(u[j])-f*complex(v[j])])
    new_u = np.matrix(new_u)
    return new_u
def first_imag_value(v):
    i = 0    
    while i < len(v)-1:
        if abs(np.imag(v[i])) > o.num_tol:
            break
        i += 1
    if i == len(v) and abs(np.imag(v[i])) < o.num_tol:
        return None
    return i    
def rotate_to_real_valued(v):           # takes vector, applies scalar phase e^{i\phi} to it in an attempt to make all entries purely real
    idx = first_imag_value(v)
    if idx == None:
        return v
    phi = np.arctan2(np.imag(v[idx]),np.real(v[idx]))
    for i in range(len(v)):        
        v[i] = v[i]*np.exp(-phi*1j)
        f = 0
        if abs(np.imag(v[i])) > o.num_tol:
            f += 1
    if f:
        print("WARNING: in rotate_to_real_valued: imaginary part remaining.")
    return v
# def rotate_to_real_valued(v):
#     phi = [np.arctan2(np.imag(v[i]),np.real(v[i])) for i in range(len(v))]
#     print(phi)
#     for i in range(len(phi)):
#         if abs(phi[i]-phi[0]) > o.num_tol:
#             # print(phi[i]-phi[0])
#             print("FAIL: in rotate_to_real_valued: no scalar phase is appropriate.")
#             return v
#     for i in range(len(v)):        
#         v[i] = v[i]*np.exp(-phi[0]*1j)
#     return v

    


    
    
