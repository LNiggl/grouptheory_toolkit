import numpy as np

from . import groups as g
from . import objects as o
from . import tools as t
from .definitions import num_tol, num_rtol,sigma


class Representation:                                       # used to access: elements, matrices, characters, ..
    def __init__(self,group,dict_matrix,name):              # initialize with dict: g -> M(g)
        self.group = group                                 
        self.elements = group.elements
        self.mult_table = group.mult_table
        self.classes = group.classes

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
                assert eps < num_tol
                
    def is_reducible(self):
        return g.inner_product(self.characters,self.characters,self.group) - 1 > num_tol
    
    def round_chars(self):
        for c,chi in self.characters.items():
            if chi - round(chi) < num_tol:
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
            if not hasattr(basis[j],"projection"):
                if basis[j].is_equal_to(acted_basis[i]):                            # "if e_i gets transformed to e_j via action A" -> ge_i = M(A)_{ji} e_j == +/- e_j
                        M[j][i] = 1                                                 # --> M(A)_{ji} == +/- 1
                elif basis[j].is_negative_of(acted_basis[i]):
                    M[j][i] = -1            
                else:
                    M[j][i] = 0
            else:
                k = basis[j].projection(acted_basis[i])
                M[j][i] = k            
    # Note: with right action: M is technically defined as the transpose of the above, but application of trafo must be by right multiplication of that matrix M.T,                          
    # therefore again by left multiplication of (M.T).T -> as implemented, the difference left/right is not noticeable in the matrices, but only in check_if_homomorphism()
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
                    new_name = t.remove_I(names[i] + names[j])
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
    return {g: np.matrix([1]) for g in group.elements}                  # returns dict of the homomorphism

def rep_determinant(rep):                                               # pass dict of homomorphism, returns new dict of homomorphism
    return {g: np.matrix(np.linalg.det(rep[g])) for g in rep.keys()}

def matrix_regular(i,group):            # formula g_i*g_j = sum g_m (delta_i)_{mj}; (delta_i)_{mj} = 1 for m = k and 0 otherwhise; (delta_i) is rep matrix of g_i in regular rep
    el = list(group.elements)
    M = np.zeros((len(el),len(el)))
    for j in range(len(el)):            # one column of regular rep
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
        
def product_rep(A, B):                                  # pass two Representation objects (of same group). Returns new Representation object
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
    assert np.linalg.norm(np.matmul(P,P) - P) < num_tol     
    return P

def project_out_irreps(rep,irreps):         #irreps: list of representation objects; manipulate rep object directly
    P = 0.0 * rep.hom["I"]
    for r in irreps:
        P += projector_irrep(rep,r)
    P = np.matrix(np.eye(len(P))) - P
    assert np.linalg.norm(np.matmul(P,P) - P) < num_tol  
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
    assert np.linalg.norm(P*P - P) < num_tol
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
    assert np.linalg.norm(P*P - P) < num_tol
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
    assert np.linalg.norm(P*P - P) < num_tol
    return P

def invert_projector(projector):
    def _proj(n):
        P = projector(n)
        P = np.matrix(np.eye(n*n)) - P
        assert np.linalg.norm(P*P - P) < num_tol
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
    assert np.linalg.norm(P*P - P) < num_tol
    return P

def apply_projectors(projectors,rep):       # input: array of projector functions, representation object which gets changed directly
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

def projected_to_zero(rep):             # returns True if the matrix of identity transformation is zero
    return (abs(rep.hom["I"]) < num_tol).all()

def list_eigvecs(M,e):                  # Matrix M and eigenvalue e -> list of eigenvectors for this eigenvalue
    result = []
    eigvals,eigvecs = np.linalg.eig(M)              
    for i in range(len(eigvals)):
        if abs(eigvals[i] - e) < num_tol:
            result.append(eigvecs[:,i])
    return result

def list_nonzero_eigvecs(M):            # Matrix M -> dict {eigenvalue: [eigenvectors]}; the vectors in np.matrix format 
    result = {}
    c = 0  
    eigvals,eigvecs = np.linalg.eig(M)          
    for i in range(len(eigvals)):
        if not abs(eigvals[i]) < num_tol:  
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
    if abs(1-mult)<num_tol:       
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

##### special functions for projecting to irreps of O_h or O or their double covers #####

def A_identify_components(P):                       # P = [proj,mult]. Projector and multiplicity of irrep                 
    components = list_nonzero_eigvecs(P[0])
    comps = {}
    for ev, evecs in components.items():
        if abs(ev-1) < o.num_tol:
            comps["1"] = gram_schmidt(*evecs)
            for i in range(len(comps["1"])):
                comps["1"][0] = rotate_to_real_valued(comps["1"][0])
                comps["1"][0] = make_first_entry_pos(comps["1"][0])        
    return comps

def T1_identify_components(P,Rep):              # P = [proj,mult] -> dict {"x" : v or [v_1,v_2, ..], "y" : ..} such that same array entries, e.g. dict["x"][0],dict["y"][0],.. form an inv. subspace
    components = {}
    M = np.matmul(Rep.hom["Rot0"],P[0])
    ev_x = list_nonzero_eigvecs(M)              # evecs of Rep(Rot0).P with ev 1 transform like x coordinate 
    for e in list(ev_x.keys()):
        if abs(e-1)< num_tol:
            components["x"] = ev_x[e]
            
            j = 0
            while j < len(components["x"]):
                for i in range(j+1,len(components["x"])):
                    components["x"][j] = subtract_for_zero_entries(components["x"][j],components["x"][i])       # makes for "nicer" vectors 
                components["x"][j] = make_first_entry_unity(components["x"][j])
                if j == 0: 
                    components["x"]= gram_schmidt(*components["x"][j:])                                         # orthonormalize
                else:
                    temp = components["x"][:]
                    components["x"] = temp[:j]
                    components["x"]+=gram_schmidt(*temp[j:])
    components["y"] = []
    components["z"] = [] 
    # defined such that basis (x,y,z) would always form a LEFT HANDED SYSTEM(same treatment of the matrices, difference only in check_if_homomorphism())
    if Rep.direction_action == "right" or  Rep.direction_action == "left":
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
        if abs(e+1)< num_tol:
            components["xx-yy"] = P_orientation[e]
            j = 0
            while j < len(components["xx-yy"]):
                for i in range(j+1,len(components["xx-yy"])):
                    components["xx-yy"][j] = subtract_for_zero_entries(components["xx-yy"][j],components["xx-yy"][i])       # makes for "nicer" vectors 
                components["xx-yy"][j] = make_first_entry_unity(components["xx-yy"][j]) 
                if j == 0: 
                    components["xx-yy"]= gram_schmidt(*components["xx-yy"][j:])                                             # orthonormalize
                else:
                    temp = components["xx-yy"][:]
                    components["xx-yy"] = temp[:j]
                    components["xx-yy"]+=gram_schmidt(*temp[j:])
                               # orthonormalize
                j += 1
    components["xx-zz"] = []
    # defined such that basis (x,y,z) would always form a LEFT HANDED SYSTEM(same treatment of the matrices, difference only in check_if_homomorphism())
    if Rep.direction_action == "right" or  Rep.direction_action == "left":
        for x in components["xx-yy"]:
            components["xx-zz"].append(np.matmul(Rep.hom["Rot0"],x))                # rotate xx-yy to xx-zz
    else:
        print("PROBLEM: in T1_identify_components: neither left nor right action applicable.")
    return components

def T2_identify_components(P,Rep):          # P = [proj,mult] -> dict {"xx-yy": .., "xx-zz": ..} such that same array entries, e.g. dict["x"][0],dict["y"][0],.. form an inv. subspace
    components = {}
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               # Eigvec with EV -1 transforms like (x_1,x_2) - (y_1,y_2) in T1m_x_T1m
    for e in list(P_orientation.keys()):
        if abs(e+1)< num_tol:
            components["xy+yx"] = P_orientation[e]
            j = 0
            while j < len(components["xy+yx"]):
                for i in range(j+1,len(components["xy+yx"])):
                    components["xy+yx"][j] = subtract_for_zero_entries(components["xy+yx"][j],components["xy+yx"][i])       # makes for "nicer" vectors 
                components["xy+yx"][j] = make_first_entry_unity(components["xy+yx"][j])
                if j == 0: 
                    components["xy+yx"]= gram_schmidt(*components["xy+yx"][j:])                                             # orthonormalize
                else:
                    temp = components["xy+yx"][:]
                    components["xy+yx"] = temp[:j]
                    components["xy+yx"] +=gram_schmidt(*temp[j:])
                j += 1
    components["xz+zx"] = []
    components["yz+zy"] = []
    # defined such that basis (x,y,z) would always form a LEFT HANDED SYSTEM(same treatment of the matrices, difference only in check_if_homomorphism())
    if Rep.direction_action == "left" or Rep.direction_action == "right":
        for a in components["xy+yx"]:
            components["xz+zx"].append(np.matmul(Rep.hom["Rot0"],a))                # rotate xy+yx to xz+zx   
        for b in components["xz+zx"]:
            components["yz+zy"].append(np.matmul(Rep.hom["Rot2"],b))                # rotate xz+zx to yz+zy
    else:
        print("PROBLEM: in T1_identify_components: neither left nor right action applicable.")
    return components

def G1_identify_components(P,Rep):
    components = {}
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               
    for e in list(P_orientation.keys()):
        a = (1/np.sqrt(2)*(1-1j))                       # Eigvec with EV a transforms like Weyl-spinor with components (1,0)
        if abs(e-a)< num_tol:
            components["s_1"] = P_orientation[e]
            j = 0
            while j < len(components["s_1"]):
                for i in range(j+1,len(components["s_1"])):
                    components["s_1"][j] = subtract_for_zero_entries(components["s_1"][j],components["s_1"][i])         # makes for "nicer" vectors 
                components["s_1"][j] = make_first_entry_unity(components["s_1"][j])
                if j == 0: 
                    components["s_1"]= gram_schmidt(*components["s_1"][j:])                                             # orthonormalize
                else:
                    temp = components["s_1"][:]
                    components["s_1"] = temp[:j]
                    components["s_1"] +=gram_schmidt(*temp[j:])
                j += 1
    components["s_2"] = []
    for x in components["s_1"]:
        components["s_2"].append(np.matmul(Rep.hom["Rot1Rot1"],x))                  # rotate to Weyl-spinor (0,1)    
    return components

def G2_identify_components(P,Rep):
    components = {}
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               
    for e in list(P_orientation.keys()):
        a = -(1/np.sqrt(2)*(1+1j))                                                  # define eigvec with EV as first basis vector t_1
        if abs(e-a)< num_tol:
            components["t_1"] = P_orientation[e]
            ##
            j = 0
            while j < len(components["t_1"]):
                for i in range(j+1,len(components["t_1"])):
                    components["t_1"][j] = subtract_for_zero_entries(components["t_1"][j],components["t_1"][i])     # makes for "nicer" vectors 
                components["t_1"][j] = make_first_entry_unity(components["t_1"][j])
                if j == 0: 
                    components["t_1"]= gram_schmidt(*components["t_1"][j:])                                         # orthonormalize
                else:
                    temp = components["t_1"][:]
                    components["t_1"] = temp[:j]
                    components["t_1"] +=gram_schmidt(*temp[j:])
                j += 1
    components["t_2"] = []
    for x in components["t_1"]:
        components["t_2"].append(np.matmul(Rep.hom["Rot1Rot1"],x))                              # rotate to Weyl-spinor (0,1)    
    return components

def H_identify_components(P,Rep):
    components = {}
    P_orientation = list_nonzero_eigvecs(np.matmul(Rep.hom["Rot2"],P[0]))               
    for e in list(P_orientation.keys()):
        a = (1/np.sqrt(2)*(1-1j))                   # Eigvec with EV a identified as h_1, with EV -a identified as h_2
        b = (1/np.sqrt(2)*(1+1j))                   # h_3 and h_4 defined via an operation from h_1 and h_2. h_3 has EV b and h_4 has EV -b 
        if abs(e-a)< num_tol:
            components["h_1"] = P_orientation[e]
            j = 0
            while j < len(components["h_1"]):
                for i in range(j+1,len(components["h_1"])):
                    components["h_1"][j] = subtract_for_zero_entries(components["h_1"][j],components["h_1"][i])     # makes for "nicer" vectors 
                components["h_1"][j] = make_first_entry_unity(components["h_1"][j])
                if j == 0: 
                    components["h_1"]= gram_schmidt(*components["h_1"][j:])                                         # orthonormalize
                else:
                    temp = components["h_1"][:]
                    components["h_1"] = temp[:j]
                    components["h_1"] +=gram_schmidt(*temp[j:])
                j += 1
        if abs(e+a)< num_tol:
            components["h_2"] = P_orientation[e]
            j = 0
            while j < len(components["h_2"]):
                for i in range(j+1,len(components["h_2"])):
                    components["h_2"][j] = subtract_for_zero_entries(components["h_2"][j],components["h_2"][i])     # makes for "nicer" vectors 
                components["h_2"][j] = make_first_entry_unity(components["h_2"][j])
                if j == 0: 
                    components["h_2"]= gram_schmidt(*components["h_2"][j:])                                         # orthonormalize
                else:
                    temp = components["h_2"][:]
                    components["h_2"] = temp[:j]
                    components["h_2"] +=gram_schmidt(*temp[j:])
                j += 1
        if abs(e-b)< num_tol:
            components["h_3"] = P_orientation[e]
            j = 0
            while j < len(components["h_3"]):
                for i in range(j+1,len(components["h_3"])):
                    components["h_3"][j] = subtract_for_zero_entries(components["h_3"][j],components["h_3"][i])     # makes for "nicer" vectors 
                components["h_3"][j] = make_first_entry_unity(components["h_3"][j])
                if j == 0:  
                    components["h_3"]= gram_schmidt(*components["h_3"][j:])                                         # orthonormalize
                else:
                    temp = components["h_3"][:]
                    components["h_3"] = temp[:j]
                    components["h_3"] +=gram_schmidt(*temp[j:])
                j += 1
        if abs(e+b)< num_tol:
            components["h_4"] = P_orientation[e]
            j = 0
            while j < len(components["h_4"]):
                for i in range(j+1,len(components["h_4"])):
                    components["h_4"][j] = subtract_for_zero_entries(components["h_4"][j],components["h_4"][i])     # makes for "nicer" vectors 
                components["h_4"][j] = make_first_entry_unity(components["h_4"][j])
                if j == 0: 
                    components["h_4"]= gram_schmidt(*components["h_4"][j:])                                         # orthonormalize
                else:
                    temp = components["h_4"][:]
                    components["h_4"] = temp[:j]
                    components["h_4"] +=gram_schmidt(*temp[j:])
                j += 1
    return components
    
def study_irreps(rep,group,filename, path_save_projectors = None):              # find all irreps via projection operators, then write results to file
                                                                                # filename: absolute or relative path ending in ".txt" 
                                                                                # returns {irrep : {basisvector_name : [basis vectors]}}
                                                                                # path_save_projectors: saves .npy file of projection matrices; if None: no save
    c_dim = 0
    P_irrep = {}
    Rep = rep.copy("Rep")                                               # to leave the input Repr. rep unchanged
    f = open(filename, "w")
    print("irreps of ", rep.name , ":")
    f.write("irreps of "+ rep.name + ":" + "\n")
    for irrep,chars in group.char_table.items():         
        temp = Rep.copy("temp")
        P = projector_irrep(temp,chars)

        if path_save_projectors != None:
            np.save(path_save_projectors+rep.name+"_"+irrep,P)          # projectors saved

        apply_projector(P,temp)
        if not projected_to_zero(temp):                
            if temp.is_reducible():
                print(irrep, "(reducible inv. subspace)-> occurrence: ", temp.characters["I"] / chars["I"])
                f.write(irrep + "(reducible inv. subspace)-> occurrence: " + str(temp.characters["I"] / chars["I"]) + "\n")
                c_dim += temp.characters["I"]
            else:
                print(irrep, "-> occurrence: 1")
                f.write(irrep+ "-> occurrence: 1" + "\n") 
                c_dim += temp.characters["I"]
            P_irrep[irrep] = [P,temp.characters["I"] / chars["I"]]
    print("sum of dims: " +  str(c_dim))
    f.write("sum of dims: " +  str(c_dim))
    
    f.write("\n\nBasis of " + rep.name + ": \n") 
    for b in rep.basis:
        f.write(b.name + "\n")
    f.write("Dimension: " + str(len(rep.basis)) + "\n")

    inv_subspaces = {}
    for key,P_n_pair in P_irrep.items(): 
        if "A1" in key: 
            vecs = A_identify_components(P_n_pair)                
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "T1" in key:
            vecs = T1_identify_components(P_n_pair,rep)    
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "A2" in key:  
            vecs = A_identify_components(P_n_pair)                    
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "T2" in key:
            vecs = T2_identify_components(P_n_pair,rep)    
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "E" in key:
            vecs = E_identify_components(P_n_pair,rep)
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "G1" in key:
            vecs = G1_identify_components(P_n_pair,rep)
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")    
        if "G2" in key:
            vecs = G2_identify_components(P_n_pair,rep)
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")
        if "H" in key:
            vecs = H_identify_components(P_n_pair,rep)
            inv_subspaces[key] = vecs
            f.write(key + " subspace:\n")
            f.write(str(vecs)+"\n")      
    f.close()
    return inv_subspaces

def scalar_prod(x,y):                           #return (complex or real) scalar product between x and y 
    x = np.array(x)
    y = np.array(y)
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

def make_first_entry_unity(v):                
    assert hasattr(v,"__len__")
    i = 0
    while abs(v[i]) < o.num_tol:
        i += 1
    norm = 0 + v[i]
    for j in range(len(v)):
        v[j] = v[j]/norm
    return v
    
def subtract_for_zero_entries(u,v):             #subtracts vector v from u such that u obtains more/most zero entries; returns u
    assert len(u) == len(v)
    new_u = []
    i = find_index_largest_entry(v)
    f = complex(u[i]/v[i])
    for j in range(len(u)):
        new_u.append([complex(u[j])-f*complex(v[j])])
    new_u = np.matrix(new_u)
    return new_u

def find_index_largest_entry(v):
    idx = 0
    for i in range(len(v)):
        if abs(v[i]) > abs(v[idx]) and abs(v[i]-v[idx] > o.num_tol):
            idx = i
    return idx

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


### double cover methods

def R_entry(A,j,k):
    return 1/2*np.trace(sigma[j]@A@sigma[k]@np.linalg.inv(A))
def R(A):                                                               # 2-to-1 homomorphism from Double-cover matrices gen. by A[0] to A[2] to single-cover matrices
    M = np.zeros((3,3),dtype = np.complex128)
    for j in range(3):
        for k in range(3):
            M[j][k] = R_entry(A,j,k)
    return M

def group_hom_via_reps(Rep1,Rep2,F = None):        #F:function relating Rep1 matrices to Rep2 matrices, e.g. Rep1: Double cover of O, Rep2: O, F: R(A); returns dict{groupelement1:groupelement2}
    group_hom = {}
    for g1,m1 in Rep1.hom.items():
        if F != None:
            n = R(m1)
        else: 
            n = m1
        matches = []
        for g2,m2 in Rep2.hom.items():
            if np.allclose(n,m2):
                matches.append(g2)
        # assert len(matches) == 1            # any-to-one homomorphism
        group_hom[g1] = matches
    return group_hom



    


    
    
