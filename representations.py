import numpy as np
import groups as g
import objects as o
class Representation(g.Group):
    # hom = {}                           #dict for homomorphism: G -> GL(n)
    # characters = {}
    def __init__(self,group,dict_matrix,name):             #initialize with dict: g -> M(g)
        self.group = group                                 # VERY ugly structure, but needed to make copy() work. RETHINK: inheritance and Group and Representation constructors
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
        # el = self.elements.copy()
        # cl = self.classes.copy()
        # mt = self.mult_table.copy()
        # inv = self.inverse_list.copy()
        hom = self.hom.copy()
        new_rep = Representation(self.group,hom,name)
        # chars = self.characters.copy()
        return new_rep
    def update(self):
        for c in self.classes:
            self.characters[c] = np.trace(self.hom[c])
        
    def check_if_homomorphism(self):
        for g in self.elements:
            for h in self.elements:
                eps = np.linalg.norm(np.matmul(self.hom[g],self.hom[h]) - self.hom[self.mult_table[g][h]])
                assert eps < 1e-13
    def is_reducible(self):
        return g.inner_product(self.characters,self.characters,self.group) - 1 > 1e-5

    

    def round_chars(self):
        for c,chi in self.characters.items():
            if chi - round(chi) < 1e-12:
                self.characters[c] = round(chi)

def matrix_from_action(A,basis):
    d = len(basis)
    M = np.zeros((d,d))
    acted_basis = []
    for b in basis:
        x = b.copy()
        x.action(A)
        acted_basis.append(x)
    for i in range(d):
        for j in range(d):
            if basis[j].is_equal_to(acted_basis[i]):               # "if e_i gets transformed to e_j via action A" -> ge_i = M(A)_{ji} e_j == +/- e_j
                M[j][i] = 1                                     # --> M(A)_{ji} == +/- 1
            if basis[j].is_negative_of(acted_basis[i]):
                M[j][i] = -1                
    return M

def hom_from_action(group,basis):                                   # returns the dict for Rep.hom
    return {g : matrix_from_action(g,basis) for g in group.elements}
def rep_from_action(group,basis,name):
    hom = hom_from_action(group,basis)
    Gamma = Representation(group,hom,name)
    return Gamma        
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
    hom = {}
    for g in A.elements:
       hom[g] = np.kron(A.hom[g],B.hom[g])
    r = Representation(A,hom,A.name + "_x_" + B.name)
    return r            

# projectors (adapted from the lecture notes: Prof. Dr. Christoph Lehner: Lattice QCD, chapter 11)
def projector_irrep(rep,irrep):                 #takes rep and EITHER Representation object OR line in char_table ( -> dict object) irrep, returns projector as matrix                            
    ret = 0.0 * rep.hom["I"]
    if isinstance(irrep,Representation):
        for c in rep.classes.keys():
            for g in rep.classes[c]:
                ret += rep.hom[g]*irrep.characters[c]
        P = np.matrix(ret * (irrep.characters["I"] / len(rep.elements)))
    # elif isinstance(irrep,dict):
    else:
        for c in rep.classes.keys():
            for g in rep.classes[c]:
                # print(rep.hom[g])
                # print(irrep)
                # print(irrep[c])
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
    # for g,v in rep.hom.items():
    #     rep.hom[g] = P*rep.hom[g]*P
    # rep.update()  
    apply_projector(P,rep)                 
def apply_projector(P,rep):
    for g,v in rep.hom.items():
        rep.hom[g] = P*rep.hom[g]*P
    rep.update()  

def projected_to_zero(rep):             #returns True if the matrix of identity transformation is zero
    return (abs(rep.hom["I"]) < 1e-8).all()
# Rot2 = np.matrix([[0,-1,0],
#         [1,0,0],
#         [0,0,1]])
def find_irreps(rep,group):             # Representation and Group objects -> dict: {irrep : [eigenvectors of invariant subspace; last entry: multiplicity of irrep occurrance]}
    c_dim = 0
    result = {}
    Rep = rep.copy("Rep")
    print("irreps of ", rep.name , ":")
    for irrep,chars in group.char_table.items():         
        temp = rep.copy("temp")
        P = projector_irrep(temp,chars)
        apply_projector(P,temp)
        if not projected_to_zero(temp):
            result[irrep] = []
            R = np.matrix(np.eye(len(P))) - P
            apply_projector(R,rep)                 
            if temp.is_reducible():
                print(irrep, "(reducible inv. subspace)-> occurrence: ", temp.characters["I"] / chars["I"])
                c_dim += temp.characters["I"]
            else:
                print(irrep, "-> occurrence: 1") 
                c_dim += temp.characters["I"]
            result[irrep] = nonzero_simult_eigvecs(Rep.hom["Rot2"],P)
            result[irrep].append(temp.characters["I"] / chars["I"])
    print("sum of dims: ",  c_dim)
    return result

def list_eigvecs(M,e):                  # Matrix M and eigenvalue e -> list of eigenvectors for this eigenvalue
    print("in list_eigvecs")
    result = []
    eigvals,eigvecs = np.linalg.eig(M)           
    for i in range(len(eigvals)):
        if abs(eigvals[i] - e) < 1e-8:
            # print("eigenvalue: ", eigvals[i], ". Eigenvector", eigvecs[:,i],  "appended")
            result.append(eigvecs[:,i])
    return result
def list_nonzero_eigvecs(M):
    print("in list_nonzero")
    result = []
    c = 0
    eigvals,eigvecs = np.linalg.eig(M)           
    for i in range(len(eigvals)):
        if not abs(eigvals[i]) < 1e-8:
            print("eigenvalue: ", eigvals[i], ". Eigenvector", eigvecs[:,i],  " appended")
            result.append(eigvecs[:,i])
            c += 1
    if not c:
        print("no nonzero EVs found")
    return result
def nonzero_simult_eigvecs(A,P):
    return list_nonzero_eigvecs(np.matmul(A,P))
             
