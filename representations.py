import numpy as np
import groups as g
import objects as o
class Representation(g.Group):
    # hom = {}                           #dict for homomorphism: G -> GL(n)
    # characters = {}
    def __init__(self,group,dict_matrix,name):             #initialize with dict: g -> M(g)
        self.elements = group.elements
        self.classes = group.classes
        self.mult_table = group.mult_table
        self.inverse_list = group.inverse_list

        self.name = name
        self.hom = dict_matrix.copy()
        # self.matrices = list(self.hom.values())
        self.characters = {}
        for c in self.classes:
            self.characters[c] = np.trace(self.hom[c])

    def check_if_homomorphism(self):
        for g in self.elements:
            for h in self.elements:
                eps = np.linalg.norm(np.matmul(self.hom[g],self.hom[h]) - self.hom[self.mult_table[g][h]])
                assert eps < 1e-13

    def update(self):
        for c in self.classes:
            self.characters[c] = np.trace(self.hom[c])

    def round_chars(self):
        for c,chi in self.characters.items():
            if chi - round(chi) < 1e-12:
                self.characters[c] = round(chi)

def matrix_from_action(A,basis):
    d = len(basis)
    M = np.zeros((d,d))
    if isinstance(A,dict):                                                  # for when A is provided as dict 
        if not any(char.isalpha() for a in A.keys() for char in a):         # and basis contains only numbers as labels   
            for i in range(d):
                for j in range(d):
                    if "{}".format(i) in A["{}".format(j)]:
                        if "-" in A["{}".format(j)]:
                            M[i][j] = -1
                        else: 
                            M[i][j] = 1          
    else:
        if callable(getattr(basis[0],"action", None)):                      # for when the basis has a class method .action(A)
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

def rep_from_action(group,basis,actions = None):
    if actions == None:
        return {g : matrix_from_action(g,basis) for g in group.elements}
    else:
        if isinstance(actions,dict):
            return {g : matrix_from_action(actions[g],basis) for g in group.elements}       #returns dict of the homomorphism
        
def rep_trivial(group):
    return {g: np.matrix([1]) for g in group.elements}                  #returns dict of the homomorphism

def rep_determinant(rep):                                               #pass dict of homomorphism, returns new dict of homomorphism
    return {g: np.matrix(np.linalg.det(rep[g])) for g in rep.keys()}

def product_rep(A, B):                                  #pass two  representation objects (of same group). Returns new representation object
    hom = {}
    for g in A.elements:
       hom[g] = np.kron(A.hom[g],B.hom[g])
    r = Representation(A,hom,A.name + "_x_" + B.name)
    return r            

# projectors (adapted from the lecture notes: Prof. Dr. Christoph Lehner: Lattice QCD, chapter 11)
def projector_irrep(rep,irrep):                                  
    ret = 0.0 * rep.hom["I"]
    for c in rep.classes.keys():
        for g in rep.classes[c]:
            ret += rep.hom[g]*irrep.characters[c]
    P = np.matrix(ret * (irrep.characters["I"] / len(rep.elements)))
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
    P = projectors[0](n)
    for p in projectors[1:]:
        P = P * p(n)
    for g,v in rep.hom.items():
        rep.hom[g] = P*rep.hom[g]*P
    rep.update()                        #updates characters 

def cut_empty_rows(rep):            #used in minimize_dimension()
    for g in rep.hom.keys():
        i = 0
        while i < len(rep.hom[g]):    
            if not rep.hom[g][i].any():                
                rep.hom[g] = np.delete(rep.hom[g],i,0)
            i += 1
    return rep

def minimize_dimension(rep):    #cut empty lines, then transpose and cut empty lines again, then transpose again
    cut_empty_rows(rep)
    for g,m in rep.hom.items():
        rep.hom[g]=rep.hom[g].T
    cut_empty_rows(rep)
    for g,m in rep.hom.items():
        rep.hom[g]=rep.hom[g].T                
