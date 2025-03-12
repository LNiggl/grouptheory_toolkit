import numpy as np
import classes as c
####################### GENERATION OF O_h #########################################

def remove_I(A):                    #remove access "I"s to make keys less ambiguous
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
    if "Rot" in O or "Inv" in O:
        while len(O)>0:
            if O[-1] == "v":
                Op = O[-3:]
                P.append(Op)
                O = O[:-3]
            else:
                Op = O[-4:]
                P.append(Op)
                O = O[:-4]
    else:
        while len(O)>0:
            n = 1
            while not O[-n].isalpha():
                n +=1
            Op = O[len(O)-n:]
            O = O[:-n]
            P.append(Op)   
    return P

def apply(A,vec,actions):  
    As = split(A)
    for i in range(len(As)):
        vec = apply_to_value(As[i],vec,actions)
    return vec

def apply_to_value(A,vec,actions):
    new_vec = vec.copy()
    for i in range(len(vec)):
        if "-" in actions[A]["{}".format(i)]:
            new_vec[int(actions[A]["{}".format(i)][1:])] = minus(vec[i])
        else:
            new_vec[int(actions[A]["{}".format(i)])] = vec[i]
    return new_vec

def find_closure(gen_actions, object):       #looks for closure of group ! regarding the provided action and applied object(usually element of some vector space)
    actions = gen_actions
    known_objects = {"I": object}       
    #plan: apply actions to all objects in known objects, if result is new, add to known_objects with name for the group element
    n = 0
    while n != len(known_objects):                          
        n = len(known_objects)
        for a in actions.copy().keys():
            for key, o in known_objects.copy().items():                
                new_o = apply(a,o,actions)
                if new_o not in known_objects.values():
                    new_name = remove_I(a + key)
                    known_objects[new_name] = new_o
                    actions[new_name] = readout_action(new_o,gen_actions)      
    return actions

def minus(name):                #makes string for mult with -1. ONLY works if "-" is first character
    if name[0] == '-':
        return name[1:]
    return '-' + name

def readout_action(vec,actions):            #gives according row in actions dict for a vector v resulting from some group action: v = gb (b "basis")
    r = {}
    for i in range(len(vec)):
        r['{}'.format(vec[i])] = '{}'.format(i)
    if len(vec) != len(actions["I"]):
        for i in range(len(vec)): 
            r[minus('{}'.format(vec[i]))] = minus('{}'.format(i))
    return r

def mult_table_from_actions(actions):                           #CAREFUL! Acting on values with A is like acting on keys with A^-1. The generating actions are self-invers, however (AB)^-1 = B^-1 A^-1  
    mt = {}
    for a in actions.keys():                             
        mt[a] = {}   
        for b in actions.keys():
            A = remove_I(a+b)
            w = apply(A,["0","1","2"],actions)                              #result is w = Ae = abe for reference vector e used for all actions, here e = ["0","1","2"]
            w = readout_action(w,actions)                                           #convert to long format {0: .. , 1: .., 2: .., -0: .., -1: .., -2: ..}            
            c = list(actions.keys())[list(actions.values()).index(w)]               # look for c such that ce = abe -> key of w      
            mt[a][b] = c                                                            # a.b = c because they have the same action
    return mt

def find_conjugacy_classes(elements,mt,inv):
    n=len(elements)
    tostudy=elements.copy()                 #copy is important, otherwhise elements gets manipulated when tostudy is manipulated!
    classes = {}
    while len(tostudy) > 0:
        x = tostudy[0]
        conjugate_to = []
        for i in range(len(elements)):
            conjugate_to.append(mt[mt[elements[i]][x]][inv[elements[i]]])
        conjugate_to = list(set(conjugate_to))
        for y in conjugate_to:
            tostudy.remove(y)
        classes[x] = conjugate_to
    return classes

def set_group_from_actions(gen_actions,obj):
    actions = find_closure(gen_actions,obj)
    el = list(actions.keys())
    mt = mult_table_from_actions(actions)
    G = c.Group(el,mt)
    return G

##################### REPRESENTATIONS ###########################

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
    r = c.Representation(A,hom,A.name + "_x_" + B.name)
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

################## DEALING WITH PARTICLES ###########################

################## PIONS ############################################

def generate_basis(pions,group):        
    O = []
    for p in pions:
        for q in p.orbit(group):
            O.append(q)
    basis = order_basis(O)
    return basis    

def order_basis(pions):          #sorts by: magnitude of momenta in ascending order
                                #           direction within same overall momentum 
    p = remove_neg(pions)
    sort_by_length(p)
    sort_by_direction(p)
    return p

def remove_neg(pions):
    new_list = []
    for p in pions:
        if p.sign == "+":
            new_list.append(p)
    return new_list

def sort_by_length(pions):
    print("sort by length")
    for i in range(len(pions)): 
        pions[i].length = np.sqrt(np.sum([x**2 for x in pions[i].momentum]))
    pions.sort(key = lambda x: x.length)

def sort_by_direction(pions):                               # needs improvement
    decorated = [(p.momentum[0],p.momentum[1],p.momentum[2],p) for p in pions]
    decorated.sort()
    pions = [p for (x,y,z,p) in decorated]

def str_sign(sign_mat):
    assert sign_mat == -1 or sign_mat == 1
    if sign_mat == 1:
        return "+"
    return "-"

def int_sign(sign):
    assert sign == "+" or sign == "-"
    if sign == "+":
        return 1
    return -1