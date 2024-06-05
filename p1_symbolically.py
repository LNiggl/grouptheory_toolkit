import numpy as np
class group:                            #includes all important things we know about the group 
    # elements = []                       #list of strings
    # mult_table = {}                     #dict of dicts: {g1: {all elements gi: result of g1gi}}
    # inverse_list = {}                   #dict: {g : g^-1} 
    # classes = {}                        #dict: {representative: list of class members}
    # char_table                        #{(dict) irrep: { (dict) class: character}}
                                        #NB: save action as dict of dicts: {g:{b_i : b_j,...}, h:...}

    def __init__(self,e,m_table):
        self.elements = e.copy()        
        self.mult_table = m_table.copy()
        self.inverse_list = self.inverses()
        self.classes = self.find_conjugacy_classes()
    
    def find_identity(self):
        e = self.elements.copy()
        mt = self.mult_table.copy()
        x,y = np.random.choice(len(e),2)                                         
        Ex = list(mt[e[x]].keys())[list(mt[e[x]].values()).index(e[x])]            
        Ey = list(mt[e[y]].keys())[list(mt[e[y]].values()).index(e[y])]
        assert Ex == Ey                                                 # I is such that gI=I for any g; evaluate from two g randomly
        if Ex != "I":
            if not any(self.elements == "I"):
                rename_key(self.elements,Ex,"I") 
            else:
                print("Name of identity: ", str(Ex))
        return Ex

    def inverses(self):
        mt = self.mult_table
        inv = {}
        I = self.find_identity()
        for g in mt.keys():
            g_inv = list(mt[g].keys())[list(mt[g].values()).index(I)]
            inv[g] = g_inv
        return inv
    
    def find_conjugacy_classes(self):
        e = self.elements
        mt = self.mult_table
        inv = self.inverse_list
        n=len(e)
        tostudy=e.copy()                 #copy is important, otherwhise elements gets manipulated when tostudy is manipulated!
        classes = {}
        while len(tostudy) > 0:
            x = tostudy[0]
            conjugate_to = []
            for i in range(len(e)):
                conjugate_to.append(mt[mt[e[i]][x]][inv[e[i]]])
            conjugate_to = list(set(conjugate_to))
            for y in conjugate_to:
                tostudy.remove(y)
            classes[x] = conjugate_to
        return classes
    
    def set_char_table(self,irreps):
        self.char_table = {}
        for r in irreps:
            self.char_table[r.name] = r.characters

    # def print_char_table(self):
    #     # print(["{}".format(c)+"\t" for c in self.classes])
    #     s = ""
    #     t = ""
    #     for i in range(len(self.classes)-2):
    #         s += "%-10s"
    #     print(s)
    #     print(s % [i for i in self.classes.keys()])
            
    def find_subgroup(self, gen_el):                            #generates the subgroup via the Cayley table
        subgroup = gen_el.copy()
        n = 0        
        while n < len(subgroup):  
            n = len(subgroup)   
            temp = subgroup.copy()          
            for a in temp:
                for b in temp:
                    subgroup.append(self.mult_table[a][b])
                    subgroup.append(self.inverse_list[self.mult_table[a][b]])
            subgroup = list(set(subgroup))
        mult_table = {}
        for s in subgroup:
            mult_table[s] = self.mult_table[s].copy()
        s = group(subgroup,mult_table)
        return s          
class representation(group):
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

####################### Generation of O_h #########################################
######### functions #####
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
    # print("in removeI: ", A , " -> ", B)
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
    # print("in apply:")    
    As = split(A)
    # print("A:", A, " -> " , As)
    for i in range(len(As)):
        vec = apply_to_value(As[i],vec,actions)
    # print("[0,1,2] -> ",vec)
    # print("apply done")
    return vec

def apply_to_value(A,vec,actions):
    # print("in apply_to_value:")
    new_vec = vec.copy()
    # print(vec)
    for i in range(len(vec)):
        if "-" in actions[A]["{}".format(i)]:
            new_vec[int(actions[A]["{}".format(i)][1:])] = minus(vec[i])
        else:
            new_vec[int(actions[A]["{}".format(i)])] = vec[i]
        # print("step ", i , vec, "->" , new_vec)
    # print(A, " " , vec, " = ",  new_vec)
    return new_vec

## outdated: old apply() and the applications to the keys; ERROR: application to the key yields action of g^-1 
# def apply(A,vec,actions):
#     n = len(A)
#     namescheme = "rot"
#     while A[-1] == "I" and n > 1:
#             A = A[:-1]
#             n = len(A)
#     if n == 1:
#         # print("action ",A)
#         return apply_any_to_key(A,vec,actions)
#     # if "I" in A:
#         # print("action ",A)
#     if A[0].isalpha():
#         if not A[1].isalpha():            
#             namescheme = "axes"
#             # print(namescheme)
#     if namescheme == "axes":
#         return apply_any_to_key(A,vec,actions)
#     if namescheme == "rot": 
#         # print("exit rot")
#         return apply_to_value(A,vec,actions)
# def apply_gen_to_key(A,vec,actions):               #applies matrix A on vector vec by changing the keys: A_1A_2vec like matrix mult acting on vec
#                                             #CAREFUL! USE ONLY FOR GENERATING ACTIONS. Acting on values with A is like acting on keys with A^-1. The generating actions are self-inverse, however (AB)^-1 = B^-1 A^-1
#     temp = {"{}".format(i) : vec[i] for i in range(len(vec))}
#     for k in temp.keys():                               #apply action on keys instead of vector values
#         rename_key(temp,k,k+"->"+actions[A][k])
#     for k,v in temp.items():                            #properly rename the keys and account for possible minus
#         if k[-2] == '-':  
#             if v[0] == "-":
#                 temp[k] = v[1:]         
#             else: 
#                 temp[k] = "-"+v 
#         rename_key(temp,k,k[-1])
#     # return [temp["0"],temp["1"],temp["2"]]
#     return [temp["{}".format(i)] for i in range(len(vec))]
# def apply_any_to_key(A,vec,actions):
#     As = split(A)
#     for a in As:
#         vec = apply_gen_to_key(a,vec,actions)
#     return vec

def find_closure(gen_actions, object):       #looks for closure of group ! regarding the provided action and applied object(usually element of some vector space) ! 
                                            #CAREFUL! Acting on values with A is like acting on keys with A^-1. Axis exchange and sign change are self-inverse, 
                                            #however (AB)^-1 = B^-1 A^-1 
                                            #maybe make more elegantly with the presentations of groups
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
                    # if new_name[-1] == 'I':
                    #     new_name = new_name[:-1]                #action on key is in inverse order compared to action on vals
                    known_objects[new_name] = new_o
                    # actions[new_name] = {"0": new_o[0], "1" : new_o[1], "2": new_o[2], "-0": minus(new_o[0]),"-1" : minus(new_o[1]), "-2" : minus(new_o[2])} 
                    actions[new_name] = readout_action(new_o,gen_actions)      
    return actions

def minus(name):                #makes string for mult with -1. ONLY works if "-" is first character
    if name[0] == '-':
        return name[1:]
    return '-' + name

def readout_action(vec,actions):            #gives according row in actions dict for a vector v resulting from some group action: v = gb (b "basis")
    r = {}
    for i in range(len(vec)):
        # r['{}'.format(i)] = '{}'.format(vec[i])
        r['{}'.format(vec[i])] = '{}'.format(i)
    if len(vec) != len(actions["I"]):
        for i in range(len(vec)):
            # r['-{}'.format(i)] = minus('{}'.format(vec[i])) 
            r[minus('{}'.format(vec[i]))] = minus('{}'.format(i))
    # print("in readout_action:", r)
    return r

#backup
# def mult_table_from_actions(actions):                           #CAREFUL! Acting on values with A is like acting on keys with A^-1. The generating actions are self-invers, however (AB)^-1 = B^-1 A^-1  
#     mt = {}
#     for a in actions.keys():                             
#         mt[a] = {}   
#         for b in actions.keys():
#             if a != "I" and b != "I":
#                 A = a+b
#             else:
#                 if a == "I":
#                     A = b
#                 if b == "I":
#                     A = b            
#             w = apply(A,["0","1","2"],actions)                              #result is w = Ae = abe for reference vector e used for all actions, here e = ["0","1","2"]
#             w = readout_action(w,actions)                                           #convert to long format {0: .. , 1: .., 2: .., -0: .., -1: .., -2: ..}
#             c = list(actions.keys())[list(actions.values()).index(w)]               # look for c such that ce = abe -> key of w      
#             mt[a][b] = c                                                            # a.b = c because they have the same action
#     return mt

def mult_table_from_actions(actions):                           #CAREFUL! Acting on values with A is like acting on keys with A^-1. The generating actions are self-invers, however (AB)^-1 = B^-1 A^-1  
    mt = {}
    
    for a in actions.keys():
        # print("action a: ", a)                             
        mt[a] = {}   
        for b in actions.keys():
            # print("action b: ", b)
            e = ["0","1","2"]
            A = remove_I(a+b)
            # v = apply(b,e,actions)
            # w = apply(a,v,actions)
            # print(a,b, ": e -> ", v , " -> " , w)
            w = apply(A,["0","1","2"],actions)                              #result is w = Ae = abe for reference vector e used for all actions, here e = ["0","1","2"]
            w = readout_action(w,actions)
           
                                                       #convert to long format {0: .. , 1: .., 2: .., -0: .., -1: .., -2: ..}            
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

########## representations #####################
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
#backup
# def matrix_from_action(A,basis):
#     d = len(basis)
#     M = np.zeros((d,d))
#     if isinstance(A,dict):                                                  # for when A is provided as dict 
#         if not any(char.isalpha() for a in A.keys() for char in a):         # and basis contains only numbers as labels   
#             for i in range(d):
#                 for j in range(d):
#                     if "{}".format(i) in A["{}".format(j)]:
#                         if "-" in A["{}".format(j)]:
#                             M[j][i] = -1
#                         else: 
#                             M[j][i] = 1          
#     else:
#         if callable(getattr(basis[0],"action", None)):                      # for when the basis has a class method .action(A)
#             acted_basis = []
#             for b in basis:
#                 x = b.copy()
#                 x.action(A)
#                 acted_basis.append(x)
#             for i in range(d):
#                 for j in range(d):
#                     if basis[j].is_equal_to(acted_basis[i]):               # "if e_i gets transformed to e_j via action A" -> ge_i = M(g)_{ji} e_j == e_j
#                         M[j][i] = 1                                     # --> M(A)_{ji} == 1
#                     if basis[j].is_negative_of(acted_basis[i]):
#                         M[j][i] = -1                
#     return M

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
    r = representation(A,hom,A.name + "_x_" + B.name)
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
        
############## elements that span group
#natural basis upon which O_h acts:
# vectors = ["0","1","2","-0","-1","-2"]
v = ["0","1","2"]
gen_actions_O_h = {
            "I" : {"0" : "0", "1" : "1", "2": "2", "-0" : "-0", "-1" : "-1", "-2" : "-2"},          #identity

           "A01": {"0" : "1", "1" : "0", "2": "2", "-0" : "-1", "-1" : "-0", "-2" : "-2"},          #axis interchange
           "A02": {"0" : "2", "1" : "1", "2": "0", "-0" : "-2", "-1" : "-1", "-2" : "-0"},
           "A12": {"0" : "0", "1" : "2", "2": "1", "-0" : "-0", "-1" : "-2", "-2" : "-1"},

           "R0": {"0" : "-0", "1" : "1", "2": "2", "-0" : "0", "-1" : "-1", "-2" : "-2"},           #axis sign change
           "R1": {"0" : "0", "1" : "-1", "2": "2", "-0" : "-0", "-1" : "1", "-2" : "-2"},
           "R2": {"0" : "0", "1" : "1", "2": "-2", "-0" : "-0", "-1" : "-1", "-2" : "2"},
           }
gen_actions_O_h_rot = {
            "I" : {"0" : "0", "1" : "1", "2": "2", "-0" : "-0", "-1" : "-1", "-2" : "-2"},          #identity

            "Rot0": {"0" : "0", "1" : "2", "2": "-1", "-0" : "-0", "-1" : "-2", "-2" : "1"},        #rotations
            "Rot1": {"0" : "-2", "1" : "1", "2": "0", "-0" : "2", "-1" : "-1", "-2" : "-0"},
            "Rot2": {"0" : "1", "1" : "-0", "2": "2", "-0" : "-1", "-1" : "0", "-2" : "-2"},

            "Inv": {"0" : "-0", "1" : "-1", "2": "-2", "-0" : "0", "-1" : "1", "-2" : "2"},         #inversion
}
################################################
# TESTS
# v = ["0","1","2"]
# for A in gen_actions_O_h_rot:
#     vec = v.copy()
#     apply(A,vec,gen_actions_O_h_rot)            #works as intended

# print("as actions:")
# print("1st:",readout_action(vv2,gen_actions_O_h_rot))
# print("2nd:",readout_action(vvv2,gen_actions_O_h))
# v = ["0","1","2"]
# rotations = find_closure(gen_actions_O_h_rot,v)
# print(rotations.keys())
#rot way
# s = split("Rot1Rot0")
# print("splitting of s: ", s)
# print(rotations)
# print(len(rotations))
# elr = list(rotations.keys())
# mtr = mult_table_from_actions(rotations)
# O_h_r = group(elr,mtr)
#previous way
# O_h_actions = find_closure(gen_actions_O_h,v)
# el = list(O_h_actions.keys())
# print("elements old fashioned way:", el)
# print("#",len(el))
# mt = mult_table_from_actions(O_h_actions)
# O_h = group(el,mt)

# v = [["0","1","2"],["1","2","0"],["2","0","1"],["0","2","1"],["1","0","2"],["2","1","0"],["-0","-1","-2"],["-1","-2","-0"],["-2","-0","-1"],["-0","-2","-1"],["-1","-0","-2"],["-2","-1","-0"],
#     ["0","-1","-2"],["1","-2","-0"],["2","-0","-1"],["0","-2","-1"],["1","-0","-2"],["2","-1","-0"]]
# for i in range(len(v)):
#     apply_to_value("Rot1",v[i],gen_actions_O_h_rot)
# print("test find_closure")
# print("rot:")
# r = find_closure(gen_actions_O_h_rot,v)
# print(r.keys())
# print("#",len(r))
# CTr = mult_table_from_actions(r)
# test = "Rot1Rot0"
# test = "I"
# print("table for ", test, ": " , CTr[test])
# print("#elements: ",len(r))
# print(list(r.keys()) == list(CTr["I"].keys()))

# print("ax:")
# a = find_closure(gen_actions_O_h,v)
# print(a.keys())
# print(len(a.keys()))
# CTa = mult_table_from_actions(a)
# test = "R1A01"
# # test = "I"
# print("table for ", test, ": " , CTa[test])
# assert False




##################################################


### Axes formalism ########
# print("axes approach")
# ########### full cubic group O_h #########################
# print("Full cubic group O_h")
# O_h_actions = find_closure(gen_actions_O_h,v)
# el = list(O_h_actions.keys())
# mt = mult_table_from_actions(O_h_actions)
# O_h = group(el,mt)
# ### irreps of O_h  ######
# A1p = representation(O_h,rep_trivial(O_h),"A1p")
# A1p.check_if_homomorphism()

# #test 
# # m = matrix_from_action("Rot2",v)
# # print("R1A01: ", m)

# Rep_T1m = rep_from_action(O_h,v,O_h_actions)
# print("R1A12: ", Rep_T1m["R1A12"])
# T1m = representation(O_h,Rep_T1m,"T1m")
# T1m.check_if_homomorphism()

# A1m = representation(O_h,rep_determinant(T1m.hom),"A1m")
# A1m.check_if_homomorphism()

# T1m_x_T1m =  product_rep(T1m,T1m)           #no irrep
# T1m_T1m_reps = T1m_x_T1m.hom.copy()
# T1m_x_T1m.check_if_homomorphism()

# T1p = representation(O_h,T1m_T1m_reps,"T1p")
# apply_projectors([antisymmetric_projector],T1p)
# T1p.check_if_homomorphism()

# T2p = representation(O_h,T1m_T1m_reps,"T2p")
# apply_projectors([symmetric_projector,invert_projector(diagonal_projector)],T2p)
# T2p.check_if_homomorphism()

# T2m = product_rep(T2p,A1m)
# T2m.name = "T2m"
# T2m.check_if_homomorphism()

# Ep = representation(O_h,T1m_T1m_reps,"Ep")
# apply_projectors([traceless_projector,diagonal_projector],Ep)
# Ep.check_if_homomorphism()
# Ep.round_chars()

# Em = product_rep(Ep,A1m)
# Em.name = "Em"
# Em.check_if_homomorphism()
# Em.round_chars()

# A2m = product_rep(T1m,product_rep(T1m,T1m))
# A2m.name = "A2m"
# project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
# A2m.check_if_homomorphism()
# A2m.round_chars()

# A2p = product_rep(A2m,A1m)
# A2p.name = "A2p"
# A2p.check_if_homomorphism()
# A2p.round_chars()

# O_h.set_char_table([A1m,A1p,T1m,T1p,T2m,T2p,Em,Ep,A2m,A2p])


# # cubic group O
# print("Cubic group O")
# gen_O = ["I","R1A12","R0A02","R1A01"]
# O = O_h.find_subgroup(gen_O)
# . . .

############### Rotations formalism #############
# print("rotational approach")
############### Full Cubic Group ################
print("Full cubic group O_h")
O_h_actions_r = find_closure(gen_actions_O_h_rot,v)
el_r = list(O_h_actions_r.keys())
# print(el_r)
# print("# ",len(el_r))
mt_r = mult_table_from_actions(O_h_actions_r)
# print(mt_r["InvRot2"])
O_h_r = group(el_r,mt_r)
### irreps of O_h  ######
A1p = representation(O_h_r,rep_trivial(O_h_r),"A1p")
A1p.check_if_homomorphism()

#test 
# M = matrix_from_action("Rot2",v)
# print("Rotation around z: ", M)

Rep_T1m = rep_from_action(O_h_r,v,O_h_actions_r)
# print("Rotation around x: ",Rep_T1m["Rot0"])
T1m = representation(O_h_r,Rep_T1m,"T1m")
T1m.check_if_homomorphism()

A1m = representation(O_h_r,rep_determinant(T1m.hom),"A1m")
A1m.check_if_homomorphism()

T1m_x_T1m =  product_rep(T1m,T1m)           #no irrep
T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

T1p = representation(O_h_r,T1m_T1m_reps,"T1p")
apply_projectors([antisymmetric_projector],T1p)
T1p.check_if_homomorphism()

T2p = representation(O_h_r,T1m_T1m_reps,"T2p")
apply_projectors([symmetric_projector,invert_projector(diagonal_projector)],T2p)
T2p.check_if_homomorphism()

T2m = product_rep(T2p,A1m)
T2m.name = "T2m"
T2m.check_if_homomorphism()

Ep = representation(O_h_r,T1m_T1m_reps,"Ep")
apply_projectors([traceless_projector,diagonal_projector],Ep)
Ep.check_if_homomorphism()
Ep.round_chars()

Em = product_rep(Ep,A1m)
Em.name = "Em"
Em.check_if_homomorphism()
Em.round_chars()

A2m = product_rep(T1m,product_rep(T1m,T1m))
A2m.name = "A2m"
project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
A2m.check_if_homomorphism()
A2m.round_chars()

A2p = product_rep(A2m,A1m)
A2p.name = "A2p"
A2p.check_if_homomorphism()
A2p.round_chars()

O_h_r.set_char_table([A1m,A1p,T1m,T1p,T2m,T2p,Em,Ep,A2m,A2p])


# cubic group O
# print("Cubic group O")
gen_O_r = ["I","Rot0","Rot1","Rot2"]
O_r = O_h_r.find_subgroup(gen_O_r)
# . . .


############### PARTICLES ######################
############### Single pions ####################

class pion:
    def __init__(self,momentum,sign):                  # important for trafo is the momentum vector and overall sign of pion operator        
        if sign != "+" and sign != "-":
            sign = str_sign(sign)
        assert sign == "+" or sign == "-"
        self.momentum = momentum
        self.sign = sign
        self.update_name()

    def update_name(self):
        self.name = self.sign + "pi({:.0f},{:.0f},{:.0f})".format(self.momentum[0],self.momentum[1],self.momentum[2])

    def copy(self):
        momentum = self.momentum.copy()
        sign = str(self.sign)
        new_p = pion(momentum, sign)
        return new_p
    
    def action(self,A):                         # make sure to apply inverse of A and det(A); maybe different archetype 
        sgn_mat = A1m.hom[A]*int_sign(self.sign)
        self.sign = str_sign(sgn_mat)
        self.momentum = np.matmul(self.momentum,T1m.hom[T1m.inverse_list[A]])           # Watch out: transforming under rep means right multiplication with matrix
        self.update_name()

    def is_equal_to(self,pion):
        if self.sign != pion.sign:
            return False
        if not np.allclose(self.momentum, pion.momentum,rtol = 1e-5, atol = 1e-10):
            return False
        return True
    
    def change_sign(self):
        if self.sign == "+":
            self.sign = "-"
        else: 
            self.sign = "+"
        self.update_name()

    def invert_momentum(self):
        for i in range(len(self.momentum)):
            self.momentum[i] = -1*self.momentum[i]
        self.update_name()

    def is_negative_of(self,pion):
        self.change_sign()
        if self.is_equal_to(pion):
            self.change_sign()
            return True
        self.change_sign()
        return False
    
    def orbit(self,group):
        orbit = [self.copy()]
        for g in group.elements:
            p = self.copy()
            p.action(g)
            if not any([q.is_equal_to(p) for q in orbit]):
                orbit.append(p)                
        return orbit
    
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

class twopion:
    def __init__(self,pi1,pi2,Momentum = [0,0,0]):                #attributes of TwoP system with capital letters
        assert isinstance(pi1,pion)
        assert isinstance(pi2,pion)
        assert all(Momentum[i]-pi1.momentum[i]-pi2.momentum[i] < 1e-10 for i in range(len(Momentum)))
        assert pi1.sign == pi2.sign
        self.pi1 = pi1.copy()
        self.pi2 = pi2.copy()
        self.Momentum = Momentum
        self.overall_sign()             # sets self.Sign
        self.update_name()  

    def update_name(self):
        self.name = self.Sign + "(" + self.pi1.name[1:] + "_x_" +  self.pi2.name[1:] + ")"

    def copy(self):
        p1 = self.pi1.copy()
        p2 = self.pi2.copy()
        M = self.Momentum.copy()
        new_twop = twopion(p1,p2,M)
        return new_twop
    
    def overall_sign(self):
        s1 = int_sign(self.pi1.sign)
        s2 = int_sign(self.pi2.sign)
        self.Sign = str_sign(s1*s2)

    def action(self,A):
        self.pi1.action(A)
        self.pi2.action(A)
        self.overall_sign()
        self.update_name()

    def is_equal_to(self, twopion):
        if not (np.allclose(self.pi1.momentum, twopion.pi1.momentum,rtol = 1e-5, atol = 1e-10) \
                and np.allclose(self.pi2.momentum, twopion.pi2.momentum,rtol = 1e-5, atol = 1e-10) and (self.Sign == twopion.Sign)):
            return False
        return True
    
    def is_negative_of(self,twopion):
        self.change_sign()
        if self.is_equal_to(twopion):
            self.change_sign()
            return True
        self.change_sign()
        return False
    
    def change_sign(self):
        if self.Sign == "+":
            self.Sign = "-"
        else: 
            self.Sign = "+"
        self.update_name()

    def orbit(self,group):
        orbit = [self.copy()]
        for g in group.elements:
            p = self.copy()
            p.action(g)
            if not any([q.is_equal_to(p) for q in orbit]):
                orbit.append(p)                
        return orbit
    
############## tests ############

P1 = pion([1,0,0],"+")
P2 = pion([-1,0,0],"+")
# P3 = pion([1,1,0],"+")
# P4 = pion([2,0,0],"-")
P5 = pion([2,1,0],"+")
# P6 = pion([2,0,1],"+")
# P7 = pion([1,0,1],"-")
# P8 = pion([0,1,1],"-")
# M = [P4,P1,P2,P3,P5,P6,P7,P8]

#single Pion
# print(P1.name)
# P1.action("Rot2")
# print("P1 after rotation: ", P1.name)
# P1.action(O_h_r.inverse_list["Rot2"])
# print("rotated back: ", P1.name)

# #orbits:
# o1 = P1.orbit(O_h_r)
# print("orbit of ",P1.name, ": ")
# for p in o1:
#     print(p.name)
# o2 = P2.orbit(O_h_r)
# print("orbit of ",P2.name, ": ")
# for p in o1:
#     print(p.name)

# #basis from orbit
# B1 = generate_basis([P1],O_h_r)
# print(len(B1))
# print("B1 generated by", P1.name, ":")
# for b in B1: 
#     print(b.name)

# A = "Rot2"
# print("Rep of Two Pi under rotation:", A )
# mat = matrix_from_action("Rot2",B1)                # todo: test also for other bases
# print(mat)

# pion_rep = rep_from_action(O_h_r,B1)
# SingleP = representation(O_h_r,pion_rep,"SingleP")
# SingleP.check_if_homomorphism()

# P1 = pion([1,0,0],"+")
# P2 = pion([-1,0,0],"+")
# # P_total = []
# # for i in range(len(P1.momentum)):
# #     P_total.append(P1.momentum[i]+P2.momentum[i])
# # print(P_total)
# print("Two Pion")
# tp = twopion(P1,P2)
# # tp.action("R0A01")
# print(tp.name)

# print("orbit")
# tporbit = tp.orbit(O_h_r)
# for p in tporbit:
#     print(p.name)
# print("#: ", len(tporbit))
# # here: orbit = basis
# m = matrix_from_action(A,tporbit)
# print("Rep of Two Pi under rotation:", A)
# print(m)
# tp_rep = rep_from_action(O_h_r,tporbit)
# TwoP = representation(O_h_r,tp_rep,"TwoP")
# TwoP.check_if_homomorphism()

# for g in gen_actions_O_h_rot:
#     print("inverse action of ", g , "on [a,b,c]: ")
#     print(apply(O_h_r.inverse_list[g],["a","b","c"],O_h_actions_r))

print("char table of O_h, rotational formalism: ")
for c in O_h_r.char_table.keys():
    print(c, ": ",O_h_r.char_table[c])

# print("char table of O_h, axes formalism: ")
# for c in O_h.char_table.keys():
#     print(c, " ",O_h.char_table[c])












###### Last work #########

# work on rotational basis formalism
# division into smaller files

####### FIX ########



### ADDITIONS ###

# add function/group method print_char_table
# add group method .isgroup() 
# latex document with all representation matrices for single and two pion trafo
# create python library

#isomorphism bewteen axes and rotation formalism

### NEXT PLAN ###

# method for extracting irreps







