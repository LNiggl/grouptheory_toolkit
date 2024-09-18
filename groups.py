import numpy as np
class Group:                            #includes all important things we know about the group 
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
        while True:
            x,y = np.random.choice(len(e),2)
            if x != y:
                break
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
        s = Group(subgroup,mult_table)
        return s
    
###################### tools for handling the names of group elements ############    

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

####################### GENERATION OF O_h ##########################

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

######################## group functions #########################

def inner_product(f,h,group,classfunction = True):           # takes (class) functions, e.g. characters, as dictionaries and returns the inner product 1/|G| * \sum_g (f(g)h*(g)) (see e.g. lecture ch. 11)
    s = 0
    if classfunction:
        for c,r in group.classes.items():
            s += f[c]*np.conjugate(h[c])*len(r)
    else:
        for e in group.elements:
            s += f[e]*np.conjugate(h[e])
    return s / len(group.elements)

def is_group_homomorphism(map,G,H):                     # Groups G,H, dictionary map: G -> H; returns true if for all elements of G, map(g1g2) = map(g1)map(g2)
    for g1 in G.elements:
        for g2 in G.elements:
            if map[G.mult_table[g1][g2]] != H.mult_table[map[g1]][map[g2]]:
                return False
    return True


