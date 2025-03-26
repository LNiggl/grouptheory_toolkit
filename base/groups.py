import numpy as np

from . import tools as t

class Group:                                                        # used to access: elements, Cayley table, classes, character table 
    def __init__(self,e,m_table,inv_list = None,class_list = None):
        self.elements = e.copy()  
        mt = m_table.copy()     
        self.mult_table = mt  
        self.find_identity()
        if inv_list == None:    
            self.inverse_list = self.inverses()
        else: 
            self.inverse_list = inv_list
        if class_list == None:
            self.classes = self.find_conjugacy_classes()
        else: 
            self.classes = class_list
    
    def find_identity(self):
        e = self.elements.copy()
        mt = self.mult_table.copy()
        while True:
            x,y = np.random.choice(len(e),2)
            if x != y:
                break
        Ex = list(mt[e[x]].keys())[list(mt[e[x]].values()).index(e[x])]            
        Ey = list(mt[e[y]].keys())[list(mt[e[y]].values()).index(e[y])]
        assert Ex == Ey                                                     # I is such that gI=I for any g; evaluate from two g randomly and compare
        if Ex != "I":
            if not any(self.elements == "I"):
                print("Name of identity: ", str(Ex))
                print(str(Ex), " is renamed to I")
                t.rename_key(self.elements,Ex,"I") 
            else:
                raise ValueError("Identity is not named 'I', and name 'I' exists for different element.")                
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
        n = len(e)
        tostudy = e.copy()                 
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
    
    def set_char_table(self,R,load_dict = False):               # R: either list of Representation objects, or nested dict in form of the final char table (used in load_group)
        if load_dict:
            self.char_table = R
        else:
            self.char_table = {}
            for r in R:
                self.char_table[r.name] = r.characters
                
    def find_subgroup(self, gen_el):                            # gen_el: generating set of subgroup. returns Group obj of subgroup (if closure is one)
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


def inner_product(f,h,group,classfunction = True):           # takes (class) functions, e.g. character(c), as dictionaries and returns the inner product 1/|G| * \sum_g (f(g)h*(g)) (see e.g. lecture ch. 11)
    s = 0
    if classfunction:
        for c,r in group.classes.items():
            s += f[c]*np.conjugate(h[c])*len(r)
    else:
        for e in group.elements:
            s += f[e]*np.conjugate(h[e])
    return s / len(group.elements)

def is_group_homomorphism(map,G,H):                     # Groups G,H, dictionary map: G -> H; returns true if for all elements of G, map(g1g2) = map(g1)map(g2)
                                                        # map could be obtained via representations.group_hom_via_rep(..)
    for g1 in G.elements:
        for g2 in G.elements:
            if map[G.mult_table[g1][g2]] != H.mult_table[map[g1]][map[g2]]:
                return False
    return True



