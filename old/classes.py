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
        x,y = np.random.choice(len(e),2)                                         
        Ex = list(mt[e[x]].keys())[list(mt[e[x]].values()).index(e[x])]            
        Ey = list(mt[e[y]].keys())[list(mt[e[y]].values()).index(e[y])]
        assert Ex == Ey                                                 # I is such that gI=I for any g; evaluate from two g randomly
        if Ex != "I":
            if not any(self.elements == "I"):
                f.rename_key(self.elements,Ex,"I") 
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
        s = Group(subgroup,mult_table)
        return s          
class Representation(Group):
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
class Pion:
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
        new_p = Pion(momentum, sign)
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
class Two_Pion:
    def __init__(self,pi1,pi2,Momentum = [0,0,0]):                #attributes of TwoP system with capital letters
        assert isinstance(pi1,Pion)
        assert isinstance(pi2,Pion)
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
        new_twop = Two_Pion(p1,p2,M)
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