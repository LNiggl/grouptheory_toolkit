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

####################### GENERATION OF O_h #########################################

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
    actions = gen_actions.copy()
    known_a_o = {"I": object.copy()}       
    #plan: apply actions to all objects in known objects, if result is new, add to known_objects with name for the group element
    n = 0
    while n != len(known_a_o):                          
        n = len(known_a_o)
        for a in actions:
            for key, obj in known_a_o.items():                
                temp = obj.copy()
                temp.action(a)
                if not all([temp.is_equal(x) for x in known_a_o.values()]):
                # if new_o not in known_a_o.values():
                    new_name = remove_I(a + key)
                    known_a_o[new_name] = temp
                    actions.append(new_name)      
    return known_a_o


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


