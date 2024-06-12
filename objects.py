import numpy as np
import numbers
import groups as g
# Structure:
# -every class should have: copy(),update(),printname(),is_equal_to(x), is_negative_of(x),action(),inverse_action()
# -action and inverse_action should require minimal group theoretical/group specific information
# -these objects are used by the functions in representations.py to create representations as Representation objects 
class Scalar:
    def __init__(self,n):
        self.num = n
        self.update()
    def __len__(self):
        return 1
    def copy(self):
        new_scalar = Scalar(self.num)
        return new_scalar
    def update(self):
        self.name = "Scalar(" + str(self.num) + ")"
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,s):
        if isinstance(s,Scalar):
            if self.num == s.num:
                return True
            if isinstance(self.num,numbers.Number) and isinstance(s.num,numbers.Number):
                return abs(self.num - s.num) < 1e-6
            return False
        if self.num == s:
            return True
        if isinstance(self.num,numbers.Number) and isinstance(s,numbers.Number):
                return abs(self.num - s) < 1e-6
        return False
    def is_negative_of(self,s):        
        if isinstance(s,Scalar):
            ms = minus(s.num)
        else:
            ms = minus(s)
        return self.is_equal_to(ms)
    def action(self,A):
        pass
    def inverse_action(self,A):
        pass

class PseudoScalar:
    def __init__(self,n):
        self.num = n
        self.update()
    def __len__(self):
        return 1
    def copy(self):
        new_scalar = PseudoScalar(self.num)
        return new_scalar
    def update(self):
        self.name = "PScalar(" + str(self.num) + ")"
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,s):
        if isinstance(s,PseudoScalar):
            if self.num == s.num:
                return True
            if isinstance(self.num,numbers.Number) and isinstance(s.num,numbers.Number):
                return abs(self.num - s.num) < 1e-6
            return False
        if self.num == s:
            return True
        if isinstance(self.num,numbers.Number) and isinstance(s,numbers.Number):
                return abs(self.num - s) < 1e-6
        return False
    def is_negative_of(self,s):        
        if isinstance(s,PseudoScalar):
            ms = minus(s.num)
        else:
            ms = minus(s)
        return self.is_equal_to(ms)
    def gen_action(self,A):
        if A == "Inv":
            self.num = minus(self.num)
            self.update()
    def action(self,A):
        As = g.split(A)
        for i in range(len(As)):
            self.gen_action(As[i])
        self.update()
    def inverse_action(self,A):
        As = g.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.gen_action(As[i])
        self.update()

class Vector:
    def __init__(self,v):
        assert len(v) == 3
        self.x = v[0]
        self.y = v[1]
        self.z = v[2]
        self.update()
    def __len__(self):
        return 3
    def __lt__(self,obj):               # ONLY use for sorting, not actual comparisons
        if self.x < obj.x and not abs(self.x - obj.x) < 1e-5:
            return True
        elif abs(self.x - obj.x) < 1e-5:
            if self.y < obj.y:
                return True 
            elif abs(self.y - obj.y) < 1e-5 and self.z < obj.z:
                return True 
        return False

    def copy(self):
        x = self.x
        y = self.y
        z = self.z
        new_v = Vector([x,y,z])
        return new_v   
    def update(self):
        self.name = "Vec(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"     
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,v):                                    
        if all([self.x == v.x,self.y == v.y, self.z == v.z]):
            return True
        if all([isinstance(self.x,numbers.Number),isinstance(self.y,numbers.Number),isinstance(self.z,numbers.Number),\
                isinstance(v.x,numbers.Number),isinstance(v.y,numbers.Number),isinstance(v.z,numbers.Number)]):
            return np.allclose([self.x,self.y,self.z],[v.x,v.y,v.z], rtol = 1e-5, atol = 1e-10)
        return False
    def is_negative_of(self,nv):
        v = nv.copy()    
        v.action("Inv")
        return self.is_equal_to(v)
    def gen_action(self,A):
        assert A in {"I","Rot0","Rot1","Rot2","Inv"}        #extend if another formalism are added
        if A in {"I","Rot0","Rot1","Rot2","Inv"}:
            #if A == "I": change nothing
            if A == "Rot0":
                temp = self.copy()
                self.x = temp.x
                self.y = minus(temp.z)
                self.z = temp.y
            if A == "Rot1":
                temp = self.copy()
                self.x = temp.z
                self.y = temp.y
                self.z = minus(temp.x)
            if A == "Rot2":
                temp = self.copy()
                self.x = minus(temp.y)
                self.y = temp.x
                self.z = temp.z
            if A == "Inv":
                temp = self.copy()
                self.x = minus(temp.x)
                self.y = minus(temp.y)
                self.z = minus(temp.z)    
        self.update()
    def action(self,A):
        As = g.split(A)
        for i in range(len(As)):
            self.gen_action(As[i])            
    def inverse_gen_action(self,A):
        assert A in {"I","Rot0","Rot1","Rot2","Inv"}        #extend if another formalism are added
        if A in {"I","Rot0","Rot1","Rot2","Inv"}:
            #if A == "I": change nothing
            if A == "Rot0":
                temp = self.copy()
                self.x = temp.x
                self.y = temp.z
                self.z = minus(temp.y)
            if A == "Rot1":
                temp = self.copy()
                self.x = minus(temp.z)
                self.y = temp.y
                self.z = temp.x
            if A == "Rot2":
                temp = self.copy()
                self.x = temp.y
                self.y = minus(temp.x)
                self.z = temp.z
            if A == "Inv":
                temp = self.copy()
                self.x = minus(temp.x)
                self.y = minus(temp.y)
                self.z = minus(temp.z)
        self.update()
    def inverse_action(self,A):
        As = g.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.inverse_gen_action(As[i]) 
    @classmethod
    def sort(cls,list_v):
        list_v.sort(key = lambda v: v.z , reverse = True)
        list_v.sort(key = lambda v: v.y , reverse = True)
        list_v.sort(key = lambda v: v.x , reverse = True) 

def vector_sum(vecs):                       # input: list of Vector objects. Output: Vector object  
    x_tot = 0
    y_tot = 0
    z_tot = 0
    for v in vecs:
        x_tot += v.x
        y_tot += v.y
        z_tot += v.z
    return Vector([x_tot,y_tot,z_tot])

class Pion:
    def __init__(self,momentum,sign):        
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if isinstance(sign,PseudoScalar):
            if not isinstance(sign.num,numbers.Number):
                sign.num = int_sign(sign.num)
            assert sign.num == 1 or sign.num == -1
            self.sign = sign
        else:
            if not isinstance(sign,numbers.Number):
                sign = int_sign(sign)
            assert sign == 1 or sign == -1    
            self.sign = PseudoScalar(sign)
        self.update()
    def __lt__(self,p):                         # ONLY use for sorting, no comparison
        if self.sign.num < p.sign.num:
            return True
        elif abs(self.sign.num - p.sign.num) < 1e-5:            
            return self.momentum < p.momentum
        return False
    def update(self):
        self.name = str_sign(self.sign.num) + "Pi(" + str(self.momentum.x) + "," + str(self.momentum.y) + "," + str(self.momentum.z) + ")"
    def copy(self):
        m = self.momentum.copy()
        s = self.sign.copy()
        new_p = Pion(m,s)
        return new_p
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,p):
        return self.momentum.is_equal_to(p.momentum) and self.sign.is_equal_to(p.sign)
    def is_negative_of(self,p):
        return self.momentum.is_equal_to(p.momentum) and self.sign.is_negative_of(p.sign)
    def action(self,A):
        self.momentum.action(A)
        self.sign.action(A)
        self.update()
    def inverse_action(self,A):
        self.momentum.inverse_action(A)
        self.sign.inverse_action(A)
        self.update()

class Two_Pion:
    def __init__(self,p1,p2, distinguishable = True):                               # initialize from two Pion objects 
            assert isinstance(p1,Pion) and isinstance(p2,Pion)
            self.pi1 = p1
            self.pi2 = p2
            self.sign = overall_sign([p1,p2])
            self.distinguishable = distinguishable
            self.update()
    def __lt__(self,tp):
        if abs(self.sign.num - tp.sign.num) < 1e-5:
            return self.pi1.momentum < tp.pi1.momentum
        return self.sign.num < tp.sign.num 
    def copy(self):
        np1 = self.pi1.copy()
        np2 = self.pi2.copy()
        dist = self.distinguishable
        new_tp = Two_Pion(np1,np2,dist)
        return new_tp
    def update(self):
        self.sign = overall_sign([self.pi1,self.pi2])
        self.name = str_sign(self.sign.num) + "[Pi(" + str(self.pi1.momentum.x) + "," + str(self.pi1.momentum.y) + "," + str(self.pi1.momentum.z) + ")"\
            + "_x_Pi(" + str(self.pi2.momentum.x) + "," + str(self.pi2.momentum.y) + "," + str(self.pi2.momentum.z) + ")]"
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,tp):
        assert (self.distinguishable and tp.distinguishable) or (not self.distinguishable and not tp.distinguishable)        
        if self.distinguishable:
            return self.pi1.momentum.is_equal_to(tp.pi1.momentum) and self.pi2.momentum.is_equal_to(tp.pi2.momentum) and self.sign.is_equal_to(tp.sign)
        else:
            return (self.pi1.momentum.is_equal_to(tp.pi1.momentum) and self.pi2.momentum.is_equal_to(tp.pi2.momentum)\
                    or self.pi2.momentum.is_equal_to(tp.pi1.momentum) and self.pi1.momentum.is_equal_to(tp.pi2.momentum)) and self.sign.is_equal_to(tp.sign)
    def is_negative_of(self,tp):
        # return self.pi1.momentum.is_equal_to(tp.pi1.momentum) and self.pi2.momentum.is_equal_to(tp.pi2.momentum) and self.sign.is_negative_of(tp.sign)
        assert (self.distinguishable and tp.distinguishable) or (not self.distinguishable and not tp.distinguishable) 
        if self.distinguishable:
            return self.pi1.momentum.is_equal_to(tp.pi1.momentum) and self.pi2.momentum.is_equal_to(tp.pi2.momentum) and self.sign.is_negative_of(tp.sign)
        else:
            return (self.pi1.momentum.is_equal_to(tp.pi1.momentum) and self.pi2.momentum.is_equal_to(tp.pi2.momentum)\
                    or self.pi2.momentum.is_equal_to(tp.pi1.momentum) and self.pi1.momentum.is_equal_to(tp.pi2.momentum)) and self.sign.is_negative_of(tp.sign)
    def action(self,A):
        self.pi1.action(A)
        self.pi2.action(A)              # CAREFUL: no action on sign object, rather: update the sign
        self.update()
    def inverse_action(self,A):
        self.pi1.inverse_action(A)
        self.pi2.inverse_action(A)              # CAREFUL: no action on sign object, rather: update the sign
        self.update()

def overall_sign(obj_list):                 # compute overall sign from all occurring Scalar or Pseudoscalars (or other) signs and return as Scalar or PseudoScalar object
    S = 1
    num_PS = 0
    for obj in obj_list:
        if hasattr(obj,"sign"):
            if isinstance(obj.sign,Scalar):
                S = S*obj.sign.num
            if isinstance(obj.sign,PseudoScalar):
                S = S*obj.sign.num
                num_PS += 1
    if num_PS % 2:              # odd number of PseudoScalars
        return PseudoScalar(S)
    return Scalar(S)

def print_all(object_list):
    for o in object_list:
        o.printname()    
def minus(n):                       # multiplies by (-1) or modifies a String according to such a multiplication
    if isinstance(n,numbers.Number):
        return (-1)*n
    n = str(n)
    if n[0] == '-':
        return n[1:]
    return '-' + n
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
def orbit(object,group):                        #applies action of all Group.elements to object and returns list of resulting objects without duplicates
    orbit = [object.copy()]
    for g in group.elements:
        o = object.copy()
        o.action(g)
        if not any([q.is_equal_to(o) for q in orbit]):
            orbit.append(o)    
    return orbit
def remove_neg(orbit):                          # receives orbit, removes negatives if antiparallel pairs exist                 
    new_list = []
    for p in orbit:
        if hasattr(p,"sign"):
            if p.sign == "+":
                new_list.append(p)
    for p in orbit:                             # NEW loop: deal with non-signed objects AND catch outlier cases where only signed object with negative sign in orbit
        if all(not p.is_negative_of(o) for o in new_list):
            new_list.append(p)
    return new_list
def generate_basis(objects,group, drop_sign = True):          # receive list of objects, generates the combination of their orbits and outputs the underlying basis as list of objects     
    O = []
    for p in objects:
        for q in orbit(p,group):
            O.append(q)
    O.sort(reverse = True)
    # could add sorting for something like Pi(-1,0,0). functions <>.ispositive()?
    if drop_sign:
        O = remove_neg(O)                   
    return O  


def find_closure(gen_actions, object):       #takes list of Strings gen_actions, some object instance, and returns the closure regarding the objects resulting from combinations of actions  
    actions = {}
    actions["I"] = object.copy()
    for a in gen_actions:
        actions[a] = object.copy().action(a)    
    known_objects = {"I": object}       
    #plan: apply actions to all objects in known objects, if result is new, add to known_objects with name for the group element
    n = 0
    while n != len(known_objects):                          
        n = len(known_objects)
        for a in actions.copy().keys():
            for key, ob in known_objects.copy().items():  
                new_o = ob.copy()  
                new_o.action(a)            
                New = True
                for x in known_objects.values():
                    if new_o.is_equal_to(x):
                        New = False
                if New:
                    new_name = g.remove_I(a + key)
                    known_objects[new_name] = new_o
                    actions[new_name] = new_o    
    return actions

def mult_table_from_actions(actions):                             
    obj = actions["I"].copy()
    mt = {}
    for a in actions.keys():                            
        mt[a] = {}   
        for b in actions.keys():
            A = g.remove_I(a+b)
            w = obj.copy()
            w.action(A)                                 #result is w = Ae = abe for reference object e = obj   
            for c in actions.keys():
                x = obj.copy()
                x.action(c)
                if x.is_equal_to(w):
                    mt[a][b] = c                        # a.b = c because they have the same action            
    return mt
def generate_group(gen_actions,obj):
    actions = find_closure(gen_actions,obj)
    mt = mult_table_from_actions(actions)
    G = g.Group(list(actions.keys()),mt)
    return G















 




