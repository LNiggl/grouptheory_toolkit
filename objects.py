import numpy as np
import numbers
import groups as g
# Structure:
# -every class should have: copy(),update(),printname(),is_equal_to(x), is_negative_of(x),action(),inverse_action()
# -action and inverse_action should require minimal group theoretical/group specific information
# -these objects are used by the functions in representations.py to create representations as Representation objects
# classes that directly know about their actions: Scalar, PseudoScalar, Vector 
class Scalar:
    def __init__(self,n):
        self.direction_action = "left"
        self.num = n
        self.update()
    def __len__(self):
        return 1
    def __lt__(self,s):
        return self.num < s.num
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
        self.direction_action = "left"
        self.num = n
        self.update()
    def __len__(self):
        return 1
    def __lt__(self,s):
        return self.num < s.num
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
        # else: do nothing
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
        self.direction_action = "left"
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
 
class PseudoVector:
    def __init__(self,v):
        self.direction_action = "left"
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
        new_v = PseudoVector([x,y,z])
        return new_v   
    def update(self):
        self.name = "PVec(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"     
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
        v.x = minus(v.x)
        v.y = minus(v.y)
        v.z = minus(v.z)
        return self.is_equal_to(v)
    def gen_action(self,A):
        assert A in {"I","Rot0","Rot1","Rot2","Inv"}        #extend if another formalism are added
        if A in {"I","Rot0","Rot1","Rot2","Inv"}:
            #if A == "I": pass
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
            # if A == "Inv":
            #     pass    
        self.update()
    def action(self,A):
        As = g.split(A)
        for i in range(len(As)):
            self.gen_action(As[i])            
    def inverse_gen_action(self,A):
        assert A in {"I","Rot0","Rot1","Rot2","Inv"}        #extend if another formalism are added
        if A in {"I","Rot0","Rot1","Rot2","Inv"}:
            #if A == "I": pass
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
            # if A == "Inv":
            #     pass
        self.update()
    def inverse_action(self,A):
        As = g.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.inverse_gen_action(As[i]) 
class Tensor:                               # (spacial components of a) Lorentz-Tensor T: T^{mu,nu} -> Lambda^{rho}_mu Lambda^{sigma}_nu T^{rho,sigma}
    def __init__(self,struc = None,sign = None, symmetric = False, antisymmetric = False):          # attributes: tuple struc from {0,..2}_x_{0,..2}, int sign
                                                                                                    # ONLY assign symm/antisymm via constructor
        if struc == None:
            self.idx1 = 0
            self.idx2 = 0
        else:
            assert len(struc) == 2
            self.idx1 = struc[0]
            self.idx2 = struc[1]
        if sign == None:
            self.sign = 1
        else:
            assert isinstance(sign,numbers.Number)
            self.sign = sign
        assert not (symmetric and antisymmetric)
        self.symmetric = symmetric
        self.antisymmetric = antisymmetric
        self.update()
    def copy(self):
        mu = self.idx1 + 0
        nu = self.idx2 + 0
        struc = [mu,nu]
        sign = self.sign
        symm = self.symmetric
        antisymm = self.antisymmetric
        new_ten = Tensor(struc,sign,symmetric = symm,antisymmetric = antisymm)
        return new_ten   
    def update(self):
        self.name = str_sign(self.sign) + "Tensor{" + str(self.idx1) + "," + str(self.idx2) + "}"  
        if self.symmetric:
            self.name += "[symm.]"
        if self.antisymmetric:
            self.name += "[antisymm.]"
    def printname(self):
        self.update()
        print(self.name)
    def assign_idx_from_vecs(self,v1,v2):
        S = 1
        if not v1.x == 0:
            if v1.x == -1:
                S = S*(-1)
            mu = 0
        if not v1.y == 0:
            if v1.y == -1:
                S = S*(-1)
            mu = 1
        if not v1.z == 0:
            if v1.z == -1:
                S = S*(-1)
            mu = 2
        
        if not v2.x == 0:
            if v2.x == -1:
                S = S*(-1)
            nu = 0
        if not v2.y == 0:
            if v2.y == -1:
                S = S*(-1)
            nu = 1
        if not v2.z == 0:
            if v2.z == -1:
                S = S*(-1)
            nu = 2
        self.sign = self.sign*S
        self.idx1 = mu
        self.idx2 = nu
    def action(self,A):
        v1 = idx_to_vec(self.idx1)
        v2 = idx_to_vec(self.idx2)
        v1.action(A)
        v2.action(A)
        self.assign_idx_from_vecs(v1,v2)
        self.update()
    def inverse_action(self,A):
        v1 = idx_to_vec(self.idx1)
        v2 = idx_to_vec(self.idx2)
        v1.inverse_action(A)
        v2.inverse_action(A)
        self.assign_idx_from_vecs(v1,v2)
        self.update()
    def is_equal_to(self,t):     
        assert isinstance(t,Tensor)     
        if self.symmetric:    
            assert t.symmetric          
            return ((self.idx1 == t.idx1 and self.idx2 == t.idx2) or (self.idx1 == t.idx2 and self.idx2 == t.idx1)) and self.sign == t.sign
        if self.antisymmetric:
            assert t.antisymmetric
            return (self.idx1 == t.idx1 and self.idx2 == t.idx2 and self.sign == t.sign) or (self.idx1 == t.idx2 and self.idx2 == t.idx1 and self.sign == -1*t.sign)
        return self.idx1 == t.idx1 and self.idx2 == t.idx2 and self.sign == t.sign 
    def is_negative_of(self,t):     
        assert isinstance(t,Tensor)   
        if self.symmetric:    
            assert t.symmetric          
            return ((self.idx1 == t.idx1 and self.idx2 == t.idx2) or (self.idx1 == t.idx2 and self.idx2 == t.idx1)) and self.sign == -1*t.sign
        if self.antisymmetric:
            assert t.antisymmetric
            return (self.idx1 == t.idx1 and self.idx2 == t.idx2 and self.sign == -1*t.sign) or (self.idx1 == t.idx2 and self.idx2 == t.idx1 and self.sign == t.sign)                            
        return self.idx1 == t.idx1 and self.idx2 == t.idx2 and self.sign == -1*t.sign   
def idx_to_vec(n):                          # takes int from 0,1,2 and creates Vector object as basis vector in that direction
    a1 = []
    i = 0
    assert n in {0,1,2}
    while i < 3:
        if i == n:
            a1.append(1)                
        else:
            a1.append(0)
        i += 1
    v = Vector(a1)
    return v

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
        self.momentum.inverse_action(A)
        self.sign.action(A)
        self.update()
    def inverse_action(self,A):
        self.momentum.action(A)
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
        if isinstance(obj,numbers.Number):
            assert obj in {-1,1}
            S = obj*S
            num = True
    if num:
        return S
    if num_PS % 2:              # odd number of PseudoScalars
        return PseudoScalar(S)
    return Scalar(S)
class Rho:
    def __init__(self,momentum,struc = None):
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.vector_structure = Vector(["rho_x","rho_y","rho_z"])
        elif isinstance(struc,Vector):
            self.vector_structure = struc
        else:
            assert len(struc) == 3
            self.vector_structure = Vector(struc)
        self.update()
    def __lt__(self,obj):               # ONLY use for sorting, not actual comparison
        if not abs(coord_to_num(self.vector_structure.x) - coord_to_num(obj.vector_structure.x)) < 1e-5: 
            return coord_to_num(self.vector_structure.x) < coord_to_num(obj.vector_structure.x)
        elif not abs(coord_to_num(self.vector_structure.y) - coord_to_num(obj.vector_structure.y)) < 1e-5: 
            return coord_to_num(self.vector_structure.y) < coord_to_num(obj.vector_structure.y)
        elif not abs(coord_to_num(self.vector_structure.z) - coord_to_num(obj.vector_structure.z)) < 1e-5: 
            return coord_to_num(self.vector_structure.z) < coord_to_num(obj.vector_structure.z)
        return self.momentum < obj.momentum
    def copy(self):
        mom = self.momentum.copy()
        struc = self.vector_structure.copy()
        new_rho = Rho(mom,struc)
        return new_rho
    def update(self):
        self.name = self.vector_structure.name + self.momentum.name[3:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.vector_structure.is_equal_to(R.vector_structure) and self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return self.vector_structure.is_negative_of(R.vector_structure) and self.momentum.is_equal_to(R.momentum)
    def action(self,A):
        self.vector_structure.action(A)
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.vector_structure.inverse_action(A)
        self.momentum.action(A)
        self.update()
class ScalarField:
    def __init__(self,momentum,struc = None):
        self.direction_action = "right"
        self.lorentz_structure = "scalar"
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.structure = Scalar(1)
        elif isinstance(struc,Scalar):
            self.structure = struc
        else:
            assert isinstance(struc,numbers.Number)
            self.structure = Scalar(struc)
        self.update()
    def __lt__(self,obj):
        return self.momentum < obj.momentum
    def copy(self):
        mom = self.momentum.copy()        
        new = ScalarField(mom)
        return new
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "ScalarField{" + self.structure.name[8:-1] + "}(p=" + self.momentum.name[4:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return False                    # impossible by construction 
    def action(self,A):
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.momentum.action(A)
        self.update()
class PseudoScalarField:                               # particle with: Lorentz pseudoscalar structure (under rotations), momentum
    def __init__(self,momentum,struc = None):
        self.direction_action = "right"
        self.lorentz_structure = "pseudoscalar"
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.structure = PseudoScalar(1)
        elif isinstance(struc,PseudoScalar):
            self.structure = struc
        else:
            assert isinstance(struc,numbers.Number)
            self.structure = PseudoScalar(struc)
        self.update()
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        new = PseudoScalarField(mom,struc)
        return new
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "PseudoScalarField{" + self.structure.name[8:-1] + "}(p=" + self.momentum.name[4:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    def action(self,A):
        self.structure.action(A)
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.structure.inverse_action(A)
        self.momentum.action(A)
        self.update()

class VectorField:                               # particle with: inner vector structure (under rotations), momentum
    def __init__(self,momentum,struc = None):
        self.direction_action = "mixed"
        self.lorentz_structure = "vector"
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.structure = Vector([1,0,0])
        elif isinstance(struc,Vector):
            self.structure = struc
        else:
            assert len(struc) == 3
            self.structure = Vector(struc)
        self.update()
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        new = VectorField(mom,struc)
        return new
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "VectorField{" + self.structure.name[3:] + "}(p=" + self.momentum.name[4:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    def action(self,A):
        self.structure.action(A)
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.structure.inverse_action(A)
        self.momentum.action(A)
        self.update()
    def set_vals(self,structure = None,momentum = None):
        if not structure == None:
            assert isinstance(structure,Vector)
            self.structure = structure
        if not momentum == None:
            assert isinstance(momentum,Vector)
            self.momentum = momentum

class PseudoVectorField:                               # particle with: Lorentz pseudovector structure (under rotations), momentum
    def __init__(self,momentum,struc = None):
        self.direction_action = "mixed"
        # for classes with "mixed": add function left_action, check for direction_action and use left_action for rep_functions as needed
        self.lorentz_structure = "pseudovector"
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.structure = PseudoVector([1,0,0])
        elif isinstance(struc,PseudoVector):
            self.structure = struc
        else:
            assert len(struc) == 3
            self.structure = PseudoVector(struc)
        self.update()
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        new = PseudoVectorField(mom,struc)
        return new
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "PseudoVectorField{" + self.structure.name[4:] + "}(p=" + self.momentum.name[4:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    def action(self,A):
        self.structure.action(A)
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.structure.inverse_action(A)
        self.momentum.action(A)
        self.update()
    def set_vals(self,structure = None,momentum = None):
        if not structure == None:
            assert isinstance(structure,PseudoVector)
            self.structure = structure
        if not momentum == None:
            assert isinstance(momentum,Vector)
            self.momentum = momentum
class TensorField:                      # Field operator(with momentum) that transforms like the spacial components of a rank-2 Lorentz-tensor T: T^{mu,nu} -> Lambda^{rho}_mu Lambda^{sigma}_nu T^{rho,sigma} 
    def __init__(self,momentum,structure = None,sign = 1,symmetric = False,antisymmetric = False):
            self.lorentz_structure = "tensor"
            if isinstance(momentum,Vector):
                self.momentum = momentum
            else:
                assert len(momentum) == 3
                self.momentum = Vector(momentum)
            if isinstance(structure,Tensor):
                self.structure = structure
            else:
                self.structure = Tensor(struc = structure,sign = sign, symmetric = symmetric,antisymmetric = antisymmetric)
            self.update()
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        new_t = TensorField(mom,struc)
        return new_t
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "TensorField" + self.structure.name[7:] + "(p=" + self.momentum.name[4:]
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    def action(self,A):
        self.structure.action(A)
        self.momentum.inverse_action(A)
        self.update()
    def inverse_action(self,A):
        self.structure.inverse_action(A)
        self.momentum.action(A)
        self.update()
class LinearCombination:    # UNFINISHED
    def __init__(self,basis,weights,linear_dep = True):                       # list of weights of each element in basis; basis: list of objects, e.g. of class <Vector>; 
                                                                            # ! all objects in basis are assumed to be linearly independent
        self.basis = basis
        self.linear_dep = linear_dep                                            # True: treat -obj and obj as linearly dependent; False: treat as linearly independent
        self.lin_comb = []
        assert len(basis) == len(weights)
        for i in range(len(basis)):
            self.lin_comb.append(TensorProduct(Scalar(weights[i]),basis[i]))
        self.update()

    def action(self,A):
        for l in self.lin_comb:
            l.action(A)
            if match_in_list(l[1],self.basis) == None:                         # means result is not in basis -> -result should be in basis 
                l.action("Inv")
                
    
class TensorProduct:       #UNFINISHED         # tensor product of arbitrary number and types of classes; each attribute of a TensorProduct instance is one object, e.g. VectorField
                                                # extend to deal nicely with Scalars and their accumulative signs
    def __init__(self,*v,distinguishable = True):  
        self.obj = {}                   
        for i in range(len(v)):
            self.obj[i] = v[i]
        self.distinguishable = distinguishable
        if not all(self.obj[i].direction_action == self.obj[0].direction_action for i in range(len(v))):
            print("action direction problem")
        self.direction_action = self.obj[0].direction_action
        self.update()
    def __lt__(self,t):                         # only for sorting purposes
        if isinstance(self.obj[0],PseudoScalarField):
            return self.obj[0].momentum < t.obj[0].momentum
        return self.obj[0] < t.obj[0]
        # i = 0
        # while i in range(len(list(self.obj.values()))):
        #     if self.obj[i] < t.obj[i]:
        #         return True
        #     if self.obj[i] > t.obj[i]:
        #         return False
        #     i += 1
        # print("TensorProduct.__lt__: Comparison failed")
        # return False

    def copy(self):
        temp = list(self.obj.values())
        new_obj = []
        for i in range(len(temp)):
            new_obj.append(temp[i].copy())
        dist = (self.distinguishable == True)
        new_t = TensorProduct(*new_obj,distinguishable = dist)
        return new_t
    def update(self):
        self.name = self.obj[0].name
        for i in range(len(self.obj)-1):
            self.name += "_x_" + self.obj[i+1].name
        if self.distinguishable:
            self.name+= "[dist.]"
        else:
            self.name += "[indist.]"
    def printname(self):
        self.update()
        print(self.name)
    def is_equal_to(self,t):                            # for DISTINGUISHABLE: True if each object (in order) is exactly the same or opposite, and n_antiparallel is even
                                                        # for INDISTINGUISHABLE: True if objects can be paired to equal or antiparallel pairs, and n_antiparallel is even 
        if len(self.obj) != len(t.obj):
            return False
        if self.distinguishable:
            assert t.distinguishable
            n_diff = 0
            for i in range(len(self.obj)):
                if not self.obj[i].is_equal_to(t.obj[i]):
                    if not self.obj[i].is_negative_of(t.obj[i]):
                        return False
                    n_diff += 1
            return n_diff % 2 == 0
        else:                       # self is indistinguishable product space
            pairs = self.pair_up(t)
            if pairs == None: 
                return False
            return pairs[1] % 2 == 0        # pairs[1]: number of pairs with differing signs. (-1)^n cancels out for n even            
    def is_negative_of(self,t):                 # for DISTINGUISHABLE: True if each object (in order) is exactly the same or opposite, and n_antiparallel is odd
                                                # for INDISTINGUISHABLE: True if objects can be paired to equal or antiparallel pairs, and n_antiparallel is odd 
        if len(self.obj) != len(t.obj):
            return False
        if self.distinguishable:
            n_diff = 0
            assert t.distinguishable
            for i in range(len(self.obj)):
                if not self.obj[i].is_equal_to(t.obj[i]):
                    if not self.obj[i].is_negative_of(t.obj[i]):
                        return False
                    n_diff += 1
            return n_diff % 2 == 1
        else:                       # self is indistinguishable product space
            pairs = self.pair_up(t)
            if pairs == None: 
                return False
            return pairs[1] % 2 == 1        # pairs[1]: number of pairs with differing signs. (-1)^n = -1 for n odd
    def action(self,A):
        for o in self.obj.values():
            o.action(A)
        self.update()
    def inverse_action(self,A):
        for o in self.obj.values():
            o.inverse_action(A)
        self.update()
    def pair_up(self,prod):                     # find pairs of objects; return [[[pair1],[pair2],..],n_different_signs]
                                                # meant for two systems of indistinguishable constituents
        assert (not self.distinguishable) and (not prod.distinguishable)
        objects1 = list(self.obj.values())
        objects2 = list(prod.obj.values())
        assert len(objects1) == len(objects2)
        pairs = []
        n_sign_diff = 0
        for o1 in objects1.copy():
            o2 = match_in_list(o1,objects2)
            if o2 == None:                      # failure to find exact match: Look for "antiparallel" pair
                neg_o1 = invert(o1.copy())
                # if hasattr(neg_o1,"lorentz_structure"):
                #     neg_o1.structure.action("Inv")
                # else:
                #     neg_o1.action("Inv")

                o2 = match_in_list(neg_o1,objects2)
                if o2 == None:
                    return None                 # failure to pair up at all: Return None
                else: 
                    n_sign_diff += 1
            pairs.append([o1,o2])
            objects2.remove(o2)
        return [pairs,n_sign_diff]
        


def print_all(object_list):
    for o in object_list:
        if isinstance(o,list):
            print("[")
            print_all(o)
            print("]")
        else:
            o.printname()    
def minus(n):                       # multiplies by (-1) or modifies a String according to such a multiplication
    if isinstance(n,numbers.Number):
        return (-1)*n
    n = str(n)
    if n[0] == '-':
        return n[1:]
    return '-' + n
def invert(object):
    obj = object.copy()                    
    if hasattr(obj,"lorentz_structure"):                # <..>Field classes
        if "vector" in obj.lorentz_structure:           # Vector and PseudoVector                
            obj.structure = invert(obj.structure.copy)
            # obj.set_vals(structure = obj.structure)               
        elif "scalar" in obj.lorentz_structure:         # Scalar and PseudoScalar
            print("lorentz (p)scalar")
            obj.structure = invert(obj.structure)
        # extend for different structures
    else:
        if hasattr(obj,"num"):              # Scalar and PseudoScalar   
            obj.num = minus(obj.num)
            return obj
        elif hasattr(obj,"x"):              # Vector and PseudoVector
            obj.x = minus(obj.x)
            obj.y = minus(obj.y)
            obj.z = minus(obj.z)
            return obj
        # extend for different structures

def negative_pair(m,n):                 #returns True if m = -n up to some precision
    return abs(m+n) < 1e-8
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
def coord_to_num(s):                # used for sorting in Rho __lt__: String s containing x,y, or z and possibly a "-" ->  corresponding int: x -> 3, y -> 2..
    s = str(s)
    n = 0
    if "x" in s:
        n = 3
    if "y" in s:
        n = 2
    if "z" in s:
        n = 1
    if s[0] == "-":
        n = n*(-1)
    if n == 0:
        return None
    return n
def match_in_list(x,l_y):                    # returns None if x is not equal up to numerical error to any value in l_y, otherwhise returns first match
    if hasattr(x, "is_equal_to"):
        for y in l_y:
            if x.is_equal_to(y):
                return y
    else:
        for y in l_y:
            if abs(x-y) < 1e-8:
                return y
    return None
def orbit(object,group):                        #applies action of all Group.elements to object and returns list of resulting objects without duplicates
    if object == None:
        return None
    if not hasattr(object,"lorentz_structure"):
        orbit = [object.copy()]
        for g in group.elements:
            o = object.copy()
            o.action(g)
            if not any([q.is_equal_to(o) for q in orbit]):
                orbit.append(o)  
    else:
        orbit = exhaustive_orbit(object,group)
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
    if not hasattr(objects,"__iter__"):
        objects = [objects]
    for p in objects:
        for q in orbit(p,group):
            if not any([q.is_equal_to(r) for r in O]):
                O.append(q)
    O.sort(reverse = True)
    # could add sorting for something like Pi(-1,0,0). functions <>.ispositive()?
    if drop_sign:
        O = remove_neg(O)                   
    return O  
def exhaustive_orbit(object,group):                   # find all possible combinations of orbits of each of the attributes; attributes currently: structure and momentum
    struc_orbit = None
    momentum_orbit = None
    if hasattr(object,"structure"):
        struc_orbit = orbit(object.structure,group)
    if hasattr(object,"momentum"):
        momentum_orbit = orbit(object.momentum,group)
    if struc_orbit == None:
        if momentum_orbit == None:
            return None
        else:
            orb = []
            for m in momentum_orbit:
                o = object.copy()
                o.momentum = m
                orb.append(o)
            return orb
    else:
        orb = []
        if momentum_orbit == None:
            for s in struc_orbit:
                o = object.copy()
                o.structure = s
                orb.append(o)
            return orb
        else:
            for s in struc_orbit:
                for m in momentum_orbit:
                    o = object.copy()
                    o.structure = s
                    o.momentum = m
                    orb.append(o)
            return orb

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















 




