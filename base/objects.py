import numbers

import numpy as np

from . import groups as g
from . import tools as t
from .definitions import num_tol, num_rtol


# Structure:
# -every class has: copy(),update(),printname(),is_equal_to(x), is_negative_of(x),action(),inverse_action()
# -attribute: direction_action: wether action is a left or right action
# -these objects are used by the functions in representations.py to create representations as Representation objects
# - some classes have: name_gpt attribute for export of data in gpt-conform format; projection(x) if comparisons is_equal_to(x) and is_negative_of(x) are insufficient for purpose 
 


class Scalar:
    def __init__(self,n):
        self.direction_action = "left"
        self.num = n
        self.update()

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
                return abs(self.num - s.num) < num_tol
            return False
        if self.num == s:
            return True
        if isinstance(self.num,numbers.Number) and isinstance(s,numbers.Number):
                return abs(self.num - s) < num_tol
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

    def __len__(self):
        return 1
    
    def __lt__(self,s):
        return self.num < s.num
    
    
class PseudoScalar:
    def __init__(self,n):
        self.direction_action = "left"
        self.num = n
        self.update()

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
                return abs(self.num - s.num) < num_tol
            return False
        if self.num == s:
            return True
        if isinstance(self.num,numbers.Number) and isinstance(s,numbers.Number):
                return abs(self.num - s) < num_tol
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
        As = t.split(A)
        for i in range(len(As)):
            self.gen_action(As[i])
        self.update()

    def inverse_action(self,A):
        As = t.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.gen_action(As[i])
        self.update()

    def __len__(self):
        return 1
    
    def __lt__(self,s):
        return self.num < s.num
        

class Vector:
    def __init__(self,v):
        self.direction_action = "left"
        assert len(v) == 3
        self.x = v[0]
        self.y = v[1]
        self.z = v[2]
        self.update()

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
        if all([self.x == v.x,self.y == v.y, self.z == v.z]):                   # for symbolic maths
            return True
        if all([isinstance(self.x,numbers.Number),isinstance(self.y,numbers.Number),isinstance(self.z,numbers.Number),\
                isinstance(v.x,numbers.Number),isinstance(v.y,numbers.Number),isinstance(v.z,numbers.Number)]):
            return np.allclose([self.x,self.y,self.z],[v.x,v.y,v.z], rtol = num_rtol, atol = num_tol)
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
        As = t.split(A)
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
        As = t.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.inverse_gen_action(As[i])

    def __len__(self):
        return 3
    
    def __lt__(self,obj):                       # only use for sorting purposes
        if all([isinstance(self.x,numbers.Number),isinstance(self.y,numbers.Number),isinstance(self.z,numbers.Number),\
                isinstance(obj.x,numbers.Number),isinstance(obj.y,numbers.Number),isinstance(obj.z,numbers.Number)]):
            if self.x < obj.x and not abs(self.x - obj.x) < num_tol:
                return True
            elif abs(self.x - obj.x) < num_tol:
                if self.y < obj.y:
                    return True 
                elif abs(self.y - obj.y) < num_tol and self.z < obj.z:
                    return True 
            return False
        else:                       #  for symbolic maths
            if str(self.x)[0] < str(obj.x)[0]:
                return True
            elif str(self.x)[0] == str(obj.x)[0]:
                if str(self.x)[0]== '-':
                    return str(self.x)[1]> str(obj.x)[1]
                else:
                    if str(self.y)[0] < str(obj.y)[0]:
                        return True
                    elif str(self.y)[0] == str(obj.y)[0]:
                        if str(self.y)[0] == '-':
                            return str(self.y)[1] > str(obj.y)[1] 
                        else: 
                            if str(self.z)[0] < str(obj.z)[0]:
                                return True
                            elif str(self.z)[0] == str(obj.z)[0]:
                                if str(self.z)[0] == '-':
                                    return str(self.z)[1] > str(obj.z)[1]           
            return False
        
        
        
 
class PseudoVector:
    def __init__(self,v):
        self.direction_action = "left"
        assert len(v) == 3
        self.x = v[0]
        self.y = v[1]
        self.z = v[2]
        self.update()

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
            return np.allclose([self.x,self.y,self.z],[v.x,v.y,v.z], rtol = num_rtol, atol = num_tol)
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
        As = t.split(A)
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
        As = t.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.inverse_gen_action(As[i]) 
    
    def __len__(self):
        return 3
    
    def __lt__(self,obj):               # only use for sorting, not actual comparisons
        if self.x < obj.x and not abs(self.x - obj.x) < num_tol:
            return True
        elif abs(self.x - obj.x) < num_tol:
            if self.y < obj.y:
                return True 
            elif abs(self.y - obj.y) < num_tol and self.z < obj.z:
                return True 
        return False
    

class Tensor:                               # (spatial components of a) Lorentz-Tensor T: T^{mu,nu} -> Lambda^{rho}_mu Lambda^{sigma}_nu T^{rho,sigma}
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


class WeylSpinor:                                     # \psi_L or \psi_R Weyl spinor(same trafo under spatial rotations(see P/S section 3.2)) as two-component, C-valued object
                                                      # note (from Wettig GT, ch. 7.3): Dirac spinor \psi = (\psi_L,\psi_R) transforms in (0,1/2)_x_(1/2,0): this is irrep iff parity is included
    def __init__(self,a,b,idx = "L",distinguish_by_type = False):
        self.direction_action = "left"
        self.idx = idx                                                     # to in principle be able to discern left from right Weyl spinors
        self.distinguish_by_type = distinguish_by_type                     # wether or not to treat R and L differently
        self.a = complex(a)
        self.b = complex(b)
        self.update()

    def copy(self):
        c = self.a
        d = self.b
        new_s = WeylSpinor(c,d)
        # print("copy result: ",type(new_s))
        return new_s
    
    def update(self):
        self.name = "WeylSpinor_" + self.idx + "{" + str(self.a) + "," + str(self.b) + "}" 

    def printname(self):
        self.update()
        print(self.name)  
    
    def is_equal_to(self,s):
        if self.distinguish_by_type:
            return self.idx == s.idx and (abs(self.a-s.a)< num_tol and abs(self.b-s.b) < num_tol)
        return (abs(self.a-s.a)< num_tol and abs(self.b-s.b) < num_tol)
    
    def is_negative_of(self,s):
        if self.distinguish_by_type:
            return self.idx == s.idx and (abs(self.a+s.a)< num_tol and abs(self.b+s.b) < num_tol)
        return (abs(self.a+s.a)< num_tol and abs(self.b+s.b) < num_tol)

    def gen_action(self,A):             
        if A == "A0":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(temp.a - 1j*temp.b)
            self.b = (1/np.sqrt(2))*(-1j*temp.a + temp.b)
        if A == "A1":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(temp.a - temp.b)
            self.b = (1/np.sqrt(2))*(temp.a + temp.b)
        if A == "A2":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(1-1j)*temp.a
            self.b = (1/np.sqrt(2))*(1+1j)*temp.b
        if A == "Rot0":
            temp = self.copy()
            self.a = np.sqrt(2)/2*(temp.a-temp.b*1j)
            self.b = np.sqrt(2)/2*(-temp.a*1j+temp.b)
        if A == "Rot1":
            temp = self.copy()
            self.a = np.sqrt(2)/2*(temp.a-temp.b)
            self.b = np.sqrt(2)/2*(temp.a+temp.b)
        if A == "Rot2":
            temp = self.copy()
            self.a = np.sqrt(2)/2*(temp.a*(1-1j))
            self.b = np.sqrt(2)/2*(temp.b*(1+1j))

    def action(self,A):
        As = t.split(A)
        for i in range(len(As)):
            AA = t.split(As[i])                             # for proper splitting of names from generating set [Inv,A0,A1,A2] 
            for j in range(len(AA)):
                self.gen_action(AA[j]) 
        self.update()

    def inverse_gen_action(self,A):
        if A == "A0":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(temp.a + 1j*temp.b)
            self.b = (1/np.sqrt(2))*(1j*temp.a + temp.b)
        if A == "A1":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(temp.a + temp.b)
            self.b = (1/np.sqrt(2))*(-temp.a + temp.b)
        if A == "A2":
            temp = self.copy()
            self.a = (1/np.sqrt(2))*(1+1j)*temp.a
            self.b = (1/np.sqrt(2))*(1-1j)*temp.b
        if A == "Rot0":
            temp = self.copy()
            self.a = np.sqrt(2)/4*(temp.a+temp.b*1j)
            self.b = np.sqrt(2)/4*(temp.a*1j+temp.b)
        if A == "Rot1":
            temp = self.copy()
            self.a = np.sqrt(2)/4*(temp.a+temp.b)
            self.b = np.sqrt(2)/4*(-temp.a+temp.b)
        if A == "Rot2":
            temp = self.copy()
            self.a = np.sqrt(2)/4*(temp.a*(1+1j))
            self.b = np.sqrt(2)/4*(temp.b*(1-1j))

    def inverse_action(self,A):
        As = t.split(A)[::-1]                                       # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            AA = t.split(As[i])[::-1]                               # for proper splitting of names from generating set [Inv,A0,A1,A2] 
            for j in range(len(AA)):
                self.inverse_gen_action(AA[j]) 
        self.update()
    
    def projection(self,s2):                # returns complex scalar product <self,s2> ; use for projection to orthogonal basis vectors
        if self.distinguish_by_type and self.idx != s2.idx:
            return 0
        return np.conj(self.a)*s2.a+np.conj(self.b)*s2.b     


class DiracSpinor:                                      # four-component, C-valued object. 
                                                        # Transformation action in chiral representation: block-diagonal for rotations, anti-blockdiagonal for parity
    def __init__(self,a,b,c,d, normalized = True, Weyl_indices = ["L","R"],LR_different = False):
        self.direction_action = "left"
        self.LR_different = LR_different
        a = complex(a)
        b = complex(b)  
        c = complex(c)
        d = complex(d)
        if normalized:
            n = np.sqrt(a*np.conj(a)+b*np.conj(b)+c*np.conj(c)+d*np.conj(d))
            a = a/n
            b = b/n
            c = c/n
            d = d/n
        self.s1 = WeylSpinor(a,b,idx = Weyl_indices[0],distinguish_by_type=self.LR_different)
        self.s2 = WeylSpinor(c,d,idx = Weyl_indices[1],distinguish_by_type=self.LR_different)
        self.update()

    def copy(self):
        new_s = DiracSpinor(self.s1.a,self.s1.b,self.s2.a,self.s2.b,Weyl_indices = [self.s1.idx,self.s2.idx])
        return new_s
    
    def update(self): 
        self.name = "DiracSpinor{" + self.s1.name + "," + self.s2.name + "}"  

    def printname(self):
        self.update()
        print(self.name)   

    def gen_action(self,A):             
        if A == "Inv":
            temp = self.copy()
            self.s1 = temp.s2
            self.s2 = temp.s1
        else:
            self.s1.action(A)
            self.s2.action(A)

    def action(self,A):
        As = t.split(A)
        for i in range(len(As)):
            self.gen_action(As[i])     

    def inverse_gen_action(self,A):
        if A == "Inv":
            temp = self.copy()
            self.s1 = temp.s2
            self.s2 = temp.s1
        else:
            self.s1.action(A)
            self.s2.action(A)

    def inverse_action(self,A):
        As = t.split(A)[::-1]               # the slicing reverses list as (AB)^-1 = B^-1 A^-1 
        for i in range(len(As)):
            self.inverse_gen_action(As[i])

    def is_equal_to(self,s):
        return self.s1.is_equal_to(s.s1) and self.s2.is_equal_to(s.s2)
    
    def is_negative_of(self,s):
        return self.s1.is_negative_of(s.s1) and self.s2.is_negative_of(s.s2)
    
    def projection(self,s2):                # projection of self-contribution in s2 via scalar product
        return self.s1.projection(s2.s1)+self.s2.projection(s2.s2)

##### classes representing fields   #####

class ScalarField:
    def __init__(self,momentum,struc = None,modified_momentum_trafo = False):
        self.direction_action = "right"
        self.modified_momentum_trafo = modified_momentum_trafo      # parameter. If False: \psi(p)->\psi(R^-1 p), if True: \psi(p)->\psi(Rp) under action R
        if self.modified_momentum_trafo:
            self.direction_action = "left"
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
        modified_mom_trafo = (self.modified_momentum_trafo == True)      
        new = ScalarField(mom,modified_momentum_trafo = modified_mom_trafo)
        return new
    
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "ScalarField{" + self.structure.name[8:-1] + "}(p=" + self.momentum.name[3:] + ")"

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,R):
        return self.momentum.is_equal_to(R.momentum)
    
    def is_negative_of(self,R):
        return False                    # impossible by construction 
    
    def action(self,A):
        if not self.modified_momentum_trafo:
            self.momentum.inverse_action(A)
        else:
            self.momentum.action(A)
        self.update()

    def inverse_action(self,A):
        if not self.modified_momentum_trafo:
            self.momentum.action(A)
        else:
            self.momentum.inverse_action(A)
        self.update()

class PseudoScalarField:                               # particle with: Lorentz pseudoscalar structure (under rotations), momentum
    def __init__(self,momentum,struc = None, modified_momentum_trafo = False):
        self.modified_momentum_trafo = modified_momentum_trafo                  # parameter. If False: \psi(p)->det(R)\psi(R^-1 p), if True: \psi(p)->det(R)\psi(Rp) under action R
        if self.modified_momentum_trafo: 
            self.direction_action = "left"
        else:
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
    
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        modified_mom_trafo = (self.modified_momentum_trafo == True)
        new = PseudoScalarField(mom,struc,modified_momentum_trafo = modified_mom_trafo)
        if hasattr(self, "name_gpt"):
            new.set_name_gpt("" + self.name_gpt)
        return new
    
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "PseudoScalarField{" + self.structure.name[8:-1] + "}(p=" + self.momentum.name[3:] + ")"
        if hasattr(self,"name_gpt"):
            self.update_name_gpt()

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)

    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def action(self,A):
        self.structure.action(A)
        if not self.modified_momentum_trafo:
            self.momentum.inverse_action(A)
        else:
            self.momentum.action(A)
        self.update()

    def inverse_action(self,A):
        self.structure.inverse_action(A)
        if not self.modified_momentum_trafo:
            self.momentum.action(A)
        else:
            self.momentum.inverse_action(A)
        self.update()
    
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
        
    def set_name_gpt(self,name):                                    # Name in accordance with gpt convention. To use e.g. in a basis of Representation 
        self.name_gpt = name

    def namestring_gpt_pion(self,type = "minus"):                   # treat this pseudoscalar object as \pi^(+/-) meson. Create according namestring for set_name_gpt
        name_gpt = "FACTOR " + str(np.real(self.structure.num)) + " " + str(np.imag(self.structure.num)) + "\n"
        if type == "minus":                                         # pi^- 
           name_gpt += "UBAR t\nGAMMA 5\nMOM ["
           name_gpt += str(self.momentum.x) + "," +  str(self.momentum.y) + "," +  str(self.momentum.z) + "] t\n" 
           name_gpt += "D t\n"
        if type == "plus":                                         # pi^- 
           name_gpt += "DBAR t\nGAMMA 5\nMOM ["
           name_gpt += str(self.momentum.x) + "," +  str(self.momentum.y) + "," +  str(self.momentum.z) + "] t\n" 
           name_gpt += "U t\n"
        #if zero pion: ...
        return name_gpt
    
    def update_name_gpt(self):                              # updates name_gpt; use in update() only if name_gpt exists
        # gpt naming in case of pions
        if "UBAR " in self.name_gpt and "D " in self.name_gpt:
            self.set_name_gpt(self.namestring_gpt_pion(type = "minus"))
        elif "DBAR " in self.name_gpt and "U " in self.name_gpt:
            self.set_name_gpt(self.namestring_gpt_pion(type = "plus"))
        # elif zero pion: ...


class VectorField:                               # particle with: inner vector structure (under rotations) and momentum as argument
    def __init__(self,momentum,struc = None, modified_momentum_trafo = False):        
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
        zero = Vector([0,0,0])
        if self.momentum.is_equal_to(zero):
            self.direction_action = "left"
        else: 
            self.direction_action = "mixed"
        self.modified_momentum_trafo = modified_momentum_trafo
        if self.modified_momentum_trafo:             # parameter. If False: V(p)->RV(R^-1 p), if True: \psi(p)->RV(Rp) under action R
            self.direction_action = "left"
        self.lorentz_structure = "vector"
        self.update()

    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        modified_mom_trafo = (self.modified_momentum_trafo == True)
        new = VectorField(mom,struc,modified_momentum_trafo=modified_mom_trafo)
        if hasattr(self,"name_gpt"):
            new.set_name_gpt(""+self.name_gpt)
        return new
    
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "VectorField{" + self.structure.name[3:] + "}(p=" + self.momentum.name[3:] + ")"
        if hasattr(self,"name_gpt"):
            self.update_name_gpt()

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def action(self,A):
        self.structure.action(A)
        if not self.modified_momentum_trafo:
            self.momentum.inverse_action(A)
        else:
            self.momentum.action(A)
        self.update()

    def inverse_action(self,A):
        self.structure.inverse_action(A)
        if not self.modified_momentum_trafo:
            self.momentum.action(A)
        else:
            self.momentum.inverse_action(A)
        self.update()

    def set_vals(self,structure = None,momentum = None):
        if not structure == None:
            assert isinstance(structure,Vector)
            self.structure = structure
        if not momentum == None:
            assert isinstance(momentum,Vector)
            self.momentum = momentum

    def structure_index(self):                      # returns nonzero index in structure vector (assumes there is only one nonzero entry)
        if abs(self.structure.x) > 0:
            return 0
        if abs(self.structure.y) > 0:
            return 1
        if abs(self.structure.z) > 0:
            return 2
        
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
        
    def set_name_gpt(self,name):
        self.name_gpt = name

    def namestring_gpt_rho(self,type = "minus"):                   # treat this vectorField object as \rho^(+/-) meson. Create according namestring for set_name_gpt
        idx = self.structure_index()
        name_gpt = "FACTOR " + str(np.real(self.structure.x+self.structure.y+self.structure.z)) + " " + str(np.imag(self.structure.x+self.structure.y+self.structure.z)) + "\n"
        if type == "minus":                                         # pi^- 
            name_gpt += "UBAR t\n"
            name_gpt += "GAMMA " + str(idx)        
            name_gpt += "\nMOM [" +str(self.momentum.x) + "," +  str(self.momentum.y) + "," +  str(self.momentum.z) + "] t\n" 
            name_gpt += "D t\n"
        if type == "plus":                                         # pi^- 
            name_gpt += "DBAR t\n"
            name_gpt += "GAMMA" + str(idx)
            name_gpt +="\nMOM [" + str(self.momentum.x) + "," +  str(self.momentum.y) + "," +  str(self.momentum.z) + "] t\n" 
            name_gpt += "U t\n"
        #if zero pion: ...
        return name_gpt
    
    def update_name_gpt(self):                              # updates name_gpt; use in update() only if name_gpt exists
        # gpt naming in case of rho^+/-
        if "UBAR " in self.name_gpt and "D " in self.name_gpt:
            self.set_name_gpt(self.namestring_gpt_rho(type = "minus"))
        elif "DBAR " in self.name_gpt and "U " in self.name_gpt:
            self.set_name_gpt(self.namestring_gpt_rho(type = "plus"))
        # elif zero rho meson: ...


class PseudoVectorField:                               # particle with: Lorentz pseudovector structure (under rotations), momentum as argument
    def __init__(self,momentum,struc = None,modified_momentum_trafo = False):
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
        zero = Vector([0,0,0])
        if self.momentum.is_equal_to(zero):
            self.direction_action = "left"
        else: 
            self.direction_action = "mixed"
        self.modified_momentum_trafo = modified_momentum_trafo
        if self.modified_momentum_trafo:             # parameter. If False: V(p)->RV(R^-1 p), if True: \psi(p)->RV(Rp) under action R
            self.direction_action = "left"
        self.update()    
        
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        modified_mom_trafo = (self.modified_momentum_trafo == True)
        new = PseudoVectorField(mom,struc,modified_momentum_trafo=modified_mom_trafo)
        return new
    
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "PseudoVectorField{" + self.structure.name[4:] + "}(p=" + self.momentum.name[3:] + ")"

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def action(self,A):
        self.structure.action(A)
        if not self.modified_momentum_trafo:
            self.momentum.inverse_action(A)
        else:
            self.momentum.action(A)
        self.update()

    def inverse_action(self,A):
        self.structure.inverse_action(A)
        if not self.modified_momentum_trafo:
            self.momentum.action(A)
        else:
            self.momentum.inverse_action(A)
        self.update()

    def set_vals(self,structure = None,momentum = None):
        if not structure == None:
            assert isinstance(structure,PseudoVector)
            self.structure = structure
        if not momentum == None:
            assert isinstance(momentum,Vector)
            self.momentum = momentum

    def structure_index(self):                      # returns nonzero index in structure vector (assumes there is only one nonzero entry)
        if abs(self.structure.x) > 0:
            return 0
        if abs(self.structure.y) > 0:
            return 1
        if abs(self.structure.z) > 0:
            return 2
        
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum
        

class TensorField:                      # Field operator(with momentum) that transforms like the spacial components of a rank-2 Lorentz-tensor T: T^{mu,nu} -> Lambda^{rho}_mu Lambda^{sigma}_nu T^{rho,sigma} 
                                        # not yet used for anything. Lacks other fields' functionalities: direction_action,modified_momentum_trafo
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


class DiracField:                               # particle with: inner DiracSpinor structure, momentum as argument
    def __init__(self,momentum,struc = None, modified_momentum_trafo = False,Weyl_indices = ["L","R"],LR_different = False):        
        if isinstance(momentum,Vector):
            self.momentum = momentum
        else:
            assert len(momentum) == 3
            self.momentum = Vector(momentum)
        if struc == None:
            self.structure = DiracSpinor(1,0,1,0,Weyl_indices=Weyl_indices,LR_different=LR_different)
        elif isinstance(struc,DiracSpinor):
            self.structure = struc
        else:
            assert len(struc) == 4
            self.structure = DiracSpinor(*struc)
        zero = Vector([0,0,0])
        if self.momentum.is_equal_to(zero):
            self.direction_action = "left"
        else: 
            self.direction_action = "mixed"
        self.modified_momentum_trafo = modified_momentum_trafo
        if self.modified_momentum_trafo:             # parameter. If False: V(p)->G1p(R)V(R^-1 p), if True: \psi(p)->G1p(R)V(Rp) under action R
            self.direction_action = "left"
        self.lorentz_structure = "spinor"
        self.update()
        
    def copy(self):
        mom = self.momentum.copy()
        struc = self.structure.copy()
        modified_mom_trafo = (self.modified_momentum_trafo == True)
        new = DiracField(mom,struc,modified_momentum_trafo=modified_mom_trafo)
        return new
    
    def update(self):
        self.structure.update()
        self.momentum.update()
        self.name = "DiracField{" + self.structure.name + "}(p=" + self.momentum.name[3:] + ")"

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,R):
        return self.structure.is_equal_to(R.structure) and self.momentum.is_equal_to(R.momentum)
    
    def is_negative_of(self,R):
        return self.structure.is_negative_of(R.structure) and self.momentum.is_equal_to(R.momentum)
        
    def action(self,A):
        self.structure.action(A)
        if not self.modified_momentum_trafo:
            self.momentum.inverse_action(A)
        else:
            self.momentum.action(A)
        self.update()

    def inverse_action(self,A):
        self.structure.inverse_action(A)
        if not self.modified_momentum_trafo:
            self.momentum.action(A)
        else:
            self.momentum.inverse_action(A)
        self.update()
    
    def projection(self,R):
        zero = Vector([0,0,0])
        if self.momentum.is_equal_to(zero) and R.momentum.is_equal_to(zero):
            return self.structure.projection(R.structure)
        else:
            if not self.momentum.is_equal_to(R.momentum):
                return 0
            return self.structure.projection(R.structure)
    
    def __lt__(self,obj):
        if self.structure < obj.structure:
            return True
        elif self.structure > obj.structure: 
            return False
        else:
            return self.momentum < obj.momentum


class TensorProduct:                            # tensor product of arbitrary number and types of classes
                                                # each attribute of a TensorProduct instance is one object, e.g. a VectorField
    
    def __init__(self,*v,distinguishable = True):                   # parameter distinguishable: for distinguishable or indistinguishable particles
        self.obj = {}                   
        for i in range(len(v)):
            self.obj[i] = v[i]
        self.distinguishable = distinguishable

        def action_direction_test(self): 
            v_excl_scalars = []
            for o in self.obj.values():
                if not (isinstance(o,Scalar) or isinstance(o,PseudoScalar)):
                    v_excl_scalars.append(o)
            if not all(v_excl_scalars[i].direction_action == v_excl_scalars[0].direction_action for i in range(len(v_excl_scalars))):
                print("WARNING: TensorProduct with mixed action_direction")
        action_direction_test(self)
        self.direction_action = self.obj[0].direction_action
        self.update()

    def copy(self):
        temp = list(self.obj.values())
        new_obj = []
        for i in range(len(temp)):
            new_obj.append(temp[i].copy())
        dist = (self.distinguishable == True)
        new_t = TensorProduct(*new_obj,distinguishable = dist)
        if hasattr(self,"name_gpt"):
            new_t.set_name_gpt(self.name_gpt)
        return new_t
    
    def update(self):
        self.obj[0].update()
        self.name = self.obj[0].name
        for i in range(len(self.obj)-1):
            self.obj[i+1].update()
            self.name += "_x_" + self.obj[i+1].name
        if self.distinguishable:
            self.name+= "[dist.]"
        else:
            self.name += "[indist.]"
        if hasattr(self,"name_gpt"):
            self.update_name_gpt()

    def printname(self):
        self.update()
        print(self.name)

    def is_equal_to(self,t):                            # for distinguishable: True if each object (in order) is exactly the same or opposite, and n_antiparallel is even
                                                        # for indistinguishable: True if objects can be paired to equal or antiparallel pairs, and n_antiparallel is even 
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
        else:                                   # self is indistinguishable product space
            pairs = self.pair_up(t)
            if pairs == None: 
                return False
            return pairs[1] % 2 == 0            # pairs[1]: number of pairs with differing signs. (-1)^n cancels out for n even         
           
    def is_negative_of(self,t):                 # for : True if each object (in order) is exactly the same or opposite, and n_antiparallel is odd
                                                # for indistinguishable: True if objects can be paired to equal or antiparallel pairs, and n_antiparallel is odd 
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
        else:                                   # self is indistinguishable product space
            pairs = self.pair_up(t)
            if pairs == None: 
                return False
            return pairs[1] % 2 == 1            # pairs[1]: number of pairs with differing signs. (-1)^n = -1 for n odd
    
    def action(self,A):
        for o in self.obj.values():
            o.action(A)
        self.update()

    def inverse_action(self,A):
        for o in self.obj.values():
            o.inverse_action(A)
        self.update()

    def projection(self,t2):
        if self.distinguishable and t2.distinguishable:
            proj = 1
            for i in range(len(self.obj)):
                if hasattr(self.obj[i],"projection") and hasattr(t2.obj[i],"projection"):
                    proj = proj*self.obj[i].projection(t2.obj[i])
                else:
                    if self.obj[i].is_equal_to(t2.obj[i]):
                        pass
                    elif self.obj[i].is_negative_of(t2.obj[i]):
                        proj = -1*proj  
                    else:
                        return 0
            return proj   
        else:
            if any(hasattr(self.obj[i],"projection") for i in range(len(self.obj))):
                raise ValueError("Necessary .projection not implemented for indistinguishable <TensorProduct>.")
            if self.is_equal_to(t2):
                return 1
            if self.is_negative_of(t2):
                return -1
            return 0            
        
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
                o2 = negative_match_in_list(o1,objects2)
                if o2 == None:
                    return None                 # failure to pair up at all: Return None
                else: 
                    n_sign_diff += 1
            pairs.append([o1,o2])
            objects2.remove(o2)
        return [pairs,n_sign_diff]
    
    def __lt__(self,t):                         # only for sorting purposes
        if isinstance(self.obj[0],PseudoScalarField):
            return self.obj[0].momentum < t.obj[0].momentum
        return self.obj[0] < t.obj[0]
    
    def set_name_gpt(self,name = None):                             # generic function to set a gpt conform name. Specific functions below.
        if name == None: 
            self.name_gpt = ""
            for i in range(len(self.obj.values())):
                assert hasattr(self.obj[i],"name_gpt")
                self.name_gpt += self.obj[i].name_gpt + "\n"
        else:
            self.name_gpt = name        

    def namestring_gpt_twopion(self):
        assert len(self.obj) == 2
        if hasattr(self.obj[0], "name_gpt") and hasattr(self.obj[1], "name_gpt"):
            cut_name_gpt_0 = self.obj[0].name_gpt.split("\n",1)[1]
            cut_name_gpt_1 = self.obj[1].name_gpt.split("\n",1)[1]
            factor_gpt = self.obj[0].structure.num*self.obj[1].structure.num
            str_factor_gpt = "FACTOR " + str(np.real(factor_gpt)) + " " + str(np.imag(factor_gpt)) + "\n"
            total_name_gpt = str_factor_gpt + cut_name_gpt_0 + cut_name_gpt_1
        return total_name_gpt
    
    def update_name_gpt(self):
        if len(self.obj) == 2:                          
            if hasattr(self.obj[0],"name_gpt") and hasattr(self.obj[1],"name_gpt"):
                if ("UBAR " in self.obj[0].name_gpt and "D " in self.obj[0].name_gpt and "DBAR " in self.obj[1].name_gpt and "U " in self.obj[1].name_gpt) or ("UBAR " in self.obj[1].name_gpt and "D " in self.obj[1].name_gpt and "DBAR " in self.obj[0].name_gpt and "U " in self.obj[0].name_gpt):
                    self.obj[0].update_name_gpt()
                    self.obj[1].update_name_gpt()
                    self.set_name_gpt(self.namestring_gpt_twopion())
    

class LinearCombination:    
    def __init__(self,basis,weights,linear_dep = True,label = None):                       # list of weights of each element in basis; basis: list of objects, e.g. of class <Vector>; 
                                                                            # ! all objects in basis are assumed to be linearly independent
        self.basis = basis
        self.linear_dep = linear_dep                                            # True: treat -obj and obj as linearly dependent; False: treat as linearly independent
        self.lin_comb = []
        if not label == None:
            self.set_label(label)
        assert len(basis) == len(weights)
        weights = np.array(weights)
        for i in range(len(basis)):
            if hasattr(weights[i],"__len__"):
                self.lin_comb.append(TensorProduct(Scalar(weights[i][0]),basis[i]))
            else:
                self.lin_comb.append(TensorProduct(Scalar(weights[i]),basis[i]))
        self.update()

    def copy(self):
        basis = []
        weights = []
        for i in range(len(self.basis)):
            basis.append(self.basis[i].copy())
            weights.append(self.lin_comb[i].obj[0].num)
        lin_dep = (self.linear_dep == True)
        new_LC = LinearCombination(basis,weights,linear_dep = lin_dep)
        if hasattr(self,"label"):
            l = "" + self.label
            new_LC.set_label(l)
        return new_LC
    
    def update(self):                   # updates order of summands and the name
        ordered_lin_comb = [0 for l in self.lin_comb]
        for l in self.lin_comb:                                         # l.obj: [0]: Scalar(weight), [1]: element of basis
            l.update()
            s = match_in_list(l.obj[1],self.basis,index = True)
            if s == None:
                assert self.linear_dep
                s = negative_match_in_list(l.obj[1],self.basis,index = True)
                l.obj[0].num = minus(l.obj[0].num)
                l.obj[1] = self.basis[s]
            assert not s == None
            ordered_lin_comb[s] = l
        self.lin_comb = ordered_lin_comb
        for l in self.lin_comb:
            l.update()
        ########## naming ###########
        relevant_factors = self.nonzero_factors()
        self.name = "Linear_Comb.{"
        for i in range(len(relevant_factors)):
            self.name += relevant_factors[i].obj[0].name[6:] + "*" + relevant_factors[i].obj[1].name + "+"
        self.name =self.name[:-1] + "}"

    def printname(self):
        self.update()
        print(self.name)
    
    def is_equal_to(self,LC):
        f = self.lin_factor(LC)
        if f == None:
            return False
        return abs(f-1) < num_tol
    
    def is_negative_of(self,LC):
        f = self.lin_factor(LC)
        if f == None:
            return False
        return abs(f+1) < num_tol
    
    def action(self,A):
        temp = self.copy()
        for l in temp.lin_comb:
            l.action(A)
        for i in range(len(self.lin_comb)):
            self.lin_comb[i] = temp.lin_comb[i]
        self.update()

    def inverse_action(self,A):
        temp = self.copy()
        for l in temp.lin_comb:
            l.inverse_action(A)
        for i in range(len(self.lin_comb)):
            self.lin_comb[i] = temp.lin_comb[i]
        self.update()
    
    def set_label(self,label):                  # for labelling, e.g. to call some linear combination "x" 
        self.label = label
    def set_name_gpt(self, name = None):
        if name == None:
            if hasattr(self.lin_comb[0].obj[1], "name_gpt"):                # concatenates names of constituents, reevaluates and rewrites first line "FACTOR .."
                name_gpt = ""
                relevant_factors = self.nonzero_factors()
                for i in range(len(relevant_factors)):
                    name_gpt += "FACTOR " + str(np.real(relevant_factors[i].obj[0].num)) + " " + str(np.imag(relevant_factors[i].obj[0].num)) + "\n"
                    TwoPi_name_part = relevant_factors[i].obj[1].name_gpt.split("\n",1)[1]
                    name_gpt += TwoPi_name_part + "\n"
                self.name_gpt = name_gpt

    def lin_factor(self,lc2):                 #returns scalar factor k=s/r such that k*self = lc2; if not linearly dependent in that way, return None
        R = []
        s = []
        n_0 = 0                                 #counter of how many wheights are zero
        k = []
        for i in range(len(self.lin_comb)):
            R.append(self.lin_comb[i].obj[0].num)
            s.append(lc2.lin_comb[i].obj[0].num)
            if abs(R[i]) < num_tol:
                if abs(s[i]) < num_tol:
                    n_0 += 1
                else:
                    return None
            else: 
                k.append(s[i]/R[i])
        if n_0 == len(self.lin_comb):
            print("LC.lin:factor: all wheights zero for both objects")#
            return 1
        for j in range(len(k)):
            if not abs(k[0]-k[j]) < num_tol:
                # print("Exit 1")
                return None
        return k[0]
    
    def nonzero_factors(self):                      #returns [self.lin_comb entries with nonzero weights]         
        nz_pairs = []   
        for x in self.lin_comb:
            if abs(x.obj[0].num) > num_tol:
                nz_pairs.append(x)
        return nz_pairs
    
    def add(self,LC2):
        for i in range(len(self.basis)):
            assert self.basis[i].is_equal_to(LC2.basis[i])
            self.lin_comb[i].obj[0].num+=LC2.lin_comb[i].obj[0].num
        self.update()
    
    def subtract(self,LC2):
        temp = LC2.negative()
        for i in range(len(self.basis)):
            assert self.basis[i].is_equal_to(temp.basis[i])
            self.lin_comb[i].obj[0].num+=temp.lin_comb[i].obj[0].num
        self.update()

    def negative(self):
        neg = self.copy()
        for i in range(len(neg.lin_comb)):
            neg.lin_comb[i].obj[0].num = (-1)* neg.lin_comb[i].obj[0].num 
        return neg
    


    

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

def negative_pair(m,n):                 #returns True if m = -n up to some precision
    return abs(m+n) < num_tol
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
def match_in_list(x,l_y, index = False):                    # returns None if x is not equal up to numerical error to any value in l_y, otherwhise returns first match
    if hasattr(x, "is_equal_to"):
        for y in l_y:
            if x.is_equal_to(y):
                if index:
                    return l_y.index(y)
                return y
    else:
        for y in l_y:
            if x == y:
                if index: 
                    return l_y.index(y)
                return y
            if abs(x-y) < num_tol:
                if index: 
                    return l_y.index(y)
                return y
    return None
def negative_match_in_list(x,l_y, index = False):                    # returns None if -x is not equal up to numerical error to any value in l_y, otherwhise returns first match
    if hasattr(x, "is_negative_of"):
        for y in l_y:
            if x.is_negative_of(y):
                if index:
                    return l_y.index(y)
                return y
    else:
        for y in l_y:
            if abs(x+y) < num_tol:
                if index:
                    return l_y.index(y)
                return y
    return None
def orbit(object,group):                        #applies action of all Group.elements to object and returns list of resulting objects without duplicates
    if object == None:
        return None
    if not hasattr(object,"lorentz_structure"):
        return true_orbit(object,group)         
    else:
        return extended_orbit(object,group)
    # return orbit
def intersection(l_x,l_y):                      # takes two lists. returns as list: intersection (set of common elements) of two lists of objects; returns None if disjoint 
    intersec = []
    for x in l_x:
        if not match_in_list(x,l_y) == None:            # since match_in_list: only adds FIRST, i.e. ONE instance of object from both sets to intersec; DONT USE this if number of occurrance is important
            intersec.append(x)
    if len(intersec) == 0:
        return None
    return intersec

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
def true_orbit(object,group):
    orbit = [object.copy()]
    for g in group.elements:
        o = object.copy()
        o.action(g)
        if not any([q.is_equal_to(o) for q in orbit]):
            orbit.append(o)
    return orbit
def extended_orbit(object,group):                   # find all possible combinations of orbits of each of the attributes; attributes currently: structure and momentum
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
        
##### functions ######

def find_closure(starting_set, object):       #takes list of Strings starting_set, some object instance, and returns dict of the closure regarding the objects resulting from combinations of actions  
    actions = {}    
    for a in starting_set:
        actions[a] = object.copy().action(a)    
    actions["I"] = object.copy()
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
                    new_name = t.remove_I(a + key)
                    known_objects[new_name] = new_o
                    actions[new_name] = new_o
    return actions

def mult_table_from_actions(actions):                             
    obj = actions["I"].copy()
    mt = {}
    for a in actions.keys():                            
        mt[a] = {}   
        for b in actions.keys():
            A = t.remove_I(a+b)
            w = obj.copy()
            w.action(A)                                 #result is w = Ae = abe for reference object e = obj   
            for c in actions.keys():
                x = obj.copy()
                x.action(c)
                if x.is_equal_to(w):
                    mt[a][b] = c                        # a.b = c because they have the same action            
    return mt

def generate_group(gen_actions,obj):
    print("Generating group.")
    actions = find_closure(gen_actions,obj)
    print("find_closure done.")
    mt = mult_table_from_actions(actions)
    print("mult_table done.")
    G = g.Group(list(actions.keys()),mt)
    print("Group generated. # elements: ", len(list(actions.keys())))
    return G

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














 




