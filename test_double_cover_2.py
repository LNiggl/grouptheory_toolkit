# import numpy as np
# from base import groups as g
# from base import objects as o
# from base import representations as r
# from base import testing as t
# from base import O_h,gen_actions_O_h                    # import group list of actions used in generation of O_h

# np.set_printoptions(precision = 8, suppress = True)

import numpy as np
from base import groups as g
from base import objects as o
from base import representations as r
from base import testing as t
from base import O_h,gen_actions_O_h                    # import group list of actions used in generation of O_h

np.set_printoptions(precision = 8, suppress = True)

gen_actions = ["Rot0","Rot1","Rot2"]
gen_actions = ["Rot0","Rot1","Rot2"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O = o.generate_group(gen_actions,v)

# Pauli matrices - convention: 0 = x, 1 = y and 2 = z
sigma = {}
sigma[0] = np.array([[0,1],[1,0]], dtype=np.complex128)
sigma[1] = np.array([[0,-1j],[1j,0]], dtype=np.complex128)
sigma[2] = np.array([[1,0],[0,-1]], dtype=np.complex128)
# generators of double cover of O -  convention: 0 = x, 1 = y and 2 = z
A = {}
A_inv = {}

A[0]= 1/np.sqrt(2)*(np.eye(2) -1j*sigma[0])
A[1] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[1])
A[2] = 1/np.sqrt(2)*(np.eye(2) -1j*sigma[2])

for i in range(3):
    A_inv[i] = np.linalg.inv(A[i])

def R_entry(A,j,k):
    return 1/2*np.trace(sigma[j]@A@sigma[k]@np.linalg.inv(A))
def R(A):
    M = np.zeros((3,3),dtype = np.complex128)
    for j in range(3):
        for k in range(3):
            M[j][k] = R_entry(A,j,k)
    return M

def group_hom_via_reps(Rep1,Rep2,F):        #F:function relating Rep1 matrices to Rep2 matrices, e.g. Rep1: Double cover of O, Rep2: O, F: R(A); returns dict{groupelement1:groupelement2}
    group_hom = {}
    for g1,m1 in Rep1.hom.items():            
        n = R(m1)
        # if len(g1) < 5:
        #     print("In group_hom_via_reps. ", g1, " for", Rep1.name, " : ", n)
        matches = []
        for g2,m2 in Rep2.hom.items():
            if np.allclose(n,m2):
                matches.append(g2)
        assert len(matches) == 1            # any-to-one homomorphism
        group_hom[g1] = matches[0]
    return group_hom

b = o.generate_basis([o.Vector([1,0,0])],O)
T1 = r.rep_from_action(O,b,"T1")
T1.check_if_homomorphism()
gen_actions = ["Rot0","Rot1","Rot2"]
gen_actions_A = ["I","A0","A1","A2"]
s_ab = o.WeylSpinor(1/np.sqrt(2),1/np.sqrt(2))
print(type(s_ab))
# O_from_spinor = o.generate_group(gen_actions,s_ab)

# s_ab = o.WeylSpinor(1/np.sqrt(2),1/np.sqrt(2))
# O_from_spinor = o.generate_group(gen_actions,s_ab)
## G1 from spinor action
# print("haha")
# s_ab = o.WeylSpinor(1/np.sqrt(2),1/np.sqrt(2))
DC_O_actions = o.find_closure(gen_actions_A,s_ab)
print("#: ",len(DC_O_actions))
i = 1
for k,v in DC_O_actions.items():
    for kk in DC_O_actions.keys():
        print("key ", kk)
        print("v: ", v)
        vv = v.action(kk)
        print("vv " , vv)
        if vv.is_equal_to(s_ab):
            print("Inverse to" , k, ": ", kk)
            print(i)
            i+= 1



# DC_O_from_spinor = o.generate_group(gen_actions_A,s_ab)
# print("# DC_O_from_spinor:", len(DC_O_from_spinor.elements))
# print("classes")
# print(DC_O_from_spinor.classes)
# s1 = o.WeylSpinor(1,0)
# s2 = o.WeylSpinor(0,1)
# G1_from_spinor = r.rep_from_action(DC_O_from_spinor,[s1,s2],"G1_from_spinor_trafo")
# G1_from_spinor.check_if_homomorphism()