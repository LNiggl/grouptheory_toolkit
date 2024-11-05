import numpy as np
import groups as g
import objects as o
import representations as r
import testing as t
gen_actions = ["Rot0","Rot1","Rot2","Inv"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O_h = o.generate_group(gen_actions,v)

# irreps of O_h

a = [o.Scalar(1)]
A1p = r.rep_from_action(O_h,a,"A1p")
A1p.check_if_homomorphism()

b = o.generate_basis([o.Vector([1,0,0])],O_h)
T1m = r.rep_from_action(O_h,b,"T1m")
T1m.check_if_homomorphism()

ps = o.PseudoScalar(1)
b_ps = o.generate_basis([ps],O_h)
A1m = r.rep_from_action(O_h,b_ps,"A1m")

T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
# T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

# T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
T1p = T1m_x_T1m.copy("T1p")
r.apply_projectors([r.antisymmetric_projector],T1p)
T1p.check_if_homomorphism()

T2p = T1m_x_T1m.copy("T2p")
r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2p)
T2p.check_if_homomorphism()

T2m = r.product_rep(T2p,A1m)
T2m.name = "T2m"
T2m.check_if_homomorphism()

Ep = T1m_x_T1m.copy("Ep")
r.apply_projectors([r.traceless_projector,r.diagonal_projector],Ep)
Ep.check_if_homomorphism()
# Ep.round_chars()

Em = r.product_rep(Ep,A1m)
Em.name = "Em"
Em.check_if_homomorphism()
# Em.round_chars()

A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
A2m.name = "A2m"
r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
A2m.check_if_homomorphism()
# A2m.round_chars()

A2p = r.product_rep(A2m,A1m)
A2p.name = "A2p"
A2p.check_if_homomorphism()


list_irreps = [A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep]
O_h.set_char_table(list_irreps)

### test trafo behavior of components of invariant subspaces 
gen_actions = ["Rot0","Rot1","Rot2","Inv"]
f = open("../tests/trafo_behavior_O_h.txt","w")
## T1m
x = o.LinearCombination(T1m.basis,[1,0,0], label = "x")
y = o.LinearCombination(T1m.basis,[0,1,0], label = "y")
z = o.LinearCombination(T1m.basis,[0,0,1], label = "z")
print("test Rot2")
print("Matrix:\n")
print(T1m.hom["Rot2"])
xx = T1m.basis[0].copy()
print(xx.name)
xx.action("Rot2")
print("->",xx.name)
LC_basis_T1m = [x,y,z]
f.write("T1m\nBasis: \n")
x_weight = [1,0,0]
Rot2_x_weight = np.matmul(T1m.hom["Rot2"],x_weight)
LC_Rot2_x = o.LinearCombination(T1m.basis,Rot2_x_weight,label = "Rot2_x")
for coord in [x,y,z]:
    print("factor between ", LC_Rot2_x.label, " and ", coord.label, ": " , LC_Rot2_x.lin_factor(coord))
for obj in LC_basis_T1m:
    f.write(obj.label + ": " + obj.name)
    f.write("\n")
trafos_T1m = {}
for A in gen_actions:
    trafos_T1m[A] = {}
    for c in LC_basis_T1m:
        temp = c.copy()       
        temp.action(A)
        res = o.match_in_list(temp,LC_basis_T1m)
        if res == None:
            res = o.negative_match_in_list(temp,LC_basis_T1m)
            if res == None: 
                print("Problem")
            trafos_T1m[A][c.label] = o.minus(res.label)
        else:
            trafos_T1m[A][c.label] = res.label
f.write(str(trafos_T1m))
f.write("\n\n")

##T1p test
b_test = o.generate_basis([o.PseudoVector([1,0,0])],O_h)
T1p_test= r.rep_from_action(O_h,b_test,"T1p_test")
T1p_test.check_if_homomorphism()

x_t = o.LinearCombination(T1p_test.basis,[1,0,0], label = "x")
y_t = o.LinearCombination(T1p_test.basis,[0,1,0], label = "y")
z_t = o.LinearCombination(T1p_test.basis,[0,0,1], label = "z")

LC_basis_T1p = [x_t,y_t,z_t]
f.write("T1p\nBasis: \n")
for obj in LC_basis_T1p:
    f.write(obj.label + ": " + obj.name)
    f.write("\n")
trafos_T1p = {}
for A in gen_actions:
    trafos_T1p[A] = {}
    for c in LC_basis_T1p:
        temp = c.copy()       
        temp.action(A)
        res = o.match_in_list(temp,LC_basis_T1p)
        if res == None:
            res = o.negative_match_in_list(temp,LC_basis_T1p)
            if res == None: 
                print("Problem")
            trafos_T1p[A][c.label] = o.minus(res.label)
        else:
            trafos_T1p[A][c.label] = res.label

trafos_T1p = t.test_trafo_behavior(LC_basis_T1p,gen_actions)
f.write(str(trafos_T1p))
f.write("\n\n")

 ## T2p

# print("T2p Basis:")
# o.print_all(T2p.basis)
T2p_P = np.array(r.symmetric_projector(3)*r.invert_projector(r.diagonal_projector)(3))
print("T2p projector:\n")
print(T2p_P)
## per reading off
T2p_weights = [T2p_P[1],T2p_P[2],T2p_P[5]]
tau_1 = o.LinearCombination(T2p.basis,T2p_weights[0], label = "tau_1")
tau_2 = o.LinearCombination(T2p.basis,T2p_weights[1], label = "tau_2")
tau_3 = o.LinearCombination(T2p.basis,T2p_weights[2], label = "tau_3")
print("Rot0 matrix:")
print(T1m_x_T1m.hom["Rot0"])
pi1 = T1m_x_T1m.basis[2].copy()
print(pi1.name)
pi1.action("Rot0")
print("Rot0 ->" , pi1.name)
LC_basis_T2p = [tau_1,tau_2,tau_3]
f.write("T2p\nBasis: \n")
for obj in LC_basis_T2p:
    f.write(obj.label + ": " + obj.name)
    f.write("\n")
trafos_T2p = {}
for A in gen_actions:
    trafos_T2p[A] = {}
    for c in LC_basis_T2p:
        temp = c.copy()    
        temp.action(A)
        res = o.match_in_list(temp,LC_basis_T2p)
        if res == None:
            res = o.negative_match_in_list(temp,LC_basis_T2p)
            if res == None: 
                print("Problem")
            trafos_T2p[A][c.label] = o.minus(res.label)
        else:
            trafos_T2p[A][c.label] = res.label
f.write(str(trafos_T2p))
f.write("\n\n")
## T2m

print("T2m Basis:")
o.print_all(T2m.basis)

# print(T2m_P)
## per reading off
T2m_weights = [T2p_P[1],T2p_P[2],T2p_P[5]]
tau_1 = o.LinearCombination(T2m.basis,T2p_weights[0], label = "tau_1")
tau_2 = o.LinearCombination(T2m.basis,T2p_weights[1], label = "tau_2")
tau_3 = o.LinearCombination(T2m.basis,T2p_weights[2], label = "tau_3")

LC_basis_T2m = [tau_1,tau_2,tau_3]
f.write("T2m\nBasis: \n")
for obj in LC_basis_T2m:
    f.write(obj.label + ": " + obj.name)
    f.write("\n")
trafos_T2m = {}
for A in gen_actions:
    trafos_T2m[A] = {}
    for c in LC_basis_T2m:
        temp = c.copy()   
        temp.action(A)
        res = o.match_in_list(temp,LC_basis_T2m)
        if res == None:
            res = o.negative_match_in_list(temp,LC_basis_T2m)
            if res == None: 
                print("Problem")
            trafos_T2m[A][c.label] = o.minus(res.label)
        else:
            trafos_T2m[A][c.label] = res.label
f.write(str(trafos_T2m))
f.write("\n\n")