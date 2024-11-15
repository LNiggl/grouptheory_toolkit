import numpy as np
import groups as g
import representations as r
import objects as o
import testing as t
np.set_printoptions(precision = 6, suppress = True)
# group O_h

gen_actions = ["Rot0","Rot1","Rot2","Inv"]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O_h = o.generate_group(gen_actions,v)

# irreps of O_h
a = [o.Scalar(1)]
A1p = r.rep_from_action(O_h,a,"A1p")

b = o.generate_basis([o.Vector([1,0,0])],O_h)
T1m = r.rep_from_action(O_h,b,"T1m")
T1m.check_if_homomorphism()

ps = o.PseudoScalar(1)
b_ps = o.generate_basis([ps],O_h)
A1m = r.rep_from_action(O_h,b_ps,"A1m")

T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
T1m_T1m_reps = T1m_x_T1m.hom.copy()
T1m_x_T1m.check_if_homomorphism()

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
Ep.round_chars()

Em = r.product_rep(Ep,A1m)
Em.name = "Em"
Em.check_if_homomorphism()
Em.round_chars()

A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
A2m.name = "A2m"
r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
A2m.check_if_homomorphism()
A2m.round_chars()

A2p = r.product_rep(A2m,A1m)
A2p.name = "A2p"
A2p.check_if_homomorphism()
A2p.round_chars()

O_h.set_char_table([A1m,A1p,T1m,T1p,T2m,T2p,Em,Ep,A2m,A2p])
print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": " , O_h.char_table[irrep])

# dummy instances
s1 = o.Scalar(1)
s2 = o.Scalar(-1)

ps1 = o.PseudoScalar(1)
ps2 = o.PseudoScalar(-1)

v1 = o.Vector([1,0,0])
v2 = o.Vector([-1,0,0])
v3 = o.Vector([2,0,0])
ten_v = o.TensorProduct(v1,v2)

pv1 = o.PseudoVector([1,0,0])
pv2 = o.PseudoVector([-1,0,0])

sf1 = o.ScalarField([0,0,0])
sf2 = o.ScalarField([1,0,0])
sf3 = o.ScalarField([-1,0,0])
ten_sf = o.TensorProduct(sf2,sf3)

psf1 = o.PseudoScalarField([0,0,0],struc = 1)
psf2 = o.PseudoScalarField([1,0,0], struc = 1)
psf3 = o.PseudoScalarField([-1,0,0],struc = 1)
ten_psf = o.TensorProduct(psf2,psf3,distinguishable=False)

vf1 = o.VectorField([0,0,0],struc = [1,0,0])
vf2 = o.VectorField([1,0,0],struc = [1,0,0])
vf3 = o.VectorField([-1,0,0],struc = [1,0,0])
ten_vf = o.TensorProduct(vf2,vf3)

pvf1 = o.PseudoVectorField([0,0,0],struc = [1,0,0])
pvf2 = o.PseudoVectorField([1,0,0],struc = [1,0,0])
pvf3 = o.PseudoVectorField([-1,0,0],struc = [1,0,0])
ten_pvf = o.TensorProduct(pvf2,pvf3)

object_list = [s1,s2,ps1,ps2,v1,v2,pv1,pv2,sf1,sf2,psf1,psf2,vf1,vf2,pvf1,pvf2]

v_basis = o.generate_basis(v1,O_h)
v_rep = r.rep_from_action(O_h,v_basis,"v_rep")
v_rep.check_if_homomorphism()
print("V_Rep Rot0:")
print(v_rep.hom["Rot0"])
print("V_Rep Rot1:")
print(v_rep.hom["Rot1"])
print("V_Rep Rot2:")
print(v_rep.hom["Rot2"])

ps_basis = o.generate_basis(ps1,O_h)
ps_rep = r.rep_from_action(O_h,ps_basis,"ps_rep")
ps_rep.check_if_homomorphism()

pv_basis = o.generate_basis(pv1,O_h)
pv_rep = r.rep_from_action(O_h,pv_basis,"pv_rep")
pv_rep.check_if_homomorphism()

print("scalarfield rep with right action")
print("action direction: ", sf2.direction_action)
sf_basis = o.generate_basis(sf2,O_h)
sf_rep = r.rep_from_action(O_h,sf_basis,"sf_rep")
sf_rep.check_if_homomorphism()

# psf_basis = o.generate_basis(psf1,O_h)
# psf_rep = r.rep_from_action(O_h,psf_basis,"psf_rep")
# psf_rep.check_if_homomorphism()

# vf_basis = o.generate_basis(vf1,O_h)
# vf_rep = r.rep_from_action(O_h,vf_basis,"vf_rep")
# vf_rep.check_if_homomorphism()

# pvf_basis = o.generate_basis(pvf1,O_h)
# pvf_rep = r.rep_from_action(O_h,pvf_basis,"pvf_rep")
# pvf_rep.check_if_homomorphism()
# print("tensor product reps")
# b_ten_v = o.generate_basis(ten_v,O_h)
# rep_ten_v = r.rep_from_action(O_h,b_ten_v,"rep_ten_v")
# print("action direction of ",rep_ten_v.name , ": " , b_ten_v[0].direction_action)
# rep_ten_v.check_if_homomorphism()

# b_ten_sf = o.generate_basis(ten_sf,O_h)
# rep_ten_sf = r.rep_from_action(O_h,b_ten_sf,"rep_ten_sf")
# print("action direction of ",rep_ten_sf.name , ": " , b_ten_sf[0].direction_action)
# rep_ten_sf.check_if_homomorphism()

# b_ten_psf = o.generate_basis(ten_psf,O_h)
# rep_ten_psf = r.rep_from_action(O_h,b_ten_psf,"rep_ten_psf")
# rep_ten_psf.check_if_homomorphism()

# b_ten_vf = o.generate_basis(ten_vf,O_h)
# rep_ten_vf = r.rep_from_action(O_h,b_ten_vf,"rep_ten_vf")
# # rep_ten_vf.check_if_homomorphism2()

# b_ten_pvf = o.generate_basis(ten_pvf,O_h)
# rep_ten_pvf = r.rep_from_action(O_h,b_ten_pvf,"rep_ten_pvf")
# print("action direction of ",rep_ten_pvf.name , ": " , b_ten_pvf[0].direction_action)
# rep_ten_pvf.check_if_homomorphism()

# test associativity axiom for group action on PseudoScalarField
A = "Rot1"
B = "Rot2"


# print("test associativity of PseudoscalarField_x_PseudoscalarField action under ", A , B, " and the product ", O_h.mult_table[A][B])
# print("subsequent actions")
# ten_psf_c = ten_psf.copy()
# ten_psf.printname()
# ten_psf.action(B)
# print("->")
# ten_psf.printname()
# ten_psf.action(A)
# print("->")
# ten_psf.printname()
# print("action of (", A,B, "):")
# ten_psf_c.printname()
# ten_psf_c.action(O_h.mult_table[A][B])
# print("->")
# ten_psf_c.printname()

# print("vector_x_vector actions to build rank-2 lorentz tensor:")
# A = "Rot0"
t0 = o.Vector([1,0,0])
t1 = o.Vector([0,1,0])
t2 = o.Vector([0,0,1])

T00 = o.TensorProduct(t0,t0)
T10 = o.TensorProduct(t1,t0)
T20 = o.TensorProduct(t2,t0)

T01 = o.TensorProduct(t0,t1)
T11 = o.TensorProduct(t1,t1)
T21 = o.TensorProduct(t2,t1)

T02 = o.TensorProduct(t0,t2)
T12 = o.TensorProduct(t1,t2)
T22 = o.TensorProduct(t2,t2)

list_t = [T00,T01,T02,T10,T11,T12,T20,T21,T22]


print("Calculations for basis trafo to book (Altmann Herzig) notation: T1m")
print("Rot0:")
print(T1m.hom["Rot0"])
print("eig:")
evR, evecsR = np.linalg.eig(T1m.hom["Rot0"])
print(evR)
for i in range(len(evecsR)):
    print(evecsR[:,i])
C_4x = np.array([[0,-1j,0],[-1j,0,0],[0,0,1]])
print("C_4x")
print(C_4x)
evC,evecsC = np.linalg.eig(C_4x)
print(evC)
for i in range(len(evecsC)):
    print(evecsC[:,i])
print("G:")
print(evecsR)
print("H:")
print(evecsC)
S = evecsR@np.linalg.inv(evecsC)
print("S = GH^-1")
print(S)
print("S^-1")
print(np.linalg.inv(S))
print("Test S^-1AS-C4x:")
print(np.linalg.inv(S)@T1m.hom["Rot0"]@S-C_4x)

print("other two rotations:")
C_4y = np.array([[1,0,0],[0,0,1],[0,-1,0]])
C_4z = np.array([[0,0,-1j],[0,1,0],[-1j,0,0]])
print("Test S^-1Rot_yS-C_4y:")
print(np.linalg.inv(S)@T1m.hom["Rot1"]@S-C_4y)

print("Test S^-1Rot_yS-C_4z:")
print(np.linalg.inv(S)@T1m.hom["Rot2"]@S-C_4z)

print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": ")
    for c,mems in O_h.classes.items():
        print(c," (" , len(mems) , "): ", O_h.char_table[irrep][c])











