import numpy as np
import groups as g
import representations as r
import objects as o
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

T1p = r.Representation(O_h,T1m_T1m_reps,"T1p")
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

vp1 = o.VectorField([1,0,0],[1,0,0])
vp2 = o.VectorField([0,1,0],[0,1,0])
vp3 = o.VectorField([0,0,1],[0,0,1])
print("Basis")
vpo = o.generate_basis([vp1,vp2,vp3],O_h,drop_sign = True)
o.print_all(vpo)
print("# ", len(vpo))
exb = o.exhaustive_orbit(vp1,O_h)
exb.sort(reverse = True)
exb2 = o.remove_neg(exb)
o.print_all(exb2)
print("# ", len(exb2))

vv = o.Pion([1,0,0],"-")
exv = o.exhaustive_orbit(vv,O_h)
exv.sort(reverse = True)
o.print_all(exv)
print("# ", len(exv))

vp1 = o.VectorField([1,0,0],[1,0,0])
vpp = o.generate_basis(vp1,O_h,drop_sign = True)
o.print_all(vpp)

vvv = o.Vector([1,0,0])
vv2 = o.Vector([-1,0,0])
print(vvv.is_negative_of(vv2))
vvx = o.generate_basis(vvv,O_h)
o.print_all(vvx)

pv = o.PseudoVectorField([1,0,0])
pv2 = o.PseudoVectorField([1,0,0],struc = [-1,0,0])
o.print_all([pv,pv2])
print(pv.structure.name, pv2.structure.name, pv.momentum.name,pv2.momentum.name)

print(pv.is_negative_of(pv2))

# dummy instances
s1 = o.Scalar(1)
s2 = o.Scalar(-1)

ps1 = o.PseudoScalar(1)
ps2 = o.PseudoScalar(-1)

v1 = o.Vector([1,0,0])
v2 = o.Vector([-1,0,0])
v3 = o.Vector([2,0,0])

pv1 = o.PseudoVector([1,0,0])
pv2 = o.PseudoVector([-1,0,0])

sf1 = o.ScalarField([1,0,0])
sf2 = o.ScalarField([-1,0,0])

psf1 = o.PseudoScalarField([1,0,0],struc = 1)
psf2 = o.PseudoScalarField([-1,0,0], struc = 1)

vf1 = o.VectorField([1,0,0],struc = [1,0,0])
vf2 = o.VectorField([1,0,0],struc = [-1,0,0])

pvf1 = o.PseudoVectorField([1,0,0],struc = [1,0,0])
pvf2 = o.PseudoVectorField([1,0,0],struc = [-1,0,0])

object_list = [s1,s2,ps1,ps2,v1,v2,pv1,pv2,sf1,sf2,psf1,psf2,vf1,vf2,pvf1,pvf2]

ten1 = o.TensorProduct(v1,v2,distinguishable=False)
ten2 = o.TensorProduct(v2,v2,distinguishable=False)
ten3 = o.TensorProduct(v2,v1)
ten4 = o.TensorProduct(v2,v1)

v_basis = o.generate_basis(v1,O_h)
v_rep = r.rep_from_action(O_h,v_basis,"testrep")
print(v_rep.hom["Rot2"])
v_rep.check_if_homomorphism()

ten1.printname()
print("test orbit")
oten = o.orbit(ten1,O_h)
o.print_all(oten)
print("orbit of ", ten3.name)
oten3 = o.orbit(ten3,O_h)
o.print_all(oten3)
b_ten3 = o.generate_basis(ten3,O_h)
print("basis of ten3")
o.print_all(b_ten3)

ten_psf = o.TensorProduct(psf1,psf2)
b_psf = o.generate_basis(ten_psf,O_h)
print("basis of ten_sf")
o.print_all(b_psf)
psf_rep = r.rep_from_action(O_h,b_psf,"pseudoscalarfield_rep")
psf_rep.check_if_homomorphism()
# print(psf_rep.hom["Rot2"])
# print(np.matmul(psf_rep.hom["Rot2"],psf_rep.hom["Rot1"]))
# print(O_h.mult_table["Rot1"]["Rot2"])
# print(psf_rep.hom["Rot2Rot0"])
# psf_rep.check_if_homomorphism2()
r.study_irreps(psf_rep,O_h,"../results/PseudoScalarField_100_irreps.txt")

psf3 = o.PseudoScalarField([1,1,0],struc = 1)
psf4 = o.PseudoScalarField([-1,-1,0], struc = 1)

# ten_psf2 = o.TensorProduct(psf3,psf4)
# b_psf2 = o.generate_basis(ten_psf2,O_h)
# print("basis of ten_sf2")
# o.print_all(b_psf2)
# psf_rep2 = r.rep_from_action(O_h,b_psf2,"pseudoscalarfield_rep")
# print(psf_rep2.hom["Inv"])
# psf_rep2.check_if_homomorphism()
# r.study_irreps(psf_rep2,O_h,"../results/PseudoScalarField_110_irreps.txt")





# print("test is_equal")
# print(ten3.is_equal_to(ten4))
# print(ten1.is_equal_to(ten2))
# print(ten1.is_negative_of(ten2))
# pairs = ten1.pair_up(ten2)
# print("pairs")
# if pairs == None:
#     print("none")
# else:
#     for p in pairs[0]:
#         o.print_all(p)
#         print("\n")
#     print("n_sign_diff ", pairs[1])


# for x in object_list:
#     x.printname()

# v = o.Vector([1,0,0])
# b_v = o.generate_basis([v],O_h) 
# V_Rep = r.rep_from_action(O_h,b_v,"V_Rep")
# # print("General Vector Basis")
# # o.print_all(V_Rep.basis)
# # print("# ",len(V_Rep.basis))
# V_x_V_Rep =r.product_rep(V_Rep,V_Rep)
# print("product basis VxV")
# o.print_all(V_x_V_Rep.basis.values())
# V_red = r.find_irreps(V_x_V_Rep,O_h)
# print("diag traceless projector:")
# P = np.matmul(r.traceless_projector(3),r.diagonal_projector(3))
# print(P)
# P_V_x_V_eigvecs = r.list_nonzero_eigvecs(np.matmul(V_x_V_Rep.hom["Rot2"],P))
# print("EIgenvectors of P:")
# print(P_V_x_V_eigvecs)
# print("after applying P to VxV:")
# # r.apply_projectors([r.traceless_projector,r.diagonal_projector],V_x_V_Rep)
# print("test: rotate A -> B and compare with result: EV of Rot1")
# print("Rotate A via Rot0")
# for k,v in P_V_x_V_eigvecs.items():
#     # print(v)
#     x = np.matmul(v[0].T,V_x_V_Rep.hom["Rot0"]).T
#     print(k,x)
# print("B from EV of Rot1P")
# P_V_x_V_eigvecs_y = r.list_nonzero_eigvecs(np.matmul(V_x_V_Rep.hom["Rot1"],P))
# print(P_V_x_V_eigvecs_y)

# print("Compare to E in two-pion system")
# p1 = o.Pion([1,0,0],"+")
# p2 = o.Pion([-1,0,0],"+")

# # Two Pions, distinguishable
# tp1 = o.Two_Pion(p1,p2)
# b_tp = o.generate_basis([tp1],O_h) 
# print("TP1 Basis: (1,0,0)-type single momenta")
# o.print_all(b_tp)
# print("# ", len(b_tp))
# TP_Rep1 = r.rep_from_action(O_h,b_tp,"TP_Rep1")

# TP_red1 = r.find_irreps(TP_Rep1,O_h)

# # vectors that span Ep:
# P_TP1_Ep = TP_red1["Ep"]
# print("projector to E of Two_Pion")
# print(P_TP1_Ep[0])
# TP1_Ep_vecs = r.list_nonzero_eigvecs(np.matmul(TP_Rep1.hom["Rot2"],P_TP1_Ep[0]))
# print("eigenvectors:")
# print(TP1_Ep_vecs)

# result = r.E_identify_components(P_TP1_Ep,TP_Rep1)
# print(result)

# #vectors that span T2:
# P_T2p = V_red["T2p"]
# o.print_all(V_x_V_Rep.basis.values())
# print("projector to T2p")
# print(P_T2p[0])
# print(V_x_V_Rep.hom["Rot2"])
# print(np.matmul(V_x_V_Rep.hom["Rot2"],P_T2p[0]))
# T2p_eig = r.list_nonzero_eigvecs(P_T2p[0])
# print(T2p_eig)
# T2p_vec_a = r.list_nonzero_eigvecs(np.matmul(V_x_V_Rep.hom["Rot2"],P_T2p[0]))
# print("vector a of T2")
# print(T2p_vec_a)
# print("Rotate a via Rot0 to b")
# for k,v in T2p_vec_a.items():
#     if abs(k+1)<1e-8:
#         x = v[0]
# x = np.matmul(V_x_V_Rep.hom["Rot0"],x)
# print(x)
# print("rotate b again to -a")
# x = np.matmul(V_x_V_Rep.hom["Rot0"],x)
# print(x)
# T2p_vec_b = r.list_nonzero_eigvecs(np.matmul(V_x_V_Rep.hom["Rot1"],P_T2p[0]))
# print("compare to eigenvalue -1 from Roty*projection")
# print(T2p_vec_b)
# vecs = r.T2_identify_components(P_T2p,V_x_V_Rep)
# print(vecs)






