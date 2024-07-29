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

O_h.set_char_table([A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep])
print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": " , O_h.char_table[irrep])

# write to file 
P_file = open("../results/Pion_irreps.txt", "w")


##(1,0,0) type single momenta

p1 = o.Pion([1,0,0],"+")
p2 = o.Pion([-1,0,0],"+")

# Two Pions, distinguishable
tp1 = o.Two_Pion(p1,p2)
b_tp = o.generate_basis([tp1],O_h) 
P_file.write("TP1 Basis: (1,0,0)-type single momenta"+"\n")
for x in b_tp:
    P_file.write(x.name + "\n")
o.print_all(b_tp)
P_file.write("# " + str(len(b_tp)))
TP_Rep1 = r.rep_from_action(O_h,b_tp,"TP_Rep1")
print(TP_Rep1.hom["Rot2"])
TP_Rep1.check_if_homomorphism()
TP_red1 = r.find_irreps(TP_Rep1,O_h)

# # A1p: 1
# # T1m: 1
# # Ep: 1
# # total dim: 6


# vectors that span A1p:
P_TP1_A1p = TP_red1["A1p"][0]
print(P_TP1_A1p)
TP1_A1p_vecs = r.list_nonzero_eigvecs(P_TP1_A1p)
P_file.write("A1p subspace:"+"\n") 
P_file.write(str(TP1_A1p_vecs)+"\n")

# vectors that span T1m:
P_TP1_T1m = TP_red1["T1m"]
TP1_vector_components = r.T1_identify_components(P_TP1_T1m,TP_Rep1)
P_file.write("(x,y,z)-like eigenvectors in TP1:"+"\n")
P_file.write(str(TP1_vector_components)+"\n")

P_file.write("Ep subspace:"+"\n")
P_TP1_Ep = TP_red1["Ep"]
TP1_Ep_components = r.E_identify_components(P_TP1_Ep,TP_Rep1)
P_file.write(str(TP1_Ep_components)+"\n")
# ## (1,1,0)-type momenta

p3 = o.Pion([1,1,0],"+")
p4 = o.Pion([-1,-1,0],"+")

# Two Pions, distinguishable
tp2 = o.Two_Pion(p3,p4)
b_tp2 = o.generate_basis([tp2],O_h) 
P_file.write("TP2 Basis: (1,1,0)-type single momenta"+"\n")
for x in b_tp2:
    P_file.write(x.name + "\n")
o.print_all(b_tp2)
P_file.write("# " + str(len(b_tp2))+"\n")
TP_Rep2 = r.rep_from_action(O_h,b_tp2,"TP_Rep2"+"\n")

TP_red2 = r.find_irreps(TP_Rep2,O_h)

# # A1p: 1
# # T1m: 1
# # T2m: 1
# # T2p: 1
# # Ep: 1

## vectors spanning subspaces:

P_TP2_A1p = TP_red2["A1p"][0]
P_TP2_T1m = TP_red2["T1m"]


TP2_A1p_vecs = r.list_nonzero_eigvecs(P_TP2_A1p)
P_file.write("A1p subspace:"+"\n") 
P_file.write(str(TP2_A1p_vecs)+"\n")

TP2_vector_components = r.T1_identify_components(P_TP2_T1m,TP_Rep2)
P_file.write("(x,y,z)-like eigenvectors in TP1:"+"\n")
P_file.write(str(TP2_vector_components)+"\n")

P_file.write("Ep subspace:"+"\n")
P_TP2_Ep = TP_red2["Ep"]
TP2_Ep_components = r.E_identify_components(P_TP2_Ep,TP_Rep2)
P_file.write(str(TP2_Ep_components)+"\n")

P_file.write("T2m subspace:" + "\n")
P_TP2_T2m = TP_red2["T2m"]
TP2_T2m_components = r.T2_identify_components(P_TP2_T2m,TP_Rep2)
P_file.write(str(TP2_T2m_components)+"\n")

P_file.write("T2p subspace:" + "\n")
P_TP2_T2p = TP_red2["T2p"]
TP2_T2p_components = r.T2_identify_components(P_TP2_T2p,TP_Rep2)
P_file.write(str(TP2_T2p_components)+"\n")

#(1,1,1) type single momenta

p5 = o.Pion([1,1,1],"+")
p6 = o.Pion([-1,-1,-1],"+")

# Two Pions, distinguishable
tp3 = o.Two_Pion(p5,p6)
b_tp3 = o.generate_basis([tp3],O_h) 
P_file.write("TP3 Basis: (1,1,1)-type single momenta"+"\n")
for x in b_tp3:
    P_file.write(x.name + "\n")
o.print_all(b_tp3)
P_file.write("# " +  str(len(b_tp3))+"\n")
TP_Rep3 = r.rep_from_action(O_h,b_tp3,"TP_Rep3"+"\n")

TP_red3 = r.find_irreps(TP_Rep3,O_h)

P_TP3_A1p = TP_red3["A1p"][0]
P_TP3_T1m = TP_red3["T1m"]

TP3_A1p_vecs = r.list_nonzero_eigvecs(P_TP3_A1p)
P_file.write("A1p subspace:"+"\n") 
P_file.write(str(TP3_A1p_vecs)+"\n")

TP3_vector_components = r.T1_identify_components(P_TP3_T1m,TP_Rep3)
P_file.write("(x,y,z)-like eigenvectors in TP3:"+"\n")
P_file.write(str(TP3_vector_components)+"\n")

P_file.write("T2p subspace:" + "\n")
P_TP3_T2p = TP_red3["T2p"]
TP3_T2p_components = r.T2_identify_components(P_TP3_T2p,TP_Rep3)
P_file.write(str(TP3_T2p_components)+"\n")

P_TP3_A2m = TP_red3["A2m"][0]
# print(P_TP1_A1p)
TP3_A2m_vecs = r.list_nonzero_eigvecs(P_TP3_A2m)
P_file.write("A2m subspace:"+"\n") 
P_file.write(str(TP3_A2m_vecs)+"\n")

#(2,1,0) type single momenta

p7 = o.Pion([2,1,0],"+")
p8 = o.Pion([-2,-1,0],"+")

# Two Pions, distinguishable
tp4 = o.Two_Pion(p7,p8)
b_tp4 = o.generate_basis([tp4],O_h) 
P_file.write("TP4 Basis: (2,1,0)-type single momenta"+"\n")
for x in b_tp4:
    P_file.write(x.name + "\n")
o.print_all(b_tp4)
P_file.write("# "+ str(len(b_tp4))+"\n")
TP_Rep4 = r.rep_from_action(O_h,b_tp4,"TP_Rep4")

TP_red4 = r.find_irreps(TP_Rep4,O_h)

# A1p : 1
# T1m: 2
# T1p : 1
# T2m : 2
# T2p : 1
# Ep : 2
# A2p : 1

P_TP4_A1p = TP_red4["A1p"][0]
P_TP4_T1m = TP_red4["T1m"]

TP4_A1p_vecs = r.list_nonzero_eigvecs(P_TP4_A1p)
P_file.write("A1p subspace:"+"\n") 
P_file.write(str(TP4_A1p_vecs)+"\n")

TP4_vector_components = r.T1_identify_components(P_TP4_T1m,TP_Rep4)
P_file.write("(x,y,z)-like eigenvectors in TP4:"+"\n")
P_file.write(str(TP4_vector_components)+"\n")

P_file.write("Ep subspace:"+"\n")
P_TP4_Ep = TP_red4["Ep"]
TP4_Ep_components = r.E_identify_components(P_TP4_Ep,TP_Rep4)
P_file.write(str(TP4_Ep_components)+"\n")

P_file.write("T2m subspace:" + "\n")
P_TP4_T2m = TP_red4["T2m"]
TP4_T2m_components = r.T2_identify_components(P_TP4_T2m,TP_Rep4)
P_file.write(str(TP4_T2m_components)+"\n")

P_file.write("T2p subspace:" + "\n")
P_TP4_T2p = TP_red4["T2p"]
TP4_T2p_components = r.T2_identify_components(P_TP4_T2p,TP_Rep4)
P_file.write(str(TP4_T2p_components)+"\n")

P_TP4_A2p = TP_red4["A2p"][0]
# print(P_TP1_A1p)
TP4_A2p_vecs = r.list_nonzero_eigvecs(P_TP4_A2p)
P_file.write("A2p subspace:"+"\n") 
P_file.write(str(TP4_A2p_vecs)+"\n")

#(2,1,1) type single momenta

p9 = o.Pion([2,1,1],"+")
p10 = o.Pion([-2,-1,-1],"+")

# Two Pions, distinguishable
tp5 = o.Two_Pion(p9,p10)
b_tp5 = o.generate_basis([tp5],O_h) 
P_file.write("TP5 Basis: (2,1,1)-type single momenta:"+"\n")
for x in b_tp5:
    P_file.write(x.name + "\n")
o.print_all(b_tp5)
P_file.write("# " + str(len(b_tp5))+"\n")
TP_Rep5 = r.rep_from_action(O_h,b_tp5,"TP_Rep5")

TP_red5 = r.find_irreps(TP_Rep5,O_h)

# A1p : 1
# T1m : 2
# T1p : 1
# T2m : 1
# T2p : 2
# Em : 1
# Ep : 1
# A2m : 1

P_TP5_A1p = TP_red5["A1p"][0]
P_TP5_T1m = TP_red5["T1m"]

TP5_A1p_vecs = r.list_nonzero_eigvecs(P_TP5_A1p)
P_file.write("A1p subspace:"+"\n") 
P_file.write(str(TP5_A1p_vecs)+"\n")

TP5_vector_components = r.T1_identify_components(P_TP5_T1m,TP_Rep5)
P_file.write("(x,y,z)-like eigenvectors in TP5:"+"\n")
P_file.write(str(TP5_vector_components)+"\n")

P_file.write("Ep subspace:"+"\n")
P_TP5_Ep = TP_red5["Ep"]
TP5_Ep_components = r.E_identify_components(P_TP5_Ep,TP_Rep5)
P_file.write(str(TP5_Ep_components)+"\n")

P_file.write("Em subspace:"+"\n")
P_TP5_Em = TP_red5["Em"]
P_file.write(str(P_TP5_Em)+"\n")
TP5_Em_components = r.E_identify_components(P_TP5_Em,TP_Rep5)
P_file.write(str(TP5_Em_components)+"\n")

P_file.write("T2m subspace:" + "\n")
P_TP5_T2m = TP_red5["T2m"]
TP5_T2m_components = r.T2_identify_components(P_TP5_T2m,TP_Rep5)
P_file.write(str(TP5_T2m_components)+"\n")

P_file.write("T2p subspace:" + "\n")
P_TP5_T2p = TP_red5["T2p"]
TP5_T2p_components = r.T2_identify_components(P_TP5_T2p,TP_Rep5)
P_file.write(str(TP5_T2p_components)+"\n")

P_TP5_A2m = TP_red5["A2m"][0]
# print(P_TP1_A1p)
TP5_A2m_vecs = r.list_nonzero_eigvecs(P_TP5_A2m)
P_file.write("A2m subspace:"+"\n") 
P_file.write(str(TP5_A2m_vecs)+"\n")


################################
P_file.close()
#test 

# sca = o.ScalarParticle([0,0,0])
# b_sca = o.generate_basis(sca,O_h)
# sca_rep = r.rep_from_action(O_h,b_sca,"scalar_rep")
# o.print_all(sca_rep.basis)
# print(sca_rep.hom["Inv"])
# sca_rep.check_if_homomorphism()
# r.study_irreps(sca_rep,O_h,"../results/ScalarParticle_irreps.txt")

psca = o.PseudoScalarParticle([0,0,0])
b_psca = o.generate_basis(psca,O_h)
psca_rep = r.rep_from_action(O_h,b_psca,"pseudoscalar_rep")
o.print_all(psca_rep.basis)
print(psca_rep.hom["Inv"])
# psca_rep.check_if_homomorphism()
# r.study_irreps(psca_rep,O_h,"../results/PseudoScalarParticle_irreps.txt")

vec = o.VectorParticle([0,0,0])
b_vec = o.generate_basis(vec,O_h)
vec_rep = r.rep_from_action(O_h,b_vec,"vector_rep")
o.print_all(vec_rep.basis)
vec_rep.check_if_homomorphism()
print(vec_rep.hom["Inv"])
r.study_irreps(vec_rep,O_h,"../results/VectorParticle_irreps.txt" )

pvec = o.PseudoVectorParticle([0,0,0])
b_pvec = o.generate_basis(pvec,O_h)
pvec_rep = r.rep_from_action(O_h,b_pvec,"pseudovector_rep")
o.print_all(pvec_rep.basis)
print(pvec_rep.hom["Inv"])
pvec_rep.check_if_homomorphism()
r.study_irreps(pvec_rep,O_h,"../results/PseudoVectorParticle_irreps.txt")






