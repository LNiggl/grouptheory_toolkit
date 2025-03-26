import numpy as np

from base import groups as g
from base import objects as o
from base import representations as r
from base import tools as t
from base.definitions import A,gen_actions_O,gen_actions_O_h                                  # generating set of matrices for G1 irrep of the double cover

##############################################################
#####   Demos of how to generate the double cover groups; 
#       for the irreps, see the functions in __init__.py  ####

#####   generation of DC_O from the G1 irrep matrices   ######

gen_matrices = {}
gen_matrices["A0"] = A[0]
gen_matrices["A1"] = A[1]
gen_matrices["A2"] = A[2]

dcover = r.find_closure_from_matrices(gen_matrices)
mtable = r.mult_table_from_matrices(dcover)
DC_O_mat = g.Group(list(dcover.keys()),mtable)
print("# DC_O from G1 matrices:", len(DC_O_mat.elements))

# assignment of G1 matrices as Representation object
G1_mat = r.Representation(DC_O_mat,dcover,"G1_mat")
G1_mat.direction_action = "left"
# assignment of basis by hand
s1 = o.WeylSpinor(1,0)
s2 = o.WeylSpinor(0,1)
G1_mat.basis = [s1,s2]
G1_mat.check_if_homomorphism()


#####   development of DC_O and G1 from action on Weyl-spinor object    #####   

WS = o.WeylSpinor(1/np.sqrt(2),1/np.sqrt(2))
DC_O_Weyl = o.generate_group(gen_actions_O,WS)

s1 = o.WeylSpinor(1,0)
s2 = o.WeylSpinor(0,1)
G1_Weyl = r.rep_from_action(DC_O_Weyl,[s1,s2],"G1_Weyl")
G1_Weyl.check_if_homomorphism()


#####   compare matrices of G1 of both approaches   #####
R = r.group_hom_via_reps(G1_Weyl,G1_mat)                    # map of group element names 
for A in O.elements:
    assert np.allclose(G1_Weyl.hom[A]-G1_mat.hom[R[A][0]],np.zeros(G1_Weyl.hom[A].shape))

#####   DC O from Dirac Spinor  #####
DiracG = o.DiracSpinor(1,2,3,4)
DC_O_Dirac = o.generate_group(gen_actions_O,DiracG)
print("# DC_O_Dirac: ", len(DC_O_Dirac.elements))

Dirac1 = o.DiracSpinor(1,0,1,0)
Dirac2 = o.DiracSpinor(0,1,0,1)
Dirac3 = o.DiracSpinor(1,0,-1,0)
Dirac4 = o.DiracSpinor(0,1,0,-1)
basis_D = [Dirac1,Dirac2,Dirac3,Dirac4]
G1_Dirac = r.rep_from_action(DC_O_Dirac,basis_D,"G1_Dirac")
G1_Dirac.check_if_homomorphism()
print("G1_Dirac irreducible: ", G1_Dirac.is_reducible())

#####   DC O_h from Dirac Spinor object   #####
DiracG = o.DiracSpinor(1,2,3,4)
DC_O_h_Dirac = o.generate_group(gen_actions_O_h,DiracG)
print("# DC_O_Dirac: ", len(DC_O_h_Dirac.elements))

Dirac1 = o.DiracSpinor(1,0,1,0)
Dirac2 = o.DiracSpinor(0,1,0,1)
Dirac3 = o.DiracSpinor(1,0,-1,0)
Dirac4 = o.DiracSpinor(0,1,0,-1)
basis_D = [Dirac1,Dirac2,Dirac3,Dirac4]
G1p_Dirac = r.rep_from_action(DC_O_h_Dirac,basis_D,"G1m_Dirac")
G1p_Dirac.check_if_homomorphism()
print("G1p_Dirac irreducible: ", G1p_Dirac.is_reducible())


#############################################
#####   test study - DiracField object  #####

from base import DC_O,DC_O_h

DF01 = o.DiracField([0,0,0],struc = [1,0,1,0])
DF02 = o.DiracField([0,0,0],struc = [0,1,0,1])
DF03 = o.DiracField([0,0,0],struc = [1,0,-1,0])
DF04 = o.DiracField([0,0,0],struc = [0,1,0,-1])
b3 = [DF01,DF02,DF03,DF04]

DiracRep = r.rep_from_action(DC_O_h,b3,"DiracRep")
DiracRep.check_if_homomorphism()

#####   test study - Pion(100)-DiracField(-100) rep #####

v = o.Vector([1,0,0])
all_vs = o.orbit(v,DC_O_h)
o.print_all(all_vs)

basis_DS = []
basis_pi = []
for vec in all_vs:
   basis_DS.append(o.DiracField(vec,struc = [1,0,1,0],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [1,0,-1,0],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [0,1,0,1],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [0,1,0,-1],modified_momentum_trafo=True))
   basis_pi.append(o.PseudoScalarField(vec,modified_momentum_trafo=True))
# o.print_all(basis_DS)

basis_DPi100 = []
for DS in basis_DS:
   for pi in basis_pi:
    if DS.momentum.is_negative_of(pi.momentum):
      basis_DPi100.append(o.TensorProduct(DS,pi))
# o.print_all(basis_DPi100)
print("# basis DiracPi: ", len(basis_DPi100))

DPi100_rep = r.rep_from_action(DC_O_h,basis_DPi100,"DPi100")
DPi100_rep.check_if_homomorphism()

filepath = "D:/Master/Masterarbeit/tests/summary_irreps/DPi100.txt"
spacesDPi100 = r.study_irreps(DPi100_rep,DC_O_h,filepath)

DPi100_spaces_ordered = t.reorder_subspaces(spacesDPi100)
trafos = t.test_trafos_matmul(DPi100_spaces_ordered,["Rot0","Rot1","Rot2","Inv"],DPi100_rep)

#####   test study - Pion(110)-DiracField(-110) rep #####

v = o.Vector([1,1,0])
all_vs = o.orbit(v,DC_O_h)
o.print_all(all_vs)

basis_DS = []
basis_pi = []
for vec in all_vs:
   basis_DS.append(o.DiracField(vec,struc = [1,0,1,0],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [1,0,-1,0],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [0,1,0,1],modified_momentum_trafo=True))
   basis_DS.append(o.DiracField(vec,struc = [0,1,0,-1],modified_momentum_trafo=True))
   basis_pi.append(o.PseudoScalarField(vec,modified_momentum_trafo=True))
# o.print_all(basis_DS)

basis_DPi110 = []
for DS in basis_DS:
   for pi in basis_pi:
    if DS.momentum.is_negative_of(pi.momentum):
      basis_DPi110.append(o.TensorProduct(DS,pi))
# o.print_all(basis_DPi110)
print("# basis DiracPi: ", len(basis_DPi110))

DPi110_rep = r.rep_from_action(DC_O_h,basis_DPi110,"DPi110")
DPi110_rep.check_if_homomorphism()

filepath = "D:/Master/Masterarbeit/tests/summary_irreps/DPi110.txt"
spacesDPi110 = r.study_irreps(DPi110_rep,DC_O_h,filepath)

DPi110_spaces_ordered = t.reorder_subspaces(spacesDPi110)
trafos = t.test_trafos_matmul(DPi110_spaces_ordered,["Rot0","Rot1","Rot2","Inv"],DPi110_rep,skip_H = True)

# Note: For Hp and Hm, r.H_identify_components() only identifies the basis vector via their eigenvalues and does not establish a consistent connection via rotations between them

