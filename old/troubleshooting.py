import numpy as np
import groups as g
import representations as r
import objects as o
import compare 

np.set_printoptions(precision = 8, suppress = True)

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
# A2p.round_chars()

O_h.set_char_table([A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep])
# chars = open("../results/troubleshoot_chars2.txt","w")
print("Character table of O_h:")
for irrep in O_h.char_table.keys():
    print(irrep, ": " , O_h.char_table[irrep])
    # chars.write(irrep + ": " + str(O_h.char_table[irrep]))
#(2,1,0) - type 2Pion    
p7 = o.PseudoScalarField([2,1,0])
p8 = o.PseudoScalarField([-2,-1,0])
pi4 = o.TensorProduct(p7,p8)

b4 = o.generate_basis(pi4,O_h)
pi210 = r.rep_from_action(O_h,b4,"pi210")

#producing the data to compare (run one after the other)
# r.study_irreps(pi210,O_h,"../results/troubleshoot1.txt")
r.study_irreps(pi210,O_h,"../results/troubleshoot2.txt")

#comparing the results
compare.nice_string("../results/troubleshoot1.txt")
compare.nice_string("../results/troubleshoot2.txt")
compare.compare_strings("../results/troubleshoot1_nice.txt","../results/troubleshoot2_nice.txt") 
#-> difference as soon as lin.comb. get "multi-valued". Here: T1m

#comparing the projectors
# for irrep in O_h.char_table.keys():
#     compare.nice_string("../results/troubleshoot1_"+ irrep + ".txt")
#     compare.nice_string("../results/troubleshoot2_" + irrep + ".txt")
# for irrep in O_h.char_table.keys():
#     compare.compare_strings("../results/troubleshoot1_"+ irrep + "_nice.txt","../results/troubleshoot2_" + irrep + "_nice.txt")    
# compare.nice_string("../results/troubleshoot1_T2p.txt")
#-> no differences

#comparing all matrices after projection
# for irrep in O_h.char_table.keys():
#     compare.nice_string("../results/troubleshoot1_"+ irrep + "_hom.txt")
#     compare.nice_string("../results/troubleshoot2_" + irrep + "_hom.txt")
#     compare.compare_strings("../results/troubleshoot1_"+ irrep + "_hom_nice.txt","../results/troubleshoot2_" + irrep + "_hom_nice.txt") 
#-> no differences

#comparing eigenvectors of P for T1p/m
# for irrep in ["T1m","T1p"]:
#     compare.nice_string("../results/troubleshoot1_"+ irrep + "_vecs.txt")
#     compare.nice_string("../results/troubleshoot2_" + irrep + "_vecs.txt")
#     compare.compare_strings("../results/troubleshoot1_" + irrep + "_vecs_nice.txt" , "../results/troubleshoot2_" + irrep + "_vecs_nice.txt")
# -> inconsistencies



# testvectors = [np.matrix([[ 0.       -0.j        ],
#         [-0.       -0.j        ],
#         [-0.       -0.j        ],
#         [ 0.       -0.j        ],
#         [ 0.       -0.j        ],
#         [ 0.3534478-0.00864032j],
#         [-0.3534478+0.00864032j],
#         [-0.       +0.j        ],
#         [-0.       -0.j        ],
#         [ 0.       -0.j        ],
#         [-0.3534478+0.00864032j],
#         [ 0.3534478-0.00864032j],
#         [-0.3534478+0.00864032j],
#         [ 0.3534478-0.00864032j],
#         [-0.       -0.j        ],
#         [ 0.       +0.j        ],
#         [-0.       +0.j        ],
#         [ 0.3534478-0.00864032j],
#         [-0.3534478+0.00864032j],
#         [ 0.       -0.j        ],
#         [-0.       +0.j        ],
#         [ 0.       +0.j        ],
#         [-0.       -0.j        ],
#         [ 0.       -0.j        ]]), np.matrix([[ 0.29178672+0.j        ],
#         [ 0.        +0.j        ],
#         [ 0.        +0.j        ],
#         [-0.29178672+0.j        ],
#         [ 0.        +0.19965097j],
#         [ 0.        +0.j        ],
#         [ 0.        +0.j        ],
#         [ 0.        -0.19965097j],
#         [ 0.        +0.19965097j],
#         [ 0.        +0.19965097j],
#         [ 0.29178672+0.j        ],
#         [ 0.29178672+0.j        ],
#         [-0.29178672+0.j        ],
#         [-0.29178672+0.j        ],
#         [ 0.        -0.19965097j],
#         [ 0.        -0.19965097j],
#         [ 0.        +0.19965097j],
#         [ 0.        +0.j        ],
#         [ 0.        +0.j        ],
#         [ 0.        -0.19965097j],
#         [ 0.29178672+0.j        ],
#         [ 0.        +0.j        ],
#         [ 0.        +0.j        ],
#         [-0.29178672+0.j        ]])]
# for i in [0,1]:
#     print(r.rotate_to_real_valued(testvectors[i]))


               
