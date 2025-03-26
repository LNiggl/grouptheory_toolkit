from . import objects as o
from . import representations as r
from . import tools as t
from .definitions import gen_actions_O, gen_actions_O_h


def generate_O():
    print("Generating O.")
    v = o.Vector(["a","b","c"])
    O = o.generate_group(gen_actions_O,v)

    # irreps of O

    a = [o.Scalar(1)]
    A1 = r.rep_from_action(O,a,"A1")
    A1.check_if_homomorphism()

    v1 = o.Vector([1,0,0])
    b = o.generate_basis(v1,O)
    T1 = r.rep_from_action(O,b,"T1")

    T1_x_T1 =  r.product_rep(T1,T1)           #no irrep
    T1_x_T1.check_if_homomorphism()


    T2 = T1_x_T1.copy("T2")
    r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2)
    T2.check_if_homomorphism()


    E = T1_x_T1.copy("E")
    r.apply_projectors([r.traceless_projector,r.diagonal_projector],E)
    E.check_if_homomorphism()
    E.round_chars()


    A2 = r.product_rep(T1_x_T1,T1)
    A2.name = "A2"
    r.project_out_irreps(A2, [A1,T1,T2,E,])
    A2.check_if_homomorphism()
    # print("A2 reducible: ", A2.is_reducible())                # irreducible
    A2.round_chars()


    list_irreps = [A1,T1,A2,T2,E]
    O.set_char_table(list_irreps)

    t.save_group(O,"./base/O/")
    return O

def generate_O_h():
    print("Generating O_h.")
    v = o.Vector(["a","b","c"])
    O_h = o.generate_group(gen_actions_O_h,v)

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
    print("shape " , (A2m.hom["I"]).shape)
    A2m.name = "A2m"
    r.project_out_irreps(A2m, [A1p,A1m,T1m,T1p,T2p,T2m,Em,Ep])
    A2m.check_if_homomorphism()
    A2m.round_chars()

    A2p = r.product_rep(A2m,A1m)
    A2p.name = "A2p"
    A2p.check_if_homomorphism()

    list_irreps = [A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep]
    O_h.set_char_table(list_irreps)

    t.save_group(O_h,"./base/O_h/")
    return O_h

def generate_DC_O():
    print("Generating DC_O.")
    DiracS = o.DiracSpinor(1,1,1,1)
    DC_O = o.generate_group(gen_actions_O,DiracS)
    a = [o.Scalar(1)]
    A1 = r.rep_from_action(DC_O,a,"A1")
    A1.check_if_homomorphism()

    s1 = o.DiracSpinor(1,0,1,0)
    s2 = o.DiracSpinor(0,1,0,1)
    G1 = r.rep_from_action(DC_O,[s1,s2],"G1")
    G1.check_if_homomorphism()

    T1 = G1.copy("T1")
    for A in T1.elements:
        T1.hom[A] = r.R(T1.hom[A])
    T1.update()   
    T1.check_if_homomorphism()
    # print(T1.characters)

    T1_x_T1 =  r.product_rep(T1,T1)           #no irrep
    T1_x_T1.check_if_homomorphism()

    T2 = T1_x_T1.copy("T2")
    r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2)
    T2.check_if_homomorphism()

    E = T1_x_T1.copy("E")
    r.apply_projectors([r.traceless_projector,r.diagonal_projector],E)
    E.check_if_homomorphism()

    A2 = r.product_rep(T1,r.product_rep(T1,T1))
    A2.name = "A2"
    r.project_out_irreps(A2, [A1,T1,T2,E])
    A2.check_if_homomorphism()

    G2 = r.product_rep(A2,G1)
    G2.name = "G2"
    G2.check_if_homomorphism()

    H = r.product_rep(G1,E)
    H.name = "H"
    H.check_if_homomorphism()
    # print("H reducible: ", H.is_reducible())

    DC_O.set_char_table([A1,T1,A2,T2,E,G1,G2,H])
    t.save_group(DC_O,"./base/DC_O/") 
    print("#elements DC_O", len(DC_O.elements))
    return DC_O

def generate_DC_O_h():
    print("Generating DC_O_h.")
    DiracS = o.DiracSpinor(1,2,3,4)
    DC_Oh = o.generate_group(gen_actions_O_h,DiracS)
    a = [o.Scalar(1)]
    A1p = r.rep_from_action(DC_Oh,a,"A1p")
    A1p.check_if_homomorphism()

    A1m = r.rep_from_action(DC_Oh,[o.PseudoScalar(1)],"A1m")
    A1m.check_if_homomorphism()
    # print("A1m candidate reducible: ",A1m.is_reducible())                             # is irreducible
    # print("characters: ",A1m.characters)                               

    Dirac1 = o.DiracSpinor(1,0,1,0)
    Dirac2 = o.DiracSpinor(0,1,0,1)
    G1p = r.rep_from_action(DC_Oh,[Dirac1,Dirac2],"G1p")
    # print("G1p candidate reducible: ",G1p.is_reducible())                             # is irreducible
    # print("characters: ",G1p.characters)                                                                           
    G1p.check_if_homomorphism()

    G1m = r.product_rep(G1p,A1m)
    G1m.name = "G1m"
    G1m.check_if_homomorphism()
    # print("G1m candidate reducible: ",G1m.is_reducible())                             # is irreducible
    # print("characters: ",G1m.characters)            

    T1p = G1p.copy("T1p")                                                               
    for A in T1p.elements:
        T1p.hom[A] = r.R(T1p.hom[A])
    T1p.update()   
    T1p.check_if_homomorphism()
    # print("T1p candidate reducible: ",T1p.is_reducible())                             # irreducible
    # print("characters: ",T1p.characters)  

    T1m = r.product_rep(T1p,A1m)
    T1m.name = "T1m"
    T1m.check_if_homomorphism()
    # print("T1m candidate reducible: ",T1m.is_reducible())                             # irreducible
    # print("characters: ",T1m.characters)

    T1m_x_T1m =  r.product_rep(T1m,T1m)           #no irrep
    T1m_x_T1m.check_if_homomorphism()

    T2p = T1m_x_T1m.copy("T2p")
    r.apply_projectors([r.symmetric_projector,r.invert_projector(r.diagonal_projector)],T2p)
    T2p.check_if_homomorphism()
    # print("T2p candidate reducible: ",T2p.is_reducible())                           # is irreducible
    # print("characters: ",T2p.characters)                               
    

    T2m = r.product_rep(T2p,A1m)
    T2m.name = "T2m"
    T2m.check_if_homomorphism()
    # print("T2m candidate reducible: ",T2m.is_reducible())                           # is irreducible
    # print("characters: ",T2m.characters)                               

    Ep = T1m_x_T1m.copy("Ep")
    r.apply_projectors([r.traceless_projector,r.diagonal_projector],Ep)
    Ep.check_if_homomorphism()
    # print("Ep candidate reducible: ",Ep.is_reducible())                           # is irreducible
    # print("characters: ",Ep.characters)                                

    Em = r.product_rep(Ep,A1m)
    Em.name = "Em"
    Em.check_if_homomorphism()
    # print("Em candidate reducible: ",Em.is_reducible())                           # is irreducible
    # print("characters: ",Em.characters)                               

    A2m = r.product_rep(T1m,r.product_rep(T1m,T1m))
    A2m.name = "A2m"
    r.project_out_irreps(A2m, [A1p,A1m,T1p,T1m,T2p,T2m,Ep,Em])
    A2m.check_if_homomorphism()
    # print("A2m candidate reducible: ",A2m.is_reducible())                           # is irreducible
    # print("characters: ",A2m.characters)                               

    A2p = r.product_rep(A2m,A1m)
    A2p.name = "A2p"
    A2p.check_if_homomorphism()
    # print("A2p candidate reducible: ",A2p.is_reducible())                           # is irreducible
    # print("characters: ",A2p.characters)                               

    G2p = r.product_rep(A2p,G1p)
    G2p.name = "G2p"
    G2p.check_if_homomorphism()
    # print("A2p candidate reducible: ",G2p.is_reducible())                           # is irreducible
    # print("characters: ",G2p.characters)                               

    G2m = r.product_rep(G2p,A1m)
    G2m.name = "G2m"
    G2m.check_if_homomorphism()
    # print("G2m candidate reducible: ",G2m.is_reducible())                           # is irreducible
    # print("characters: ",G2m.characters)                               

    Hp = r.product_rep(G1p,Ep)
    Hp.name = "Hp"
    Hp.check_if_homomorphism()
    # print("Hp candidate reducible: ",Hp.is_reducible())                           # is irreducible
    # print("characters: ",Hp.characters)                               

    Hm = r.product_rep(Hp,A1m)
    Hm.name = "Hm"
    Hm.check_if_homomorphism()
    irreps_DC_O_h = [A1m,A1p,T1m,T1p,A2m,A2p,T2m,T2p,Em,Ep,G1m,G1p,G2m,G2p,Hm,Hp]
    # print("# irreps: ", len(irreps_DC_Oh), ". # classes: ", len(DC_Oh.classes))
    DC_Oh.set_char_table(irreps_DC_O_h)
    t.save_group(DC_Oh,"./base/DC_O_h/") 
    return DC_Oh

generate_O()

import os
if not os.path.isdir("./base/O_h/"):
    generate_O_h()
if not os.path.isdir("./base/O/"):
    generate_O()
if not os.path.isdir("./base/DC_O_h/"):
    generate_DC_O_h()
if not os.path.isdir("./base/DC_O/"):
    generate_DC_O()

O_h = t.load_group("./base/O_h/")
O = t.load_group("./base/O/")
DC_O_h = t.load_group("./base/DC_O_h/")
DC_O = t.load_group("./base/DC_O/")
