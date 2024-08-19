import numpy as np
import groups as g
import representations as r
import objects as o

np.set_printoptions(precision = 8, suppress = True)

gen_actions = ["Rot0","Rot1","Rot2",]
v = o.Vector(["a","b","c"])
ac = o.find_closure(gen_actions,v)
mt = o.mult_table_from_actions(ac)
O = o.generate_group(gen_actions,v)

print("# elements of O: ", len(O.elements))