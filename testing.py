import numpy as np
import groups as g
import objects as o
import representations as r
import numbers
def transform_by_explicit_action(A,Rep,lin_comb):               # takes String A, e.g. "Rot0", Representation object, vector according to the basis Representation.basis
                                                                # returns lin_comb after acting on each component with A