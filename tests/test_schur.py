import numpy as np
import representations as r
import compare as c
from scipy.linalg import schur,eigvals,solve
P = np.load("../results/projectors/Rep_T2m_t1.npy")
np.save("./projector",P)
T,Z = schur(P)
ev = eigvals(T)
# t,z = np.linalg.eig(P)
# for i in range(len(t)):
#     if abs(t[i])>0.1:
#         print(z[i])
print(ev)
ev_nonzero = []
for i in range(len(ev)):
    if abs(ev[i])>0.1:
        ev_nonzero.append(ev[i])
# print(ev_nonzero)
evecs = []
for i in range(len(ev_nonzero)):
    ev_i = [ev_nonzero[i] for j in range(len(P))]
    print(np.diag(ev_i)-P)
    evec = solve(np.diag(ev_i)-P,np.zeros(len(P)))
    evecs.append(evec)
print(evecs)
f = open("../tests/test_list_nonzero_eigvals2.txt","w")
d = r.list_nonzero_eigvecs(P)
f.write(str(d))
f.close()

c.compare_strings("../tests/test_list_nonzero_eigvals.txt","../tests/test_list_nonzero_eigvals2.txt")




