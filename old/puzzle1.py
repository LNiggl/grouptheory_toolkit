import numpy as np
#create group by studying its defining action on basis
basis = [[1,0,0],[0,1,0],[0,0,1]]
def interchange(x,a,b):
    xp=[xi for xi in x]
    xp[a]=x[b]
    xp[b]=x[a]
    return xp
def reflection(x,a):
    xp=[xi for xi in x]
    xp[a]=-x[a]
    return xp
def rotateZ90(x):
    return [-x[1],x[0],x[2]]
def rotateX90(x):
    return [x[0],-x[2],x[1]]
def find_closure(actions,basis):        #give actions on a suited basis for a vectorspace        
    elements = [basis]
    sequence = { str(basis) : [] }
    n = 0
    while n != len(elements):
        n = len(elements)
        for a in actions:
            for e in elements:
                ep = a(e)
                if ep not in elements:
                    elements.append(ep)
                    # print([a])
                    sequence[str(ep)] = sequence[str(e)] + [a]
    return elements, sequence

# First only study group of rotations, i.e., cubic group O
elements = find_closure([
    #lambda e: [interchange(b,0,1) for b in e],
    #lambda e: [interchange(b,0,2) for b in e],
    #lambda e: [interchange(b,1,2) for b in e],
    #lambda e: [reflection(b,0) for b in e],
    #lambda e: [reflection(b,1) for b in e],
    #lambda e: [reflection(b,2) for b in e],
    lambda e: [rotateZ90(b) for b in e],
    lambda e: [rotateX90(b) for b in e]
],basis)
# print(elements[1])
# print(type(elements[1]['[[0, -1, 0], [0, 0, -1], [1, 0, 0]]']))
# key = [k for k,v in elements[1].items() if v  [<function <lambda> at 0x0000018F0C9215E0>, <function <lambda> at 0x0000018F0C9215E0>]]
# print(key)
# print(elements[1].index("[<function <lambda> at 0x0000018F0C9215E0>, <function <lambda> at 0x0000018F0C9215E0>]"))
def apply(el, e0, idx):
    for acc in el[1][str(el[0][idx])]:          #applies actions of entire sequence to create group element el[0][idx]
        e0 = acc(e0)
    return e0                                   #return vector where this element has acted on basis
print(apply(elements,[[1,0,0]],1))
def multiplication_table(el, ba):
    n=len(el[0])
    M=[[-1 for i in range(n)] for j in range(n)]
    for i in range(n):                          #loop through all elements i of group
        ei = apply(el, ba, i)                   #apply action of ith element to basis -> ei
        for j in range(n):                      #apply action of jth element to ei
            eji = apply(el, ei, j)
            M[i][j] = el[0].index(eji)          #read index -> in M[i][j] stands result of ej*ei
    return M

mt = multiplication_table(elements, basis)
def find_identity(el, ba):              #identity is element that leaves basis invariant
    return el[0].index(ba)
ident = find_identity(elements, basis)      
assert ident == 0 # for convenience always have first element be the identity
def inverse_list(mt, e):
    return [x.index(e) for x in mt]     #list of elements y such that  yx = xy = e
inv = inverse_list(mt, ident)
# find conjugation classes / invariant subspaces
def find_conjugation_classes(mt, inv, ident):
    n=len(inv)
    tostudy=list(range(n))
    classes = []
    while len(tostudy) > 0:
        x = tostudy[0]
        conjugate_to = []
        for i in range(n):
            conjugate_to.append(mt[mt[i][x]][inv[i]])
        conjugate_to = set(conjugate_to)
        for y in conjugate_to:
            tostudy.remove(y)
        classes.append(conjugate_to)
    return classes
classes = find_conjugation_classes(mt, inv, ident)
print("classes: ",classes)
# find a matrix representation of individual elements
def matrix_vector(l):
    # transform a lattice vector
    return np.matrix(elements[0][l]).T

def matrix_trivial(l):
    return np.matrix([[1]])

def matrix_determinant(l):
    return np.matrix([[np.linalg.det(matrix_vector(l))]], dtype=np.int32)

def matrix_regular(l):              #creates regular rep matrix for l-th item in group
    M = np.zeros(shape=(len(mt),len(mt)), dtype=np.int32)
    for i in range(len(mt)):
        M[i,mt[l][i]] = 1           #sets a 1 at index j in row i where g_l*g_j = g_i
    return np.matrix(M)

def test_representation(matrix):
    print(matrix.__name__)
    for i in range(len(mt)):
        for j in range(len(mt)):
            eps = np.linalg.norm(matrix(j) * matrix(i) - matrix(mt[i][j]))
            assert eps < 1e-13

############# methods for action on pion and representation
def pion(p):                        #returns "important symmetry characteristics" of pion; here: momentum                      
    return p
def P(p):
    p = np.array(p)
    return -1*p
    # return [-1*pi for pi in p]
def pion_action(g,group,inv,basis):          #NEEDS global mt and inv; insert group element and basis for momenta p of pion
    g_inv = inv[group.index(g)]
    return [pion(np.matmul(matrix_vector(g_inv),p)) for p in basis]

def build_rep_matrix(g,action,basis):                  # takes group element g, function for the action and appropriate basis; returns matrix of action of g on basis
    # for b in basis:
    #     b = np.matrix(b)
    t_basis = action(g,basis)
    M = np.zeros(shape=(len(basis),len(basis)), dtype=np.int32)
    for i in range(len(basis)):
        for j in range(len(basis)):
            # n = np.matmul(t_basis[i], basis[j])
            if (t_basis[i] == basis[j]).all() :         #! expand to check for proportionality rather than only equality
                M[i][j] = 1
    return M

def matrix_rep(group,action,basis):
    # list_matrices = {ident : np.identity(len(basis))}
    list_matrices = {}
    for i in range(len(group)): 
        M = build_rep_matrix(group[i],action,basis)
        list_matrices[i] = M
    return list_matrices

def test_matrix_rep(mt, rep):                           #takes mult table and checks if rep is homomorphism from G -> GL(n)
    for i in range(len(mt)):
        for j in range(len(mt)):
            print("------------------")
            print(rep[j] * rep[i])
            print(rep[mt[i][j]])
            eps = np.linalg.norm(rep[j] * rep[i] - rep[mt[i][j]])
            # assert eps < 1e-13
            print("difference: ",eps)

def pion_product_rep(group,action,basis):
    b1 = basis
    b2 = [P(p) for p in p_basis]
    M1 = list(matrix_rep(group,action,b1).values())    
    M2 = list(matrix_rep(group,action,b2).values())
    return np.kron(M1,M2)


    

print("test pion action")

p_basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
#print(pion_action(elements[0][1],p_basis))
# print(build_rep_matrix(elements[0][1],pion_action,p_basis))
pion_rep = matrix_rep(elements[0],pion_action,p_basis)
# print(pion_rep[6])
# test_matrix_rep(mt,pion_rep)
# print(P(p_basis))
product_rep = pion_product_rep(elements[0],pion_action,p_basis)
print(product_rep[0],product_rep[5])
test_matrix_rep(mt,product_rep)


O_h = find_closure([
    lambda e: [interchange(b,0,1) for b in e],
    lambda e: [interchange(b,0,2) for b in e],
    lambda e: [interchange(b,1,2) for b in e],
    lambda e: [reflection(b,0) for b in e],
    lambda e: [reflection(b,1) for b in e],
    lambda e: [reflection(b,2) for b in e]
], basis)

mt = multiplication_table(elements, basis)
ident = find_identity(elements, basis)      
assert ident == 0 # for convenience always have first element be the identity
inv = inverse_list(mt, ident)
pion_O_h = matrix_rep(O_h[0],pion_action,p_basis)
product_p_O_h = pion_product_rep([0],pion_action,p_basis)











######################    

# p_basis = [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]] 
def pion_rotateZ90(x):
    return [x[1],x[3],x[2],x[4],x[0],x[5]]
def pion_rotateX90(x):
    return [x[0],x[2],x[4],x[3],x[5],x[1]]
# basis = p_basis
# p_elements = find_closure([
#     lambda e: [pion_rotateZ90(b) for b in e],
#     lambda e: [pion_rotateX90(b) for b in e]
#     ],basis)
# print(len(p_elements)) 
# for e in range(len(p_elements)):
#     print(p_elements[0][e])

# print(pion_rotateZ90(pion_rotateZ90(pion_rotateZ90(pion_rotateZ90(p_basis[0])))))