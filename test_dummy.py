import numpy as np
class Vector:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def copy(self):
        x = self.x
        y = self.y
        z = self.z
        new_v = Vector(x,y,z)
        return new_v
    def is_equal_to(self,v):
        return np.allclose([self.x,self.y,self.z],[v.x,v.y,v.z], rtol = 1e-6, atol = 1e-10)

test = Vector(1,2,3)
test2 = test.copy()
test2.x = 1.000001
print(test2.x)
print(test.x)
print(test.is_equal_to(test2))