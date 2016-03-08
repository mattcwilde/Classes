import bst
import numpy as np

randos = np.random.randint(0,100,100)

b = bst.BinarySearchTree()

for i in randos:
    b.insert(i)
    
    
b.traverse()
