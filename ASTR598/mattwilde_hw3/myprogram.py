import numpy as np
import stack

randos = np.random.randint(10, size=10)

s = stack.Stack()
for r in randos:
	s.push(r)
	print r


for i in range(11):
	print s.pop()
