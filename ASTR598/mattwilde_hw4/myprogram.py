import queue

# create a list of 10 'random' numbers
randos = [1,3,26,9,19,22,17,15,8,10]

# create an instance of a queue
q = queue.Queue()

# put the random integers on the queue
for r in randos:
	q.put(r)

for i in range(10):
	print q.get()

print q.isEmpty()