import priorityqueue

randos = [68, 29, 94, 25, 69, 62, 43, 49, 1, 90, 
			23, 69, 25, 44, 41, 4, 75, 84, 41, 39]


#create instance of priority queue
pq = priorityqueue.PriorityQueue(30)

# insert 20 integers into priority queue
for r in randos:
    pq.insert(r)

print "My priorty queue is: \n",pq.a
print "\nnow lets delMin():"

# delMin() from the priority queue 20 times and print
for i in range(20):
    print "min = ",pq.delMin()