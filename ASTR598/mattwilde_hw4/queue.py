class Node(object):
	""" 
	An element in a linked list stack with a data attribute
	and a pointer attribute. Initializes to None.
	"""
	def __init__(self, data=None, next=None):
		self.data = data
		self.next = next

class Queue(object):
	""" A class representing a Queue.
	
	First in First Out.

	Composed of a Linked List.
	"""
	def __init__(self, head=None, tail=None, N=0):
		self.tail = tail
		self.head = head
		self.N = N

	def put(self,i):
		""" Add a node to the tail of the Queue.

			input = int i

			returns None
		"""
		# create an instance of a Node with data i
		n = Node(data=i, next=None)
		
		# if empty
		if self.head == None:
			self.head = n
			self.tail = n
		# if length 1
		elif self.tail == self.head:
			self.tail = n # add to the tail
			self.head.next = n # set the the nodes pointer to the 
		# generally
		else:
			self.tail.next = n # tails pointer to n. should be None?
			self.tail = n # add to the tail

		#increment size
		self.N += 1
		return None

	def get(self):
		# empty
		if self.head == None:
			print "cannot get. queue empty"
			return None
		
		else:
			tmp = self.head.data
			self.head = self.head.next
			self.N -= 1
			return tmp

	# like get but does not remove head
	def front(self):
		# empty
		if self.head == None:
			print "cannot front, empty queue"
			return None
		else:
			return self.head.data

	def size(self):
		return self.N

	def isEmpty(self):
		if self.N == 0:
			return True
		else:
			return False


