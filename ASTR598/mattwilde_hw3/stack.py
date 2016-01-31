class Node:
    """ 
    An element in a linked list stack with a data attribute
    and a pointer attribute. Initializes to None.
    """
    def __init__(self, data=None, next=None):
       self.data = data
       self.next = next

class Stack():
	""" A class representing a Stack """
	def  __init__(self, size=0, head=None):
		self.sizes = size
		self.head = head

	def push(self, i):
		self.sizes += 1
		if self.head == None:
			self.head = Node(i,None)
		else:
			n = Node()
			n.data = i
			n.next = self.head
			self.head = n
		return None
   
	def pop(self):
		if self.head is None:
			print "stack is empty"  
		else:
			self.sizes -= 1
			self.old_head = self.head
			self.head = self.head.next
			return self.old_head.data

	def top(self):
		if self.head is None:
			print "stack is empty"
		else:
			self.s = self.head.data
			return self.s

	def size(self):
		return self.sizes

	def isEmpty(self):
		if self.head is None:
			self.boo = True
		else:
			self.boo = False
		return self.boo


