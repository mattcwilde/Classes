class PriorityQueue:
	
	def __init__(self, MAX=0):
		self.a = []
		self.N = 0
		self.max = MAX

	def min(self):
		"""Error code"""
		self.m = self.a[1]
		return self.m

	def size(self):
		return self.N

	def isEmpty(self):
		if len(self.a) == 0
			return True
		else return False

	def isFull(self):
		if 
		# ...
		return False

	def insert(self, i):
		self.i = i
		self.temp = 0
		self.k = 0
		self.N += 1 
		if self.N > self.max:
			a[self.N] = self.i
			return Error
		self.k = self.N
		while (self.k > 1) and (self.a[self.k/2] > self.a[self.k]):
			""" 
			as long as not at root and as long as parents value is > child
			do something 
				
			for child of index k:
			parent index is k/2 

			!!!!! integer division super important
			"""
			
			# exchange parent and child
			self.temp = self.a[self.k/2]
			self.a[self.k/2] = self.a[self.k]
			self.a[self.k] = self.temp
			self.k = self.k/2
		return None