class PriorityQueue(object):

    def __init__(self, MAX=0):
        self.a = [0]*MAX
        self.N = 0
        self.max = MAX
        return None

    def min(self):
        # check if empty
        if self.isEmpty():
            print "empty PriorityQueue"
            return None
        
        m = self.a[1]
        return m

    def size(self):
        return self.N

    def isEmpty(self):
        if len(self.a) == 0:
            return True
        else:
            return False

    def isFull(self):
        if self.N == self.max+1:
            return True
        else:
            return False

    def insert(self, i):
        self.N = self.N + 1
        if self.N > self.max:
            self.k = self.N
            a[self.N] = i
            print "Error, N > max"
            return None
            k = self.N
        else:
            k = self.N
            self.a[self.N] = i
            while (k > 1) and (self.a[k/2] > self.a[k]):
                """
                heapsort O(NlogN)

                as long as not at root and as long as parents value is 
                > child do something 

                for child of index k:
                parent index is k/2 

                !!!!! integer division super important
                """

                # exchange parent and child
                temp = self.a[k/2]
                self.a[k/2] = self.a[k]
                self.a[k] = temp
                k = k/2
        return None
    
    def delMin(self):
        
        
        if self.isEmpty():
            print "Error, PriorityQueue already empty"
            return None
        else:
            # minimum value
            m = self.a[1]
            
            # put last value to top
            self.a[1] = self.a[self.N]
            self.a[self.N] = 0

            # loop backwards
            self.N = self.N - 1
            
            # start loop at root
            k = 1 
            while 2*k <= self.N:
                # Note: right child is 2k+1, left is 2k
                
                # left child
                if (2*k == self.N) or (self.a[2*k] < self.a[2*k+1]):
                    j = 2*k # lesser child or only child
                else:
                    j = 2*k + 1 # right child and is the smaller child
                    
                # See if heap has been ordered by checking if 
                # parent > child
                if (self.a[k] > self.a[j]):
                    temp = self.a[k]
                    self.a[k] = self.a[j]
                    self.a[j] = temp
                    k = j
                else:
                    # binary heap condition satisfied
                    # parent < child
                    break
            return m