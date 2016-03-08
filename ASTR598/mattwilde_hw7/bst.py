class Node(object):
    def __init__(self,data=None,left=None,right=None):
        self.data = data # key->critical for deciding how BST behaves
        self.left = left # leaves
        self.right = right # leaves

class BinarySearchTree(object):
    """
    A binary search tree: left leaves are always greater than corresponding right leaf

        10 (root)
      /    \
     5      14  (leaves)
    /  \    /  \
   1   7 12   18
              /
             15

    
    """
    def __init__(self,N=0,root=None):
        """ Initialize with Length N and root = emtpy node. """
        self.N = N
        self.root = root
    
    # Given a local root, insert n in correct place recursively
    # By traversing down tree until the correct place is found
    def insert_rec(self,local_root,n):
        
        if local_root == None:
            print "Insert_rec error: local_root==None"
            return -1
        
        if n.data < local_root.data:
            if local_root.left == None: # No leaf, insert here!
                local_root.left = n
            else: # Left child exists, move down the tree
                self.insert_rec(local_root.left,n)
        else:
            if local_root.right == None: # No leaf, insert here (if n.data > local_root.data)
                local_root.right = n
            else: # Right child exists, move down the tree
                self.insert_rec(local_root.right,n)
    
    
    
    def insert(self, value):
        """ Function to insert nodes into BST.
            
            input: node value
            
            returns: None
        """
        
        # Create new node
        n = Node(data=value,left=None,right=None) # explicit with instantiation for note taking purposes
        
        # Case: Empty BST
        if self.root == None:
            self.root = n
        else:
            self.insert_rec(self.root,n) # recursive insert function
        
        # Increase size of BST
        self.N = self.N + 1
        
    
    # Check if value is in the BST
    def find(self,value):
        return self.find_rec(self.root,value) # Call recursive find
    
    
    # Recursive find method that traverses tree looking for value
    def find_rec(self,local_root,value):
        
        # Reached the end, didn't find value
        if local_root == None:
            return False
        # Found value!
        elif local_root.data == value:
            return True
        # Recursive step: pick right leaf to look at
        # Left value is greater than value, look there
        elif value < local_root.data:
            return self.find_rec(local_root.left,value)
        # Right value is greater, look there
        else:
            return self.find_rec(local_root.right,value)
        
    
    # Go down the BST in order printing everything along the way
    def traverse(self):
        self.traverse_in_order(self.root)
        
    
    # More complicated recursive traverse
    # First process the parent, the left subtree, then right subtree
    def traverse_pre_order(self,local_root):
        # Check for null (base) case
        if local_root == None:
            return
        else:
            print local_root.data
            # Now go through rest of the tree recursively
            self.traverse_pre_order(local_root.left) # Take care of left subtree
            self.traverse_pre_order(local_root.right) # Take care of right subtree
            
            
    # Traverse tree in order
    # Leverages the fact that left.data < local_root.data < right.data
    def traverse_in_order(self,local_root):
        # Check for null (base) case
        if local_root == None:
            return
        else:
            # Process left subtree
            self.traverse_in_order(local_root.left) # Take care of left subtree
            # Print parent's value
            print local_root.data
            # Now process right subtree
            self.traverse_in_order(local_root.right) # Take care of right subtree
            
    
    # Print current node after left, right 
    def traverse_post_order(self, local_root):
        # Check for null (base) case
        if local_root == None:
            return
        else:
            # Process left subtree
            self.traverse_post_order(local_root.left) # Take care of left subtree
            # Now process right subtree
            self.traverse_post_order(local_root.right) # Take care of right subtree
            # Print parent's value
            print local_root.data
    
    
    # Return minimum value of BST
    def mininum(self):
        # Base case: empty tree
        if self.root == None:
            return None
        
        local_root = self.root
        i = 0
        
        while local_root != None:
            i = local_root.data
            local_root = local_root.left
        
        return i
    
    
    # Return maximum value of BST
    def maximum(self):
        if(self.root == None):
            return None
        
        local_root = self.root
        i = 0
        
        while local_root != None:
            i = local_root.data
            local_root = local_root.right
        
        return i
            
    
    def isEmpty(self):
        if self.root == None or self.N == 0:
            return True
        else:
            return False
    
    def size(self):
        return self.N
