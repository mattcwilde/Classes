# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 17:38:10 2016

@author: mwilde
"""

class Node:
    """ 
    An element in a linked list with a data attribute
    and a pointer attribute. Initializes to None.
    """
    def __init__(self, data=None, next=None):
       self.data = data
       self.next = next
      
      
class LinkedList:
    """ Create a class where the input are two nodes.

        Initializes to an empty LinkedList.
    """
    def __init__(self, head=None, tail=None):
        self.head = head
        self.tail = tail
    
    def insert_at_head(self,i):
        if self.head == None:
            self.head = Node(i,None)
            self.tail = self.head
        elif self.tail == None:
            n = Node()
            n.data = i
            n.next = self.head
            self.head = n
            self.head.next = self.tail
        else:
            n = Node()
            n.data = i
            n.next = self.head
            self.head = n
       
    def delete_at_head(self):
        """Make cases for empty list and len = 1 list"""
        
        if self.head is None:
            pass
            # linked list is empty. cannot delete.
        elif self.head == self.tail:
            # Length of linked list is one. Deleting list.
            self.old_head = self.head
            self.head = None
            self.tail = None     
        else:
            self.old_head = self.head
            self.head = self.head.next
            return self.old_head

    def insert_at_tail(self, i):
        
        n = Node(i, None)
        if self.head is None:
            self.head = n
            self.tail = n
        elif self.head == self.tail:
            self.tail.next = n
            self.tail = n
        else:
            self.tail.next = n
            self.tail = n