# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 18:26:17 2016

@author: mwilde


Program to test the linkedlist module

Module contains a Node class and LinkedList class.
"""

import linkedlist as ll

# create empty list
list = ll.LinkedList()
print "created empty list \n"


# add node to head
list.insert_at_head('foo')
print "added a node to the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# add a second node
list.insert_at_head('bar')
print "insert a second node to the head of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# add a third node
list.insert_at_head('spam')
print "insert a third node to the head of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# delete a node at the head
list.delete_at_head()
print "delete a node from the head of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# add a node to tail
list.insert_at_tail('green')
print "insert a node to the tail of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# add a second node to tail
list.insert_at_tail('eggs')
print "insert a 2nd node to the tail of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {} \n".format(list.tail.data,list.tail.next)

# add a third node to tail
list.insert_at_tail('ham')
print "insert a 3rd node to the tail of the linkedlist." 
print "head data is {} and head pointer = {}".format(list.head.data,list.head.next)
print "tail data is {} and tail pointer = {}".format(list.tail.data,list.tail.next)



