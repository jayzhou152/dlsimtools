# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 11:17:07 2020

@author: Jay
"""


def get_num(file = "outputf"):
    
    with open (file,'r') as fh:
        line = fh.readline()
        num = int(line)
        
    return num

def edit_num(num, file = "inputf"):
    
    with open (file,'w') as fh:
        fh.write(num)

def working (infile = "inputf", outfile ="outputf"):
    
    with open (infile,'r') as fh:
        line = fh.readline()
        num = int(line)
        
    sqr = num**2
    
    with open (outfile,'w') as fw:
        fw.write(str(sqr))
        
