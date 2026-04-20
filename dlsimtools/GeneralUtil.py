# -*- coding: utf-8 -*-
"""
@author: Jay
"""


class GeneralUtil():
    
    def __init__(self):
        pass
    
    def edit_anything (self, key, val, file):
        
        cont = []
        val = str(val)
        with open (file, 'r') as fh:
            for line in fh:
                if key in line:
                    curval = line.split()[1]
                    line = line.replace(curval, val)
                cont.append(line)
        
        with open (file,'w') as fw:
            for i in cont:
                fw.write(i)
            
    def ev_to_K (self,ev):
        
        return ev*11604.525
    
    def get_anything (self,key,file):
        
        with open (file,'r') as fh:
            
            for line in fh:
                if key in line:
                    curval = line.split()[1]
                    return curval
    
    def lj_out (self, lj, temp, e_unit):
        
        Kelvins = lj * temp
        
        if e_unit == "ev":
            ev = Kelvins * 8.61732814974056E-05
            return ev
        
        elif e_unit == "J":
            joules = Kelvins * 1.38064878066922E-23
            return joules 
        
        elif e_unit == "K":
            return Kelvins
        
        else:
            print ("Incorrect e_unit.")
            raise ValueError
    