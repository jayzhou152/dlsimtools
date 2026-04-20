#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:48:50 2020

@author: jay
"""
import subprocess
import numpy as np
import os

class MetaUtil():
    
    def __init__(self):
        pass
    
    def run_metadise(self, meta_com = "metadise"):

        FNULL = open(os.devnull,'w')
        subprocess.run(meta_com,stdout = FNULL, stderr = subprocess.STDOUT)
        
    def get_title_vec (self, file = "config__o0001sl.dlp"):
        
        with open (file, 'r') as fh:
            title = fh.readline()
        title = title.split()
        title = np.asarray (title[1:4], dtype = float)
        
        return title
    
    def edit_grow(self, new_growth, file = "input.txt"):
        
        cont = []
        
        with open (file,'r') as fh:
            
            for line in fh:
                if "grow" in line:
                    old_growth = line.split()[1:4]
                    line = line.split()[0] + " " + new_growth + "\n"
                    
                cont.append(line)
        
        with open (file,'w') as fw:
            for i in cont:
                fw.write(i)
        
        old_gstr = "{} {} {}".format(old_growth[0],old_growth[1],old_growth[2])
        
        return old_gstr
    
    def edit_molrad (self, newmolrad, file = "input.txt"):
        
        cont = []
        
        with open (file,'r') as fh:
            
            for line in fh:
                if "molecule" in line:
                    line = "molecule radall {}\n".format(newmolrad)
                    
                cont.append(line)
        
        with open (file,'w') as fw:
            for i in cont:
                fw.write(i)
        
    
    def fres_check (self, atmnum, file = "mol_o0001sl_res.out"):
        
        mol_count = 0
        flag = True
        
        with open (file, 'r') as fh:
            
            for line in fh:
                if "basi" in line:
                    break
                
            for line in fh:
                if "#molecule" in line:
                    atmcount = int (line.split()[2])
                    mol_count += 1
                    
                    if atmcount != atmnum:
                        flag = False
        
        return flag, mol_count, atmcount
    
    def molrad_check (self, molrad,atmnum):
        
        self.edit_molrad(molrad)
        self.run_metadise()
        flag, mol_count, atmcount = self.fres_check(atmnum)
        subprocess.run ("rm *000*", shell=True)
        
        return flag, mol_count, atmcount
    
    def molrad_optimiser (self, atmnum, rrange = [0.2,1], tol = 0.001):
        
        #with open ('input.txt', 'r') as fh:
        #    line = fh.readline()
        #if line.split()[0] != "dlpoly":
        #    print ("Currently only accepts dlpoly config files.")
        #    return None
        
        if len(rrange) != 2:
            print ("Specified optimiser range must contain only two entries.")
            return None
        
        try:
            rrange = list(map(float,rrange))
            if rrange[0] >= rrange[1]:
                print ("Specified optimiser range has its starting value \
                       larger than its finishing value.")
                return None
        
        except ValueError:
            print ("Specified optimiser range \
                   must contain numbers only.")
            return -1
        
        print ("Value check complete, looking to optimise molecular seek distance \
               for optimal molecular separation.")
        flag1, mol_count1, atmcount1 = self.molrad_check(rrange[0], atmnum)
        if flag1:
            return rrange[0]
        else:
            if atmcount1 > atmnum:
                print ("Lower bound gave low molecular numbers, try lower range.")
                return None
        flag2, mol_count2, atmcount1 = self.molrad_check(rrange[1], atmnum)
        if flag2:
            return rrange[1]
        else:
            if atmcount1 < atmnum:
                print ("Upper bound gave too many molecules, try higher range.")
                return None
        print ("Range check complete, convergence within bound, calcs starting...")
        

        while (rrange[1] - rrange[0]) > tol:
            val = (rrange[0] + rrange[1]) /2
            flag, mol_count, atmcount = self.molrad_check(val, atmnum)
            if flag:
                print ("Sucessfully converged to a molecule radius range of {} Angstrom.".format((rrange[0] + rrange[1]) /2))
                return val
            else:
                if atmcount < atmnum:
                    rrange[0] = val
                else:
                    rrange[1] = val
        
        print ("Converged to a molecule radius range of {} Angstrom. Did not meet tolerance.".format((rrange[0] + rrange[1]) /2))

        return (rrange[0] + rrange[1]) /2
    
    def add_line(self, line, trigger = "check", file = "input.txt"):
        
        with open(file,'r') as fh:
            cont = []
            for i in fh:
                cont.append(i)
                if trigger in i:
                    cont.append(line + "\n")
        
        with open(file, 'w') as fw:
            for i in cont:
                fw.write(i)
                

        
        
        
        
        
        