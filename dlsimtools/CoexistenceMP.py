#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .PolyCore import PolyCore
import os

"""
Class for creating coexistence systems.

"""

class CoexistenceMP ():
    
    dlpc = PolyCore()
    
    def __init__(self, k, tether_pot="harm",t_ratio = 0.5):
        pass
        #self.tether (k,tether_pot = tether_pot, t_ratio = t_ratio)
        #print ("{}% of supercell now tethered with {} type potential (k = {}).".format(int(t_ratio*100),tether_pot,k))
        
    def tether(self, k = 20, tether_pot = "harm", t_ratio = 0.5, filename = "FIELD"):
        
        tethered = []
        moltemp = []
        vdwtemp = []
        
        with open (filename,'r') as fh:
            for line in fh:
                if "Molecular types" in line:
                    # check if single component
                    moltnum = float(line.split()[2])
                    if moltnum != 1:
                        print ("More than 1 moltype in target FIELD file or \
                               field file is already tethered")
                        return 0
                    else:
                        line = line.replace(line.split()[2],"2")
                        tethered.append(line)
                        break
                tethered.append(line)
            for line in fh:
                moltemp.append(line)
                if "finish" in line:
                    break
            for line in fh:
                vdwtemp.append(line)
            
        for i in moltemp:
            if "Molecule name" in i:
                mname = i.split()[2]
                tname = i.split()[2] + "_teth"
                i = i.replace(mname,tname)
            if "nummols" in i:
                indh = moltemp.index(i)
                totnummol = float(i.split()[1])
                free_nummol = totnummol*(1-t_ratio)
                if not free_nummol.is_integer():
                    print("Total number of molecules is not wholly divisible by tether ratio")
                    return 0
                tethered_nummol = totnummol - free_nummol
                moltemp[indh] = i.replace(i.split()[1],str(int(free_nummol)))
                i = i.replace(i.split()[1],str(int(tethered_nummol)))
            if "atoms" in i:
                atmnum = int (i.split()[1])
            if "finish" in i:
                tethered_block = self.tether_block(1,atmnum,k,tether_pot)
                moltemp.extend(tethered_block)
                moltemp.append(i)
                break
            moltemp.append(i)
        
        tethered.extend(moltemp)
        tethered.extend(vdwtemp)
        
        with open ("FIELD_tethered",'w') as fw:
            for i in tethered:
                fw.write(i)
    
    def tether_block (self, start_ind, end_ind, k, tether_pot, rc=1):
        
        block = ["teth {}\n".format(end_ind-start_ind+1)]
        
        for i in range(start_ind,end_ind+1):
            
            teth_str = "{}        {}        {}\n".format(tether_pot,i,k)
            block.append (teth_str)
        
        return block
    
    def construct (self, temp = 800, pressure = 0.001, steps = 500000):
        
        print ("Constructing coexistence model...")
        
        ndir = self.dlpc.get_new_run("Coexi",files = ["CONFIG","CONTROL","FIELD_tethered"])
        os.chdir (ndir)
        self.dlpc.edit_control("temperature", temp)
        self.dlpc.edit_control("pressure", pressure)
        self.dlpc.edit_control("ensemble", "nvt")
        self.dlpc.edit_control("steps", steps)
        self.dlpc.run_poly(mode="mpi",np=2)
        
        print ("Please wait for model to reach coexistence...")
    
    def scan_temp (self, start_t, end_t, interval, tol = 1e-10):
        temps = list (range (start_t,end_t+1,interval))
        self.dlpc.edit_control("ensemble", "nst")
        self.dlpc.paralooper(temps, "temperature")