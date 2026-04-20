# -*- coding: utf-8 -*-

from .MonteCore import MonteCore
from .MonteData import MonteData
import numpy as np
import os
import shutil
import time
import subprocess

tar_mol = 114
mol_name = "MEOH"

def chempot_change(new_val):
    
    """
    Changes CONTROL file chemical potential value.
    
    Parameters
    ----------
    new_val : float,int
        New value of switch bias

    Returns
    -------
    None.

    """

    new_val = np.exp(new_val)

    with open('CONTROL','r') as fh:
        
        cont = []
        flag = 0
        for line in fh:
            if "move gcinsertmol" in line:
                flag = 1
            if mol_name in line and flag == 1:
                cur_val = line.split()[1]
                line = line.replace(cur_val,str(new_val))
                flag = 0
            cont.append(line)
        
    with open('CONTROL','w') as fw:
        for i in cont:
            fw.write(i)

def get_nummol (end_mode = "average"):
    md = MonteData()

    zero_count = 0
    same_count = 0
    high_count = 0
    low_count = 0
    prev_nummol = 0
    prev_n = 0 
    prev_n_count = 0

    print ("trying to seek chemical potential for {} {} molecules.".format(mol_name,tar_mol))

    with open ("jobid",'r') as fh:
        jobid = fh.readline()
    time.sleep(30)


    while True:
        
        time.sleep(2)
        nummol = md.last_yaml_data()[-1]
        nummol_array = md.yaml_one(-1)
        
        if len(nummol_array) <= 1000:
            continue

        if nummol == prev_nummol:
            same_count += 1
        else:
            same_count = 0
        if nummol == 0:
            zero_count += 1
        else:
            zero_count = 0
        if nummol > tar_mol:
            high_count += 1
        else:
            high_count = 0
        if nummol < tar_mol:
            low_count += 1
        else:
            low_count = 0
        
        prev_nummol = nummol
        #if nummol > (tar_mol*1.02):
        #    n_count += 1
        #else:
        #    n_count = 0
        
        
        av_nummol = np.average(nummol_array[1000:])
        
        if prev_n == len(nummol_array):
            prev_n_count += 1
        else:
            prev_n_count = 0

        prev_n = len(nummol_array)

        if (av_nummol < tar_mol*0.7 or av_nummol > tar_mol*1.3):
            if len(nummol_array) < 4000:
                continue
            else:
                print("N too high/low.")
                subprocess.run("kill {}".format(jobid),shell=True)
                return av_nummol

        if prev_n_count == 100:
            print("Sim ended.")
            subprocess.run("kill {}".format(jobid),shell=True)
            if end_mode == "average":
                return av_nummol 
            elif end_mode == "final":
                return nummol
            else:
                print("nummol output mode not recognised, outputting average value...")
                return av_nummol

        if (zero_count >= 10 or same_count == 1000 or high_count == 1500 or low_count == 1500 ):
            print("too many zeroes, same or high/low count")
            subprocess.run("kill {}".format(jobid),shell=True)
            return nummol

    
    
def decision(nummol):

    if nummol < tar_mol:
        return 1
    else:
        return 2


def run_test():
    time.sleep(1)

    subprocess.run("rm *.000",shell=True)
    proc = subprocess.Popen([r"DLMONTE-SRL.X","&"])

    with open ("jobid",'w') as fw:
        fw.write(str(proc.pid))
    


def chempot_range_check(init_range):

    if len(init_range) != 2:
        print ("Range must be a list/array of two values.")
        return 0 
    
    elif init_range[1] < init_range[0]:
        print ("Second value in range is larger than first, please check values were entered correctly.")
    
    else:
        chempot_change(init_range[0])
        tar_dir = run_test()
        nummol = get_nummol()
        
        if nummol < tar_mol:
            chempot_change(init_range[1])
            tar_dir = run_test()
            nummol2 = get_nummol()
            if nummol2 > tar_mol:
                return True
            else:
                print("Chempot range should be to the right of the current range, redefine range.")
                return False
            
        else:
            print("Chempot range should be to the left of the current range, redefine range.")
            return False
