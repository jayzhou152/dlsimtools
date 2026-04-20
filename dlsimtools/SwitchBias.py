# -*- coding: utf-8 -*-

from .MonteCore import MonteCore
import os
import shutil
import time


mc = MonteCore()
    
def bias_change(new_val, sig = 1):
    
    """
    Changes CONTROL file switch bias value.
    
    Parameters
    ----------
    new_val : float,int
        New value of switch bias

    Returns
    -------
    None.

    """
    
    if sig == 1:
        new_val = -new_val
    
    with open('CONTROL','r') as fh:
        
        cont = []
        for line in fh:
            if "switchbias" in line:
                cur_val = line.split()[1]
                line = line.replace(cur_val,str(new_val))
            
            cont.append(line)
        
    with open('CONTROL','w') as fw:
        for i in cont:
            fw.write(i)



def bias_mute (pval):
    
    if pval == 2:
        return 2
    elif pval == 1:
        return 1
    else:
        print ("Swtich_bias phase is not recognised")
        return 0
    
    
def get_phase(tar_dir):
    
    """
    Obtain phase at first PSDATA data point.

    Returns
    -------
    val : INT
        Phase for first PSDATA data point.

    """

    
    os.chdir(tar_dir)
    
    with open ('PSDATA.000','r') as fh:
        for i in range(1):
            fh.readline()
            
        line = fh.readline()
        val = int(line.split()[1])
    
    os.chdir("..")
    
    shutil.rmtree(tar_dir)
    
    return val

def run_sbtest():
    
    dir_name = mc.get_new_run(key="sbias")
    
    if dir_name == 0:
        return 0
    else:
        os.chdir(dir_name)
        mc.run_dlm(mode ="short")    
        os.chdir("..")
        return dir_name 

def range_check(init_range):
    
    if len(init_range) != 2:
        print ("Range must be a list/array of two values.")
        return 0 
    
    elif init_range[1] < init_range[0]:
        print ("Second value in range is larger than first, please check values were entered correctly.")
    
    else:
        bias_change(init_range[0])
        tar_dir = run_sbtest()
        phase1 = get_phase(tar_dir)
        
        if phase1 == 1:
            bias_change(init_range[1])
            tar_dir = run_sbtest()
            phase2 = get_phase(tar_dir)
            if phase2 == 2:
                return True
            else:
                print("Switch bias range should be to the right of the current range, redefine range.")
                return False
            
        else:
            print("Switch bias range should be to the left of the current range, redefine range.")
            return False
        
        

