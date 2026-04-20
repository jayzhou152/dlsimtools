# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:53:20 2020

@author: Jay
"""

import numpy as np
from . import SwitchBias
from random import randint


class GeneralOptimizer():

    def __init__(self, gettor, mutator, worker, range_checker, close_in_rule):
        
        self.get = gettor
        self.mut = mutator
        self.work = worker
        self.range_check = range_checker
        self.close_in_rule = close_in_rule
        
        
    def optimize(self, vrange, inc_method = "linear_half", inc_factor = 1, tol = 0.05, max_iter = 200):
        
        if self.range_check(vrange):
            # incur loop here after checking range is correct, applying incrementation method.
            
            for i in range(max_iter):
                
                if inc_method == "linear_half":
                    new_val = (vrange[0] + vrange[1])/2
                elif inc_method == "linear_quarter":
                    vals = [0.25,0.75]
                    new_val = vrange[0] + (vrange[1]-vrange[0])*vals[randint(0,1)]
                else:
                    print ("Incremental method does not exist.")
                
                # proceed to work out new range.
                
                self.mut(new_val)
                tar_dir = self.work()
                if not tar_dir:
                    out = self.get()
                else:
                    out = self.get(tar_dir)
                
                close_flag = self.close_in_rule(out)
                
                if close_flag == 1:
                    vrange[0] = vrange[0] + (new_val-vrange[0])*inc_factor
                elif close_flag == 2:
                    vrange[1] = vrange[1] - (vrange[1]-new_val)*inc_factor
                else:
                    print ("Range close flag is not recognised.")
                
                print(vrange)
                if vrange[1] == 0:
                    ref = vrange[1] + 0.0000000000000001
                else:
                    ref = vrange[1]
                if abs((vrange[0]-vrange[1])/ref) < tol:
                    print("successfully converged")
                    print(vrange)
                    return vrange
            return vrange
        else:
            print("range_check failed")
            return 0
