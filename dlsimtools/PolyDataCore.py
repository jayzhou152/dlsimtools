import subprocess
import re
"""
DLPoly Utility for treating output data sets.

Available functionalities: concatenate sequel run results.

"""

import os
import numpy as np
from .PolyOutput import PolyOutput
from .STATIS import Statis
import pandas as pd
import matplotlib.pyplot as plt

class PolyDataCore:
    
    dlo= PolyOutput()
    sta = Statis()
    
    def __init__ (self):
        print ("dlpoly data core initiated, you are now able to operate on dlpoly datasets")
        pass
    
    def get_av_all (self, key, mode = "seq"):

        """
        Collect all single value averages from all runs.
        
        Arguments:
            key (str) ==> the key of the desired variable
            mode (str) ==> the key of the runs, defaults to seq (sequel runs)
        
        Returns:
            ?
        """
        if mode == "seq":
            folderkey = "dlprunseq"
        else:
            folderkey = "dlprun"
        
        out = []
        
        for i in os.listdir():
            if folderkey in i:
                os.chdir(i)
                av_val = self.dlo.get_av(key)
                regstr = r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
                val = re.findall(regstr,i.split('_')[1])
                o_string = "{}    {}\n".format(val[0],av_val)
                out.append(o_string)
                os.chdir("..")
        
        with open ("DLPOUT",'w') as fw:
            
            for i in sorted(out):
                fw.write(i)
    
    def concatenate (self, folderkey ="dlprunseq"):
        
        """
        Concatenate multiple sequel run STATIS, HISTORY and rolling average results.
        
        Arguments:
            folderkey (str) ==> the key to seek the DLP runs, defaults to "dlprunseq", which
            is the default name for sequel run looper written in "dlpcore".
        
        Returns:
            none.
        """
        
        STATIS_content = []
        HIS_content = []
        ROLL_content = []
        cwdir = os.getcwd()
        
        for i in sorted(os.listdir()):
            if folderkey in i:
                os.chdir(i)
                with open ("STATIS",'r') as fh:
                    for line in fh:
                        STATIS_content.append(line)
                with open ("HISTORY",'r') as fh:
                    for line in fh:
                        HIS_content.append(line)
                a,b = self.dlo.get_rolling()
                ROLL_content.extend(a)
                os.chdir(cwdir)
        
        with open ("STATISCON",'w') as fw:
            for i in STATIS_content:
                fw.write(i)
                
        with open ("HISTCON",'w') as fw:
            for i in HIS_content:
                fw.write(i)
        print(ROLL_content)
        with open("ROLLINGCON",'w') as fw:
            for i in ROLL_content:
                for j in i:
                    fw.write(j)
                    fw.write(" ")
                fw.write("\n")
    
    def read_rolling(self, file = "ROLLINGCON", keys=dlo.keys):
        
        """
        Reads the "ROLLINGCON" file and returns the results as a Pandas DataFrame.
        
        Arguments:
            file (str) ==> defaults to source file name "ROLLINGCON".
            keys (list(str)) ==> list of keys as column names for the DataFrame.
            
        Returns:
            data (pd.DataFrame) ==> Pandas dataframe containing rolling average data.
        
        """
        
        data = []
        
        with open("ROLLINGCON",'r') as fh:
            for line in fh:
                line = line.replace("\n","")
                data.append(line.split())
        
        index = [i*250 for i in range(len(data))]
        data = pd.DataFrame(data,columns=keys,index=index)
        
        return data
    
    def plot_rolling (self, data,keyone,keytwo):
        
        """
        plots rolling av from data produced by read_rolling which reads from 
        ROLLINGCON file.
        """
        dataone = np.asarray(list(data[keyone]),dtype=float)
        datatwo = np.asarray(list(data[keytwo]),dtype=float)
        plt.plot(dataone,datatwo,'x')
        plt.show()