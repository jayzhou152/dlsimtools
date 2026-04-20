import numpy as np
import pandas as pd

"""
Contains DLPOLY OUTPUT file data methods. Mainly aimed at obtaining rolling averages.
STATIS.py would deal with discrete datapoint data. Which obviously requires a STATIS file to work.

OUTPUT KEYS:
    [
    "engtot","temptot","engcfg","engvdw","engcou","engbnd","engang","engdih","engtet",
     "engpv","temprot","vircfg","virvdw","vircou","virbnd", "virang","vircon","virtet",
     "volume","tempshl","engshl","virshl","alpha","beta","gamma", "virpmf","press"
    ]
    
"""

class PolyOutput():
    
    # 
    keys = ["engtot","temptot","engcfg","engvdw","engcou","engbnd","engang","engdih","engtet", 
            "engpv","temprot","vircfg","virvdw","vircou","virbnd", "virang","vircon","virtet", 
            "volume","tempshl","engshl","virshl","alpha","beta","gamma", "virpmf","press"]
    
    def __init__ (self):
        pass
    
    def get_rolling (self,file="OUTPUT",outputsteps = 250):
        
        keys = self.keys
        
        with open(file,'r') as fh:
            
            flag = 0
            block = []
            data = []
            
            for line in fh:
                if "rolling" in line:
                    flag = 1
                if flag == 1:
                    block.append(line)
                if len(block) == 3:
                    data.append(block)
                    block = []
                    flag = 0
                    
        reformat = []
        row = []
        for i in data:
            i[0]=i[0].replace("rolling","")
            i[1]=i[1].replace("averages","")
            for j in i:
                row.extend(j.split())
            reformat.append(row)
            row = []
        
        steps = [i*outputsteps for i in range(len(reformat))]
        data = pd.DataFrame(reformat,columns=keys,index= steps)
        return reformat, data
    
    def get_one (self, step, key, file = "OUTPUT"):
        try:
            temp,data = self.get_rolling ()
            if step == "end":
                return data.tail(1)[key]
            else:
                return data[key][step]
        except IndexError:
            return -1000000
    
    def get_av (self, key, file = "OUTPUT"):
        
        temp,data = self.get_rolling ()
        tar_s = np.asarray(data[key].values,dtype=float)
        av = np.mean(tar_s)
        
        return av