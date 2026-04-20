""" DLPOLY Utility for handling STATIS files

STATIS FILE => DLPOLY output file containing high precision data of key variables

The STATIS class creates an object for transferring data in STATIS
file into an accessable format. 

In particular, this class is able to create one single dataframe 
containing all STATIS data with labelled column names with the 
loadall() function. If file size is large, the load() function 
is capable of fetching one data variable at a time to lower memory usage.
For visualisation, loadplot() method allows 2D plots of given variables.

valid variable keys are shown below, base keys are always present, stress tensor(st)
and npt(NPT ensemble) keys are valid depending on the DLPOLY setup.

base_keys = ["nstep","time","nument","engcns","temp","engcfg","engsrp","engcpe",
             engbnd","engang","engdih","engtet","enthal","tmprot","vir","virsrp","vircpe",
             "virbnd","virang","vircon","virtet","volume","tmpshl","engshl","virshl","alpha",
             beta","gamma","virpmf","press"]
st_keys = ["stress1","stress2","stress3","stress4","stress5",
           "stress6","stress7","stress8","stress9"]
npt_keys = ["cell1","cell2","cell3","cell4","cell5","cell6",
            "cell7","cell8","cell9"]
mean squared displacement for each atom type are keyed as "atomtypeN" where N is 
the abstract atomtype number. At the current state, to key in specific atom type names
is not possible since this info is not included in STATIS file, however, one could refer
to other output files as instructed by manual to relate atomtype numbers to atomtype names.

Currently defaults to current directory to locate STATIS file, however directory can be passed as
an argument to all three public functions in this class to locate this STATIS file. Similarly if
STATIS file was renamed, file name could also be passed to these functions as an input argument.

"""


import os
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt



class Statis:
    
    # Does not initiate any data, instead fetches data when called.
    def __init__(self):
        pass
    
    def keygen(self,n_ele,st,npt):
        
        """ returns list of keys based on the number of variable detected and 
            whether or not st and npt was enabled
            
            Arguments:
                n_ele (int) => number of variables detected
                st (int): 1 => enabled, otherwise => disabled
                npt (int): 1 => enabled, otherwise => disabled
            
            Returns:
                base_keys (list of strings): finalized list of keys
        """
        
        base_keys = ["nstep","time","nument","engcns","temp","engcfg","engsrp","engcpe",
                "engbnd","engang","engdih","engtet","enthal","tmprot","vir","virsrp","vircpe",
                "virbnd","virang","vircon","virtet","volume","tmpshl","engshl","virshl","alpha",
                "beta","gamma","virpmf","press"]
        st_keys = ["stress1","stress2","stress3","stress4","stress5",
                  "stress6","stress7","stress8","stress9"]
        npt_keys = ["cell1","cell2","cell3","cell4","cell5","cell6",
                    "cell7","cell8","cell9"]
        
        # natom => number of unique atom types
        natom = n_ele - 27 - st*9 - npt*9
        atom_keys = []
        
        # Check for the number of elements, if less than expected, terminate.
        if natom < 1:
            print ("data format incomplete")
            sys.exit(1)
        
        for atom in range(natom):
            atom_key = "atomtype" + str(atom+1)
            atom_keys.append(atom_key)
            
        base_keys.extend(atom_keys)
        
        # Append st and npt keys if enabled.
        if st == 1:
            base_keys.extend(st_keys)
        if npt == 1:
            base_keys.extend(npt_keys)
            
        return base_keys
    
    
    def loadall(self, directory=os.curdir, file="STATIS",st=0,npt=0):
        
        """ 
            Returns all STATIS data in a panda dataframe with labels
            
            Arguments:
                directory (string) => directory of STATIS file, defaults to curdir
                file (string) => file name, defaults to "STATIS"
            
            Returns:
                data (panda dataframe object, dtype=float64) => all data labeled
        """
        
        filename = os.path.join(directory, file)
        
        # Check file is in directory
        try:
            open(filename, "r")
        except FileNotFoundError:
            print("Cannot locate file " + file + " in current directory")
            sys.exit(1)
        
        # Fetch number of elements, block size (one block is the data in one timestep)
        # and keys.
        with open(filename, "r") as fh:
            for n in range(3):
                if n == 2:
                    conf=fh.readline().split()
                    n_ele = int(conf[2])
                    block_size = int(n_ele/5)+2
                    keys = self.keygen(n_ele,st,npt)
                    if len(conf)!=3:
                        print("bad file format")
                        sys.exit(1)
                else:
                    fh.readline()

        # Proceed to read line by line, creating block data lists and dump blocks at end of each timestep
        # block are appened to storage list before dumped
        # final block is made into dataframe and labelled.
        with open(filename, "r") as fh:
            
            global title 
            title = fh.readline()
            global e_u 
            e_u = fh.readline()
            
            count = 2
            block = []
            blocks = []

            for line in fh:
                count += 1
                block.extend(line.split())
                if count % block_size == 2:
                    blocks.append(block)
                    block = []
            data = pd.DataFrame(np.asarray(blocks,dtype=float))
            data.columns = keys
        return data
    
    def load(self, target, directory=os.curdir,file="STATIS",st=0,npt=0):
        """ Analogous to loadall(), except this class return one single column vector
            containing data for one specified variable.
            
            Arguments:
                target (string) => string of the target variable key (availabe keys shown in header)
            
            Returns:
                target_data (numpy array float64) => data of target variable
        """
        
        filename = os.path.join(directory, file)

        try:
            open(filename, "r")
        except FileNotFoundError:
            print("Cannot locate file " + file + " in current directory")
            sys.exit(1)
            
        with open(filename, "r") as fh:
            for n in range(3):
                if n == 2:
                    conf=fh.readline().split()
                    n_ele = int(conf[2])
                    block_size = int(n_ele/5)+2
                    keys = self.keygen(n_ele,st,npt)
                    if len(conf)!=3:
                        print("bad file format")
                        sys.exit(1)
                else:
                    fh.readline()
        
        # Check if key is valid.
        if target not in keys:
            print("targete variable is an invalid key")
            sys.exit(1)
        
        # Instead of gathering all block data, find data based on index of target key
        with open(filename, "r") as fh:
            global title 
            title = fh.readline()
            global e_u 
            e_u = fh.readline()
            
            count = 2
            block = []
            target_data = []
            target_index = keys.index(target)
            
            for line in fh:
                count += 1
                block.extend(line.split())
                if count % block_size == 2:
                    target_data.append(block[target_index])
                    block = []
            print(target_data)
            target_data = np.asarray(target_data, dtype=float)
            
        return target_data
    
    def loadplot(self,arg1,arg2,directory=os.curdir,file="STATIS",st=0,npt=0):
        
        """ Plots data for input keys. arg1 on x axis, arg 2 on y axis
        
            Arguments:
                arg1 (string) => string of x axis key
                arg2 (string) => string of y axis key
            
            Returns:
                none            
        """
        
        arg1_data = self.load(arg1,directory,file,st,npt)
        arg2_data = self.load(arg2,directory,file,st,npt)
        
        plt.plot(arg1_data,arg2_data)
        plt.title(title + e_u.strip() + "\n", loc='left')
        plt.xlabel(arg1)
        plt.ylabel(arg2)
        plt.show()

    