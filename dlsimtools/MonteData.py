# -*- coding: utf-8 -*-
"""
@author: Jay
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from .GeneralUtil import GeneralUtil
import math
from dlmontepython.simtask import analysis
from decimal import Decimal
import shutil
import time

class MonteData():
    
    def __init__ (self, filename = "YAMLDATA.000"):

        self.filename = filename
    
    def yaml_data (self):
        
        """

        Parameters
        ----------
        file : STR, optional
            file to fetch yaml data from, must be in yaml data format. The default is "YAMLDATA.000".

        Returns
        -------
        Array
            Contains the YAML data as a numpy array.

        """
        
        data = []
        
        with open (self.filename, 'r') as fh:
            
            for i in range(1000):
                line = fh.readline()
                if "---" in line:
                    fh.readline()
                    break
            
            block = []
            for line in fh:
                
                if line.split()[0] == "-":
                    data.append(block)
                    block =[]
                    continue
                if '[' in line:
                    line = line.replace('[','')
                    line = line.replace(']','')
                if ',' in line:
                    line = line.replace(',',' ')
                vals = line.split(":")[1]
                if len(vals.split())== 1:
                    block.append(float(vals.split()[0]))
                else:
                    block.append([float(i) for i in vals.split()])
        
        return np.array(data,dtype=object)

    def get_bin_num(self):
        with open ('CONTROL','r') as fh:
            for line in fh:
                if "fed order param" in line:
                    lines = line.split()
                    return int(float(lines[4]))
    
    def last_yaml_data (self):
        
        """
        Get last yaml data record.

        Returns
        -------
        Array
            Array of the final yaml data record.

        """
        
        data = self.yaml_data()
        
        return data[-1]
    
    def last_yaml_ten (self):
        
        """
        Get last ten yaml data records.
        
        Returns
        -------
        Array
            Last ten yaml data records.

        """

        data = self.yaml_data()
        
        return data[-10:]
    
    def plot_yaml (self,ind1,ind2):
        
        """
        Get specific rows of data from yaml data and plot 2D graphs, e.g Time vs Orderparamter (time
        corresponds to the first row and order parameter corresponds to the third row, therefore
        use plot_yaml (0,2) to plot.

        Parameters
        ----------
        ind1 : INT
            Indices of the first desired yaml data value.
        ind2 : TYPE
            Indices of the second desired yaml data value.

        Returns
        -------
        None.

        """
        

        data = self.yaml_data()
        xd = data[:,ind1]
        yd = data[:,ind2]
        
        plt.plot(xd,yd,"x")

    def last_yaml_ten_percent (self):
        
        data = self.yaml_data()
        length = int(len(data)/2)
        
        return data[-length:]
    

    def last_yaml_n_percent (self,n):

        data = self.yaml_data()
        length = int(len(data)*n)
        
        return data[-length:]

    def get_tmat (self,file = "TMATRX.000",quick=True):
        
        """
        Fetch transition matrix data from TMATRX.000 and outputs it as a numpy array.

        Parameters
        ----------
        file : STR, optional
            File name to fetch tmatrix from, defaults to TMATRX.000, must be in transition matrix form. The default is "TMATRX.000".

        Returns
        -------
        Array
            Transition matrix info.

        """
        mat = []
        try:
            with open (file, 'r') as fh:
                fh.readline()
                fh.readline()
                
                for i in fh:
                    try:
                        row = list(map(float, i.split()))
                    except ValueError:
                        row = []
                        for j in i.split():
                            try:
                                row.append(float(j))
                            except ValueError:
                                row.append(float(0))
                    mat.append(row)
        except (FileNotFoundError,UnicodeDecodeError) as e:
            print ("Tmat not found in {}...".format(os.getcwd()))
        return np.array(mat,dtype=object)  
    
    def yaml_one (self, ind):
        
        data = self.yaml_data()
        one = data[:,ind]
        return one

    def ps_data (self, file = "PSDATA.000"):
        cont = []
        
        with open (file,'r') as fh:
            for line in fh:
                line = line.split()
                line = [float(i) for i in line]
                cont.append(line)
        cont = np.array(cont)
        cont.astype(float)
        
        return cont
    
    def ps_data_col (self, col, plot=False):
        
        data = self.ps_data()
        
        if plot:
            plt.plot(data[:,col])
            
        return data[:,col]


    def calc_pratio(self,tmin,tmax):
        
        data = self.yaml_data()
        
        sumone = 0
        sumtwo = 0

        for i in data:
            if i[0] < tmin or i[0] > tmax:
                continue
            else:
                if i[3] == 1:
                    sumone += Decimal(i[1]).exp()
                else:
                    sumtwo += Decimal(i[1]).exp()
        print(sumone/sumtwo)
        return sumone/sumtwo
    
    def calc_free(self,temperature,tmin,tmax):

        kb = 8.617333262145E-5 #ev/K
        pratio = self.calc_pratio(tmin,tmax)
        fe = -kb*float(temperature)*float(Decimal.ln(pratio))
        return fe

    def softedge(self, file = 'CONTROL'):
        try:
            with open (file,'r') as fh:
                for line in fh:
                    if "softedges" in line:
                        if len(line.split()) > 11:
                            return True
                        else:
                            return False
                return False  
        except FileNotFoundError:
            return False


    def fe_analysis(self, prodn, tcut, detect_se = False, write_out = True):
        
        if prodn == "auto":
            prod_list = []
            for i in os.listdir():
                if "lsmc_prod" in i:
                    os.chdir(i)
                    print(i)
                    if detect_se:
                        if self.softedge():
                            if self.check_prod():
                                prod_list.append(i)
                    else:
                        if self.check_prod:
                            prod_list.append(i)
                    os.chdir("..")
            prod_list = sorted(prod_list)
            print(prod_list)
        else:
            if not isinstance(prodn,(list,str)):
                print("prodn is not a list or 'auto' (automated folder fetch). prodn must \
                contain two numbers, for the start and end of production run folder number \
                respectively.")
                return 0
            
            if prodn[0] > prodn[1]:
                print("First number in prodn must be smaller than the second numner \
                in prodn.")
                return 0
            

            prod_list = []
            for i in range(prodn[0],prodn[1]+1):
                tdir = "lsmc_prod{}".format(i)
                os.chdir(tdir)
                if self.check_prod():
                    prod_list.append(tdir)
                os.chdir("..")

        gu = GeneralUtil()
        temperature = gu.get_anything("temperature","CONTROL")
        fes = []

        for i in prod_list:
            print("currently working on folder {}".format(i))
            os.chdir(i)
            time_series = self.yaml_one(0)
            m_series = self.yaml_one(2)

            if tcut == "auto":
                tmin_pos = int(analysis.autocorrelation_time(m_series))
                tmin = time_series[tmin_pos-1]
            else:
                tmin = int(max(time_series)*tcut)
            tmax = max(time_series)
            print(tmin,tmax)
            fe = self.calc_free(temperature, tmin, tmax)
            fes.append(fe)
            print(fe)
            os.chdir("..")

        se = np.std(fes)/math.sqrt(len(fes))
        mean = np.average(fes)

        if write_out:
            with open ("fe_report.dat",'w') as fw:
                fw.write("Free energy report.\n")
                for i in range(len(fes)):
                    fw.write(prod_list[i])
                    fw.write(":    ") 
                    fw.write(str(fes[i])+"eV\n")
                fw.write("Free energy mean:    {}eV\n".format(mean))
                fw.write("Free energy standard error:    {}eV\n".format(se))
        
        return fes,mean,se

    def check_prod(self, cut_ratio = 0.5, pratio = 0.98):
        
        from MonteCore import MonteCore
        mc = MonteCore()
        bgap, abound, wbound = mc.get_bandgap()

        try:
            data = np.loadtxt("PSDATA.000")
        except OSError:
            return False
        
        data = data[(int(len(data)*cut_ratio)):,3]
        data_new = [x for x in data if x <= abound[1]]

        if max(data_new) - min(data_new) >= pratio * bgap:
            return True
        else:
            return False

    def tri_check (self, file= "CONTROL"):
        
        with open(file, 'r') as fh:
            for line in fh:
                if "fed method" in line:
                    if "tri" in line:
                        return True
                    else:
                        return False
        
        return False

    def master_report_plot(self, plot = True):
        data_col = []
        with open('master_report','r') as fh:
            fh.readline()
            for i in fh:
                data = [float(x) for x in i.split()]
                data_col.append(data)
        data_col = np.array(data_col)
        plt.errorbar(data_col[:,0], data_col[:,1],yerr=data_col[:,2])
        return data_col

    def get_ncount (self):
        count = 0

        for i in os.listdir():
            if "tmmc" in i:
                count += 1
        return count

    def fep_all(self,keyone,keytwo,keythree,ncount):

        for i in os.listdir():
            if keyone in i :

                os.chdir(i)
                for j in os.listdir():
                    if keytwo in j:
                        os.chdir(j)
                        ncount = self.get_ncount()
                        if "fe_prof" in os.listdir():
                            os.rename("fe_prof","fe_prof_back")
                        self.tmat_profile(ncount, folderkey=keythree, valwrite=True)
                        os.chdir("..")
                os.chdir("..")
    
    def tmat_profile (self, nwin, folderkey = "win", valplot = False, valwrite = False, make_back= False):
        
        """
        Fetches transition matrix data from windowed runs and calculates the resulting bias profile.

        Parameters
        ----------
        nwin : INT
            Number of windows to fetch transition matrix info.
        folderkey : STR, optional
            Folder key to search for, refer to get_newrun function in core module. The default is "lsmc".
        valplot : BOOLEAN, optional
            True enables plots of bias versus bin number. The default is False.

        Returns
        -------
        Array
            Tmatrix bias vector.

        """

        
        for i in np.arange(nwin):
            dirname = "{}{}".format(folderkey,i+1)
            os.chdir(dirname)

            if i == 0:
                
                tmat = self.get_tmat()
                bin_n = self.get_bin_num()
                while len(tmat) != bin_n:
                    time.sleep(10)
                    tmat = self.get_tmat()

                if make_back:
                    shutil.copy("TMATRX.000","TMATRX_back.000") 
                print("TMATRX SIZE:")
                print(len(tmat))
            else:
                count = 0
                while True:
                    if count >= 180 and "TMATRX_back.000" in os.listdir():
                        tmat_cur = self.get_tmat(file="TMATRX_back.000")
                        print("TMAT size mismatch in window{}, using TMAT_BACK INSTEAD".format(i+1))
                        break
                    else:
                        if count > 600:
                            raise RuntimeError("TMAT writing complete in window {} and correct bin number  not detected in 10 minutes. Terminating checker...".format(i+1))
                            break
                    tmat_cur =  self.get_tmat()
                    if len(tmat_cur) == bin_n:
                        if make_back:
                            shutil.copy("TMATRX.000","TMATRX_back.000") 
                        break
                    time.sleep(1)
                    count += 1
                
                tmat += tmat_cur


            os.chdir("..")
        bias_vec = np.zeros(len(tmat))
        if self.tri_check():        
            for i in np.arange(len(bias_vec)-1):
                bias_vec[i+1] = bias_vec[i] + np.log(((tmat[i,2]+1)/(tmat[i+1,0]+1)) *\
                ((sum(tmat[i+1])+1)/(sum(tmat[i])+1)))
        else:
            for i in np.arange(len(bias_vec)-1):
                bias_vec[i+1] = bias_vec[i] + np.log(((tmat[i,i+1]+1)/(tmat[i+1,i]+1)) *\
                ((sum(tmat[i+1])+1)/(sum(tmat[i])+1)))            
        bias_vec = bias_vec - max(bias_vec)
        f_e = -bias_vec
        
        if valplot:
            plt.plot(np.arange(len(tmat)),f_e)
        
        if valwrite:
            with open("fe_prof",'w') as fw:
                for i in f_e:
                    fw.write(str(i)+"\n")
                
        
        return bias_vec
         
