# -*- coding: utf-8 -*-
"""
@author: Jay
"""
from .MetaSurf import MetaSurf
from .MetaUtil import MetaUtil
from .MonteCore import MonteCore
from .MonteData import MonteData
from .GeneralUtil import GeneralUtil
from .InputConverter import InputConverter
import numpy as np
import shutil
import os
import subprocess



class LSMC():
    
    
    def __init__ (self, dlm_com=None, dlm_com_par=None):

        """
        CONSTRUCTOR to create function objects.

        Returns
        -------
        None.

        """

        mc_kwargs = {}
        if dlm_com is not None:
            mc_kwargs["dlm_com"] = dlm_com
        if dlm_com_par is not None:
            mc_kwargs["dlm_com_par"] = dlm_com_par
        self.mc = MonteCore(**mc_kwargs)
        self.ms = MetaSurf()
        self.mcd = MonteData()

    def run_metagen (self, surflist, atmnum, growth=[1,1,1], vac_space = 20, code_s = "auto", **kwargs):
        """Generates DL_MONTE inputs using miller indices.

        Args:
            surflist (list): list of miller indices 
            atmnum (int): number of atoms per molecule for solid.
            growth (list, optional): growth multiple in 3d space. Defaults to [1,1,1].
            vac_space (int, optional): vacuum space length. Defaults to 20.
            code_s (str, optional): unit cell stack code, "auto" defaults to DL_POLY based selections. Defaults to "auto".
            **kwargs:
            steps (int, optional): steps for energy relaxation during for "relaxed" mode. Defaults to 500000.
            cut_off (list, optional): minimum allowed dimension length [xmin,ymin]. Defaults to [20,20].
            mode (str, optional): single point energy("singp") for lowest energy unit cell, or "relaxed" energy \
            for lowest energy unit cell. Defaults to "singp".
            ff (str, optional): force field for solid. Defaults to "opls2005".
            solvent (str, optional): force field for solvent. Defaults to "none".
        """
        mu = MetaUtil()
        gu = GeneralUtil()
        ic = InputConverter()
        stacks = {}
        
        for i in surflist:
            if i in os.listdir():
                shutil.rmtree(i)
            tar_code = self.ms.prod_meta(i, atmnum, **kwargs)
            with open ("misc",'w') as fw:
                fw.write("tar_code: {}\n".format(tar_code))
            if code_s == "auto":
                os.chdir(tar_code)
            else:
                os.chdir(code_s)
            stack,growth,molcount = self.stack_info()
            stacks[i] = (stack,growth,molcount)
            nummol = int(float(stack[0])*float(growth[0])*float(growth[1])*molcount)
            bulk_vec = mu.get_title_vec()
            self.slab2bulk(bulk_vec)
            ic.cut_field(nummol)
            self.mc.convert_pm(uc=True)
            self.bns_config(vac_space=vac_space)
            gu.edit_anything("NCONFIG","2","FIELD_DLM")
            print ("Metagen finished and DLMONTE files for PSMC is now prepared")
            os.chdir("../..")
    
    def change_cells (self, mat1, mat2, turn_frac = True):
        
        self.frac_carti_swap(file = "CONFIG.1")
        self.frac_carti_swap(file = "CONFIG.2")
        self.mc.change_cellmat(mat1,file = "CONFIG.1")
        self.mc.change_cellmat(mat2,file = "CONFIG.2")
        
        with open ("CONFIG",'w') as fw:
            for i in ["CONFIG.1","CONFIG.2"]:
                with open(i) as fh:
                    for j in fh:
                        fw.write(j)
                
    def frac_carti_swap (self, file = "CONFIG", flag = 0):
        
        if os.path.isfile(file):
            pass
        else:
            raise FileNotFoundError("Frac2Carti Target file not found.")
            
        if flag == 0 :
            val = 1
        else:
            val = 11
            
        met_inp_str = "dlmonte\n\
                       @include {}\n\
                       ends\n\
                       nopoten\n\
                       nosort\n\
                       noshift\n\
                       print dlmonte {}\n\
                       check\n\
                       stop\n\
                       start\n\
                       stop\n".format(file,val)
        
        with open ("input.txt",'w') as fw:
            fw.write(met_inp_str)
        
        subprocess.run("metadise")
        os.remove(file)
        os.rename("config__o0001no.dlp",file)
        os.remove("input.txt")
        for i in os.listdir():
            if "0001" in i:
                os.remove(i)
                
            
    def run_rangeseek (self, surflist, temp = 300, pressure = 0):
        
        """
        Run two simulations to find the lower and upper bound of the PSMC order parameter for specified surfaces.
        

        Parameters
        ----------
        surflist : STR LIST
            List of surface Miller indices to run rangeseek on. Each char must be
            separated by single spaces.
        temp : FLOAT, optional
            System temperature. The default is 300.
        pressure : FLOAT, optional
            System pressure. The default is 0.

        Returns
        -------
        None.

        """
        
        for i in surflist:
            os.chdir("surface_{}".format(i))
            self.ms.direct_tar()
            os.chdir(self.mc.get_new_run(key="rangeone"))
            self.mc.create_psmc_control("-1E20","1E20", temp=temp, pressure = pressure,cartdisp=True)
            self.mc.run_dlm()
            print ("Simulation looking for order parameter bound in phase 1...")
            os.chdir("..")
            os.chdir(self.mc.get_new_run(key="rangetwo"))
            self.mc.create_psmc_control("-1E20","1E20", temp=temp, pressure = pressure,init=2,cartdisp=True)
            self.mc.run_dlm()
            print ("Simulation looking for order parameter bound in phase 2...")
    
    def npt_relax(self, update = True, mpi=True, n = 8, **kwargs):

        nptfname1 = self.mc.get_new_run(key = "vrelax")
        os.chdir(nptfname1)
        self.mc.create_psmc_control("-1E10","1E10", init=1, **kwargs)
        if not mpi:
            self.mc.run_dlm()
        else:
            self.mc.run_dlm("mpi", n = min(n,32))
        os.chdir("..")
        #nptfname2 = self.mc.get_new_run(key = "vrelax")
        #os.chdir(nptfname2)
        #self.mc.create_psmc_control("-1E10","1E10", init=2, **kwargs)
        #if not mpi:
        #    self.mc.run_dlm()
        #else:
        #    self.mc.run_dlm("mpi", n = n)
        #os.chdir("..")
        self.mc.check_runs_terminate([nptfname1],20)

        os.chdir(nptfname1)
        matdat = np.vstack(self.mcd.last_yaml_n_percent(0.3)[:,-1])
        bulkz = np.average(matdat[:,0])
        bulkx = np.average(matdat[:,1])
        bulky = np.average(matdat[:,2])
        bulkarea =  bulkx * bulky
        bulkvol = bulkz * bulkx * bulky
        bulkmat = np.array([[bulkz,0,0],[0,bulkx,0],[0,0,bulky]])
        os.chdir("..")

        #os.chdir(nptfname2)
        #matdat = np.vstack(self.mcd.last_yaml_n_percent(0.3)[:,-1])
        #surfz = np.average(matdat[:,0])
        #surfx = np.average(matdat[:,1])
        #surfy = np.average(matdat[:,2])
        #surfarea =  surfx * surfy
        #surfmat = np.array([[surfz,0,0],[0,surfx,0],[0,0,surfy]])
        #os.chdir("..")

        if update:
            
            confonecont = self.mc.change_cellmat(bulkmat, file = "CONFIG.1")
            conftwocont = self.mc.change_cellmat(bulkmat, file = "CONFIG.2", mode = "surf")
            
            with open ('CONFIG',"w") as fw:
                for i in confonecont:
                    fw.write(i)
                fw.write('\n')
                for i in conftwocont:
                    fw.write(i)
        
        return bulkz,bulkx,bulky,bulkarea,bulkvol


    def rangeseek(self, temperature = 300, pressure = 0, steps = 100000, mpi = True, n = 8, rigid = False):
        
        print("Initiating range seek...")
        
        fnames = []
        nameone = self.mc.get_new_run(key = "rangeone")
        os.chdir(nameone)
        fnames.append(nameone)
        if rigid == True:
            self.mc.create_psmc_control("-1E20","1E20", temp=temperature, pressure = pressure, steps = steps, moveatm=False, movemol=True,)
        else:
            self.mc.create_psmc_control("-1E20","1E20", temp=temperature, pressure = pressure, steps = steps)
        
        if not mpi:
            self.mc.run_dlm()
        else:
            self.mc.run_dlm("mpi", n = min(int(n/2),32))
        
        print("Lower bound simulation is running.")
        
        os.chdir("..")
        nametwo = self.mc.get_new_run(key = "rangetwo")
        os.chdir(nametwo)
        fnames.append(nametwo)
        if rigid == True:
            self.mc.create_psmc_control("-1E20","1E20", temp=temperature, pressure = pressure, steps = steps, init =2, moveatm=False, movemol=True,)
        else:
            self.mc.create_psmc_control("-1E20","1E20", temp=temperature, pressure = pressure, steps = steps, init =2)
        if not mpi:
            self.mc.run_dlm()
        else:
            self.mc.run_dlm("mpi", n = min(int(n/2),32))
        
        os.chdir("..")
        print("Upper bound simulation is running.")
        
        return fnames
    
    def run_tmmc_bias (self, lbound, ubound, nbin, nwin, surf, temp = 300, pressure = 0):
        
        """
        Runs transition matrix biasing method with windows for specified surface after metagen.
        If none metagen input files are present, please name it CONFIG_MASTER, CONFIG.1, CONFIG.2, FIELD_DLM respectively
        and specify surface string as "here".
        
        
        Parameters
        ----------
        lbound : FLOAT
            Lower bound of the transition matrix order parameter.
        ubound : FLOAT
            Upper bound of the transition matrix order parameter.
        nbin : INT
            Number of total bins in the simulation.
        nwin : INT
            Number of total windows in the simulation.
        surf : STR
            Surface miller index "h k l" to run the simulation, must run metagen beforehand to
            genereate relevant files (or type in "here" to use files in current folder)
        temp : FLOAT, optional
            System temperature. The default is 300.
        pressure : FLOAT, optional
            System pressure. The default is 0.

        Returns
        -------
        None.

        """
        
        if surf != "here":
            os.chdir("surface_{}".format(surf))
            self.ms.direct_tar()
        
        self.mc.tmmc(lbound,ubound,nbin,nwin,temp,pressure)
        os.chdir("../..")
        
    def run_bias (self, lbound, ubound, surf, temp = 300, pressure = 0, mode = "ee"):
        """
        Runs none windowed biasing methods for specified surface after metagen.
        If none metagen input files are present, please name it CONFIG_MASTER, CONFIG.1, CONFIG.2, FIELD_DLM respectively
        and specify surface string as "here".       

        Parameters
        ----------
        lbound : FLOAT
            Lower bound of the FE order parameter.
        ubound : FLOAT
            Upper bound of the FE order parameter.
        surf : STR
            Surface miller index "h k l" to run the simulation, must run metagen beforehand to
            genereate relevant files (or type in "here" to use files in current folder)
        temp : FLOAT, optional
            System temperature. The default is 300.
        pressure : FLOAT, optional
            System pressure. The default is 0.
        mode : STR, optional
            FE method to use, currently only accepts "ee". The default is "ee".

        Returns
        -------
        None.

        """
        if surf != "here":
            os.chdir("surface_{}".format(surf))
            self.ms.direct_tar()
        
        if mode == "ee":
            ndir = self.mc.get_new_run(key = "ee")
            os.chdir(ndir)
            self.mc.create_psmc_control(lbound,ubound, temp = temp, pressure=pressure, step = 10000000, freq=500)
            self.mc.run_dlm()
            os.chdir("..")
        
    def fep_to_feddat (self, data):
        
        if data == "fep":
            data = []
            try:
                with open ("fe_prof",'r') as fh:
                    for line in fh:
                        dataf = round(float(line),10)
                        data.append(dataf)
                        
            except FileNotFoundError:
                print ("fe_prof file not found in current directory.")
            
        newfed = []
        try:
            
            with open ("FEDDAT.000_001",'r') as fh:
                
                for line in fh:
                    if "#" not in line:
                        break
                    newfed.append(line)
                datas = line.split()
                s = "{:>20}{:>20}{:>20}\n".format(datas[0],"%10f"%data[0],datas[2])
                newfed.append(s)
                
                for i in range(len(data)-1):
                    s = fh.readline()
                    datas = s.split()
                    s = "{:>20}{:>20}{:>20}\n".format(datas[0],"%10f"%data[i+1],datas[2])
                    newfed.append(s)
        
            with open ("FEDDAT.000",'w') as fw:
                for i in newfed:
                    fw.write(i)
        
        except FileNotFoundError:
            print ("FEDDAT.000_001 file not found in current directory.")
            
    def bias_to_prod_control(self, mode = "ee"):
        cont = []

        with open ("CONTROL",'r') as fh:
            for i in fh:
                if "fed method" in i:
                    if mode == "ee":
                        i = "fed method ee 1.0 1.0 1000000000\n"
                    elif mode == "tm_update":
                        i = "fed method tm 10000000 100\n"
                    else:
                        raise ValueError("Wrong mode input.")
                
                if "fed order param" in i:
                    spl = i.split()
                    spl = spl[:7]
                    spl.append("-1")
                    i = "{} {} {} {} {} {} {} {}\n".format(*spl)
                cont.append(i)
        with open ("CONTROL",'w') as fw:
            for i in cont:
                fw.write(i)

    
    def run_production (self,lbound,ubound):

        ndir = self.mc.get_new_run()
        shutil.copy(os.path.join(os.curdir,"FEDDAT.000_002"),os.path.join(ndir,"FEDDAT.000"))
        os.chdir(ndir)
        self.mc.create_psmc_control(lbound,ubound,biasf=-1,freq=500)
        
            
    def slab2bulk (self, bulk_vec, file = "CONFIG"):
        """
        Turn slab CONFIG to bulk CONFIG.

        Parameters
        ----------
        bulk_vec : ARRAY
            Stacking direction cell vector of the slab.
        file : STR, optional
            Name of file to edit. The default is "CONFIG".

        Returns
        -------
        None.

        """
        
        print("Looking for bulk length from metadise...")
        cont = []
        with open (file, 'r') as fh:
            for i in fh:
                cont.append(i)

        vecstr = "{:>19} {:>19} {:>19}\n".format("%.10f"%bulk_vec[0],"%.10f"%bulk_vec[1],"%.10f"%bulk_vec[2])
        cont[2] = vecstr
        
        with open(file, 'w') as fw:
            for i in cont:
                fw.write(i)
                
    def stack_info (self):
        """
        Fetch stacking related info from METADISE input file (Staco only).

        Returns
        -------
        stack : INT LIST
            Z direction info, includes slice stack and vaccumn region.
        growth : INT LIST
            growth vector.
        molcount : TYPE
            Number of mols read off of the staco file, low mol count suggests error in staco prep, should turn to use reorder algorithm.

        """
        
        file = "input.txt"
        molcount = 0
        with open (file,'r') as fh:
            for i in fh:
                if "region" in i:
                    temp = i.split()
                    stack = temp[1:5]
                if "MOLECULE" in i:
                    molcount += 1
                if "grow" in i:
                    temp = i.split()
                    growth = temp[1:3]
                    return stack,growth,molcount
        print ("Lattice super cell info obtained.")
    
    def bns_config (self, vac_space = 20, file = "CONFIG_DLM"):

        """
        Prepares all three CONFIG files required for LSMC, CONFIG_MASTER, CONFIG.1 and CONFIG.2.
        Requires a single slab CONFIG file in DLMONTE format named CONFIG_DLM.

        Parameters
        ----------
        vac_space : FLOAT, optional
            Amount of vaccumn space to add to bulk. The default is 20.
        file : STR, optional
            File name of the source CONFIG. The default is "CONFIG_DLM".

        Returns
        -------
        None.
        """
        
        cont = []
        
        with open (file) as fh:
            for i in fh:
                cont.append(i)
        
        with open ("CONFIG.1",'w') as fw1:
            for i in cont:
                fw1.write(i)
        
        vec_len = cont[2].split()[0]
        vac_len = float(vec_len) + vac_space
        cont_copy = cont[:]
        cont[2] = cont[2].replace(vec_len,"%.10f"%(vac_len))
        
        with open ("CONFIG.2",'w') as fw2:
            for i in cont:
                fw2.write(i)
                
        with open ("CONFIG_MASTER",'w') as fw3:
            for i in cont_copy:
                fw3.write(i)
            for i in cont:
                fw3.write(i)
        
        print ("Generated slab and bulk CONFIG files")
    
    def get_opt_bin(self, temp, input_range, nw, gradco, opt_lj_per_bin = 2.5, max_bin = 300000):

        max_grad = (-gradco[0]*temp + gradco[1])/temp
        m_per_bin = opt_lj_per_bin / max_grad
        bin_n = (input_range / m_per_bin)
        
        if bin_n > max_bin:
            bin_n = max_bin

        mul = int(bin_n/(nw+1))
        bins = mul*(nw+1)

        return bins

        
    def smooth_fe_prof(self, mul = 6):

        data = np.loadtxt("fe_prof")
        grad = data[1:] - data[:-1]
        sec_grad = grad[1:] - grad[:-1]
        import matplotlib.pyplot as plt
        abs_sec = abs(sec_grad)
        slice = grad[110:150]
        #for i in range(len(sec_grad)-1):
        #    i += 1
        #    if abs(sec_grad[i]-sec_grad[i-1]) > np.average(abs_sec):
        #        if abs(sec_grad[i])!=max(abs_sec) and abs(sec_grad[i-1])!= max(abs_sec):
        #            sec_grad[i] = sec_grad[i-1]
        #grad = [grad[0]]
        #for i in sec_grad:
        #    grad.append(grad[-1]+i)
        plt.plot(slice)
        plt.show()

    def get_max_atm_disp(self, col):
        disps = []
        with open("OUTPUT.000",'r') as fh:
            for line in fh:
                if "displacement (Angstroms) for" in line:
                    linestr = line.split()
                    if col == 1:
                        val = float(linestr[-2])
                    else:
                        val = float(linestr[-1])
                    disps.append(val)
        
        return max(disps)

    def get_range(self, surf, cont_fac = [0.3,0.3], get_energy = False, get_max_disp = False, mode = "lowavrange"):
        """
        Fetches OP range from rangeseek simulations.
        Rangeseek must be completed before using this function.

        Parameters
        ----------
        surf : STR
            Surface Miller Index to operate on.
        cont_fac : FLOAT, optional
            Contigency room added onto both sides of the range. The default is 0.2.

        Returns
        -------
        cont_opone : FLOAT
            Lower bound OP.
        cont_optwo : FLOAT
            Upper bound OP.

        """
        if surf != "/.":
            os.chdir("surface_{}".format(surf))
            self.ms.direct_tar()
            
        dirlist = sorted(os.listdir())
        dirone = [x for x in dirlist if "rangeone" in x]
        dirtwo = [x for x in dirlist if "rangetwo" in x]
        
        os.chdir(dirone[-1])
        dh = self.mcd.last_yaml_ten_percent()
        opone = round(np.average(dh[:,2]))
        if get_energy:
            ener_one = np.average(dh[:,4])
            from math import sqrt
            se_one = np.std(dh[:,4])/sqrt(len(dh[:,4]))

        os.chdir("..")
        os.chdir(dirtwo[-1])
        dh = self.mcd.last_yaml_ten_percent()
        
        if mode == "lowavrange":
            optwo = np.average(sorted(dh[:,2])[0:int(len(dh[:,2])/5)])
        elif mode == "lowmaxrange":
            optwo = max(sorted(dh[:,2])[0:int(len(dh[:,2])/5)])
        else:
            optwo = np.average(dh[:,2])
            
        if get_energy:
            ener_two = np.average(dh[:,4])
            se_two = np.std(dh[:,4])/sqrt(len(dh[:,4]))
        if get_max_disp:
            max_disp = self.get_max_atm_disp(2)
        
        os.chdir("..")
        cont_opone = opone - abs(cont_fac[0]*opone)
        cont_optwo = optwo + abs(cont_fac[1]*optwo)
        
        with open("misc", 'a') as fa:
            fa.write("Actual lbound: {}\n".format(opone))
            fa.write("Recommended lbound: {}\n".format(cont_opone))
            fa.write("Actual ubound: {}\n".format(optwo))
            fa.write("Recommended ubound: {}\n".format(cont_optwo))
        
        if surf != "/.":
            os.chdir("..")
            
        if get_energy:
            if get_max_disp:
                return cont_opone,cont_optwo,ener_one,ener_two,se_one,se_two,max_disp
            else:
                return cont_opone,cont_optwo,ener_one,ener_two,se_one,se_two
        else:
            if get_max_disp:
                return cont_opone,cont_optwo,get_max_disp
            else:
                return cont_opone,cont_optwo

    def get_configs(self, lbound, ubound, nbin, np = 16, mode = "inter"):

        incre = (ubound - lbound)/nbin
        intervals = []
        mid_points = []
        mc = MonteCore()
        gu = GeneralUtil()
        runs = []

        for i in range(nbin):
            wlbound = lbound + (i+1)*incre
            wubound = lbound + (i+2)*incre
            intervals.append([wlbound,wubound])
            mid_points.append((wlbound+wlbound)/2)
        
        if mode == "inter":
            for i in intervals:
                if len(runs) >= np:
                    mc.check_runs_terminate(runs,60)
                    runs = []
                tdir = mc.get_new_run ("configs_i")
                os.chdir(tdir)
                runs.append(tdir)
                mc.edit_windowed_control(i[0],i[1],mode="wrange")
                gu.edit_anything("steps","100000","CONTROL")
                mc.run_dlm("bg")
                os.chdir("..")
            if len(runs) > 0:
                mc.check_runs_terminate(runs,60)

        elif mode == "mid":
            for i in mid_points:
                if len(runs) >= np:
                    mc.check_runs_terminate(runs,60)
                    runs = []
                tdir = self.mc.get_new_run ("configs_m")
                runs.append(tdir)
                os.chdir(tdir)
                lb = i- incre/100
                ub = i+ incre/100
                mc.edit_windowed_control(lb,ub,mode="wrange")
                gu.edit_anything("steps","100000","CONTROL")
                mc.run_dlm("bg")
                os.chdir("..")
            if len(runs) > 0:
                mc.check_runs_terminate(runs,60)
                
        else:
            print("mode not recognised, available mode is mid or inter")
    
    
    
    
    