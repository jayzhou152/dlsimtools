# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:04:52 2020

@author: Jay
"""
import numpy as np
import os
import shutil
import subprocess
import sys
import time
import json
import csv

from .InputConverter import InputConverter
from .GeneralUtil import GeneralUtil
from .MonteCon import MonteCon

class MonteCore():
    
    """
    
    """
    
    
    def __init__(self, dlm_com = "/work/e05/e05/jz662/dlmonte_new/DL_MONTE-2/bin/DLMONTE-SRL.X", dlm_com_par = "/work/e05/e05/jz662/dlmonte_par/DL_MONTE-2/bin/DLMONTE-PRL.X"):
        
        """
        CONSTRUCTOR
        
        Parameters
        ----------
        dlm_com : STRING, optional
            DL_Monte call command in user PC. The default is "DLMONTE-SRL_psmcv.x".

        Returns
        -------
        None.

        """
        self.dlm_com = dlm_com
        self.dlm_com_par = dlm_com_par

    
    def run_dlm (self, mode = "bg", n = 8):
        """
        Function to run DLMONTE in shell command/other platforms.
        
        Parameters
        ----------
        mode : STRING, optional
            Mode for DLM commands, bg > running in background, fg > running in foreground. The default is "bg".

        Returns
        -------
        int
            error.

        """

        
        if mode == "bg":
            subprocess.run(self.dlm_com + " &", shell = True)
        elif mode == "fg":
            subprocess.run (self.dlm_com)
        elif mode == "short":
            process = subprocess.Popen(self.dlm_com)
            while "YAMLDATA.000" not in os.listdir():
                time.sleep(5)
            time.sleep(10)
            while not self.check_initial():
                time.sleep(5)
            process.terminate()
            time.sleep(120)  
        elif mode == "mpi":
            if n <= 16:
                mpi_com = "mpirun --oversubscribe -n {} {} &".format(n,self.dlm_com_par)
            else:
                subprocess.Popen(["export OMP_NUM_THREADS=1"],shell=True)
                subprocess.Popen(["export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"],shell=True)
                mpi_com = "srun --nodes={} --ntasks={} --hint=nomultithread --distribution=block:block {} &".format(int(n/129)+1,n,self.dlm_com_par)
            process = subprocess.Popen(mpi_com.split())           
        else:
            print ("DLMONTE mode not recognised")
            return 0
    
    def check_initial (self):
        
        with open ("YAMLDATA.000",'r') as fh:
            for line in fh:
                if line == "-\n":
                    return True 
        return False
        
    
    def convert_pm(self, uc=True):
        """
        To convert DLMONTE input to DLPOLY inputs, requires all three DLPOLY input files.
        

        Parameters
        ----------
        uc : BOOLEAN, optional
            Conversion to unique charges format. The default is True.

        Returns
        -------
        int
            error.

        """
        
        print ("Converting from DLPOLY to DLMONTE input.")
        
        files = ["CONFIG","CONTROL","FIELD"]
        
        if any(x not in os.listdir() for x in files):
            print ("Files required for conversion missing")
            return 0
        
        ic = InputConverter()
        ic.write_config()
        if "FIELD_sol" in os.listdir() and "CONFIG_sol" in os.listdir():
            ic = InputConverter(cfile = "CONFIG_sol", ffile = "FIELD_sol")
            ic.write_field(uc=True,lr_label=False)
        else:
            ic = InputConverter()
            ic.write_field(uc=True,lr_label=False)
    
    def create_psmc_control (self, lbound, ubound, biasf=1, temp=0.1, pressure=0, steps=1000000, 
                             freq=0, init = 1, fed_method = "ee", win = [], nbin = 100, cartdisp=False,
                             moveatm = True, movemol = False, movevol = False):
        """
        Creates a psmc CONTROL file, enabling run automation in current directory without the need for different 
        versions of CONTROL FILES.

        Parameters
        ----------
        lbound : FLOAT/INT
            Lower bound of the psmc order parameter.
        ubound : FLOAT/INT
            Upper bound of the psmc order parameter.
        biasf : INT, optional
            The bias scaling factor. The default is 1.
        temp : FLOAT, optional
            System temperature. The default is 0.1.
        pressure : FLOAT, optional
            System pressure. The default is 0.
        steps : INT, optional
            Number of steps in the run. The default is 1000000.
        freq : INT, optional
            Switch frequency for psmc. The default is 0.
        init : INT, optional
            Phase to start in, can only be 1 or 2. The default is 1.
        fed_method : STR, optional
            Fed method type. The default is "ee".
        win : FLOAT LIST, optional
            Window size for psmc (if empty defaults to no windows). The default is [].
        nbin : INT, optional
            Number of total bins in the simulation. The default is 100.
        cartdisp : BOOLEAN, optional
            Enables conservation of cartisian displacement (alternative is scaled displacement). The default is False.
        moveatm : BOOLEAN, optional
            Enables automatic fetching of atom move set from FIELD file. The default is True.
        movemol : BOOLEAN, optional
            Enables construction of molecule move list (creates a default XYZ molecule move). The default is False.
        movevol : BOOLEAN, optional
            Enables a vplume move set,. The default is False.

        Returns
        -------
        int
            error.

        """
        title = "PSMC simulation\n"
        
        mis = "use ortho\n use clists\n"
        
        if fed_method == "ee":
            method_str = "fed method ee 1.0 1.0 1000000000"
        elif fed_method == "tm":
            method_str = "fed method tm 1000000 100"
        else:
            print ("Input fed method not recognised.")
            return 0
        
        order_str = "fed order psmc {} {} {} {}".format(nbin,lbound,ubound,biasf)
        if not win:
            pass
        else:
            if len (win) == 2 and win[1] > win[0]:
                order_str += " win {} {}".format(*win)
            else:
                print ("Check input fed order parameter window.")
                return 0
        
        if cartdisp:
            disp_str = "     switchcartdisp\n"
        else:
            disp_str = ""
        
        fed = "  use fed psmc\n\n      {} switchfreq {}\n      initactive {}\n\n    psmc done\n\n    {}\n    {}\n\n  fed done\n\nfinish\n\n".format(disp_str, freq, init, method_str, order_str)
        
        seed = "ranseed\n"
        
        steps = "steps {}\n".format(steps)
        
        temperature = "temperature {}\n".format(temp)
        pressure = "pressure {}\n".format(pressure)
        
        yamlfreq = "yamldata 3000\n"
        
        with open("CONTROL", 'w') as fw:
            
            fw.write(title)
            fw.write(mis)
            fw.write(fed)
            fw.write(seed)
            fw.write(steps)
            fw.write("ewald prec 1e-7\n")
            fw.write(temperature)
            fw.write(pressure)
            fw.write(yamlfreq)
            fw.write("sample coordinate  10000\n")
            fw.write("archiveformat dlpoly4\n")
            
            if moveatm == True:
                atmlist = self.get_atm_list(file="FIELD")
                moveatm_str = "move atoms {} {}\n".format(len(atmlist), 20)
                fw.write(moveatm_str)
                for i in atmlist:
                    fw.write(i[0]+"    "+i[1]+"\n")
                
                
            if movemol == True:
                fw.write("move molecules 1 100\n")
                mollist = "XYZ\n"
                fw.write(mollist)
                
            if movevol == True:
                movevol_str = "move volume orthoani 1\n"
                fw.write(movevol_str)
            
            fw.write("start")
    
    def get_atm_list (self, file = "FIELD"):
        
        """
        Fetch atom list from DLMONTE FIELD file.

        Parameters
        ----------
        file : STR, optional
            File name to fetch atom list from, file must be in DLMONTE file format. The default is "FIELD".

        Returns
        -------
        atmlist : STR LIST
            List of atoms.

        """
        
        atmlist = []
        
        with open (file, 'r') as fh:
            
            for line in fh:
                if "ATOMS" in line:
                    for i in np.arange(int(line.split()[1])):
                        temp = fh.readline()
                        atmstr = temp.split()[0:2]
                        atmlist.append(atmstr)
                    return atmlist

    def get_new_run (self, key="", mode="LSMC", get_cont = True, get_feddat = False):
    
        
        """
        Automated generation of separate and new folders for DLMONTE runs.

        Parameters
        ----------
        key : STR, optional
            Key name for . The default is "".
        mode : STR, optional
            Name for new run modes, different modes moves different files. The default is "LSMC".

        Returns
        -------
        INT/STR
            Int 0 if error, string name for new folder if creation is successful.

        """
        
        
        for i in range(10000):
            
            if i > 1000:
                print("you have too many simulations in current folder, consider cleaning up")
            
            if mode == "LSMC":
                ndir_name = "mc_{}{}"
            else:
                ndir_name = "mc_{}{}"
            
            ndir_name = ndir_name.format(key,str(i+1))
            
            if ndir_name not in os.listdir():
                os.mkdir(ndir_name)
                
                try:

                    if mode == "LSMC":
                        shutil.copy(os.path.join(os.curdir,"CONFIG"),os.path.join(ndir_name,"CONFIG"))
                        shutil.copy(os.path.join(os.curdir,"CONFIG.1"),os.path.join(ndir_name,"CONFIG.1"))
                        shutil.copy(os.path.join(os.curdir,"CONFIG.2"),os.path.join(ndir_name,"CONFIG.2"))
                        shutil.copy(os.path.join(os.curdir,"FIELD"),os.path.join(ndir_name,"FIELD"))
                    elif mode == "normal":
                        shutil.copy(os.path.join(os.curdir,"CONFIG"),os.path.join(ndir_name,"CONFIG"))
                        shutil.copy(os.path.join(os.curdir,"FIELD"),os.path.join(ndir_name,"FIELD"))
                    else:
                        return 0
                    
                    
                    if get_cont:
                        shutil.copy(os.path.join(os.curdir,"CONTROL"),os.path.join(ndir_name,"CONTROL"))
                    if get_feddat:
                        shutil.copy(os.path.join(os.curdir,"FEDDAT.000"),os.path.join(ndir_name,"FEDDAT.000"))
                except FileNotFoundError:
                    print("Missing file, please check run file is in current directory.")
                
                return ndir_name
            
        return 0
    
    def edit_windowed_control(self, winlb, winub, file = "CONTROL", mode="wrange", softedge = False):
        
        flag = 0
        cont = []
        try:
            with open (file, mode = 'r') as fh:
                for i in range(1000):
                    
                    line = fh.readline()
                    if "fed method" in line:
                        flag = 1
                        cont.append(line)
                    elif "fed order param" in line and flag == 1 and line.split()[3] != "done":
                        line = line.split()
                        new_line = ""
                        for words in line:
                            if words == "win":
                                break
                            new_line += words + " "
                        if mode == "wrange":
                            new_line += "win {} {}".format(winlb,winub)
                            if softedge:
                                new_line += " softedges"
                            new_line += "\n"
                        elif mode == "allrange":
                            new_list = new_line.split()
                            new_list[-3] = str(winlb)
                            new_list[-2] = str(winub)
                            new_line = " ".join(new_list)
                            if softedge:
                                new_line += " softedges"
                            new_line += "\n"
                        elif "bins" in mode:
                            bin_str = mode.split()
                            if len(bin_str) != 2:
                                raise ValueError("Bin string is not given correctly, unable to fetch bin value.")
                            print(new_line)
                            new_list = new_line.split()
                            new_list[4] = bin_str[1]
                            new_line = " ".join(new_list)
                            new_line += "\n"
                        else:
                            print ("Editing mode is not recognised.")
                            return 0
                        cont.append(new_line)
                    else:
                        cont.append(line)
                    
                    if "start" in line:
                        break
            
            if flag == 0:
                print("fed method directive is not found.")
                sys.exit()

            with open (file, mode = 'w') as fw:
                for x in cont:
                    fw.write(x)
                    
        except FileNotFoundError:
            print("Control file does not exist.")

    
    def tmmc (self, mmin, mmax, nbin, nwin, temp, pressure, mode = "loose", ls=True, cont = True, upper_softedge = False, lower_softedge = False, bubbleless = False, nwin_crit=127):
        """
        

        Parameters
        ----------
        mmin : FLOAT
            Lower bound of the order parameter range.
        mmax : FLOAT
            Upper bound of the order parameter range.

            Number of total bins in the simulation.
        nwin : INT
            Number of total windows in the simulation.
        temp : FLOAT
            System temperature.
        pressure : FLOAT
            System pressure.
        mode : STR, optional
            Mode for order parameter checks (strict disables the use of floating point order parameter range)
            For example, if the order parameter range is atom numbers floating point order parameter should be disabled.
            The default is "loose".
        
        Returns
        -------
        int
            error.

        """
        run_folders = []

        if bubbleless:
            con_h = MonteCon()
            com = np.copy(con_h.cents)
            com[:,2] = abs(com[:,2])
            com = com[com[:,2].argsort()]

            slab_cent = np.vstack([i for i in com if i[0]=="XYZ"])
            solv_cent = np.vstack([i for i in com if i[0]!="XYZ"])

            solvent_list = solv_cent[:,1]
            n_slab = len(slab_cent)

        if mmax < mmin:
            print ("mmax must be larger than mmin")
            return 0
        
        if mode == "strict":
            ph = (mmax - mmin)/nbin
            
            if not ph.is_integer():
                raise ValueError("strict mode checks requires fed param to be divisible by nwin*nbin.")
                return 0
        
#        if mode == "gcmc":
#            ph = (mmax - mmin)/nbin
#            
#            if ph != 1:
#                raise ValueError("gcmc checks requires Nmin-Nmax to equal to Nbin")
#                return 0

        if nwin >= nwin_crit:
            print("Invoking cross node functionality.")
            proc_str = subprocess.check_output(['scontrol','show','hostnames'],universal_newlines=True)
            nodelist = proc_str.split()
            with open('nodes_list.txt','w') as fw:
                for i in nodelist:
                    fw.write(i+'\n')
            
        width = (mmax-mmin)/(nwin+1)
        self.edit_windowed_control(0,0,mode="bins {}".format(nbin))
        for i in np.arange(nwin):
            winmin = round((mmin + i * width),5)
            winmax = round ((mmin + (i+2) * width), 5)
            
            
            if upper_softedge == True and i == nwin-1:
                winmax = 1E10
            if lower_softedge == True and i == nwin-1:
                winmin = -1E10
            if ls:
                fname = self.get_new_run(key="tmmc",get_cont = cont)
            else:
                fname = self.get_new_run(key="tmmc", mode="normal",get_cont = cont)
            os.chdir(fname)

            run_folders.append(fname)

            if bubbleless:
                con_h_s = MonteCon()
                mid_point = int((winmin+winmax)/2)
                n_solv = mid_point - n_slab - 1
                rm_list = solvent_list[n_solv:]

                for j in rm_list:
                    con_h_s.del_mol(j)
                    con_h_s.write_out(filename="CONFIG")

            if not cont:
                self.create_psmc_control(mmin, mmax, temp=temp, pressure=pressure, steps=50000000,freq=500, \
                                         fed_method = "tm",win = [winmin,winmax],nbin = nbin, cartdisp = False )
            else:
                if upper_softedge == False and lower_softedge == False:
                    self.edit_windowed_control(winmin,winmax)
                else:
                    if i < nwin-1:
                        self.edit_windowed_control(winmin,winmax)
                    else:
                        self.edit_windowed_control(winmin,winmax,softedge = True)
            
            print(self.get_bandgap())

            if nwin >= nwin_crit:
                node_name = nodelist[int(i/nwin_crit)]
                print(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                    '--exact', '--mem=1500M', self.dlm_com, '&'])
                subprocess.Popen(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                    '--exact', '--mem=1500M', self.dlm_com, '&'])
            else:
                self.run_dlm()
            print ("tmmc no.{} running..".format(i+1))
            os.chdir("..")
        
        return run_folders
    
    
    def change_cellmat(self, mat, file = "CONFIG", mode = "normal"):
        
        mat = np.array(mat)
        if mat.shape != (3,3):
            raise ValueError("Input matrix is not a 3 by 3 array.")
        
        cellmat_str = ["{:>14} {:>14} {:>14}\n".format("%.7f"%mat[0,0],"%.7f"%mat[0,1],"%.7f"%mat[0,2]),
                       "{:>14} {:>14} {:>14}\n".format("%.7f"%mat[1,0],"%.7f"%mat[1,1],"%.7f"%mat[1,2]),
                       "{:>14} {:>14} {:>14}\n".format("%.7f"%mat[2,0],"%.7f"%mat[2,1],"%.7f"%mat[2,2])]
        
        cont = []
        with open (file,'r') as fh:
            
            cont.append(fh.readline())
            cont.append(fh.readline())
            
            if mode == "surf":
                cont.append(fh.readline())
                for i in range(2):
                    fh.readline()
                    cont.append(cellmat_str[i+1])
            else:
                for i in range(3):
                    fh.readline()
                    cont.append(cellmat_str[i])
            
            for line in fh:
                cont.append(line)
                
        with open(file,'w') as fw:
            for i in cont:
                fw.write(i)
        
        return cont

    def get_bandgap(self):

        with open ('CONTROL','r') as fh:
            for line in fh:
                if "fed order" in line:
                    lines = line.split()
                    if "param" not in line:
                        lines.insert(2,"param")
                    if len(lines) < 11:
                        print("Windowed bands not detected, fetched whole range.")
                        return abs(float(lines[6]) - float(lines[5])),[float(lines[5]),float(lines[6])],[0,0]
                    else:
                        if "softedges" in line:
                            if float(lines[10]) == 1e10 and float(lines[9])!= -1e10:
                                return abs(float(lines[9]) - float(lines[6])), [float(lines[5]),float(lines[6])], [float(lines[9]),float(lines[6])]
                            elif float(lines[10]) != 1e10 and float(lines[9]) == -1e10:
                                return abs(float(lines[5]) - float(lines[10])), [float(lines[5]),float(lines[6])], [float(lines[5]),float(lines[10])]
                        else:
                            return abs(float(lines[9]) - float(lines[10])), [float(lines[5]),float(lines[6])], [float(lines[9]),float(lines[10])]
    


    def check_in_window(self, runs, interval=1, threshold = 0.9):
        
        from .MonteData import MonteData
        md = MonteData()
        same_length_count = 1
        prev_length = 1

        print("initiating process of moving sims into window ranges")

        for i in runs[:int(len(runs)*threshold)]:

            os.chdir(i)
            window = self.get_bandgap()[-1]
            
            with open ("clue",'w') as fw:
                fw.write("loop1")
            while True:
                time.sleep(1)

                try:
                    same_length_count, prev_length = self.check_alive(600,same_length_count,prev_length)
                except TypeError:
                    return False

                try: 
                    cont = []
                    with open("YAMLDATA.000",'r') as fh:
                        for line in fh:
                            cont.append(line)
                except OSError:
                    continue
                if len(cont) > 20:
                    break

            with open ("clue",'w') as fw:
                fw.write("loop2")

            same_length_count = 1
            prev_length = 1

            while True:
                time.sleep(interval)

                try:
                    same_length_count, prev_length = self.check_alive(600,same_length_count,prev_length)
                except TypeError:
                    return False

                try:    
                    last = md.last_yaml_data()
                except ValueError:
                    continue
                print("Window lower bound is {}, upper bound is {}, current at {}.".format(window[0],window[1],last[2]))
                if last[2] >= window[0] and last[2] <= window[1]:
                    print("range reached, exiting...")
                    break
            os.chdir("..")

        return True
                
    def check_alive(self, n_count , same_length_count, prev_length):
        
        while True:
            if "OUTPUT.000" in os.listdir():
                if os.stat("OUTPUT.000").st_size > 1000:
                    break
        with open ("OUTPUT.000", 'rb') as fh:
                
            fh.seek(-2,2)
            while fh.read(1) != b"\n":
                fh.seek(-2,1)
                
            last_line = fh.read().decode()

        if last_line == " normal exit \n" or "maximum number of molecules reached" in last_line:
            print("Sim ended with normal exit")
            return False
            
        if "error" in last_line and "WARNING" not in last_line:
            raise KeyError("Process failure detected in {}, terminating run loop...".format(os.getcwd()))
            
            
        with open("OUTPUT.000", 'r') as f:
            try:
                for i,l in enumerate(f):
                    pass
            except OSError:
                i = prev_length

            if i == prev_length:
                same_length_count += 1
            else:
                same_length_count = 1


        if same_length_count >= n_count:
            print("sim ended due to OUTPUT not updating")
            return False
        
        else:
            prev_length = i
            return same_length_count , prev_length


    def check_terminate(self, interval, threshold = 600, sleep = 5, check_fe = False, tol = .1):
        
        same_length_count = 1
        prev_length = 1
        time.sleep(sleep)
        from .MonteData import MonteData
        md = MonteData()
        ncount = int(10800/interval)

        if check_fe: 
            print("checking for free energy convergence.")
            curdir = os.getcwd()
            os.chdir("..")
            run_count = 0 
            for i in os.listdir():
                if "tmmc" in i:
                    try:
                        os.chdir(i)
                        run_count += 1
                        os.chdir("..")
                    except NotADirectoryError:
                        continue
            cur_dir_n = curdir.split("/")[-1]
            prefix = cur_dir_n.split("tmmc")[0]
            last_dir = curdir.replace(cur_dir_n, prefix + "tmmc" + str(run_count))
            os.chdir(last_dir)

            fe_series_data = {}

            count = 0

            while True:
                time.sleep(interval)
                try:
                    same_length_count, prev_length = self.check_alive(ncount,same_length_count,prev_length)
                except TypeError:
                    return True
                
                if "TMATRX.000" in os.listdir():
                    time.sleep(30)
                    os.chdir("..")
                    print("TMATRX FOUND IN LAST DIR, calculating initial fe and transitioning to window 1 based fe calculations.")
                    initial_fe_series = md.tmat_profile(run_count,folderkey="mc_tmmc",valwrite= True, make_back = True)
                    while isinstance(initial_fe_series,int):
                        initial_fe_series = md.tmat_profile(run_count,folderkey="mc_tmmc",valwrite= True, make_back = True)
                    
                    if "fe_history.csv" not in os.listdir():
                        fe_series_data = np.linspace(0,len(initial_fe_series)-1,len(initial_fe_series))
                    else:
                        with open("fe_history.csv") as fp:
                            reader = csv.reader(fp, delimiter=",", quotechar='"')
                            # next(reader, None)  # skip the headers
                            data_read = [row for row in reader]
                        fe_series_data = np.array(data_read,dtype=float).T
                    
                    fe_series_data = np.vstack([fe_series_data,-initial_fe_series])
                    
                    os.chdir(curdir)
                    break

            with open ("TMATRX.000",'r') as fh:
                initial_tmat = fh.readlines()

            while len(initial_tmat) < 10:
                with open ("TMATRX.000",'r') as fh:
                    initial_tmat = fh.readlines()
            
            same_length_count = 1
            prev_length = 1
            
            while True:

                time.sleep(interval)
                try:
                    same_length_count, prev_length = self.check_alive(ncount,same_length_count,prev_length)
                except TypeError:
                    return True


                with open ("TMATRX.000",'r') as fh:
                    current_tmat = fh.readlines()

                if current_tmat != initial_tmat and len(current_tmat) > 10:

                    initial_tmat = current_tmat

                    
                    os.chdir("..")
                    current_fe_series = md.tmat_profile(run_count,folderkey="mc_tmmc",valwrite= True, make_back = True)
                    while isinstance(current_fe_series,int):
                        current_fe_series = md.tmat_profile(run_count,folderkey="mc_tmmc",valwrite= True, make_back = True)
                    fe_series_data= np.vstack([fe_series_data,-current_fe_series])

                    np.savetxt("fe_history.csv",fe_series_data.T,delimiter=",")
                        
                    max_er = max(abs(fe_series_data[-1] - fe_series_data[-2]))
                    with open("err_prog","a") as fw:
                        fw.write("{}\n".format(max_er))

                    if  max_er < tol:
                        print("converged to max error tolerance {} KT, current max error = {}".format(tol,max_er))
                        return True

                    os.chdir(curdir)


        else:
            print("checking for normal process behaviour")
            while True:
                time.sleep(interval)
                try:
                    same_length_count, prev_length = self.check_alive(ncount,same_length_count,prev_length)
                except TypeError:
                    print(self.check_alive(ncount,same_length_count,prev_length))
                    return True
            
    
    def gcexclude (self,filename= "CONFIG", mol_name = "XYZ", mode = "centre"):

        a = MonteCon()
        if a.nummol == 0:
            exclusion_factor = 0
        else:
            thickness = a.get_slab_thickness(mol_name = mol_name, mode = mode)
            exclusion_factor = thickness/a.cell_mat[a.stackd,a.stackd]/2

        cont = []
        flag = 0
        with open ("CONTROL",'r') as fh:
            for line in fh:
                if "use gcexcludeslab" in line:
                     line = "use gcexcludeslab {} -{} {} \n".format(a.stackd+1, exclusion_factor, exclusion_factor)
                     flag = 1
                if "finish" in line and flag == 0:
                    cont.append("use gcexcludeslab {} -{} {} \n".format(a.stackd+1, exclusion_factor, exclusion_factor))
                
                if "move gcinsertmol" in line:
                    mol_ins_line =  fh.readline()
                    gcmc_chempot = np.log(float(mol_ins_line.split()[1]))
                    new_chempot = np.exp(gcmc_chempot + np.log(1-2*exclusion_factor))
                    mol_ins_line = mol_ins_line.replace(mol_ins_line.split()[1],str(new_chempot))
                    cont.append(line)
                    cont.append(mol_ins_line)
                    continue
                
                cont.append(line)


        with open ("CONTROL",'w') as fw:
            for i in cont:
                fw.write(i)


    def tm_backup(self):
        
        for i in os.listdir():
            if "tmmc" in list:
                try:
                    os.chdir(i)
                    shutil.copy("TMATRX.000","TMATRX_back.000")
                    os.chdir("..")
                except OSError:
                    continue

 
    def check_runs_terminate(self, runs, interval, threshold = 600, check_fe= False, tol = 0.1):

        #same_line_count = np.zeros(len(runs))
        #prev_line = ["" for i in range(len(runs))]

        count = 0

        for i in runs:
            os.chdir(i)

            if self.check_terminate(interval,threshold = threshold, check_fe= check_fe, tol = tol):
                count +=1

            os.chdir("..")

            if check_fe:
                if count == 1:
                    return True
            else:
                if count == len(runs):
                    return True              
        