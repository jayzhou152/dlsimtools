
import numpy as np
from .LSMC import LSMC
from .MonteCore import MonteCore
from . import SwitchBias
import time
import os
import shutil
from .GeneralOptimizer import GeneralOptimizer
from .GeneralUtil import GeneralUtil
from .HPCworker import HPCWorker
from .MonteData import MonteData
from math import sqrt
from .MonteCon import MonteCon
import natsort
import subprocess


class Controller():
    
    def __init__(self, dlm_exec="default" , dlm_exec_par ="default" , rs_steps = 20000000, tmmc_steps = 2000000000, psteps = 2000000000, rlims = [-4e7,4e7], lsmc_mpi = False, use_sb= "se", rel_steps = 100000, cont_fac = [0.5,2], bin_fac = 1, added_sb = 0, rmode = "normal", nodes=1, nw=16, gradco = [5.089E-5,0.6136], prodmode = "ee", gcexclude = True, slab_molname = "XYZ", thickness_mode="edge", check_m = False, check_fe=True, tol=0.1, get_max_disp = True, upper_softedge = True, lower_softedge = False, rigid = False, volrel = "no_rel", volrel_param = [0,0,0,0,0,0], pull_cont_dis=1, nmol_strict= False,setup_cp = 1E100, setup_steps=10000000, new_control = False, sched="default"):
        """_summary_

        Args:
            rs_steps (int, optional): _description_. Defaults to 20000000.
            tmmc_steps (int, optional): _description_. Defaults to 2000000000.
            psteps (int, optional): _description_. Defaults to 2000000000.
            rlims (list, optional): _description_. Defaults to [-4e7,4e7].
            lsmc_mpi (bool, optional): _description_. Defaults to False.
            use_sb (str, optional): _description_. Defaults to "se".
            rel_steps (int, optional): _description_. Defaults to 100000.
            cont_fac (list, optional): _description_. Defaults to [0.5,2].
            bin_fac (int, optional): _description_. Defaults to 1.
            added_sb (int, optional): _description_. Defaults to 0.
            rmode (str, optional): _description_. Defaults to "normal".
            nodes (int, optional): _description_. Defaults to 1.
            nw (int, optional): _description_. Defaults to 16.
            gradco (list, optional): _description_. Defaults to [5.089E-5,0.6136].
            prodmode (str, optional): _description_. Defaults to "ee".
            gcexclude (bool, optional): _description_. Defaults to True.
            slab_molname (str, optional): _description_. Defaults to "XYZ".
            thickness_mode (str, optional): _description_. Defaults to "edge".
            check_m (bool, optional): _description_. Defaults to False.
            check_fe (bool, optional): _description_. Defaults to True.
            tol (float, optional): _description_. Defaults to 0.1.
            get_max_disp (bool, optional): _description_. Defaults to True.
            upper_softedge (bool, optional): _description_. Defaults to True.
            lower_softedge (bool, optional): _description_. Defaults to False.
            rigid (bool, optional): _description_. Defaults to False.
            volrel (str, optional): _description_. Defaults to "no_rel".
            volrel_param (list, optional): _description_. Defaults to [0,0,0,0,0,0].
            pull_cont_dis (int, optional): _description_. Defaults to 1.

        Raises:
            ValueError: _description_
        """
        self.dlm_exec = dlm_exec
        self.dlm_exec_par = dlm_exec_par
        self.rs_steps = rs_steps
        self.tmmc_steps = tmmc_steps
        self.psteps = psteps
        self.rlims = rlims
        self.lsmc_mpi = lsmc_mpi
        self.use_sb = use_sb
        self.rel_steps = rel_steps
        self.cont_fac = cont_fac
        self.bin_fac = bin_fac
        self.added_sb = added_sb
        self.rmode = rmode
        self.nodes = nodes
        self.nw = nw
        self.gradco = gradco
        self.prodmode = prodmode
        self.gcexclude = gcexclude
        self.slab_molname = slab_molname
        self.thickness_mode = thickness_mode
        self.check_m = check_m
        self.check_fe = check_fe
        self.tol = tol
        self.get_max_disp = get_max_disp
        self.upper_softedge = upper_softedge
        self.lower_softedge = lower_softedge
        self.rigid = rigid
        self.volrel = volrel
        self.volrel_param = volrel_param
        self.pull_cont_dis = pull_cont_dis
        self.nmol_strict = nmol_strict
        self.setup_cp = setup_cp
        self.setup_steps = setup_steps
        self.new_control = new_control
        self.sched = sched
        self.nwin_crit = 143 if sched == "isambard3" else 127
        if dlm_exec == "default" and dlm_exec_par == "default":
            self.mc = MonteCore()
        else:
            if dlm_exec == "default":
                self.mc = MonteCore(dlm_com_par = dlm_exec_par)
            elif dlm_exec_par == "default" :
                self.mc = MonteCore(dlm_com = dlm_exec)
            else:
                self.mc = MonteCore(dlm_com = dlm_exec, dlm_com_par = dlm_exec_par)
        
        if self.check_m and self.check_fe:
            raise ValueError("check_m and check_fe are mutually exclusive runtime checkers, please only enable one of them.")

    def master_surf_tmmc(self,temp):

        ls = LSMC(dlm_com=None if self.dlm_exec == "default" else self.dlm_exec,
                  dlm_com_par=None if self.dlm_exec_par == "default" else self.dlm_exec_par)
        mc = self.mc
        gu = GeneralUtil()

        instance = GeneralOptimizer(SwitchBias.get_phase,SwitchBias.bias_change,SwitchBias.run_sbtest,SwitchBias.range_check,SwitchBias.bias_mute)
        
        os.chdir(mc.get_new_run("temp{}_".format(temp)))
        #os.chdir("lsmc_temp{}_1".format(temp))
        gu.edit_anything("temperature",temp,"CONTROL")
        gu.edit_anything("steps",self.tmmc_steps,"CONTROL")
        if self.volrel == "sim_rel":
            print("Volume relaxation is caried out using simulations.")
            bulkz,bulkx,bulky,bulkarea,bulkvol = ls.npt_relax(temp = temp,pressure = 1, mpi = self.lsmc_mpi, n = self.nw, movevol = True, steps = self.rel_steps)
            with open("vol_info","a") as fw:
                fw.write(str(temp) + " {} {} {} {} {}\n".format(bulkz,bulkx,bulky,bulkarea,bulkvol))
        elif self.volrel == "no_rel":
            print("Volume relaxation is disallowed")
            a = MonteCon("CONFIG.1")
            with open("vol_info","a") as fw:
                fw.write(str(temp) + " {} {} {} {} {}\n".format(a.cell_mat[0,0],a.cell_mat[1,1],a.cell_mat[2,2],a.cell_mat[1,1]*a.cell_mat[2,2],a.cell_mat[0,0]*a.cell_mat[1,1]*a.cell_mat[2,2]))
        elif self.volrel == "param_rel":

            if self.volrel_param == [0,0,0,0,0,0]:
                raise ValueError("Volume relaxation using parameters enabled, however parameters not given in volrel_param. Exiting program.")

            print("Volume relaxtion is carried out using input functional parameters")
            z = self.volrel_param[0] * temp + self.volrel_param[1]
            x= self.volrel_param[2] * temp + self.volrel_param[3]
            y= self.volrel_param[4] * temp + self.volrel_param[5]

            bulkmat = [[z,0,0],[0,x,0],[0,0,y]]
            confonecont = mc.change_cellmat(bulkmat,file="CONFIG.1")
            conftwocont = mc.change_cellmat(bulkmat,file="CONFIG.2",mode="surf")

            with open ('CONFIG',"w") as fw:
                for i in confonecont:
                    fw.write(i)
                fw.write('\n')
                for i in conftwocont:
                    fw.write(i)

            with open("vol_info","a") as fw:
                fw.write(str(temp) + " {} {} {} {} {}".format(z,x,y,x*y,x*y*z))
        else:
            raise ValueError("Volume relaxation mode not recognised")

        rsnames = ls.rangeseek(temperature = temp, steps = self.rs_steps, mpi = self.lsmc_mpi, n = self.nw, rigid = self.rigid)
        time.sleep(2)
        mc.check_runs_terminate(rsnames, 10)

        if self.get_max_disp:
            rone, rtwo, eone, etwo, ph1, ph2, md = ls.get_range("/.",cont_fac = self.cont_fac,get_energy=True, get_max_disp = self.get_max_disp, mode = self.rmode)
            gu.edit_anything ("maxatmdist", md,"CONTROL")        
        else:
            rone, rtwo, eone, etwo, ph1, ph2 = ls.get_range("/.",cont_fac = self.cont_fac,get_energy=True, get_max_disp = self.get_max_disp, mode = self.rmode)
        rone = max(self.rlims[0],rone)
        rtwo = min(self.rlims[1],rtwo)
        with open ("range", 'w') as fw:
            fw.write("{}\n".format(rone))
            fw.write("{}".format(rtwo))
        mc.edit_windowed_control(round(rone,1), round(rtwo,1), mode = "allrange")
        bins =ls.get_opt_bin(temp, abs(rtwo-rone), self.nw, self.gradco)
        bins =int(bins * self.bin_fac)
        mc.edit_windowed_control(0,0,mode = "bins {}".format(bins))
        sb = (gu.ev_to_K(max([eone,etwo]) - min([eone,etwo]))/temp)
        se = abs(eone-etwo)
        with open("se",'w') as fw:
            fw.write(str(se) + "     ")
            fw.write(str(sqrt(ph1**2+ph2**2)))
            fw.write("    " + str(eone)) 
            fw.write("    " + str(etwo))
        sb_range = [sb*-5,sb*5]
        if self.use_sb == "optimised":
            result = instance.optimize(sb_range, tol = 0.00001)
            mid_p = (result[0] + result[1])/2 - self.added_sb
            SwitchBias.bias_change(mid_p)
        elif self.use_sb == "se":
            SwitchBias.bias_change(sb)
        elif self.use_sb == "zero":
            SwitchBias.bias_change(0)
        else:
            raise ValueError("sb mode must be 'zero', 'se' or 'optimised'.")
        tmruns = mc.tmmc (rone, rtwo, bins, self.nw, temp, 1, mode = "loose", cont = True, upper_softedge=self.upper_softedge, lower_softedge=self.lower_softedge, nwin_crit=self.nwin_crit)
        mc.check_runs_terminate(tmruns, 60 , check_m = self.check_m,check_fe=self.check_fe,tol=self.tol)
        os.chdir("..")

    def tmmc(self,temp,oprange):
        mc = self.mc
        os.chdir("lsmc_temp{}_1".format(temp))
        if oprange =="auto":
            with open ("range",'r') as fh:
                rone = float(fh.readline())
                rtwo = float(fh.readline())
        oprange = [rone,rtwo]
        tmruns = mc.tmmc(oprange[0],oprange[1],0,self.nw,temp,1,cont=True, upper_softedge=self.upper_softedge, lower_softedge=self.lower_softedge, nwin_crit=self.nwin_crit)
        time.sleep(120)
        mc.check_runs_terminate(tmruns,60,check_m=self.check_m,check_fe=self.check_fe)
        os.chdir("..")
    
    def write_tmmc(self, temp, process = "bias", contdir = "lsmc",bubbleless= False):
        
        with open("tmmc_master{}.py".format(temp),'w') as fw:

            fw.write("from dlsimtools.tmmc_master_control import Controller\n")

            f_str = "pc = Controller(dlm_exec = '{}', dlm_exec_par ='{}', rs_steps = {},tmmc_steps = {}," +\
                "psteps = {}, rlims = [{},{}], lsmc_mpi = {}, use_sb ='{}', rel_steps = {}, cont_fac = [{},{}], bin_fac = {}, added_sb = {}, rmode = '{}', nodes = {}, nw = {}, gradco = {}," +\
                "prodmode = '{}', gcexclude = {}, slab_molname = '{}', thickness_mode='{}', check_m = {}, check_fe = {}, tol = {}, get_max_disp = {}, upper_softedge = {}, lower_softedge = {}, rigid = {}, volrel = '{}', volrel_param=[{},{},{},{},{},{}], pull_cont_dis = {}, nmol_strict = {}, setup_cp = {}, setup_steps = {}, new_control = {}, sched = '{}')\n"
            f_str = f_str.format(self.dlm_exec, self.dlm_exec_par, self.rs_steps, self.tmmc_steps, self.psteps, self.rlims[0], self.rlims[1], self.lsmc_mpi, self.use_sb, self.rel_steps, self.cont_fac[0],
            self.cont_fac[1], self.bin_fac, self.added_sb, self.rmode, self.nodes, self.nw, self.gradco,
            self.prodmode, self.gcexclude, self.slab_molname, self.thickness_mode, self.check_m, self.check_fe, self.tol, self.get_max_disp, self.upper_softedge, self.lower_softedge, self.rigid, self.volrel, self.volrel_param[0],self.volrel_param[1],self.volrel_param[2],self.volrel_param[3],self.volrel_param[4],self.volrel_param[5],self.pull_cont_dis,self.nmol_strict,self.setup_cp,self.setup_steps,self.new_control,self.sched)

            fw.write (f_str)

            if "bias" in process:
                fw.write("pc.master_surf_tmmc({})".format(temp))
            elif "prod" in process:
                fw.write("pc.run_prods({})".format(temp))
            elif "tmmc_only" in process:
                fw.write("pc.tmmc({},'auto')".format(temp))
            elif "tmmc_cont" in process:
                dir_list = []
                for i in os.listdir():
                    if "{}_temp{}".format(contdir,temp) in i:
                        dir_list.append(int(i.split("_")[-1]))
                fw.write("pc.tmmc_continue(tdir='{}_temp{}_{}')".format(contdir,temp,max(dir_list)))
            elif "tmmc_wipe_cont" in process:
                dir_list = []
                for i in os.listdir():
                    if "{}_temp{}".format(contdir,temp) in i:
                        dir_list.append(int(i.split("_")[-1]))
                fw.write("pc.tmmc_continue(tdir='{}_temp{}_{}',wipe_prog=True)".format(contdir,temp,max(dir_list)))
            elif "tmmc_pull" in process:
                fw.write("pc.tmmc_distance({})".format(temp))
            elif "tmmc_pcont" in process:
                fw.write("pc.tmmc_continue(tdir='mc_dis_temp{}_1')".format(temp))
            elif "tmmc_nmols" in process:
                fw.write("pc.tmmc_nmols({},bubbleless={})".format(temp,bubbleless))
            elif "tmmc_nspec" in process:
                fw.write("pc.tmmc_nmols({},nspec=True)".format(temp))
            else:
                raise ValueError ("Wrong process chosen, available processes \
                    tmmc_only, tmmc_cont, tmmc_pull, tmmc_pcont, tmmc_nmols.")
    
    def run_folder(self, folder, timec=24, nodes=1, mem=0, prt="standard", qos="standard",
                   Qtype="premium", pcode=None, env=None):
        """Submit a single HPC job that runs all DL_MONTE simulations in `folder` sequentially.

        Each immediate subdirectory of `folder` that contains a CONTROL file is treated as one
        simulation. The simulations are executed one after another inside the job using the serial
        DL_MONTE executable (run_dlm mode='fg').

        Args:
            folder (str): Path to the parent directory whose subdirectories are simulations.
            timec (int): Wall-time in hours. Defaults to 24.
            nodes (int): Number of HPC nodes. Defaults to 1.
            mem (int): Memory in GB (0 = scheduler default). Defaults to 0.
            prt (str): Partition name. Defaults to 'standard'.
            qos (str): QoS string (archer2 only). Defaults to 'standard'.
            Qtype (str): Queue type ('premium' or 'free'). Defaults to 'premium'.
            pcode (str): Premium account code. Defaults to None.
            env (str): Conda/virtual environment name. Defaults to scheduler-appropriate default.

        Returns:
            str: Submitted job ID.
        """
        if self.sched not in ("balena", "archer2", "isambard3"):
            raise ValueError(
                "run_folder requires an HPC scheduler. "
                "Set sched='archer2', 'isambard3', or 'balena' when creating Controller."
            )

        if env is None:
            env = "myenv" if self.sched == "isambard3" else "py_env"

        folder = os.path.abspath(folder)
        sims = natsort.natsorted([
            d for d in os.listdir(folder)
            if os.path.isdir(os.path.join(folder, d))
            and "CONTROL" in os.listdir(os.path.join(folder, d))
        ])

        if not sims:
            raise ValueError(
                "No simulation directories (containing CONTROL) found in '{}'.".format(folder)
            )

        print("Found {} simulation(s) in '{}': {}".format(len(sims), folder, sims))

        cwdir = os.getcwd()
        os.chdir(folder)

        runner_name = "folder_runner.py"
        with open(runner_name, "w") as fw:
            fw.write("import os\n")
            fw.write("from dlsimtools.MonteCore import MonteCore\n\n")
            fw.write("mc = MonteCore(dlm_com='{}', dlm_com_par='{}')\n".format(
                self.mc.dlm_com, self.mc.dlm_com_par))
            fw.write("base = os.path.dirname(os.path.abspath(__file__))\n")
            fw.write("sims = {}\n\n".format(repr(sims)))
            fw.write("for sim in sims:\n")
            fw.write("    print('Running: ' + sim, flush=True)\n")
            fw.write("    os.chdir(os.path.join(base, sim))\n")
            fw.write("    mc.run_dlm(mode='fg')\n")
            fw.write("    print('Completed: ' + sim, flush=True)\n\n")
            fw.write("print('All simulations complete.')\n")

        exe_name = "python" if self.sched == "archer2" else "python3"
        exe = "{} {}".format(exe_name, runner_name)

        if self.sched == "balena":
            hpw = HPCWorker("slurm")
            sname = hpw.write_jobscript("slurm", env=env, exe=exe)
            runcom = hpw.get_runcom(nodes, timec, sname, prt=prt, mem=mem, qos=qos, Qtype=Qtype)
        elif self.sched == "archer2":
            if pcode is None:
                pcode = "e05-surfin-par"
            hpw = HPCWorker("archer2")
            sname = hpw.write_jobscript("archer2", env=env, exe=exe)
            runcom = hpw.get_runcom(nodes, timec, sname, prt=prt, mem=mem, qos=qos, Qtype=Qtype, premiumcode=pcode)
        elif self.sched == "isambard3":
            hpw = HPCWorker("isambard3")
            sname = hpw.write_jobscript("isambard3", env=env, exe=exe)
            runcom = hpw.get_runcom(nodes, timec, sname, prt=prt, mem=mem, Qtype=Qtype, premiumcode=pcode)

        jobid = hpw.submit_job(runcom)
        print("Submitted job {} for folder '{}'.".format(jobid, folder))

        os.chdir(cwdir)
        return jobid

    def tmmc_looper(self, temps, timec = 96, mem = 0, prt = "standard", qos ="standard", contdir = "mc_nmol", bubbleless=False, Qtype = "premium", pcode = None, process = "bias", bulk_image=False, env=None):

        if env is None:
            env = "myenv" if self.sched == "isambard3" else "py_env"

        print ("worker {}".format(self.sched))


        if bulk_image:
            #do some checks
            if all((process != "tmmc_nmols",process!= "tmmc_cont",process!="tmmc_wipe_cont")):
                raise RuntimeError("Process must be nmols, tmmc_cont or tmmc_wipe_cont to invoke bulk_image mode.")
            print ("invoking bulk image mode, creating addtional gcmc sim sets for the bulk solvent image.")

            mc = self.mc

        if self.sched == "balena":
            hpw = HPCWorker("slurm")

            for temp in temps:
                self.write_tmmc(temp, process=process, contdir=contdir, bubbleless = bubbleless)
                time.sleep(1)
                sname = hpw.write_jobscript("slurm", env="work", exe = "python3 tmmc_master{}.py".format(temp))
                runcom = hpw.get_runcom(self.nodes, timec, sname, prt = prt, mem = mem, qos=qos, Qtype=Qtype)
                jobid = hpw.submit_job(runcom)
                time.sleep(1)
                print ("Tmmc simulation for temperature {}K is running under job id {}".format(temp,jobid))

        elif self.sched == "archer2":
            if pcode is None:
                pcode = "e05-surfin-par"
            hpw = HPCWorker("archer2")

            for temp in temps:        
                if bulk_image:
                    if "cont" not in process:
                        dirone, dirbulk, dirslab = self.bulk_image_setup(temp)
                        os.chdir(dirone)
                        os.chdir(dirslab)
                    else:
                        dir_list = []
                        for i in os.listdir():
                            if "mc_gcmc_duo_temp{}".format(temp) in i:
                                dir_list.append(int(i.split("_")[-1]))
                        os.chdir("mc_gcmc_duo_temp{}_{}".format(temp,max(dir_list)))
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")
                            shutil.copy("../CONTROL_new","mc_slab1/CONTROL_new")
                        os.chdir("mc_slab1")
                
                self.write_tmmc(temp, process=process, contdir=contdir, bubbleless = bubbleless)
                time.sleep(1)
                sname = hpw.write_jobscript("archer2", env=env, exe = "python tmmc_master{}.py".format(temp))
                runcom = hpw.get_runcom(self.nodes, timec, sname, prt = prt, mem = mem, qos=qos, Qtype=Qtype, premiumcode=pcode)
                jobid = hpw.submit_job(runcom)
                time.sleep(1)
                print ("Tmmc simulation for temperature {}K is running under job id {}".format(temp,jobid))
                if bulk_image:
                    os.chdir("..")
                    self.gcexclude=False
                    if "cont" not in process:
                        os.chdir(dirbulk)
                    else:
                        os.chdir("mc_bulk_image1")
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")

                    self.write_tmmc(temp, process=process, contdir=contdir, bubbleless = bubbleless)
                    time.sleep(1)
                    sname = hpw.write_jobscript("archer2", env=env, exe = "python tmmc_master{}.py".format(temp))
                    runcom = hpw.get_runcom(self.nodes, timec, sname, prt = prt, mem = mem, Qtype=Qtype, premiumcode=pcode)
                    jobid = hpw.submit_job(runcom)
                    print ("Tmmc bulk image simulation for temperature {}K is running under job id {}".format(temp,jobid))
                    time.sleep(1)

        elif self.sched == "normal":
            if len (temps) > 20:
                raise ValueError("Too many processes for normal mode to handle at once.")
                return 0
            for temp in temps:

                if bulk_image:
                    if "cont" not in process:
                        dirone, dirbulk, dirslab = self.bulk_image_setup(temp)
                        os.chdir(dirone)
                        os.chdir(dirslab)
                    else:
                        dir_list = []
                        for i in os.listdir():
                            if "mc_gcmc_duo_temp{}".format(temp) in i:
                                dir_list.append(int(i.split("_")[-1]))
                        os.chdir("mc_gcmc_duo_temp{}_{}".format(temp,max(dir_list)))
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")
                            shutil.copy("../CONTROL_new","mc_slab1/CONTROL_new")
                        os.chdir("mc_slab1")

                if "bias" in process:
                    self.master_surf_tmmc(temp)
                elif "prod" in process:
                    self.run_prods(temp)
                elif "tmmc_only" in process:
                    self.tmmc(temp,"auto")
                elif "tmmc_cont" in process:
                    dir_list = []
                    for i in os.listdir():
                        if "{}_temp{}".format(contdir,temp) in i:
                            dir_list.append(int(i.split("_")[-1]))
                    self.tmmc_continue(tdir="{}_temp{}_{}".format(contdir,temp,max(dir_list)))
                elif "tmmc_wipe_cont" in process:
                    dir_list = []
                    for i in os.listdir():
                        if "{}_temp{}".format(contdir,temp) in i:
                            dir_list.append(int(i.split("_")[-1]))
                    self.tmmc_continue(tdir="{}_temp{}_{}".format(contdir,temp,max(dir_list)),wipe_prog=True)
                elif "tmmc_pull" in process:
                    self.tmmc_distance(temp)
                elif "tmmc_pcont" in process:
                    self.tmmc_continue(tdir="mc_dis_temp{}_1".format(temp))
                elif "tmmc_nmols" in process:
                    self.tmmc_nmols(temp,bubbleless=bubbleless)
                elif "tmmc_nspec" in process:
                    self.tmmc_nmols(temp,nspec=True)
                else:
                    raise ValueError ("Wrong process chosen, choose either bias or production.")
                
                if bulk_image:
                    os.chdir("..")
                    self.gcexclude=False
                    if "cont" not in process:
                        os.chdir(dirbulk)
                    else:
                        os.chdir("mc_bulk_image1")
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")

                    if "bias" in process:
                        self.master_surf_tmmc(temp)
                    elif "prod" in process:
                        self.run_prods(temp)
                    elif "tmmc_only" in process:
                        self.tmmc(temp,"auto")
                    elif "tmmc_cont" in process:
                        dir_list = []
                        for i in os.listdir():
                            if "{}_temp{}".format(contdir,temp) in i:
                                dir_list.append(int(i.split("_")[-1]))
                        self.tmmc_continue(tdir="{}_temp{}_{}".format(contdir,temp,max(dir_list)))
                    elif "tmmc_wipe_cont" in process:
                        dir_list = []
                        for i in os.listdir():
                            if "{}_temp{}".format(contdir,temp) in i:
                                dir_list.append(int(i.split("_")[-1]))
                        self.tmmc_continue(tdir="{}_temp{}_{}".format(contdir,temp,max(dir_list)),wipe_prog=True)
                    elif "tmmc_pull" in process:
                        self.tmmc_distance(temp)
                    elif "tmmc_pcont" in process:
                        self.tmmc_continue(tdir="mc_dis_temp{}_1".format(temp))
                    elif "tmmc_nmols" in process:
                        self.tmmc_nmols(temp,bubbleless=bubbleless)
                    elif "tmmc_nspec" in process:
                        self.tmmc_nmols(temp,nspec=True)
                    else:
                        raise ValueError ("Wrong process chosen, choose either bias or production.")
                
        elif self.sched == "isambard3":
            hpw = HPCWorker("isambard3")

            for temp in temps:
                if bulk_image:
                    if "cont" not in process:
                        dirone, dirbulk, dirslab = self.bulk_image_setup(temp)
                        os.chdir(dirone)
                        os.chdir(dirslab)
                    else:
                        dir_list = []
                        for i in os.listdir():
                            if "mc_gcmc_duo_temp{}".format(temp) in i:
                                dir_list.append(int(i.split("_")[-1]))
                        os.chdir("mc_gcmc_duo_temp{}_{}".format(temp,max(dir_list)))
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")
                            shutil.copy("../CONTROL_new","mc_slab1/CONTROL_new")
                        os.chdir("mc_slab1")

                self.write_tmmc(temp, process=process, contdir=contdir, bubbleless=bubbleless)
                time.sleep(1)
                sname = hpw.write_jobscript("isambard3", env=env, exe="python tmmc_master{}.py".format(temp))
                runcom = hpw.get_runcom(self.nodes, timec, sname, prt=prt, mem=mem, Qtype=Qtype, premiumcode=pcode)
                jobid = hpw.submit_job(runcom)
                time.sleep(1)
                print("Tmmc simulation for temperature {}K is running under job id {}".format(temp, jobid))
                if bulk_image:
                    os.chdir("..")
                    self.gcexclude = False
                    if "cont" not in process:
                        os.chdir(dirbulk)
                    else:
                        os.chdir("mc_bulk_image1")
                        if self.new_control:
                            shutil.copy("../CONTROL_new","CONTROL_new")

                    self.write_tmmc(temp, process=process, contdir=contdir, bubbleless=bubbleless)
                    time.sleep(1)
                    sname = hpw.write_jobscript("isambard3", env=env, exe="python tmmc_master{}.py".format(temp))
                    runcom = hpw.get_runcom(self.nodes, timec, sname, prt=prt, mem=mem, Qtype=Qtype, premiumcode=pcode)
                    jobid = hpw.submit_job(runcom)
                    print("Tmmc bulk image simulation for temperature {}K is running under job id {}".format(temp, jobid))
                    time.sleep(1)

        else:
            raise ValueError("Inappropriate sched value on Controller. Use either 'balena', 'archer2', 'isambard3' or 'normal'.")
    
    def bulk_image_setup (self,temp):
        mc = self.mc
        dirone = mc.get_new_run(key="gcmc_duo_temp{}_".format(temp),mode="normal")
        os.chdir(dirone)
        dirbulk = mc.get_new_run(key = "bulk_image",mode="normal")
        dirslab = mc.get_new_run(key = "slab",mode="normal")

        os.chdir(dirbulk)
        con_cont = []
        a = MonteCon()
        thickness = a.get_slab_thickness(mol_name = self.slab_molname,mode = self.thickness_mode)
        exclusion_factor = thickness/a.cell_mat[a.stackd,a.stackd]/2
        a.nummol = 0
        a.cell_mat[a.stackd,a.stackd] =  a.cell_mat[a.stackd,a.stackd] * (1-2*exclusion_factor)
        a.write_out(filename="CONFIG")

        bg = mc.get_bandgap()
        lb = -0.5
        ub = lb + bg[0]
        mc.edit_windowed_control(lb,ub,mode="allrange")         
        cont = []
        with open ("CONTROL",'r') as fh:
            for line in fh:
                if "move atoms" in line:
                    n_read = int(line.split()[2])
                    for i in range(n_read):
                        fh.readline()
                elif "maxatmdist" in line:
                    continue
                elif "acceptatmmoveupdate" in line:
                    continue
                else:    
                    cont.append(line)
        with open("CONTROL", 'w') as fw:
            for i in cont:
                fw.write(i)

        os.chdir("..")            
        os.chdir("..")
        return dirone, dirbulk, dirslab


    def tmmc_get_fe_vals(self, update = True):
        
        gu = GeneralUtil()
        md = MonteData()
        sb = float(gu.get_anything("switchbias", "CONTROL"))
        dirs = os.listdir()
        
        count = 0
        for i in dirs:
            if "lsmc_tmmc" in i:
                count +=1
        if update:
            bias_vec = md.tmat_profile(count,folderkey="lsmc_tmmc",valwrite=True)
            fe = -bias_vec
        else:
            try:
                fe = np.loadtxt("fe_prof")
            except OSError:
                print("fe_prof not found, turn on update to produce fe_prof")
                return 0
        fe = fe[int(len(fe)/50):-int(len(fe)/50)]
        grad = fe[1:] - fe[:-1]
        secgrad = abs(grad[1:] - grad[:-1])
        midpoint = np.where(secgrad == max(secgrad))[0][0]
        
        bulk_fe = fe[:midpoint]
        surf_fe = fe[midpoint:]
        fed = min(bulk_fe) - min(surf_fe)
        print(sb,fed, sb + fed, midpoint)
        
        fed_corrected = sb + fed
        
        return -fed_corrected
    
    def produce_se(self, temps):

        ls = LSMC(dlm_com=None if self.dlm_exec == "default" else self.dlm_exec,
                  dlm_com_par=None if self.dlm_exec_par == "default" else self.dlm_exec_par)

        #if temps == "auto":
        #    temps = []
        #    dirs = os.listdir()
        #    for i in dirs:
        #        if "lsmc_temp" in i:
        #            x = i.split("lsmc_temp")[1]
        #            y = int(x[:-2])
        #            temps.append(y)
        #    temps = sorted(temps)

        dirs = os.listdir()
        tdirs=[]
        for i in dirs:
            if "mc_temp" in i:
                if temps == "auto":
                    tdirs.append(i)
                else:
                    x = i.split("mc_temp")[1]
                    y = int(x[:-2])
                    if y in temps:
                        tdirs.append(i)
        tdirs = sorted(tdirs)


        for temp in tdirs:
            print("Grabbing surface energy info for temp {}K.".format(temp))
            os.chdir(temp)
            rone,rtwo,eone,etwo,se_one,se_two = ls.get_range("/.",get_energy=True)
            from math import sqrt
            se = abs(eone-etwo)
            print(temp,eone,etwo,se)
            standard_error = sqrt(se_one**2 + se_two**2)
            with open("se",'w') as fw:
                fw.write(str(se))
                fw.write(str("    "))
                fw.write(str(standard_error))
                fw.write("    ")
                fw.write(str(eone))
                fw.write("    ")
                fw.write(str(etwo))            
            os.chdir("..")
        print("Surface energy info gathered and written to 'se' file in corresponding directories.")
    
    def get_fe_val_looper(self, temps, write = True, unit = "ev", update= True, get_dim = True):
        
        feds= []
        entro = []
        ses = []
        area = []
        xlength = []
        ylength = []
        volume = []
        stes = []
        slab_es = []
        bulk_es = []
        dirs = os.listdir()
        tdirs = []
        temps_ph = []
            #for i in dirs:
            #    if "lsmc_temp" in i:
            #        x = i.split("lsmc_temp")[1]
            #        y = int(x[:-2])
            #        temps.append(y)
            #temps = sorted(temps)

        for i in dirs:
            if "mc_temp" in i:
                if temps == "auto":
                    tdirs.append(i)
                else:
                    x = i.split("mc_temp")[1]
                    y = int(x[:-2])
                    if y in temps:
                        tdirs.append(i)
        tdirs = natsort.natsorted(tdirs)

        
        for tdir in tdirs:
            print(tdir)
            os.chdir(tdir)
            fed = self.tmmc_get_fe_vals(update = update)
            x = tdir.split("mc_temp")[1]
            temp = int(x[:-2])
            if unit != "lj":
                gu = GeneralUtil()
                fed = gu.lj_out(fed,temp,unit)
            feds.append(fed)
            with open("se",'r') as fh:
                dat_line = fh.readline()
                se = float(dat_line.split()[0])
                ste = float(dat_line.split()[1])
                bulk_e = float(dat_line.split()[2])
                slab_e = float(dat_line.split()[3])
            ses.append(se) 
            stes.append(ste)
            bulk_es.append(bulk_e)
            slab_es.append(slab_e)
            entro.append((se-fed)/temp)
            temps_ph.append(temp)
            with open("vol_info",'r') as fh:
                data = fh.readline()
                area_t = float(data.split()[2]) * float(data.split()[3]) *2
                xlength.append(float(data.split()[2]))
                ylength.append(float(data.split()[3]))
                volume.append(float(data.split()[5]))
                area.append(area_t)
            os.chdir("..")
        
        fe_pa = [x*16.02/y for x,y in zip(feds,area)]
        se_pa = [x*16.02/y for x,y in zip(ses,area)]
        se_se_pa = [x*16.02/y for x,y in zip(stes,area)]


        if write:
            with open("feds",'w') as fw:
                fw.write ("temp(K)/ FreeEnergy(eV) /SurfaceEnergy(EV) / Entropy(EV)/ X(A)/Y(A)/Area(A^2)/volume(A^3)/ FE (J/m^2) / SE (J/m^2)/ SE_STE (J/m^2)\n")
                for i,j,k,l,m,n,o,p,q,r,s,t,u in zip(temps_ph,feds,ses,bulk_es,slab_es,entro,xlength,ylength,area,volume,fe_pa,se_pa,se_se_pa):
                    fw.write(str(i) + "    " + str(j) +  "    " + str(k) + "    " + str(l) + "    " + str(m) + "    " + str(n) + "    " + str(o) + "    " + str(p) + "    " + str(q) + "    " + str(r) + "    " + str(s) + "    " + str(t) + "    " + str(u) + " \n")
                    
        return temps_ph,feds,entro

    def tmmc_continue(self, tdir="mc_nmols_temp300_1", mckey = "mc", wipe_prog= False, active=False):

        mc = self.mc
        md = MonteData()
        os.chdir(tdir)
        count = 0
        tmruns = []
        bin_n = md.get_bin_num()

        if active:
            jobid = os.environ['SLURM_JOB_ID']
            steps_str = subprocess.check_output(['sstat','-o','JOBID',jobid],universal_newlines=True)
            sstats = steps_str.split()[2:]

        if wipe_prog:
            try:
                shutil.move("fe_history.csv","fe_history_back.csv")
                print("Wipe_progress, fe_history is backed up, proceeding to wipe initial TMAT progress.")
            except FileNotFoundError:
                print("Wipe_progress, fe_history does not exist, however proceeding to wipe initial TMAT progress.")
                pass
        
        else:
            if self.new_control:
                raise ValueError("new_control keyword turned on when progress wipe is turned off, this functionality can only be used when progress wipe is turned on ")
        
        if self.new_control:
            try:
                shutil.copy("../CONTROL_new","CONTROL")
            except FileNotFoundError:
                raise FileNotFoundError ("CONTROL_new must exist in the parent directory when new_control is turned on.")

        if self.nw >= 127:
            print("Invoking cross node functionality.")
            proc_str = subprocess.check_output(['scontrol','show','hostnames'],universal_newlines=True)
            nodelist = proc_str.split()
            with open('nodes_list.txt','w') as fw:
                for i in nodelist:
                    fw.write(i+'\n')

        for i in os.listdir():

            if "{}_tmmc".format(mckey) in i:
                os.chdir(i)
                tmruns.append(i)
                sim_n = int(i.split("tmmc")[1]) -1
                if active:
                    subprocess.Popen(['scancel','--signal=TERM', sstats[count]])
                    count += 1
                    time.sleep(10)
                if "REVCON.000" in os.listdir():
                    shutil.move("REVCON.000","CONFIG")
                
                if not wipe_prog:

                    if "TMATRX.000" in os.listdir():
                        if len(md.get_tmat()) == bin_n:
                            shutil.copy("TMATRX.000","TMATRX")
                        else:
                            shutil.copy("TMATRX_back.000","TMATRX")

                    if "TMATRX" not in os.listdir():
                        raise FileNotFoundError("TMATRX does not exist in this folder.")

                    cont = []
                    with open ("CONTROL",'r') as fh:
                        for line in fh:
                            if "fed method" in line:
                                line = line[:-1] + " resume\n"
                            cont.append(line)
                    with open ("CONTROL",'w') as fw:
                        for j in cont:
                            fw.write(j)
                else:
                    if self.new_control:
                        try:
                            bd = mc.get_bandgap()[2]
                        except IndexError:
                            raise ValueError("Cannot find window range in folder {}.".format(i+1))
                        shutil.copy ("../CONTROL", "CONTROL")
                        mc.edit_windowed_control(bd[0],bd[1],mode="wrange")

                    if "TMATRX.000" in os.listdir():
                        os.remove("TMATRX.000")
                        
                if self.nw >= 127:
                    node_name = nodelist[int(sim_n/127)]
                    print(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                        '--exact', '--mem=1500M', mc.dlm_com, '&'])
                    subprocess.Popen(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                        '--exact', '--mem=1500M', mc.dlm_com, '&'])
                else:
                    mc.run_dlm()
                os.chdir("..")
        
        mc.check_runs_terminate(tmruns, 60 ,threshold = 6000,check_fe = self.check_fe, check_m = self.check_m, tol = self.tol)
        os.chdir("..")

    def tmmc_distance(self, temp):

        mc = self.mc
        gu = GeneralUtil()

        with open ("CONTROL",'r') as fh:
            for line in fh:
                if "fed order parameter com2" in line:
                    drangeone = float(line.split()[5])
                    drangetwo = float(line.split()[6])
                    bin = int(line.split()[4])
                
        with open ("CONTROL",'r') as fh:
            if "move gcinsertmol" in fh.read():
                mode = "solvent"
                print("solvent control file detected, invoking solvent molecule pulling mode")
            else:
                mode = "vacuum"
                print("solvent control file not detected, invoking vacuum molecule pulling mode")

        
        if mode == "solvent":
            with open ("densities",'r') as fh:
                sol_mol_weight = float(fh.readline())
                for line in fh:
                    if int(line.split()[0]) == temp:
                        density = float(line.split()[1])
                        break
            a = MonteCon()
            if a.nummol == 0:
                thickness = 0
            else:
                thickness = a.get_slab_thickness(mol_name = self.slab_molname, mode = self.thickness_mode)
                n_layers, layers_label, layers_dis = a.find_layers()
                layer_thick = thickness/n_layers
                full_nummol_per_layer = len(layers_label[int(len(layers_label)/2)])
                full_con_n = n_layers * full_nummol_per_layer
                missing_n = full_con_n - a.nummol
                missing_v = (a.cell_mat[1,1] * a.cell_mat[2,2] * layer_thick) * (missing_n/full_nummol_per_layer) * 1e-30
                volume = ((a.cell_mat[0,0] - thickness) * a.cell_mat[1,1] * a.cell_mat[2,2] )*1e-30 + missing_v
                mol_dens = density/(sol_mol_weight/1000)*6.02214e23
                target_n = int(volume * mol_dens)
            print(target_n)
            a.mol_max[-1] = target_n
            
            
        tdir = mc.get_new_run(key = "dis_temp{}_".format(temp), mode = "normal")
        os.chdir(tdir)
        if mode == "solvent":
            a.write_out(filename = "CONFIG")
        gu.edit_anything("temperature",temp,"CONTROL")
        gu.edit_anything("steps",0,"CONTROL")
        


        if bin <= self.nw + 1:
            bin = self.nw + 1
        else:
            bin = int((bin/(self.nw + 1))) * (self.nw+1)

        runs = mc.tmmc(drangeone,drangetwo,bin,self.nw,temp,0,ls=False,nwin_crit=self.nwin_crit)
        mc.check_runs_terminate(runs, 5, check_fe = False, check_m = False)

        #mc.check_in_window(runs)
        #subprocess.run("pkill DLMONTE-SRL.X",shell=True)
        #time.sleep(60)

 

        for i in runs:
            os.chdir(i)

            bdgap = mc.get_bandgap()
            print(bdgap)
            wind = bdgap[-1]
            tot_range = bdgap[-2]

            mid = (wind[0] + wind[1]) / 2 
            print(mid)
            dif = mid - tot_range[0]

            con = MonteCon(filename="CONFIG")
            
            for i in con.data:
                if i[0] == "TAR_MOL":
                    for j in i[2]:
                        j[1] = round(j[1] + dif - self.pull_cont_dis,5)
        
            con.write_out(filename = "CONFIG")

            os.chdir("..")

        if self.nw >= 127:
            print("Invoking cross node functionality.")
            proc_str = subprocess.check_output(['scontrol','show','hostnames'],universal_newlines=True)
            nodelist = proc_str.split()
            with open('nodes_list.txt','w') as fw:
                for i in nodelist:
                    fw.write(i+'\n')

        if mode == "solvent":
            for i in runs:
                cont = []
                os.chdir(i)

                with open ("CONTROL",'r') as fh:
                    for line in fh:
                        if "TAR_MOL" in line:
                            continue
                        elif "TP3O" in line:
                            line = "TP3O\n"
                            cont.append(line)
                        elif "move gcinsertmol" in line:
                            line = "move gcinsertmol {} {}\n".format(line.split()[2],"1E4")
                            cont.append(line)
                            x = fh.readline()
                            #chempot = float(x.split()[1])
                            x = x.replace(x.split()[1],"1E30") 
                            x = x.replace("#","")
                            cont.append(x)
                        else:
                            cont.append(line)
                
                with open ("CONTROL",'w') as fw:
                    for j in cont:
                        fw.write(j)
                
                
                #os.remove("CONFIG")
                #os.rename("REVCON.000","CONFIG")
                sim_n = int(i.split("tmmc")[1]) - 1
                gu.edit_anything("steps",1000000000,"CONTROL")

                if self.nw >= 127:
                    node_name = nodelist[int(sim_n/127)]
                    print(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                        '--exact', '--mem=1500M', mc.dlm_com, '&'])
                    subprocess.Popen(['srun','--nodelist={}'.format(node_name),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
                        '--exact', '--mem=1500M', mc.dlm_com, '&'])
                else:
                    mc.run_dlm()

                #mc.run_dlm()
                os.chdir("..")

            mc.check_runs_terminate(runs, 5, check_fe = False, check_m = False)

            #try:
            #    mc.check_runs_terminate(runs, 5, check_fe = False, check_m = False)
            #except KeyError:
            #    os.chdir("..")
            #    os.chdir(runs[0])
            #    b = MonteCon("REVCON.000")
            #    if b.nummol != target_n + a.nummol:
            #        raise ValueError("Molecule filling sims failed without reaching target number of molecules.")
            #    os.chdir("..")
            
            if self.nw >= 127:
                jobid = os.environ['SLURM_JOB_ID']
                subprocess.Popen(['scancel','--signal=TERM',jobid])
            else:
                subprocess.Popen(["pkill","DLMONTE-SRL.X"])

        for i in runs:
            
            os.chdir(i)

            if mode == "solvent":
                cont = []
                
                with open("CONTROL",'r') as fh:
                    
                    for line in fh:
                        #if "move gcinsertmol" in line:
                        #    x = fh.readline()
                        #    x = x.replace(x.split()[1],str(chempot))
                        #    line = "move gcinsertmol {} {}\n".format(line.split()[2],ins_freq)
                        #    cont.append(line)
                        #    cont.append(x)

                        if "move molecule" in line or "move rotatemol" in line:
                            line = line.split()
                            new_line = "{} {} {} {}\n".format(line[0],line[1],int(line[2])+1,line[3])
                            cont.append(new_line)
                            cont.append("TAR_MOL\n")
                        elif "move gcinsertmol" in line:
                            for i in range(int(line.split()[2])):
                                fh.readline()
                        else:
                            cont.append(line)

                with open ("CONTROL",'w') as fw:
                    for j in cont:
                        fw.write(j)
            
                os.remove("CONFIG")
                os.rename("REVCON.000","CONFIG")
            else:                         
                os.remove("REVCON.000")

            gu.edit_anything("steps",self.tmmc_steps,"CONTROL")
            
            os.chdir("..")            
        os.chdir("..")
        
        self.tmmc_continue(tdir = tdir, mckey="mc", wipe_prog=True)

            
    def removal_order(self, temp, indices, limit = False, ndir="find_order", mode = "relaxed"):
        
        if "order" in os.listdir():
            with open ("order",'r') as fh:
                o_str = fh.readline()
                rm_order = list(map(float,o_str.split()))
                return rm_order

        gu = GeneralUtil()  
        mcr = self.mc
        md = MonteData()
        shutil.rmtree(ndir)
        os.mkdir(ndir)
        shutil.copy ("CONFIG",ndir)
        shutil.copy ("CONTROL",ndir)
        shutil.copy ("FIELD",ndir)
        os.chdir(ndir)
        shutil.move("CONFIG","CONFIG.O")
        if mode == "singp":
            gu.edit_anything("steps","0","CONTROL")
        elif mode == "relaxed":
            gu.edit_anything("steps","100000","CONTROL")
        gu.edit_anything("temperature",str(temp),"CONTROL")
        rm_order = []
        

        if limit:
            stopper = int(len(indices)/4)*3
        else:
            stopper = 0

        while len(indices) > stopper:
            
            energies = []
            run_dirs = []

            for i in indices:
                mc = MonteCon(filename="CONFIG.O")
                for mols in rm_order:
                    mc.del_mol(mols)
                mc.del_mol(i)
                mc.write_out(filename = "CONFIG")
                de_dir = mcr.get_new_run(key="defeng",mode = "normal", get_cont=True)
                os.chdir(de_dir)
                mcr.run_dlm(mode="bg")
                os.chdir("..")
                run_dirs.append(de_dir)

            mcr.check_runs_terminate(run_dirs,5)

            for i in run_dirs:
                os.chdir(i)
                data = md.yaml_one(1)
                data = data[int(0.4*len(data)):]
                eng = np.average(data)
                energies.append(eng)
                os.chdir("..")
            
            print(energies)
            print(indices)

    
            low = energies.index(min(energies))
            ind = indices.pop(low)
            rm_order.append(ind)
            #mc = MonteCon(filename = "CONFIG.O")
            #mc.del_mol(ind)
            #mc.write_out(filename = "CONFIG.O")
        os.chdir("..")    
        
        with open ("order","w") as fw:
            for i in rm_order:
                fw.write(str(i))
                fw.write("  ")
        
        return rm_order

    def multiple_removal_set_up(self, temps, preemble, indices, pull_dis, tmmc_cont = False, **kwargs):
        
        from .MonteCon import MonteCon
        #checks for preemble inside indices.
        if any(np.in1d(preemble,indices)):
            raise ValueError("Preemble molecules are pre-removed molecules and is therefore not allowed to be removed.")

        layer = len(preemble) + 1

        for i in indices:
            if not tmmc_cont:
                a = MonteCon()
                for j in preemble:
                    a.del_mol(j)
                tar_i_n = i - sum(k < i for k in preemble)
                a.change_mollabel(tar_i_n,"TAR_MOL")

            ndir = "step_layer{}_mol{}".format(layer,i)    
            for j in preemble:
                ndir += "_" + str(j)

            if not tmmc_cont:
                try:

                    os.mkdir(ndir)
                    shutil.copy("CONTROL",ndir)
                    shutil.copy("FIELD", ndir)
                    shutil.copy("densities",ndir)
                except FileExistsError:

                    shutil.copy("CONTROL",ndir)
                    shutil.copy("FIELD", ndir)
                    shutil.copy("densities",ndir)
                

            os.chdir(ndir)
            if not tmmc_cont:
                a.write_out(filename = "CONFIG")
                origin_ind = int(a.nummol/2)
                ori_dis = round(abs(a.cents[tar_i_n][2]-a.cents[int(a.nummol/2)][2]),2) - self.pull_cont_dis
                end_dis = round(ori_dis + pull_dis,2)
                cont = []
                with open ("CONTROL",'r') as fh:
                    for line in fh: 
                        if "fed order parameter com2" in line:
                            line_sp = line.split()
                            line_sp[-2] = str(end_dis)
                            line_sp[-3] = str(ori_dis)
                            cont.append(" ".join(line_sp) + "\n")
                            fh.readline()
                            fh.readline()
                            cont.append("com2 molecule {}\n".format(tar_i_n+1))
                            cont.append("com2 molecule {}\n".format(origin_ind+1))
                        else:
                            cont.append(line)
                with open ("CONTROL",'w') as fw:
                    for l in cont:
                        fw.write(l)
            if not tmmc_cont:
                self.tmmc_looper(temps, process="tmmc_pull",  **kwargs)
            else:
                self.tmmc_looper(temps, process="tmmc_cont",contdir = "mc_dis", **kwargs)
            os.chdir("..")

    def one_by_one_remove (self, temps, filename = "CONFIG", limit = True, **kwargs):

        a = MonteCon(filename=filename)
        layers,layer_indices,layer_height = a.find_layers()
        print("removing top layer, layer number {} of height {} \
            angstrom containing {} molecules.".format(layers, layer_height[-1],len(layer_indices[-1])))
        
        for temp in temps:
            order = self.removal_order(temp,layer_indices[-1],limit = limit, mode = "singp")
            
            #if limit:
            #    limit_num = int(len(order)/4)
            #    order = order[:limit_num]
            
            for i in order:

                x = a.seek_ind(i)
                y = a.seek_ind(order[-1])
                a.del_mol(i,write_current=True,Tar_mol_label=True)
                cont = []

                with open ("CONTROL",'r') as fh:
                    for line in fh:
                        if "fed order parameter com2" in line:
                            cont.append(line)
                            fh.readline()
                            fh.readline()
                            if i == order[-1]:
                                cont.append("com1 molecule {}\n".format("65"))
                            else:
                                cont.append("com1 molecule {}\n".format("65"))
                            cont.append("com2 molecule {}\n".format(x+1))
                            continue
                        cont.append(line)
                with open("CONTROL",'w') as fw:
                    for j in cont:
                        fw.write(j)
                ndir = "step_mr_{}".format(i)
                os.mkdir(ndir)
                shutil.copy("CONFIG",ndir)
                shutil.copy("CONTROL",ndir)
                shutil.copy("FIELD",ndir)
                os.chdir(ndir)
                self.tmmc_looper([temp], process="tmmc_pull",  **kwargs)
                os.chdir("..")

    def get_pull_results(self,temps):

        record = {}

        for i in temps:

            temp_block = []

            for j in natsort.natsorted(os.listdir()):
                if "step" in j:
                    j_new = j.split("_")
                    os.chdir(j)
                    os.chdir("mc_dis_temp{}_1".format(i))
                    data = np.loadtxt("fe_prof")
                    vac_pull_energy = np.average(data[int(len(data)*0.85):int(len(data)*0.95)])
                    temp_block.append("{} {} {}\n".format(j_new[1],j_new[2],vac_pull_energy))
                    os.chdir("../..")
            
            record[i] = temp_block

        with open ("pull_res", 'w') as fw:

            for i in temps:
                fw.write(str(i) + "\n")
                for j in record[i]:
                    fw.write(j)
            

    def tmmc_nmols(self, temp, nspec = False, bubbleless = False):

        mc = self.mc
        gu = GeneralUtil()
        if not nspec:
            if not isinstance(self.setup_cp,int) and not isinstance(self.setup_cp,float):
                raise ValueError("setup_cp must be an non iterable integer or float for normal nmols.")
        else:
            if not isinstance(self.setup_cp, list):
                raise ValueError("setup_cp must be a list for nspec mc")

        tdir = mc.get_new_run(key = "nmol_temp{}_".format(temp), mode = "normal")
        os.chdir(tdir)
        if nspec:
            self.gcexclude = False
            print ("nspec turns off gcexclude automatically.")

        if self.gcexclude: 
            mc.gcexclude(mol_name = self.slab_molname,mode = self.thickness_mode)
        else:
            cont = []
            with open ("CONTROL",'r') as fh:
                for line in fh:
                    if "use gcexcludeslab" in line:
                        continue
                    else:
                        cont.append(line)
            
            with open("CONTROL",'w') as fw:
                for i in cont:
                    fw.write(i)
                    
        if bubbleless:
            steps = 0
        else:
            steps = self.setup_steps
        gu.edit_anything("temperature",temp,"CONTROL")
        gu.edit_anything("steps",steps,"CONTROL")
        with open ("CONTROL",'r') as fh:
            for line in fh:
                if "fed order" in line:
                    if "nmols" in line:
                        if nspec == True:
                            raise ValueError("nspec invoked but nmols fed order detected in CONTROL file")
                        x = line.split("nmols")
                        bin = float(x[1].split()[0])
                        drangeone = float(x[1].split()[1])
                        drangetwo = float(x[1].split()[2])
                        break
                    elif "nspec" in line:
                        if nspec == False:
                            raise ValueError("nspec not invoked but nspec detected in CONTROL file.")
                        x = line.split("nspec")
                        bin = float(x[1].split()[1])
                        drangeone = float(x[1].split()[2])
                        drangetwo = float(x[1].split()[3])
                        tar_spec_name = fh.readline().split()[0]
                    else:
                        raise ValueError("either nmols or nspecmol must be in the fed directive for this workflow (process).")

        cont = []
        with open("CONTROL",'r') as fh:
                
            for line in fh:

                if line[0] == "#":
                    cont.append(line)
                    continue
                if "fed order" in line:
                    if "param" not in line:
                        lines = line.split()
                        lines.insert(2,"param")
                        line = " ".join(lines) + "\n"

                if nspec:
                    if "move semi" in line:
                        n_semi = int(float(line.split()[2]))
                        if len(self.setup_cp) != n_semi:
                            raise ValueError("setup_cp does not include the same number of chemical potentials as the number of semigrand moves.")
                        semi_freq = int(float(line.split()[3]))
                        x = line.split()
                        x[3] = "10000"
                        x = " ".join(x) + "\n"
                        cont.append(x)
                        nd_count = 0
                        ori_cp = []
                        for i in range(n_semi):
                            x = fh.readline().split()
                            ori_cp.append(" ".join(x) + "\n")
                            if x[0] == tar_spec_name:
                                x[-1] = str(-self.setup_cp[i])
                            elif x[1] == tar_spec_name:
                                x[-1] = str(self.setup_cp[i])
                            else:
                                nd_count += 1
                            cont.append(" ".join(x) + "\n")
                        
                        if nd_count == n_semi:
                            raise ValueError("target nspec species not detected in semi grand directives.")
                    else:
                        cont.append(line)

                else:
                    if "move gcinsertmol" in line:
                        x = fh.readline()
                        chempot = float(x.split()[1])
                        ins_freq = int(line.split()[3])
                        x = x.replace(x.split()[1],str(self.setup_cp))
                        line = line.replace(line.split()[3],"10000")
                        cont.append(line)
                        cont.append(x)
                    else:
                        cont.append(line)

        with open ("CONTROL",'w') as fw:
            for j in cont:
                fw.write(j)

        if self.nmol_strict:
            bin = (int(bin/(self.nw))+1)*(self.nw + 1)
            bg = mc.get_bandgap()[1]
            drangeone = bg[0]
            drangetwo = bg[0] + bin
            mc.edit_windowed_control(drangeone,drangetwo,mode="allrange")

        runs = mc.tmmc(drangeone,drangetwo,bin,self.nw,temp,1,mode="gcmc",ls=False,bubbleless=bubbleless,nwin_crit=self.nwin_crit)
        if not mc.check_in_window(runs):
            print ("Window set up sim failed, terminating run loop...")
            return False
        #subprocess.Popen(['kill' '$()'])
        #jobid = os.environ["SLURM_JOB_ID"]
        #nodelist = os.environ["SLURM_JOB_NODELIST"]
        #nodelist = nodelist[5:-1].split(',')
        
        #with open ("kill.py",'w') as fw:
        #    fw.write("import subprocess\n")
        #    fw.write("subprocess.Popen(['pkill','DLMONTE-SRL.X'])\n")
            
        #for i in nodelist:
        #    subprocess.Popen(['srun','--nodelist={}'.format('nid'+i),'--nodes=1','--ntasks=1','--tasks-per-node=1',\
        #            '--exact', '--mem=1500M', 'python', 'kill.py', '&'])
        if self.nw >= 127:
            jobid = os.environ['SLURM_JOB_ID']
            subprocess.Popen(['scancel','--signal=TERM',jobid])
        else:
            subprocess.Popen(["pkill","DLMONTE-SRL.X"])
        #steps_str = subprocess.check_output(['sstat','-o','JOBID',jobid],universal_newlines=True)
        #for i in steps_str.split()[2:]:
        #    subprocess.Popen(['scancel',i])   
        
        time.sleep(100)

        for i in os.listdir():
            
            if "mc_tmmc" in i:
                os.chdir(i)
                gu.edit_anything("steps",self.tmmc_steps,"CONTROL")
                cont = []
                with open("CONTROL",'r') as fh:
                
                    for line in fh:
                        if line[0] == "#":
                            cont.append(line)
                            continue
                        if nspec:
                            if "move semi"in line:
                                x = line.split()
                                x[3] = str(semi_freq)
                                x[2] = str(n_semi)
                                x = " ".join(x) + "\n"
                                cont.append(x)
                                for i in range(n_semi):
                                    x = fh.readline()
                                    cont.append(ori_cp[i])
                            else:
                                cont.append(line)
                        else:
                            if "move gcinsertmol" in line:
                                x = fh.readline()
                                x = x.replace(x.split()[1],str(chempot))
                                line = line.replace(line.split()[3],str(ins_freq))
                                cont.append(line)
                                cont.append(x)
                            else:
                                cont.append(line)

                with open ("CONTROL",'w') as fw:
                    for j in cont:
                        fw.write(j)
                os.chdir("..")
        os.chdir("..")
        self.tmmc_continue(tdir = tdir, mckey="mc", wipe_prog=True)

    def prod_analysis_looper(self, temps, get_confeng = True, detect_se = True):
        md = MonteData()
        means = []
        ses = []
        if get_confeng:
            self.produce_se("auto")
            nf_se_mean = []
            nf_se_std = []
        for i in temps:
            os.chdir("lsmc_temp{}_1".format(str(i)))
            fes,mean,se = md.fe_analysis("auto", 0.3, detect_se = detect_se)
            if get_confeng:
                with open("se",'r') as fh:
                    dat_line = fh.readline()
                    nf_se = float(dat_line.split()[0])
                    nf_se_ste = float(dat_line.split()[1])
                    nf_se_mean.append(nf_se)
                    nf_se_std.append(nf_se_ste)
            means.append(mean)
            ses.append(se)
            
            os.chdir("..")
        with open ("master_report",'w') as fw:
             
            fw.write("Temperature/free energy(eV)/free energy standard error(eV)(surface energy eV/surface energy standard deviation eV)\n")

            for i in range(len(temps)):
                if get_confeng:
                    fw.write("{} {} {} {} {}\n".format(temps[i],means[i],ses[i],nf_se_mean[i],nf_se_std[i]))
                else:
                    fw.write("{} {} {}\n".format(temps[i],means[i],ses[i]))


    def run_prods (self, temp):
        os.chdir("lsmc_temp{}_1".format(str(temp)))
        
        shutil.copy("lsmc_tmmc1/FEDDAT.000_001","FEDDAT.000_001")
        ls = LSMC(dlm_com=None if self.dlm_exec == "default" else self.dlm_exec,
                  dlm_com_par=None if self.dlm_exec_par == "default" else self.dlm_exec_par)
        mc = self.mc
        gu = GeneralUtil()

        ls.fep_to_feddat("fep")
        dir_lists = []
        for i in range(self.nw):
            tdir = mc.get_new_run(key = "prod", get_feddat=True)
            os.chdir(tdir)
            dir_lists.append(tdir)
            ls.bias_to_prod_control(mode = self.prodmode)
            gu.edit_anything("steps",self.psteps,"CONTROL")
            if self.upper_softedge == True or self.lower_softedge == True:
                bd,vals,windbnd = mc.get_bandgap()
                if self.upper_softedge:
                    vals[1] = 1e10
                if self.lower_softedge:
                    vals[0] = -1e10
                mc.edit_windowed_control(vals[0],vals[1],mode="wrange",softedge=True)
            mc.run_dlm()
            os.chdir("..")
        mc.check_runs_terminate(dir_lists,120, check_m = self.check_m)
    