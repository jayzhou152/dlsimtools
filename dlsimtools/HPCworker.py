import subprocess
from .PolyCore import PolyCore
import time
import shutil
import os

class HPCWorker():
    
    sched = ""
    
    def __init__(self,scheduler):
        
        self.sched = scheduler
        self.write_jobscript (self.sched)
            
    def write_jobscript (self, sched, mods=[], env="", exe="python3 tmmc_master300.py", mpi=1):

        cont =["#!/bin/bash\n"]
        loadcom = "module load {}\n"

        if mpi > 1:
            if sched == "isambard3":
                com = "srun --mpi=cray_shasta -n {} {}\n".format(mpi, exe)
            else:
                com = "mpirun -np {} {}\n".format(mpi, exe)
        else:
            com = exe + "\n"

        if sched == "slurm":
            cont.append("module purge\n")
            cont.append(loadcom.format("gcc"))
            cont.append(loadcom.format("slurm"))
            cont.append(loadcom.format("python3"))
            cont.append(loadcom.format("openmpi"))

        if sched == "archer2":
            cont.append(loadcom.format("cray-python"))

        if sched == "isambard3":
            cont.insert(1, "#SBATCH --output=job_%j.out\n")
            cont.insert(1, "#SBATCH --error=job_%j.err\n")
            cont.append("module purge\n")
            cont.append(loadcom.format("brics/default"))
            cont.append(loadcom.format("brics/userenv"))
            cont.append(loadcom.format("PrgEnv-gnu"))

        for i in mods:
            cont.append(loadcom.format(i))

        if env != "":
            if self.sched == "slurm":
                cont.append("source activate {}\n".format(env))
            elif self.sched == "archer2":
                env_path = os.environ["VIRTUAL_ENV"]
                cont.append("source {}/bin/activate\n".format(env_path))
            elif self.sched == "isambard3":
                cont.append("source ~/miniforge3/bin/activate\n")
                cont.append("conda activate {}\n".format(env))       
        
        cont.append (com)
        sname = "{}script".format(sched)
        
        with open (sname,'w') as fw:
            for i in cont:
                fw.write(i)
        
        return sname
    
    def get_runcom (self, nodes, wtime, sname, mem = 0, Qtype= "free", qos = "standard", prt = "batch-all", premiumcode = "rc-ch1297"):
        
        
        if self.sched == "slurm":
            runcom = "{} {} {} {} {}"

            if not float(nodes).is_integer():
                print ("Number of nodes must be specified as a whole number")
                return ''
            
            runstr = "sbatch"
            if Qtype == "premium":
                runstr += " --account={}".format(premiumcode)
            
            prtstr = "--partition={}".format(prt)
            timestr = "--time={}:00:00".format(wtime)
            nodestr = "--nodes={}".format(nodes)
            
            tarstr = sname
            runcom = runcom.format(runstr,prtstr,nodestr,timestr,tarstr)

        if self.sched == "archer2":

            runcom = "{} {} {} {} {} {} {} {} {}"
            if not float(nodes).is_integer():
                print ("Number of nodes must be specified as a whole number")
                return ''

            runstr = "sbatch"

            prtstr = "--partition={}".format(prt)
            timestr = "--time={}:00:00".format(wtime)
            nodestr = "--nodes={}".format(nodes)
            tpnstr = "--tasks-per-node=1"
            cptstr = "--cpus-per-task=1"
            jobstr = "--job-name=tmmc_jz662"
            qosstr = "--qos={}".format(qos)

            if Qtype == "premium":
                timestr += " --account={}".format(premiumcode)
            if mem != 0:
                timestr += " --mem={}G".format(mem)
            tarstr = sname
            runcom = runcom.format(runstr, jobstr, nodestr, tpnstr, cptstr, timestr, prtstr, qosstr, tarstr)

        if self.sched == "isambard3":

            runcom = "{} {} {} {} {} {} {}"
            if not float(nodes).is_integer():
                print("Number of nodes must be specified as a whole number")
                return ''

            runstr = "sbatch"
            jobstr = "--job-name=tmmc_jz662"
            nodestr = "--nodes={}".format(nodes)
            timestr = "--time={}:00:00".format(wtime)
            # grace partition has 144 cores per node
            ntasksstr = "--ntasks={}".format(int(nodes) * 144)
            prtstr = "--partition=grace"
            qosstr = "--qos=grace_qos"

            if Qtype == "premium":
                runstr += " --account={}".format(premiumcode)
            if mem != 0:
                ntasksstr += " --mem={}G".format(mem)
            tarstr = sname
            runcom = runcom.format(runstr, jobstr, nodestr, ntasksstr, timestr, prtstr, qosstr, tarstr)

        return runcom
    
    def move_script(self, sname, tardir):
        
        shutil.copy(os.path.join(os.curdir,sname),os.path.join(tardir,sname))
    
    def submit_job (self, runcom, scancel = False):

        print("Submitting job with command:\n  {}".format(runcom))
        try:
            ostr = subprocess.check_output(runcom.split(), stderr=subprocess.STDOUT, timeout=30)
        except subprocess.TimeoutExpired:
            print("ERROR: sbatch timed out after 30s — check HPC connectivity or command")
            raise
        except subprocess.CalledProcessError as e:
            print("ERROR: sbatch failed (exit {}): {}".format(e.returncode, e.output.decode("ascii", errors="replace")))
            raise
        ostr = ostr.decode("ascii")
        
        if self.sched in ("slurm", "archer2", "isambard3"):
            jobid = ostr.split()[len(ostr.split())-1]
        
        if scancel:
            subprocess.run (["scancel",jobid])
        
        return jobid
    
    def status_check (self, jobid, userid= "jz662"):
        
        qccom = "squeue --user={}".format(userid) 
        ostr = subprocess.check_output (qccom.split(), stderr = subprocess.STDOUT)
        ostr = ostr.decode ("ascii")
        ostr = ostr.split("\n")
        for i in ostr:
            if not i.split():
                continue
            if i.split()[0] == jobid:
                return True
            
        return False
    
    
    
    def HPCseqlooper(self, vals, kw, step, wtime = 6, nodes=6, prt = "batch-all", restart_key= -1):
        
        if prt == "batch-sky":
            mpi = nodes*24
        else:
            mpi = nodes*16
        
        sname = self.write_jobscript (self.sched,mpi=mpi)
        runcom = self.get_runcom (wtime = wtime, nodes=nodes, prt = prt)
        dlp = PolyCore()
        cwdir = os.getcwd()
        
        for i in range(len(vals)):
            
            val = vals[i]
            print("(sequel) running {} {} simulation on DL_Poly, please wait...".format(kw,val))
            
            if i == 0:
                val = str(val)
                if restart_key <= 0:
                    
                    nwdir = dlp.get_new_run(kw[0]+val, mode ="seq")
                    self.move_script (sname,nwdir)
                    os.chdir(nwdir)
                    
                    self.dlp.edit_control("steps",str(step))
                    self.dlp.edit_control(kw,val)
                else:
                    nwdir = "dlprunseq_t{}_{}".format(val,restart_key)
                    os.chdir (nwdir)
                    removal_list = ["HISTORY","OUTPUT","REVCON","REVIVE","STATIS"]
                    for j in removal_list:
                        if os.path.exists(j):
                            os.remove(j)
                    
                jobid = self.submit_job(runcom)
                
                while True:
                    time.sleep(60)
                    if self.status_check(jobid):
                        continue
                    else:
                        break
                os.chdir(cwdir)
            
            else:
                
                val = str(val)
                nwdir = dlp.get_next_run(kw[0]+val,nwdir)
                self.move_script (sname,nwdir)
                os.chdir(nwdir)
                
                self.dlp.edit_control(kw,val)
                cur_step = dlp.get_steps()
                new_step = int(cur_step) + step
                dlp.edit_control("steps",str(new_step))
                dlp.restart()
                
                jobid = self.submit_job(runcom)
                
                while True:
                    time.sleep(60)
                    if self.status_check(jobid):
                        continue
                    else:
                        break
                
                os.chdir(cwdir)
