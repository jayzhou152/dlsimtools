"""
Core functionalities for editing, handling  and ultimately automating dl_poly runs

The dlpcore class is a fundamental class containing multiple small ultilities useful 
for any programs involving DL_POLY(Molecular Dynamics MD) calcs.

"""

import os
import shutil
import subprocess
import time
import socket

class PolyCore:
        
    # Instantiate class with a message
    def __init__(self):
        pass
    
    
    def run_poly(self, mode = "bg", np = 4, runcom = "DLPOLY_MPI.X"):
        
        """
        Runs DLPOLY.X on command line, which initiates one single MD run. stdout
        and stderr are currently supressed for readable output. could enable them 
        simply by deteleting second and third arugument of suprocess.run() call.
        
        Arguments:
            mode (str) => mode for running command, bg correspond to background
            minimized calculation, fg corresponds to foreground shell occupying
            runs, fg mode is mainly used for procedural logic where the MD run and
            result collection are bound together. bg mode is mainly used for lengthy sims
            where users would prefer shell is not occupied during the time of the run
            and result collection is detached from the MD run.
            
        Returns:
            n/a
        """
        
        files = ["CONTROL", "CONFIG", "FIELD"]
        FNULL = open(os.devnull,'w')
        if "balena" in socket.gethostname():
            mode = "hpc"
            
        for i in files:
            if i not in os.listdir(os.curdir):
                print("{} file missing".format(i))
                return -1
        try:
            if mode == "bg":
                subprocess.run("DLPOLY.X &", stdout = FNULL, stderr = subprocess.STDOUT, shell=True)
            elif mode == "fg":
                subprocess.run("DLPOLY.X", stdout = FNULL, stderr = subprocess.STDOUT)
            elif mode == "mpi":
                runstr = "mpirun -np {} {} &".format(np, runcom)
                subprocess.run(runstr, stdout = FNULL, stderr = subprocess.STDOUT, shell=True)
            elif mode == "mpipara":
                runstr = "mpirun -np {} {}".format(np, runcom)
                subprocess.Popen(runstr, stdout = FNULL, stderr = subprocess.STDOUT)
            elif mode == "hpc":
                hw = HPCworker()
                runcom = hw.get_runcom()
                jobID = hw.submit_job(runcom)
                return jobID
            else:
                print("Mode not recognised, available modes: 'bg','fg','mpi' or 'mpipara'.")
    
        except OSError:
            print("DLPOLY.X is not a valid command, check installation and PATH setup")
            raise
    
    def edit_control (self, kw, newval, filename = "CONTROL", directory = os.curdir):
        
        """
        Edits DLPOLY control file using keyword and values. NOTE this function 
        only changes the first value that comes after the keyword, keywords with
        multiple values assigned will need to be changed manually.
        
        Arguments:
            kw (str) ==> keyword(directive) of the target variable, check an
            example control file for what type of keyword is available to change
            
            newval (str) ==> new value to change to.
            
        Returns:
            n/a
        
        """

        cont = []
        file = os.path.join(directory,filename)
        
        with open(file,'r') as fh:
            for i in range(10000):
                temp = fh.readline()
                if len(temp.split())>0:
                    if temp.split()[0] == kw:
                        temp = temp.replace(temp.split()[1],str(newval))
                    cont.append(temp)
        
        with open(file,'w') as fw:
            for i in cont:
                fw.write(i)
    
    
    def move_input(self, tar_dir, cur_dir = os.curdir, files = ["CONTROL", "CONFIG", "FIELD"]):
    
        """
        Move inputs from one directory to another.
        
        Arguments:
            tar_dir (str,pathlike) ==> target directory
            cur_dir (str,pathlike) ==> source directory, defaults to the current 
            working directory.
            files (list(str)) ==> file names to move, defaults to the standard
            DLPOLY input files.
            
        Returns:
            n/a
        """
        
        for i in files:
            # REVCON always moved as CONFIG, saves time for renaming.
            if i == "REVCON":
                tar_filename = "CONFIG"
            # REVIVE always moved as REVOLD.
            elif i == "REVIVE":
                tar_filename = "REVOLD"
            elif "FIELD" in i:
                tar_filename = "FIELD"
            else:
                tar_filename = i
            
            try:
                shutil.copy(os.path.join(cur_dir,i),os.path.join(tar_dir,tar_filename))
            
            except FileNotFoundError:
                print("{} file not in current directory".format(i))
    
    
    def get_new_run(self,key,directory = os.curdir,mode="para",files = ["CONTROL", "CONFIG", "FIELD"]):
        
        """
        Create a directory/folder for new DLPOLY runs, with a specific folder
        naming system.
        
        This new_run function defaults to creating a folder called dlprun_@key_@runID.
        When called during sequel DLPOLY runs to intiate first in the sequence of runs,
        this function creates a folder called dlprunseq_@key_@runID to keep consistency
        with follow_up runs created by @get_next_run.
        
        runID is created based on the number of runs currently exist with the same key,
        for example if dlprun_temp_1 exists, this function would create a folder for
        dlp_run_temp_2.
        
        
        Arguments:
            key (str) ==> a short keyword description for the DLPOLY run.
            directory (str,pathlike) ==> the folder to create this run in.
            mode (str) ==> whether it is a sequel run or parallel run. "seq" for sequel,
            defaults to parallel mode "para".
        
        Returns:
            ndir (str,pathlike) ==> directory of the created run.
            0 if no run was created.
            -1 if input arg error.
        """
        
        for i in range(10000):
            
            if i > 20:
                print("you have too many simulations in current folder, consider using cleanup()")
            
            if mode == "para":
                ndir_name = "dlprun_{}_{}"
            elif mode == "seq":
                ndir_name = "dlprunseq_{}_{}"
            else:
                print ("Run creation only accepts mode 'seq' or mode 'para'")
                return -1
            
            ndir_name = ndir_name.format(key,str(i+1))
            ndir = os.path.join(directory,ndir_name)
            
            if not os.path.exists(ndir):
                os.mkdir(ndir)
                self.move_input(ndir,files=files)
                return ndir
            
        return 0
    
    
    def get_next_run (self,key,source_dir,directory = os.curdir):
        
        """
        Create new runs for sequel runs, defaults to folder name dlprunseq_@key_@runID.
        Works with the restart directive in DLPOLY, utilises the REVCON and REVIVE files.
        
        Arguments:
            key (str) ==> a short keyword description for the DLPOLY run.
            source_dir (str,pathlike) ==> folder to fetch the continuation run files from.
            directory (str,pathlike) ==> the folder to create this run in.
            
        Returns:
            ndir ==> directory of the created run.
            0 if nothing is created.
        """
        
        for i in range(10000):
            
            if i > 20:
                print("you have too many simulations in current folder, consider using cleanup()")
            
            ndir_name = "dlprunseq_{}_{}"
            ndir_name = ndir_name.format(key,str(i+1))
            
            ndir = os.path.join(directory,ndir_name)
            
            if not os.path.exists(ndir):
                os.mkdir(ndir)
                self.move_input(ndir,cur_dir = source_dir, files = ["CONTROL","FIELD","REVCON","REVIVE"])
                return ndir
            
        return 0
    
    
    def get_steps (self, file = "CONTROL"):
        
        """
        Fetch the current number of steps, this is used for sequel runs for incremental increase
        in the number of timesteps.
        
        Arguments:
            file (str)(optional) ==> file name defaults to "CONTROL".
            
        Returns:
            number of steps (str).
        """
        
        with open(file,'r') as fh:
            for line in fh:
                if "steps" in line:
                    return line.split()[1]

    def paralooper (self, vals, kw , mpi=False, npcount = 2):
        
        """
        Run multiple simultaneous DLPOLY simulations with a set of values for a 
        chosen variable.
        
        Arguments:
            vals (list(str)) ==> list of string values to loop through.
            kw (str) ==> keyword (variable) that correpsonds to the list of values
            provided.
            
        Returns:
            none.
        """
        
        cwdir = os.getcwd()
        
        for val in vals:
            val = str(val)
            nwdir = self.get_new_run(kw[0]+val)
            os.chdir(nwdir)
            self.edit_control(kw,val)
            print("trying to run {} {} simulation on DL_Poly...".format(kw,val))
            if mpi:
                # mpi runs here is usually not favoured, only use if doesn't overload core thread count.
                self.run_poly (mode = "mpipara",np=npcount)
            else:
                self.run_poly()
            os.chdir(cwdir)
    
    
    def restart (self,scale=False):
        
        """
        Used to simply insert restart directive to the CONTROL file in CWD, does
        not insert if restart is already present.
        
        Arguments:
            scale(boolean) ==> temperature scaling (starts a new sim)
        
        Returns:
            0 if restart is already a directive in the CONTROL file.
        """
        
        cont = []
        with open ("CONTROL", 'r') as fh:
            for line in fh:
                if line.split()[0] == "restart":
                    continue
                else:
                    cont.append(line)

        if scale:
            cont.insert(len(cont)-1,"restart scale\n")
        else:
            cont.insert(len(cont)-1,"restart\n")
            
        with open ("CONTROL",'w') as fw:
            for i in cont:
                fw.write(i)
            
    def seqlooper (self, vals, step, kw, mpi=False, npcount=8):
        
        """
        Run consecutive DLPOLY runs with the restart directive, mainly used to simulate
        real time changes of a certain variable.
        
        Arguments:
            vals (list(str)) ==> list of string values to loop through.
            step (int) ==> number of steps to increment each run, currently only
            constant increments are allowed.
            kw (str) ==> keyword (variable) that correpsonds to the list of values
            provided.
            
        Returns:
            none
            
        """
        
        cwdir = os.getcwd()
        
        for i in range(len(vals)):
            val = vals[i]
            print("(sequel) running {} {} simulation on DL_Poly, please wait...".format(kw,val))
            if i == 0:
                val = str(val)
                nwdir = self.get_new_run(kw[0]+val,mode ="seq")
                os.chdir(nwdir)
                self.edit_control("steps",str(step))
                self.edit_control(kw,val)
                if mpi:
                    self.run_poly (mode = "mpi", np=npcount)
                else:
                    self.run_poly(mode = "fg")
                print("hi first step here")
                os.chdir(cwdir)
            else:
                val = str(val)
                nwdir = self.get_next_run(kw[0]+val,nwdir)
                os.chdir(nwdir)
                self.edit_control(kw,val)
                cur_step = self.get_steps()
                new_step = int(cur_step) + step
                self.edit_control("steps",str(new_step))
                self.restart()
                if mpi:
                    self.run_poly (mode = "mpi", np=npcount)
                else:
                    self.run_poly(mode = "fg")
                os.chdir(cwdir)
        print ("sequel runs complete")
    
    def get_molnum (self,atomnum,file =  "FIELD" ):
        
        """
        Used to obtain the number of molecules contained in current DLPOLY run.
        Regardless of whether the FIELD is constructed as multiple molecules defined
        as a single molecule or they are defined as a single molecule.
        
        Arguments:
            atomnum (int) ==> number of atoms in a molecule, needs to be user defined
            as occasionally FIELD file could defined multiple molecules as one and this
            parameter is neccesary to work out the correct nummol.
            file (str)(optional) ==> defaults to FIELD file.
        
        Returns:
            nummols (int) ==> number of molecules in FIELD file
        
        """
                
        with open (file,'r') as fh:
            for i in range(4):
                fh.readline()
            nummols = fh.readline().split()[1]
            if nummols == "1":  # still checks for nummol when nummol given in FIELD is 1.
                atoms = int(fh.readline().split()[1])
                return atoms/atomnum
            else:
                return int(nummols)
            
    def surf_a(self,file = "CONFIG"):
        
        """
        Calculates the surface area of a slab.
        
        Arguments:
            file (str) (optional) ==> defaults to seek info in CONFIG file
        
        Returns:
            surf_a (float) ==> surface area of the slab (top and bottom area)
        """
        
        with open (file,'r') as fh:
            for i in range(3):
                fh.readline()
            temp = fh.readline()
            x = float(temp.split()[1])
            temp = fh.readline()
            y = float(temp.split()[2])
        
        surf_a = 2*x*y
        return surf_a
    
    def check_term_error (self, mode = "ini", wait = False):
        
        """
        utility to check for unexpected run terminations.
        
        """
        if wait:
            time.sleep(60)
        
        with open ("OUTPUT",'r') as fh:
            for line in fh:
                if "error" in line:
                    err_code = line.split()[-1]
                    for moreline in fh:
                        if "error" in moreline:
                            print(moreline)
                    return err_code
        
        return ''
    
    def check_term (self, file = "OUTPUT", gettime=False):
        
        term_str = "time elapsed"

        while True:
            time.sleep(60)
            count = 0
            with open ("OUTPUT",'r') as fh:
                for line in fh:
                    if "error" in line:
                        for moreline in fh:
                            if "error" in moreline:
                                print(moreline)
                        return False
                    if term_str in line:
                        time = float(line.split()[6])
                        count += 1
                        if count == 2:
                            if gettime:
                                return time
                            else:
                                return True
    
    def check_term_looper (self, folders):

        nfolders = len(folders)
        succount = 0

        for i in folders:
            os.chdir(i)
            if self.check_term:
                succount += 1
            else:
                print ("Simulation in folder {} has encountered an error.".format(i))
            os.chdir("../..")

        if succount == nfolders:
            print("All dl_poly process check completed. Currently in folder{}.".format(os.getcwd()))
            return

    
    def cleanup(self, check, directory = os.curdir):
        
        """
        Utility to clean up DLPOLY run folders.
        
        Arguments:
            check (str) ==> only clean if passed the string "yes", act as a safeguard
            check.
            directory (str,pathlike) (optional) ==> defaults to clean runs in current
            directory.
            
        Returns:
            none.
        """
        if check == 'yes':
            for i in os.listdir(directory):
                if "dlprun" in i:
                    shutil.rmtree(i)
                    
    def killall (self,check):
        
        """
        Utility to kill all DLPOLY runs running in shell background. Shell independent.
        
        Arguments:
            check (str) ==> only clean if passed the string "yes", act as a safeguard
            check.
            
        Returns:
            none
        """
        
        if check == 'yes':
            subprocess.run ("ps -A | grep 'DLP'| awk '{print $1}' | xargs kill".split())
    
