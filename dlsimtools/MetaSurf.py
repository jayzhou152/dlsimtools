#!/usr/bin/env python
# coding: utf-8
""" Utility for automation of integrated surface energy collection with metadise and dl_poly
update: surface length now reboxed
slice relaxation currently has termination checks and dl_field +5 corrections
, does not rebox.

"""
from .FieldTools import FieldTools
from .PolyCore import PolyCore
from .PolyOutput import PolyOutput
from .InputConverter import InputConverter
from pymatgen.io import cif
from pymatgen.core.surface import Structure,Lattice
from .CifTools import CifTools
from .MetaUtil import MetaUtil
import subprocess
import numpy as np
import os
import shutil

class MetaSurf:
    
    dlf = FieldTools()
    dlp = PolyCore()
    dlo = PolyOutput()
    dlpf = InputConverter()
    mu = MetaUtil()

    def __init__ (self):
        pass

    def change_cellmat(self, mat, file = "CONFIG", mode = "normal"):
        
        mat = np.array(mat)
        if mat.shape != (3,3):
            raise ValueError("Input matrix is not a 3 by 3 array.")
        
        cellmat_str = ["{:>20}{:>20}{:>20}\n".format("%.10f"%mat[0,0],"%.10f"%mat[0,1],"%.10f"%mat[0,2]),
                       "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[1,0],"%.10f"%mat[1,1],"%.10f"%mat[1,2]),
                       "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[2,0],"%.10f"%mat[2,1],"%.10f"%mat[2,2])]
        
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
    def prep_inp (self, mi, stack, atmnum, file = "input.txt", clean=False, search = True):
        
        col = []
        p1 = "check \nnosort\nnoshift\nmolecule\nslab\nmiller {}\nstart\nstackgen {} \nnewt 0 \nbfgs 1 \nstart \nstop".format(mi,stack)
        FNULL = open(os.devnull,'w')
        
        if clean:
            subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT, shell = True)

        with open (file,'r') as fh:
            for line in fh:
                if "check" in line:
                    break
                else:
                    col.append(line)
        col.append(p1)
        
        with open ("input.txt",'w') as fw:
            for i in col:
                fw.write(i)
                
        self.mu.run_metadise()
        if search:
            flag, mc, ac = self.mu.fres_check(atmnum)

            if not flag:
                subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT, shell = True)
                self.mu.molrad_optimiser(atmnum)
                self.mu.run_metadise()
    
    def prod_meta_looper (self, surfaces, atmnum, nsim =1,  **kwargs):
        
        dirs_a = []
        count = 0
        for i in surfaces:
            dirs = self.prod_meta(i,atmnum, mode = "relaxed",**kwargs)
            dirs_a += dirs
            count += 1
            print(count)
            if count == nsim:
                self.dlp.check_term_looper(dirs_a)
                dirs_a = []
                count = 0
            os.chdir("..")


    def prep_staco (self, clean=False, sliceprep = False, file="input.txt", cutoff = [20,20]):
        
        col = []
        p1 = "ends\nnopoten\n{}slab\nstart\nminimise\nconj 0\nbfgs 1\nnewt 0\nmaxu 20\nnokill\nstart\nanal\nstart\nstop\n" 
        c_vector=[]
        struc = ""
        FNULL = open(os.devnull,'w')
        
        if clean:
            subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT, shell = True)
        
        with open (file, 'r') as fh:
            
            for i in range (1000000):
        
                if i == 1:
                    struc = fh.readline()
                    col.append(struc)
                    struc = struc.split()
                
                elif i == 2:
                    info = fh.readline()
                    col.append(info)
                    info = info.split()[1:5]
                
                elif (4 < i < 8):
                    itemp = fh.readline()
                    c_vector.append(itemp.split())
                    col.append(itemp)
                
                else:
                    itemp = fh.readline()
                    if "ends" in itemp:
                        break
                    else:
                        col.append(itemp)

        c_vector = np.asarray(c_vector,dtype=float)
        c_vector_new = c_vector.copy()
        growth = [1,1,1]
        
        for i in range (3):
            
            mul = 1  
            x, y = i , i+1

            if i == 2:
                if sliceprep:
                    x, y = 2, 0
                else:
                    break
            if i < 2:
                while c_vector_new[x,y] < cutoff[i]:
                    mul += 1
                    c_vector_new[i] += c_vector[i]
                    
            if i == 2 and mul<3:
                growth [2] = 3
                c_vector_new[2] = c_vector[2] * 3
            else:
                growth [i] = mul
        
        c_vector = c_vector_new
        

        if sliceprep:
            if growth[2] <3:
                growth[2] = 3
            
            growth_string = "check\nprint dlpoly 1\ngrow {} {} 1\nstart\nmeta\ncheck\ngrow 1 1 {}\nstart\n".format(growth[0],growth[1],growth[2])
            
        else:
            growth_string = "check\nprint dlpoly 1\ngrow {} {} {}\n".format(growth[0],growth[1],growth[2])
        
        p1 = p1.format(growth_string)
        
        c_vector[[0,2]] = c_vector[[2,0]]
        c_vector[[1,2]] = c_vector[[2,1]]
        
        col.append(p1)
        
        with open("input.txt",'w') as fw:
            for i in col:
                fw.write(i)
        FNULL = open(os.devnull,'w')
        self.mu.run_metadise() 
        
        if sliceprep:
            return growth[2], c_vector
        else:
            return c_vector
    
    def prod_meta (self, surface, atmnum, steps=500000, cut_off = [20,20], mode = "singp", ff="opls2005", solvent="none"):
        """_summary_

        Args:
            surface (_type_): _description_
            atmnum (_type_): _description_
            steps (int, optional): steps for energy relaxation during for "relaxed" mode. Defaults to 500000.
            cut_off (list, optional): minimum allowed dimension length [xmin,ymin]. Defaults to [20,20].
            mode (str, optional): single point energy("singp") for lowest energy unit cell, or "relaxed" energy \
            for lowest energy unit cell. Defaults to "singp".
            ff (str, optional): force field for solid. Defaults to "opls2005".
            solvent (str, optional): force field for solvent. Defaults to "none".

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        import warnings
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF")
        
        ndir = "surface_" + surface
        
        os.mkdir(ndir)
        shutil.copy("input.txt",os.path.join(ndir,"input.txt"))
        shutil.copy("CONTROL",os.path.join(ndir,"CONTROL"))
        os.chdir(ndir)
        self.dlp.edit_control("steps",str(steps))
        self.prep_inp(surface,"systematic", atmnum, clean=True)
        self.dlp.edit_control("ensemble","npt")
        codes = self.get_codes()
        eng = []
        if mode == "relaxed":
            dirs = []
        print ("Starting code search for surface {}.".format(surface))
        for j in codes:
            
            self.prep_inp(surface, j, atmnum, clean=True)
            os.mkdir(j)
            shutil.copy("CONTROL", os.path.join(j,"CONTROL"))
            shutil.copy("staco0001sl.out", os.path.join(j,"input.txt"))
            os.chdir(j)
            c_vector = self.prep_staco(clean=True, cutoff = cut_off)
            surf_a= c_vector[1,1] * c_vector[2,2]
            mat = self.get_surf_mat()
            self.dlf.change_vector(mat)
            self.dlf.get_1mol_xyz(atmnum,filename = "bef_o0001sl.xyz")
            self.dlf.copyTofield(filename="1mol.xyz")
            self.dlf.edit_dlfc(10,"1mol.xyz")
            self.dlf.edit_dlfc(3,ff)
            self.dlf.edit_dlfc(35,solvent)
            print ("Running force field generation.")
            self.dlf.run_dlf()
            self.dlf.get_output("config")
            self.dlf.get_output("field")
            os.rename("CONFIG","CONFIG_sol")
            os.rename("FIELD","FIELD_sol")
            self.dlf.edit_dlfc(35,"none")
            self.dlf.run_dlf()
            self.dlf.get_output("field")
            os.rename("input.txt","input.backup")
            bcif = CifTools("bef_o0001sl.cif")
            a,b,c = bcif.get_abc()
            al,be,ga = bcif.get_angles()
            cell_str = "          {}   {}   {}   {}   {}   {}".format(a,b,c,al,be,ga)
            self.dlpf.xyz_to_polycon("bef_o0001sl.xyz",cell_str)
            os.rename("config__o0002no.dlp","CONFIG")
            self.dlpf.map_config()
            self.change_cellmat(mat)
            os.rename("input.txt","input.xyz")
            os.rename("input.backup","input.txt")

            run_dir = self.dlp.get_new_run("surf",os.curdir)
            if mode == "singp":
                os.chdir(run_dir)
                self.dlp.edit_control("steps",str(1))
                print ("Collecting single point energy for surface {} code {} ...".format(surface,j))
                self.dlp.run_poly(mode="fg")
                sing_p = float(self.dlo.get_one(0,"engcfg"))/surf_a
                eng.append(sing_p)                
                os.chdir("../..")
            elif mode == "relaxed":
                os.chdir(run_dir)
                print ("Running surface{} code{} relaxtions in the background...".format(surface,j))
                self.dlp.run_poly(mode="bg")
                dirs.append(os.getcwd())
                os.chdir("../..")
            else:
                raise ValueError("Mode not recognised..")
        
        if mode == "relaxed":
            return dirs
        tar_code = codes [eng.index(min(eng))]
        print ("Energy minimum is found in code {}.".format(tar_code))

        with open ("SurfEng", 'w') as fw:
            for i,j in zip(codes,eng):
                fw.write("EngUnit = eV/A^2 \n")
                fw.write("code " +i + "    ")
                fw.write(str(j) +"\n")
        
        return tar_code

    def surflooper(self,surfaces,atmnum,steps=500000, skipmeta=[], cut_off=[20,20], dosurface=True, doslice = True, clean = False):
        
        import warnings
        warnings.filterwarnings("ignore", message="Issues encountered while parsing CIF")

        if clean:
            self.surf_clean()

        for i in surfaces:

            if len(skipmeta) == 0:
                target_code = self.prod_meta(i, atmnum,cut_off=cut_off)
                os.chdir(target_code)
            else:
                os.chdir("surface_" + i)
                ind = surfaces.index(i)
                target_code = skipmeta[ind]
                os.chdir(target_code)

            if dosurface:
                print ("Surface energy sim for surface {} code {} is running in the background".format(i,target_code))
                run_dir = self.dlp.get_new_run("tar_surf",os.curdir)
                os.chdir(run_dir)
                self.dlp.edit_control("steps", str(steps))
                self.dlp.run_poly(mode="bg")
                os.chdir("..")

            
            if doslice:
                print ("Starting slice relaxation sim for surface {}...".format(i))
                
                stack_num, c_vector = self.prep_staco(sliceprep=True, clean=True, cutoff = cut_off)
                print(c_vector)
                self.dlf.change_vector(c_vector)
                self.dlf.prep_xyz("af_co0001no.xyz")
                self.dlf.copyTofield(filename="af_co0001no.xyz")
                self.dlf.edit_dlfc(10,"af_co0001no.xyz")
                self.dlf.run_dlf()
                self.dlf.get_output("config")
                self.dlf.get_output("field")
                run_dir = self.dlp.get_new_run("slice",os.curdir)
                os.chdir(run_dir)
                with open ("stackn" , 'w') as fw:
                    fw.write(str(stack_num)+"\n")
                    fw.write(str(atmnum))
                self.dlp.edit_control("steps",str(steps))
                self.dlp.edit_control("ensemble","nst")
                self.dlp.run_poly(mode="bg")
                if self.dlp.check_term_error(wait=True) == "95":
                    c_vector = np.tril(np.triu(c_vector))
                    c_vector[0,0] += 5
                    print ("Correcting slice config by adding contingency vacuumn space...")
                    self.dlf.change_vector(c_vector)
                    self.dlf.run_dlf()
                    self.dlf.get_output("config")
                    self.dlf.get_output("field")
                    self.dlp.run_poly(mode="bg")
                os.chdir("..")
                print ("Slice relaxation sim has started in the background for surface {}".format(i))
            
            os.chdir("../..")
            

    def surf_clean (self):
        for i in os.listdir():
            if "surface" in i:
                shutil.rmtree(i)

    def unit_surfa (self):
        
        with open ("input.txt",'r') as fh:
            
            for i in range (5):
                fh.readline()
            
            temp = fh.readline()
            x = float(temp.split()[1])
            temp = fh.readline()
            y = float(temp.split()[2])
            
            unit_a = x*y
        return unit_a
    
    def calc_eng (self,bulk_eng,atmnum,dosurface=True,doslice=True):
        output = ["Energy unit: eV (surface eV/A^2,attachment eV/molecule) \n"]
        
        for i in os.listdir():
            if "surface" in i:
                try:
                    os.chdir(i)
                except OSError:
                    continue
                
                mi = i.split('_')[1]
                self.direct()
                dsp = self.get_dsp()
                o_string = "surface {}     {}".format(mi,dsp)
                
                if dosurface:
                    print(os.getcwd())
                    os.chdir("dlprun_tar_surf_1")
                    surf_tot = self.dlo.get_av("engcfg")
                    nummol = self.dlp.get_molnum(atmnum)
                    bulk_tot = bulk_eng*nummol
                    surf_a = self.dlp.surf_a()
                    surf_eng = (surf_tot - bulk_tot)/surf_a
                    
                    os.chdir("..")
                    o_string += "        {}".format(str(surf_a))
                    unit_surfa = self.unit_surfa() 
                    o_string += "        {}".format(str(unit_surfa))
                    o_string += "        {}".format(str(surf_eng))
                    
                if doslice:
                    os.chdir("dlprun_slice_1")
                    attach_eng = self.attach_eng()
                    o_string += "        {}".format(str(attach_eng))
                    os.chdir("..")
                    
                os.chdir("..")
                os.chdir("..")
                output.append(o_string + "\n")
                #print(os.getcwd())
                
        with open("SURFENG", 'w') as fw:
            for i in output:
                fw.write(i)
    
    def attach_eng(self):
        
        bulk_eng = self.dlo.get_av("engcfg")
        with open ('stackn','r') as fh:
            stack = int(fh.readline())
            atmnum = int(fh.readline())
        bulk_slice = bulk_eng/stack
        
        if "singp" in os.listdir():
            shutil.rmtree("singp")
        
        os.mkdir("singp")
        self.dlp.move_input("singp",files = ["CONTROL","FIELD","REVCON"])
        os.chdir("singp")
        slice_atm_num = self.dlpf.cut_config(stack)
        nummol_slice = slice_atm_num/atmnum
        self.dlpf.cut_field(stack, retain = True)
        self.dlp.edit_control("steps","1")
        self.dlp.edit_control("ensemble","npt")
        self.dlp.run_poly(mode="fg")
        vacumn_slice = float(self.dlo.get_one(0,"engcfg"))
        
        atta_eng = (vacumn_slice - bulk_slice)/nummol_slice
        os.chdir("..")
        
        return atta_eng
        
    def direct (self):
        retrace = os.getcwd()
        for i in os.listdir():
            try:
                os.chdir(i)
                for j in os.listdir():
                    if "tar_surf" in j or "slice" in j:
                        return 1
                os.chdir(retrace)
                
            except OSError:
                continue
        return 0
    
    def direct_tar (self):
        
        with open ("misc", 'r') as fh:
            tar_code = fh.readline()
            
        x = tar_code.split()[1]
        
        print ("misc tracking found target code {}".format(x))
        os.chdir(x)
        
        return x
    
    def get_codes(self):
        
        file= "summ_o0001.out"
        codes = []
        
        with open (file,'r') as fh:
            for i in range(5):
                fh.readline()
            for line in fh:
                codes.append(line.split()[1])
        return codes
    
    def get_dsp (self,file="input.txt",directory = os.curdir):
        
        os.chdir(directory)
        with open (file,'r') as fh:
            for i in range(7):
                fh.readline()
            temp = fh.readline()
            dsp = float (temp.split()[0])
        
        return dsp
    
    def read_SURFENG (self,mode="attach",filename="SURFENG"):
        surfinfo = {}
        
        with open (filename,'r') as fh:
            fh.readline()
            if mode == "attach":
                ind = 8
            if mode == "surf":
                ind = 7
            for line in fh:
                temp = line.split()
                surfinfo[tuple(map(int,temp[1:4]))] = float(temp[ind])

        return surfinfo

    def get_matrix (self,filename ="bef_o0001sl.cif"):
        
        cif_parse = cif.CifParser(filename)
        struct = cif_parse.get_structures(primitive=False)[0]
        mat = struct.lattice.matrix
        tol = 1e-10
        mat.setflags(write=1)
        mat[abs(mat) < tol] = 0.00
        return mat

    def get_surf_mat (self, filename = "bef_o0001sl_res.out"):
        
        
        with open (filename,'r') as fh:
            fh.readline()
            surf_mat = []
            for i in range(3):
                vec = fh.readline()
                vec = list(map(float,vec.split()[1:]))
                surf_mat.append(vec)
            surf_mat = np.array(surf_mat)
        
        return surf_mat
    
    def check_runs (self):
        
        print("Checking if all surface MD runs terminated properly..")
        count = 0
        
        for i in os.listdir():
            
            if "surface" in i:
                try:
                    os.chdir(i)
                except OSError:
                    continue
            
                mi = i.split('_')[1]
                self.direct()
                
                for j in ["dlprun_tar_surf_1","dlprun_slice_1"]:
                    
                    try:
                        os.chdir(j)
                    except OSError:
                        continue
                    
                    if j == "dlprun_tar_surf_1":
                        os_m = "surface"
                    if j == "dlprun_slice_1":
                        os_m = "slice"
                    
                    
                    if self.dlp.check_term():
                        
                        O_str = "Surface {} {} simulation has terminated successfully.".format(mi, os_m)
                    
                    else:
                        O_str = "Surface {} {} simulation has failed due to the error above.".format(mi, os_m)
                        count += 1
                    
                    print (O_str)
                    os.chdir ("..")
                os.chdir("../..")
        if count == 0:
            print ("All surface runs terminated successfully.")
        else:
            print ("{} runs terminated with an error.".format(count))
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                