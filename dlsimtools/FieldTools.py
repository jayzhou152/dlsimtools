import os
import shutil
import subprocess
import numpy as np
import random
from .MetaUtil import MetaUtil
from .CifTools import CifTools
"""
Utility for automating DL_FIELD output generation, this scripts could also solve 
problems related to molecular definition in DL_FIELD output file.


keys (dl_field control)

1 => whether to construct dlp output
2 => unused slot
3* => force field: charmm, amber, opls2005, dreiding, pcff, charmm19, cvff, inorganic, multiple
4 => Energy unit: kcal/mol, kJ/mol, eV, or K
5 => Conversion criteria (strict, normal, loose)
6 => Bond type (0=default, 1=harmonic , 2=Morse)
7 => Angle type (0=default, 1=harmonic, 2=harmonic cos)
8 => Include user-defined information. Put 'none' or a .udff filename
9 => Verbosity mode: 1 = on, 0 = off   
10* => Configuration file name
11 => Output file in PDB. Put 'none' if not needed
12 => Solution Maker: on/off, density, unit, cutoff
13 => Optimise FIELD output size, if possible? 1=yes  0=no
14 => Atom display: 1 = DL_FIELD format. 2 = Standard format
15 => Vdw display format: 1 = 12-6 format   2 = LJ format
16 => Epsilon mixing rule (organic FF only) : default, or 1 = geometric, 2 = arithmatic 
17 => Sigma mixing rule (organic FF only) : default, or 1 = geometric, 2 = arithmatic 
18 => Epsilon mixing rule (inorganic FF only) : 1 = geometric   2 = arithmatic
19 => Sigma mixing rule (inorganic FF only) : 1 = geometric   2 = arithmatic 
20 => Epsilon mixing rule (BETWEEN different FF) : 1 = geometric   2 = arithmatic
21 => Sigma mixing rule (BETWEEN different FF): 1 = geometric 2 = arithmatic
22 => Display additional info. for protein 1=Yes  0=No
23 => Freeze atoms? 1 = Yes (see below)  0 = No
24 => Tether atoms? 1 = Yes (see below)  0 = No
25 => Constrain bonds? 1 = Yes (see below) 0 = No
26 => Apply rigid body? 1 = Yes (see below) 0 = No
27 => Periodic condition ? 0=no, other number = type of box (see below)
28 =>    22.1544000000        0.0000000000        0.0000000000* Cell vector a (x, y, z) 
29 =>    0.0000000000       17.3213000000        0.0000000000* Cell vector b (x, y, z) 
30 =>    0.0000000000        0.0000000000       23.8724000000* Cell vector c (x, y, z)    
31 => 1-4 scaling for coulombic (put default or x for scaling=x)
32 => 1-4 scaling for vdw (put default or x for scaling=x)
33 => Include velocity? 1=yes, 0=no and scaling temperature.
34 => Position solute at origin? 1 = yes, 0=no 
35 => Solvate model? none or specify solvent (see below) and distance criteria.
36 => Add counter ions? 1=yes, 0=no,  minimum distance from solute  
37 => MM energy calculation. 1=Yes, 0=No
38 => Cut off for electrostatic energy calculation (angstrom)
39 => Cut off for vdw energy calculation (angstrom) 


# Water
TIP3P_O - Original TIP3P - water_tip3p_o
TIP3P_C - Charmm-specific TIP3P model - water_tip3p_c
TIP3P_E - TIP3P water model optimised for ewald sum - water_tip3p_e
TIP4P   - TIP4P water model (one pseudo point) - water_tip4p
TIP5P   - Five-site water model (two pseudo points) - water_tip5p
TIP5P_E - Five site water model optimised for ewald sum - water_tip5p_e
SPC     - SPC water model - water_spc
SPCE    - SPC/E water model - water_spce
WATER_PCFF - Water model specific to PCFF - water_pcff
# Alcohol
MeOH    - Methanol, obtained from CHARMM FF - methanol
EtOH    - Ethanol, from OPLS2005 - ethanol
IPA     - Isopropyl alcohol (propan-2-ol), from CHARMM FF - ipa
# Misc.
hexane  - linear hexane, from CHARMM36 CGenFF - hexane
THF     - tetrahydrofuran, from CHARMM FF - thf 
EA      - ethyl acetate, EtOAc, from OPLS2005 - ethylacetate
CCL4    - carbon tetrachloride, CCl4, from OPLS2005 - ccl4
HCPE    - hexachloropropene, from OPLS2005 - hcpe
acetone - propanone, acetone, OPLS2005 - acetone
CH3CN   - acetonitrile, ethanenitrile, OPLS2005 - acetonitrile
dmso    - dimethylsulphoxide, OPLS2005 - dmso
toluene - methylbenzene, from OPLS2005 - toluene

"""

class FieldTools:
    # installation directory for dl_field, MUST change to appropriate directory.
    dlf_dir = "/home/jay/dl_field/dl_f_4.5"
    metaut = MetaUtil()

    def __init__(self):
        pass
      
    def copyTofield (self, filename="1mol.xyz", directory=os.curdir):
        
        """
        Copies target file to DL_FIELD folder
        
        Arguments:
            filename (str) ==> name of file to copy to DL_FIELD.
            directory (str) (optional) ==> defaults to copying file from current directory.
            
        Returns:
            none.
            
        """
        
        loc = os.path.join(directory,filename)
        dlf_loc = os.path.join(self.dlf_dir,filename)
        shutil.copy(loc,dlf_loc)
        
    
    def edit_dlfdir (self, ndir):
        
        """
        DL_FIELD installation directory modifier.
        
        Arguments:
            ndir (str,pathlike) ==> new DL_FIELD directory
            
        Returns:
            none.
    
        """
        self.dlf_dir = ndir
    
    def run_dlf(self):
        
        """
        Runs DL_FIELD, stdout and stderr both suppressed.
        
        Arguments:
            none.
            
        Returns:
            none
        """
        
        FNULL = open(os.devnull,'w')
        cur_dir = os.getcwd()
        os.chdir(self.dlf_dir)
        subprocess.run("./dl_field", stdout=FNULL,stderr = subprocess.STDOUT)
        os.chdir(cur_dir)
        
    def edit_dlfc(self, n, nw, filename="dl_field.control"):
        """
        Edits DL_FIELD control file when given a numerical key and the new value to
        change to. Numerical keys and their corresponding variables are given in the 
        code header. Previous control file is backed up as dl_field_backup.control.
        
        Arguments:
            n (int) ==> Numerical key for a variable, a list is shown in the header.
            nw (str) ==> value to change to, only accepts txt like strings.
            filename (str,pathlike) (optional) ==> defaults to changing dl_field.control file.
            
        Returns:
            none.
        
        """
        
        cur_dir = os.getcwd()
        os.chdir(self.dlf_dir)
        cont = []
        count = 0
        with open (filename,'r') as fh:
            for line in fh:
                if line.isspace():
                    continue
                else:
                    if count == n:
                        words = line.split("*")
                        old_w = words[0]
                        n_line = line.replace(old_w,nw+" ")
                        cont.append(n_line)
                    else:
                        cont.append(line)
                    count += 1
                    
        os.rename("dl_field.control","dl_field_backup.control")
        
        with open (filename,'w') as fw:
            for data in cont:
                fw.write(data)
                
        os.chdir(cur_dir)
    
    def get_output(self,kw):
        
        """
        Grabs outputs from DL_FIELD output folder.
        Keywords ==> "all" signals all files, however this keyword does not remove
        the dl_field prefix in the filenames. "field" signals FIELD file and "CONFIG" signals
        CONTROL file.
        
        Arguments:
            kw (str) ==> signals which file(s) to grab from the output folder.
            
        Returns:
            none.
        """
        
        out_dir = os.path.join(self.dlf_dir,"output")
        tar_dir = os.curdir
        
        try:
            os.mkdir(tar_dir)
        except FileExistsError:
            pass
        
        if kw == "all":
            for i in os.listdir(out_dir):
                if i != 'readme.txt':
                    shutil.copy(os.path.join(out_dir,i),os.path.join(tar_dir,i))
                    
        elif kw == "field":
            shutil.copy(os.path.join(out_dir,"dl_poly.FIELD"),os.path.join(tar_dir,"FIELD"))
            
        elif kw == "config":
            shutil.copy(os.path.join(out_dir,"dl_poly.CONFIG"),os.path.join(tar_dir,"CONFIG"))
        
        else:
            print("target keyword not recognised")
            
    
    def reOrder_config(self,f_block,directory=os.curdir,filename="CONFIG"):
        
        """
        Reorders Config file according to the atom name order in the first molecule. First molecule
        data block obtained from 1mol.xyz is inputed.
        
        Arguments:
            f_block(list(str)) ==> first molecule data block, entered as raw txt format.
            directory (str,pathlike) ==> directory of the config file to reorder.
            filename (str) ==> file name of the config file, defaults to "CONFIG".
            
        Returns:
            none
        """
        data = []
        numatm = len(f_block)
        atoms = []
        
        for i in f_block:
            atoms.append(i.split()[0])
            
        atomnames = list(set(atoms)) # Create unique atom names
        
        atom_indices = {}
        
        # Get indices of each atoms
        for atomname in atomnames:
            atom_indices[atomname] = []
            for i, j in enumerate(f_block):
                if atomname == j.split()[0]:
                    atom_indices[atomname].append(i)
        print(atom_indices)
                    
        block = []
        molecule = []
        count = 0
        
        with open(filename,'r') as fh:
            
            for n in range(5):
                data.append(fh.readline())
            
            for line in fh:
                block.append(line)
                count += 1
                
                if count % 2 == 0:
                    molecule.append(block)
                    block = []
                    
                    if count % (2*numatm) == 0:
                        move_count = 0
                        new_molecule = list(np.zeros(numatm))
                        for atomname in atom_indices.keys():
                            for i in molecule:
                                if i[0].split()[0]==atomname:
                                    new_line = i
                                    numbering = new_line[0].split()[1]
                                    new_numbering = str(int(count/2 - (numatm- atom_indices[atomname][move_count]-1)))
                                    new_line[0] = new_line[0].replace(numbering, new_numbering)
                                    new_molecule[atom_indices[atomname][move_count]] = new_line
                                    move_count += 1
                            move_count = 0    
                        data.extend(new_molecule)
                        molecule = []
                        
        with open('reordered.config','w') as fw:
            for i in data:
                print (i)
                if len(i)>2 :
                    fw.write(i)
                else:
                    for j in i:
                        for k in j:
                            fw.write(k)


    def get_field(self,nummol,directory=os.curdir,filename="FIELD"):
        
        """
        Redefine multimols FIELD file as singmol FIELD file (adds the number of
        molecules to the field file).
        
        Arguments:
            nummol (int) ==> number of molecules in multimols FIELD file.
            directory (str,pathlike) ==> directory of FILED file, defaults to CWD.
            filename (str) ==> defaults to filename "FIELD".
            
        Returns:
            none.
        """
        
        cont = []
        block = []
        
        with open(filename,'r') as fh:
            
            for n in range(10000):
                temp = fh.readline()
                
                if "nummols" in temp:
                    temp = temp.replace('1',str(nummol))
                cont.append(temp)
                
                if "atoms" in temp:
                    numatm = int(temp.split()[1])
                    break
            
            for n in range(numatm):
                block.append(fh.readline())
                
            cont.extend(block)
            
            for line in fh:
                cont.append(line)
                
        with open("FIELD",'w') as fw:
            for i in cont:
                fw.write(i)
        
    def update_vector(self, directory = os.curdir, filename ="mol_o0001no_res.out"):
        
        """
        Grabs cell vector from molres.out file and updates the DL_FILED control file
        with the obtained cell_vector (includes the determination and update of periodic 
        BC type).
        
        Arguments:
            directory (str,pathlike) (optional) ==> directory of the molres.out file.
        
        Returns:
            vector_a (arraylike) ==> 3x3 cell vector
        """
        vector = []
        test_ph = []
        
        try:
            open(filename, 'r')
        except FileNotFoundError:
            print ("unable to locate source (mol_res) file")
            return
    
        with open(filename, 'r') as fh:
            check = fh.readline()

            if check == "latt\n":
                for i in range(3):
                    vector.append(fh.readline())
            else:
                print ("incorrect mol file")
                return
        
        if len(vector) == 3:
            for i in range(3):
                test = vector[i].split()
                test_ph.append(test)
                vector_a = np.asarray(test_ph,dtype=float)
        else:
            print ("vector incomplete, check data in mol.out file")
            return
        
        self.change_vector(vector_a)
        return vector_a
    
    def box_type(self,vector_a):
        
        """
        Determines the periodic boundary condition of a given cell vector.
        
        Arguments:
            vector_a (arraylike) ==> input cell vector.
        
        Returns:
            per_c （int) ==> periodic boundary condition type, 1:cubic, 2:orthorhombic
            3:paralellopiped.
        """
        
        if (vector_a/np.diag(vector_a)==np.identity(3)).all():
            diagonals = list(np.diag(vector_a))
            if diagonals.count(diagonals[0]) == len(vector_a):
                per_c = 1
            else:
                per_c = 2
        else:
            per_c = 3
            
        return per_c
    
    def change_vector (self,vector_a):
        
        """
        Update DL_FIELD control file with a given input cell vector.
        
        Arguments:
            vector_a (arraylike) ==> input cell vector.
        
        Returns:
            none.
        """
        vector_a = self.remake_cell(vector_a)
        
        for i in range(3):
            vec_string = "  {}  {}  {} ".format(*vector_a[i])
            self.edit_dlfc(i+28,vec_string)
            
        per_c = self.box_type(vector_a)
        self.edit_dlfc(27,str(per_c))

    
    def meta_sym_coords(self, sym_n, cif_ratio= 1, mol_file = "af_co0001no.xyz"):

        print ("Calculating sym coords...")
        
        grow_holder = self.metaut.edit_grow("1 1 1")
        
        FNULL =open(os.devnull,'w')
        subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT,shell=True)
        subprocess.run("metadise",stdout = FNULL, stderr = subprocess.STDOUT)
        cell_vector = self.update_vector()
        
        count = 0
        coords = {}
        onemol = []
        titles = []
        block = []
        
        
        with open (mol_file,'r') as fh:
            
            titles.append(fh.readline())
            vector_para = fh.readline()
            titles.append("CRYST1   {}    {}    {}    {}    {}    {}\n".format(*vector_para.split()[-6:]))

            for line in fh:
                
                coord = list(map(float,line.split()[1:4]))
                block.append (coord)
                
                count += 1
                
                if count != 0 and count % sym_n == 0:
                    coords[int((count/sym_n)-1)] = block
                    block = []
                if count == 0 or count % sym_n == 1 or sym_n == 1:
                    onemol.append(line)
        if cif_ratio != 1:
            print ("CIF cut detected, please make sure cif atom coords are seprerated molecule-wise \
                   and are ordered correctly")
            
            newatmnum = int(len(onemol)/cif_ratio)
            onemol = onemol [0:newatmnum]
            keys = list(coords.keys())
            for i in keys:
                if i > newatmnum-1:
                    coords[i- (cif_ratio - 1)*newatmnum].extend(coords[i])
                    del coords[i]
                    
        atmnum = len(onemol)
        nummol = int(sym_n * np.prod(list(map(int,grow_holder.split())))*cif_ratio)
        print ("atom number  =  {}".format(atmnum))
        print ("molecule number  = {}".format(nummol))

        self.metaut.edit_grow (grow_holder)
        subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT,shell=True)
        subprocess.run("metadise",stdout = FNULL, stderr = subprocess.STDOUT)
        return atmnum, nummol, coords, onemol, cell_vector
        
    
    def reOrder_xyz (self,onemol,coords,cell_vector,filename="mol_o0001no.xyz"):
        
        """
        Reorders xyz file according to initial mol coords and its symmetry coords. This 
        would provide the correct order for DL_POLY to read off of the singmol FIELD file.
        
        Arguments:
            onemol (list(str)) ==> initial molecule xyz as raw txt string block.
            cell_vector (arraylike) ==> cell vector of grown system.
            coords (dict) ==> symmetry coordinates of each atom in the molecule.
        
        Returns:
            none.
        """
        
        print ("Reordering xyz...")
        data = []
        numatm = len(onemol)
        atoms = []
        print(onemol)
        for i in onemol:
            atoms.append(i.split()[0])
            
        atomnames = list(set(atoms))
        
        atom_indices = {}
        
        for atomname in atomnames:
            atom_indices[atomname] = []
            for i, j in enumerate(onemol):
                if atomname == j.split()[0]:
                    atom_indices[atomname].append(i)                    
        block = []
        count = 0
        
        with open(filename,'r') as fh:
            
            for n in range(2):
                data.append(fh.readline())
                        
            for line in fh:
                print(line)
                block.append(line)
                count += 1
                
                if count % numatm == 0:
                    new_block = list(np.zeros(numatm))
                    for atom in atom_indices.keys():
                        for indices in atom_indices[atom]:
                            poss_coords = coords[indices]
                            for atom_ in block:
                                atom_coord = list(np.array(atom_.split()[1:4],dtype=float))
                                if atom == atom_.split()[0]:
                                    for poss_coord in poss_coords:
                                        if self.check_coordinate(atom_coord,poss_coord,cell_vector):
                                            new_block[indices]=atom_
                                            break
                    data.extend(new_block)
                    block = []

        with open('target.xyz','w') as fw:
            for i in data:
                fw.write(i)
            
    def check_coordinate(self,coordone,coordtwo,cell_vector):
        
        """
        Checks if one set of coordinates is a scalar growth of another given the 
        scalar growth vector.
        
        Arguments:
            coordone (array(float)) ==> first set of coordinates.
            coordtwo (array(float)) ==> second set of coordinates.
        
        Returns:
            Boolean result for the check.
        """
        
        rel_tol = 0.001 # Relative tolerance to account for floating point uncertainties.

        cgrow = (coordone[2]-coordtwo[2])/cell_vector[2,2]
        bgrow = (coordone[1]-cgrow * cell_vector[2,1] - coordtwo[1])/cell_vector[1,1]
        agrow = (coordone[0]-bgrow * cell_vector[1,0] - cgrow * cell_vector[2,0] - coordtwo[0])/cell_vector[0,0]
        
            
        
        #if coordone[0] == 7.1081508874:
            #print (cgrow)
            #print (bgrow)
            #print (agrow)
            #print((abs(cgrow) % 1 + abs(bgrow) % 1 + abs(agrow) % 1))
        if (round(cgrow,3) % 1 + round(bgrow,3) % 1 + round(agrow,3) % 1) < rel_tol:
            retrace = np.asarray(coordone) - cgrow * cell_vector[2] - bgrow * cell_vector[1] - agrow * cell_vector[0]
            if sum(abs(retrace-coordtwo)) < rel_tol:
                return True
        return False
    
    def get_dlf(self, cif_ratio=1, cleanmeta=True):
        
        """
        Main logic for obtaining singmol FIELD file and reorder CONFIG file in order
        for singmol FIELD file to work
        
        Arguements:
            numatm (int) ==> number of atoms in one target molecule.
            nummol (int) ==> number of mols in the singmol FIELD file.
        """
        
        #logic
        #get 1 mol xyz
        #run on dl_field
        #get output field file
        #reorder xyz
        #run dl_f
        #get output config file
        FNULL =open(os.devnull,'w')
        
        print ("Trying to construct model...")
        cf = CifTools("input.txt")
        symn = cf.get_symn()
        atmnum, nummol, coords, onemol, cell_vector = self.meta_sym_coords(symn,cif_ratio = cif_ratio)
        self.reOrder_xyz(onemol,coords,cell_vector)
        self.get_1mol_xyz(atmnum)
        self.copyTofield()
        self.edit_dlfc(10,"1mol.xyz")
        self.change_vector(cell_vector)
        self.run_dlf()
        self.get_output("field")
        self.get_field(nummol)
        self.copyTofield(filename="target.xyz")
        self.edit_dlfc(10,"target.xyz")
        self.update_vector()
        self.run_dlf()
        self.get_output("config")
        print ("Model generated successfully.")
        if cleanmeta:
            subprocess.run("rm *000*",stdout=FNULL, stderr = subprocess.STDOUT,shell=True)
            print ("Residual files cleared.")
            
    def get_1mol_xyz(self,numatm,directory=os.curdir,filename="target.xyz"):
        
        """
        Cut mol.xyz from to one molecule only. writes it in 1mol.xyz file and outputs
        the one molecule xyz data as a string block.
        
        Auguments:
            numatm (int) ==> number of atoms in the first molecule.
            directory (str, pathlike) ==> directory of the mol.xyz file.
            filename (str) ==> file name of the mol.xyz file.
            
        Returns:
            block (str) ==> a raw txt string list for the one molecule xyz file.
        """
        cont = []
        block = []
        count = 0
        
        with open(filename,'r') as fh:  
            
            fh.readline()
            cont.append(("      "+str(numatm)+'\n'))
            cont.append(fh.readline())
            
            for line in fh:
                block.append(line)
                count += 1
                if count == numatm:
                    cont.extend(block)
                    break

        with open('1mol.xyz','w') as fw:
            for i in cont:
                fw.write(i)
        
        return block
       
    def get_nosort(self, vector_a,  file = "target.xyz"):
        
        """
        Main logic for none reordered DL_FILED operations.
        """
        self.change_vector (vector_a)
        self.copyTofield(filename=file)
        self.edit_dlfc(10,file)
        self.run_dlf()
        self.get_output("all")
    
    def remake_cell(self, c_vector):
        """        
        Cell vector remake for cells with too much angle (which may cause 
        DL_FIELD to terminate)
        
        Parameters
        ----------
        c_vector : Array
            cell vector array

        Returns
        -------
        c_vector : Array
            Modified cell vector array

        """
        
        of_diag =  abs(c_vector[0,1]) + abs(c_vector[0,2])
        th = 5
        if of_diag/abs(c_vector[0,0]) > th:
            c_vector = self.box_optimiser(c_vector)
            return c_vector
        else:
            return c_vector
    
    def box_optimiser (self, c_vector, convergence = "loose"):
        
        if convergence == "loose":
            quo = c_vector[0,2] // c_vector[2,2]
            c_vector[0,2] %= c_vector[2,2]
            c_vector[0,1] -= quo * c_vector[2,1]
            c_vector[0,1] %= c_vector[1,1]
        
        return c_vector                    
        
    def prep_xyz (self,file):
        
        """
        Prepares metadise default output xyz file format to work for DL_FIELD (applying CRYST1 statement).
        
        Arguments:
            file (str) ==> xyz file name.
            
        Returns:
            none.
        """
        
        col = []
        
        with open (file, 'r') as fh:
            col.append(fh.readline())
            temp = fh.readline().split()
            new = " CRYST1 {} {} {} {} {} {}\n".format(temp[2],temp[3],temp[4],temp[5],temp[6],temp[7])
            col.append(new)
            
            for line in fh:
                col.append(line)
        
        with open (file,'w') as fw:
            for i in col:
                fw.write(i)
    
    def cut_1m_field(self, nummol, file = "FIELD"):
        """
        must be ordered by molecules.
        """
        cont = []
        itera = ["atoms","bonds","angles","dihedral"]
        with open (file, 'r') as fh:
            
            while True:
                line = fh.readline()
                
                if any(x in line for x in itera):
                    nf = float(line.split()[1])/nummol
                    
                    if not nf.is_integer():
                        return 0
                    
                    line = line.replace(line.split()[1],str(int(nf)))
                    cont.append(line)
                    for i in range(int(nf)):
                        inf = fh.readline()
                        cont.append(inf)
                    for i in range (int(nf)*(nummol-1)):
                        fh.readline()
                else:
                    if "nummols" in line:
                        line = "nummols {}\n".format(nummol)
                    cont.append(line)

                if line.split()[0] == "close":
                    break
        with open("FIELD", 'w') as fw:
            
            for i in cont:
                fw.write(i)
                        




























