import numpy as np
import subprocess
from .MetaUtil import MetaUtil
"""
Input file converter from DLPOLY to DLMONTE.
Includes a few field utilities related to other programmes.
"""

class InputConverter:
    
    def __init__(self, ffile = "FIELD", cfile = "CONFIG", contfile = "CONTROL"):
        self.ffile = ffile
        self.cfile = cfile
        self.contfile = contfile
    
    def ffetch(self, keyword,):
        
        """
        Fetches data blocks from FIELD file directly, these data blocks includes atoms,
        bonds,angles etc...
        
        Arguments:
            keyword (str) ==> the data block to obtain.
            
        Returns:
            blocktile (str) ==> the title that intiates the block.
            block (list(str)) ==> data block (split no space).
        """

        block = []
        blocktitle = ""
        
        with open (self.ffile,'r') as fh:
            for line in fh:
                if line.split()[0] == keyword:
                    blocktitle = line
                    break
                
            for i in range(int(line.split()[1])):
                line = fh.readline()
                block.append(line.split())
    
        return blocktitle,block
    
    def bfetch(self, keyword, data):

        """
        Fetches data block from a raw txt string block of the entire FIELD file.
        
        Arguments:
            keyword (str) ==> the data block to obtain.
            data (list(str)) ==> raw txt string data.
            
        Returns:
            blocktile (str) ==> the title that intiates the block.
            block (list(str)) ==> data block (split no space).
        
        """
        
        block = []
        blocktitle = ""
        start = 0
        end = 0
        
        for linenum in range(len(data)):            
            if data[linenum].split()[0] == keyword:
                blocktitle = data[linenum]
                start = linenum + 1
        end = int(blocktitle.split()[1]) + start
        
        tempblock = data[start:end]
        
        block = [i.split() for i in tempblock]
            
        return blocktitle, block

    def get_atoms(self,data):
        
        """
        Get Atoms block and edit it to the format of DLMONTE input.
        
        Argument:
            data (list(str)) ==> raw txt string data for the entire FIELD file.
        
        Returns:
            atoms (list(str)) ==> atoms block for DLMONTE input.
        """
        
        title,block = self.bfetch("atoms",data)
        charges = [i[2] for i in block]
        atom_names = [i[0] for i in block]
        atoms = np.array(block)
        uniq = set(list(atoms[:,0]))
        atoms= []
        for i in uniq:
            for j in block:
                if i in j:
                    atoms.append(j)
                    break
        atoms = [i[0:3] for i in atoms]
        atoms = [np.insert(i,1,"CORE") for i in atoms]
        atoms = self.spacecreator(atoms)
        
        return atoms, charges, atom_names
    
    def get_pos(self,atmnum,atmnames, uc=[]):
        
        """
        Get atom positions block and edit it to the format of DLMONTE input.
        
        Argument:
            atmnum (int) ==> the number of atom coordinates to grab from CONFIG file
        
        Returns:
            pos (list(str)) ==> atom position block.
        """
        pos = []
        
        with open (self.cfile,'r') as fh:
            count = 0
            for i in range(10000000):
                temp = fh.readline()
                
                if not temp:
                    continue
                else:
                    count += 1
                
                if count == 5:
                    count = 0
                    break
            
            
            while True:

                name = fh.readline().split()[0]
                coords = fh.readline().split()
                if name == atmnames[count]:
                    pos.append([name]+coords)
                    count += 1
                else:
                    pos = []
                    count = 0
                
                if count == atmnum:
                    break
                

            #for line in fh:
            #    
            #    if count % 2 == 0:
            #        
            #        if count / (2*atmnum) == 1:
            #            break
            #        else:
            #            block.append(line.split()[0])
            #        
            #    else:
            #        block.extend(line.split())
            #        pos.append(block)
            #        block = []
            #    
            #    if not line == False:
            #        count += 1


            for i in pos:
                i.insert(1,"CORE")
            
            if len(uc) > 0:
                if len(uc) == len(pos):
                    for i in range(len(uc)):
                        pos[i][-1] = pos[i][-1] + " charge " + uc[i]
                else:
                    print(len(uc),len(pos))
                    print("unique charges profile doesn't match the molecule.")
            pos = self.spacecreator(pos)                
        return pos
    
    def get_vdw(self, lr_label = True, low_e_filter = True):
        
        """
        Get vdw force field block and convert into DLMONTE format.
        
        Arguments:
            none
        Returns:
            vdw (list(str)) ==> vdw data block.
        """
        
        title,block = self.ffetch("vdw")
        
        vdw = []
        vdw.append(title)
        for i in block:
            temp = ""
            for j in range(len(i)):
                temp += i[j] + "    "
                if j < 2:
                    temp += "CORE    "
                if lr_label and j == 1:
                    temp += "+"
                if j == len(i) -1 :
                    temp += "\n"
            if low_e_filter and all([float(i[-2])<1e-5, float(i[-1])<1e-4]):
                continue
            else:
                vdw.append(temp)
        vdw[0] = "vdw {}\n".format(len(vdw)-1)
        return vdw
    
    def sort_block(self, indend, block, ins):
        
        """
        Sort the interation force field blocks into DLMONTE format, which consists of
        an bonding block and a force field constants block
        
        Arguments:
            indend (int) ==> number of columns to grab to obtain bonding info.
            block (list(str)) ==> data block to sort
            ins (int) ==> starting indices to start the numbering of the bonding interactions.
        
        Returns:
            targets (list(str)) ==> target interaction block.
            targettypes (list(str)) ==> target interaction block containing the constants.
        """
        
        constants = []
        targets = []
        types = []
        targettypes = []
        
        for i in block:
            if i[indend:] not in constants:
                constants.append(i[indend:])
                types.append(i[0])
            temp = i[1:indend]
            temp.append(str(constants.index(i[indend:])+1+ins))
            targets.append(temp)
        
        for i in range(len(constants)):
            temp = constants[i]
            temp.insert(0,types[i])
            targettypes.append(temp)
        
        return targets,targettypes
    
    def spacecreator (self,block):
        """
        Generates spaces and next line from a block of space splitted string lists, to form a raw txt data format.
        
        Arguments:
            block （list(str)) ==> data to create space for.
    
        Returns:
            lines (list(str)) ==> raw txt format data.
        """
        lines = []
        
        for i in block:
            line = ""
            for j in range(len(i)):
                if j +1 == len(i):
                    line += i[j] + "\n"
                else:
                    line += i[j] + "    "
            lines.append(line)
            
        return lines
    
    
    """
    Methods for obtaining DLMONTE Dataformat for bonds,angles and dihedrals, theses
    methods are analogous to eachother.
    """
    
    def get_bonds(self, moldata, ins):
        
        title, block = self.bfetch ("bonds", moldata)
        
        bonds,bondtypes = self.sort_block(3,block,ins)
        
        bondblock = self.spacecreator(bonds)
        bondtypeblock = self.spacecreator(bondtypes)
        
        return title,bondblock,bondtypeblock

    def get_angles(self, moldata, ins):
        
        title, block = self.bfetch("angles", moldata)
        
        angles, angletypes = self.sort_block(4,block,ins)
        
        angleblock = self.spacecreator(angles)
        angletypeblock = self.spacecreator(angletypes)
        
        return title,angleblock, angletypeblock
    
    def get_dihedrals(self, moldata, ins):
        
        title, block = self.bfetch("dihedral", moldata)
        
        dihedrals, dihedraltypes = self.sort_block(5,block,ins)
        
        dihedralblock = self.spacecreator(dihedrals)
        dihedraltypeblock = self.spacecreator(dihedraltypes)
        
        return title, dihedralblock, dihedraltypeblock
    
    
    def get_info (self):
        
        """
        Get extra none block info from the FIELD file and converts them into DLM format.
        Fetch raw txt block info from DLFIELD file.
        
        Arguments:
            field (str) ==> default FIELD file name
            control (str) ==> default CONTROL file name
            
        Returns:
            misl (list(str)) ==> contains title, cutoff, unit and NCONFIG info.
            Molecular types (str) ==> contains the number of molecular types info.
            MolInfo(dict) ==> raw txt block info for all molecular types.
        """
        
        cutoff = "CUTOFF    "
        unit = "UNITS    "
        title = ""
        NCONFIG = "NCONFIG    1\n"
        Moleculartypes = "MOLTYPES    "
        molnum = 0
        MolInfo = {}
        
        with open(self.contfile,'r') as fh:
            for i in range(10000000):
                temp = fh.readline()
                if temp.split()[0] == "cutoff":
                    cutoff += temp.split()[1] + "\n"
                    break
                
        with open(self.ffile,'r') as fh:
            
            for i in range(10000000):
                temp = fh.readline()
                if i == 0:
                    title = temp
                if temp.split()[0] == "Units":
                    unit += temp.split()[1] + "\n"
                if temp.split()[0] == "Molecular" and temp.split()[1] == "types":
                    Moleculartypes += temp.split()[2] + "\n"
                    molnum = int(temp.split()[2])
                    break
                
            for i in range(molnum):
                molcont = []
                molname = ""
                for linenum in range(10000000):
                    temp = fh.readline()
                    if not temp:
                        continue
                    elif "Molecule name" in temp:
                        molname = temp.split()[2]
                    elif "finish" in temp:
                        break
                    else:
                        molcont.append(temp)
                MolInfo[molname] = molcont
        
        misl = [title,cutoff,unit,NCONFIG]
        
        return misl, Moleculartypes, MolInfo
    
    def get_atomnum (self,data):
        
        """
        Calculates the number of atoms per molecule from the data block.
        
        Arguments:
            data (list(str)) ==> raw txt data block.
            
        Returns:
            atomnum (int) ==> number of atoms per molecule.
        """
        
        title,block = self.bfetch("atoms",data)
        atomnum = int(title.split()[1])
        return atomnum
        
    def write_field(self, uc= True, lr_label = True, low_e_filter = True):
        
        """
        Method to collect all necessary info and write the FIELD file in DLMONTE format.
        
        """
        filecont = []
        AtomInfo = []
        BondInfo = []
        AngleInfo = []
        DihedralInfo = []
        uc_info = {}
        names_info = {}

        misl, Moleculartypes, MolInfo = self.get_info() # Grabs info
        
        filecont.extend(misl) 
        
        for i in MolInfo.keys():
            inf,uc,names = self.get_atoms(MolInfo[i])
            AtomInfo.extend(inf) # Concat Atom info if multiple mol types
            uc_info[i] = uc
            names_info[i] = names
        AtomTitle = "ATOMS    " + str(len(AtomInfo)) + "\n"
        filecont.append(AtomTitle)
        filecont.extend(AtomInfo)
        
        print(names_info)

        filecont.append(Moleculartypes)
        
        for i in MolInfo.keys():
            filecont.append(i+"\n")
            atomnum = self.get_atomnum(MolInfo[i])
            if uc:
                pos = self.get_pos(atomnum,names_info[i],uc=uc_info[i])
            else:
                pos = self.get_pos(atomnum,names_info[i])
            filecont.append("ATOMS    " + str(len(pos)) + "    " + str(len(pos)) + "\n")
            filecont.extend(pos)
            if any(['rigid' in j for j in MolInfo[i]]):
                filecont.append("RIGID\n")
                filecont.append("EXCLUDE\n")
            else:
                title,bonds,bondtypes = self.get_bonds(MolInfo[i],len(BondInfo))
                BondInfo.extend(bondtypes)
                filecont.append(title)
                filecont.extend(bonds)
                title,angles,angletypes = self.get_angles(MolInfo[i],len(AngleInfo))
                AngleInfo.extend(angletypes)
                filecont.append(title)
                filecont.extend(angles)
                try:
                    title,dihedrals,dihedraltypes = self.get_dihedrals(MolInfo[i],len(DihedralInfo))
                    DihedralInfo.extend(dihedraltypes)
                    filecont.append(title)
                    filecont.extend(dihedrals) 
                except IndexError:
                    print("No dihedrals found for {}.".format(i))
            filecont.append("FINISH \n")
            
        filecont.append("BOND    "+ str(len(BondInfo)) + "\n")
        filecont.extend(BondInfo)
        filecont.append("ANGLE    "+ str(len(AngleInfo)) + "\n")
        filecont.extend(AngleInfo)
        filecont.append("DIHED    "+ str(len(DihedralInfo)) + "\n")
        filecont.extend(DihedralInfo)
        
        vdw = self.get_vdw(lr_label=lr_label,low_e_filter=low_e_filter)
        filecont.extend(vdw)
        filecont.append("CLOSE")
        
        with open("FIELD_DLM",'w') as fw:
            for i in filecont:
                fw.write(i)
    
    def get_nummol(self):
        
        """
        Get necessary molecular info from FIELD file.
        
        Returns:
            molnum (list(str)) ==> molecular info containing molecular name, number of molecules
            and number of atoms in the molecule.
        """
        
        with open(self.ffile, 'r') as fh:
            molname = []
            nummol = []
            molatm = []
            for line in fh:
                if "Molecule name" in line:
                    molname.append(line.split()[2])
                elif "nummols" in line:
                    nummol.append(line.split()[1])
                elif "atoms" in line:
                    molatm.append(line.split()[1])
                else:
                    continue
        molnum = list(zip(molname,nummol,molatm))
        
        return molnum
    
    def unique_atomise (self, output_filename = "FIELD_DLM_UA", molecule_name = "MEOH"):

        unique_temp = []
        atom_names = []
        with open (self.ffile, 'r') as fh:
            for line in fh:

                if molecule_name in line:
                    x = fh.readline()
                    x = int(x.split()[1])
                    for i in range(x):
                        atom_line = fh.readline()
                        st_data = atom_line.split()[0] + " " + atom_line.split()[-1]
                        if any([st_data == i for i in unique_temp]):
                            continue
                        else:
                            unique_temp.append(st_data)

        
        cont = []

        with open (self.ffile,'r') as fh:
            flag = 0
            atoms_cont = []
            st_data_cont = []
            for line in fh:

                if "MOLTYPES" in line:
                    flag = 1
                    cont.append(line)

                elif ("ATOMS" in line and flag == 0):
                    atomnum = int(line.split()[1])
                    cont.append("ATOMS {}\n".format(atomnum))

                    for atom_c in range(atomnum):
                        atom_line = fh.readline()
                        atoms_cont.append(atom_line)
                        st_data_cont.append(atom_line.split()[0] +" "+ atom_line.split()[-1])
                    
                    other_target = [i.split()[0] for i in st_data_cont if i not in unique_temp]
                    print(other_target)
                    for i in np.arange(len(atoms_cont)):
                        st_data = st_data_cont[i]
                        atom_line = atoms_cont[i]
                        if st_data in unique_temp and st_data.split()[0] in other_target:
                            atom_list = atom_line.split()
                            atom_names.append(atom_list[0])
                            atom_list[0] = atom_list[0] + "_" + molecule_name[0]
                            cont.append("    ".join(atom_list) + "\n")
                        else:
                            cont.append(atom_line)


                elif molecule_name in line:
                    cont.append(line)
                    atom_info = fh.readline()
                    cont.append(atom_info)
                    atom_info = atom_info.split()
                    for atoms in range(int(atom_info[1])):
                        atom_line = fh.readline()
                        atom_name = atom_line.split()[0]
                        atom_list = atom_line.split()
                        if atom_list[0] in atom_names:
                            atom_list[0] = atom_name + "_" + molecule_name[0]
                            atom_line = "    ".join(atom_list) + "\n"
                        cont.append(atom_line)

                elif "vdw" in line:
                    count = 0
                    vdw_info = line
                    cont.append(line)
                    vdw_line_num = int(line.split()[1])
                    for vdws in range (vdw_line_num):
                        vdw_line = fh.readline()
                        vdw_list = vdw_line.split()
                        vdw_names = [vdw_list[0],vdw_list[2]]
                        dup_cont = [x in atom_names for x in vdw_names]
                        print(dup_cont)
                        if all(dup_cont):
                            cont.append(vdw_line)
                            if vdw_names[0] == vdw_names[1]:
                                vdw_list[0] = vdw_names[0] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                                vdw_list[2] = vdw_names[0] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                                count += 2
                            else:
                                vdw_list[0] = vdw_names[0] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                                vdw_list[2] = vdw_names[1] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                                vdw_list[0] = vdw_names[0]
                                cont.append("    ".join(vdw_list) + "\n")
                                count += 3
                        elif any(dup_cont):
                            cont.append(vdw_line)
                            if dup_cont[0]:
                                vdw_list[0] = vdw_names[0] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                            else:
                                vdw_list[2] = vdw_names[1] + "_" + molecule_name[0]
                                cont.append("    ".join(vdw_list) + "\n")
                            count += 1
                        else:
                            cont.append(vdw_line)

                else:
                    cont.append(line)
        
        cont[cont.index(vdw_info)] = "vdw {}\n".format(count + int(vdw_info.split()[1]))

        with open (output_filename, 'w') as fw:
            for i in cont:
                fw.write(i)
        
        return atom_names
    
    def lambda_field (self, molecule_name = "XYZ", output_filename="FIELD_lambda", coeff = "0.5 1 2 6", lpath = [0.0,0.8,0.8,1.0]):

        cont = []
        mol_cont = []

        with open(self.ffile, 'r') as fh:

            for line in fh:

                if "MOLTYPES" in line:
                    moltype_num = str(int(line.split()[1]) + 1)
                    cont.append(line.replace(line.split()[1],moltype_num))
                elif molecule_name in line:
                    cont.append(line)
                    mol_cont.append(line.replace(molecule_name,molecule_name + "_L"))
                    for i in range(1000):
                        atom_c_line = fh.readline()
                        cont.append(atom_c_line)
                        mol_cont.append(atom_c_line)
                        if "FINISH" in atom_c_line:
                            break
                    for i in mol_cont:
                        cont.append(i)
                    
                    cont.append("lambdacharge\n")
                else:
                    cont.append(line)

        with open (output_filename , 'w') as fw:
            for i in cont:
                fw.write(i)
        
        atom_names = self.unique_atomise(molecule_name=molecule_name + "_L")
        atom_names = [x + "_" + molecule_name[0] for x in atom_names]

        with open (output_filename, 'r') as fh:
            cont = []

            for line in fh:

                if "vdw" in line:
                    cont.append(line)
                    for i in range(int(line.split()[1])):
                        vdw_line = fh.readline()
                        vdw_list = vdw_line.split()
                        vdw_names = [vdw_list[0],vdw_list[2]]
                        if any([x in atom_names for x in vdw_names]):
                            vdw_list[4] = "+lambdalj"

                            lj_cos = list(map(float,[vdw_list[5],vdw_list[6]]))
                            delta = (lj_cos[0]/ lj_cos[1])**(1/6)
                            epsilon = lj_cos[1]**2/(4*lj_cos[0])
                            vdw_list[5] = f'{epsilon:.9f}'
                            vdw_list[6] = f'{delta:.9f}'
                            vdw_list.append(coeff)
                            cont.append("    ".join(vdw_list) + "\n")
                        else:
                            cont.append(vdw_line)
                    cont.append("LAMBDAPATH\n")
                    cont.append("vdwstart {}\n".format(lpath[0]))
                    cont.append("vdwend {}\n".format(lpath[1]))
                    cont.append("coulstart {}\n".format(lpath[2]))
                    cont.append("coulend {}\n".format(lpath[3]))
                    cont.append("FINISH\n")
                else:
                    cont.append(line)

        with open (output_filename, 'w') as fw:
            for i in cont:
                fw.write(i)

        

    def cut_field(self,cut_rat, retain = False):
        
        """
        Reduce the system by a specific ratio. provided that the system consists of repeatable units
        and the number of repeated units is a multiple of the cutting ratio.
        
        Arguments:
            cut_rat (int) ==> cutting ratio, e.g a cutting ratio of 2 cuts the system by a half.
        
        Returns:
            none.
        
        """
        
        con = []
        with open (self.ffile,'r') as fh:
            for i in range (10000):
                temp= fh.readline()
                if "nummols" in temp:
                    if retain:
                        con.append(temp)
                        continue
                    else:
                        temp = "nummols {}\n".format(cut_rat)
                if "atoms" in temp:
                    atomnum = int(temp.split()[1])
                    if atomnum % cut_rat != 0:
                        print("total atom number is no wholely dividable by inputed cut ratio")
                        return 0
                    cut_num = int(atomnum/cut_rat)
                    con.append("atoms {}\n".format(str(cut_num)))
                    for i in range(cut_num):
                        con.append(fh.readline())
                    break
                
                con.append(temp)
        bondtitle, bonds = self.ffetch ("bonds")
        angletitle, angles = self.ffetch ("angles")
        dihedraltitle, dihedrals = self.ffetch("dihedral")
        
        bondnum = str(int(int(bondtitle.split()[1])/cut_rat))
        con.append("bonds {}\n".format(bondnum))
        bonds = self.block_filter(bonds,cut_num,2)
        bonds = self.spacecreator(bonds)
        con.extend(bonds)
        
        anglenum = str(int(int(angletitle.split()[1])/cut_rat))      
        con.append("angles {}\n".format(anglenum))
        angles = self.block_filter(angles,cut_num,3)
        angles = self.spacecreator(angles)
        con.extend(angles)
        
        dihednum = str(int(int(dihedraltitle.split()[1])/cut_rat))      
        con.append("dihedral {}\n".format(dihednum))
        dihedrals = self.block_filter(dihedrals,cut_num,2)
        dihedrals = self.spacecreator(dihedrals)
        con.extend(dihedrals)
        con.append("finish\n")
        
        vdwtitle,vdw = self.ffetch("vdw")
        con.append(vdwtitle)
        con.extend(self.spacecreator(vdw))
        con.append("close")
        
        with open("FIELD",'w') as fw:
            for i in con:
                fw.write(i)
    
    


    def cut_config (self, cut_rat, mode = "REVCON"):
        
        """
        Analogous to cut_field, currently operatable on REVCON and CONFIG files only.
        Does not work with DLM inputs.
        """
        
        con = []
        
        if mode == "REVCON":
            blocksize = 4
        elif mode == "CONFIG":
            blocksize = 2
        else:
            print ("CONFIG file reduction mode is not allowed, available mode: CONFIG or REVCON")
            return 0
        
        with open("CONFIG",'r') as fh:
            if mode == "REVCON":
                con.append(fh.readline())
                temp = fh.readline().split()
                new_string = "         {}         {}\n".format (temp[0],temp[1])
                con.append(new_string)
            for line in fh:
                con.append(line)
        
        if (len(con)-5) % blocksize != 0:
            print ("CONFIG file missing some content...")
            return 0
        
        atomnum = (len(con)-5) / blocksize
        cut_num = int(atomnum/cut_rat)
        con = con[0:cut_num*blocksize+5]
        
        with open ("CONFIG",'w') as fw:
            for i in con:
                fw.write(i)
                
        return cut_num
    
    def block_filter (self,block, filnum, filcol):
        
        """
        Filters a specific interaction block, e.g if filnum is 260, any interaction involving
        atom number larger than 260 would be filtered out.
        
        Argumments:
            block (list(str)) ==> data block in space split format.
            filnum (int) ==> atom number to filter by.
            filcol (int) ==> number of columns to look for interation atom numbering.
            
        Returns:
            blockh (list(str)) ==> filtered data block.
        """
        
        blockh = []
        
        for i in block:
            if int(max(i[1:filcol+1])) <= filnum:
                blockh.append(i)
            else:
                continue
        return blockh
    
    def xyz_to_polycon(self, file, cell_str):

        data,dump = self.read_xyz(file)
        tail_string = "ends\nnopoten\nnosort\ncheck\nprint dlpoly 1\nstop\nstart\nstop\n"
        with open ("input.txt",'w') as fw:
            fw.write("cell\n")
            fw.write(cell_str+"\n")
            fw.write("basi\n")
            for i in data:
                fw.write(i)
            fw.write(tail_string)
        mu = MetaUtil()
        mu.run_metadise()


    def read_xyz(self,file):
        str_data = []
        array_data = []
        with open(file,'r') as fh:
            for line in fh:
                if len(line.split()) == 4:
                    str_data.append(line)
                    array_data.append(line.split())
        return str_data,array_data

    def map_config (self):
        title, data = self.ffetch("atoms")
        data = np.array(data)
        atmnum = len(data)
        cont = []
        count = 0
        with open (self.cfile,'r') as fh:
            for i in range(5):
                cont.append(fh.readline())
            for line in fh:
                work_str = line.split()
                if len(work_str) == 2:
                    ind = int(work_str[1])
                    line = line.replace(work_str[0],data[(ind-1)%atmnum][0])
                    count +=1
                cont.append(line)
                
        mol_num = int(count/atmnum)
        with open(self.cfile,'w') as fw:
            for i in cont:
                fw.write(i)
        cont = []

        with open(self.ffile,'r') as fh:
            for line in fh:
                if "nummols" in line:
                    line = line.replace(line.split()[1],str(mol_num))
                cont.append(line)

        with open(self.ffile,'w') as fw:
            for i in cont:
                fw.write(i)


    def write_config(self):
        """
        Writes DLM config file from DLP config file.
        
        Arguments:
            configfile (str) (optional) ==> defaults to "CONFIG" file name.
            
        Returns:
            none.
        """
        
        configcont = []
        molnum = self.get_nummol()
        molstr = ""
        moltot = 0
        
        for i in molnum:
            moltot += int(i[1])
            molstr += str(2*moltot)
            molstr += "    "
    
        nummolstr = "NUMMOL    " + str(moltot) + "    " + molstr + "\n" 
        
        with open(self.cfile,'r') as fh:
            
            count = 0

            for i in range(1000000):
                if count == 5:
                    break
                else:
                    temp = fh.readline()
                    count += 1

                if not temp:
                    count -= 1
                    continue
                
                elif count == 2:
                    configcont.append("    0    1 \n")
                
                else:
                    configcont.append(temp)
            
            configcont.append(nummolstr)
            
            for i in molnum:
                for j in range (int(i[1])):
                    configcont.append("molecule " +i[0] + "    " + i[2] + "    " + i[2] +"\n")
                    for k in range (2*int(i[2])):
                        temp = fh.readline()
                        if k % 2 == 0:
                            temp = temp.split()[0] + "    CORE \n"
                        configcont.append(temp)
            
            with open("CONFIG_DLM",'w') as fw:
                for i in configcont:
                    fw.write(i)

