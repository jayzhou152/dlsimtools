# small cif processing tool
import re

class CifTools:
    
    file = []
    
    def __init__(self, filename):
        with open (filename, 'r') as fh:
            for line in fh:
                self.file.append(line)
        
    def show (self):
        for i in self.file:
            print (i)
    
    def sort (self):
        
        revcon = []

        for line in self.file:
            line = line.lstrip()
            if line:
                if line[0] == "#":
                    continue
                else:
                    line = [s.strip() for s in line.split('  ') if s.split()]
                    revcon.append(line)
        return revcon
    
    def get_coords (self, labelsite, atomnum, siteini = '_atom_site_fract_z'):
        
        con = self.sort()
        
        start = con.index([siteini])
        end = start + atomnum
        info = con[start+1:end]
        
        labels = []
        coords = []
        
        for i in info:
            i = i[0].split()
            labels.append(i[labelsite-1])
            print (i)
            coords.append (i[labelsite:labelsite+3])
        coords = [map(float,coord) for coord in coords]
        
        return labels,coords
    
    def get_abc (self):
        
        con = self.sort()
        a = 0
        b = 0
        c = 0
        for i in con:
            if i[0] == '_cell_length_a':
                a = float(re.sub(r"\([^)]+\)", "", i[1]))
            if i[0] == '_cell_length_b':
                b = float(re.sub(r"\([^)]+\)", "", i[1]))
            if i[0] == '_cell_length_c':
                c = float(re.sub(r"\([^)]+\)", "", i[1]))
        return a,b,c
        
    def get_angles (self):
        con = self.sort()
        alpha,beta,gamma = 0,0,0
        for i in con:
            if i[0] == '_cell_angle_alpha':
                alpha = float(re.sub(r"\([^)]+\)", "", i[1]))
            if i[0] == '_cell_angle_beta':
                beta = float(re.sub(r"\([^)]+\)", "", i[1]))
            if i[0] == '_cell_angle_gamma':
                gamma = float(re.sub(r"\([^)]+\)", "", i[1]))
        return alpha,beta,gamma

    def get_symn (self):
        con = self.sort()
        flag = 0
        count = 0
        for i in con:
            if flag == 1:
                if "_" in i[0]:
                    print ("symmetries =  {}".format(count))
                    return count
                else:
                    count += 1
            if i[0] == "_symmetry_equiv_pos_as_xyz":
                flag = 1
        return 0

