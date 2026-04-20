from multiprocessing.sharedctypes import Value
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from molmass import ELEMENTS

class MonteCon():

    def __init__(self, filename = "CONFIG",stackd = "x", mol_divide = 2):
        
        self.ph =  self.read_in(filename=filename)
        self.data = self.ph[4]
        self.nummol = self.ph[3][0]
        self.mol_max = self.ph[3][1:]
        self.cell_mat = self.ph[2]
        self.cents = self.calc_centre()
        self.mol_type_num = self.ph[5]
        self.atm_count = self.ph[6]

        if stackd == "x":
            self.stackd = 0
        elif stackd == "y":
            self.stackd = 1
        elif stackd == "z":
            self.stackd = 2
        else:
            raise ValueError("stack direction must be x,y or z")

        self.mol_divide = mol_divide

        if len(self.data) != self.nummol:
            print("nummol in header does not match actual data, adjusting data length")
            self.data = self.data[:int(self.nummol)]

        else:
            print("data successfully imported")

    def read_in(self, filename = "CONFIG"):
        
        cell_mat = []

        with open (filename , 'r') as fh:
            title = fh.readline()
            mis = fh.readline()
            for i in range(3):
                line = fh.readline()
                cell_mat.append(list(map(float,line.split())))
            mol_info = fh.readline()
            nummol = [int(float(i)) for i in mol_info.split()[1:]]

            mol_type_num = len(mol_info.split()) - 2

            data = []
            mol_count = -1
            atm_count = 0
            for line in fh:
                if "MOLECULE" in line or "molecule" in line:
                    mol_count += 1
                    mol_name = line.split()[1]
                    coords_block = []
                    atom_num = line.split()[2]
                    atm_count += int(float(atom_num))
                    for i in range(int(float(atom_num))):
                        name_l = fh.readline()
                        name = name_l.split()[0]
                        coords_l = fh.readline()
                        coords = list(map(float,coords_l.split()[:3]))
                        atom_line = [name] + coords
                        coords_block.append(atom_line)
                    data.append([mol_name,mol_count,coords_block])

        return title, mis, np.vstack(cell_mat), nummol, data, mol_type_num, atm_count
    
    def get_slab_thickness (self,mol_name = "XYZ", mode = "centre"):
        if mode == "edge":
            data_col = [[j[1] for j in i[2]] for i in self.data if i[0] == mol_name]
            data_flat = [j for sub in data_col for j in sub]
            return max(data_flat) - min(data_flat)

        elif mode == "centre":
            data_flat = [i[2] for i in self.cents if i[0] == mol_name]
            return max(data_flat) - min(data_flat)
        else:
            raise ValueError("Solvent layer thickness mode must be 'edge' or 'centre'.")
    
    
    def density_profile(self, slab_thickness, mol = "TP3O", ngrid=100, mode = "edge_in"):

        excess_ = self.cell_mat[self.stackd,self.stackd] - slab_thickness

        ad_data = self.cents.copy()
        
        z_data = ad_data[:,2]

        n,x = np.histogram(z_data, bins = ngrid)
        n_prob, x = np.histogram(z_data, bins = ngrid, density = True)
        x_cent = (x[:-1] + x[1:])/2

        if mode == "normal":
            plt.plot(x_cent,n)
            plt.xlabel("z length (A)")
            plt.ylabel("molecule count (-)")
        elif mode == "probability":
            plt.plot(x_cent,n_prob)
            plt.xlabel("z length (A)")
            plt.ylabel("probability density (-)")

        elif mode == "edge_in":
            z_data_edge = [i + self.cell_mat[self.stackd,self.stackd] if i<0 else i for i in z_data]
            n, x = np.histogram(z_data_edge, bins = ngrid)
            x_cent = (x[:-1] + x[1:])/2

            plt.plot(x_cent,n)
            plt.xlim((0,self.cell_mat[self.stackd,self.stackd]))
            plt.xlabel("z length (A)")
            plt.ylabel("molecule count (-)")
        elif mode == "density":
            if mol == "TP3O":
                mol_weight = 18
            new_n = (n*mol_weight/1000/6.0221408e23)/(self.cell_mat[0,0]*self.cell_mat[1,1]*self.cell_mat[2,2]*1e-30/ngrid)
            plt.plot(x_cent,new_n)
            plt.ylabel("density (kg/m3)")
            plt.xlabel("z length(A)")


    def plot_crystal(self, mode = "3d", layers="all", numbering = True, no_ticks = False):

        if layers == "all":
            new_cents = self.cents
            ind = np.arange(len(new_cents))
        else:
            if isinstance(layers,list):
                ind, new_cents = self.cent_layers(layers)
            else:
                raise ValueError("if 'all' was not specified for layers, layers must be given as a list.")

        if mode == "3d":
            fig = plt.figure(figsize = (12,12))
            ax = fig.add_subplot(projection='3d')
            ax.scatter(new_cents[:,2],new_cents[:,3],new_cents[:,4],marker="o",s=100)
            if numbering:
                for i in range(len(ind)):
                    ax.text(new_cents[i,2],new_cents[i,3],new_cents[i,4],'%s' % (str(new_cents[i,1])),size=20,zorder=1,color='k')
            ax.set_xlabel('Z (A)')
            ax.set_ylabel('X (A)')
            ax.set_zlabel('Y (A)')
            ax.set_xlim(self.cell_mat[0,0]/2,-self.cell_mat[0,0]/2)
            ax.set_ylim(self.cell_mat[1,1]/ 2,-self.cell_mat[1,1]/2)
            ax.set_zlim(self.cell_mat[2,2]/2,-self.cell_mat[2,2]/2)
            ax.set_box_aspect((self.cell_mat[0,0],self.cell_mat[1,1],self.cell_mat[2,2]))
        elif mode == "2d":
            plt.scatter(new_cents[:,3],new_cents[:,4], marker = "x",color = 'red', s=100)
            if numbering:
                for i in range(len(ind)):
                    plt.text(new_cents[i,3],new_cents[i,4],str(new_cents[i,1]),size=10,zorder=1,color='k')
            plt.xlabel(r'$x (\mathring{A})$',fontsize = 16)
            plt.ylabel(r'$y (\mathring{A})$',fontsize = 16)
            plt.xlim(-self.cell_mat[1,1]/2,self.cell_mat[1,1]/2)
            plt.ylim(-self.cell_mat[2,2]/2,self.cell_mat[2,2]/2)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            if no_ticks:
                plt.xticks([])
                plt.yticks([])
            plt.gca().set_aspect(1)
            plt.gcf().set_size_inches(10,10)
        else:
            raise ValueError("mode is not recognised, must be '3d' or '2d'.")

    def calc_centre (self, atomname = "all", molname = "all"):

        cents = []

        for i in self.data:
            if molname != "all":
                if i[0] != molname:
                    continue

            if atomname == "all":
                coords = np.vstack([j[1:] for j in i[2]])
                atm_names = np.array([k[0] for k in i[2]])

                weights = []
                for h in atm_names:
                    if  h[0] == "H" or h[0] == "h":
                        weights.append(1)
                    elif h[0] == "O" or h[0] == "o":
                        weights.append(16)
                    elif h[0] == "C" or h[0] == "c":
                        weights.append(12)
                    elif h[0] == "N" or h[0] == "n":
                        weights.append(14)
                    else:
                        try:
                            elem = ELEMENTS[h[0].upper()]
                            weights.append(elem.mass)
                        except KeyError:
                            try:
                                elem = ELEMENTS[h[0].upper() + h[1].lower()]
                                weights.append(elem.mass)
                            except KeyError:
                                try:
                                    elem = ELEMENTS[h[0].upper() + h[1].lower() + h[2].lower()]
                                    weights.append(elem.mass)
                                except (KeyError, IndexError) as err:
                                    print(h)
                                    print(err)
                                    raise ValueError("ATOM NAME NOT RECOGNISED IN CONFIG, WEIGHTS NOT OBTAINTABLE")

            else:
                coords = np.vstack([j[1:] for j in i[2] if j[0] == atomname])
                atm_names = np.array([k[0] for k in i[2] if k[0] == atomname])

                weights = []
                for h in atm_names:
                    if  h[0] == "H" or h[0] == "h":
                        weights.append(1)
                    elif h[0] == "O" or h[0] == "o":
                        weights.append(16)
                    elif h[0] == "C" or h[0] == "c":
                        weights.append(12)
                    elif h[0] == "N" or h[0] == "n":
                        weights.append(14)
                    else:
                        try:
                            elem = ELEMENTS[h[0].upper()]
                            weights.append(elem.mass)
                        except KeyError:
                            try:
                                elem = ELEMENTS[h[0].upper() + h[1].lower()]
                                weights.append(elem.mass)
                            except KeyError:
                                try:
                                    elem = ELEMENTS[h[0].upper() + h[1].lower() + h[2].lower()]
                                    weights.append(elem.mass)
                                except (KeyError, IndexError) as err:
                                    print(err)
                                    raise ValueError("ATOM NAME NOT RECOGNISED IN CONFIG, WEIGHTS NOT OBTAINTABLE")

        
            weights = np.array(weights,dtype="float")
            weights = weights/np.sum(weights)
             
            result = [i[0],i[1],np.sum(coords[:,0]*weights),np.sum(coords[:,1]*weights),np.sum(coords[:,2]*weights)]
            cents.append(result)
        
        return np.array(cents,dtype ='object')
    
    def inter_slice_distance (self, layers = "all", mode = "2d"):

        if layers == "all":
            new_cents = self.cents
            ind = np.arange(len(new_cents))
        else:
            if isinstance(layers,list):
                ind, new_cents = self.cent_layers(layers)
            else:
                raise ValueError("if 'all' was not specified for layers, layers must be given as a list.")


        distances = {}

        xlen = self.cell_mat[1,1]
        ylen = self.cell_mat[2,2]

        for i in combinations(np.arange(len(new_cents)),2):
            one_co = new_cents[i[0],-2:]
            two_co = new_cents[i[1],-2:]

            xdis = abs(one_co[0]-two_co[0])
            ydis = abs(one_co[1]-two_co[1])
            xdis = min(xdis,xlen-xdis)
            ydis = min(ydis,ylen-ydis)
            dis = np.sqrt(xdis**2+ydis**2)
            distances[i] = dis
        return ind, distances 
    
    def removal_order (self,**kwargs):
        
        ind, distances = self.inter_slice_distance(**kwargs)
        ind_indices = list(np.arange(len(ind)))
        removal_order = []
        
        for atom_num in range(len(ind)):
            av_distances = []

            for i in ind_indices:
                total_dis = 0
                for j in ind_indices:
                    if i == j :
                        continue
                    else:
                        total_dis += distances[(min(i,j),max(i,j))]

                av_distances.append(total_dis)

            removal_order .append(ind_indices.pop(av_distances.index(max(av_distances))))
        
        return [ind[i] for i in removal_order]



    def cent_layers (self,layers_in):

        layers,layer_indices,layer_height = self.find_layers()

        if any([i > layers for i in layers_in]):
            raise ValueError ("layer provided is larger than the amount of layers in this crystal.")
        elif any([i < 1 for i in layers_in]):
            raise ValueError ("layer value cannot be smaller than 1.")
        else:
            pass
        
        layer_indices_flat = []
        for i in layers_in:
            layer_indices_flat.extend(layer_indices[i-1])

        layers_cents = [self.cents[i] for i in layer_indices_flat]
        
        return layer_indices_flat,np.vstack(layers_cents)

    def find_layers (self):

        layers = 0
        layer_indices = []
        cents = np.array([i[self.stackd+2] for i in self.cents])

        diff =  abs(cents[:-1] - cents[1:])

        layer_block=[0]
        count = 1

        for i in diff:


            if i >= self.mol_divide:
                layers += 1 
                layer_indices.append(layer_block)
                layer_block = []

            layer_block.append(count)
            count += 1
        
        layers += 1 
        layer_indices.append(layer_block)

        layer_height = []

        for i in layer_indices:
            layer_height.append(cents[i[0]])

        return layers, layer_indices, layer_height

    def change_mollabel (self, ind, new_label):
                        
        old_label = self.data[ind][0]

        self.data[ind][0] = new_label
        self.cents[ind,0] = new_label

        return old_label

    def seek_ind (self,ind):

        for i in self.data:
            if i[1] == ind:
                ac_ind = self.data.index(i)

        return ac_ind

    def del_mol (self, mol, write = False, plot = False, write_current = False, Tar_mol_label = False):

        mol = self.seek_ind(mol)

        if Tar_mol_label:
            self.change_mollabel(mol,"TAR_MOL")
        
        if write_current:
            self.write_out(filename = "CONFIG")
    
        self.data.pop(mol)
        self.cents = np.delete(self.cents,mol,axis=0)
        self.nummol = len(self.data)
            
        if plot:
            self.plot_crystal(mode="3d")

        if write:
            self.write_out(filename = "CONFIG_Final")
        

    def write_out (self, filename = "CONFIG_new", fformat = "DLMONTE"):
        if fformat == "DLMONTE":
            with open(filename, 'w') as fw:
                
                fw.write(self.ph[0])
                fw.write(self.ph[1])
                mat = np.vstack(self.ph[2])
                cellmat_str = ["{:>20}{:>20}{:>20}\n".format("%.10f"%mat[0,0],"%.10f"%mat[0,1],"%.10f"%mat[0,2]),
                        "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[1,0],"%.10f"%mat[1,1],"%.10f"%mat[1,2]),
                        "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[2,0],"%.10f"%mat[2,1],"%.10f"%mat[2,2])]
                        
                for i in cellmat_str:
                    fw.write(i)

                nummol_str = "NUMMOL        {}".format(self.nummol)
                for i in self.mol_max:
                    nummol_str += "      {}".format(i)
                nummol_str += "\n"
                fw.write (nummol_str)

                for i in self.data:
                    fw.write ("MOLECULE {}    {}    {}\n".format(i[0],len(i[2]),len(i[2])))
                    for j in i[2]:
                        fw.write(" {}      c \n".format(j[0]))
                        fw.write("{:>15}{:>15}{:>15}\n".format(j[1],j[2],j[3]))
    
        elif fformat == "xyz":
            with open(filename+".xyz",'w') as fw:
                
                fw.write(str(self.atm_count) + "\n")
                lat_str = "Lattice=\"{} {} {} {} {} {} {} {} {}\"".format(*self.cell_mat.flatten()) + " Origin=\"{} {} {}\"\n".format(*-0.5*np.diagonal(self.cell_mat))
                fw.write(lat_str)
                for i in self.data:
                    for j in i[2]:
                        fw.write("{} {} {} {}\n".format(*j))
                

    def remove_top_layer(self):

        info = self.find_layers()
        layer_info = info[1][-1]
        target_index = layer_info[int(len(layer_info)/2)]
        layer_info.pop(target_index)
        self.del_mol(layer_info)


    
