import numpy as np
import subprocess
import matplotlib.pyplot as plt

class traj_analysis():

    def __init__ (self,filename="XYZ.000",format="xyz",n_points = 10000):
        self.filename = filename 
        self.format = format
        self.vectors, self.data = self.read_xyz_traj(n_points)

    def read_xyz_traj(self, n_points):
        
        data = {}
        count = 0

        with open (self.filename,'r') as fh:
            for line in fh:
                data_log = []

                for i in range(int(line.split()[0])+1):
                    atom_line = fh.readline()
                    if i == 0:
                        atom_line = atom_line.split("\"")[1]
                        vectors = [atom_line.split()[0],atom_line.split()[4],atom_line.split()[8]]
                        vectors = np.array(list(map(float,vectors)))
                        continue
                    else:
                        atom_line = atom_line.split()
                        atom_line = [atom_line[0], float(atom_line[1]),float(atom_line[2]),float(atom_line[3])]
                        data_log.append(atom_line)
                data[count] = np.array(data_log,dtype=object)
                count += 1
                if count >= n_points:
                    return vectors, data
        return vectors, data
    
    def av_pos (self):

        total = np.array([])

        for i in self.data.keys():
            if i == 0:
                total = self.data[i][:,1:4]
            else:
                total += self.data[i][:,1:4]
        inc = len(self.data)
        total = total/inc
        total = np.array(total,dtype = float)
        total = np.round(total,6)
        return np.concatenate((self.data[0][:,0].reshape(len(total),1),total),axis=1)
    
    def xyz_density_profile(self, atomfilter = "OT3", ngrid = 100, mode = "normal"):

        collection = []
        his_range = (-self.vectors[0]/2,self.vectors[0]/2)
        for i in self.data.keys():
            filtered_data = np.array([j[1] for j in self.data[i] if j[0] == atomfilter],dtype=float)
            n,x = np.histogram(filtered_data, bins = ngrid, range = his_range)
            if i == 0:
                collection = n 
            else:
                collection += n

        x_cent = (x[:-1] + x[1:])/2

        if mode == "normal":
            plt.plot(x_cent,collection)

            with open ("density_prof.txt",'w') as fw:
                for i in range(len(x_cent)):
                    fw.write(str(x_cent[i]) + "    " + str(collection[i]) + "\n" )
        elif mode == "rescaled":
            data = np.vstack((x_cent,collection)).T
            data = np.array([i for i in data if i[1]>0.01* max(collection)])
            neg_x = np.array([i for i in data[:,0] if i < 0])
            neg_x = -neg_x + max(neg_x)
            pos_x = np.array([i for i in data[:,0] if i >= 0]) 
            pos_x = -pos_x + max(pos_x) + max(neg_x)
            data[:,0] = np.concatenate([neg_x,pos_x],axis=0)
            data = data[data[:, 0].argsort()]
            with open ("density_prof.txt", 'w') as fw:
                for i in data:
                    fw.write(str(i[0]) + "    " + str(i[1]) + "\n")
            plt.plot(data[:,0],data[:,1])
        else:
            raise ValueError("Plot mode must be 'normal' or 'rescaled'.")

    def out_put_av_con(self, atmdict,filename="av_con",title ="traj_created_average_config", mis = "    0    1 "):
        

        with open(filename, 'w') as fw:
            
            fw.write(title + "\n")
            fw.write(mis + "\n")
            mat = np.array([[self.vectors[0],0,0],[0,self.vectors[1],0],[0,0,self.vectors[2]]])
            cellmat_str = ["{:>20}{:>20}{:>20}\n".format("%.10f"%mat[0,0],"%.10f"%mat[0,1],"%.10f"%mat[0,2]),
                       "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[1,0],"%.10f"%mat[1,1],"%.10f"%mat[1,2]),
                       "{:>20}{:>20}{:>20}\n".format("%.10f"%mat[2,0],"%.10f"%mat[2,1],"%.10f"%mat[2,2])]
                    
            for i in cellmat_str:
                fw.write(i)
            
            data = self.av_pos()
            nummol = 0
            natms = 0
            for i in atmdict:
                if i[1] != "n":
                    nummol += i[1]
                    natms += i[1]*i[2]
                else:
                    atmdict[-1][1] = (len(data) - natms)/i[2]
                    if (atmdict[-1][1]).is_integer():
                        atmdict[-1][1] = int(atmdict[-1][1])
                    else:
                        raise ValueError("Molecule number for solvent is not a whole number.")
                    
                    nummol += atmdict[-1][1]

            nummol_str = "NUMMOL        {}".format(nummol)
            for i in range(len(atmdict)):
                nummol_str += "      10000"

            nummol_str += "\n"
            fw.write (nummol_str)
            cur_atm_count = 0
            for i in range(len(atmdict)):
                for j in range(atmdict[i][1]):
                    if j == 0:
                        block = np.array([atmdict[i][0],j,data[cur_atm_count + j*atmdict[i][2]: cur_atm_count + (j+1)*atmdict[i][2]]],dtype=object)
                    else:
                        block = np.vstack((block,np.array([atmdict[i][0],j,data[cur_atm_count + j*atmdict[i][2]: cur_atm_count + (j+1)*atmdict[i][2]]],dtype=object)))

                cur_atm_count += atmdict[i][1]*atmdict[i][2]
                for i in block:
                    fw.write ("MOLECULE {}    {}    {}\n".format(i[0],len(i[2]),len(i[2])))
                    for j in i[2]:
                        fw.write(" {}      c \n".format(j[0]))
                        fw.write("{:>15}{:>15}{:>15}\n".format("%.10f"%j[1],"%.10f"%j[2],"%.10f"%j[3]))
    

    def plot_atoms(self,atom_name, mode = "xy", zfilter = 0, plot_mode = "heatmap"):

        atom_data = []

        for i in range(len(self.data)):
            for j in self.data[i]:
                if zfilter > 0:
                    if j[0] in atom_name and j[1] > zfilter:
                        atom_data.append(j)
                elif zfilter < 0:
                    if j[0] in atom_name and j[1] < zfilter:
                        atom_data.append(j)               
                else:
                    if j[0] in atom_name:
                        atom_data.append(j)           
        atom_data = np.array(atom_data,dtype=object)
        if mode == "xy":
            x = np.array(atom_data[:,2],dtype=float)
            y = np.array(atom_data[:,3],dtype=float)
            xmin = -self.vectors[1]/2
            xmax = -xmin
            ymin = -self.vectors[2]/2
            ymax = -ymin
        elif mode == "xz":
            x = np.array(atom_data[:,1],dtype=float)
            y = np.array(atom_data[:,2],dtype=float)
            xmin = -self.vectors[0]/2
            xmax = -xmin
            ymin = -self.vectors[1]/2
            ymax = -ymin
        elif mode == "yz":
            x = np.array(atom_data[:,1],dtype=float)
            y = np.array(atom_data[:,3],dtype=float)
            xmin = -self.vectors[0]/2
            xmax = -xmin
            ymin = -self.vectors[2]/2
            ymax = -ymin
        else:
            raise ValueError("mode must be 'xy','xz','yz'.")
        #Z,xedges,yedges = np.histogram2d(atom_data[:,2],atom_data[:,3])
        #ax = fig.add_subplot(projection='3d')
        #plt.pcolormesh(xedges,yedges,Z.T)
        if plot_mode == "heatmap":
            from scipy.stats.kde import gaussian_kde

            k = gaussian_kde(np.vstack([x, y]))

            xi, yi = np.mgrid[xmin:xmax:x.size**0.5*1j,ymin:ymax:y.size**0.5*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))

            
            fig = plt.figure(figsize=(10,20*(abs(ymin/xmin))))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)

            # alpha=0.5 will make the plots semitransparent
            ax1.pcolormesh(xi, yi, zi.reshape(xi.shape), alpha=1)
            ax2.contourf(xi, yi, zi.reshape(xi.shape), alpha=1)

            ax1.set_xlim(x.min(), x.max())
            ax1.set_ylim(y.min(), y.max())
            ax2.set_xlim(xmin, xmax)
            ax2.set_ylim(ymin, ymax)
            fig.savefig("heatmap.png",dpi = 600)
        elif plot_mode == "scatter":
            cont = []
            for i in atom_data[:,0]:
                if i == "OT3":
                    cont.append("blue")
                elif i == "O":
                    cont.append("yellow")
                elif i == "OHP":
                    cont.append("orange")
                else:
                    cont.append("black")
            
            plt.scatter(x,y,c=cont,s=0.3)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)


