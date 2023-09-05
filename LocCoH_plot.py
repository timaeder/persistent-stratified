import numpy as np

#graphic plotting
import matplotlib.pyplot as plt

#plotting high-dim data
from sklearn.decomposition import PCA

#plotting calculation progress
from tqdm import tqdm

#density filtration
from sklearn.neighbors import NearestNeighbors


def LocCoH_read(
              dataname, pairname, 
              path_data = "",
              path_output = ""):

    X = np.genfromtxt(path_data + dataname, delimiter=',')
    
    diags1 = []
    diags2 = []
    ClusterKeys1 = []
    ClusterKeys2 = []
    
    for i in range(0, len(X)):
    
        dat = np.genfromtxt(path_output + pairname + str(i) + ".txt",
                            delimiter=',')
        
        diag1 = []
        diag2 = []
    
        if(type(dat[0]) is np.ndarray):
    
            for j in range(0, len(dat)):
                if (dat[j][2] == 0):
                    if (dat[j][1] > 0.00001):
                        diag1.append([dat[j][0], dat[j][1]])
                        
            for j in range(0,len(dat)):
                if (dat[j][2] == 1):
                    if (dat[j][1] > 0.00001):
                        diag2.append([dat[j][0],dat[j][1]])
    
    
        if(len(diag1)):
            diags1.append(diag1)
            ClusterKeys1.append(1)
        else:
            diags1.append([[0.0, 0.0]])
            ClusterKeys1.append(0)
            
        if(len(diag2)):
            diags2.append(diag2)
            ClusterKeys2.append(1)
        else:
            diags2.append([[0.0, 0.0]])
            ClusterKeys2.append(0)
        
    return diags1, diags2, ClusterKeys1, ClusterKeys2
    

def PhiPLH(diags, CK, k = 1):
    
    HeatList1 = []
    HeatList2 = []
    HeatList3 = []
    
    HeatListA = []
    M = 0.0
    
    for i in range(0,len(diags)):
        
        maxsofar = 0.0
        plus2 = 0.0
        plus3 = 0.0
        plusA = 0.0
        if(CK[i]):
            for j in range(0,len(diags[i])):    
                
                plusA += diags[i][j][1]
                
                if(maxsofar<diags[i][j][1]):
                    maxsofar = diags[i][j][1]
                else:
                    if(plus2 < diags[i][j][1]):
                        plus2 = diags[i][j][1]
                        
                    else:
                        if(plus3 < diags[i][j][1]):
                            plus3 = diags[i][j][1]
                        
                
        
            HeatList1.append(maxsofar)
            HeatList2.append(plus2)
            HeatList3.append(plus3)
            HeatListA.append(plusA)
        else:
            HeatList1.append(0.0)
            HeatList2.append(0.0)
            HeatListA.append(0.0)
            
        if maxsofar > M:
            M = maxsofar
    
    print(M)
    
    HeatList = []
    
    for i in range(0,len(HeatList1)):
        
         HeatList.append(np.max([M - HeatList1[i],HeatList2[i]/2]))
        
    H = [HeatList, HeatList1, HeatList2, HeatListA]
        
    return H

def plot_data(data, HL=[], pca=False):
    
    dim = data.shape[1]
    
    if dim > 3:
        pca3 = PCA(n_components=3)
        pca3.fit(data)
        xyz = pca3.transform(data)

    else:
        if(dim == 2):
            xyz = np.zeros((len(data[:,0]),3))
            xyz[:,0] = data[:,0]
            xyz[:,1] = data[:,1]
            xyz[:,2] = np.zeros((len(data[:,0]),))
        else:
            xyz = data

    
    fig = plt.figure(figsize = (16,10))
    ax = fig.add_subplot(111,projection="3d")
    ax.view_init(elev=50, azim=110)

    if len(HL):
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],depthshade=False,picker=True,
                       c=HL)
    
    else:
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],depthshade=False,picker=True)

    plt.grid(False)
    plt.axis('off')



def main():

    #Insert name of the point cloud data as dataname and the output of the local homology computations as pairname.
    #If the data is not in the same repository as this python script, please specify path_data and path_output.

    DC = LocCoH_read(dataname="data.txt",pairname = "LocCoH_output")
    X = np.genfromtxt("data.txt",delimiter=",")

    Phi = PhiPLH(DC[0],DC[2])

    plot_data(X,HL=Phi[0])


if __name__ == "__main__":

    main()

