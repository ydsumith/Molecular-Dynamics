import numpy as np
import matplotlib.pyplot as plt
import math
import re


#%%
def read_raw_data(filename):
    with open(filename) as file_in:
        for line in file_in:
            temp1 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
            # An expected line should be 5117 atoms [N atoms]
            if len(temp1) == 1:
                N = int(temp1[0])
                print("No: of atoms read = %d " % N)
                break;
    galax = np.zeros((N,4))
    i = 0
    with open(filename) as file_in:
        for line in file_in:
            temp1 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
            # An expected line should be 1 2 -4.9725e+02 -2.475e+01 0.0e+00 0 0 0 [ID type x y z nx ny nz]
            if len(temp1) == 8:  # for a valid entry
                galax[i][0] = float(temp1[2]) # x value
                galax[i][1] = float(temp1[3]) # y value
                galax[i][2] = float(temp1[4]) # z value
                galax[i][3] = int(temp1[1]) # type
                i += 1
    N = i
    return galax, N

#%%
def write_lammps_data(galaxy, xlonew, xhinew,
                      ylonew, yhinew, zlonew, zhinew,N):
    with open("raw_data_minimized.data", "w") as file_out:
        file_out.write("\n%d atoms\n" % N)
        file_out.write("2 atom types\n\n")
        file_out.write("%f %f xlo xhi\n" % (xlonew,xhinew) )
        file_out.write("%f %f ylo yhi\n" % (ylonew,yhinew) )
        file_out.write("%f %f zlo zhi\n\n" % (zlonew,zhinew) )
        file_out.write("Masses\n\n")
        file_out.write("%d %f\n" % (1,72.0611) )
        file_out.write("%d %f\n\n" % (2,10) )
        file_out.write("Atoms\n\n")
        for i in range(N):
            file_out.write("%d %d %.4f %.4f %.4f\n" % (i+1, galaxy[i,3], galaxy[i,0],galaxy[i,1],galaxy[i,2]) )

#%%
def send_to_c(galaxy):
    N = len(galaxy)
    with open("data_to_C.data", "w") as file_out:
        for i in range(N):
            file_out.write("%.4f %.4f %.4f\n" % (galaxy[i,0],galaxy[i,1],galaxy[i,2]) )
    print("data written to data_to_C.data")
#-------------------------------------
#-------------------------------------
def read_from_c(N):
    galax = np.zeros((N,1))
    galay = np.zeros((N,1))
    galaz = np.zeros((N,1))
    i = 0
    with open("data_from_C.data") as file_in:
        for line in file_in:
            temp1 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
            if len(temp1) == 3:  # for a valid entry
                galax[i][0] = float(temp1[0]) # x value
                galay[i][1] = float(temp1[1]) # y value
                galaz[i][2] = float(temp1[2]) # z value
                i += 1
    return galax,galay,galaz
# %%
def main():
    print("Current version doesn't support the following features\n\
          1) Periodic Boundary check\n\
          2) Selective atom types")
    galaxy, N= read_raw_data("raw_data.data")
    print("N = %d" % (N))
    galaxy = galaxy[galaxy[:,3].argsort()] # sort to separate water and plat

    send_to_c(galaxy)
    galax,galay,galaz = read_from_c(N)


    #%% find the limits, atoms etc
    zhinew = max(galaxy[:,2])
    zlonew = min(galaxy[:,2])

    ylonew = min(galaxy[:,1])
    yhinew = max(galaxy[:,1])

    xlonew = min(galaxy[:,0])
    xhinew = max(galaxy[:,0])

    write_lammps_data(galaxy, xlonew, xhinew,
                      ylonew, yhinew, zlonew, zhinew,N)

    print("Program finished")

# %%
if __name__ == "__main__":
    main()
