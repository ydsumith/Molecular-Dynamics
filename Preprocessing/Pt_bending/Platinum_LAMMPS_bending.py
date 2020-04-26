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
    NPt = 0
    NWat = 0
    with open(filename) as file_in:
        for line in file_in:
            temp1 = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
            # An expected line should be 1 2 -4.9725e+02 -2.475e+01 0.0e+00 0 0 0 [ID type x y z nx ny nz]
            if len(temp1) == 8:  # for a valid entry
                galax[i][0] = float(temp1[2]) # x value
                galax[i][1] = float(temp1[3]) # y value
                galax[i][2] = float(temp1[4]) # z value
                galax[i][3] = int(temp1[1]) # type
                if galax[i][3] == 1:
                    NWat += 1
                else:
                    NPt += 1
                i += 1
    return galax , NPt, NWat, N

#%%
def write_lammps_data(xnew, ynew, znew, xlonew, xhinew, ylonew, yhinew, zlonew, zhinew, Natoms,typer):
    with open("system.data", "w") as file_out:
        file_out.write("\n%d atoms\n" % Natoms)
        file_out.write("2 atom types\n\n")
        file_out.write("%f %f xlo xhi\n" % (xlonew,xhinew) )
        file_out.write("%f %f ylo yhi\n" % (ylonew,yhinew) )
        file_out.write("%f %f zlo zhi\n\n" % (zlonew,zhinew) )
        file_out.write("Masses\n\n")
        file_out.write("%d %f\n" % (1,72.0611) )
        file_out.write("%d %f\n\n" % (2,10) )
        file_out.write("Atoms\n\n")
        for i in range(Natoms):
            file_out.write("%d %d %.4f %.4f %.4f\n" % (i+1, typer[i], xnew[i],ynew[i],znew[i]) )

# %%
def main():
    xlo = -500 # Ang
    xhi = 500 # Ang
    dx = 10 # Ang
    ylo = -20
    yhi = 20
    chanl_heigh = 500 # 50 nm

    galaxy, NPt, NWat, N= read_raw_data("raw_data.data")

    print("NPt = %d, NWat = %d, N = %d" % (NPt, NWat, N))

    galaxy = galaxy[galaxy[:,3].argsort()] # sort to separate water and plat

    xplat = galaxy[NWat:N-1,0] # x coord of platinum
    yplat = galaxy[NWat:N-1,1] # y coord of platinum
    zplat = galaxy[NWat:N-1,2] # z coord of platinum

    xdata = np.copy(xplat)
    ydata = np.copy(zplat)

    NX = int((xhi-xlo)/dx)+1
    print("NX = %.2f" % NX)

    xbend = np.linspace(xlo, xhi, NX)
    ybend = np.random.uniform(low=ylo, high=yhi, size=(NX,))
    R = np.zeros((2, 2))

#%%  get the end point
    xmax = -100000
    xmaxloc = -1
    i = 0
    for x in xdata:
        if x > xmax:
            xmax = x
            xmaxloc = i
        i+=1
    print("xmax = %.2f, xmaxloc = %d" % (xmax, xmaxloc))

#%% make the shape
    ydata = ydata + ybend[0]
    j = 0
    for xb in xbend:
        yb = ybend[j]
        theta1 = math.atan( (ydata[xmaxloc]-yb) / (xdata[xmaxloc]-xb))
        if j == NX-1:
            theta2 = 0
        else:
            theta2 = math.atan( (ybend[j+1]-yb) / (dx))
        theta = theta2-theta1
        R[0][0] = math.cos(theta)
        R[0][1] = -math.sin(theta)
        R[1][0] = math.sin(theta)
        R[1][1] = math.cos(theta)
        i=0
        for x in xdata:
            if x >= xb:
                x_x = x - xb
                y_y = ydata[i] - yb
                x_p = x_x * R[0][0] + y_y * R[0][1]
                y_p = x_x * R[1][0] + y_y * R[1][1]
                xdata[i] = xb + x_p
                ydata[i] = yb + y_p
            i += 1
        j+=1


    plt.plot(xdata, ydata, 'o', color='red', markersize=2)
    plt.plot(xbend, ybend, 'o', color='blue', markersize=2)

#%% find the limits, atoms etc
    zhinew = max(max(ybend), abs(min(ybend)))
    zlonew = -zhinew

    ylonew = min(galaxy[NWat:N-1,1])
    yhinew = max(galaxy[NWat:N-1,1])

    xlonew = xbend[0]
    xhinew = xbend[NX-1]

    xnew = []
    ynew = []
    znew = []
    typer = []
    i = 0
    j = 0
    for x in xdata:
        if x > xlonew and x < xhinew:
            xnew.append(x)
            ynew.append(galaxy[j+NWat,1])
            znew.append(ydata[j])
            typer.append(2) # platinum
            i += 1
        j += 1
    Natoms = i
    # add water
    for i in range(NWat):
        if galaxy[i,0] > xlonew and galaxy[i,0] < xhinew:
            xnew.append(galaxy[i,0])
            ynew.append(galaxy[i,1])
            znew.append(galaxy[i,2] + zhinew)
            typer.append(1) # water
            Natoms += 1
    # add top platinum plate
    j = 0
    for x in xplat:
        if x > xlonew and x < xhinew:
            xnew.append(x)
            ynew.append(yplat[j])
            znew.append(zplat[j]+chanl_heigh)
            typer.append(2) # platinum
            Natoms += 1
        j += 1

    zhinew = max(znew) + 1.0
    zlonew = zlonew - 1.0

    write_lammps_data(xnew, ynew, znew, xlonew, xhinew, ylonew, yhinew, zlonew, zhinew, Natoms,typer)

    print("Program finished")


# %%

if __name__ == "__main__":
    main()
