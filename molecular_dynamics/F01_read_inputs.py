# read variables from input file
import csv

def read_inputs():
    with open('RUN_PARAMETERS.cfg') as f:
        for line in f:
            a = line.split()
            first = a[0][0]
            if first != ';':
                vartype = a[0][:]
                varvalue = a[2][:]
                if vartype == 'dt':
                    global DT
                    DT = float(varvalue)
                elif vartype == 'nsteps':
                    global TOT_STEPS
                    TOT_STEPS = int(varvalue)
                elif vartype == 'PRESS_INTVL':
                    global PRESS_INTVL
                    PRESS_INTVL = int(varvalue)
                elif vartype == 'SLAB_THICK':
                    global SLAB_THICK
                    SLAB_THICK = float(varvalue)
                elif vartype == 'THERMOSTAT':
                    global THERMOSTAT
                    THERMOSTAT = int(varvalue)
                elif vartype == 'T_REF':
                    global T_REF
                    T_REF = float(varvalue)
                elif vartype == 'sigma_ARPT':
                    global sigma_ARPT
                    sigma_ARPT = float(varvalue)
                elif vartype == 'eps_ARPT':
                    global eps_ARPT
                    eps_ARPT = float(varvalue)
                elif vartype == 'rc_wall':
                    global rc_wall
                    rc_wall = float(varvalue)
                elif vartype == 'RCRIT':
                    global RCRIT
                    RCRIT = float(varvalue)
                elif vartype == 'INIT_MOMENT':
                    global INIT_MOMENT
                    INIT_MOMENT = int(varvalue)
                elif vartype == 'RCUT':
                    global RCUT
                    RCUT = float(varvalue)
                elif vartype == 'NST_COMM':
                    global NST_COMM
                    NST_COMM = int(varvalue)
                elif vartype == 'nstxout':
                    global nstxout
                    nstxout = int(varvalue)
                elif vartype == 'PBC':
                    global PBC
                    PBC = int(varvalue)
                elif vartype == 'WALL_MODEL':
                    global WALL_MODEL
                    WALL_MODEL = int(varvalue)
                elif vartype == 'top_w':
                    global top_w
                    top_w = float(varvalue)
                elif vartype == 'bot_w':
                    global bot_w
                    bot_w = float(varvalue)
                elif vartype == 'TAU':
                    global TAU
                    TAU = float(varvalue)
                else:
                    print "[ERROR 100]: Input file read error"
                    print "parameter no found: ",vartype
            #end of if ";"
        #end of for loop
        if NST_COMM == 0: NST_COMM = TOT_STEPS
        print "File read successfully.."
        
def allot_memory(ARG):
    global ACCEL, FORCE #, VELO
    ACCEL = FORCE = [] # VELO = []
    for i in range(ARG):
        ACCEL.append([0,0,0])
        FORCE.append([0,0,0])
        #VELO.append([0,0,0])
def read_argon():
    with open('argon.csv','r') as f:
        reader = csv.reader(f, delimiter=",")
        data = list(reader)
        global ARGON, VELO
        ARGON = []
        VELO = []
        global NARGON
        NARGON = int(data[0][0])
        print "N_ARGON = ", NARGON
        for i in range(NARGON):
            b = []
            c = []
            b.append(float(data[i+1][3]))
            b.append(float(data[i+1][4]))
            b.append(float(data[i+1][5]))
            c.append(float(data[i+1][6]))
            c.append(float(data[i+1][7]))
            c.append(float(data[i+1][8]))
            ARGON.append(b)
            VELO.append(c)
        global box_L,box_B,box_H
        box_L = float(data[-1][0])
        box_B = float(data[-1][1])
        box_H = float(data[-1][2])
        print "Read argon successfully"
        allot_memory(NARGON)
        print "Allotted initial memory"
        #print data
