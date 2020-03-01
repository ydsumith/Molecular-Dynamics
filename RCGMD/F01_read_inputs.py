# read variables from input file
import F07_memory as mem
import csv
import sys
from sys import platform as _platform
import os


# ----------------------------------------------------------------
def read_inputs():
    with open('out_log.log', 'w') as f:
        f.write('Log file \n')
    f.close()
    with open('out_traj_pyth.mol', 'w') as f:
        f.close()
    with open('out_energy.res', 'w') as f:
        f.write('Step[ps]\tE_Total[kJ/mol]\tT[K]\tKE[kJ/mol]\tLJ[kJ/mol]\tBOND[kJ/mol]\tANGLE[kJ/mol]\n')
    f.close()
    with open('out_pcg_log.log', 'w') as f:
        f.write('step\tnE4\tnD4\n')
    f.close()

    if _platform == "linux" or _platform == "linux2":
        # linux
        mem.path = os.getcwd() + "/input_files/"
        print mem.path
    elif _platform == "darwin":
        # MAC OS X
        mem.path = os.getcwd() + "/input_files/"
        print mem.path
    elif _platform == "win32":
        # Windows
        mem.path = os.getcwd() + "\input_files\\"
        print mem.path

    with open(mem.path + 'RUN_PARAMETERS.cfg') as f:
        for line in f:
            a = line.split()
            first = a[0][0]
            if first != ';':
                vartype = a[0][:]
                varvalue = a[2][:]
                if vartype == 'dt':
                    mem.DT = float(varvalue)
                elif vartype == 'nsteps':
                    mem.TOT_STEPS = int(varvalue)
                elif vartype == 'THERMOSTAT':
                    mem.THERMOSTAT = int(varvalue)
                elif vartype == 'T_REF':
                    mem.T_REF = float(varvalue)
                elif vartype == 'INIT_MOMENT':
                    mem.INIT_MOMENT = int(varvalue)
                elif vartype == 'RCUT':
                    mem.RCUT = float(varvalue)
                    mem.rcut_sqr = mem.RCUT * mem.RCUT
                elif vartype == 'NST_COMM':
                    mem.NST_COMM = int(varvalue)
                elif vartype == 'nstxout':
                    mem.nstxout = int(varvalue)
                elif vartype == 'pcg_morse_lj':
                    mem.pcg_morse_lj = int(varvalue)
                elif vartype == 'pcg_thro_start':
                    mem.pcg_thro_start = int(varvalue)
                elif vartype == 'pcg_poly_start':
                    mem.pcg_poly_start = int(varvalue)
                elif vartype == 'TAU':
                    mem.TAU = float(varvalue)
                else:
                    print "[ERROR 100]: Input file read error"
                    print "parameter no found: ", vartype
                    mem.error_flag = 1
                    # end of if ";"
        # end of for loop
        if mem.pcg_thro_start > mem.TOT_STEPS or mem.pcg_poly_start > mem.TOT_STEPS or mem.pcg_thro_start > mem.pcg_poly_start:
            print "[Error 101]: pcg_thro_start < pcg_poly_start < nsteps"
            mem.error_flag = 1
        if mem.NST_COMM == 0: mem.NST_COMM = mem.TOT_STEPS
        print "cfg File read successfully.."


# ----------------------------------------------------------------
def read_data():
    with open(mem.path + 'inp_1pcgmd.data', 'r') as f:
        massstart = 0
        linenum = 1
        for line in f:
            if linenum > 1:
                content = [i for i in line.split()]
                length = len(content)
                if length == 1:
                    massstart = 1
                if length == 2:
                    varvalue = content[0]
                    vartype = content[1]
                    if massstart == 0:
                        if vartype == 'atoms':
                            mem.N = int(varvalue)
                        elif vartype == 'bonds':
                            mem.NB = int(varvalue)
                        elif vartype == 'angles':
                            mem.NA = int(varvalue)
                        elif vartype == 'atom-types': # atom types
                            mem.NTYPES = int(varvalue)
                        elif vartype == 'bond-types': # bond types
                            mem.NBTYPES = int(varvalue)
                        elif vartype == 'angle-types': # angle types
                            mem.NATYPES = int(varvalue)
                        elif vartype == 'xlo':
                            mem.xlo = float(varvalue) / 10.0
                        elif vartype == 'xhi':
                            mem.xhi = float(varvalue) / 10.0
                        elif vartype == 'ylo':
                            mem.ylo = float(varvalue) / 10.0
                        elif vartype == 'yhi':
                            mem.yhi = float(varvalue) / 10.0
                        elif vartype == 'zlo':
                            mem.zlo = float(varvalue) / 10.0
                        elif vartype == 'zhi':
                            mem.zhi = float(varvalue) / 10.0
                        else:
                            print "Error 0Xss"
                            mem.error_flag = 1
                    else:
                        mem.massvalue.append(float(vartype))
            linenum += 1
        print "inp_1pcgmd.data file read successfully.."

        mem.pcg_neigh_list = [0 for i in range(mem.NTYPES + 1)]

        mem.L = mem.xhi - mem.xlo
        mem.B = mem.yhi - mem.ylo
        mem.H = mem.zhi - mem.zlo
        mem.half_L = mem.L/2.0
        mem.half_B = mem.B/2.0
        mem.half_H = mem.H/2.0

        read_atoms(mem.N)
        calc_mass(mem.N)
        read_bonds(mem.NB)
        read_angles(mem.NA)

        create_coeffs()
        create_pcg_sets()

        read_coeff()
        build_lj_params()



# ----------------------------------------------------------------
def read_atoms(N):
    with open(mem.path + 'inp_2atoms.data', 'r') as f:
        try:
            linenum = 1
            reader = csv.reader(f, delimiter=",")
            data = list(reader)
            for i in range(N):
                mem.MOLID.append(int(data[i][1]))
                atom_type = int(data[i][2])
                mem.TYPE.append(atom_type)
                mem.CHARGE.append(float(data[i][3]))

                # mem.pcg_type_list[atom_type] += 1  # counting stars for each type

                mem.X.append(float(data[i][4]) / 10.0)
                mem.Y.append(float(data[i][5]) / 10.0)
                mem.Z.append(float(data[i][6]) / 10.0)

                mem.ACCX.append(0)
                mem.ACCY.append(0)
                mem.ACCZ.append(0)
                mem.FX.append(0)
                mem.FY.append(0)
                mem.FZ.append(0)
                mem.VX.append(0)
                mem.VY.append(0)
                mem.VZ.append(0)
                mem.MASS.append(0)
                linenum += 1
        except:
            print "unknown error read_atoms"
            mem.error_flag = 1

    print "%d Atoms file read successfully" % (linenum - 1)


# ---------------create_coeffs()-------------------------------------------------
def read_bonds(NB):
    with open(mem.path + 'inp_3bonds.data', 'r') as f:
        try:
            linenum = 1
            reader = csv.reader(f, delimiter=",")
            data = list(reader)
            for i in range(NB):
                mem.BTYPE.append(int(data[i][1]))
                mem.B1.append(int(data[i][2]))
                mem.B2.append(int(data[i][3]))
                linenum += 1
        except:
            print "unknown error read_bonds"
            mem.error_flag = 1

    print "%d Bonds file read successfully" % (linenum - 1)


# ----------------------------------------------------------------
def read_angles(NA):
    with open(mem.path + 'inp_4angles.data', 'r') as f:
        try:
            linenum = 1
            reader = csv.reader(f, delimiter=",")
            data = list(reader)
            for i in range(NA):
                mem.ATYPE.append(int(data[i][1]))
                mem.A1.append(int(data[i][2]))
                mem.A2.append(int(data[i][3]))
                mem.A3.append(int(data[i][4]))
                # print "%f %f %f %f" % (ATYPE[i], A1[i], A2[i], A3[i])
                linenum += 1
        except:
            print "unknown error read_angles"
            mem.error_flag = 1

    print "%d angles file read successfully" % (linenum - 1)


# ----------------------------------------------------------------
def calc_mass(N):
    for i in range(N):
        mem.MASS[i] = mem.massvalue[mem.TYPE[i]]


# ----------------------------------------------------------------
def read_coeff():
    with open(mem.path + 'inp_5coeff.data', 'r') as f:
        for line in f:
            content = [i for i in line.split()]
            length = len(content)

            if length > 1:
                vartype = content[0]

                if vartype == 'pair_coeff':
                    if content[1] == '*':
                        itype = content[1]
                    else:
                        itype = int(content[1])
                    if content[2] == '*':
                        jtype = content[2]
                    else:
                        jtype = int(content[2])
                    eps = float(content[3])
                    sig = float(content[4])

                    set_pair_coeff(itype, jtype, eps, sig)

                elif vartype == 'bond_coeff':
                    btype = int(content[1])
                    k = float(content[2])
                    d0 = float(content[3])

                    mem.bond_k[btype] = k
                    mem.bond_d0[btype] = d0

                elif vartype == 'angle_coeff':
                    atype = int(content[1])
                    k = float(content[2])
                    d0 = float(content[3])

                    mem.ang_k[atype] = k
                    mem.ang_t0[atype] = d0 * 0.01745329251994329576923690768489 #converting deg to rad

                elif vartype == 'pcgbond_coeff':
                    itype = int(content[1])
                    jtype = int(content[2])
                    k1 = float(content[3])
                    rbond = float(content[4])

                    mem.r_bond.append(rbond)

                    mem.pcg_bond[itype][jtype] = k1 * rbond * rbond
                    mem.pcg_bond[jtype][itype] = k1 * rbond * rbond

                elif vartype == 'pcgpair_coeff':
                    itype = int(content[1])
                    jtype = int(content[2])
                    k2 = float(content[3])
                    d0 = float(content[4])

                    mem.pcg_pair_A[itype][jtype] = k2 * pow(d0, 8.0)
                    mem.pcg_pair_A[jtype][itype] = k2 * pow(d0, 8.0)
                    mem.pcg_pair_B[itype][jtype] = k2 * pow(d0, 4.0)
                    mem.pcg_pair_B[jtype][itype] = k2 * pow(d0, 4.0)

                    mem.pcg_d0[itype][jtype] = d0
                    mem.pcg_d0[jtype][itype] = d0
                    mem.pcg_k2[itype][jtype] = k2
                    mem.pcg_k2[jtype][itype] = k2

                elif vartype == 'pcg_val':
                    itype = int(content[1])
                    jtype = int(content[2])
                    valency = int(content[3]) #currently not used

                    mem.pcg_val[itype][jtype] = valency  #currently not used
                    mem.ipcg.append(itype)
                    mem.jpcg.append(jtype)

                else:
                    print "unknown vartype in inp_5coeff.data", vartype
                    mem.error_flag = 1

        check_coeff_errors()

    print "coefficients file read successfully"


# ----------------------------------------------------------------
def create_coeffs():
    ntypes = mem.NTYPES
    btypes = mem.NBTYPES
    atypes = mem.NATYPES
    # ----------[[0 for i in range(cols_count)] for i in range(rows_count)]
    mem.pair_eps = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.pair_sig = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.bond_k = [0 for i in range(btypes + 1)]
    mem.bond_d0 = [0 for i in range(btypes + 1)]
    mem.ang_k = [0 for i in range(atypes + 1)]
    mem.ang_t0 = [0 for i in range(atypes + 1)]
    mem.lj_A = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.lj_B = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    print("coeff memory created")


# ----------------------------------------------------------------
def set_pair_coeff(itype, jtype, eps, sig):
    ntypes = mem.NTYPES
    if itype == '*' and jtype == '*':
        for i in range(1, ntypes + 1):
            for j in range(1, ntypes + 1):
                mem.pair_eps[i][j] = eps
                mem.pair_sig[i][j] = sig

    elif itype == '*':
        for i in range(1, ntypes + 1):
            mem.pair_eps[i][jtype] = eps
            mem.pair_sig[i][jtype] = sig
        for i in range(1, ntypes + 1):
            mem.pair_eps[jtype][i] = eps
            mem.pair_sig[jtype][i] = sig

    elif jtype == '*':
        for j in range(1, ntypes + 1):
            mem.pair_eps[itype][j] = eps
            mem.pair_sig[itype][j] = sig
        for j in range(1, ntypes + 1):
            mem.pair_eps[j][itype] = eps
            mem.pair_sig[j][itype] = sig

    else:
        mem.pair_eps[itype][jtype] = eps
        mem.pair_sig[itype][jtype] = sig
        mem.pair_eps[jtype][itype] = eps
        mem.pair_sig[jtype][itype] = sig

# ----------------------------------------------------------------
def check_coeff_errors():
    ntypes = mem.NTYPES

    for i in range(1, ntypes + 1):
        for j in range(1, ntypes + 1):
            if mem.pair_eps[i][j] == 0 or mem.pair_sig[i][j] == 0:
                print "[Warning]: type %s %s pair interaction is zero!" % (i, j)
    print "coeff checking complete."

    if mem.RCUT > mem.xhi or mem.RCUT > mem.yhi or mem.RCUT > mem.zhi:
        mem.error_flag = 1
        print "ERROR: rcut (%s) > box size (%s, %s, %s)" % (mem.RCUT, mem.xhi, mem.yhi, mem.zhi)


# ----------------------------------------------------------------
def build_lj_params():
    ntypes = mem.NTYPES
    for i in range(1, ntypes + 1):
        for j in range(1, ntypes + 1):
            mem.lj_A[i][j] = 4 * mem.pair_eps[i][j] * pow(mem.pair_sig[i][j], 12)
            mem.lj_B[i][j] = 4 * mem.pair_eps[i][j] * pow(mem.pair_sig[i][j], 6)


# ----------------------------------------------------------------
def create_pcg_sets():
    N = mem.N
    ntypes = mem.NTYPES
    types = mem.TYPE

    # ----------[[0 for i in range(cols_count)] for i in range(rows_count)]
    mem.pcg_neigh_list = [[0 for i in range(ntypes + 1)] for j in range(N + 1)]
    mem.pcg_bond = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.pcg_val = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.pcg_pair_A = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.pcg_pair_B = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]

    mem.pcg_d0 = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    mem.pcg_k2 = [[0 for i in range(ntypes + 1)] for j in range(ntypes + 1)]
    print "pcg related memory created"

    ntype2 = 0
    ntype3 = 0
    ntype4 = 0

    for i in range(N):
        if types[i] == 2:
            ntype2 += 1
        elif types[i] == 3:
            ntype3 += 1
        elif types[i] == 4:
            ntype4 += 1
        else:
            pass

    mem.ntype2 = ntype2
    mem.ntype3 = ntype3
    mem.ntype4 = ntype4

    mem.type2list = [0 for j in range(ntype2)]
    mem.type3list = [0 for j in range(ntype3)]
    mem.type4list = [0 for j in range(ntype4)]

    ntype2 = 0
    ntype3 = 0
    ntype4 = 0

    for i in range(N):
        if types[i] == 2:
            mem.type2list[ntype2] = i
            ntype2 += 1
        elif types[i] == 3:
            mem.type3list[ntype3] = i
            ntype3 += 1
        elif types[i] == 4:
            mem.type4list[ntype4] = i
            ntype4 += 1
        else:
            pass