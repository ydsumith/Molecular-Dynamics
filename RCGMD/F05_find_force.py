import F07_memory as mem
from math import sqrt, exp
from math import acos, cos, sin
import multiprocessing
import math
import F06_output as out


#-------------------------------------------------------------------
def lj_worker(istart, iend, L, B, H, half_L, half_B, half_H,
                      RCUT2, rbond2, X, Y, Z, TYPE, N, lj_A, lj_B):
    try:
        outdict1 = []
        outdict2 = []
        outdict3 = []

        ndict1 = []
        ndict2 = []
        ndict3 = []
        ndict4 = []
        e_lj = 0

        print "istart = %s, iend = %s, N = %s" % (istart, iend, N)

        for i in range(istart, iend, 1):
            xi = X[i]
            yi = Y[i]
            zi = Z[i]
            itype = TYPE[i]

            sumFX = 0
            sumFY = 0
            sumFZ = 0
            neigh1 = 0
            neigh2 = 0
            neigh3 = 0
            neigh4 = 0

            for j in range(N):
                if i != j:
                    xj = X[j]
                    yj = Y[j]
                    zj = Z[j]
                    jtype = TYPE[j]

                    xij = xi - xj
                    yij = yi - yj
                    zij = zi - zj

                    if (abs(xij) > half_L):
                        if (xij < 0):
                            xij = xij + L
                        else:
                            xij = xij - L

                    if (abs(yij) > half_B):
                        if (yij < 0):
                            yij = yij + B
                        else:
                            yij = yij - B

                    if (abs(zij) > half_H):
                        if (zij < 0):
                            zij = zij + H
                        else:
                            zij = zij - H

                    rsqr = xij * xij + yij * yij + zij * zij

                    if (rsqr < RCUT2):

                        if rsqr < rbond2:  # found a neighbor for pcg
                            if jtype == 1:
                                neigh1 += 1
                            elif jtype == 2:
                                neigh2 += 1
                            elif jtype == 3:
                                neigh3 += 1
                            elif jtype == 4:
                                neigh4 += 1

                        irr = 1 / rsqr
                        ir6 = irr * irr * irr
                        ir12 = ir6 * ir6

                        e_lj = e_lj + (lj_A[itype][jtype] * ir12 - lj_B[itype][jtype] * ir6)
                        FORCE_MAG = (12 * lj_A[itype][jtype] * ir12 - 6 * lj_B[itype][jtype] * ir6) * irr

                        sumFX = sumFX + FORCE_MAG * xij
                        sumFY = sumFY + FORCE_MAG * yij
                        sumFZ = sumFZ + FORCE_MAG * zij
            outdict1.append(sumFX)
            outdict2.append(sumFY)
            outdict3.append(sumFZ)
            ndict1.append(neigh1)
            ndict2.append(neigh2)
            ndict3.append(neigh3)
            ndict4.append(neigh4)

        return outdict1, outdict2, outdict3, ndict1, ndict2, ndict3, ndict4, e_lj
    except Exception, e:
        print "ERROR: " , str(e)
        out.loggit("ERROR: %s " % str(e))
        mem.error_flag = 1
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
def find_parall_lj_force():
    N = mem.N
    X = mem.X
    Y = mem.Y
    Z = mem.Z
    TYPE = mem.TYPE
    # FX = mem.FX (dont do this! sync issues)

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    RCUT2 = mem.rcut_sqr

    rbond2 = mem.r_bond[0] * mem.r_bond[0]  # the index 0 has to be generalized
    lj_A = mem.lj_A
    lj_B = mem.lj_B

    # Virial = [0, 0, 0]
    e_lj = 0

    nprocs = mem.nprocs

    chunksize = int(math.ceil(N / float(nprocs)))
    procs = []

    pool = multiprocessing.Pool(nprocs)

    for i in range(nprocs):
        istart = chunksize * i
        iend = chunksize * (i + 1)
        if iend > N:
            iend = N
        if istart >= iend:
            out.loggit("ERROR: istart %s >= iend %s" % (istart, iend))

        p = pool.apply_async(lj_worker,(istart, iend, L, B, H, half_L, half_B, half_H,
                      RCUT2, rbond2, X, Y, Z, TYPE, N, lj_A, lj_B))
        procs.append(p)

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    resultdict1 = []
    resultdict2 = []
    resultdict3 = []
    neigh_dict1 = []
    neigh_dict2 = []
    neigh_dict3 = []
    neigh_dict4 = []
    e_lj = 0
    for p in procs:
        t1, t2, t3, t4, t5, t6, t7, t8 = p.get()
        resultdict1.extend(t1)
        resultdict2.extend(t2)
        resultdict3.extend(t3)
        neigh_dict1.extend(t4)
        neigh_dict2.extend(t5)
        neigh_dict3.extend(t6)
        neigh_dict4.extend(t7)
        e_lj += t8

    pool.close()
    pool.join()

    e_lj = e_lj / 2.0 # because of the double loop
    mem.e_lj = e_lj/N
    for i in range(N):
        mem.FX[i] = resultdict1[i]
        mem.FY[i] = resultdict2[i]
        mem.FZ[i] = resultdict3[i]
        mem.pcg_neigh_list[i][1] = neigh_dict1[i]
        mem.pcg_neigh_list[i][2] = neigh_dict2[i]
        mem.pcg_neigh_list[i][3] = neigh_dict3[i]
        mem.pcg_neigh_list[i][4] = neigh_dict4[i]

    return 0

# -----------------------------------------------------------------------------------------
def find_lj_force():
    N = mem.N
    X = mem.X
    Y = mem.Y
    Z = mem.Z
    TYPE = mem.TYPE
    # FX = mem.FX (dont do this! sync issues)

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    RCUT = mem.RCUT
    RCUT2 = mem.rcut_sqr

    rbond2 = mem.r_bond[0] * mem.r_bond[0]  # the index 0 has to be generalized

    # Virial = [0, 0, 0]
    e_lj = 0

    for i in range(N - 1):

        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        itype = TYPE[i]

        for j in range(i + 1, N):

            xj = X[j]
            yj = Y[j]
            zj = Z[j]
            jtype = TYPE[j]

            xij = xi - xj
            yij = yi - yj
            zij = zi - zj

            if (abs(xij) > half_L):
                if (xij < 0):
                    xij = xij + L
                else:
                    xij = xij - L

            if abs(xij) > RCUT:
                continue

            if (abs(yij) > half_B):
                if (yij < 0):
                    yij = yij + B
                else:
                    yij = yij - B

            if abs(yij) > RCUT:
                continue

            if (abs(zij) > half_H):
                if (zij < 0):
                    zij = zij + H
                else:
                    zij = zij - H

            if abs(zij) > RCUT:
                continue

            rsqr = xij * xij + yij * yij + zij * zij

            if (rsqr < RCUT2):

                if rsqr < rbond2:  # found a neighbor for pcg
                    mem.pcg_neigh_list[i][jtype] += 1
                    mem.pcg_neigh_list[j][itype] += 1

                irr = 1 / rsqr
                ir6 = irr * irr * irr
                ir12 = ir6 * ir6

                e_lj = e_lj + (mem.lj_A[itype][jtype] * ir12 - mem.lj_B[itype][jtype] * ir6)
                FORCE_MAG = (12 * mem.lj_A[itype][jtype] * ir12 - 6 * mem.lj_B[itype][jtype] * ir6) * irr

                mem.FX[i] = mem.FX[i] + FORCE_MAG * xij
                mem.FY[i] = mem.FY[i] + FORCE_MAG * yij
                mem.FZ[i] = mem.FZ[i] + FORCE_MAG * zij

                mem.FX[j] = mem.FX[j] - FORCE_MAG * xij
                mem.FY[j] = mem.FY[j] - FORCE_MAG * yij
                mem.FZ[j] = mem.FZ[j] - FORCE_MAG * zij

                # Virial[0] = Virial[0] - 0.5 * FORCE_MAG * xij *xij
                # Virial[1] = Virial[1] - 0.5 * FORCE_MAG * yij *yij
                # Virial[2] = Virial[2] - 0.5 * FORCE_MAG * zij *zij

    mem.e_lj = e_lj / mem.N


#    return Virial
# -----------------------------------------------------------------------------------------
def bond_force():
    nbonds = mem.NB
    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    e_bond = 0

    for i in range(nbonds):
        k = mem.bond_k[mem.BTYPE[i]]
        d0 = mem.bond_d0[mem.BTYPE[i]]

        ith_atom = mem.B1[i] - 1
        jth_atom = mem.B2[i] - 1

        delx = mem.X[ith_atom] - mem.X[jth_atom]
        dely = mem.Y[ith_atom] - mem.Y[jth_atom]
        delz = mem.Z[ith_atom] - mem.Z[jth_atom]

        if (abs(delx) > half_L):
            if (delx < 0):
                delx = delx + L
            else:
                delx = delx - L

        if (abs(dely) > half_B):
            if (dely < 0):
                dely = dely + B
            else:
                dely = dely - B

        if (abs(delz) > half_H):
            if (delz < 0):
                delz = delz + H
            else:
                delz = delz - H

        rsq = delx * delx + dely * dely + delz * delz
        r = sqrt(rsq)
        dr = r - d0
        rk = k * dr
        fbond = -2.0 * rk / r
        e_bond += rk * dr

        mem.FX[ith_atom] += fbond * delx
        mem.FY[ith_atom] += fbond * dely
        mem.FZ[ith_atom] += fbond * delz

        mem.FX[jth_atom] -= fbond * delx
        mem.FY[jth_atom] -= fbond * dely
        mem.FZ[jth_atom] -= fbond * delz
    mem.e_bond = e_bond / mem.N


# -----------------------------------------------------------------------------------------
def angle_force():
    nangles = mem.NA
    X = mem.X
    Y = mem.Y
    Z = mem.Z

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    e_angle = 0

    for i in range(nangles):

        i1 = mem.A1[i] - 1
        i2 = mem.A2[i] - 1
        i3 = mem.A3[i] - 1

        atype = mem.ATYPE[i]
        k = mem.ang_k[atype]
        t0 = mem.ang_t0[atype]

        # 1st bond

        delx1 = X[i1] - X[i2]
        dely1 = Y[i1] - Y[i2]
        delz1 = Z[i1] - Z[i2]

        if (abs(delx1) > half_L):
            if (delx1 < 0):
                delx1 = delx1 + L
            else:
                delx1 = delx1 - L

        if (abs(dely1) > half_B):
            if (dely1 < 0):
                dely1 = dely1 + B
            else:
                dely1 = dely1 - B

        if (abs(delz1) > half_H):
            if (delz1 < 0):
                delz1 = delz1 + H
            else:
                delz1 = delz1 - H

        rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1
        #r1 = sqrt(rsq1)

        # 2nd bond

        delx2 = X[i3] - X[i2]
        dely2 = Y[i3] - Y[i2]
        delz2 = Z[i3] - Z[i2]

        if (abs(delx2) > half_L):
            if (delx2 < 0):
                delx2 = delx2 + L
            else:
                delx2 = delx2 - L

        if (abs(dely2) > half_B):
            if (dely2 < 0):
                dely2 = dely2 + B
            else:
                dely2 = dely2 - B

        if (abs(delz2) > half_H):
            if (delz2 < 0):
                delz2 = delz2 + H
            else:
                delz2 = delz2 - H

        rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2
        #r2 = sqrt(rsq2)

        r1r2 = sqrt(rsq1 * rsq2)

        # angle (cos and sin)

        c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2
        #c = c / (r1 * r2)
        c /= r1r2

        if c > 1.0: c = 1.0
        if c < -1.0: c = -1.0

        s = sqrt(1.0 - c * c)
        if s < 0.001: s = 0.001
        s = 1.0 / s

        # ------ force & energy

        dtheta = acos(c) - t0
        tk = k * dtheta

        e_angle = e_angle + tk * dtheta

        a = -2.0 * tk * s
        a11 = a * c / rsq1
        #a12 = -a / (r1 * r2)
        a12 = -a / r1r2
        a22 = a * c / rsq2

        f1x = a11 * delx1 + a12 * delx2
        f1y = a11 * dely1 + a12 * dely2
        f1z = a11 * delz1 + a12 * delz2
        f3x = a22 * delx2 + a12 * delx1
        f3y = a22 * dely2 + a12 * dely1
        f3z = a22 * delz2 + a12 * delz1

        # apply force to each of 3 atoms

        mem.FX[i1] += f1x
        mem.FY[i1] += f1y
        mem.FZ[i1] += f1z

        mem.FX[i2] -= f1x + f3x
        mem.FY[i2] -= f1y + f3y
        mem.FZ[i2] -= f1z + f3z

        mem.FX[i3] += f3x
        mem.FY[i3] += f3y
        mem.FZ[i3] += f3z

    mem.e_angle = e_angle / mem.N


# ----------------------Thrombin -> E region interaction-----------------------------------
# -----------------------------------------------------------------------------------------
def bonding_force1():
    X = mem.X
    Y = mem.Y
    Z = mem.Z

    RCUT2 = mem.rcut_sqr

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    #dont call below function as it already called in lj_force module ;)
    #calc_pcg_neighbors()  # --------call the resetting function

    rbond2 = mem.r_bond[0] * mem.r_bond[0]  # the index 0 has to be generalized

    e_pcg = 0

    type3list = mem.type3list
    type4list = mem.type4list

    for i in type3list:

        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        itype = 3 #TYPE[i]

        for j in type4list:

            xj = X[j]
            yj = Y[j]
            zj = Z[j]
            jtype = 4

            Nn1 = mem.pcg_neigh_list[i][jtype]
            Nn2 = mem.pcg_neigh_list[j][itype]
            Nn = max(Nn1,Nn2)

            xij = xi - xj
            yij = yi - yj
            zij = zi - zj

            if abs(zij) > half_H:
                if zij < 0:
                    zij = zij + H
                else:
                    zij = zij - H

            if abs(xij) > half_L:
                if xij < 0:
                    xij = xij + L
                else:
                    xij = xij - L

            if abs(yij) > half_B:
                if yij < 0:
                    yij = yij + B
                else:
                    yij = yij - B

            rsqr = xij * xij + yij * yij + zij * zij

            if (rsqr < RCUT2):
                irr = 1 / rsqr
                r = sqrt(rsqr)
                ir = 1/r

                E_new1 = -mem.pcg_bond[itype][jtype] * ir
                Fnew1 = -mem.pcg_bond[itype][jtype] * irr *ir  # Attractive force

                if mem.pcg_morse_lj == 0:
                    ir4 = irr * irr
                    ir8 = ir4 * ir4
                    Fnew2 = (8 * mem.pcg_pair_A[itype][jtype] * ir8 - 4 * mem.pcg_pair_B[itype][jtype] * ir4) * irr # bonding force
                    E_new2 = (mem.pcg_pair_A[itype][jtype] * ir8 - mem.pcg_pair_B[itype][jtype] * ir4)
                elif mem.pcg_morse_lj == 1:
                    exponential = exp(-mem.alpha*(r - mem.pcg_d0[itype][jtype]))
                    oneminusexp2 = (1-exponential)**2
                    Fnew2 = 2*mem.alpha * mem.pcg_k2[itype][jtype] * exponential * ir * oneminusexp2
                    E_new2 = mem.pcg_k2[itype][jtype] * oneminusexp2
                else:
                    print "Error: unknown pcg_morse_lj style"
                    mem.error_flag = 1

                #Fnew3 = mem.pcg_bond[itype][jtype] * irr *ir   # Repulsive force
                E_new3 = -2 * E_new1
                Fnew3 = -2 * Fnew1  # Repulsive force

                theta1 = heavyside(0.9 - Nn)
                theta2 = heavyside((Nn - 0.9) / (1.1 - Nn))
                theta3 = heavyside(Nn - 1.9)
                theta4 = heavyside(rsqr - rbond2)
                theta5 = heavyside(rbond2 - rsqr)

                theta_block1 = (theta1 + theta2 * theta5 + theta3 * theta5)
                theta_block2 = theta2 * theta5
                theta_block3 = (theta3 + theta2 * theta4)

                e_pcg += theta_block1 * E_new1 + theta_block2 * E_new2 + theta_block3 * E_new3
                Fbonding = theta_block1 * Fnew1 + theta_block2 * Fnew2 + theta_block3 * Fnew3

                mem.FX[i] += Fbonding * xij
                mem.FY[i] += Fbonding * yij
                mem.FZ[i] += Fbonding * zij

                mem.FX[j] -= Fbonding * xij
                mem.FY[j] -= Fbonding * yij
                mem.FZ[j] -= Fbonding * zij
    mem.e_pcg += (e_pcg / mem.N)

# -----------------------------------------------------------------------------------------
# ------------Thrombin -> D region interaction  -----------------------------------
def bonding_force2():
    X = mem.X
    Y = mem.Y
    Z = mem.Z

    RCUT2 = mem.rcut_sqr

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    rbond2 = mem.r_bond[1] * mem.r_bond[1]  # the index 0 has to be generalized

    e_pcg = 0

    type2list = mem.type2list
    type4list = mem.type4list

    for i in type2list:

        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        itype = 2

        for j in type4list:

            xj = X[j]
            yj = Y[j]
            zj = Z[j]
            jtype = 4

            Nn3 = mem.pcg_neigh_list[j][3]  # type 3 neighbors around j (4)

            if Nn3 == 1: # work only thrombins (4) attached to E region (3)

                xij = xi - xj
                yij = yi - yj
                zij = zi - zj

                if abs(zij) > half_H:
                    if zij < 0:
                        zij = zij + H
                    else:
                        zij = zij - H

                if abs(xij) > half_L:
                    if xij < 0:
                        xij = xij + L
                    else:
                        xij = xij - L

                if abs(yij) > half_B:
                    if yij < 0:
                        yij = yij + B
                    else:
                        yij = yij - B

                rsqr = xij * xij + yij * yij + zij * zij

                if rsqr < RCUT2:
                    Nn24 = mem.pcg_neigh_list[i][4]  # type 4 neighbors around i (2)
                    Nn42 = mem.pcg_neigh_list[j][2]  # type 2 neighbors around j (4)
                    irr = 1 / rsqr
                    r = sqrt(rsqr)
                    ir = 1/r
                    E_new1 = -mem.pcg_bond[itype][jtype] * ir
                    Fnew1 = -mem.pcg_bond[itype][jtype] * irr *ir  # Attractive force

                    if mem.pcg_morse_lj == 0:
                        ir4 = irr * irr
                        ir8 = ir4 * ir4
                        E_new2 = (mem.pcg_pair_A[itype][jtype] * ir8 - mem.pcg_pair_B[itype][jtype] * ir4)
                        Fnew2 = (8 * mem.pcg_pair_A[itype][jtype] * ir8 - 4 * mem.pcg_pair_B[itype][jtype] * ir4) * irr # bonding force
                    elif mem.pcg_morse_lj == 1:
                        exponential = exp(-mem.alpha*(r - mem.pcg_d0[itype][jtype]))
                        oneminusexp2 = (1 - exponential) ** 2
                        Fnew2 = 2 * mem.alpha * mem.pcg_k2[itype][jtype] * exponential * ir * oneminusexp2
                        E_new2 = mem.pcg_k2[itype][jtype] * oneminusexp2

                    #Fnew3 = mem.pcg_bond[itype][jtype] * irr *ir   # Repulsive force
                    E_new3 = -2 * E_new1
                    Fnew3 = -2 * Fnew1  # Repulsive force

                    theta1 = heavyside(0.9 - Nn42) * heavyside(0.9 - Nn24)
                    theta2 = (heavyside((Nn42 - 0.9) / (1.1 - Nn42)) + heavyside((Nn42 - 1.9) / (2.1 - Nn42))) * heavyside((Nn24 - 0.9) / (1.1 - Nn24))
                    theta3 = max(heavyside(Nn24 - 1.9), heavyside(Nn42 - 2.9))
                    theta4 = heavyside(rsqr - rbond2)
                    theta5 = heavyside(rbond2 - rsqr)

                    theta_block1 = (theta1 + theta2 * theta5 + theta3 * theta5)
                    theta_block2 = theta2 * theta5
                    theta_block3 = (theta3 + theta2 * theta4)

                    e_pcg += theta_block1 * E_new1 + theta_block2 * E_new2 + theta_block3 * E_new3
                    Fbonding = theta_block1 * Fnew1 + theta_block2 * Fnew2 + theta_block3 * Fnew3

                    mem.FX[i] += Fbonding * xij
                    mem.FY[i] += Fbonding * yij
                    mem.FZ[i] += Fbonding * zij

                    mem.FX[j] -= Fbonding * xij
                    mem.FY[j] -= Fbonding * yij
                    mem.FZ[j] -= Fbonding * zij
    mem.e_pcg += (e_pcg / mem.N)

# -----------------------------------------------------------------------------------------
# ------------E -> D region interaction  -----------------------------------
def bonding_force3():
    rbond2 = mem.r_bond[1] * mem.r_bond[1]  # the index 0 has to be generalized

    e_pcg = 0

    ntype3 = mem.ntype3
    type3list = mem.type3list
    i1list = [-1, -1]
    i2list = [-1, -1]
    j1list = -1
    j2list = -1
    r1sqr = [-1, -1]
    r2sqr = [-1, -1]

    for m1 in range(0, ntype3-1, 1):
        m = type3list[m1]
        i1list[0] = m - 4
        i1list[1] = m + 4
        j1list = m
        for n1 in range(m1 + 1, ntype3, 1):
            n = type3list[n1]
            i2list[0] = n - 4
            i2list[1] = n + 4
            j2list = n

            connection = 0 # do these molecules have a bond?
            connection1 = 0
            connection2 = 0

            j = j2list
            tmp = 0
            for i in i1list:
                rsqr = get_rsqr(i,j)
                r1sqr[tmp] = rsqr
                tmp += 1
                if rsqr < 2*rbond2:
                    connection += 1
                    connection1 += 1

            j = j1list
            tmp = 0
            for i in i2list:
                rsqr = get_rsqr(i,j)
                r2sqr[tmp] = rsqr
                tmp += 1
                if rsqr < 2*rbond2:
                    connection += 1
                    connection2 += 1

            if connection == 1:
                if connection1 == 1:
                    j = j1list
                    if r2sqr[0] > r2sqr[1]:
                        i = i2list[1]
                    else:
                        i = i2list[0]

                    (rsqr,xij,yij,zij) = get_rsqr_all(i,j)

                    if rsqr > rbond2: # ask them to attract
                        irr = 1 / rsqr
                        r = sqrt(rsqr)
                        ir = 1 / r
                        E_new1 = -mem.pcg_bond[3][4] * ir
                        Fnew1 = -mem.pcg_bond[3][4] * irr * ir  # Attractive force

                        e_pcg += E_new1
                        Fbonding = Fnew1

                        mem.FX[i] += Fbonding * xij
                        mem.FY[i] += Fbonding * yij
                        mem.FZ[i] += Fbonding * zij

                        mem.FX[j] -= Fbonding * xij
                        mem.FY[j] -= Fbonding * yij
                        mem.FZ[j] -= Fbonding * zij
                else: #connection2 == 1
                    j = j2list
                    if r1sqr[0] > r1sqr[1]:
                        i = i1list[1]
                    else:
                        i = i1list[0]

                    (rsqr, xij, yij, zij) = get_rsqr_all(i, j)

                    if rsqr > rbond2:  # ask them to attract
                        irr = 1 / rsqr
                        r = sqrt(rsqr)
                        ir = 1 / r
                        E_new1 = -mem.pcg_bond[3][4] * ir
                        Fnew1 = -mem.pcg_bond[3][4] * irr * ir  # Attractive force

                        e_pcg += E_new1
                        Fbonding = Fnew1

                        mem.FX[i] += Fbonding * xij
                        mem.FY[i] += Fbonding * yij
                        mem.FZ[i] += Fbonding * zij

                        mem.FX[j] -= Fbonding * xij
                        mem.FY[j] -= Fbonding * yij
                        mem.FZ[j] -= Fbonding * zij

    mem.e_pcg += (e_pcg / mem.N)

# -----------------------------------------------------------------------------------------
def calc_pcg_neighbors():
    N = mem.N
    ntypes = mem.NTYPES
    X = mem.X
    Y = mem.Y
    Z = mem.Z
    TYPE = mem.TYPE

    xlo = mem.xlo
    xhi = mem.xhi
    ylo = mem.ylo
    yhi = mem.yhi
    zlo = mem.zlo
    zhi = mem.zhi

    L = xhi - xlo
    B = yhi - ylo
    H = zhi - zlo
    half_L = L / 2.0
    half_B = B / 2.0
    half_H = H / 2.0

    rbond2 = mem.r_bond[0] * mem.r_bond[0]  # the index 0 has to be generalized

    for i in range(N):
        for j in range(ntypes):
            mem.pcg_neigh_list[i][j] = 0

    for i in range(N - 1):
        xi = X[i]
        yi = Y[i]
        zi = Z[i]
        itype = TYPE[i]

        for j in range(i + 1, N):
            xj = X[j]
            yj = Y[j]
            zj = Z[j]
            jtype = TYPE[j]

            xij = xi - xj
            yij = yi - yj
            zij = zi - zj

            if abs(zij) > half_H:
                if zij < 0:
                    zij = zij + H
                else:
                    zij = zij - H

            if abs(xij) > half_L:
                if xij < 0:
                    xij = xij + L
                else:
                    xij = xij - L

            if abs(yij) > half_B:
                if yij < 0:
                    yij = yij + B
                else:
                    yij = yij - B

            rsqr = xij * xij + yij * yij + zij * zij

            if rsqr < rbond2:  # found a neighbor
                mem.pcg_neigh_list[i][jtype] += 1
                mem.pcg_neigh_list[j][itype] += 1


# -----------------------------------------------------------------------------------------
def heavyside(val):
    if val < 0:
        return 0
    elif val == 0:
        return 0.5
    else:
        return 1

# -----------------------------------------------------------------------------------------

def get_rsqr(i, j):

    xij = mem.X[i] - mem.X[j]
    yij = mem.Y[i] - mem.Y[j]
    zij = mem.Z[i] - mem.Z[j]

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    if abs(zij) > half_H:
        if zij < 0:
            zij = zij + H
        else:
            zij = zij - H

    if abs(xij) > half_L:
        if xij < 0:
            xij = xij + L
        else:
            xij = xij - L

    if abs(yij) > half_B:
        if yij < 0:
            yij = yij + B
        else:
            yij = yij - B

    rsqr = xij * xij + yij * yij + zij * zij
    return rsqr

def get_rsqr_all(i, j):

    xij = mem.X[i] - mem.X[j]
    yij = mem.Y[i] - mem.Y[j]
    zij = mem.Z[i] - mem.Z[j]

    L = mem.L
    B = mem.B
    H = mem.H
    half_L = mem.half_L
    half_B = mem.half_B
    half_H = mem.half_H

    if abs(zij) > half_H:
        if zij < 0:
            zij = zij + H
        else:
            zij = zij - H

    if abs(xij) > half_L:
        if xij < 0:
            xij = xij + L
        else:
            xij = xij - L

    if abs(yij) > half_B:
        if yij < 0:
            yij = yij + B
        else:
            yij = yij - B

    rsqr = xij * xij + yij * yij + zij * zij
    return rsqr, xij, yij, zij