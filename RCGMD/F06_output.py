import F07_memory as mem

def pcg_out(content):
    with open('out_pcg_log.log', 'a') as f:
        f.write('%s\n' % (content))
    f.close()

def WRITE_VIRIAL(virial, kinetic, box_VOL, STEP):
    N = mem.N
    if STEP == 1:
        with open('out_Pressure.res', 'w') as f:
            f.write('VirXX\tVirYY\tVirZZ\tPXX\tPYY\tPZZ\tPavg\n')
        f.close()
    else:
        with open('out_Pressure.res', 'a') as f:
            # print "VirXX = %s, VirYY = %s, VirZZ = %s" % (virial[0],virial[1],virial[2])
            TOBAR = 16.6054
            PXX = 2 * TOBAR * (kinetic[0] - virial[0]) / box_VOL
            PYY = 2 * TOBAR * (kinetic[1] - virial[1]) / box_VOL
            PZZ = 2 * TOBAR * (kinetic[2] - virial[2]) / box_VOL
            P = (PXX + PYY + PZZ) / 3.0
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (virial[0], virial[1], virial[2], PXX, PYY, PZZ, P))
            # print "PXX = %s, PYY = %s, PZZ = %s, Pavg = %s" % (PXX,PYY,PZZ,P)
        f.close()


def WRITE_XYZ(STEP):
    N = mem.N
    DT = mem.DT
    curr = DT * float(STEP)
    with open('out_traj_pyth.mol', 'a') as f:
        f.write('ITEM: TIMESTEP\n')
        f.write('%s\n' % STEP)
        f.write('ITEM: NUMBER OF ATOMS\n')
        f.write('%s\n' % N)
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('%s %s\n' % (mem.xlo, mem.xhi))
        f.write('%s %s\n' % (mem.ylo, mem.yhi))
        f.write('%s %s\n' % (mem.zlo, mem.zhi))
        f.write('ITEM: ATOMS id mol type x y z\n')
        for i in range(N):
            f.write('%d %d %d %.6f %.6f %.6f\n' % (i + 1, mem.MOLID[i], mem.TYPE[i], mem.X[i], mem.Y[i], mem.Z[i]))
    f.close()


def WRITE_RES(step):
    with open('out_energy.res', 'a') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (mem.DT*step, mem.e_total, mem.T_BULK, mem.e_ke, mem.e_lj, mem.e_bond, mem.e_angle))
    f.close()


def loggit(content):
    with open('out_log.log', 'a') as f:
        f.write('%s\n' % (content))
    f.close()

def write_force(fileext):
    N = mem.N
    FX = mem.FX
    FY = mem.FY
    FZ = mem.FZ
    filename = "out_" + fileext + ".res"
    with open(filename, 'w') as f:
        for i in range(N):
            f.write('%s, %s, %s, %s\n' % (i+1, FX[i], FY[i], FZ[i]))
    f.close()