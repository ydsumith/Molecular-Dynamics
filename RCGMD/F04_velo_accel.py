import F07_memory as mem


def FIND_VELOCITY():
    DT = mem.DT
    HALF_DT = 0.5 * DT
    N = mem.N
    for i in range(N):
        mem.VX[i] = mem.VX[i] + HALF_DT * mem.ACCX[i]
        mem.VY[i] = mem.VY[i] + HALF_DT * mem.ACCY[i]
        mem.VZ[i] = mem.VZ[i] + HALF_DT * mem.ACCZ[i]


def CLEAR_FORCE():
    N = mem.N
    for i in range(N):
        mem.FX[i] = 0
        mem.FY[i] = 0
        mem.FZ[i] = 0
    mem.e_angle = 0
    mem.e_bond = 0
    mem.e_lj = 0
    mem.e_ke = 0
    mem.e_total = 0
    mem.e_pcg = 0

    mem.pcg_n_D_4 = 0
    mem.pcg_n_E_4 = 0

    ntypes = mem.NTYPES
    for i in range(N):
        for j in range(ntypes):
            mem.pcg_neigh_list[i][j] = 0


def FIND_ACCELERATION():
    N = mem.N
    mass = mem.MASS

    for i in range(N):
        mem.ACCX[i] = mem.FX[i] / mass[i]
        mem.ACCY[i] = mem.FY[i] / mass[i]
        mem.ACCZ[i] = mem.FZ[i] / mass[i]
