import F07_memory as mem


def FIND_DISPLACEMENT():
    DT = mem.DT
    DT_SQR = DT * DT
    HALF_DT_SQR = 0.5 * DT_SQR
    N = mem.N

    for i in range(N):
        mem.X[i] = mem.X[i] + DT * mem.VX[i] + HALF_DT_SQR * mem.ACCX[i]
        mem.Y[i] = mem.Y[i] + DT * mem.VY[i] + HALF_DT_SQR * mem.ACCY[i]
        mem.Z[i] = mem.Z[i] + DT * mem.VZ[i] + HALF_DT_SQR * mem.ACCZ[i]
