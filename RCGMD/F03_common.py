import F07_memory as mem
import F06_output as out

from math import sqrt
from random import *


# ---------------------------------------------------------------------------------
def CHECK_AND_UPDATE():
    N = mem.N
    xlo = mem.xlo
    xhi = mem.xhi
    ylo = mem.ylo
    yhi = mem.yhi
    zlo = mem.zlo
    zhi = mem.zhi
    L = xhi - xlo
    B = yhi - ylo
    H = zhi - zlo

    for i in range(N):
        if (mem.X[i] < xlo): mem.X[i] = mem.X[i] + L
        if (mem.X[i] > xhi): mem.X[i] = mem.X[i] - L
        if (mem.Y[i] < ylo): mem.Y[i] = mem.Y[i] + B
        if (mem.Y[i] > yhi): mem.Y[i] = mem.Y[i] - B
        if (mem.Z[i] < zlo): mem.Z[i] = mem.Z[i] + H
        if (mem.Z[i] > zhi): mem.Z[i] = mem.Z[i] - H


# ---------------------------------------------------------------------------------
def COM_REMOVAL():
    VXSUM = VYSUM = VZSUM = 0
    N = mem.N
    VXSUM = sum(mem.VX[:]) / N
    VYSUM = sum(mem.VY[:]) / N
    VZSUM = sum(mem.VZ[:]) / N
    mem.VX[:] = [i - VXSUM for i in mem.VX[:]]
    mem.VY[:] = [i - VYSUM for i in mem.VY[:]]
    mem.VZ[:] = [i - VZSUM for i in mem.VZ[:]]


# ---------------------------------------------------------------------------------
def FIND_KE_TEMP():
    kinetic_pressure_X = kinetic_pressure_Y = kinetic_pressure_Z = 0
    BOLTZ = 0.008314462175  # (kJ/mol/K)
    # BOLTZ = 0.001985898889 #(kcal/mol/K)
    N = mem.N
    KE = 0
    MASS = mem.MASS
    # log.loggit('writing velocities---')
    for i in range(N):
        vxx = mem.VX[i] * mem.VX[i]
        vyy = mem.VY[i] * mem.VY[i]
        vzz = mem.VZ[i] * mem.VZ[i]
        # log.loggit('%s, %s, %s, %s' % (mem.VX[i],mem.VY[i],mem.VZ[i],mem.MASS[i]))
        KE = KE + (vxx + vyy + vzz) * MASS[i]
        # kinetic_pressure_X = kinetic_pressure_X + 0.5 * vxx* MASS[i]
        # kinetic_pressure_Y = kinetic_pressure_Y + 0.5 * vyy* MASS[i]
        # kinetic_pressure_Z = kinetic_pressure_Z + 0.5 * vzz* MASS[i]

    KE = KE * 0.5
    mem.T_BULK = KE * 2 / (3 * N * BOLTZ)

    mem.e_ke = KE / N
    # kinetic = [kinetic_pressure_X, kinetic_pressure_Y, kinetic_pressure_Z]


# ---------------------------------------------------------------------------------
def VELOCITY_SCALING():
    if (mem.THERMOSTAT == 1): scale = sqrt(1 + mem.DT * ((mem.T_REF / mem.T_BULK) - 1) / mem.TAU)  # Berendsen
    if (mem.THERMOSTAT == 2):
        scale = sqrt(mem.T_REF / mem.T_BULK)  # Velocity scaling
    N = mem.N
    for i in range(N):
        mem.VX[i] = mem.VX[i] * scale
        mem.VY[i] = mem.VY[i] * scale
        mem.VZ[i] = mem.VZ[i] * scale


# ---------------------------------------------------------------------------------

def INIT_VELO():
    N = mem.N
    minvel = -0.05
    maxvel = 0.05
    for i in range(N):
        mem.VX[i] = uniform(minvel, maxvel)
        mem.VY[i] = uniform(minvel, maxvel)
        mem.VZ[i] = uniform(minvel, maxvel)

# ---------------------------------------------------------------------------------


def happy_ending(step):

    N = mem.N
    TYPE = mem.TYPE

    for i in range(N):
        itype = TYPE[i]
        if itype == 4:
            mem.pcg_n_E_4 += mem.pcg_neigh_list[i][3]
            mem.pcg_n_D_4 += mem.pcg_neigh_list[i][2]
    out.pcg_out("%s\t%s\t%s" % (step,mem.pcg_n_E_4,mem.pcg_n_D_4))