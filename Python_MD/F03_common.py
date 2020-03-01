import F01_read_inputs as var1
from math import sqrt

def CHECK_AND_UPDATE():
    NARG = var1.NARGON
    L = var1.box_L
    B = var1.box_B
    H = var1.box_H
    PBC = var1.PBC
    
    for i in range(NARG):
        if(var1.ARGON[i][0] < 0): var1.ARGON[i][0] = var1.ARGON[i][0] + L
        if(var1.ARGON[i][0] > L): var1.ARGON[i][0] = var1.ARGON[i][0] - L
        if(var1.ARGON[i][1] < 0): var1.ARGON[i][1] = var1.ARGON[i][1] + B;
        if(var1.ARGON[i][1] > B): var1.ARGON[i][1] = var1.ARGON[i][1] - B;

    if PBC == 1:
        for i in range(NARG):
            if(var1.ARGON[i][2] < 0): var1.ARGON[i][2] = var1.ARGON[i][2] + H
            if(var1.ARGON[i][2] > H): var1.ARGON[i][2] = var1.ARGON[i][2] - H
def COM_REMOVAL():
    VXSUM = VYSUM = VZSUM = 0
    NARG = var1.NARGON
    VXSUM = sum(var1.VELO[:][0])/NARG
    VYSUM = sum(var1.VELO[:][1])/NARG
    VZSUM = sum(var1.VELO[:][2])/NARG

    var1.VELO[:][0] = [x-VXSUM for x in var1.VELO[:][0]]
    var1.VELO[:][1] = [x-VYSUM for x in var1.VELO[:][1]]

    if(var1.PBC == 1):
        var1.VELO[:][2] = [x-VZSUM for x in var1.VELO[:][2]]

def FIND_KE_TEMP(mass):
    kinetic_pressure_X = kinetic_pressure_Y = kinetic_pressure_Z = 0
    BOLTZ = 0.008314462175 #(kJ/mol/K)
    NARG = var1.NARGON
    KE = 0
    for i in range(NARG):
        vxx = var1.VELO[i][0] * var1.VELO[i][0]
        vyy = var1.VELO[i][1] * var1.VELO[i][1]
        vzz = var1.VELO[i][2] * var1.VELO[i][2]

        KE = KE + (vxx + vyy + vzz)
        kinetic_pressure_X = kinetic_pressure_X + 0.5 * vxx* mass
        kinetic_pressure_Y = kinetic_pressure_Y + 0.5 * vyy* mass
        kinetic_pressure_Z = kinetic_pressure_Z + 0.5 * vzz* mass
        
    KE = KE * 0.5 * mass
    T_BULK = KE * 2 / (3 * NARG * BOLTZ)
    KE = KE / NARG
    kinetic = [kinetic_pressure_X, kinetic_pressure_Y, kinetic_pressure_Z]
    return [T_BULK, kinetic,KE]
def VELOCITY_SCALING(T_BULK):
    if(var1.THERMOSTAT == 1): scale = sqrt(1+var1.DT*((var1.T_REF/T_BULK)-1)/var1.TAU) #Berendsen
    if(var1.THERMOSTAT == 2):
        scale =  sqrt(var1.T_REF/T_BULK) #Velocity scaling
    NARG = var1.NARGON
    for i in range(NARG):
        var1.VELO[i][0] = var1.VELO[i][0] * scale
        var1.VELO[i][1] = var1.VELO[i][1] * scale
        var1.VELO[i][2] = var1.VELO[i][2] * scale
