import F01_read_inputs as var1

def FIND_VELOCITY():
    DT = var1.DT
    HALF_DT = 0.5 * DT
    NARGON = var1.NARGON
    
    for i in range(NARGON):
        var1.VELO[i][0] = var1.VELO[i][0] + HALF_DT * var1.ACCEL[i][0]
        var1.VELO[i][1] = var1.VELO[i][1] + HALF_DT * var1.ACCEL[i][1]
        var1.VELO[i][2] = var1.VELO[i][2] + HALF_DT * var1.ACCEL[i][2]

def CLEAR_FORCE():
    NARGON = var1.NARGON
    
    for i in range(NARGON):
        var1.FORCE[i][:] = [0,0,0]

def FIND_ACCELERATION(mass):
    NARGON = var1.NARGON
    
    for i in range(NARGON):
        var1.ACCEL[i][0] = var1.FORCE[i][0]/mass
        var1.ACCEL[i][1] = var1.FORCE[i][1]/mass
        var1.ACCEL[i][2] = var1.FORCE[i][2]/mass
