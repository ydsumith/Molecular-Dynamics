import F01_read_inputs as var1

def FIND_DISPLACEMENT():    
    DT = var1.DT
    DT_SQR = DT*DT
    HALF_DT_SQR = 0.5 * DT_SQR
    NARGON = var1.NARGON
    
    for i in range(NARGON):
        var1.ARGON[i][0] = var1.ARGON[i][0] + DT * var1.VELO[i][0] + HALF_DT_SQR * var1.ACCEL[i][0]
        var1.ARGON[i][1] = var1.ARGON[i][1] + DT * var1.VELO[i][1] + HALF_DT_SQR * var1.ACCEL[i][1]
        var1.ARGON[i][2] = var1.ARGON[i][2] + DT * var1.VELO[i][2] + HALF_DT_SQR * var1.ACCEL[i][2]
