"""
Molecular Dynamics program for estimating Virial Pressure for
Argon system in a cube
"""
import F01_read_inputs as init
import F02_displacement as disp
import F03_common as misc
import F04_velo_accel as va
import F05_find_force as fo
import F06_output as out

global pressure
init.read_inputs()
init.read_argon()
mass = 39.9
PRESS_COUNT = 0
box_VOL = init.box_L * init.box_B * init.box_H

for STEP_VAL in range(1,init.TOT_STEPS):
    pressure = False
    disp.FIND_DISPLACEMENT()
    misc.CHECK_AND_UPDATE()
    va.FIND_VELOCITY() # 1st pass
    if(STEP_VAL % 50 == 0):
        pressure = True
        PRESS_COUNT = PRESS_COUNT + 1
    va.CLEAR_FORCE()
    RET = fo.FIND_FORCE()
    EWW = RET[0]
    EWPT = RET[1]
    virial = RET[2]
    va.FIND_ACCELERATION(mass)
    va.FIND_VELOCITY() # 2nd pass
    if(STEP_VAL % init.NST_COMM == 0): misc.COM_REMOVAL()
    RET1 = misc.FIND_KE_TEMP(mass)
    T_BULK = RET1[0]
    kinetic = RET1[1]
    KE = RET1[2]
    if(init.THERMOSTAT != 0):        misc.VELOCITY_SCALING(T_BULK)
    if(STEP_VAL % init.PRESS_INTVL == 0): out.WRITE_VIRIAL(virial,kinetic,box_VOL,STEP_VAL)
    E_TOTAL = EWW + EWPT
    out.WRITE_RES(KE,E_TOTAL,T_BULK,STEP_VAL)
    print "Exec =",float(STEP_VAL)*100/init.TOT_STEPS," %,PE = ", EWW, "T = ",T_BULK
    if(STEP_VAL % init.nstxout == 0):  out.WRITE_XYZ(STEP_VAL)
#End of main time stepping
print "Program completed"
