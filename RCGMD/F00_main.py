"""
PCG-MD (Polymer Coarse Grain MD)
copyright 2017
Author: Dr. Sumith Yesudasan
Molecular Dynamics program for polymerization of Fibrinogen

internal units: L(nm), energy (kJ/mol), time (ps)
"""
import F01_read_inputs as init
import F02_displacement as disp
import F03_common as misc
import F04_velo_accel as va
import F05_find_force as fo
import F06_output as out
import F07_memory as mem
import time

if __name__ == '__main__':
    start_time = time.time()

    mem.nprocs = 2 #-------write code to fetch this automatically
    print "using %s processors" % mem.nprocs
    mem.error_flag = 0

    init.read_inputs()
    init.read_data()
    PRESS_COUNT = 0
    box_VOL = (mem.xhi - mem.xlo) * (mem.yhi - mem.ylo) * (mem.zhi - mem.zlo)

    if mem.INIT_MOMENT == 1:    misc.INIT_VELO()  # initialize with velocities

    if mem.r_bond[0] > mem.RCUT:
        print "Error: rbond > rcut"
        mem.error_flag = 1

    for STEP_VAL in range(1, mem.TOT_STEPS):
        if mem.error_flag == 0:

            disp.FIND_DISPLACEMENT()

            misc.CHECK_AND_UPDATE()

            va.FIND_VELOCITY()  # 1st pass

            if (STEP_VAL % 50 == 0):
                pressure = True
                PRESS_COUNT = PRESS_COUNT + 1

            va.CLEAR_FORCE()

            #fo.find_lj_force() # serial version
            ret_val = fo.find_parall_lj_force() # parallel version

            #out.write_force("1_lj_parallel")

            fo.bond_force()
            fo.angle_force()

            #out.write_force("2_bef_phs1")

            if STEP_VAL >= mem.pcg_thro_start:
                fo.bonding_force1() # Thrombin - E region
            #out.write_force("3_aft_phs1")

            if STEP_VAL >= mem.pcg_poly_start:
                fo.bonding_force2() # Thrombin - D region
                fo.bonding_force3()
            #out.write_force("4_aft_phs2")

            va.FIND_ACCELERATION()
            va.FIND_VELOCITY()  # 2nd pass

            if STEP_VAL % mem.NST_COMM == 0: misc.COM_REMOVAL()

            misc.FIND_KE_TEMP()

            if (mem.THERMOSTAT != 0):        misc.VELOCITY_SCALING()

            # if(STEP_VAL % mem.PRESS_INTVL == 0): out.WRITE_VIRIAL(virial,kinetic,box_VOL,STEP_VAL)
            mem.e_total = mem.e_lj + mem.e_bond + mem.e_angle + mem.e_ke + mem.e_pcg

            out.loggit("Exec = %s , Tot_E = %s, T= %s" % (float(STEP_VAL) * 100 / mem.TOT_STEPS, mem.e_total, mem.T_BULK))

            if (STEP_VAL % mem.nstxout == 0):
                out.WRITE_XYZ(STEP_VAL)
                out.WRITE_RES(STEP_VAL)
                misc.happy_ending(STEP_VAL)
        else:
            print "Exiting due to error, check log file for details"
            break
    # End of main time stepping
    print "Program completed"
    print("--- %s seconds ---" % (time.time() - start_time))


