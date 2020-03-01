import F01_read_inputs as init

def WRITE_VIRIAL(virial,kinetic,box_VOL,STEP):
    N = init.NARGON
    if STEP == 1:
        with open('Pressure.res','w') as f:
            f.write('VirXX\tVirYY\tVirZZ\tPXX\tPYY\tPZZ\tPavg\n')
            f.close()
    else:
        with open('Pressure.res','a') as f:
            #print "VirXX = %s, VirYY = %s, VirZZ = %s" % (virial[0],virial[1],virial[2])
            TOBAR = 16.6054
            PXX = 2*TOBAR*(kinetic[0] - virial[0])/box_VOL
            PYY = 2*TOBAR*(kinetic[1] - virial[1])/box_VOL
            PZZ = 2*TOBAR*(kinetic[2] - virial[2])/box_VOL
            P = (PXX + PYY + PZZ)/3.0
            f.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (virial[0],virial[1],virial[2],PXX,PYY,PZZ,P))
            #print "PXX = %s, PYY = %s, PZZ = %s, Pavg = %s" % (PXX,PYY,PZZ,P)
        f.close()
    
def WRITE_XYZ(STEP):
    N = init.NARGON
    DT = init.DT
    curr = DT*float(STEP)
    if STEP == 1:
        with open('traj_pyth.xyz','w') as f:     f.close()
    with open('traj_pyth.xyz','a') as f:
        f.write('%s\n' % N)
        f.write('@ %s ps\n' % curr)
        for i in range(N):
            f.write('%-2s %15.6f %15.6f %15.6f\n' % ('Ar',init.ARGON[i][0]*10,init.ARGON[i][1]*10,init.ARGON[i][2]*10) )
    f.close()
def WRITE_RES(KE,PE,T,STEP):
    if STEP == 1:
        with open('energy.res','w') as f:
            f.write('Step\tT[K]\tKE\tPE\tTotalE\n')
        f.close()
    else:
        with open('energy.res','a') as f:
            TE = float(KE)+float(PE)
            f.write( '%s\t%s\t%s\t%s\t%s\n' % (STEP,T,KE,PE,TE))
        f.close()
