import F01_read_inputs as init

def FIND_FORCE():
    N = init.NARGON
    L = init.box_L
    B = init.box_B
    H = init.box_H
    half_L = L/2.0
    half_B = B/2.0
    half_H = H/2.0
    
    sig12 = pow(0.34,12)
    sig6 = pow(0.34,6)
    eps = 1.005841
    RCUT = init.RCUT
    RCUT2 = init.RCUT*init.RCUT
    PBC = init.PBC
    LJ1 = (6*sig12*pow(RCUT,-12) - 3*sig6*pow(RCUT,-6))/RCUT2
    LJ2 = -7*sig12*pow(RCUT,-12) + 4*sig6*pow(RCUT,-6)
    LJ3 = (48*eps*sig12*pow(RCUT,-12) - 24*eps*sig6*pow(RCUT,-6))/RCUT2
    EWW = 0
    Virial = [0,0,0]
    for i in range(N-1):
        xi = init.ARGON[i][0]
        yi = init.ARGON[i][1]
        zi = init.ARGON[i][2]
        for j in range(i+1,N):
            xj = init.ARGON[j][0]
            yj = init.ARGON[j][1]
            zj = init.ARGON[j][2]

            xij = xi - xj
            yij = yi - yj
            zij = zi - zj
            
            if(PBC == 1):
                if(abs(zij) > half_H):
                    if(zij < 0): zij = zij + H
                    else: zij = zij - H
            
            if(abs(xij) > half_L):
                if(xij < 0): xij = xij + L
                else: xij = xij - L
            
            if(abs(yij) > half_B):
                if(yij < 0): yij = yij + B
                else: yij = yij - B

            rsqr = xij*xij + yij*yij + zij*zij
            if(rsqr < RCUT2):
                irr = 1/rsqr
                ir6 = irr * irr * irr
                ir12 = ir6 * ir6
                EWW = EWW + 4*eps*(sig12*ir12 - sig6*ir6 + LJ1*rsqr + LJ2)
                FORCE_MAG = (48*eps*sig12*ir12 - 24*eps*sig6*ir6)*irr - LJ3
                init.FORCE[i][0] = init.FORCE[i][0] + FORCE_MAG * xij
                init.FORCE[i][1] = init.FORCE[i][1] + FORCE_MAG * yij
                init.FORCE[i][2] = init.FORCE[i][2] + FORCE_MAG * zij

                init.FORCE[j][0] = init.FORCE[j][0] - FORCE_MAG * xij
                init.FORCE[j][1] = init.FORCE[j][1] - FORCE_MAG * yij
                init.FORCE[j][2] = init.FORCE[j][2] - FORCE_MAG * zij

                Virial[0] = Virial[0] - 0.5 * FORCE_MAG * xij *xij
                Virial[1] = Virial[1] - 0.5 * FORCE_MAG * yij *yij
                Virial[2] = Virial[2] - 0.5 * FORCE_MAG * zij *zij
    EWPT = 0
    RET = [EWW/N, EWPT, Virial]
    return RET
