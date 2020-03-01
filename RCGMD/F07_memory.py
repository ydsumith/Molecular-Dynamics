# memory holder for global variables

nprocs = 0

error_flag = 0

path = ''
DT = 0
TOT_STEPS = 0
PRESS_INTVL = 0
SLAB_THICK = 0
THERMOSTAT = 0
T_REF = 0
sigma_ARPT = 0
eps_ARPT = 0
rc_wall = 0
RCRIT = 0
INIT_MOMENT = 0
RCUT = 0
NST_COMM = 0
nstxout = 0
WALL_MODEL = 0
top_w = 0
bot_w = 0
TAU = 0
massvalue = []
massvalue.append(0)

T_BULK = 0

N = 0  # No: of Atoms
NB = 0  # No: of Bonds
NA = 0  # No: of Angles
NTYPES = 0  # No: of atom types
NBTYPES = 0  # No: of bond types
NATYPES = 0  # No: of angle types
xlo = xhi = ylo = yhi = zlo = zhi = 0

L = 0
B = 0
H = 0
half_L = 0
half_B = 0
half_H = 0
rcut_sqr = 0

MOLID = []  # Matrix storing molids
TYPE = []  # Matrix with type of atoms
CHARGE = []
X = []
Y = []
Z = []

BTYPE = []
B1 = []
B2 = []

ATYPE = []
A1 = []
A2 = []
A3 = []

ACCX = []
ACCY = []
ACCZ = []
FX = []
FY = []
FZ = []
VX = []
VY = []
VZ = []
MASS = []

pair_eps = []
pair_sig = []
bond_k = []
bond_d0 = []
ang_k = []
ang_t0 = []

lj_A = []
lj_B = []

e_angle = 0
e_bond = 0
e_lj = 0
e_ke = 0
e_total = 0
e_pcg = 0

pcg_bond = []   # bond coeffs for Fnew1 and Fnew3
pcg_pair_A = []  # pair coeff for Fnew2
pcg_pair_B = []  # pair coeff for Fnew2

pcg_morse_lj = -1 # put 0 for LJ and 1 for morse potential in the run_parameters file

pcg_d0 = []
pcg_k2 = []
alpha = 3.0 # spread of the morse potential

pcg_thro_start = 0 # starting thrombin cleaving
pcg_poly_start = 0 # starting factor XIIIa attach

pcg_n_E_4 = 0 # no: of E (type 3) around type 4
pcg_n_D_4 = 0 # no: of D (type 2) around type 4

pcg_neigh_list = []  # store neighbors (Nn) of type i

r_bond = []  # r_bond store 3.5,future2, future3..
pcg_val = []  # valency 1,future2, future3...
ipcg = []   # ipcg store 3,future2, future3..
jpcg = []   # jpcg store 4,future2, future3..

ntype2 = 0
ntype3 = 0
ntype4 = 0

type2list = []
type3list = []
type4list = []
