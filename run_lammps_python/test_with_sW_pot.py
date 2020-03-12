import numpy as np
import subprocess


def create_sw_file(m):
    fileID = open("run_data/force_mW.sw","w");
    fileID.write("# comment\n# comment\n\n");
    fileID.write("mW mW mW %.4f %.4f %.4f %.4f %.4f %.4f\n" % (m[0],m[1],m[2],m[3],m[4],m[5]));
    fileID.write("         %.4f %.4f %.4f %.4f %.4f\n" % (m[6],m[7],m[8],m[9],m[10]));
    fileID.close();

def create_batch_file():
    fileID = open("run_data/runme.bat","w");
    fileID.write("echo 0 > program_status.txt\n");
    fileID.write("lmp_serial -in film_run.in\n");
    fileID.write("echo 1 > program_status.txt\nexit\n");
    fileID.close();
    
def create_vapor_run_script(dt,Tref,equ_run,prod_run):
    fileID = open("run_data/vapor_run.in","w");
    fileID.write("units real\n");
    fileID.write("boundary p p p\n");
    fileID.write("atom_style atomic\n\n");
    fileID.write("pair_style sw\n");
    fileID.write("read_data vapor.data\n");
    fileID.write("pair_coeff * * force_mW.sw mW\n");
    fileID.write("timestep %.3f\n" % (dt));
    fileID.write("neighbor 2.0 bin\n");
    fileID.write("neigh_modify delay 0 every 1\n\n");
    fileID.write("fix 2 all nvt temp %.3f %.3f %.3f\n" % (Tref,Tref,dt*100));
    fileID.write("thermo 100\n");
    fileID.write("thermo_style one\n");
    fileID.write("\nlog dontuselog_vapor.log \n");
    fileID.write("\nmin_style cg\n");
    fileID.write("minimize 1.0e-4 1.0e-6 100 1000\n");
    fileID.write("\nvelocity all create %.3f %d\n" % (Tref,np.round_(np.random.rand(1)*36000)));
    fileID.write("run 10\n");
    fileID.write("velocity all scale %.3f\n" % Tref);
    fileID.write("\nrun %d\n" % equ_run);
    fileID.write("\nthermo 1\n");
    fileID.write("thermo_style custom step temp ke pe etotal press density enthalpy\n");
    fileID.write("thermo_modify norm yes flush yes temp thermo_temp press thermo_press\n");
    fileID.write("\nlog log_vapor.log\n");
    fileID.write("\nrun %d\n" % prod_run);
    fileID.close();
    
def create_liquid_run_script(dt,Tref,equ_run,prod_run):
    fileID = open("run_data/liquid_run.in","w");
    fileID.write("units real\n");
    fileID.write("boundary p p p\n");
    fileID.write("atom_style atomic\n\n");
    fileID.write("pair_style sw\n");
    fileID.write("read_data liquid.data\n");
    fileID.write("pair_coeff * * force_mW.sw mW\n");
    fileID.write("timestep %.3f\n" % (dt));
    fileID.write("neighbor 2.0 bin\n");
    fileID.write("neigh_modify delay 0 every 1\n\n");
    fileID.write("fix 2 all npt temp %.3f %.3f %.3f iso 1.0009 1.0009 %.3f\n" % (Tref,Tref,dt*100,dt*500));
    fileID.write("thermo 100\n");
    fileID.write("thermo_style one\n");
    fileID.write("\nlog dontuselog_liquid.log \n");
    fileID.write("\nmin_style cg\n");
    fileID.write("minimize 1.0e-4 1.0e-6 100 1000\n");
    fileID.write("\nvelocity all create %.3f %d\n" % (Tref,np.round_(np.random.rand(1)*36000)));
    fileID.write("run 10\n");
    fileID.write("velocity all scale %.3f\n" % Tref);
    fileID.write("\nrun %d\n" % equ_run);
    fileID.write("\nthermo 1\n");
    fileID.write("thermo_style custom step temp ke pe etotal press density enthalpy\n");
    fileID.write("thermo_modify norm yes flush yes temp thermo_temp press thermo_press\n");
    fileID.write("\nlog log_liquid.log\n");
    fileID.write("\nrun %d\n" % prod_run);
    fileID.close();
    
def create_film_run_script(dt,Tref,equ_run,prod_run):
    print("Module for creating lammps run script");
    fileID = open("run_data/film_run.in","w");
    fileID.write("units real\n");
    fileID.write("boundary p p p\n");
    fileID.write("atom_style atomic\n\n");
    fileID.write("pair_style sw\n");
    fileID.write("read_data film.data\n");
    fileID.write("pair_coeff * * force_mW.sw mW\n");
    fileID.write("timestep %.3f\n" % (dt));
    fileID.write("neighbor 2.0 bin\n");
    fileID.write("neigh_modify delay 0 every 1\n\n");
    fileID.write("fix 2 all nvt temp %.3f %.3f %.3f\n" % (Tref,Tref,dt*100));
    fileID.write("variable gamm equal (101325*1e-10*1000)*lz*0.5*(pzz-(pxx+pyy)/2.0)\n");
    fileID.write("thermo 100\n");
    fileID.write("thermo_style one\n");
    fileID.write("\nlog dontuselog_film.log \n");
    fileID.write("\nmin_style cg\n");
    fileID.write("minimize 1.0e-4 1.0e-6 100 1000\n");
    fileID.write("\nvelocity all create %.3f %d\n" % (Tref,np.round_(np.random.rand(1)*36000)));
    fileID.write("run 10\n");
    fileID.write("velocity all scale %.3f\n" % Tref);
    fileID.write("\nrun %d\n" % equ_run);
    fileID.write("\nthermo 1\n");
    fileID.write("thermo_style custom step temp ke pe etotal press v_gamm\n");
    fileID.write("thermo_modify norm yes flush yes\n");
    fileID.write("\nlog log_film.log\n");
    fileID.write("\nrun %d\n" % prod_run);
    fileID.close();

def run_all_lammps():
    cmd = 'call C:/Sumith_Projects/Python_codes/run_lammps_with_python/sssss.bat'
    failure = subprocess.call(cmd, shell=True)
    if failure:
        print ('Execution of "%s" failed!\n' % cmd)

def main():
    print("main function");
    dt = 10.0;
    Tref = 373.15;
    equ_run = 500; #50000;
    prod_run = 500; #75000;
    parameters = [8.385138, 2.751924, 2.495835, 10.649491, 1.192579, -0.302891, 4.190227, 1.981567, 5.004732, 0.132740, 0.175862];

    create_batch_file();
    
    create_sw_file(parameters);
    create_film_run_script(dt,Tref,equ_run,prod_run);
    create_liquid_run_script(dt,Tref,equ_run,prod_run);
    create_vapor_run_script(dt,Tref,equ_run,prod_run);
    run_all_lammps();
    
    
if __name__ == "__main__":
    main()
