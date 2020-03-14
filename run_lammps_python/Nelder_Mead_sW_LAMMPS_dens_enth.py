# Author - Prof. Sumith Yesudasan
# Year 2020

import numpy as np
import subprocess
import math
import time

def create_sw_file(m):
    fileID = open("run_data/force_mW.sw","w");
    fileID.write("# comment\n# comment\n\n");
    fileID.write("mW mW mW %.4f %.4f %.4f %.4f %.4f %.4f\n" % (m[0],m[1],m[2],m[3],m[4],m[5]));
    fileID.write("         %.4f %.4f %.4f %.4f %.4f\n" % (m[6],m[7],m[8],m[9],m[10]));
    fileID.close();

def create_batch_file():
    fileID = open("run_data/runme.bat","w");
    fileID.write("echo 0 > program_status.txt\n");
    fileID.write("mpiexec -np 6 lmp_mpi -in liquid_run_273.in\n");
    fileID.write("mpiexec -np 6 lmp_mpi -in liquid_run_323.in\n");
    fileID.write("mpiexec -np 6 lmp_mpi -in liquid_run_373.in\n");
    fileID.write("mpiexec -np 6 lmp_mpi -in liquid_run_423.in\n");
    fileID.write("mpiexec -np 2 lmp_mpi -in vapor_run_273.in\n");
    fileID.write("mpiexec -np 2 lmp_mpi -in vapor_run_323.in\n");
    fileID.write("mpiexec -np 2 lmp_mpi -in vapor_run_373.in\n");
    fileID.write("mpiexec -np 2 lmp_mpi -in vapor_run_423.in\n");
    fileID.write("echo 1 > program_status.txt\nexit\n");
    fileID.close();
    #---
    fileID = open("master.bat","w");
    fileID.write("cd C:/Users/Sumith.Yesudasan/PROJECTS/LAMMPSrun/for_paper_6/Automated_Nelder_mW/run_data\n");
    fileID.write("start C:/Users/Sumith.Yesudasan/PROJECTS/LAMMPSrun/for_paper_6/Automated_Nelder_mW/run_data/runme.bat\n");
    fileID.close();
    #---
    fileID = open("results/results.txt","w");
    fileID.close();
    #---
    fileID = open("results/parameters.txt","w");    
    fileID.close();
    
def create_vapor_run_script(dt,Tref,equ_run,prod_run):
    file_to_open = "run_data/vapor_run_" + str(round(Tref)) + ".in";
    fileID = open(file_to_open,"w");
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
    fileID.write("\nlog log_vapor_%d.log\n" % round(Tref));
    fileID.write("\nrun %d\n" % prod_run);
    fileID.close();
    
def create_liquid_run_script(dt,Tref,equ_run,prod_run):
    file_to_open = "run_data/liquid_run_" + str(round(Tref)) + ".in";
    fileID = open(file_to_open,"w");
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
    fileID.write("\nlog log_liquid_%d.log\n" % round(Tref));
    fileID.write("\nrun %d\n" % prod_run);
    fileID.close();


def run_all_lammps():
    cmd = 'call C:/Users/Sumith.Yesudasan/PROJECTS/LAMMPSrun/for_paper_6/Automated_Nelder_mW/master.bat'
    failure = subprocess.call(cmd, shell=True)
    if failure:
        print ('Execution of "%s" failed!\n' % cmd)
    #---
    #print("Running the lammps programs, wait...");
    waitforme = 0;
    while waitforme < 1:
        time.sleep(1.0);
        with open("run_data/program_status.txt") as fileID:
            for line in fileID:
                temp = line.split(" ");
                waitforme = int(temp[0]);
    #print("ran all lammps");

def vapor_loger(T,prod_run):
    file_to_open = "run_data/log_vapor_" + str(round(T)) + ".log";
    with open(file_to_open) as fileID:
        iteration = 0; ncuttan = 0; Hvap = 0;
        for line in fileID:
            iteration += 1;
            if iteration > 5 and iteration < prod_run:
                ncuttan += 1;
                temp = line.split();
                Hvap = Hvap + float(temp[7]);
        Hvap = Hvap / ncuttan;
    return Hvap;

def liquid_loger(T,prod_run):
    file_to_open = "run_data/log_liquid_" + str(round(T)) + ".log";
    with open(file_to_open) as fileID:
        iteration = 0; ncuttan = 0; Hliq = 0; density = 0;
        for line in fileID:
            iteration += 1;
            if iteration > 5 and iteration < prod_run:
                ncuttan += 1;
                temp = line.split();
                Hliq = Hliq + float(temp[7]);
                density = density + float(temp[6]);
        Hliq = Hliq / ncuttan;
        density = 1000*density / ncuttan;
    return Hliq, density;
    
def estimate_costs(prod_run,T):
    Hvap = np.zeros((4,1));
    Hliq = np.zeros((4,1));
    density = np.zeros((4,1));
    for i in range(4):
        Hvap[i] = vapor_loger(T[i],prod_run);
        Hliq[i], density[i] = liquid_loger(T[i],prod_run);
    #--estimate cost value
    ref_density = [999.79, 987.99, 958.35, 917];
    ref_enthalpy = [10.768, 10.256, 9.715, 9.101];
    dH = Hvap - Hliq;
    dens_error = 0; enth_error = 0;
    for i in range(4):
        dens_error = dens_error + math.sqrt(((density[i]-ref_density[i])**2)/(ref_density[i]**2));
        enth_error = enth_error + math.sqrt(((dH[i] - ref_enthalpy[i])**2)/(ref_enthalpy[i]**2));
    #--weights for each error will be assigned
    curr_error = 0.50*dens_error + 0.50*enth_error;
    curr_error = curr_error *100;
    fileID = open("results/results.txt","a");
    fileID.write("\nCurrent Error: %.2f\n" % (curr_error));
    for i in range(4):
        fileID.write("T = %.2f, Density = %.2f, Enthalpy = %.2f\n" % (T[i],density[i],dH[i]));
    fileID.close();
    return curr_error;
   
    
def cost_estimate(x,n,dt,equ_run,prod_run,typee,T):
    create_sw_file(x);
    for Tref in T:
        create_liquid_run_script(dt,Tref,equ_run,prod_run);
        create_vapor_run_script(dt,Tref,equ_run,prod_run);
    
    run_all_lammps();
    fileID = open("results/results.txt","a");
    if typee == 1:  fileID.write("Initial simplex:\n");
    if typee == 2:  fileID.write("Reflection:\n");
    if typee == 3:  fileID.write("Expansion:\n");
    if typee == 4:  fileID.write("Contraction:\n");
    if typee == 5:  fileID.write("Shrink:\n");
    fileID.close();
    curr_error = estimate_costs(prod_run,T);
    return curr_error;

def create_initial_simplex(x,n,dh):
    x_mat = np.zeros((n+1,n));
    for i in range(n):
        x_mat[0][i] = x[i];
    for j in range(1,n+1): # row navigation
        for i in range(n): # column nav
            if i == j-1:
                x_mat[j][i] = x[i] * dh;
            else:
                x_mat[j][i] = x[i];
    return x_mat

def sort_mats(cost_mat,x_mat,n):
    ind = np.argsort(cost_mat, axis=0);
    cost_mat = np.take_along_axis(cost_mat, ind, axis=0);
    x_mat = np.take_along_axis(x_mat, ind, axis=0);
    return cost_mat,x_mat;

def write_parameters(x_mat,iteration,n,cur_cost):
    fileID = open("results/parameters.txt","a");  
    fileID.write("\nCurrent iteration = %d, cost = %.3f,  writing best parameters below\n" % (iteration,cur_cost));
    for j in range(n):
        fileID.write("%.6f " % (x_mat[0][j]));
    fileID.write("\n");
    fileID.close();


def main():
    start_time = time.time()
    n = 11; dh = 1.5; # 50 %
    alpha = 1.0; gamma = 2.0; rho = 0.5; sigma = 0.5;
    term_tol = 1e-4; MAX_ITER = 200; cur_iter = 0;
    stddev = 1000; # initial val
    T = [273.15,323.15, 373.15, 423.15];

    cost_mat = np.zeros((n+1,1));
    x_0 = np.zeros((n,1));

    dt = 5.0;
    equ_run = 50; #25000;
    prod_run = 50; #25000;
    x = [8.412576, 2.743976, 2.442725, 10.653902, 1.225309, -0.240518, 4.334005, 1.990876, 5.006381, 0.192317, 0.207676];

    create_batch_file();
    
    # step 1: estimate the cost of the current simplex
    x_mat = create_initial_simplex(x,n,dh);
    print("Here is the initial simplex\n",x_mat);
    for i in range(n+1): # from x_i to x_n+1
        cost_mat[i] = cost_estimate(x_mat[i][:],n,dt,equ_run,prod_run,1,T);
        
    while stddev > term_tol and cur_iter < MAX_ITER: # termination criteria
        print("Current iteration = ", cur_iter);
        # sorting
        (cost_mat,x_mat) = sort_mats(cost_mat,x_mat,n);
        write_parameters(x_mat,cur_iter,n,cost_mat[0]);
        stddev = np.std(cost_mat);        
        # step 2: centroid
        x_0 = np.sum(x_mat[0:n][:],axis=0)/n;
        # step 3: reflection
        x_n1 = x_mat[n][:]; # x_(n+1)
        x_r = x_0 + alpha*(x_0-x_n1);
        cost_r = cost_estimate(x_r,n,dt,equ_run,prod_run,2,T);
        if cost_mat[0] <= cost_r and cost_r < cost_mat[n]:
            x_mat[n][:] = x_r; # Obtaining a new simplex
            cost_mat[n] = cost_r;
        else:
            # step 4: Expansion
            if cost_r < cost_mat[0]:
                x_e = x_0 + gamma*(x_r - x_0);
                cost_e = cost_estimate(x_e,n,dt,equ_run,prod_run,3,T);
                if cost_e < cost_r:
                    x_mat[n][:] = x_e; # Obtaining a new simplex
                    cost_mat[n] = cost_e;
                    cur_iter += 1;
                    continue;
                else:
                    x_mat[n][:] = x_r;  # Obtaining a new simplex
                    cost_mat[n] = cost_r;
                    cur_iter += 1;
                    continue;
            #step 5 Contraction
            else: 
                x_c = x_0 + rho*(x_mat[n][:] - x_0);
                cost_c = cost_estimate(x_c,n,dt,equ_run,prod_run,4,T);
                if cost_c < cost_mat[n]:
                    x_mat[n][:] = x_c;  # Obtaining a new simplex
                    cost_mat[n] = cost_c;
                    cur_iter += 1;
                    continue;
                #step 6 Shrink
                else:
                    for i in range(1,n+1): # from x_2 to x_n+1
                        x_mat[i][:] = x_mat[0][:] + sigma*(x_mat[i][:] - x_mat[0][:]);
                        cost_mat[i] = cost_estimate(x_mat[i][:],n,dt,equ_run,prod_run,5,T);
        cur_iter += 1;
    # while ends here, terminate the whole business
    print("\nTerminating..\nstddev = ",stddev,", tolerance = ",term_tol);
    print("Total iterations = ",cur_iter);
    print("Best x values = ", x_mat[0][:]);
    print("Cost = \n", cost_mat);
    end_time = time.time();
    print("Elapsed time: ",end_time - start_time);
    
if __name__ == "__main__":
    main()
