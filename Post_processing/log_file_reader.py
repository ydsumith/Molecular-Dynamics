# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:05:38 2020

@author: ydsum
"""

def get_current_folder():
    cur_dir1 = input("Enter the location of this python file : ")
    cur_dir2 = cur_dir1.replace("\\","/");
    return cur_dir1,cur_dir2;

def vapor_loger(T,prod_run,dir2,r):
    file_to_open = dir2 + "/log_" + str(T) + "_" + str(r) + ".log";
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

def liquid_loger(T,prod_run,dir2,r):
    file_to_open = dir2 + "/log_" + str(T) + "_" + str(r) + ".log";
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

def main():
    typer = 0; # 0 for liquid, 1 for vapor
    prod_run = 50000;
    
    dir1,dir2 = get_current_folder();
    Tref = [273.15, 283.15,293.15,303.15,313.15,323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,413.15,423.15,433.15,443.15];
    rnd = [28831, 28777, 33365, 4923, 31803, 9001, 5157, 1598, 9629, 16290, 17240, 30198, 6627, 31351, 18654, 18544, 15907, 26478, 21025, 15585];
    print (dir1)
    print (dir2)
    if typer == 0:
        print("Now running liquid module");
        for T in Tref:
            for r in rnd:
                Hliq, density = liquid_loger(T,prod_run,dir2,r);
                print("T = %.2f, Hliq = %.3f, density = %.3f" % (T,Hliq, density));
    else:
        print("Now running vapor module");
        for T in Tref:
            for r in rnd:
                Hvap = vapor_loger(T,prod_run,dir2,r);
                print("T = %.2f, Hvap = %.3f" % (T, Hvap));
        
if __name__ == "__main__":
    main()
