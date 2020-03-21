import tkinter as tk
from tkinter import filedialog
import numpy as np
import re

def get_current_folder():
    root = tk.Tk()
    root.withdraw()

    in_file = filedialog.askopenfilename()
    out_file = filedialog.asksaveasfilename(initialdir = in_file ,\
               title = "Select file to save",\
               filetypes = (("Text files","*.txt"),("all files","*.*")))
    return in_file, out_file;

def main():    
    in_file, out_file = get_current_folder();
    print ("File to read is \n", in_file);
    print ("File to write is \n", out_file);
    Temp = np.array([273.15, 283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15,\
            353.15, 363.15, 373.15, 383.15, 393.15, 403.15, 413.15, 423.15,\
            433.15,443.15]);
    lTemp = len(Temp);
    n = 5;
    res = np.zeros([lTemp,n]); # 6th coumn reserved for counter
    
    with open(in_file) as file_in:        
        for line in file_in:
            temp1 = re.findall( r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line );
            l1 = len(temp1);
            cur_temperature = float(temp1[0]);
            for j in range(lTemp): # finds the location of the data in the main
                if cur_temperature == Temp[j]:
                    cur_index = j;
                    break;
            for i in range(l1):
                res[cur_index][i] = res[cur_index][i] + float(temp1[i]);
                
            res[cur_index][n-1] += 1; # nth or last column                
                
    for i in range(lTemp):
        if res[i][n-1] != 0:
            for j in range(n-1):
                res[i][j] = res[i][j] / res[i][n-1];
    
    with open(out_file,"w") as file_out:
        for i in range(lTemp):
            for j in range(n-1):
                file_out.write("%.5f\t" % res[i][j]);
            file_out.write("\n");
    

        
if __name__ == "__main__":
    main()
