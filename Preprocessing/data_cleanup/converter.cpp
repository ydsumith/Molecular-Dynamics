#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

FILE *dumpfile;

int main()
{
    string linebuffer, filename, delim;
    ifstream in;

    double exclude_distance = 1.0; // Ang

    int length, pos, counter;
    int error_flag, N = 0;
    double *xdata;
    double *ydata;
    double *zdata;

    filename = "data_to_C.data";

    in.open(filename.c_str());
    if(!in) {
        cout << "Error: unable to open file " <<filename<< endl;
        error_flag = 1;
        return -1;
    } else cout << "Reading file " <<filename <<"  pass 1..."<< endl;

    // PASS 1
    while ( getline(in, linebuffer)) {
        N = N + 1;
    }
    in.close();
    cout << "Total lines: " << N;
    //PASS 2
    xdata = (double *)malloc((N+1)*sizeof(double));
    ydata = (double *)malloc((N+1)*sizeof(double));
    zdata = (double *)malloc((N+1)*sizeof(double));

    counter = 0;
    in.open(filename.c_str());
    while ( getline(in, linebuffer)) {
        if (linebuffer.find (" ") != string::npos)
            delim = " ";
        else if (linebuffer.find ("\t") != string::npos)
            delim = "\t";
        else {
            error_flag = 1;
            return -1;
        }

        length = linebuffer.length();
        pos = linebuffer.find(delim);
        xdata[counter] = atof(linebuffer.substr(0, pos).c_str());
        linebuffer.erase(0, pos+1);

        pos = linebuffer.find(delim);
        ydata[counter] = atof(linebuffer.substr(0, pos).c_str());
        linebuffer.erase(0, pos+1);

        length = linebuffer.length();
        zdata[counter] = atof(linebuffer.substr(0,length).c_str());

        counter++;
    }
    in.close();

    //#%% - Minimize
    cout<<"\ncurrently processing the minimization, please wait...\n";
    double rx, ry, rz, r, rsp, nx, ny, nz;
    counter = 0;
    float percent;
    int moder = round(N / 200);
    for(int i = 0; i<N; i++) {
        if(i % moder == 0 ){
            percent = 100*(static_cast<double>(i)/static_cast<double>(N);
            cout << percent << " % done, " << i << " out of "<<N << "\n";
        }
        int j = i+1;
        while (j < N) {
            rx = xdata[i]-xdata[j];
            ry = ydata[i]-ydata[j];
            rz = zdata[i]-zdata[j];
            if (abs(rx) <= exclude_distance) {
                if (abs(rz) <= exclude_distance) {
                    if (abs(ry) <= exclude_distance) {
                        r = sqrt(rx*rx + ry*ry + rz*rz);
                        if (r <= exclude_distance) {
                            nx = rx/r;
                            ny = ry/r;
                            nz = rz/r;
                            rsp = (exclude_distance - r) / 2;
                            xdata[i] = xdata[i] + rsp * nx;
                            ydata[i] = ydata[i] + rsp * ny;
                            zdata[i] = zdata[i] + rsp * nz;

                            xdata[j] = xdata[j] - rsp * nx;
                            ydata[j] = ydata[j] - rsp * ny;
                            zdata[j] = zdata[j] - rsp * nz;
                            counter += 1;
                        }
                    }
                }
            }
            j += 1;
        }
    }
    cout<<"\nTotal "<< counter <<" modifications made\n";

    dumpfile = fopen("data_from_C.data","w");
    for (int i=0; i<N ; i++){
        fprintf(dumpfile,"%f %f %f\n", xdata[i], ydata[i], zdata[i]);
    }

    fflush(dumpfile);
    fclose(dumpfile);

    free(xdata);
    free(ydata);
    free(zdata);
    cout << "\nProgram completed\n";
         return 0;
}
