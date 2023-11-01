#include <iostream>
#include <fstream>
#include <thread>
#include <math.h>
#include <string>
#include <chrono>

//#define outname "out(5)_0.002z_0.1m_s10^5.csv"
#define outname "out.csv"

using namespace std;


/*
int main() {
    double f[502], w[502], fN[502], wN[502], z[502];
    double S, t, dt, dz;
    int i;
    dz = 0.002;
    dt = 0.00000001;
    S = 1000000;
    i = 0;
    t = 0;

    ofstream fstr;
    fstr.open("1.txt");
    while (i < 501)
    {
        z[i] = i * dz;
        f[i] = 0;
        w[i] = 0;
        fN[i] = 0;
        wN[i] = 0;

        i = i + 1;
    }


beginning: i = 1;
    while (i < 500)
    {
        wN[i] = w[i] + dt * (((w[i + 1] - 2 * w[i] + w[i - 1]) / (dz * dz)) - f[i] * ((w[i + 1] - w[i - 1]) / (2 * dz)) + S * (z[i] - 0.5));
        fN[i] = f[i] + 100 * dt * (((f[i + 1] - 2 * f[i] + f[i - 1]) / (dz * dz)) - w[i]);
        i = i + 1;
    }

    fN[0] = 0;
    fN[500] = 0;
    wN[0] = 2 * fN[1] / (dz * dz);
    wN[500] = 2 * fN[499] / (dz * dz);
    i = 0;
    while (i < 501)
    {
        w[i] = wN[i];
        f[i] = fN[i];
        i = i + 1;
    }

    //fstr << t << " " << f[10] << "\n";


    t = t + dt;
    if (t < 1) { goto beginning; }



    i = 0;
    while (i < 501)
    {
        fstr << z[i] << "  " << f[i] << "\n";
        i = i + 1;
    }

}
*/

void threadPrint(long int &size_mass_x, ofstream &ofile, long double &t, long double &dx, long double* F, long double* W, bool &exit){
    while (true) {
        string input;
        cin >> input;
        
        if (input == "p") {
            for (long int k = 0; k < size_mass_x; k+=1)
              ofile << t << ";" << k * dx << ";" << F[k] << ";" << W[k] << endl;
            ofile << endl;
            cout << "Print!" << endl;
            continue;
        }

        if (input == "e") {
            exit = 1;
            break;
        }
    }
}

long int main(){
    cout << "Start" << endl << endl;
    
    // настройки
    long double dx = 0.001;
    long double dt = powl(10, -8);
    long double print_step = 0.001;
    long double s = powl(10, 0);       // !
    long double m = 1.00;

    long double print_pos = 0.0;
    long double t = 0;
    long double delta = 0.0;
    long double delta_max = 0.0;
    long double max_x = 1.0;
    chrono::steady_clock::time_point time1 = chrono::steady_clock::now();
    chrono::steady_clock::time_point time2;
    int mod = 2;
    bool exit = 0;

    // инициализация основных переменных
    long int size_mass_x = ceil(fabs(max_x / dx)) + 1;

    // инициализация массивов
    long double* F_1 = new long double[size_mass_x];
    long double* F_2 = new long double[size_mass_x];

    long double* W_1 = new long double[size_mass_x];
    long double* W_2 = new long double[size_mass_x];
    
    // ввод S
    //cout << "Input S = "; cin >> s;

    // вывод 
    ofstream ofile(outname);
    ofile << fixed;
    ofile.precision(12);
    ofile << "S =;" << s << ";;" << endl << "dx =;" << dx << ";;" << endl << "dt =;" << dt << ";;" << endl << "m =;" << m << ";;" << endl << "t;x;F;W" << endl;

    // вывод 
    ofstream ofile_time("out_time_0.2_0.4.csv");
    ofile_time << fixed;
    ofile_time.precision(12);
    ofile_time << "t;F_0.2;F_0.4" << endl;


    // начальные условия
    for (long int i = 0; i < size_mass_x; i++) {
        F_1[i] = 0.0;
        W_1[i] = 0.0;
    }

    // поток для управления
    thread thr(threadPrint, ref(size_mass_x), ref(ofile), ref(t), ref(dx), ref(F_2), ref(W_2), ref(exit));

    //расчет
    while (exit != 1) {
        delta_max = 0;

        if (mod == 2) {
            for (long int i = 1; i < size_mass_x - 1; i++) {
                F_2[i] = F_1[i] + dt/m*((F_1[i+1] - 2*F_1[i] + F_1[i-1])/dx/dx - W_1[i]);
                W_2[i] = W_1[i] + dt*((W_1[i-1] - 2*W_1[i] + W_1[i+1])/dx/dx - (W_1[i+1] - W_1[i-1])/(2*dx)*F_1[i] + s*(i*dx - 0.5));

                delta = F_1[i] - F_2[i];
                if (abs(delta) > abs(delta_max)) delta_max = delta;
            }
            W_2[0] = 2*F_2[1]/dx/dx;
            W_2[size_mass_x-1] = 2*F_2[size_mass_x-2]/dx/dx;

            mod = 1;
        }
        else {
            for (long int i = 1; i < size_mass_x - 1; i++) {
                F_1[i] = F_2[i] + dt/m*((F_2[i + 1] - 2 * F_2[i] + F_2[i - 1]) / dx / dx - W_2[i]);
                W_1[i] = W_2[i] + dt * ((W_2[i - 1] - 2 * W_2[i] + W_2[i + 1]) / dx / dx - (W_2[i + 1] - W_2[i - 1]) / (dx*2) * F_2[i] + s * (i * dx - 0.5));
                
                delta = F_1[i] - F_2[i];
                if (abs(delta) > abs(delta_max)) delta_max = delta;
            }
            W_1[0] = 2*F_1[1] / dx / dx;
            W_1[size_mass_x - 1] = 2*F_1[size_mass_x - 2] / dx / dx;

            mod = 2;
        }

        if (t >= print_pos) {
            time2 = chrono::steady_clock::now();
            cout << "delta = " << delta_max << "\tt = " << t << "\tn = " << chrono::duration_cast<chrono::milliseconds>(time2-time1).count() << " ms" << endl;
            time1 = time2;
            print_pos += print_step;


            // вывод 
            ofstream ofile("out_0.01z_1m_s10^9_" + std::to_string(t) + ".csv");
            ofile << fixed;
            ofile.precision(12);
            ofile << "S =;" << s << ";;" << endl << "dx =;" << dx << ";;" << endl << "dt =;" << dt << ";;" << endl << "m =;" << m << ";;" << endl << "t;x;F;W" << endl;

            for (long int k = 0; k < size_mass_x; k += 1)
                ofile << t << ";" << k * dx << ";" << F_2[k] << ";" << W_2[k] << endl;
            ofile.close();

            ofile_time << t << ";" << F_1[(int)(0.2 / dx)] << ";" << F_1[(int)(0.4 / dx)] << endl;
        }

        t += dt;
    }

    thr.join();

    ofile.close();
    ofile_time.close();

    cout << endl << "Ende." << endl;
    return 0;
}

/*
 void threadPrint(long int &size_mass_x, ofstream &ofile, bigD &t, long double &dx, bigD* F, bigD* W, bool &exit){
    while (true) {
        string input;
        cin >> input;
        
        if (input == "p") {
            for (long int k = 0; k < size_mass_x; k++)
              ofile << t << ";" << k * dx << ";" << F[k] << ";" << W[k] << endl;
            ofile << endl;
            cout << "PRint!" << endl;
            continue;
        }

        if (input == "e") {
            exit = 1;
            break;
        }
    }
}

long int main(){
    cout << "Start" << endl << endl;
    
    // настройки
    long double dx = 0.002;
    bigD dt = powl(10, -10);
    long double print_step = 0.0001;
    long double s = 10000.0;       // !
    long double m = 0.01;

    bigD print_pos = 0.0;
    bigD t = 0;
    bigD error = 0.0;
    bigD error_max = 0.0;
    long double max_x = 1.0;
    int mod = 2;
    bool exit = 0;

    // инициализация основных переменных
    long int size_mass_x = ceil(fabs(max_x / dx)) + 1;

    // инициализация массивов
    bigD* F_1 = new bigD[size_mass_x];
    bigD* F_2 = new bigD[size_mass_x];

    bigD* W_1 = new bigD[size_mass_x];
    bigD* W_2 = new bigD[size_mass_x];
    
    // вывод 
    ofstream ofile(outname);
    ofile << fixed;
    ofile.precision(12);
    ofile << "t;x;F;W" << endl;

    // начальные условия
    for (long int i = 0; i < size_mass_x; i++) {
        F_1[i] = 0.0;
        W_1[i] = 0.0;
    }

    // поток для управления
    thread thr(threadPrint, ref(size_mass_x), ref(ofile), ref(t), ref(dx), ref(F_2), ref(W_2), ref(exit));

    //расчет
    while (exit != 1) {
        error_max = 0;

        if (mod == 2) {
            for (long int i = 1; i < size_mass_x - 1; i++) {
                F_2[i] = F_1[i] + dt/m*((F_1[i+1] - 2*F_1[i] + F_1[i-1])/dx/dx - W_1[i]);
                W_2[i] = W_1[i] + dt*((W_1[i-1] - 2*W_1[i] + W_1[i+1])/dx/dx - (W_1[i+1] - W_1[i-1])/dx*F_1[i] + s*(i*dx - 0.5));

                error = abs(F_1[i] - F_2[i]);
                if (error > error_max) error_max = error;
            }
            W_2[0] = 2*F_2[1]/dx/dx;
            W_2[size_mass_x-1] = 2*F_2[size_mass_x-2]/dx/dx;

            mod = 1;
        }
        else {
            for (long int i = 1; i < size_mass_x - 1; i++) {
                F_1[i] = F_2[i] + dt/m*((F_2[i + 1] - 2 * F_2[i] + F_2[i - 1]) / dx / dx - W_2[i]);
                W_1[i] = W_2[i] + dt * ((W_2[i - 1] - 2 * W_2[i] + W_2[i + 1]) / dx / dx - (W_2[i + 1] - W_2[i - 1]) / dx * F_2[i] + s * (i * dx - 0.5));
                
                error = abs(F_1[i] - F_2[i]);
                if (error > error_max) error_max = error;
            }
            W_1[0] = 2*F_1[1] / dx / dx;
            W_1[size_mass_x - 1] = 2*F_1[size_mass_x - 2] / dx / dx;

            mod = 2;
        }

        if (t >= print_pos) {
            cout << "error = " << error_max << "\t t = " << t << endl;
            print_pos += print_step;
        }

        t += dt;
    }

    thr.join();

    cout << endl << "Ende." << endl;
    return 0;
}
*/
