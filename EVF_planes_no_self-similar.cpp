#include <iostream>
#include <sstream>
#include <iomanip>
#include "fstream"
#include "cmath"

using namespace std;

int main() {

    // начальные параметры
    const long double dz = 0.02;
    const long double dr = 0.05;
    const long double R = 10;
    const long double S = powl(10, 0);

    // счет на установление
    long double dt = powl(10, -5);
    long double t_W = 0;
    long double t_PSI = 0;
    long double W_max = 0, PSI_max = 0;

    // вывод
    const unsigned int filename_precision = 5;
    const long double time_step = powl(10, -5);
    long double last_time_print = 0;

    // точность
    const long double eps_W = powl(10, -9);
    const long double eps_PSI = powl(10, -9);
    long double eps_max_W = eps_W + 1;
    long double eps_max_PSI = eps_PSI + 1;

    // сетка
    const long long N_R = ceil(R / dr) + 1;
    const long long N_Z = ceil(1.0 / dz) + 1;

    // инициализация массивов
    long double **W = new long double *[N_R];
    for (long long i = 0; i < N_R; i++)
        W[i] = new long double[N_Z];

    long double **W_new = new long double *[N_R];
    for (long long i = 0; i < N_R; i++)
        W_new[i] = new long double[N_Z];

    long double **PSI = new long double *[N_R];
    for (long long i = 0; i < N_R; i++)
        PSI[i] = new long double[N_Z];

    long double **PSI_new = new long double *[N_R];
    for (long long i = 0; i < N_R; i++)
        PSI_new[i] = new long double[N_Z];

    // начальные значения
    for (long long i = 0; i < N_R; i++){
        for (long long k = 0; k < N_Z; k++) {
            W[i][k] = 0;
            W_new[i][k] = 0;

            PSI[i][k] = 0;
            PSI_new[i][k] = 0;
        }
    }

    // установление PSI
    long double eps;
    while (eps_max_PSI > eps_PSI){
        // начальные условия W
        eps_max_PSI = 0.0;
        for (long long k = 0; k < N_Z; k++){
            W[0][k] = 2.0/dr/dr * PSI[1][k];
            W[N_R-1][k] = 2.0/dr/dr * PSI[N_R-2][k];

            W_new[0][k] = W[0][k];
            W_new[N_R-1][k] = W[N_R-1][k];
        }

        for (long long i = 0; i < N_R; i++){
            W[i][0] = 2.0/dz/dr * PSI[i][1];
            W[i][N_Z-1] = 2.0/dz/dr * PSI[i][N_Z-2];

            W_new[i][0] = W[i][0];
            W_new[i][N_Z-1] = W[i][N_Z-1];
        }

        // установление W
        while (eps_max_W > eps_W/* | t_W < 0.001*/){
            // обновление слоя W
            eps_max_W = 0.0;
            for (long long i = 1; i < N_R-1; i++){
                for (long long k = 1; k < N_Z-1; k++) {
                    W_new[i][k] = W[i][k] + dt*(
                            (W[i-1][k] - 2*W[i][k] + W[i+1][k])/dr/dr +
                            (W[i][k-1] - 2*W[i][k] + W[i][k+1])/dz/dz +
                            1.0/(i*dr)*(W[i+1][k] - W[i-1][k])/(2.0*dr) -
                            W[i][k]*W[i][k]/(i*dr)/(i*dr) +
                            S*(i*dr)*(k*dz - 0.5)
                    );

                    if (W[i][k] == 0)
                        eps = eps_W + 1;
                    else
                        eps = abs((W_new[i][k] - W[i][k])/W[i][k]);

                    if (eps > eps_max_W)
                        eps_max_W = eps;

                    if (abs(W[i][k]) > W_max)
                        W_max = W[i][k];
                }
            }

            for (long long i = 1; i < N_R-1; i++)
                for (long long k = 1; k < N_Z-1; k++)
                    W[i][k] = W_new[i][k];

            t_W += dt;

            //cout << "eps_max_W = " << eps_max_W << endl;
        }
        t_W = 0;
        eps_max_W = eps_W + 1;


        // обновление слоя PSI
        for (long long i = 1; i < N_R-1; i++){
            for (long long k = 1; k < N_Z-1; k++) {
                PSI_new[i][k] = PSI[i][k] + dt*(
                        (PSI[i-1][k] - 2*PSI[i][k] + PSI[i+1][k])/dr/dr +
                        (PSI[i][k-1] - 2*PSI[i][k] + PSI[i][k+1])/dz/dz +
                        1.0/(i*dr)*(PSI[i+1][k] - PSI[i-1][k])/(2.0*dr) -
                        PSI[i][k]*PSI[i][k]/(i*dr)/(i*dr) +
                        W[i][k]
                );

                if (PSI[i][k] == 0)
                    eps = eps_PSI + 1;
                else
                    eps = abs((PSI_new[i][k] - PSI[i][k])/PSI[i][k]);

                if (eps > eps_max_PSI)
                    eps_max_PSI = eps;

                if (abs(PSI[i][k]) > PSI_max)
                    PSI_max = PSI[i][k];
            }
        }

        for (long long i = 1; i < N_R-1; i++)
            for (long long k = 1; k < N_Z-1; k++)
                PSI[i][k] = PSI_new[i][k];

        t_PSI += dt;

        // вывод промежуточных результатов в файл
        if (t_PSI - last_time_print > time_step){
            last_time_print = t_PSI;
            cout << "t_PSI = " << t_PSI << "\teps_max_PSI = " << eps_max_PSI << "\t|PSI_max| = " << PSI_max << "\t|W_max| = " << W_max << endl;

            stringstream ss;
            ss << std::fixed << setprecision(filename_precision);
            ss << last_time_print;

            ofstream outFile_csv("out" + ss.str() + ".csv");
            for (long long i = 0; i < N_R; i++)
                for (long long k = 0; k < N_Z; k++)
                    outFile_csv << i * dr << ";" << k * dz << ";" << PSI[i][k] << ";" << W[i][k] << endl;
            outFile_csv.close();


            ofstream outFile_txt("out" + ss.str() + ".txt");
            for (long long i = 0; i < N_R; i++)
                for (long long k = 0; k < N_Z; k++)
                    outFile_txt << i * dr << "\t" << k * dz << "\t" << PSI[i][k] << "\t" << W[i][k] << endl;
            outFile_txt.close();

        }
    }

    // вывод в файл последний слой PSI и W
    ofstream outFile_csv("out_last.csv");
    for (long long i = 0; i < N_R; i++)
        for (long long k = 0; k < N_Z; k++)
            outFile_csv << i * dr << ";" << k * dz << ";" << PSI[i][k] << ";" << W[i][k] << endl;
    outFile_csv.close();


    ofstream outFile_txt("out_last.txt");
    for (long long i = 0; i < N_R; i++)
        for (long long k = 0; k < N_Z; k++)
            outFile_txt << i * dr << "\t" << k * dz << "\t" << PSI[i][k] << "\t" << W[i][k] << endl;
    outFile_txt.close();


    return 0;
}
