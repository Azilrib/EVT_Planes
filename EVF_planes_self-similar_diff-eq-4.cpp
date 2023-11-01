#define PRECISION 20

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <filesystem>

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return std::move(out).str();
}

long double getStepTimeForOut(long double *time){
    return pow(*time, 1.15)/10.0;
}

int main(){
    //params
    int N = 500; // число интервалов по X
    long double dx = 1.0/N;
    long double S = pow(10.0, 5.0);
    long double dt = pow(dx, 5.0);
    long double T = 1.0;
    string folderName = "data_500/";
    long double stepTimeForOut = 0.0000001;

    MPI_Init(nullptr, nullptr);

    int rank;
    int rankCount;

    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    int *pointsCount = new int[rankCount];
    int *displs = new int[rankCount];
    for (int i = 0; i < rankCount; i++){
        div_t div_mod = div(N+1, rankCount);
        pointsCount[i] = div_mod.quot + (i < div_mod.rem ? 1 : 0);
    }
    displs[0] = 0;
    for (int i = 1; i < rankCount; i++)
        displs[i] = displs[i-1] + pointsCount[i-1];

    int sizeMass = pointsCount[rank] + (rank == 0 || rank == rankCount-1 ? 2 : 4);
    long double a = displs[rank] * dx;
    long double b = (displs[rank] + pointsCount[rank] - 1) * dx;
    long double time = 0.0;
    long double lastTime = time;

    long double *F = new long double [sizeMass];
    long double *newF = new long double [sizeMass];
    long double *finalF = new long double [N + 1];

    for (int i = 0; i < sizeMass; i++)
        F[i] = 0;

    filesystem::create_directories(folderName);
    while (time < T){
        for (int i = 2; i < sizeMass-2; i++)
            newF[i] = F[i] + dt*(
                    -(F[i-2] - 4.0*F[i-1] + 6.0*F[i] - 4.0*F[i+1] + F[i+2])/(dx*dx*dx*dx)
                    - F[i]*(-F[i-2] + 2.0*F[i-1] - 2.0*F[i+1] + F[i+2])/(2.0*dx*dx*dx)
                    + S*(a + (i-2.0)*dx - 0.5)
                    );

        if (rank == 0){
            newF[0] = 0.0;
            newF[1] = 0.0;

            MPI_Sendrecv(&newF[sizeMass-4], 2, MPI_LONG_DOUBLE, 1, 0,
                         &newF[sizeMass-2], 2, MPI_LONG_DOUBLE, 1, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
        }else if (rank > 0 && rank < rankCount - 1){
            MPI_Sendrecv(&newF[2], 2, MPI_LONG_DOUBLE, rank-1, 0,
                         &newF[0], 2, MPI_LONG_DOUBLE, rank-1, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
            MPI_Sendrecv(&newF[sizeMass-4], 2, MPI_LONG_DOUBLE, rank+1, 0,
                         &newF[sizeMass-2], 2, MPI_LONG_DOUBLE, rank+1, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
        }else{
            MPI_Sendrecv(&newF[2], 2, MPI_LONG_DOUBLE, rank-1, 0,
                         &newF[0], 2, MPI_LONG_DOUBLE, rank-1, MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);

            newF[sizeMass-2] = 0.0;
            newF[sizeMass-1] = 0.0;
        }

        for (int i = 0; i < sizeMass; i++)
            F[i] = newF[i];

        if (time - lastTime > getStepTimeForOut(&time)){
            if (rank == 0) {
                MPI_Gatherv(&F[0], sizeMass - 2, MPI_LONG_DOUBLE,
                            finalF, pointsCount, displs, MPI_LONG_DOUBLE,
                            0, MPI_COMM_WORLD);
            }else if (rank > 0 && rank < rankCount - 1){
                MPI_Gatherv(&F[2], sizeMass - 4, MPI_LONG_DOUBLE,
                            nullptr, nullptr, nullptr, MPI_LONG_DOUBLE,
                            0, MPI_COMM_WORLD);
            }else{
                MPI_Gatherv(&F[2], sizeMass - 2, MPI_LONG_DOUBLE,
                            nullptr, nullptr, nullptr, MPI_LONG_DOUBLE,
                            0, MPI_COMM_WORLD);
            }

            if (rank == 0){
                ofstream ofile(folderName + "out_" + to_string_with_precision(time, PRECISION) + ".csv");
                //ofile.precision(PRECISION);
                ofile << fixed;
                for (int i = 0; i < N+1; i++)
                    ofile << i*dx << ";" << finalF[i] << endl;
                ofile.close();
            }

            lastTime = time;
        }

        time += dt;
    }

    //for (int i = 0; i < sizeMass; i++)
        //cout << rank << "\time" << i << "\time" << F[i] << endl;

    if (rank == 0) {
        MPI_Gatherv(&F[0], sizeMass - 2, MPI_LONG_DOUBLE,
                    finalF, pointsCount, displs, MPI_LONG_DOUBLE,
                    0, MPI_COMM_WORLD);
    }else if (rank > 0 && rank < rankCount - 1){
        MPI_Gatherv(&F[2], sizeMass - 4, MPI_LONG_DOUBLE,
                    nullptr, nullptr, nullptr, MPI_LONG_DOUBLE,
                    0, MPI_COMM_WORLD);
    }else{
        MPI_Gatherv(&F[2], sizeMass - 2, MPI_LONG_DOUBLE,
                    nullptr, nullptr, nullptr, MPI_LONG_DOUBLE,
                    0, MPI_COMM_WORLD);
    }

    if (rank == 0){
        ofstream ofile(folderName + "out.csv");
        //ofile.precision(PRECISION);
        ofile << fixed;
        for (int i = 0; i < N+1; i++)
            ofile << i*dx << ";" << finalF[i] << endl;
        ofile.close();
    }

    MPI_Finalize();
    return 0;
}
