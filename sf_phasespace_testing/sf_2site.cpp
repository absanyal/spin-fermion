#include <iostream>
#include <complex>
#include <vector>
#include "Matrix.h"
#include "diag.h"
#include <cmath>
#include <fstream>

using std::cout;
using std::endl;
using std::pair;
using std::vector;

int Lx, Ly, numw, Mx, My, nx, ny;
float minw, maxw, JH, G;
int k(int x, int y) { return x + Lx * y; }
pair<int, int> kinv(int);
double lorentzian(double, double);
cd phasex(int);
cd phasey(int);

double d_theta_x = 0.0;
double d_theta_y = 0.0;

cd tx0 = cd(1.0, 0.0);
cd ty0 = cd(1.0, 0.0);
cd imagi = cd(0, 1);
cd t = cd(1.0, 0.0);

cd phasex(int nx)
{
    return exp(imagi * M_2_PI * (nx * 1.0) / (1.0 * Mx));
}

cd phasey(int ny)
{
    return exp(imagi * M_2_PI * (ny * 1.0) / (1.0 * My));
}

pair<int, int> kinv(int M)
{
    int x = (int)M / Lx;
    int y = M % Lx;
    return std::make_pair(x, y);
}

double lorentzian(double x, double x0)
{
    return (1 / M_PI) * 0.5 * G / (pow(x - x0, 2) + pow(0.5 * G, 2));
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "JH, num of pts"
                  << endl;
        exit(1);
    }
    // Lx = std::atoi(argv[1]);
    // Ly = std::atoi(argv[2]);
    // Mx = std::atoi(argv[3]);
    // My = std::atoi(argv[4]);
    // nx = std::atoi(argv[5]);
    // ny = std::atoi(argv[6]);
    // minw = std::atof(argv[7]);
    // maxw = std::atof(argv[8]);
    JH = std::atof(argv[1]);
    numw = std::atoi(argv[2]);
    // G = std::atof(argv[11]);

    double theta, phi;

    // if (Lx < 2 or Ly < 2)
    // {
    //     std::cerr << "Lattice cannot be 1D." << endl;
    //     exit(1);
    // }

    std::ofstream outfile;
    outfile.open("2sitedata.dat");

    // Create a 2*N*N x 2*N*N matrix
    Matrix<cd> H;
    H.resize(4, 4);

    for (int i = 0; i < numw; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            theta = i * 2 * 3.1415 / (numw * 1.0);
            phi = j * 2 * 3.1415 / (numw * 1.0);
            cout << theta << "\t" << phi << endl;

            H(0, 0) = -0.5 * JH;
            H(1, 1) = -0.5 * JH * cos(theta);
            H(2, 2) = 0.5 * JH;
            H(3, 3) = 0.5 * JH * cos(theta);
            H(0, 1) = t;
            H(2, 3) = t;
            H(1, 3) = 0.5 * JH * sin(theta) * exp(-imagi * phi);

            vector<double> eigs_;
            Diagonalize('V', H, eigs_);

            outfile << theta << "\t" << phi << "\t" << eigs_[0] << "\t" << eigs_[1] << "\t" << eigs_[2] << "\t" << eigs_[3] << endl;
        }
        // outfile << endl;
    }

    outfile.close();

    return 0;
}