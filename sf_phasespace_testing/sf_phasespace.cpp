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

int Lx, Ly, numw, Mx, My, nx, ny, thetanum, numparticles;
float minw, maxw, JH, G;
int k(int x, int y) { return x + Lx * y; }
pair<int, int> kinv(int);
double lorentzian(double, double);
cd phasex(int);
cd phasey(int);

cd tx0 = cd(1.0, 0.0);
cd ty0 = cd(1.0, 0.0);
cd t = cd(1.0, 0.0);
cd imagi = cd(0, 1);

cd phasex(int nx)
{
    return exp(imagi * M_2_PI * (nx * (Mx - 1) * 1.0) / (1.0 * Mx));
}

cd phasey(int ny)
{
    return exp(imagi * M_2_PI * (ny * (My - 1) * 1.0) / (1.0 * My));
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
    if (argc != 6)
    {
        std::cerr << "Enter Lx, Ly, JH, thetanum, numparticles"
                  << endl;
        exit(1);
    }
    Lx = std::atoi(argv[1]);
    Ly = std::atoi(argv[2]);
    JH = std::atof(argv[3]);
    thetanum = std::atoi(argv[4]);
    numparticles = std::atoi(argv[5]);

    if (Lx < 2 or Ly < 2)
    {
        std::cerr << "Lattice cannot be 1D." << endl;
        exit(1);
    }

    std::ofstream outfile;
    outfile.open("2D_phasepace.dat");
    Matrix<cd> H;
    H.resize(2 * Lx * Ly, 2 * Lx * Ly);

    Matrix<cd> theta;
    theta.resize(1, Lx * Ly);
    Matrix<cd> phi;
    phi.resize(1, Lx * Ly);

    for (int theta_c = 0; theta_c < thetanum; theta_c++)
    {

        // Theta and phi arrays
        double d_theta_x = theta_c * 2.0 * 3.1415 / (thetanum * 1.0);
        double d_theta_y = theta_c * 2.0 * 3.1415 / (thetanum * 1.0);
        double d_phi_x = 0.0;
        double d_phi_y = 0.0;

        for (int i = 0; i < Lx; i++)
        {
            for (int j = 0; j < Ly; j++)
            {
                int pos = k(i, j);
                theta(0, pos) = i * d_theta_x + j * d_theta_y;
                phi(0, pos) = i * d_phi_x + j * d_phi_y;
            }
        }

        // Populate TB part
        for (int j = 0; j < Ly; j++)
        {
            for (int i = 0; i < Lx; i++)
            {
                // Coordinate mapping
                pair<int, int> h1 = std::make_pair(i, (j + 1) % Ly);      //y+1
                pair<int, int> h2 = std::make_pair(i, (j - 1 + Ly) % Ly); //y-1
                pair<int, int> h3 = std::make_pair((i + 1) % Lx, j);      // x+1
                pair<int, int> h4 = std::make_pair((i - 1 + Lx) % Lx, j); //x-1

                int m = k(i, j);
                int n1 = k(h1.first, h1.second);
                int n2 = k(h2.first, h2.second);
                int n3 = k(h3.first, h3.second);
                int n4 = k(h4.first, h4.second);

                // Fill TB elemnts
                H(n1, m) = t;
                H(n2, m) = t;
                H(n3, m) = t;
                H(n4, m) = t;

                H(n1 + Lx * Ly, m + Lx * Ly) = t;
                H(n2 + Lx * Ly, m + Lx * Ly) = t;
                H(n3 + Lx * Ly, m + Lx * Ly) = t;
                H(n4 + Lx * Ly, m + Lx * Ly) = t;
            }
        }

        // Hund's terms
        for (int i = 0; i < Lx * Ly; i++)
        {
            H(i, i) = -0.5 * JH * cos(theta(0, i));
            H(i + Lx * Ly, i + Lx * Ly) = 0.5 * JH * cos(theta(0, i));
            H(i, i + Lx * Ly) =
                -0.5 * JH * sin(theta(0, i)) * exp(-imagi * phi(0, i));
        }

        // H.print();

        vector<double> eigs_;
        Diagonalize('V', H, eigs_);
        // cout << "Finished diagonalization." << endl;

        outfile << d_theta_x;

        for (int i = 0; i < eigs_.size(); i++)
        {
            outfile << "\t" << eigs_[i];
        }

        // int filledlevels = (int)numparticles / 2 + numparticles % 2;
        double totale = 0;
        for (int i = 0; i < numparticles; i++){
            totale += eigs_[i];
        }

        outfile << "\t" << totale;
        outfile << endl;
    }

    outfile.close();
}