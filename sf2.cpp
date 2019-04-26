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
    if (argc != 12)
    {
        std::cerr << "Enter Lx, Ly, Mx, My, nx, ny, minw, "
                     "maxw, numw, JH, Gamma"
                  << endl;
        exit(1);
    }
    Lx = std::atoi(argv[1]);
    Ly = std::atoi(argv[2]);
    Mx = std::atoi(argv[3]);
    My = std::atoi(argv[4]);
    nx = std::atoi(argv[5]);
    ny = std::atoi(argv[6]);
    minw = std::atof(argv[7]);
    maxw = std::atof(argv[8]);
    numw = std::atoi(argv[9]);
    JH = std::atof(argv[10]);
    G = std::atof(argv[11]);

    if (Lx < 2 or Ly < 2)
    {
        std::cerr << "Lattice cannot be 1D." << endl;
        exit(1);
    }

    // Create a 2*N*N x 2*N*N matrix
    Matrix<cd> H;
    H.resize(2 * Lx * Ly, 2 * Lx * Ly);
    // Theta and phi arrays
    Matrix<cd> theta;
    theta.resize(1, Lx * Ly);
    Matrix<cd> phi;
    phi.resize(1, Lx * Ly);

    for (int i = 0; i < Lx; i++)
    {
        for (int j = 0; j < Ly; j++)
        {
            int pos = k(i, j);
            theta(0, pos) = i * d_theta_x + j * d_theta_y;
        }
    }

    // H(1, 5) = cd(1.0, 0.0);
    // H.print();
    cd tempphase1;
    cd tempphase2;
    cd tempphase3;
    cd tempphase4;

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

            // 1
            if (h1.second - j == 1 - Ly)
            {
                tempphase1 = phasey(ny);
            }
            else
            {
                assert(h1.second - j == 1);
                tempphase1 = cd(1.0, 0);
            }
            // 2
            if (h2.second - j == -1 + Ly)
            {
                tempphase2 = std::conj(phasey(ny));
            }
            else
            {
                assert(h2.second - j == -1);
                tempphase2 = cd(1.0, 0);
            }
            // 3
            if (h3.first - i == 1 - Lx)
            {
                tempphase3 = phasex(nx);
            }
            else
            {
                assert(h3.first - i == 1);
                tempphase3 = cd(1.0, 0);
            }
            // 4
            if (h4.first - i == -1 + Lx)
            {
                tempphase4 = std::conj(phasex(nx));
            }
            else
            {
                assert(h4.first - i == -1);
                tempphase4 = cd(1.0, 0);
            }

            // Fill TB elemnts
            H(n1, m) = ty0 * tempphase1;
            H(n2, m) = ty0 * tempphase2;
            H(n3, m) = tx0 * tempphase3;
            H(n4, m) = tx0 * tempphase4;

            H(n1 + Lx * Ly, m + Lx * Ly) = ty0 * tempphase1;
            H(n2 + Lx * Ly, m + Lx * Ly) = ty0 * tempphase2;
            H(n3 + Lx * Ly, m + Lx * Ly) = tx0 * tempphase3;
            H(n4 + Lx * Ly, m + Lx * Ly) = tx0 * tempphase4;
        }
    }

    // H.print();

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
    cout << "Finished diagonalization." << endl;

    std::ofstream outfile;
    // outfile.open("eigenvalues.dat");
    // for (int i = 0; i < eigs_.size(); i++)
    // {
    //     outfile << eigs_[i] << endl;
    // }
    // outfile.close();

    // std::ofstream outfile;

    outfile.open(std::to_string(nx) + '_' + std::to_string(ny) + "_dos.dat");

    int counter = 0;
    for (int xi = 0; xi <= numw; xi++)
    {
        double x = minw + xi * ((maxw - minw) / numw);

        double fullrho = 0;
        double rhoup = 0;
        double rhodown = 0;
        for (int i = 0; i < eigs_.size(); i++)
        {
            fullrho += (1.0 / eigs_.size()) * lorentzian(x, eigs_[i]);
            for (int q = 0; q < Lx * Ly; q++)
            {
                rhoup += std::real(H(q, i) * std::conj(H(q, i))) * (1.0 / eigs_.size()) * lorentzian(x, eigs_[i]);
                rhodown += std::real(H(q + Lx * Ly, i) * std::conj(H(q + Lx * Ly, i))) * (1.0 / eigs_.size()) * lorentzian(x, eigs_[i]);
            }
        }

        outfile << x << '\t' << fullrho << '\t' << rhoup << '\t' << rhodown << endl;
    }

    outfile.close();

    // cout << "Finished calculating full density." << endl;
    cout << 0.5 * (eigs_[Lx * Ly - 1] + eigs_[Lx * Ly]) << endl;
    cout << "Finished" << endl;
    return 0;
}