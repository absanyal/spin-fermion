#include <iostream>
#include <complex>
#include <vector>
#include "Matrix.h"
#include "diag.h"
#include <cmath>
#include <fstream>

using std::cout;
using std::endl;
using std::ifstream;
using std::pair;
using std::vector;

int mx, my;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Enter Mx, My"
                  << endl;
        exit(1);
    }
    mx = std::atoi(argv[1]);
    my = std::atoi(argv[2]);

    vector<double> w_list;
    vector<double> dos_list;
    vector<double> temp_dos_list;

    ifstream file("0.000000_0.000000_0_0_dos.dat");
    double w, dos, dummy;

    while (file >> w >> dos >> dummy >> dummy)
    {
        w_list.push_back(w);
        dos_list.push_back(0);
    }

    file.close();

    for (int nx = 0; nx < mx; nx++)
    {
        for (int ny = 0; ny < my; ny++)
        {
            ifstream file("0.000000_0.000000_" + std::to_string(nx) + "_" + std::to_string(ny) + "_dos.dat");
            while (file >> w >> dos >> dummy >> dummy)
            {
                temp_dos_list.push_back(dos);
            }
            for (int i = 0; i < dos_list.size(); i++)
            {
                dos_list[i] += temp_dos_list[i] / (mx * my * 1.0);
            }
            temp_dos_list.clear();
        }
    }
    file.close();

    std::ofstream outfile;
    outfile.open("avgdos.dat");

    for (int i = 0; i < w_list.size(); i++)
    {
        outfile << w_list[i] << "\t" << dos_list[i] << endl;
    }

    outfile.close();
}