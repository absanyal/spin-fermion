#include "matworks.hpp"

int N;
int k(int x, int y) { return x + N * y; }
std::pair<int, int> kinv(int);
int t = 1;
double JH = 1.0;
cd iota = cd(0, 1);

std::pair<int, int> kinv(int M)
{
    int x = (int)M / N;
    int y = M % N;
    return std::make_pair(x, y);
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Enter N" << endl;
        exit(1);
    }

    N = std::atoi(argv[1]);
    MatrixXcd htbup = MatrixXcd::Zero(N * N, N * N);
    MatrixXcd htbdown = MatrixXcd::Zero(N * N, N * N);
    MatrixXcd spinflip = MatrixXcd::Zero(N * N, N * N);
    MatrixXcd H = MatrixXcd::Zero(2 * N * N, 2 * N * N);
    VectorXd theta = VectorXd::Random(N * N);
    VectorXd phi = VectorXd::Random(N * N);

    // TB for both spin sectors
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            std::pair<int, int> h1 = std::make_pair(i, (j + 1) % N);
            std::pair<int, int> h2 = std::make_pair(i, (j - 1 + N) % N);
            std::pair<int, int> h3 = std::make_pair((i + 1) % N, j);
            std::pair<int, int> h4 = std::make_pair((i - 1 + N) % N, j);

            int m = k(i, j);
            int n1 = k(h1.first, h1.second);
            int n2 = k(h2.first, h2.second);
            int n3 = k(h3.first, h3.second);
            int n4 = k(h4.first, h4.second);

            htbup(m, n1) = t;
            htbup(m, n2) = t;
            htbup(m, n3) = t;
            htbup(m, n4) = t;

            htbdown(m, n1) = t;
            htbdown(m, n2) = t;
            htbdown(m, n3) = t;
            htbdown(m, n4) = t;
        }
    }

    // Hund's coupling z component
    for (int i = 0; i < N * N; i++)
    {
        htbup(i, i) = 0.5 * JH * cos(theta(i));
        htbdown(i, i) = -0.5 * JH * cos(theta(i));
    }

    // Hund's coupling spin flip
    for (int i = 0; i < N * N; i++)
    {
        spinflip(i, i) = 0.5 * JH * sin(theta(i)) * exp(-iota * phi(i));
    }

    H.block(0, 0, N * N, N * N) = htbup;
    H.block(N * N, N * N, N * N, N * N) = htbdown;
    H.block(0, N * N, N * N, N * N) = spinflip;
    H.block(N * N, 0, N * N, N * N) = spinflip.adjoint();

    cout << htbdown << endl;
}