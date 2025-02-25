#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
#include <Eigen/Dense> // install eigen3 library first
#include <cmath>

using namespace Eigen;
using namespace std;

// bool comp(vector<std::complex<double>> a, vector<std::complex<double>> b)
// {

//     return b[0].real() >  a[0].real();
// } 

int main()
{

    int n = 500, lwork, info;
    double t, hbar = 1, m = 0.5, dx, l = 5, e, x[n], v[n];
    vector<vector<double>> evec(n); // vector for storing eigenvalues and eigenvectors
    MatrixXd h(n, n);           // hamiltonian matrix
    MatrixXcd D, o;                 // complex matrices for eigenvalues and eigenvectors

    dx = l / real(n - 1);
    t = pow(hbar, 2) / (2 * m * pow(dx, 2));

    fstream file1, file2, file3;
    file1.open("square_well_hamiltonian.txt", ios::out);
    file2.open("square_well_eigenvalue.txt", ios::out);
    file3.open("square_well_eigenvector.txt", ios::out);

    for (int i = 0; i < n; i++)
    {
        x[i] = i * dx;
        if((i * dx <= 0.2*l) || (i * dx >=0.8*l))
        // if((i * dx >=4))sq   
        {
            v[i] = 8;
        } else{

            v[i] = 0;
        }
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                h(i, j) = 2 * t + v[i];
            else if (abs(i - j) == 1)
                h(i, j) = -t;
            else
                h(i, j) = 0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            file1 << h(i, j) << " ";
        }
        file1 << endl;
    }
    file1 << endl;

    EigenSolver<MatrixXd> es(h);

    D = es.eigenvalues();
    o = es.eigenvectors();

    for (int i = 0; i < n; i++)
    {
        evec[i].push_back(D(i, 0).real()); // adding eigenvalues in first column
        for (int j = 0; j < n; j++)
        {
            evec[i].push_back(o(j, i).real()); // adding eigenvectors in n columns (transpose of o)
        }
    }
    sort(evec.begin(), evec.end()); // sorting increaing eigenvalues and respective eigenvectors thats why transposed and saved in one row
    for (int i = 0; i < n; i++)
    {
        e = pow((i + 1) * 2 * M_PI, 2) / (8 * m * pow(l, 2));    // calculating eigenvalues
        file2 << i + 1 << " " << evec[i][0] << " " << e << endl; // value of n, eigenvalue using eigensolver, calculated eigenvalue
        file3 << x[i] << " ";                                    // grid points
        for (int j = 0; j < n; j++)
        {
            file3 << evec[j][i + 1] << " "; // printing transpose so that each column contains eigenvectos for each n 
        }
        file3 << endl;
    }

    for (int i = 0; i < n; i++)
    {
        cout << v[i] << " ";
    }
    return 0;
}
