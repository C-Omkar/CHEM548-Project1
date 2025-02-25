#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
#include <Eigen/Dense> // install eigen3 library first
#include <cmath>

using namespace Eigen;
using namespace std;

const double hbar {1.05457180013e-34};
const double electron_mass {9.10938356e-31};
const double evtoJ {1.60217662e-19};
const double jtoEv {6.241509e18};
const double h {hbar * 2 * M_PI};
const double mtoAng {1e10};
const double angtoM {1e-10};

const double L {10};
const int num_points {100};



int main()
{
    double kinetic, unit_length, analytical_energy, x_grid[num_points], potential[num_points];
    // vector<double> eigenvals(num_points);
    vector<vector<double>> val_and_vec(num_points, std::vector<double>(num_points+1));

    MatrixXd hamiltonian(num_points, num_points);
    hamiltonian.setZero();
    MatrixXcd D, o;

    unit_length = 2*L/(num_points-1);
    kinetic = pow(hbar, 2) / (2 * electron_mass * pow(unit_length, 2)) * jtoEv * pow(mtoAng, 2);

    fstream file1, file2, file3;
    file1.open("finite_potential_hamiltonian.txt", ios::out);
    file2.open("finite_potential_eigenvalue.txt", ios::out);
    file3.open("finite_potential_eigenvector.txt", ios::out);

    for(int i = 0; i < num_points; i++)
    {   
        x_grid[i] = -L + i * unit_length;
        potential[i] = 0;
        if((x_grid[i] <= -0.5*L) || (x_grid[i] >= 0.5*L))
        {
            potential[i] = 8;
        }
        
        hamiltonian(i, i) = 2*kinetic + potential[i];
        if(i != 0)
        {
            hamiltonian(i, i-1) = -kinetic;
            hamiltonian(i-1, i) = -kinetic;
        }
    }

    for (int i = 0; i < num_points; i++)
    {
        for (int j = 0; j < num_points; j++)
        {
            file1 << hamiltonian(i, j) << " ";
        }
        file1 << endl;
    }
    file1 << endl;

    EigenSolver<MatrixXd> eigensolver(hamiltonian);

    D = eigensolver.eigenvalues();
    o = eigensolver.eigenvectors();

    for (int i = 0; i < num_points; i++)
    {
        val_and_vec[i].push_back(D(i, 0).real()); // adding eigenvalues in first column
        for (int j = 0; j < num_points; j++)
        {
            val_and_vec[i].push_back(o(j, i).real()); // adding eigenvectors in n columns (transpose of o)
        }
    }
    sort(val_and_vec.begin(), val_and_vec.end()); // sorting increaing eigenvalues and respective eigenvectors thats why transposed and saved in one row
    for (int i = 0; i < num_points; i++)
    {
        analytical_energy = pow((i+1) * h, 2) / (8 * electron_mass * pow(2*L, 2)) * jtoEv * pow(mtoAng, 2);    // calculating eigenvalues
        file2 << i + 1 << " " << val_and_vec[i][0] << " " << analytical_energy << endl; // value of n, eigenvalue using eigensolver, calculated eigenvalue
        // file3 << x_grid[i] << " ";                                    // grid points
        for (int j = 0; j < num_points; j++)
        {
            file3 << val_and_vec[j][i + 1] << " "; // printing transpose so that each column contains eigenvectos for each n 
        }
        file3 << endl;
    }



    return 0;
}
