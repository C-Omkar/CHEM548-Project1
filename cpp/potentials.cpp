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
const double omega {1e15}; // CHANGE THIS VALUE
const double dissoc_energy {4.747}; // this is in eV
const double morse_width {0.0153};
const double displacement {0};

const double L {10};
const int num_points {300};


int main()
{
    double kinetic, unit_length, analytical_energy, x_grid[num_points];
    unit_length = 2*L/(num_points-1);
    kinetic = pow(hbar, 2) / (2 * electron_mass * pow(unit_length, 2)) * jtoEv * pow(mtoAng, 2);

    unordered_map<std::string, vector<double>> potentials;

    vector<std::string> pot_types = {"infinite_well", "finite_well", "rectangle_barrier", "harmonic_well", "morse_well"};
    fstream file4;
    file4.open("grid.txt", ios::out);

    for (auto& x : pot_types)
    {
        potentials[x] = vector<double>(num_points, 0.0);
    }

    for (int i  = 0; i < num_points; i++)
    {   
        x_grid[i] = -L + i * unit_length;
        file4 << x_grid[i] << " ";                                    // grid points
        for (auto& x : pot_types)
        {
            if (x == "infinite_well")
            {
                continue;
            }
            else if (x == "finite_well")
            {
                if((x_grid[i] <= -0.5*L) || (x_grid[i] >= 0.5*L))
                {
                    potentials[x][i] = 8;
                }
            }
            else if (x == "rectangle_barrier")
            {
                if((x_grid[i] >= -0.5*L) && (x_grid[i] <= 0.5*L))
                {
                    potentials[x][i] = 8;
                }
            }
            else if (x == "harmonic_well")
            {
                potentials[x][i] = (0.5) * electron_mass * pow(omega*x_grid[i], 2);
            }
            else if (x == "morse_well")
            {
                potentials[x][i] = dissoc_energy * pow(1 - exp(-morse_width * (x_grid[i] - displacement)), 2);
            }
        }
    }

    for (auto& pot : pot_types)
    {
        vector<vector<double>> val_and_vec(num_points, std::vector<double>(num_points+1));

        MatrixXd hamiltonian(num_points, num_points);
        hamiltonian.setZero();
        MatrixXcd D, o;

        fstream file1, file2, file3;
        file1.open(pot + "_hamiltonian.txt", ios::out);
        file2.open(pot + "_eigenvalue.txt", ios::out);
        file3.open(pot + "_eigenvector.txt", ios::out);

        for(int i = 0; i < num_points; i++)
        {            
            hamiltonian(i, i) = 2*kinetic + potentials[pot][i];
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
            val_and_vec[i][0] = (D(i, 0).real()); // adding eigenvalues in first column
            for (int j = 0; j < num_points; j++)
            {
                val_and_vec[i][j+1] = (o(j, i).real()); // adding eigenvectors in n columns (transpose of o)
            }
        }
        sort(val_and_vec.begin(), val_and_vec.end()); // sorting increaing eigenvalues and respective eigenvectors thats why transposed and saved in one row
        for (int i = 0; i < num_points; i++)
        {
            analytical_energy = pow((i+1) * h, 2) / (8 * electron_mass * pow(2*L, 2)) * jtoEv * pow(mtoAng, 2);    // calculating eigenvalues
            file2 << i + 1 << " " << val_and_vec[i][0] << " " << analytical_energy << endl; // value of n, eigenvalue using eigensolver, calculated eigenvalue
            for (int j = 0; j < num_points; j++)
            {
                file3 << val_and_vec[j][i + 1] << " "; // printing transpose so that each column contains eigenvectos for each n 
            }
            file3 << endl;
        }
        file1.close();
        file2.close();
        file3.close();

    }
    

    file4.close();

    return 0;
}
