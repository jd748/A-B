#include <random>
#include "Eigen/Eigen"
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "DP.hh"

int main()
{
    int p = 5;
    int n = 100;

    std::mt19937 generator1 (61245);
    std::mt19937 generator2 (16746);
    std::mt19937 generator3 (27351);
    std::mt19937 generator4 (66459);
    std::mt19937 generator5 (12612);
    std::normal_distribution <double> nd(0.0, 1.0);
    std::chi_squared_distribution <double> xs(p - 2);
    std::uniform_real_distribution <double> u(0.0, 1.0);

    //s is Cholesky factorization of Covar matrix
    //mu is vector of patient attribute means
    double s [p][p];
    double zee [n][p];
    Eigen::MatrixXf Z2(n,p);
    double Z[n][p];
    double mu [p];

    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){
            zee[i][j] = nd(generator1);
        }
    }

    for (int i = 0; i < p; i++){
        for (int j = 0; j < p; j++){
            if (i==j)
            {
                s[i][j] = 1.0;
            }
            else
            {
                s[i][j] = 0.1;
            }
        }
    }

    double temp = 0;

    //Z is matrix of patient attributes (n x p)
    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){ 
            temp = 0;

            for (int k = 0; k < p; k++){ 
                temp += s[k][j]*zee[i][k];
            }
                Z[i][j] = temp + mu[j];
                Z2(i,j) = Z[i][j];
        } 
    }

    //Eta + Xi computation
    
    double eta [n];
    double xi [n];

    for (int i = 0; i < n; i++){
        eta[i] = nd(generator3);
        xi[i] = xs(generator4);
    }

    //Solving 2-D DP using mesh
    double M = 2000.0;
    double delta = 0.2;
    int steps = ceil(M/delta)+1;
    int N = 10000;

    DP table(n, delta, M);

    double minus;
    double plus;
    double lambda;
    int m;
    int r;


    //Boundary condition
    for (int i = 0; i < 2*n + 1; i++){
        for (int j = 0; j < steps;j++){
            m = i-n;
            lambda = delta*j;
            table.set(0, i, j, pow(m, 2) + lambda);
        }
    }
    
    //Solving mesh
    for (int l = 1;l < n; l++){
        std::cout << l << "\n";
        for (int i = l; i < 2*n+1-l; i++){
            for (int j = 0; j < steps; j++){
                temp = 0;
                table.set(l,i,j,i*l*j);
                std::cerr<<table.at(l,i,j);
            }
        }
    }

    temp = 0;

    //Vs Naive random allocation
    //Eff = x^T P_Z x where x is allocations, P_Z = I - Z(Z^T Z)^(-1) Z^T
    
    Eigen::MatrixXf PZ;
    Eigen::MatrixXf I = Eigen::MatrixXf::Identity(n,n);
    Eigen::MatrixXf Z2I;
    Z2I = Z2.transpose()*Z2;
    PZ = I - Z2*(Z2I.inverse())*Z2.transpose();

    int rand_x [n]; 
    for (int i = 0; i < n; i++){
        rand_x[i] = -1 + 2*floor(2*i/n);
    }

    std::vector<int> myvector (rand_x, rand_x + n);

    Eigen::VectorXf random_x(n); 
    for (int i = 0; i < N; i++){
        std::random_shuffle ( myvector.begin(), myvector.end());
        for (int j = 0; j < n; j++){
            random_x(j) = myvector[j];
        }
        temp += random_x.transpose() * PZ * random_x;
    }
    double eff_r = temp/N;
    std::cout << eff_r;
}

