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
    // Patients (n), attributes (p), iterations in solving DP (N1), iterations in estimation of ratio (N2);
    int p = 5;
    int n = 100;
    int N1 = 1000;
    int N2 = 1000;

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
    double mu [p];

    //Generates matrix of standard normals, then uses cholesky factorization to get correct covar matrix
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
            Z2(i,j) = temp + mu[j];
        } 
    }

    //Eta + Xi computation
    
    double eta [N1];
    double xi [N1];

    for (int i = 0; i < N1; i++){
        eta[i] = nd(generator3);
        xi[i] = xs(generator4);
    }

    //Solving 2-D DP using mesh
    double M = 2000.0;
    double delta = 0.2;
    int steps = ceil(M/delta)+1;

    DP table(n, delta, M);

    //syntax helpers
    double lambda = 0;
    int m = 0;

    double lambda_up = 0;
    double lambda_down = 0;
    double roundup_plus = 0;
    double rounddown_plus = 0;
    double roundup_minus = 0;
    double rounddown_minus = 0;

    double distdown = 0;
    double distup = 0;

    double plus = 0;
    double minus = 0;

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
            m = i-n;
            for (int j = 0; j < steps; j++){
                lambda = delta*j;
                temp = 0;
                for (int iter = 0; iter < N1; iter++){ 
                    lambda_up = pow(sqrt(lambda)+eta[iter],2)+xi[iter];
                    roundup_plus = ceil(lambda_up/delta);
                    rounddown_plus = floor(lambda_up/delta); 
                    distdown = lambda_up - rounddown_plus;
                    distup = rounddown_plus - lambda_down;
                    plus = (1/delta)*(distdown*table.at(l-1, i+1, rounddown_plus) + distup*table.at(l-1, i+1, roundup_plus));

                    lambda_down = pow(sqrt(lambda)+eta[iter],2)+xi[iter];
                    roundup_minus = ceil(lambda_down/delta);
                    rounddown_minus = floor(lambda_down/delta);
                    distdown = lambda_up - rounddown_minus;
                    distup = roundup_minus - lambda_up;
                    minus = (1/delta)*(distdown*table.at(l-1,i-1,rounddown_minus) + distup*table.at(l-1,i-1,roundup_minus));

                    temp = temp*(iter-1)/iter + std::min(plus,minus)/iter;
                }

                table.set(l,i,j,temp);          
                
            }
        }
    }

    //Vs Naive random allocation
    //Eff = x^T P_Z x where x is allocations, P_Z = I - Z(Z^T Z)^(-1) Z^T

    Eigen::MatrixXf PZ;
    Eigen::MatrixXf I = Eigen::MatrixXf::Identity(n,n);
    Eigen::MatrixXf Z2I;
    
    //Generates matrix of standard normals, then uses cholesky factorization to get correct covar matrix
     
    for (int i = 0; i < p; i++){
        for (int j = 0; j < p; j++){
            if (i==j) {
                s[i][j] = 1.0; 
            } else {
                s[i][j] = 0.1;
            }
        }
    }

    //Vector of n/2 1's and -1's
    int rand_x[n];
    for (int i = 0; i < n; i++){
        rand_x[i] = -1 + 2*floor(2*i/n);
    }

    std::vector <int> myvector (rand_x, rand_x+n);
    Eigen::VectorXf random_x(n);
    Eigen::VectorXf dp_x(n);

    double eff_r = 0;
    double eff_dp = 0;

    //Large loop to test policies
    for (int asdf = 0; asdf < N2; asdf++){

        for (int i = 0; i < n; i++){
            for (int j = 0; j < p; j++){
                zee[i][j] = nd(generator1);
            }
        }    
                
        //Z is matrix of patient attributes (n x p)  
        for (int i = 0; i < n; i++){
            for (int j = 0; j < p; j++){
                temp = 0;
                for (int k = 0; k < p; k++){
                    temp += s[k][j]*zee[i][k];
                }
                Z2(i,j) = temp + mu[j];
            }
        }    
    
        Z2I = Z2.transpose()*Z2;
        PZ = I - Z2*(Z2I.inverse())*Z2.transpose(); 
    
        temp = 0;
        for (int i = 0; i < N1; i++){
            std::random_shuffle ( myvector.begin(), myvector.end());
            for (int j = 0; j < n; j++){
                random_x(j) = myvector[j];
            }
            temp += random_x.transpose() * PZ * random_x;
        }

        eff_r = temp/N1;

        for (int i = 0; i < n; i++){
        //Need to do argmin procedure
        }

    }
}

