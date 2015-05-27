//Generates instance of random customers
#include <iostream>
#include <random>
#include "Eigen/Eigen"
#include <math.h>
#include <algorithm>
#include <vector>

int main()
{
    int p;
    std::cout << "\n Please enter dim of covariates \n";
    std::cin >> p;

    int n;
    std::cout << "\n Please enter number of patients n \n";
    std::cin >> n;

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
            //s[i][j] = 2*nd(generator2);
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
    int N = 10000;

    double table [2*n+1][N][n];
    double minus;
    double plus;
    double lambda;
    double temp;
    double round;

    //Boundary condition
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 2*n+1;j++){
            table[j][i][0] = (-n + i)^2 + (0 + delta*j);
        }
    }
    
    for (int l = 1;l < n; k++){
        for (int i = 0; i < N; i++){
            for (int j = 0; j < 2*(n-l) + 1; j++){
                temp = 0;
                for (int w = 0; w <n; w++){
                    round = floor(((sqrt(i)-eta[w])^2 + xi[w]+ delta/2)/delta)+N/2;
                    if (round > M)
                        round = M;
                    minus = table[j-1][round][l];
                    round = floor(((sqrt(i)+eta[w])^2 + xi[w] + delta/2)/delta) + N/2;
                    if (round > M)
                        round = M;
                    plus = table[j+1][round][l];
                    if(minus<plus){
                        temp = temp + minus;
                    } 
                    else{
                        temp = temp + plus;
                    }          
                }
            }
        }
    }

    
    //Vs Naive random allocation
    //Eff = x^T P_Z x where x is allocations, P_Z = I - Z(Z^T Z)^(-1) Z^T
    double e_rand = n/ 
    

}

