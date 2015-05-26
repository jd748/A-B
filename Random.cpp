//Generates instance of random customers
#include <iostream>
#include <random>

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
    std::normal_distribution <double> nd(0.0, 1.0);
    std::chi_squared_distribution <double> xs(p - 2);

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
            s[i][j] = 2*nd(generator2);
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
    M = 2000.0;
    delta = 0.2;
    N = 10000;

    double table [2*n+1][N][n];
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 2*n+1;j++){
            table[j][i][0] = (-n + i)^2 + (0 + delta*j);
        }
    }
    
    for (int k = 1;k < n; k++){
        for (int i = 0; i < N; i++){
            for (int j = 0; j < 2*n + 1; j++){
                                    
            }
        }
    }
}

