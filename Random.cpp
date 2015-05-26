//Generates instance of random customers
#include <iostream>
#include <random>
#include "Eigen/Eigen"

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
    std::normal_distribution <double> nd(0.0, 1.0);

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


}

