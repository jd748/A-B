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
    std::normal_distribution <double> nd(0.0, 1.0);

    double s [n][n];
    double zee [n][p];

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


       
}
