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

    double x [n][p];
    std::default_random_engine generator;
    std::normal_distribution <float> nd(0.0, 2.0);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){
            x[i][j] = nd(generator);
        }
    }
       
}
